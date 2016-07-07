!> Provides the suite of mixing parameterizations via GOTM
module MOM_GOTM

! License goes here?

use MOM_coms,          only : max_across_PEs
use MOM_checksums,     only : hchksum, is_NaN
use MOM_diag_mediator, only : time_type, diag_ctrl, safe_alloc_ptr, post_data
use MOM_diag_mediator, only : query_averaging_enabled, register_diag_field
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_PE
use MOM_EOS,           only : EOS_type, calculate_density
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_file_parser,   only : openParameterBlock, closeParameterBlock
use MOM_grid,          only : ocean_grid_type, isPointInCell
use MOM_verticalGrid,  only : verticalGrid_type
use MOM_wave_interface, only : Wave_parameters_CS
use turbulence, only : init_turbulence, do_turbulence
use turbulence, only : G_Kdm=>num, G_Kdt=>nuh, G_Kds=>nus, G_TKE=>tke, G_EPS=>eps, G_TKEo =>tkeo
use turbulence, only : G_l => l, G_cde=>cde, g_k_min=>k_min, g_eps_min=>eps_min
use mtridiagonal, only: init_tridiagonal
implicit none ; private

#include "MOM_memory.h"

public :: GOTM_init
public :: GOTM_calculate
!public :: GOTM_end

! Enumerated constants
integer, private, parameter :: NLT_SHAPE_CVMIX     = 0 !< Use the CVmix profile

!> Control structure for containing KPP parameters/data
type, public :: GOTM_CS ; private

  ! Parameters
  real    :: Ri_crit                   !< Critical bulk Richardson number (defines OBL depth)
  real    :: vonKarman                 !< von Karman constant (dimensionless)
  logical :: debug

  real ALLOCABLE_, dimension(NIMEM_,NJMEM_,NKMEM_) :: TKE, TKEo
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_,NKMEM_) :: EPS, LenScale
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_,NKMEM_) :: KV
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_,NKMEM_) :: KT
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_,NKMEM_) :: KS


end type GOTM_CS

! Module data used for debugging only
logical, parameter :: verbose = .False.
#define __DO_SAFETY_CHECKS__

contains

!> Initialize the GOTM module and set up diagnostics
!! Returns True if GOTM is to be used, False otherwise.
logical function GOTM_init(paramFile, G, diag, Time, CS, passive)

  ! Arguments
  type(param_file_type),   intent(in)    :: paramFile !< File parser
  type(ocean_grid_type),   intent(in)    :: G         !< Ocean grid
  type(diag_ctrl), target, intent(in)    :: diag      !< Diagnostics
  type(time_type),         intent(in)    :: Time      !< Time
  type(GOTM_CS),            pointer       :: CS        !< Control structure
  logical, optional,       intent(out)   :: passive   !< Copy of %passiveMode

  ! Local variables
  integer :: is, ie, js, je, isd, ied, jsd, jed, nz

#include "version_variable.h"
  character(len=40) :: mod = 'MOM_GOTM' ! name of this module
  character(len=20) :: string          ! local temporary string
  integer :: namlst=10
  if (associated(CS)) call MOM_error(FATAL, 'MOM_KPP, KPP_init: '// &
           'Control structure has already been initialized')
  allocate(CS)

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  ! Read parameters
  call log_version(paramFile, mod, version, 'This is the MOM wrapper to GOTM\n' // &
            'See http://gotm.net')
  call get_param(paramFile, mod, "USE_GOTM", GOTM_init, &
                 "If true, turns on the GOTM wrapper\n"// &
                 "to calculate diffusivities. ",     &
                 default=.false.)
  ! Forego remainder of initialization if not using this scheme
  if (.not. GOTM_init) return
  call init_turbulence (namlst,'gotmturb.nml',G%Ke)
  call init_tridiagonal (G%Ke)
  !allocation on grid center of computational domain indices
  ALLOC_ (CS%TKE(is:Ie,js:je,nz+1)) ; CS%TKE(:,:,:) = g_k_min
  ALLOC_ (CS%TKEo(is:Ie,js:je,nz+1)) ; CS%TKEo(:,:,:) = g_k_min
  ALLOC_ (CS%EPS(is:Ie,js:je,nz+1)) ; CS%EPS(:,:,:) = g_eps_min
  ALLOC_ (CS%LenScale(is:Ie,js:je,nz+1)) ; CS%LenScale(:,:,:) = g_cde*g_k_min**1.5/g_eps_min
  
  ALLOC_ (CS%KV(is:Ie,js:je,nz+1)) ; CS%KV(:,:,:) = 1.E-6
  ALLOC_ (CS%KS(is:Ie,js:je,nz+1)) ; CS%KS(:,:,:) = 1.E-6
  ALLOC_ (CS%KT(is:Ie,js:je,nz+1)) ; CS%KT(:,:,:) = 1.E-6

end function GOTM_init



!> KPP vertical diffusivity/viscosity and non-local tracer transport
subroutine GOTM_calculate(CS, G, GV, DT, h, Temp, Salt, u, v, EOS, uStar,&
                          Kt, Ks, Kv, Waves)

  ! Arguments
  type(GOTM_CS),                           pointer       :: CS             !< Control structure
  type(ocean_grid_type),                  intent(in)    :: G              !< Ocean grid
  type(verticalGrid_type),                intent(in)    :: GV             !< Ocean vertical grid
  real, intent(in)                                         :: Dt             !< Time step of MOM6 [s] for GOTM turbulence solver
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: h              !< Layer/level thicknesses (units of H)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: Temp           !< potential/cons temp (deg C)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: Salt           !< Salinity (ppt)
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)    :: u              !< Velocity i-component (m/s)
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)    :: v              !< Velocity j-component (m/s)
  type(EOS_type),                         pointer       :: EOS            !< Equation of state
  real, dimension(SZI_(G),SZJ_(G)),         intent(in)    :: uStar          !< Surface friction velocity (m/s)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(inout) :: Kt       !< (in)  Vertical diffusivity of heat w/o KPP (m2/s)
                                                                        !< (out) Vertical diffusivity including KPP (m2/s)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(inout) :: Ks       !< (in)  Vertical diffusivity of salt w/o KPP (m2/s)
                                                                        !< (out) Vertical diffusivity including KPP (m2/s)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(inout) :: Kv       !< (in)  Vertical viscosity w/o KPP (m2/s)
                                                                        !< (out) Vertical viscosity including KPP (m2/s)
  type(Wave_parameters_CS), pointer, optional              :: Waves  !< Surface wave related control structure.

  ! Local variables
  integer :: i, j, k, km1, kgotm                 ! Loop indices
  real, dimension( 0:G%ke )   :: H_1d            ! 1-d level thickness in [m] for GOTM
  real, dimension( 0:G%ke )   :: N2_1d           ! Brunt-Vaisala frequency squared, at interfaces (1/s2)
  real, dimension( 0:G%ke )   :: S2_1d           ! Shear frequency at interfaces (1/s2)
  real, dimension( 0:G%ke )   :: S2x_1d          ! X Shear frequency at interfaces (1/s2)
  real, dimension( 0:G%ke )   :: S2y_1d          ! Y Shear frequency at interfaces (1/s2)
  real, dimension( 0:G%ke )   :: S2xStk_1d       ! X Stokes Shear frequency at interfaces (1/s2)
  real, dimension( 0:G%ke )   :: S2yStk_1d       ! Y Stokes Shear frequency at interfaces (1/s2)
  real, dimension( 0:G%ke, 2) :: Kdiffusivity    ! Vertical diffusivity at interfaces (m2/s)
  real, dimension( 0:G%ke )   :: Kviscosity      ! Vertical viscosity at interfaces (m2/s)


  ! for EOS calculation
  real, dimension( 3*G%ke )   :: rho_1D
  real, dimension( 3*G%ke )   :: pres_1D
  real, dimension( 3*G%ke )   :: Temp_1D
  real, dimension( 3*G%ke )   :: Salt_1D

  real :: GoRho, pRef, Uabove,Ubelow,Vabove,Vbelow,&
          uSabove, uSbelow, vSabove, vSbelow
  real :: delH                 ! Thickness of a layer (m)
  real :: Dpt, DZ
  real :: Z0S
  integer :: kk, ksfc, ktmp
  logical, save :: first=.true.

#ifdef __DO_SAFETY_CHECKS__
  if (CS%debug) then
!    call hchksum(h*GV%H_to_m, "KPP in: h",G,haloshift=0)
  endif
#endif

  ! Initializations:
  H_1d(:) = 0.0
  N2_1d(:) = 0.0
  S2_1d(:) = 0.0

 
  ! some constants
  GoRho = G%g_Earth / GV%Rho0

!$OMP parallel do default(none) shared(G,GV,CS,EOS,uStar,Temp,Salt,u,v,h,GoRho,       &
!$OMP                                  Kt, Ks, KV)
!$OMP                     firstprivate(nonLocalTrans)                                 &
!$OMP                          private(Coriolis,surfFricVel,SLdepth_0d,hTot,surfTemp, &
!$OMP                                  surfHtemp,surfSalt,surfHsalt,surfU,            &
!$OMP                                  surfHu,surfV,surfHv,iFaceHeight,               &
!$OMP                                  pRef,km1,cellHeight,Uk,Vk,deltaU2,             &
!$OMP                                  rho1,rhoK,rhoKm1,deltaRho,N2_1d,N_1d,delH,     &
!$OMP                                  surfBuoyFlux,Ws_1d,Vt2_1d,BulkRi_1d,           &
!$OMP                                  OBLdepth_0d,zBottomMinusOffset,Kdiffusivity,   &
!$OMP                                  Kviscosity,sigma,kOBL,kk,pres_1D,Temp_1D,      &
!$OMP                                  Salt_1D,rho_1D,surfBuoyFlux2,ksfc)

  ! loop over horizontal points on processor
  do j = G%jsc, G%jec
     do i = G%isc, G%iec
        Dpt = 0.0
        PRef = 0.0
        !compute NN, compute SS
        do k=1,G%ke
           kgotm = G%ke-k+1
           delH = h(i,j,k)*GV%H_to_m
           H_1d(k) = delH
           Dpt = Dpt + delH
           
           ! pressure, temp, and saln for EOS
           ! kk+1 = surface fields
           ! kk+2 = k fields
           ! kk+3 = km1 fields
           km1  = max(1, k-1)
           kk   = 3*(k-1)
           pres_1D(kk+1) = pRef
           pres_1D(kk+2) = pRef
           pres_1D(kk+3) = pRef
           Temp_1D(kk+1) = Temp(i,j,1)
           Temp_1D(kk+2) = Temp(i,j,k)
           Temp_1D(kk+3) = Temp(i,j,km1)
           Salt_1D(kk+1) = Salt(i,j,1)
           Salt_1D(kk+2) = Salt(i,j,k)
           Salt_1D(kk+3) = Salt(i,j,km1)
           
           ! pRef is pressure at interface between k and km1.
           ! iterate pRef for next pass through k-loop.
           pRef = pRef + GV%H_to_Pa * h(i,j,k)
           
        enddo ! k-loop finishes
        
        ! compute in-situ density
        call calculate_density(Temp_1D, Salt_1D, pres_1D, rho_1D, 1, 3*G%ke, EOS)
        
        do k = 1, G%ke
           km1 = max(1, k-1)
           kk = 3*(k-1)
           kgotm = G%ke-k+1
           DZ = ((0.5*(h(i,j,km1) + h(i,j,k))+GV%H_subroundoff)*GV%H_to_m)
           
           N2_1d(kgotm)    = (GoRho * (rho_1D(kk+2) - rho_1D(kk+3)) ) / dz
                
           ! C-grid average to get Uk and Vk on T-points.
           ! k is above, k+1 is below (relative to interface k)
           Uabove  = 0.5*(u(i,j,k)+u(i-1,j,k))
           Ubelow  = 0.5*(u(i,j,km1)+u(i-1,j,km1))
           Vabove  = 0.5*(v(i,j,k)+v(i,j-1,k))
           Vbelow  = 0.5*(v(i,j,km1)+v(i,j-1,km1))
           if (present(Waves).and.associated(Waves)) then
              uSabove  = 0.5*(waves%us_x(i,j,k)+waves%us_x(i-1,j,k))
              uSbelow  = 0.5*(waves%us_x(i,j,km1)+waves%us_x(i-1,j,km1))
              vSabove  = 0.5*(waves%us_y(i,j,k)+waves%us_y(i,j-1,k))
              vSbelow  = 0.5*(waves%us_y(i,j,km1)+waves%us_y(i,j-1,km1))
           else
              uSabove = 0.0
              uSbelow = 0.0
              vSabove = 0.0
              vSbelow = 0.0
           endif
           S2_1d(kgotm) = ((Uabove-Ubelow)*(Uabove-Ubelow)  + &
                (Vabove-Vbelow)*(Vabove-Vbelow)) / (dz*dz)
           S2x_1d(kgotm) = (Uabove-Ubelow)*(Uabove-Ubelow) / (dz*dz)
           S2y_1d(kgotm) = (Vabove-Vbelow)*(Vabove-Vbelow) / (dz*dz)
           S2xStk_1d(kgotm) = (uSabove-uSbelow)*(uSabove-uSbelow) / (dz*dz)
           S2yStk_1d(kgotm) = (vSabove-vSbelow)*(vSabove-vSbelow) / (dz*dz)
        enddo
        N2_1d(0) = 0.0
        S2_1d(0) = 0.0
        S2x_1d(0) = 0.0
        S2y_1d(0) = 0.0
        S2xStk_1d(0) = 0.0
        S2yStk_1d(0) = 0.0
        
      !   GOTM do_turbulence format:
      ! nlev - number of levels
      ! dt   - time step [s]
      ! depth - distance between surface and bottom [m]
      ! u_taus - surface friction velocity [m/s]
      ! u_taub - bottom friction velocity [m/s]
      ! z0s - surface roughness [m]
      ! z0b - bottom roughness [m]
      ! h - surface thickness array [m]
      ! NN - buoyancy frequency array [1/s2]
      ! SS - shear freuqnecy array [1/s2]
      ! SSu - shear frequency in x
      ! SSv - shear frequency in y
      ! SSuS - Stokes shear frequency in x
      ! SSvS - Stokes shear frequency in y
      ! xP - TKE production due to [perhaps] seagrass [m2/s3]
      !   subroutine do_turbulence(nlev,dt,depth,u_taus,u_taub,z0s,z0b,h,      &
      !                            NN,SS,xP)
      !N2_1d=0.0;S2_1d=0.0;
        do k=1,G%ke+1
           kgotm=G%ke-k+1
           !Fill GOTM arrays for time stepping
           g_TKE(kgotm) = CS%TKE(i,j,k)
           g_TKEo(kgotm) = CS%TKEo(i,j,k)
           g_EPS(kgotm) = CS%EPS(i,j,k)
           g_Kdm(kgotm) = CS%Kv(i,j,k)
           g_Kdt(kgotm) = CS%Kt(i,j,k)
           g_Kds(kgotm) = CS%Ks(i,j,k)
           g_l(kgotm)   = CS%LenScale(i,j,k)
           !if (k.lt.10) then
           ! print*,k,'in----------'
           ! print*,g_kdm(kgotm),g_Kdt(kgotm),g_Kds(kgotm)
           ! print*,S2_1d(kgotm),N2_1d(kgotm)
           !endif
        enddo
      !Original GOTM do_turbulence
      !call do_turbulence(G%ke,dt,dpt,ustar(i,j),0.,0.0,0.0,H_1d,N2_1d,S2_1d)
      !Modified GOTM do_turbulence for including wave impacts.
        z0s=0.02
        call do_turbulence(G%ke,dt,dpt,ustar(i,j),0.,z0s,0.0,H_1d,N2_1d,S2_1d,&
                           S2x_1d,S2y_1d,S2xStk_1d,S2yStk_1d) 
      do k=1,G%ke+1
         kgotm=G%ke-k+1
         !Fill MOM arrays from GOTM
         if (k.lt.10) then
            !print*,k,'out----------'
            !print*,g_kdm(kgotm),g_Kdt(kgotm),g_Kds(kgotm)
            !print*,S2_1d(kgotm),N2_1d(kgotm)
         endif
         CS%Kv(i,j,k) = g_Kdm(kgotm)
         CS%Kt(i,j,k) = g_Kdt(kgotm)
         CS%Ks(i,j,k) = g_Kds(kgotm)
         Kv(i,j,k) = g_Kdm(kgotm) 
         Kt(i,j,k) = g_Kdt(kgotm) 
         Ks(i,j,k) = g_Kds(kgotm)
         CS%TKEo(i,j,k) = CS%TKE(i,j,k)
         CS%TKE(i,j,k) = g_TKE(kgotm)
         CS%EPS(i,j,k) = g_EPS(kgotm)
         CS%LenScale(i,j,k) = g_L(kgotm)
      enddo
    enddo ! i
  enddo ! j
  
  ! do i=3,4
  !    do j=3,4
  !       print*,i,j,'-----'
  !       print*,u(i,j,1),u(i,j,2)
  !    enddo
  ! enddo
#ifdef __DO_SAFETY_CHECKS__
  if (CS%debug) then
!    call hchksum(Ks, "KPP out: Ks",G,haloshift=0)
  endif
#endif

  first=.false.

end subroutine GOTM_calculate

end module MOM_GOTM
