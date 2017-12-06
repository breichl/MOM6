!> Provides the suite of mixing parameterizations via GOTM
module MOM_GOTM

! License goes here?

! Interface to MOM modules
use MOM_checksums,     only : hchksum
use MOM_diag_mediator, only : time_type, diag_ctrl, safe_alloc_ptr, post_data
use MOM_diag_mediator, only : query_averaging_enabled, register_diag_field
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING
use MOM_EOS,           only : EOS_type, calculate_density
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_file_parser,   only : openParameterBlock, closeParameterBlock
use MOM_grid,          only : ocean_grid_type, isPointInCell
use MOM_verticalGrid,  only : verticalGrid_type
! Interface to GOTM modules
use turbulence, only : init_turbulence, do_turbulence
! Point the MOM6 variables to the GOTM variables.
use turbulence, only : G_Kdm=>num, G_Kdt=>nuh, G_Kds=>nus, G_TKE=>tke, G_EPS=>eps
use turbulence, only : G_l => l, G_cde=>cde
use mtridiagonal, only: init_tridiagonal

implicit none ; private

#include "MOM_memory.h"

public :: GOTM_init      ! Allocate the control structure and read inputs.
public :: GOTM_calculate ! Interface MOM to GOTM turbulence.
public :: GOTM_end       ! Gracefully deallocate the control structure.

!/
!> Control structure for containing GOTM parameters/data
type, public :: GOTM_CS ; private

  ! Parameters needed saved for time-step.
  logical :: debug

  real ALLOCABLE_, dimension(NIMEM_,NJMEM_,NKMEM_) :: TKE !< Turbulent kinetic energy
                                                          !  [m2 s-2]
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_,NKMEM_) :: EPS !< Turbulent dissipation rate
                                                          !  [m2 s-3]
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_,NKMEM_) :: L   !< Turbulent mixing length
                                                          !  [m]
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_,NKMEM_) :: KV  !< Turbulent mixing of momentum
                                                          !  [m2 s-1]
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_,NKMEM_) :: KT  !< Turbulent mixing of heat
                                                          !  [m2 s-1]
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_,NKMEM_) :: KS  !< Turbulent mixing of salt
                                                          !  [m2 s-1]
  real ALLOCABLE_, dimension(NIMEMB_,NJMEM_,NKMEM_) :: UO
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_,NKMEM_) :: VO


  ! Diagnostics
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_,NKMEM_) :: N2 !< Buoyancy frequency used by GOTM [s-2]
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_,NKMEM_) :: S2 !< Shear frequency used by GOTM [s-2]


  ! Diagnostic handles and pointers
  type(diag_ctrl), pointer :: diag => NULL()
  integer :: id_TKE = -1
  integer :: id_Dissipation = -1
  integer :: id_LengthScale = -1
  integer :: id_GOTM_KM = -1
  integer :: id_GOTM_KH = -1
  integer :: id_GOTM_KS = -1
  integer :: id_GOTM_S2 = -1
  integer :: id_GOTM_N2 = -1

end type GOTM_CS
!\

! Module data used for debugging only
logical, parameter :: verbose = .False.
#define __DO_SAFETY_CHECKS__

contains

!/
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
  integer :: i, j, k, kgotm

#include "version_variable.h"
  character(len=40) :: mod = 'MOM_GOTM' ! name of this module
  character(len=20) :: string          ! local temporary string
  integer :: namlst=10
  if (associated(CS)) call MOM_error(FATAL, 'MOM_KPP, KPP_init: '// &
           'Control structure has already been initialized')
  allocate(CS)

  ! Read parameters
  call log_version(paramFile, mod, version, 'This is the MOM wrapper to GOTM\n' // &
            'See http://gotm.net')
  call get_param(paramFile, mod, "USE_GOTM", GOTM_init, &
       "If true, turns on the GOTM wrapper to calculate\n"// &
       " turbulent diffusivities. ", default=.false.)
    ! Forego remainder of initialization if not using this scheme
   if (.not. GOTM_init) return

  !Allocate on grid center of computational domain indices
  ! Note the allocated values will not matter, on first time step
  ! completion all values are set to GOTM values.  GOTM handles initialization.
  ALLOC_ (CS%TKE(SZI_(G), SZJ_(G), SZK_(G)+1)) ; CS%TKE(:,:,:) = 0.0
  ALLOC_ (CS%EPS(SZI_(G), SZJ_(G), SZK_(G)+1)) ; CS%EPS(:,:,:) = 0.0
  ALLOC_ (CS%L(SZI_(G), SZJ_(G), SZK_(G)+1))  ; CS%L(:,:,:)    = 0.0
  ALLOC_ (CS%KV(SZI_(G), SZJ_(G), SZK_(G)+1)) ; CS%KV(:,:,:) = 0.0
  ALLOC_ (CS%KS(SZI_(G), SZJ_(G), SZK_(G)+1)) ; CS%KS(:,:,:) = 0.0
  ALLOC_ (CS%KT(SZI_(G), SZJ_(G), SZK_(G)+1)) ; CS%KT(:,:,:) = 0.0
  ! To use KE conserving shear we store old u/v
  ALLOC_ (CS%UO(SZIB_(G), SZJ_(G), SZK_(G))) ; CS%Uo(:,:,:) = 0.0
  ALLOC_ (CS%VO(SZI_(G), SZJB_(G), SZK_(G))) ; CS%Vo(:,:,:) = 0.0

  ! Call GOTM initialization, which reads the gotmturb.nml namelist file
  call init_turbulence (namlst,'gotmturb.nml',G%Ke)
  ! Initialize the GOTM tridiagonal solver, which is used within the GOTM
  ! turbulence evolution and therefore kept distinct from the MOM routine.
  call init_tridiagonal (G%Ke)
  ! Fill the GOTM initial condition to all grid points
  do j = G%jsc, G%jec
    do i = G%isc, G%iec
      do k=1,G%ke+1
        kgotm=G%ke-k+1
        CS%Kv(i,j,k)  = g_Kdm(kgotm)
        CS%Kt(i,j,k)  = g_Kdt(kgotm)
        CS%Ks(i,j,k)  = g_Kds(kgotm)
        CS%TKE(i,j,k) = g_TKE(kgotm)
        CS%EPS(i,j,k) = g_EPS(kgotm)
        CS%L(i,j,k)   = g_l(kgotm)
      enddo
    enddo
  enddo

!/ Prep GOTM related Diagnostics
  CS%diag => diag
  CS%id_TKE = &
       register_diag_field('ocean_model','GOTM_TKE', diag%axesTi, Time, &
                           'Turbulent kinetic energy from GOTM','m2 s-2')
  CS%id_Dissipation = &
       register_diag_field('ocean_model','GOTM_Dissipation', diag%axesTi, Time, &
                           'Turbulent dissipation rate from GOTM','m2 s-3')
  CS%id_LengthScale = &
       register_diag_field('ocean_model','GOTM_L', diag%axesTi, Time, &
                           'Turbulent mixing length from GOTM','m')
  CS%id_GOTM_KM = &
       register_diag_field('ocean_model','GOTM_KM', diag%axesTi, Time, &
                           'Turbulent viscosity of momentum from GOTM','m2 s-1')
  CS%id_GOTM_KH = &
       register_diag_field('ocean_model','GOTM_KH', diag%axesTi, Time, &
                           'Turbulent viscosity of heat from GOTM','m2 s-1')
  CS%id_GOTM_KS = &
       register_diag_field('ocean_model','GOTM_KS', diag%axesTi, Time, &
                           'Turbulent viscosity of salt from GOTM','m2 s-1')
  CS%id_GOTM_S2 = &
       register_diag_field('ocean_model','GOTM_S2', diag%axesTi, Time, &
                           'Shear frequency squared','s-2')
  if (CS%id_GOTM_S2.gt.0) then
    ALLOC_ (CS%S2(SZI_(G), SZJ_(G), SZK_(G)+1)) ; CS%S2(:,:,:) = 0.0
  endif
  CS%id_GOTM_N2 = &
       register_diag_field('ocean_model','GOTM_N2', diag%axesTi, Time, &
                           'Buoyancy frequency squared','s-2')
  if (CS%id_GOTM_N2.gt.0) then
    ALLOC_ (CS%N2(SZI_(G), SZJ_(G), SZK_(G)+1)) ; CS%N2(:,:,:) = 0.0
  endif
  !\

end function GOTM_init
!\

!/
! This subroutine is called to deallocate the control structure
subroutine GOTM_end(CS)
  type(GOTM_CS), pointer :: CS
  deallocate(CS)
end subroutine GOTM_end
!\

!/
!> This subroutine performs the actual stepping of the GOTM turbulence scheme.
!>  NOTE: With GOTM, many flavors of TKE-based vertical mixing parameterizations
!>        can be chosen.  See the gotm documentation at gotm.net for a very
!>        detailed description of the parameterizations.
subroutine GOTM_calculate(CS, G, GV, DT, h, Temp, Salt, u, v, EOS, uStar,&
                          Kt, Ks, Kv)

  ! Arguments
  type(GOTM_CS),                           pointer          :: CS     !< Control structure
  type(ocean_grid_type),                  intent(in)        :: G      !< Ocean grid
  type(verticalGrid_type),                intent(in)        :: GV     !< Ocean vertical grid
  real, intent(in)                                          :: Dt     !< Time step of MOM6 [s] for GOTM turbulence solver.
                                                                      !<   We could subcycle this.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)     :: h      !< Layer/level thicknesses [units of H]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)     :: Temp   !< potential/cons temp at h points [deg C]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)     :: Salt   !< Salinity (ppt) at h points
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)     :: u      !< U velocity at u points [m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)     :: v      !< V velocity at v points [m s-1]
  type(EOS_type),                         pointer           :: EOS    !< Equation of state
  real, dimension(SZI_(G),SZJ_(G)),         intent(in)      :: uStar  !< Surface friction velocity [m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(inout) :: Kt     !< (in)  Vertical diffusivity of heat w/o GOTM [m2 s-1]
                                                                      !< (out) Vertical diffusivity including GOTM [m2 s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(inout) :: Ks     !< (in)  Vertical diffusivity of salt w/o GOTM [m2 s-1]
                                                                      !< (out) Vertical diffusivity including GOTM [m2 s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(inout) :: Kv     !< (in)  Vertical viscosity w/o GOTM [m2 s-1]
                                                                      !< (out) Vertical viscosity including GOTM [m2 s-1]

  ! Local variables
  integer :: i, j, k, km1, kgotm                 ! Loop indices
  real, dimension( 0:G%ke )   :: H_1d            ! 1-d level thickness in [m] for GOTM
  real, dimension( 0:G%ke )   :: N2_1d           ! Brunt-Vaisala frequency squared, at interfaces [s-2]
  real, dimension( 0:G%ke )   :: S2_1d           ! Shear frequency at interfaces [s-2]

  ! for EOS calculation
  !  *These profiles are size 3x the profile, specifically formulated for computing the N2 profile.
  real, dimension( 3*G%ke )   :: rho_1D  ! Density [kg m-3]
  real, dimension( 3*G%ke )   :: pres_1D ! Pressure [Pa]
  real, dimension( 3*G%ke )   :: Temp_1D ! Temperature [deg C]
  real, dimension( 3*G%ke )   :: Salt_1D ! Salinity [PPT]

  real :: GoRho, pRef, Uabove,Ubelow,Vabove,Vbelow
  real :: uoabove, uobelow,voabove,vobelow
  real :: delH                 ! Thickness of a layer [h units]
  real :: Dpt, DZ              ! In-situ depth and dz [h units]
  integer :: kk, ksfc, ktmp

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
      !compute NN and SS on the GOTM vertical grid
      do k=1,G%ke
        ! Note gotm uses a flipped vertical grid compared to MOM.
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
        UOabove  = 0.5*(CS%uo(i,j,k)+CS%uo(i-1,j,k))
        UObelow  = 0.5*(CS%uo(i,j,km1)+CS%uo(i-1,j,km1))
        VOabove  = 0.5*(CS%vo(i,j,k)+CS%vo(i,j-1,k))
        VObelow  = 0.5*(CS%vo(i,j,km1)+CS%vo(i,j-1,km1))

        ! Use KE conserving shear following Burchard et al. (2002)
        ! -> need to check that old/present values are correctly identified.
        S2_1d(kgotm) = .5*( (Uabove-Ubelow)*(UOabove-Ubelow)  + &
                            (Uabove-Ubelow)*(Uabove-UObelow)  + &
                            (Vabove-Vbelow)*(VOabove-Vbelow)  + &
                            (Vabove-Vbelow)*(Vabove-VObelow)  ) &
                            / (dz*dz)

      enddo
      N2_1d(0) = 0.0
      S2_1d(0) = 0.0

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
      ! xP - TKE production due to [perhaps] seagrass [m2/s3]
      !   subroutine do_turbulence(nlev,dt,depth,u_taus,u_taub,z0s,z0b,h,      &
      !                            NN,SS,xP)
      !N2_1d=0.0;S2_1d=0.0;
      do k=1,G%ke+1
        kgotm=G%ke-k+1
        !Fill GOTM arrays for time stepping
        g_Kdm(kgotm) = CS%Kv(i,j,k)
        g_Kdt(kgotm) = CS%Kt(i,j,k)
        g_Kds(kgotm) = CS%Ks(i,j,k)
        g_TKE(kgotm) = CS%TKE(i,j,k)
        g_EPS(kgotm) = CS%EPS(i,j,k)
        g_l(kgotm)   = CS%L(i,j,k)
      enddo
      call do_turbulence(G%ke,&! Number of levels
                           dt,&! time step
                          dpt,&! depth
                   ustar(i,j),&! ustar surface
                           0.,&! ustar bottom
                         0.02,&! Z0 surface ! This could be made to compute from Charnock, or run-time
                         0.02,&! Z0 bottom  ! ditto
                         H_1d,&! Thickness (e.g., see calculation above)
                        N2_1d,&! Buoyancy Frequency
                        S2_1d &! Shear squared
                            ) ! An additional TKE production source is optional input, neglected here
      do k=1,G%ke+1
        kgotm=G%ke-k+1
        !Fill MOM arrays from GOTM
        CS%Kv(i,j,k) = g_Kdm(kgotm)
        CS%Kt(i,j,k) = g_Kdt(kgotm)
        CS%Ks(i,j,k) = g_Kds(kgotm)
        Kv(i,j,k) = g_Kdm(kgotm)
        Kt(i,j,k) = g_Kdt(kgotm)
        Ks(i,j,k) = g_Kds(kgotm)
        CS%TKE(i,j,k) = g_TKE(kgotm)
        CS%EPS(i,j,k) = g_EPS(kgotm)
        CS%L(i,j,k) = g_l(kgotm)
        ! Writing for Output the shear/buoyancy frequencies
        if (CS%id_GOTM_S2.gt.0) then
          CS%S2(i,j,k) = S2_1d(kgotm)
        endif
        if (CS%id_GOTM_N2.gt.0) then
          CS%N2(i,j,k) = N2_1d(kgotm)
        endif
      enddo
    enddo ! i
  enddo ! j

  ! Set old values (needed for KE conserving S2 calculation)
  CS%UO(:,:,:) = U(:,:,:)
  CS%VO(:,:,:) = V(:,:,:)

!/ Diagnostics
  if (CS%id_TKE.gt.0)         call post_data(CS%id_TKE        , CS%TKE, CS%diag)
  if (CS%id_Dissipation.gt.0) call post_data(CS%id_Dissipation, CS%EPS, CS%diag)
  if (CS%id_LengthScale.gt.0) call post_data(CS%id_LengthScale, CS%L  , CS%diag)
  if (CS%id_GOTM_KM.gt.0)     call post_data(CS%id_GOTM_KM    , CS%Kv , CS%diag)
  if (CS%id_GOTM_KH.gt.0)     call post_data(CS%id_GOTM_KH    , CS%Kt , CS%diag)
  if (CS%id_GOTM_KS.gt.0)     call post_data(CS%id_GOTM_KS    , CS%Ks , CS%diag)
  if (CS%id_GOTM_N2.gt.0)     call post_data(CS%id_GOTM_N2    , CS%N2 , CS%diag)
  if (CS%id_GOTM_S2.gt.0)     call post_data(CS%id_GOTM_S2    , CS%S2 , CS%diag)
!\

#ifdef __DO_SAFETY_CHECKS__
  if (CS%debug) then
!    call hchksum(Ks, "KPP out: Ks",G,haloshift=0)
  endif
#endif

end subroutine GOTM_calculate

end module MOM_GOTM
