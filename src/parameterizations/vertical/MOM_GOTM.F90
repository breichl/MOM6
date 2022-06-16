!> Provides the suite of mixing parameterizations via GOTM
module MOM_GOTM

! License goes here?
! Interface to MOM modules
use MOM_checksums,     only : hchksum
use MOM_diag_mediator, only : time_type, diag_ctrl, safe_alloc_alloc, post_data
use MOM_diag_mediator, only : query_averaging_enabled, register_diag_field
use MOM_domains,       only : pass_var
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
use turbulence, only : G_B => B, G_P=>P
use mtridiagonal, only: init_tridiagonal

implicit none ; private

#include "MOM_memory.h"

public :: GOTM_init      ! Allocate the control structure and read inputs.
public :: GOTM_calculate ! Interface MOM to GOTM turbulence.
public :: GOTM_calculate_vertex ! Interface MOM to GOTM turbulence.
public :: GOTM_at_vertex
public :: GOTM_end       ! Gracefully deallocate the control structure.

!/
!> Control structure for containing GOTM parameters/data
type, public :: GOTM_CS ; private

  ! Parameters needed saved for time-step.
  logical :: debug

  real, allocatable, dimension(:,:,:) :: TKE, & !< Turbulent kinetic energy [m2 s-2]
                                         EPS, & !< Turbulent dissipation rate [m2 s-3]
                                         L, &   !< Turbulent mixing length [m]
                                         KV, &  !< Turbulent mixing of momentum [m2 s-1]
                                         KT, &  !< Turbulent mixing of heat [m2 s-1]
                                         KS, &  !< Turbulent mixing of salt [m2 s-1]
                                         UO, &  !< U velocity component, previous step [m s-1]
                                         VO     !< V velocity component, previous step [m s-1]

  ! Diagnostics
  real, allocatable, dimension(:,:,:) :: N2, & !< Buoyancy frequency used by GOTM [s-2]
                                         S2, & !< Shear frequency used by GOTM [s-2]
                                         P,  & !< Shear Production by GOTM [m2 s-3]
                                         B     !< Buoyancy Production by GOTM [m2 s-3]

  real :: minK = 1.e-4

  integer :: NSubCycle

  ! Diagnostic handles and pointers
  type(diag_ctrl), pointer :: diag => NULL()
  integer :: id_TKE = -1
  !integer :: id_Dissipation = -1
  integer :: id_LengthScale = -1
  integer :: id_GOTM_KM = -1
  integer :: id_GOTM_KH = -1
  integer :: id_GOTM_KS = -1
  integer :: id_GOTM_S2 = -1
  integer :: id_GOTM_N2 = -1
  integer :: id_GOTM_P = -1
  integer :: id_GOTM_B = -1
  integer :: id_GOTM_Eps = -1


end type GOTM_CS
!\

! Module data used for debugging only
logical, parameter :: verbose = .False.
#define __DO_SAFETY_CHECKS__

contains

!/
!> Initialize the GOTM module and set up diagnostics
!! Returns True if GOTM is to be used, False otherwise.
logical function GOTM_init(paramFile, G, GV, diag, Time, CS, passive)

  ! Arguments
  type(param_file_type),   intent(in)    :: paramFile !< File parser
  type(ocean_grid_type),   intent(in)    :: G         !< Ocean grid
  type(verticalGrid_type), intent(in)    :: GV        !< The ocean's vertical grid structure.
  type(diag_ctrl), target, intent(in)    :: diag      !< Diagnostics
  type(time_type),         intent(in)    :: Time      !< Time
  type(GOTM_CS),            pointer       :: CS        !< Control structure
  logical, optional,       intent(out)   :: passive   !< Copy of %passiveMode

  ! Local variables
  integer :: isd, ied, jsd, jed, isdB, iedB, jsdB, jedB
  integer :: isc, iec, jsc, jec, iscB, iecB, jscB, jecB
  integer :: nz
  integer :: i, j, k, kgotm

#include "version_variable.h"
  character(len=40) :: mod = 'MOM_GOTM' ! name of this module
  character(len=20) :: string          ! local temporary string
  integer :: namlst=10
  logical :: VS

  nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  isdB = G%isdB ; iedB = G%iedB ; jsdB = G%jsdB ; jedB = G%jedB
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  iscB = G%iscB ; iecB = G%iecB ; jscB = G%jscB ; jecB = G%jecB

  if (associated(CS)) call MOM_error(FATAL, 'MOM_GOTM_mixing, GOTM_init: '// &
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
  call get_param(paramFile, mod, "GOTM_NSUBCYCLE", CS%NSubCycle, &
       "Number of times to cycle GOTM per timestep. ", default=1)
  call get_param(paramFile, mod, "GOTM_VERTEX_SHEAR", VS, &
       "Checking Vertex Shear, not logged. ", &
       default=.false.,do_not_log=.true.)

  ! Call GOTM initialization, which reads the gotmturb.nml namelist file
  ! Can we silence this on all but one PE???
  call init_turbulence (namlst,'gotmturb.nml',G%Ke)
  ! Initialize the GOTM tridiagonal solver, which is used within the GOTM
  ! turbulence evolution and therefore kept distinct from the MOM routine.
  call init_tridiagonal (G%Ke)

  !Allocate on grid center of computational domain indices
  ! Note the allocated values will not matter, on first time step
  ! completion all values are set to GOTM values.  GOTM handles initialization.
  if (VS) then
    call safe_alloc_alloc(CS%TKE, isdB, iedB, jsdB, jedB, nz+1) ; CS%TKE(:,:,:) = 0.0
    call safe_alloc_alloc(CS%EPS, isdB, iedB, jsdB, jedB, nz+1) ; CS%EPS(:,:,:) = 0.0
    call safe_alloc_alloc(CS%L, isdB, iedB, jsdB, jedB, nz+1) ; CS%L(:,:,:) = 0.0
    call safe_alloc_alloc(CS%Kv, isdB, iedB, jsdB, jedB, nz+1) ; CS%Kv(:,:,:) = 0.0
    call safe_alloc_alloc(CS%Ks, isdB, iedB, jsdB, jedB, nz+1) ; CS%Ks(:,:,:) = 0.0
    call safe_alloc_alloc(CS%Kt, isdB, iedB, jsdB, jedB, nz+1) ; CS%Kt(:,:,:) = 0.0
    ! To use KE conserving shear we store old u/v (not allocated, not used, needs work)
    !call safe_alloc_alloc(CS%UO, isdB, iedB, jsd, jed, nz) ; CS%UO(:,:,:) = 0.0
    !call safe_alloc_alloc(CS%VO, isd, ied, jsdB, jedB, nz) ; CS%VO(:,:,:) = 0.0
    ! Fill the GOTM initial condition to all grid points
    do j = jsc-1, jecB ; do i = isc-1, iecB
      do k=1,nz+1
        kgotm=nz-k+1
        CS%Kv(i,j,k)  = g_Kdm(kgotm)
        CS%Kt(i,j,k)  = g_Kdt(kgotm)
        CS%Ks(i,j,k)  = g_Kds(kgotm)
        CS%TKE(i,j,k) = g_TKE(kgotm)
        CS%EPS(i,j,k) = g_EPS(kgotm)
        CS%L(i,j,k)   = g_l(kgotm)
      enddo
    enddo ; enddo
  else
    call safe_alloc_alloc(CS%TKE, isd, ied, jsd, jed, nz+1) ; CS%TKE(:,:,:) = 0.0
    call safe_alloc_alloc(CS%EPS, isd, ied, jsd, jed, nz+1) ; CS%EPS(:,:,:) = 0.0
    call safe_alloc_alloc(CS%L, isd, ied, jsd, jed, nz+1) ; CS%L(:,:,:) = 0.0
    call safe_alloc_alloc(CS%Kv, isd, ied, jsd, jed, nz+1) ; CS%Kv(:,:,:) = 0.0
    call safe_alloc_alloc(CS%Ks, isd, ied, jsd, jed, nz+1) ; CS%Ks(:,:,:) = 0.0
    call safe_alloc_alloc(CS%Kt, isd, ied, jsd, jed, nz+1) ; CS%Kt(:,:,:) = 0.0
    ! To use KE conserving shear we store old u/v (not allocated, not used, needs work)
    !call safe_alloc_alloc(CS%UO, isdB, iedB, jsd, jed, nz) ; CS%UO(:,:,:) = 0.0
    !call safe_alloc_alloc(CS%VO, isd, ied, jsdB, jedB, nz) ; CS%VO(:,:,:) = 0.0
    ! Fill the GOTM initial condition to all grid points
    do j = jsc, jec ; do i = isc, iec
      do k=1,nz+1
        kgotm=nz-k+1
        CS%Kv(i,j,k)  = g_Kdm(kgotm)
        CS%Kt(i,j,k)  = g_Kdt(kgotm)
        CS%Ks(i,j,k)  = g_Kds(kgotm)
        CS%TKE(i,j,k) = g_TKE(kgotm)
        CS%EPS(i,j,k) = g_EPS(kgotm)
        CS%L(i,j,k)   = g_l(kgotm)
      enddo
    enddo ; enddo
  endif


!/ Prep GOTM related Diagnostics
  CS%diag => diag
  if (VS) then
    CS%id_TKE = &
         register_diag_field('ocean_model','GOTM_TKE', diag%axesBi, Time, &
                             'Turbulent kinetic energy from GOTM','m2 s-2')
!    CS%id_Dissipation = &
!         register_diag_field('ocean_model','GOTM_Dissipation', diag%axesBi, Time, &
!                             'Turbulent dissipation rate from GOTM','m2 s-3')
    CS%id_LengthScale = &
         register_diag_field('ocean_model','GOTM_L', diag%axesBi, Time, &
                             'Turbulent mixing length from GOTM','m')
    CS%id_GOTM_KM = &
         register_diag_field('ocean_model','GOTM_KM', diag%axesBi, Time, &
                             'Turbulent viscosity of momentum from GOTM','m2 s-1')
    CS%id_GOTM_KH = &
         register_diag_field('ocean_model','GOTM_KH', diag%axesBi, Time, &
                             'Turbulent viscosity of heat from GOTM','m2 s-1')
    CS%id_GOTM_KS = &
         register_diag_field('ocean_model','GOTM_KS', diag%axesBi, Time, &
                             'Turbulent viscosity of salt from GOTM','m2 s-1')
    CS%id_GOTM_S2 = &
         register_diag_field('ocean_model','GOTM_S2', diag%axesBi, Time, &
                             'Shear frequency squared','s-2')
    CS%id_GOTM_N2 = &
         register_diag_field('ocean_model','GOTM_N2', diag%axesBi, Time, &
                             'Buoyancy frequency squared','s-2')
    CS%id_GOTM_B = &
         register_diag_field('ocean_model','GOTM_B', diag%axesBi, Time, &
                             'Buoyancy Production of TKE from GOTM','m2 s-3')
    CS%id_GOTM_P = &
         register_diag_field('ocean_model','GOTM_P', diag%axesBi, Time, &
                             'Shear Production of TKE from GOTM','m2 s-3')
    CS%id_GOTM_Eps = &
         register_diag_field('ocean_model','GOTM_EPS', diag%axesBi, Time, &
                             'TKE dissipation from GOTM','m2 s-3')
    if (CS%id_GOTM_N2.gt.0) &
         ALLOC_ (CS%N2(SZIB_(G), SZJB_(G), SZK_(G)+1)) ; CS%N2(:,:,:) = 0.0
    if (CS%id_GOTM_S2.gt.0) &
         ALLOC_ (CS%S2(SZIB_(G), SZJB_(G), SZK_(G)+1)) ; CS%S2(:,:,:) = 0.0
    if (CS%id_GOTM_P.gt.0) &
         ALLOC_ (CS%P(SZIB_(G), SZJB_(G), SZK_(G)+1)) ; CS%P(:,:,:) = 0.0
    if (CS%id_GOTM_B.gt.0) &
         ALLOC_ (CS%B(SZIB_(G), SZJB_(G), SZK_(G)+1)) ; CS%B(:,:,:) = 0.0
  else
    CS%id_TKE = &
         register_diag_field('ocean_model','GOTM_TKE', diag%axesTi, Time, &
                             'Turbulent kinetic energy from GOTM','m2 s-2')
!    CS%id_Dissipation = &
!         register_diag_field('ocean_model','GOTM_Dissipation', diag%axesTi, Time, &
!                             'Turbulent dissipation rate from GOTM','m2 s-3')
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
    CS%id_GOTM_N2 = &
         register_diag_field('ocean_model','GOTM_N2', diag%axesTi, Time, &
                             'Buoyancy frequency squared','s-2')
    CS%id_GOTM_B = &
         register_diag_field('ocean_model','GOTM_B', diag%axesTi, Time, &
                             'Buoyancy Production of TKE from GOTM','m2 s-3')
    CS%id_GOTM_P = &
         register_diag_field('ocean_model','GOTM_P', diag%axesTi, Time, &
                             'Shear Production of TKE from GOTM','m2 s-3')
    CS%id_GOTM_Eps = &
         register_diag_field('ocean_model','GOTM_EPS', diag%axesTi, Time, &
                             'TKE dissipation from GOTM','m2 s-3')
    if (CS%id_GOTM_N2.gt.0) &
         ALLOC_ (CS%N2(SZI_(G), SZJ_(G), SZK_(G)+1)) ; CS%N2(:,:,:) = 0.0
    if (CS%id_GOTM_S2.gt.0) &
         ALLOC_ (CS%S2(SZI_(G), SZJ_(G), SZK_(G)+1)) ; CS%S2(:,:,:) = 0.0
    if (CS%id_GOTM_P.gt.0) &
         ALLOC_ (CS%P(SZI_(G), SZJ_(G), SZK_(G)+1)) ; CS%P(:,:,:) = 0.0
    if (CS%id_GOTM_B.gt.0) &
         ALLOC_ (CS%B(SZI_(G), SZJ_(G), SZK_(G)+1)) ; CS%B(:,:,:) = 0.0
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

!> This subroutine performs the actual stepping of the GOTM turbulence scheme.
!!  NOTE: With GOTM, many flavors of TKE-based vertical mixing parameterizations
!!        can be chosen.  See the gotm documentation at gotm.net for a very
!!        detailed description of the parameterizations.
subroutine GOTM_calculate(CS, G, GV, DT, h, Temp, Salt, u, v, EOS, uStar,&
                          TKE, EPS, L, Kt, Ks, Kv)

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
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(inout)   :: TKE    !< TKE (h points) [m2 s-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(inout)   :: EPS    !< TKE (h points) [m2 s-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(inout)   :: L      !< TKE (h points) [m2 s-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(inout) :: Kt     !< (in)  Vertical diffusivity of heat w/o GOTM [m2 s-1]
                                                                      !< (out) Vertical diffusivity including GOTM [m2 s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(inout) :: Ks     !< (in)  Vertical diffusivity of salt w/o GOTM [m2 s-1]
                                                                      !< (out) Vertical diffusivity including GOTM [m2 s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(inout) :: Kv     !< (in)  Vertical viscosity w/o GOTM [m2 s-1]
                                                                      !< (out) Vertical viscosity including GOTM [m2 s-1]

  ! Local variables
  integer :: i, j, k, klo, kup, kgotm, nk_gotm      ! Loop indices
  real :: ustar_b, z0_b
  real, dimension( G%ke )     :: H_1d ! 1-d level thickness in [m]
  real, dimension( 0:G%ke )   :: g_h  ! 1-d level thickness in [m] for GOTM
  real, dimension( 0:G%ke )   :: g_N2 ! Brunt-Vaisala frequency squared, at interfaces [s-2]
  real, dimension( 0:G%ke )   :: g_S2 ! Shear frequency at interfaces [s-2]

  ! for EOS calculation
  !  *These profiles are size 2x the profile, specifically formulated for computing the N2 profile.
  real, dimension( 2*G%ke )   :: rho_1D  ! Density [kg m-3]
  real, dimension( 2*G%ke )   :: pres_1D ! Pressure [Pa]
  real, dimension( 2*G%ke )   :: Temp_1D ! Temperature [deg C]
  real, dimension( 2*G%ke )   :: Salt_1D ! Salinity [PPT]

  real :: GoRho, pRef, Uabove,Ubelow,Vabove,Vbelow,Habove,Hbelow
  real :: uoabove, uobelow,voabove,vobelow
  real :: delH                 ! Thickness of a layer [h units]
  real :: Dpt, DZ, DU, DV              ! In-situ depth and dz [h units]
  integer :: kk, ksfc, ktmp

  integer :: sCycle

#ifdef __DO_SAFETY_CHECKS__
  if (CS%debug) then
!    call hchksum(h*GV%H_to_m, "KPP in: h",G,haloshift=0)
  endif
#endif

  ! Initializations:
  g_h(:) = GV%H_subroundoff*GV%H_to_m
  g_N2(:) = 0.0
  g_S2(:) = 0.0

  ! some constants
  GoRho = GV%g_Earth / GV%Rho0

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
    do i = G%isc, G%iec; if (G%mask2dT(i,j) > 0.5) then

      pRef = 0.0 ; Dpt = 0.0
      do k=1,G%ke
        H_1d(k) = h(i,j,k)*GV%H_to_m
      enddo

      !compute NN and SS on the GOTM vertical grid
      do K=2,G%ke
        ! Set vertical indices
        kup=K-1
        klo=K
        ! pRef is pressure at interface K
        pRef = pRef + GV%H_to_Pa * H_1d(kup)
        ! pressure, temp, and saln for EOS
        ! kk set for the bounding density of the interfaces, with the first pair for the 1st bounded interface (K=2)
        ! kk+1 = kcenup fields
        ! kk+2 = kcenlo fields
        kk   = 2*(K-2) !(see note above for K=2)
        pres_1D(kk+1) = pRef
        pres_1D(kk+2) = pRef
        Temp_1D(kk+1) = Temp(i,j,kup)
        Temp_1D(kk+2) = Temp(i,j,klo)
        Salt_1D(kk+1) = Salt(i,j,kup)
        Salt_1D(kk+2) = Salt(i,j,klo)
      enddo

      ! compute in-situ density
      call calculate_density(Temp_1D, Salt_1D, pres_1D, rho_1D, 1, 2*(G%ke-1), EOS)

      !Find the cells that contain volume above a certain threshold
      do k=1,G%ke
        if (h_1d(k)>0.1) then
          nk_gotm = k
        endif
      enddo
      ustar_b = sqrt(0.001*((0.5*(u(i,j,kup)+u(i-1,j,kup)))**2+(0.5*(v(i,j,kup)+v(i,j-1,kup)))**2))
      z0_b = 0.001
      !Sum those cells to get a depth
      g_h(:)=GV%H_subroundoff*GV%H_to_m
      Dpt = 0.0
      do k=1,nk_gotm
        kgotm=nk_gotm-k+1
        g_h(kgotm) = H_1d(k)
        Dpt = Dpt + H_1d(k)
      enddo

      
      ! Zero out arrays
      g_N2(:) = 0.0
      g_S2(:) = 0.0
      do K = 2, nk_gotm !Run the loop from the bottom interface of the top cell to the top interface of the bottom cell.
        kup = K-1
        klo = K
        kk = 2*(K-2)
        Kgotm = nk_gotm-K+1 !top is mom's K=1, and gotm's K=nk, botm is mom's K=nk+1 and gotm's K=0

        ! 1. Compute N2
        DZ = 0.5*( h_1d(kup) + h_1d(klo) + GV%H_subroundoff*GV%H_to_m )
        g_N2(kgotm) = GoRho * (rho_1D(kk+2) - rho_1D(kk+1)) / dz

        ! 2. Compute S2

        ! C-grid average to get Uk and Vk on T-points.
        ! k is above, k+1 is below (relative to interface k)
        Uabove  = 0.5*(u(i,j,kup)+u(i-1,j,kup))
        Ubelow  = 0.5*(u(i,j,klo)+u(i-1,j,klo))
        Vabove  = 0.5*(v(i,j,kup)+v(i,j-1,kup))
        Vbelow  = 0.5*(v(i,j,klo)+v(i,j-1,klo))        ! Needs rewritten for vertex points.  Not presently used.
        !UOabove = 0.5*(CS%uo(i,j,k)+CS%uo(i-1,j,k))
        !UObelow = 0.5*(CS%uo(i,j,kp1)+CS%uo(i-1,j,kp1))
        !VOabove = 0.5*(CS%vo(i,j,k)+CS%vo(i,j-1,k))
        !VObelow = 0.5*(CS%vo(i,j,kp1)+CS%vo(i,j-1,kp1))
        Habove  = h_1d(kup)
        Hbelow  = h_1d(klo)

        ! This code is for KE conserving shear following Burchard et al. (2002)
        ! -> need to check that old/present values are correctly identified.
        ! Not obvious that this gives positive-definite shear?
        !DU = 0.5 * ( (Uabove-Ubelow)*(Uabove-UObelow) /&
        !             (0.5*(Habove+Hbelow)*Hbelow)&
        !            +(Uabove-Ubelow)*(UOabove-Ubelow) /&
        !             (0.5*(Habove+Hbelow)*Habove))
        !DV = 0.5 * ( (Vabove-Vbelow)*(Vabove-VObelow) /&
        !             (0.5*(Habove+Hbelow)*Hbelow)&
        !            +(Vabove-Vbelow)*(VOabove-Vbelow) /&
        !             (0.5*(Habove+Hbelow)*Habove))
        ! This code computes shear from U/V at present time-step
        DU = (Uabove-Ubelow)*(Uabove-Ubelow) /&
             (0.5*(Habove+Hbelow+GV%H_subroundoff*GV%H_to_m))**2
        DV = (Vabove-Vbelow)*(Vabove-Vbelow) /&
             (0.5*(Habove+Hbelow+GV%H_subroundoff*GV%H_to_m))**2
        g_S2(kgotm) = DU + DV
      enddo

      ! The next few lines of code should only affect output (boundary
      !  conditions are otherwise set).
      g_S2(0) = g_S2(1)
      g_S2(nk_gotm) = g_S2(nk_gotm-1)
      g_N2(0) = g_N2(1)
      g_N2(nk_gotm) = g_N2(nk_gotm-1)

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
      g_Kdm(:)=0.
      g_Kdt(:)=0.
      g_Kds(:)=0.
      g_EPS(:)=0.
      g_TKE(:)=0.
      g_l(:)=0.
      do K=1,nk_gotm
        Kgotm=nk_gotm-K+1
        !Fill GOTM arrays for time stepping
        g_Kdm(kgotm) = Kv(i,j,k)
        g_Kdt(kgotm) = Kt(i,j,k)
        g_Kds(kgotm) = Ks(i,j,k)
        g_TKE(kgotm) = TKE(i,j,k)
        g_EPS(kgotm) = EPS(i,j,k)
        g_l(kgotm)   = L(i,j,k)
      enddo


      do sCycle=1,CS%NSubCycle
        call do_turbulence(nk_gotm,&! Number of levels
                          dt/real(CS%NSubCycle),&! time step
                          dpt,&! depth
                          ustar(i,j),&! ustar surface
                          ustar_b,&! ustar bottom
                          0.02,&! Z0 surface ! This could be made to compute from Charnock, or run-time
                          z0_b,&! Z0 bottom  ! ditto
                          g_h(0:nk_gotm),&! Thickness (e.g., see calculation above)
                          g_N2(0:nk_gotm),&! Buoyancy Frequency
                          g_S2(0:nk_gotm) &! Shear squared
                          ) ! An additional TKE production source is optional input, neglected here
      enddo
      !--------------------------------------------------------------------------------------------------
      ! Here is what actually happens in this call (first_order):
      ! 1. production
      !    Takes the values of initial diffusivity, intial N2, and initial S2 and computes the production terms (P, B, and Pb)
      !     -> production terms are therefore from the initial state
      ! 2. alpha_mnb
      !    Comes "alpha" values, which include
      !      as = k^2/eps^2 * S^2
      !      an = k^2/eps^2 * N^2
      !      at = k/eps * kb/eps
      !     -> All from the initial state
      ! 3. stabilityfunctions
      !    Various options to compute cmue1 and cmue2
      !     -> All from the initial state
      ! 4. do_tke
      !    Evaluates the time step of tke w/ an implicit diffusion equation
      !   UPDATES: TKE
      !   USES: production terms, epsilon, diffusivity
      !    -> TKE projected to future state from production terms from initial state
      ! 5. do_lengthscale
      !    Evaluates the time step of epsilon w/ an implicit diffusion equation
      !   UPDATES: epsilon, length scale
      !   USES: production terms, TKE, diffusivity
      !    -> length scale projected to future state from projected TKE and initial state diffusivity and production terms
      ! 6. kolpran
      !    Updates the diffusivities
      !    -> diffusivity projected to future state from projected TKE and length-scale
      !   UPDATES: diffusivity
      !   USES: TKE, length scale, cmue1 and cmue2
      !----------------------------------------------------------------------------------------------------

      ! Writing everything w/ the value from the interface above the bottom.
      ! We will fill the interfaces above that point below.
      Kv(I,J,:) = g_Kdm(1)
      Ks(I,J,:) = g_Kds(1)
      Kt(I,J,:) = g_Kdt(1)
      TKE(I,J,:) = g_TKE(1)
      EPS(I,J,:) = g_EPS(1)
      L(I,J,:) = g_l(1)
      CS%N2(i,j,:) = 0.0
      CS%S2(i,j,:) = 0.0
      CS%B(i,j,:) = 0.0
      CS%P(i,j,:) = 0.0
      do K=1,nk_gotm
        Kgotm=nk_gotm-K+1
        !Fill MOM arrays from GOTM
        CS%Kv(i,j,k) = g_Kdm(kgotm)
        CS%Kt(i,j,k) = g_Kdt(kgotm)
        CS%Ks(i,j,k) = g_Kds(kgotm)
        Kv(i,j,k) = g_Kdm(kgotm)
        Kt(i,j,k) = g_Kdt(kgotm)
        Ks(i,j,k) = g_Kds(kgotm)
        TKE(i,j,k) = g_TKE(kgotm)
        EPS(i,j,k) = g_EPS(kgotm)
        L(i,j,k) = g_l(kgotm)
        ! Writing for Output the shear/buoyancy frequencies
        if (CS%id_GOTM_S2.gt.0) &
             CS%S2(i,j,k) = g_S2(kgotm)
        if (CS%id_GOTM_N2.gt.0) &
             CS%N2(i,j,k) = g_N2(kgotm)
        if (CS%id_GOTM_B.gt.0) &
             CS%B(i,j,k) = g_B(kgotm)
        if (CS%id_GOTM_P.gt.0) &
             CS%P(i,j,k) = g_P(kgotm)
      enddo
    else
      CS%Kv(i,j,:) = CS%minK
      CS%Kt(i,j,:) = CS%minK
      CS%Ks(i,j,:) = CS%minK
      CS%N2(i,j,:) = 0.0
      CS%S2(i,j,:) = 0.0
    endif;enddo ! i
  enddo ! j

  ! Set old values (needed for KE conserving S2 calculation)
  CS%UO(:,:,:) = U(:,:,:)
  CS%VO(:,:,:) = V(:,:,:)

!/ Diagnostics
  if (CS%id_TKE.gt.0)         call post_data(CS%id_TKE        , TKE, CS%diag)
!  if (CS%id_Dissipation.gt.0) call post_data(CS%id_Dissipation, EPS, CS%diag)
  if (CS%id_LengthScale.gt.0) call post_data(CS%id_LengthScale, L  , CS%diag)
  if (CS%id_GOTM_KM.gt.0)     call post_data(CS%id_GOTM_KM    , Kv , CS%diag)
  if (CS%id_GOTM_KH.gt.0)     call post_data(CS%id_GOTM_KH    , Kt , CS%diag)
  if (CS%id_GOTM_KS.gt.0)     call post_data(CS%id_GOTM_KS    , Ks , CS%diag)
  if (CS%id_GOTM_N2.gt.0)     call post_data(CS%id_GOTM_N2    , CS%N2 , CS%diag)
  if (CS%id_GOTM_S2.gt.0)     call post_data(CS%id_GOTM_S2    , CS%S2 , CS%diag)
  if (CS%id_GOTM_B.gt.0)      call post_data(CS%id_GOTM_B     , CS%B  , CS%diag)
  if (CS%id_GOTM_P.gt.0)      call post_data(CS%id_GOTM_P     , CS%P  , CS%diag)
  if (CS%id_GOTM_Eps.gt.0)    call post_data(CS%id_GOTM_Eps   , Eps, CS%diag)
!\

#ifdef __DO_SAFETY_CHECKS__
  if (CS%debug) then
!    call hchksum(Ks, "KPP out: Ks",G,haloshift=0)
  endif
#endif

end subroutine GOTM_calculate

!> This subroutine performs the actual stepping of the GOTM turbulence scheme.
!!  NOTE: With GOTM, many flavors of TKE-based vertical mixing parameterizations
!!        can be chosen.  See the gotm documentation at gotm.net for a very
!!        detailed description of the parameterizations.
subroutine GOTM_calculate_vertex(CS, G, GV, DT, h, Temp, Salt, u, v, EOS, uStar,&
                                 TKE, EPS, L, Kt, Ks, Kv_Bu)

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
  real, dimension(SZI_(G),SZJ_(G)),         intent(inout)      :: uStar  !< Surface friction velocity [m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(inout) :: TKE    !< TKE (vertex points) [m2 s-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(inout) :: EPS    !< dissipation (vertex points) [m2 s-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(inout) :: L      !< mixing length (vertex points) [m2 s-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(inout) :: Kt     !< (in)  Vertical diffusivity of heat w/o GOTM [m2 s-1]
                                                                      !< (out) Vertical diffusivity including GOTM [m2 s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(inout) :: Ks     !< (in)  Vertical diffusivity of salt w/o GOTM [m2 s-1]
                                                                      !< (out) Vertical diffusivity including GOTM [m2 s-1]
  real, dimension(SZIB_(G),SZJB_(G),SZK_(G)+1), intent(inout) :: Kv_Bu !< (in)  Vertical viscosity w/o GOTM [m2 s-1]
                                                                      !< (out) Vertical viscosity including GOTM [m2 s-1]

  real, dimension(SZIB_(G),SZK_(GV)) :: &
    h_2d, &             ! A 2-D version of h, but converted to [Z ~> m].
    u_2d, v_2d, &       ! 2-D versions of u_in and v_in, converted to [L T-1 ~> m s-1].
    T_2d, S_2d          ! 2-D versions of T [degC], S [ppt], and rho [R ~> kg m-3].
  real, dimension(SZIB_(G),SZK_(GV)+1) :: &
    Ks_2d, Kt_2d, TKE_2d, EPS_2d, L_2d        ! A 2-D version of Kx, but converted to [].
  real, dimension(SZIB_(G)) :: ust_1d ! A 1-D version of ustar

  ! Local variables
  integer :: i, j, k, kup, klo, kgotm, nk_gotm      ! Loop indices
  real, dimension( G%ke )     :: H_1d ! 1-d level thickness in [m]
  real, dimension( 0:G%ke )   :: g_h  ! 1-d level thickness in [m] for GOTM
  real, dimension( 0:G%ke )   :: g_N2 ! Brunt-Vaisala frequency squared, at interfaces [s-2]
  real, dimension( 0:G%ke )   :: g_S2 ! Shear frequency at interfaces [s-2]

  ! for EOS calculation
  !  *These profiles are size 2x the number of bounded interfaces, neglecting the top and bottom
  !   specifically formulated for computing the N2 profile.
  real, dimension( 2*(G%ke-1) )   :: rho_1D  ! Density [kg m-3]
  real, dimension( 2*(G%ke-1) )   :: pres_1D ! Pressure [Pa]
  real, dimension( 2*(G%ke-1) )   :: Temp_1D ! Temperature [deg C]
  real, dimension( 2*(G%ke-1) )   :: Salt_1D ! Salinity [PPT]

  real :: GoRho, pRef, Uabove,Ubelow,Vabove,Vbelow,Habove,Hbelow
  real :: uoabove, uobelow,voabove,vobelow
  real :: ustar_b, z0_b
  real :: delH                 ! Thickness of a layer [h units]
  real :: Dpt, DZ, DU, DV              ! In-situ depth and dz [h units]
  integer :: kk, ksfc, ktmp

  integer :: sCycle

  real :: i_hwt, i_intwt, mask_weight
  integer :: IsB, IeB, JsB, JeB, nz, nzc

  isB = G%isc-1 ; ieB = G%iecB ; jsB = G%jsc-1 ; jeB = G%jecB ; nz = GV%ke

#ifdef __DO_SAFETY_CHECKS__
  if (CS%debug) then
!    call hchksum(h*GV%H_to_m, "KPP in: h",G,haloshift=0)
  endif
#endif

  ! Initializations:
  g_h(:) = GV%H_subroundoff*GV%H_to_m
  g_N2(:) = 0.0
  g_S2(:) = 0.0

  ! The input diffusivities need to be updated on Halos somewhere before computing on the vertices.
  call pass_var(Ks,G%Domain)
  call pass_var(Kt,G%Domain)
  call pass_var(TKE,G%Domain)
  call pass_var(EPS,G%Domain)
  call pass_var(L,G%Domain)
  !Might this needed?  I don't think we want this to be intent in/out.
  call pass_var(ustar,G%Domain)
  ! some constants
  GoRho = GV%g_Earth / GV%Rho0

  do J=JsB,JeB
    do k=1,nz+1 ; do I=IsB,IeB
      mask_weight = G%mask2dT(i,j) + G%mask2dT(i+1,j+1) + G%mask2dT(i+1,j) + G%mask2dT(i,j+1)
      if (mask_weight>0.) then
        ! Interpolate the various interface quantities that are on cell centers to the corners, using masks.
        i_intwt = 1./mask_weight
        if (k==1) then
          ust_1d(I) = (G%mask2dT(i,j)     * ustar(i,j)     + &
                      G%mask2dT(i+1,j+1) * ustar(i+1,j+1) + &
                      G%mask2dT(i+1,j)   * ustar(i+1,j)   + &
                      G%mask2dT(i,j+1)   * ustar(i,j+1) ) * i_intwt
        endif
        Ks_2d(I,k) = (G%mask2dT(i,j)     * Ks(i,j,k)     + &
                      G%mask2dT(i+1,j+1) * Ks(i+1,j+1,k) + &
                      G%mask2dT(i+1,j)   * Ks(i+1,j,k)   + &
                      G%mask2dT(i,j+1)   * Ks(i,j+1,k) ) * i_intwt
        Kt_2d(I,k) = (G%mask2dT(i,j)     * Kt(i,j,k)     + &
                      G%mask2dT(i+1,j+1) * Kt(i+1,j+1,k) + &
                      G%mask2dT(i+1,j)   * Kt(i+1,j,k)   + &
                      G%mask2dT(i,j+1)   * Kt(i,j+1,k) ) * i_intwt
        EPS_2d(I,k) = (G%mask2dT(i,j)    * EPS(i,j,k)     + &
                      G%mask2dT(i+1,j+1) * EPS(i+1,j+1,k) + &
                      G%mask2dT(i+1,j)   * EPS(i+1,j,k)   + &
                      G%mask2dT(i,j+1)   * EPS(i,j+1,k) ) * i_intwt
        TKE_2d(I,k) = (G%mask2dT(i,j)    * TKE(i,j,k)     + &
                      G%mask2dT(i+1,j+1) * TKE(i+1,j+1,k) + &
                      G%mask2dT(i+1,j)   * TKE(i+1,j,k)   + &
                      G%mask2dT(i,j+1)   * TKE(i,j+1,k) ) * i_intwt
        L_2d(I,k) =  (G%mask2dT(i,j)     * L(i,j,k)     + &
                      G%mask2dT(i+1,j+1) * L(i+1,j+1,k) + &
                      G%mask2dT(i+1,j)   * L(i+1,j,k)   + &
                      G%mask2dT(i,j+1)   * L(i,j+1,k) ) * i_intwt
      else !if no surrounding cell centers are open the corner should also be land.
        Ks_2d(I,k) = CS%minK
        Kt_2d(I,k) = CS%minK
      endif
    enddo; enddo
    ! Interpolate the various layer quantities that are on c-grid to the corners, using masks.
    do k=1,nz ; do I=IsB,IeB
      u_2d(I,k) = (u(I,j,k)   * (G%mask2dCu(I,j)   * (h(i,j,k)   + h(i+1,j,k))) + &
                   u(I,j+1,k) * (G%mask2dCu(I,j+1) * (h(i,j+1,k) + h(i+1,j+1,k))) ) / &
                  ((G%mask2dCu(I,j)   * (h(i,j,k)   + h(i+1,j,k)) + &
                    G%mask2dCu(I,j+1) * (h(i,j+1,k) + h(i+1,j+1,k))) + GV%H_subroundoff)
      v_2d(I,k) = (v(i,J,k)   * (G%mask2dCv(i,J)   * (h(i,j,k)   + h(i,j+1,k))) + &
                   v(i+1,J,k) * (G%mask2dCv(i+1,J) * (h(i+1,j,k) + h(i+1,j+1,k))) ) / &
                  ((G%mask2dCv(i,J)   * (h(i,j,k)   + h(i,j+1,k)) + &
                    G%mask2dCv(i+1,J) * (h(i+1,j,k) + h(i+1,j+1,k))) + GV%H_subroundoff)
      i_hwt = 1.0 / (((G%mask2dT(i,j) * h(i,j,k) + G%mask2dT(i+1,j+1) * h(i+1,j+1,k)) + &
                      (G%mask2dT(i+1,j) * h(i+1,j,k) + G%mask2dT(i,j+1) * h(i,j+1,k))) + &
                     GV%H_subroundoff)
      T_2d(I,k) = ( ((G%mask2dT(i,j) * h(i,j,k)) * Temp(i,j,k) + &
                     (G%mask2dT(i+1,j+1) * h(i+1,j+1,k)) * Temp(i+1,j+1,k)) + &
                    ((G%mask2dT(i+1,j) * h(i+1,j,k)) * Temp(i+1,j,k) + &
                    (G%mask2dT(i,j+1) * h(i,j+1,k)) * Temp(i,j+1,k)) ) * i_hwt
      S_2d(I,k) = ( ((G%mask2dT(i,j) * h(i,j,k)) * Salt(i,j,k) + &
                     (G%mask2dT(i+1,j+1) * h(i+1,j+1,k)) * Salt(i+1,j+1,k)) + &
                    ((G%mask2dT(i+1,j) * h(i+1,j,k)) * Salt(i+1,j,k) + &
                     (G%mask2dT(i,j+1) * h(i,j+1,k)) * Salt(i,j+1,k)) ) * i_hwt
      h_2d(I,k) = GV%H_to_Z * ((G%mask2dT(i,j) * h(i,j,k) + G%mask2dT(i+1,j+1) * h(i+1,j+1,k)) + &
                               (G%mask2dT(i+1,j) * h(i+1,j,k) + G%mask2dT(i,j+1) * h(i,j+1,k)) ) / &
                              ((G%mask2dT(i,j) + G%mask2dT(i+1,j+1)) + &
                              (G%mask2dT(i+1,j) + G%mask2dT(i,j+1)) + 1.0e-36 )
    enddo ; enddo;
    !---------------------------------------
    ! Work on each column.
    !---------------------------------------
    do I=IsB,IeB ; if ((G%mask2dCu(I,j) + G%mask2dCu(I,j+1)) + &
                       (G%mask2dCv(i,J) + G%mask2dCv(i+1,J)) > 0.0) then
      pRef = 0.0 ; Dpt = 0.0
      do k=1,G%ke
        H_1d(k) = h_2d(I,k)*GV%H_to_m
      enddo

      !compute NN and SS on the GOTM vertical grid
      do K=2,G%ke
        ! Set vertical indices
        kup=K-1
        klo=K
        ! pRef is pressure at interface K
        pRef = pRef + GV%H_to_Pa * H_1d(kup)
        ! pressure, temp, and saln for EOS
        ! kk set for the bounding density of the interfaces, with the first pair for the 1st bounded interface (K=2)
        ! kk+1 = kcenup fields
        ! kk+2 = kcenlo fields
        kk   = 2*(K-2) !(see note above for K=2)
        pres_1D(kk+1) = pRef
        pres_1D(kk+2) = pRef
        Temp_1D(kk+1) = T_2d(I,kup)
        Temp_1D(kk+2) = T_2d(I,klo)
        Salt_1D(kk+1) = S_2d(I,kup)
        Salt_1D(kk+2) = S_2d(I,klo)
      enddo

      ! compute in-situ density
      call calculate_density(Temp_1D, Salt_1D, pres_1D, rho_1D, 1, 2*(G%ke-1), EOS)

      !Find the cells that contain volume above a certain threshold
      do k=1,G%ke
        if (h_1d(k)>0.1) then
          nk_gotm = k
        endif
      enddo
      ustar_b = sqrt(0.001*(u_2d(I,nk_gotm)**2+v_2d(I,nk_gotm)**2))
      z0_b = 0.001
      !Sum those cells to get a depth
      g_h(:)=GV%H_subroundoff*GV%H_to_m
      Dpt = 0.0
      do k=1,nk_gotm
        kgotm=nk_gotm-k+1
        g_h(kgotm) = H_1d(k)
        Dpt = Dpt + H_1d(k)
      enddo

      ! Zero out arrays
      g_N2(:) = 0.0
      g_S2(:) = 0.0
      do K = 2, nk_gotm !Run the loop from the bottom interface of the top cell to the top interface of the bottom cell.
        kup = K-1
        klo = K
        kk = 2*(K-2)
        Kgotm = nk_gotm-K+1 !top is mom's K=1, and gotm's K=nk, botm is mom's K=nk+1 and gotm's K=0

        ! 1. Compute N2
        DZ = 0.5*( h_1d(kup) + h_1d(klo) + GV%H_subroundoff*GV%H_to_m )
        g_N2(kgotm) = GoRho * (rho_1D(kk+2) - rho_1D(kk+1)) / dz

        ! 2. Compute S2

        ! C-grid average to get Uk and Vk on T-points.
        ! k is above, k+1 is below (relative to interface k)
        Uabove  = u_2d(I,kup)
        Ubelow  = u_2d(I,klo)
        Vabove  = v_2d(I,kup)
        Vbelow  = v_2d(I,klo)
        ! Needs rewritten for vertex points.  Not presently used.
        !UOabove = 0.5*(CS%uo(i,j,k)+CS%uo(i-1,j,k))
        !UObelow = 0.5*(CS%uo(i,j,kp1)+CS%uo(i-1,j,kp1))
        !VOabove = 0.5*(CS%vo(i,j,k)+CS%vo(i,j-1,k))
        !VObelow = 0.5*(CS%vo(i,j,kp1)+CS%vo(i,j-1,kp1))
        Habove  = h_1d(kup)
        Hbelow  = h_1d(klo)

        ! This code is for KE conserving shear following Burchard et al. (2002)
        ! -> need to check that old/present values are correctly identified.
        ! Not obvious that this gives positive-definite shear?
        !DU = 0.5 * ( (Uabove-Ubelow)*(Uabove-UObelow) /&
        !             (0.5*(Habove+Hbelow)*Hbelow)&
        !            +(Uabove-Ubelow)*(UOabove-Ubelow) /&
        !             (0.5*(Habove+Hbelow)*Habove))
        !DV = 0.5 * ( (Vabove-Vbelow)*(Vabove-VObelow) /&
        !             (0.5*(Habove+Hbelow)*Hbelow)&
        !            +(Vabove-Vbelow)*(VOabove-Vbelow) /&
        !             (0.5*(Habove+Hbelow)*Habove))
        ! This code computes shear from U/V at present time-step
        DU = (Uabove-Ubelow)*(Uabove-Ubelow) /&
             (0.5*(Habove+Hbelow+GV%H_subroundoff*GV%H_to_m))**2
        DV = (Vabove-Vbelow)*(Vabove-Vbelow) /&
             (0.5*(Habove+Hbelow+GV%H_subroundoff*GV%H_to_m))**2
        g_S2(kgotm) = DU + DV
      enddo

      ! The next few lines of code should only affect output (boundary
      !  conditions are otherwise set).
      g_S2(0) = g_S2(1)
      g_S2(nk_gotm) = g_S2(nk_gotm-1)
      g_N2(0) = g_N2(1)
      g_N2(nk_gotm) = g_N2(nk_gotm-1)

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
      do K=1,nk_gotm
        Kgotm=nk_gotm-K+1
        !Fill GOTM arrays for time stepping
        g_Kdm(Kgotm) = Kv_Bu(i,j,K)
        g_Kdt(Kgotm) = Kt_2d(i,K)
        g_Kds(Kgotm) = Ks_2d(i,K)
        g_TKE(Kgotm) = TKE_2d(i,K)
        g_EPS(Kgotm) = EPS_2d(i,K)
        g_l(Kgotm)   = L_2d(i,K)
      enddo

      do sCycle=1,CS%NSubCycle
        call do_turbulence(nk_gotm,&! Number of levels
                          dt/real(CS%NSubCycle),&! time step
                          dpt,&! depth
                          ust_1d(i),&! ustar surface
                          ustar_b,&! ustar bottom
                          0.02,&! Z0 surface ! This could be made to compute from Charnock, or run-time
                          z0_b,&! Z0 bottom  ! ditto
                          g_h(0:nk_gotm),&! Thickness (e.g., see calculation above)
                          g_N2(0:nk_gotm),&! Buoyancy Frequency
                          g_S2(0:nk_gotm) &! Shear squared
                          ) ! An additional TKE production source is optional input, neglected here
      enddo

      ! Writing everything w/ the value from the interface above the bottom.
      ! We will fill the interfaces above that point below.
      CS%Kv(I,J,:) = g_Kdm(1)
      CS%Ks(I,J,:) = g_Kds(1)
      CS%Kt(I,J,:) = g_Kdt(1)
      CS%TKE(I,J,:) = g_TKE(1)
      CS%EPS(I,J,:) = g_EPS(1)
      CS%L(I,J,:) = g_l(1)
      CS%N2(i,j,:) = 0.0
      CS%S2(i,j,:) = 0.0
      CS%B(i,j,:) = 0.0
      CS%P(i,j,:) = 0.0
      do K=1,nk_gotm
        Kgotm=nk_gotm-K+1
        !Fill MOM arrays from GOTM
        CS%Kv(I,J,k) = g_Kdm(kgotm) !for output on corners
        Kv_Bu(I,J,k) = g_Kdm(kgotm) !for passing back to model
        CS%Kt(I,J,k) = g_Kdt(kgotm) !for outputs on corners and interpolated horizontally to pass back
        CS%Ks(I,J,k) = g_Kds(kgotm) !for outputs on corners and interpolated horizontally to pass back
        CS%TKE(I,J,k) = g_TKE(kgotm)
        CS%EPS(I,J,k) = g_EPS(kgotm)
        CS%L(I,J,k) = g_l(kgotm)
        !if (CS%TKE(I,J,k)>0.01) then
        !   print*,k,G%geolatt(I,J),G%geolont(I,J)
        !   print*,k,ust_1d(i),g_Kdm
        !   print*,k,dpt,g_h
        !   print*,k,'tke',g_TKE
        !   print*,k,'diss',g_EPS
        !   print*,k,'n2',g_N2
        !   print*,k,'s2',g_S2
        !   stop
        !endif
        ! Writing for Output the shear/buoyancy frequencies
        if (CS%id_GOTM_S2.gt.0) &
             CS%S2(i,j,k) = g_S2(kgotm)
        if (CS%id_GOTM_N2.gt.0) &
             CS%N2(i,j,k) = g_N2(kgotm)
        if (CS%id_GOTM_B.gt.0) &
             CS%B(i,j,k) = g_B(kgotm)
        if (CS%id_GOTM_P.gt.0) &
             CS%P(i,j,k) = g_P(kgotm)
      enddo
    else
     CS%Kv(i,j,:) = CS%minK
     CS%Kt(i,j,:) = CS%minK
     CS%Ks(i,j,:) = CS%minK
     CS%N2(i,j,:) = 0.0
     CS%S2(i,j,:) = 0.0
     CS%TKE(i,j,:) = 1.e-10 !When smoothing, creates a floor
     CS%EPS(i,j,:) = 1.e-12 !When smoothing, creates a floor
     CS%L(I,J,:) = 1.e-8 !When smoothing, creates a floor
   endif;enddo ! i
  enddo ! j

  ! Interpolate Kt/Ks from Q to H
  do j=G%jsc,G%jec
    do i=G%isc,G%iec
      do k=2,G%ke
        if ((G%mask2dBu(I,J) + G%mask2dBu(I-1,J) + &
              G%mask2dBu(I,J-1) + G%mask2dBu(I-1,J-1)) > 1.) then
          Kt(i,j,K) = (CS%Kt(I,J,K)*G%mask2dBu(I,J) + &
                       CS%Kt(I-1,J,K)*G%mask2dBu(I-1,J) + &
                       CS%Kt(I,J-1,K)*G%mask2dBu(I,J-1) + &
                       CS%Kt(I-1,J-1,K)*G%mask2dBu(I-1,J-1)) / &
                      (G%mask2dBu(I,J) + G%mask2dBu(I-1,J) + &
                       G%mask2dBu(I,J-1) + G%mask2dBu(I-1,J-1))
          Ks(i,j,K) = (CS%Ks(I,J,K)*G%mask2dBu(I,J) + &
                       CS%Ks(I-1,J,K)*G%mask2dBu(I-1,J) + &
                       CS%Ks(I,J-1,K)*G%mask2dBu(I,J-1) + &
                       CS%Ks(I-1,J-1,K)*G%mask2dBu(I-1,J-1)) / &
                      (G%mask2dBu(I,J) + G%mask2dBu(I-1,J) + &
                      G%mask2dBu(I,J-1) + G%mask2dBu(I-1,J-1))
          TKE(i,j,K) = (CS%TKE(I,J,K)*G%mask2dBu(I,J) + &
                       CS%TKE(I-1,J,K)*G%mask2dBu(I-1,J) + &
                       CS%TKE(I,J-1,K)*G%mask2dBu(I,J-1) + &
                       CS%TKE(I-1,J-1,K)*G%mask2dBu(I-1,J-1)) / &
                      (G%mask2dBu(I,J) + G%mask2dBu(I-1,J) + &
                      G%mask2dBu(I,J-1) + G%mask2dBu(I-1,J-1))
          EPS(i,j,K) = (CS%EPS(I,J,K)*G%mask2dBu(I,J) + &
                       CS%EPS(I-1,J,K)*G%mask2dBu(I-1,J) + &
                       CS%EPS(I,J-1,K)*G%mask2dBu(I,J-1) + &
                       CS%EPS(I-1,J-1,K)*G%mask2dBu(I-1,J-1)) / &
                      (G%mask2dBu(I,J) + G%mask2dBu(I-1,J) + &
                      G%mask2dBu(I,J-1) + G%mask2dBu(I-1,J-1))
          L(i,j,K) = (CS%L(I,J,K)*G%mask2dBu(I,J) + &
                       CS%L(I-1,J,K)*G%mask2dBu(I-1,J) + &
                       CS%L(I,J-1,K)*G%mask2dBu(I,J-1) + &
                       CS%L(I-1,J-1,K)*G%mask2dBu(I-1,J-1)) / &
                      (G%mask2dBu(I,J) + G%mask2dBu(I-1,J) + &
                      G%mask2dBu(I,J-1) + G%mask2dBu(I-1,J-1))
        else
          Kt(i,j,K) = 0.0
          Ks(i,j,K) = 0.0
          TKE(i,j,K) = 0.0
          EPS(i,j,K) = 0.0
          L(i,j,K) = 0.0
        endif
      enddo
    enddo
  enddo
  ! Set old values (needed for KE conserving S2 calculation)
  CS%UO(:,:,:) = U(:,:,:)
  CS%VO(:,:,:) = V(:,:,:)

  !Smooth TKE/Eps?
  !call pass_var(TKE_gotm,G%Domain)
  !call pass_var(EPS_gotm,G%Domain)
  !do J=JsB,JeB; do I=IsB,IeB
  !  mask_weight =    G%mask2dBu(I-1,J-1) + 2.*G%mask2dBu(I-1,J) +    G%mask2dBu(I-1,J+1) &
!                   +2.*G%mask2dBu(I,J-1)   + 4.*G%mask2dBu(I,J)   + 2.*G%mask2dBu(I,J+1) &
  !                 +2.*G%mask2dBu(I,J-1)   + 4.                   + 2.*G%mask2dBu(I,J+1) &
  !                 +   G%mask2dBu(I+1,J-1) + 2.*G%mask2dBu(I+1,J) +    G%mask2dBu(I+1,J+1)
  !  do k=1,G%ke+1
  !    TKE(I,J,k) = (G%mask2dBu(I-1,J-1)*TKE_gotm(I-1,J-1,k) +&
  !                  2.*G%mask2dBu(I-1,J)*TKE_gotm(I-1,J,k) +&
  !                  G%mask2dBu(I-1,J+1)*TKE_gotm(I-1,J+1,k) +&
  !                  2.*G%mask2dBu(I,J-1)*TKE_gotm(I,J-1,k) +&
!                    4.*G%mask2dBu(I,J)*TKE_gotm(I,J,k) +&
  !                  4.*TKE_gotm(I,J,k) +&
  !                  2.*G%mask2dBu(I,J+1)*TKE_gotm(I,J+1,k) +&
  !                  G%mask2dBu(I+1,J-1)*TKE_gotm(I+1,J-1,k) +&
  !                  2.*G%mask2dBu(I+1,J)*TKE_gotm(I+1,J,k) +&
  !                  G%mask2dBu(I+1,J+1)*TKE_gotm(I+1,J+1,k))/mask_weight
  !     EPS(I,J,k) = (G%mask2dBu(I-1,J-1)*EPS_gotm(I-1,J-1,k) +&
  !                   2.*G%mask2dBu(I-1,J)*EPS_gotm(I-1,J,k) +&
  !                   G%mask2dBu(I-1,J+1)*EPS_gotm(I-1,J+1,k) +&
  !                   2.*G%mask2dBu(I,J-1)*EPS_gotm(I,J-1,k) +&
!                     4.*G%mask2dBu(I,J)*EPS_gotm(I,J,k) +&
  !                   4.*EPS_gotm(I,J,k) +&
  !                   2.*G%mask2dBu(I,J+1)*EPS_gotm(I,J+1,k) +&
  !                   G%mask2dBu(I+1,J-1)*EPS_gotm(I+1,J-1,k) +&
  !                   2.*G%mask2dBu(I+1,J)*EPS_gotm(I+1,J,k) +&
  !                   G%mask2dBu(I+1,J+1)*EPS_gotm(I+1,J+1,k))/mask_weight
  !enddo ; enddo ; enddo

  ! do J=JsB,JeB; do I=IsB,IeB ;do k=2,G%ke
  !  if (EPS(I,J,k)==0. .and. EPS_gotm(I,J,k).ne.0.) then
  !      print*,I,J,k
  !      print*,i+j+k,EPS(I-1,J-1,k),EPS_gotm(I-1,J-1,k),G%mask2dBu(I-1,J-1)
  !      print*,i+j+k,EPS(I,J-1,k),  EPS_gotm(I,J-1,k),G%mask2dBu(I,J-1)
  !      print*,i+j+k,EPS(I+1,J-1,k),EPS_gotm(I+1,J-1,k),G%mask2dBu(I+1,J-1)
  !      print*,i+j+k,EPS(I-1,J,k),EPS_gotm(I-1,J,k),G%mask2dBu(I-1,J)
  !      print*,i+j+k,EPS(I,J,k),  EPS_gotm(I,J,k),G%mask2dBu(I,J)
  !      print*,i+j+k,EPS(I+1,J,k),EPS_gotm(I+1,J,k),G%mask2dBu(I+1,J)
  !      print*,i+j+k,EPS(I-1,J+1,k),EPS_gotm(I-1,J+1,k),G%mask2dBu(I-1,J+1)
  !      print*,i+j+k,EPS(I,J+1,k),  EPS_gotm(I,J+1,k),G%mask2dBu(I,J+1)
  !      print*,i+j+k,EPS(I+1,J+1,k),EPS_gotm(I+1,J+1,k),G%mask2dBu(I+1,J+1)
  !      stop
  !   endif
  !   enddo ; enddo ; enddo

!/ Diagnostics
  if (CS%id_TKE.gt.0)         call post_data(CS%id_TKE        , CS%TKE, CS%diag)
  !if (CS%id_Dissipation.gt.0) call post_data(CS%id_Dissipation, EPS, CS%diag)
  if (CS%id_GOTM_Eps.gt.0)    call post_data(CS%id_GOTM_Eps   , CS%EPS, CS%diag)
  if (CS%id_LengthScale.gt.0) call post_data(CS%id_LengthScale, CS%L  , CS%diag)
  if (CS%id_GOTM_KM.gt.0)     call post_data(CS%id_GOTM_KM    , CS%Kv , CS%diag)
  if (CS%id_GOTM_KH.gt.0)     call post_data(CS%id_GOTM_KH    , CS%Kt , CS%diag)
  if (CS%id_GOTM_KS.gt.0)     call post_data(CS%id_GOTM_KS    , CS%Ks , CS%diag)
  if (CS%id_GOTM_N2.gt.0)     call post_data(CS%id_GOTM_N2    , CS%N2 , CS%diag)
  if (CS%id_GOTM_S2.gt.0)     call post_data(CS%id_GOTM_S2    , CS%S2 , CS%diag)
  if (CS%id_GOTM_B.gt.0)      call post_data(CS%id_GOTM_B     , CS%B  , CS%diag)
  if (CS%id_GOTM_P.gt.0)      call post_data(CS%id_GOTM_P     , CS%P  , CS%diag)
!\

#ifdef __DO_SAFETY_CHECKS__
  if (CS%debug) then
!    call hchksum(Ks, "KPP out: Ks",G,haloshift=0)
  endif
#endif

end subroutine GOTM_calculate_vertex

logical function GOTM_at_vertex(param_file)
  type(param_file_type), intent(in) :: param_file !< A structure to parse for run-time parameters

  ! Local variables
  character(len=40)  :: mdl = "MOM_GOTM"  ! This module's name.
  logical :: do_GOTM
  ! This function returns true only if the parameters "USE_GOTM" and "GOTM_VERTEX_SHEAR" are both true.

  GOTM_at_vertex = .false.

  call get_param(param_file, mdl, "USE_GOTM", do_GOTM, &
                 default=.false., do_not_log=.true.)
  if (do_GOTM) &
    call get_param(param_file, mdl, "GOTM_VERTEX_SHEAR", GOTM_at_vertex, &
                 "If true, do the calculations of the GOTM mixing "//&
                 "at the cell vertices (i.e., the vorticity points).", &
                 default=.false., do_not_log=.false.)

end function GOTM_at_vertex

end module MOM_GOTM
