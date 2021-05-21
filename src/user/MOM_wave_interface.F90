!> Interface for surface waves
module MOM_wave_interface

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_alloc
use MOM_diag_mediator, only : diag_ctrl
use MOM_domains,       only : pass_var, pass_vector, AGRID
use MOM_domains,       only : To_South, To_West, To_All
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_forcing_type, only : mech_forcing
use MOM_grid,          only : ocean_grid_type
use MOM_safe_alloc,    only : safe_alloc_ptr
use MOM_time_manager,  only : time_type, operator(+), operator(/)
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs, surface
use MOM_verticalgrid,  only : verticalGrid_type
use data_override_mod, only : data_override_init, data_override

implicit none ; private

#include <MOM_memory.h>

public MOM_wave_interface_init, MOM_wave_interface_init_lite
public Update_Surface_Waves, Update_Stokes_Drift
public get_Langmuir_Number
public Stokes_PGF
public StokesMixing, CoriolisStokes
public Waves_end

integer, private, parameter :: NULL_WAVEMETHOD = -99
integer, private, parameter :: TESTPROF = 0
integer, private, parameter :: SURFBANDS = 1
integer, private, parameter :: DHH85 = 2
integer, private, parameter :: LF17 = 3
integer, private, parameter :: WIND = 4
integer, private, parameter :: DATAOVR = 1
integer, private, parameter :: COUPLER = 2
integer, private, parameter :: INPUT = 3

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Control structure for all surface wave related parameters
type, public :: wave_parameters_CS ; private

  logical, public :: UseWaves !< Master Flag for module

  ! Primary "Process" Flags
  logical, public :: LagrangianMixing !< Developmental:
                                      !! True if Stokes drift is present and mixing
                                      !! should be applied to Lagrangian current
                                      !! (mean current + Stokes drift).
  logical, public :: StokesMixing !< Developmental:
                                  !! True if vertical mixing of momentum
                                      !! should be applied directly to Stokes current
                                      !! (with separate mixing parameter for Eulerian
                                      !! mixing contribution).
  logical, public :: Stokes_VF     !< Developmental:
                                   !! True if Stokes vortex force is used
  logical, public :: Stokes_PGF    !< Developmental:
                                   !! True if Stokes shear pressure Gradient force is used
  logical, public :: Passive_Stokes_PGF !< Keeps Stokes_PGF on, but doesn't affect dynamics
  logical, public :: Stokes_DDT    !< Developmental:
                                   !! True if Stokes d/dt is used

  
  ! Primary source flags
  integer :: DataSource !< Integer that specifies where the Model Looks for Data
                      !! Valid choices are:
                      !! 1 - FMS DataOverride Routine
                      !! 2 - Reserved For Coupler
                      !! 3 - User input (fixed values, useful for 1d testing)

  ! Other Integer/Logical Flags
  integer, public :: StkLevelMode=1   !< Flags if Stokes drift is defined at mid-points or layer
                                      !! averaged.  Set to 0 if mid-point and set to 1 if average
                                      !! value of Stokes drift over level. If advecting with
                                      !! Stokes transport, 1 is the correct approach.

  ! Surface Wave Dependent 1d/2d/3d vars
  real, allocatable, dimension(:), public :: &
       WaveNum_Cen        !< Wavenumber bands for read/coupled [m-1]
  real, allocatable, dimension(:), public :: &
       Freq_Cen           !< Frequency bands for read/coupled [s-1]
  real, allocatable, dimension(:), public :: &
       PrescribedSurfStkX !< Surface Stokes drift if prescribed [m s-1]
  real, allocatable, dimension(:), public :: &
       PrescribedSurfStkY !< Surface Stokes drift if prescribed [m s-1]
  real, allocatable, dimension(:,:,:), public :: &
       Us_x               !< 3d zonal Stokes drift profile [m s-1]
                          !! Horizontal -> U points
                          !! Vertical -> Mid-points
  real, allocatable, dimension(:,:,:), public :: &
       Us_y               !< 3d meridional Stokes drift profile [m s-1]
                          !! Horizontal -> V points
                          !! Vertical -> Mid-points
    real, allocatable, dimension(:,:,:), public :: &
       ddt_Us_x           !< 3d zonal Stokes drift profile [m s-1]
                          !! Horizontal -> U points
                          !! Vertical -> Mid-points
  real, allocatable, dimension(:,:,:), public :: &
       ddt_Us_y           !< 3d meridional Stokes drift profile [m s-1]
                          !! Horizontal -> V points
                          !! Vertical -> Mid-points
  real, allocatable, dimension(:,:), public :: &
       La_SL,&            !< SL Langmuir number (directionality factored later)
                          !! Horizontal -> H points
       La_Turb            !< Aligned Turbulent Langmuir number
                          !! Horizontal -> H points
  real, allocatable, dimension(:,:), public :: &
       US0_x              !< Surface Stokes Drift (zonal, m/s)
                          !! Horizontal -> U points
  real, allocatable, dimension(:,:), public :: &
       US0_y              !< Surface Stokes Drift (meridional, m/s)
                          !! Horizontal -> V points
  real, allocatable, dimension(:,:,:), public :: &
       STKx0              !< Stokes Drift spectrum (zonal, m/s)
                          !! Horizontal -> U points
                          !! 3rd dimension -> Freq/Wavenumber
  real, allocatable, dimension(:,:,:), public :: &
       STKy0              !< Stokes Drift spectrum (meridional, m/s)
                          !! Horizontal -> V points
                          !! 3rd dimension -> Freq/Wavenumber
  real, allocatable, dimension(:,:,:), public :: &
       KvS                !< Viscosity for Stokes Drift shear [Z2 T-1 ~> m2 s-1]

  ! Pointers to auxiliary fields
  type(diag_ctrl), pointer, public :: diag !< A structure that is used to regulate the
                                           !! timing of diagnostic output.

  !> An arbitrary lower-bound on the Langmuir number.  Run-time parameter.
  !! Not physical, but a limit against numerical artifacts.
  real :: La_min = 0.05

  ! Options related to specific Stokes drift choices.
  integer :: WaveMethod = -99 !< Options for including wave information
                          !! Valid (tested) choices are:
                          !!   0 - Test Profile
                          !!   1 - Surface Stokes Drift Bands
                          !!   2 - DHH85
                          !!   3 - LF17
                          !!   4 - WIND
                          !! -99 - No waves computed, but empirical Langmuir number used.

  ! Options For Test Prof
  Real    :: TP_STKX0, &! X-direction surface Stokes drift in test profile
             TP_STKY0, &! Y-direction surface Stokes drift in test profile
             TP_WVL ! Decay wavelength in Stokes drift test profile

  ! Options for DHH85 Profile
  logical :: WaveAgePeakFreq ! Flag to use Wave Age to get peak frequency
  logical :: StaticWaves ! Flag to use a temporally fixed Stokes drift profile
                         ! (e.g., to avoid recomputing the profile every time step in column mode)
  logical :: DHH85_Is_Set ! Flag to indicate that the profile is set in StaticWaves mode
  real    :: WaveAge ! Fixed value of WaveAge when using the DHH85 method.
  logical :: Wind_from_Tau ! Get wind from Tau or:
  real    :: WaveWindX ! Fixed value of Wind (optional) when using the DHH85 method
  real    :: WaveWindY ! Fixed value of Wind (optional) when using the DHH85 method

  ! Options if WaveMethod is Surface Stokes Drift Bands (1)
  integer, public :: NumBands =0 !< Number of wavenumber/frequency partitions to receive
                               !! This needs to match the number of bands provided
                               !! via either coupling or file.
  integer, public :: PartitionMode !< Method for partition mode (meant to check input)
                                   !! 0 - wavenumbers
                                   !! 1 - frequencies
                                   !! \todo Module variable! Move into a control structure.
  character(len=40)  :: SurfBandFileName !< Filename if using DataOverride
  logical :: dataoverrideisinitialized !< Flag for DataOverride Initialization

  ! Mathematical constant that may be defined elsewhere, but is set here from trig identities.
  real    :: PI


  ! Options for computing Langmuir number
  real :: LA_FracHBL !< Fraction of OSBL for averaging Langmuir number
  logical :: LA_Misalignment = .false. !< Flag to use misalignment in Langmuir number
                                     !! \todo Module variable! Move into a control structure.

  ! Diagnostic handles
  integer, public :: id_surfacestokes_x = -1 , id_surfacestokes_y = -1
  integer, public :: id_3dstokes_x = -1 , id_3dstokes_y = -1
  integer, public :: id_ddt_3dstokes_x = -1 , id_ddt_3dstokes_y = -1
  integer, public :: id_La_turb = -1
  integer, public :: id_PFu_Stokes = -1 , id_PFv_Stokes = -1



end type wave_parameters_CS

! Options not needed outside of this module

! This include declares and sets the variable "version".
#include "version_variable.h"

character(len=40)  :: mdl = "MOM_wave_interface" !< This module's name.

contains

!> Initializes parameters related to MOM_wave_interface
subroutine MOM_wave_interface_init(time, G, GV, US, param_file, CS, diag )
  type(time_type), target, intent(in)    :: Time       !< Model time
  type(ocean_grid_type),   intent(inout) :: G          !< Grid structure
  type(verticalGrid_type), intent(in)    :: GV         !< Vertical grid structure
  type(unit_scale_type),   intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Input parameter structure
  type(wave_parameters_CS), pointer      :: CS         !< Wave parameter control structure
  type(diag_ctrl), target, intent(inout) :: diag       !< Diagnostic Pointer
  ! Local variables
  ! I/O
  character*(13) :: TMPSTRING1,TMPSTRING2

  ! Dummy Check
  if (associated(CS)) then
     call MOM_error(FATAL, "wave_interface_init called with an associated"//&
                             "control structure.")
     return
  endif

  ! Allocate CS
  allocate(CS)

  ! Set pointers
  CS%diag => diag

  ! Initialize module pi constant
  CS%PI = 4.0*atan(1.0)

  ! Add any initializations needed here
  CS%DataOverrideIsInitialized = .false.
  CS%UseWaves = .true. ! The only way to get here is with UseWaves enabled.

  call log_version(param_file, mdl, version)

  ! Wave modified physics
  !  Presently these are all in research mode
  call get_param(param_file, mdl, "LAGRANGIAN_MIXING", CS%LagrangianMixing, &
       "Flag to use Lagrangian Mixing of momentum", units="", &
       Default=.false.)
  if (CS%LagrangianMixing) then
    ! Force Code Intervention
    call MOM_error(FATAL,"Should you be enabling Lagrangian Mixing? Code not ready.")
  endif
  call get_param(param_file, mdl, "STOKES_MIXING", CS%StokesMixing, &
       "Flag to use Stokes Mixing of momentum", units="", &
       Default=.false.)
  if (CS%StokesMixing) then
    ! Force Code Intervention
    call MOM_error(FATAL,"Should you be enabling Stokes Mixing? Code not ready.")
  endif
  call get_param(param_file, mdl, "STOKES_VF", CS%Stokes_VF, &
       "Flag to use Stokes vortex force", units="", &
       Default=.false.)
  call get_param(param_file, mdl, "STOKES_PGF", CS%Stokes_PGF, &
       "Flag to use Stokes pressure gradient force", units="", &
       Default=.false.)
  call get_param(param_file, mdl, "PASSIVE_STOKES_PGF", CS%Passive_Stokes_PGF, &
       "Flag to make Stokes pressure gradient force diagnostic only.", units="", &
       Default=.false.)
  call get_param(param_file, mdl, "STOKES_DDT", CS%Stokes_DDT, &
       "Flag to use Stokes d/dt", units="", &
       Default=.false.)

  ! Get Wave Method and write to integer WaveMethod
  call get_param(param_file,mdl,"WAVE_METHOD",TMPSTRING1,&
       "Choice of wave method, valid options include: \n"//&
       " TEST_PROFILE  - Prescribed from surface Stokes drift and a decay wavelength.\n"//&
       " SURFACE_BANDS - Computed from surface values and decay wavelengths.\n"//&
       " DHH85         - Uses Donelan et al. 1985 empirical wave spectrum. \n"//  &
       " LF17          - Infers Stokes drift profile from wind speed following Li and Fox-Kemper 2017.\n"//&
       " WIND          - A Simple relationship from the wind vector for testing.",&
       units='', default="EMPTY")
  select case (TRIM(TMPSTRING1))
  case ("EMPTY")! No Waves
    call MOM_error(FATAL, "wave_interface_init called with no specified WAVE_METHOD.")
  case ("TEST_PROFILE")! Test Profile
    CS%WaveMethod = TESTPROF
    call get_param(param_file,mdl,"TP_STKX_SURF",CS%TP_STKX0,&
         'Surface Stokes (x) for test profile',&
         units='m/s',default=0.1)
    call get_param(param_file,mdl,"TP_STKY_SURF",CS%TP_STKY0,&
         'Surface Stokes (y) for test profile',&
         units='m/s',default=0.0)
    call get_param(param_file,mdl,"TP_WVL",CS%TP_WVL,&
         units='m', default=50.0, scale=US%m_to_Z)
  case ("SURFACE_BANDS")! Surface Stokes Drift Bands
    CS%WaveMethod = SURFBANDS
    call get_param(param_file, mdl, "SURFBAND_SOURCE",TMPSTRING2,       &
       "Choice of SURFACE_BANDS data mode, valid options include: \n"// &
       " DATAOVERRIDE  - Read from NetCDF using FMS DataOverride. \n"// &
       " COUPLER       - Look for variables from coupler pass \n"//     &
       " INPUT         - Testing with fixed values.",                   &
       units='', default="EMPTY")
    select case (TRIM(TMPSTRING2))
    case ("EMPTY")! Default
      call MOM_error(FATAL, "wave_interface_init called with SURFACE_BANDS"//&
                           " but no SURFBAND_SOURCE.")
    case ("DATAOVERRIDE")! Using Data Override
      CS%DataSource = DATAOVR
      call get_param(param_file, mdl, "SURFBAND_FILENAME", CS%SurfBandFileName,&
           "Filename of surface Stokes drift input band data.", default="StkSpec.nc")
    case ("COUPLER")! Reserved for coupling
      CS%DataSource = COUPLER
    case ("INPUT")! A method to input the Stokes band (globally uniform)
      CS%DataSource = Input
      call get_param(param_file,mdl,"SURFBAND_NB",CS%NumBands,&
         "Prescribe number of wavenumber bands for Stokes drift. "//&
         "Make sure this is consistnet w/ WAVENUMBERS, STOKES_X, and "//&
         "STOKES_Y, there are no safety checks in the code.",&
         units='', default=1)
      allocate( CS%WaveNum_Cen(1:CS%NumBands) )
      CS%WaveNum_Cen(:) = 0.0
      allocate( CS%PrescribedSurfStkX(1:CS%NumBands))
      CS%PrescribedSurfStkX(:) = 0.0
      allocate( CS%PrescribedSurfStkY(1:CS%NumBands))
      CS%PrescribedSurfStkY(:) = 0.0
      allocate( CS%STKx0(G%isdB:G%iedB,G%jsd:G%jed,1:CS%NumBands))
      CS%STKx0(:,:,:) = 0.0
      allocate( CS%STKy0(G%isd:G%ied,G%jsdB:G%jedB,1:CS%NumBands))
      CS%STKy0(:,:,:) = 0.0
      CS%partitionmode=0
      call get_param(param_file,mdl,"SURFBAND_WAVENUMBERS",CS%WaveNum_Cen,      &
           "Central wavenumbers for surface Stokes drift bands.",units='rad/m', &
           default=0.12566)
      call get_param(param_file,mdl,"SURFBAND_STOKES_X",CS%PrescribedSurfStkX,      &
           "X-direction surface Stokes drift for bands.",units='m/s', &
           default=0.15)
      call get_param(param_file,mdl,"SURFBAND_STOKES_Y",CS%PrescribedSurfStkY,      &
           "Y-direction surface Stokes drift for bands.",units='m/s', &
           default=0.0)
   case default! No method provided
     call MOM_error(FATAL,'Check WAVE_METHOD.')
   end select

  case ("DHH85")!Donelan et al., 1985 spectrum
    CS%WaveMethod = DHH85
    call get_param(param_file,mdl,"DHH85_WAVEAGE_FP",CS%WaveAgePeakFreq,   &
         "Choose true to use waveage in peak frequency.", &
         units='', default=.false.)
    if (CS%WaveAgePeakFreq) then
      call get_param(param_file,mdl,"DHH85_AGE",CS%WaveAge,   &
           "Wave Age for DHH85 spectrum.", &
           units='', default=1.2)
    endif
    call get_param(param_file,mdl,"DHH85_WIND_FROM_TAU",CS%Wind_from_Tau,&
         "Choose true to estimate wind from tau.", &
         units='', default=.false.)
    if (.not.CS%Wind_from_Tau) then
      call get_param(param_file,mdl,"DHH85_WIND_X",CS%WaveWindx,&
           "Wind speed for DHH85 spectrum.",&
           units='', default=10.0)
      call get_param(param_file,mdl,"DHH85_WIND_Y",CS%WaveWindy,&
           "Wind speed for DHH85 spectrum.", &
           units='', default=0.0)
    endif
    call get_param(param_file,mdl,"STATIC_DHH85",CS%StaticWaves,&
         "Flag to disable updating DHH85 Stokes drift.",&
          default=.false.)
  case ("LF17")!Li and Fox-Kemper 17 wind-sea Langmuir number
    CS%WaveMethod = LF17
  case ("WIND")!Wind speed based
    CS%WaveMethod = WIND
  case default
      call MOM_error(FATAL,'Check WAVE_METHOD.')
  end select

  ! Langmuir number Options
  call get_param(param_file, mdl, "LA_DEPTH_RATIO", CS%LA_FracHBL,&
         "The depth (normalized by BLD) to average Stokes drift over in "//&
         "Langmuir number calculation, where La = sqrt(ust/Stokes).",&
         units="nondim",default=0.04)
  call get_param(param_file, mdl, "LA_MISALIGNMENT", CS%LA_Misalignment,&
         "Flag (logical) if using misalignment bt shear and waves in LA",&
         default=.false.)
  call get_param(param_file, mdl, "MIN_LANGMUIR", CS%La_min,&
         "A minimum value for all Langmuir numbers that is not physical, "//&
         "but is likely only encountered when the wind is very small and "//&
         "therefore its effects should be mostly benign.",units="nondim",&
         default=0.05)

  ! Allocate and initialize
  ! a. Stokes driftProfiles
  allocate(CS%Us_x(G%isdB:G%IedB,G%jsd:G%jed,G%ke))
  CS%Us_x(:,:,:) = 0.0
  allocate(CS%Us_y(G%isd:G%Ied,G%jsdB:G%jedB,G%ke))
  CS%Us_y(:,:,:) = 0.0
    allocate(CS%ddt_Us_x(G%isdB:G%IedB,G%jsd:G%jed,G%ke))
  CS%ddt_Us_x(:,:,:) = 0.0
  allocate(CS%ddt_Us_y(G%isd:G%Ied,G%jsdB:G%jedB,G%ke))
  CS%ddt_Us_y(:,:,:) = 0.0
  ! b. Surface Values
  allocate(CS%US0_x(G%isdB:G%iedB,G%jsd:G%jed))
  CS%US0_x(:,:) = 0.0
  allocate(CS%US0_y(G%isd:G%ied,G%jsdB:G%jedB))
  CS%US0_y(:,:) = 0.0
  ! c. Langmuir number
  allocate(CS%La_SL(G%isc:G%iec,G%jsc:G%jec))
  allocate(CS%La_turb(G%isc:G%iec,G%jsc:G%jec))
  CS%La_SL(:,:) = 0.0
  CS%La_turb (:,:) = 0.0
  ! d. Viscosity for Stokes drift
  if (CS%StokesMixing) then
    allocate(CS%KvS(G%isd:G%Ied,G%jsd:G%jed,G%ke))
    CS%KvS(:,:,:) = 0.0
  endif

  ! Initialize Wave related outputs
  CS%id_surfacestokes_x = register_diag_field('ocean_model','uos_Stokes', &
       CS%diag%axesCu1,Time,'Surface Stokes drift (meridional)','m s-1')
  CS%id_surfacestokes_y = register_diag_field('ocean_model','vos_Stokes', &
       CS%diag%axesCv1,Time,'Surface Stokes drift (zonal)','m s-1')
  CS%id_3dstokes_y = register_diag_field('ocean_model','v_Stokes', &
       CS%diag%axesCvL,Time,'Stokes drift (meridional)','m s-1')
  CS%id_3dstokes_x = register_diag_field('ocean_model','u_Stokes', &
       CS%diag%axesCuL,Time,'Stokes drift (zonal)','m s-1')
  CS%id_ddt_3dstokes_y = register_diag_field('ocean_model','dvdt_Stokes', &
       CS%diag%axesCvL,Time,'d/dt Stokes drift (meridional)','m s-2')
  CS%id_ddt_3dstokes_x = register_diag_field('ocean_model','dudt_Stokes', &
       CS%diag%axesCuL,Time,'d/dt Stokes drift (zonal)','m s-2')
  CS%id_PFv_Stokes = register_diag_field('ocean_model','PFv_Stokes', &
       CS%diag%axesCvL,Time,'PF from Stokes drift (meridional)','m s-2')
  CS%id_PFu_Stokes = register_diag_field('ocean_model','PFu_Stokes', &
       CS%diag%axesCuL,Time,'PF from Stokes drift (zonal)','m s-2')
  CS%id_La_turb = register_diag_field('ocean_model','La_turbulent',&
       CS%diag%axesT1,Time,'Surface (turbulent) Langmuir number','nondim')

  return
end subroutine MOM_wave_interface_init

!> A 'lite' init subroutine to initialize a few inputs needed if using wave information
!! with the wind-speed dependent Stokes drift formulation of LF17
subroutine MOM_wave_interface_init_lite(CS, diag, param_file)
  type(wave_parameters_CS), pointer :: CS !< Wave parameter control structure
  type(diag_ctrl), target, intent(inout) :: diag !< Diagnostic Pointer
  type(param_file_type), intent(in) :: param_file !< Input parameter structure
  !/
  character*(13) :: TMPSTRING1
  logical :: StatisticalWaves

  ! Dummy Check
  if (associated(CS)) then
    call MOM_error(FATAL,"wave_interface_init_lite called with an associated"//&
                         "control structure.")
    return
  endif

  ! Allocate CS
  allocate(CS)

  ! Set pointers
  CS%diag => diag

  ! Initialize module pi constant
  CS%PI = 4.0*atan(1.0)

  ! Add any initializations needed here
  CS%UseWaves = .false.

  ! Langmuir number Options
  call get_param(param_file, mdl, "LA_DEPTH_RATIO", CS%LA_FracHBL,&
       "The depth (normalized by BLD) to average Stokes drift over in "//&
       "Langmuir number calculation, where La = sqrt(ust/Stokes).",&
       units="nondim",default=0.04)

  ! Check if using LA_LI2016 (UseWaves remains false)
  call get_param(param_file,mdl,"USE_LA_LI2016",StatisticalWaves,&
                 do_not_log=.true.,default=.false.)
  if (StatisticalWaves) then
    CS%WaveMethod = LF17
    CS%PI=4.0*atan(1.0)
  else
    CS%WaveMethod = NULL_WaveMethod
  end if

  return
end subroutine MOM_wave_interface_init_lite

!> Subroutine that handles updating of surface wave/Stokes drift related properties
!!  at the driver level, particularly for handling coupler and file read updates.
subroutine Update_Surface_Waves(G, GV, US, Day, dt, CS)
  type(wave_parameters_CS), pointer    :: CS  !< Wave parameter Control structure
  type(ocean_grid_type), intent(inout) :: G   !< Grid structure
  type(verticalGrid_type), intent(in)  :: GV  !< Vertical grid structure
  type(unit_scale_type),   intent(in)  :: US   !< A dimensional unit scaling type
  type(time_type),         intent(in)  :: Day !< Current model time
  type(time_type),         intent(in)  :: dt  !< Timestep as a time-type
  ! Local variables
  integer :: ii, jj, kk, b
  type(time_type) :: Day_Center

  ! Computing central time of time step
  Day_Center = Day + DT/2

  if (CS%WaveMethod==SURFBANDS) then
    if (CS%DataSource==DATAOVR) then
      call Surface_Bands_by_data_override(day_center, G, GV, US, CS)
    elseif (CS%DataSource==Coupler) then
      ! Reserve for coupler hooks
    elseif (CS%DataSource==Input) then
      do b=1,CS%NumBands
        do II=G%isdB,G%iedB ; do jj=G%jsd,G%jed
          CS%STKx0(II,jj,b) = CS%PrescribedSurfStkX(b)
        enddo ; enddo
        do ii=G%isd,G%ied ; do JJ=G%jsdB, G%jedB
          CS%STKY0(ii,JJ,b) = CS%PrescribedSurfStkY(b)
        enddo ; enddo
      enddo
    endif
  endif

  return
end subroutine Update_Surface_Waves

!> Constructs the Stokes Drift profile on the model grid based on
!! desired coupling options
subroutine Update_Stokes_Drift(G, GV, US, CS, h, forces, dt)
  type(wave_parameters_CS),  pointer     :: CS    !< Wave parameter Control structure
  type(ocean_grid_type),   intent(inout) :: G     !< Grid structure
  type(verticalGrid_type), intent(in)    :: GV    !< Vertical grid structure
  type(unit_scale_type),   intent(in)    :: US    !< A dimensional unit scaling type
  type(mech_forcing), intent(in)         :: forces !< Container with mechanical forcing terms
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
       intent(in)    :: h     !< Thickness [H ~> m or kg m-2]
  real :: dt
  ! Local Variables
  real    :: Top, MidPoint, Bottom, one_cm
  real    :: DecayScale
  real    :: CMN_FAC, WN, UStokes
  real    :: La
  integer :: ii, jj, kk, b, iim1, jjm1
  real :: ustary, ustarx, tauy, taux, u10, ustar
  real :: idt


  idt = 1./dt
  one_cm = 0.01*US%m_to_Z

  CS%ddt_us_x(:,:,:) = CS%US_x(:,:,:)
  CS%ddt_us_y(:,:,:) = CS%US_y(:,:,:)
  
  ! 1. If Test Profile Option is chosen
  !    Computing mid-point value from surface value and decay wavelength
  if (CS%WaveMethod==TESTPROF) then
    DecayScale = 4.*CS%PI / CS%TP_WVL
    do II = G%iscB,G%iecB ; do jj = G%jsc,G%jec
      Bottom = 0.0
      MidPoint = 0.0
      do kk = 1,G%ke
        Top = Bottom
        MidPoint = Bottom - GV%H_to_Z*0.25*(h(II,jj,kk)+h(IIm1,jj,kk))
        Bottom = Bottom - GV%H_to_Z*0.5*(h(II,jj,kk)+h(IIm1,jj,kk))
        CS%Us_x(II,jj,kk) = CS%TP_STKX0*exp(MidPoint*DecayScale)
      enddo
    enddo ; enddo
    do ii = G%isc,G%iec ; do JJ = G%jscB,G%jecB
      Bottom = 0.0
      MidPoint = 0.0
      do kk = 1,G%ke
        Top = Bottom
        MidPoint = Bottom - GV%H_to_Z*0.25*(h(ii,JJ,kk)+h(ii,JJm1,kk))
        Bottom = Bottom - GV%H_to_Z*0.5*(h(ii,JJ,kk)+h(ii,JJm1,kk))
        CS%Us_y(ii,JJ,kk) = CS%TP_STKY0*exp(MidPoint*DecayScale)
      enddo
    enddo ; enddo
  ! 2. If Surface Bands is chosen
  !    In wavenumber mode compute integral for layer averaged Stokes drift.
  !    In frequency mode compuate value at midpoint.
  elseif (CS%WaveMethod==SURFBANDS) then
    ! Computing Surface Stokes drift
    do b = 1,CS%NumBands
      if (CS%PartitionMode==0) then
        ! In wavenumber we are averaging over (small) level
        CMN_FAC = (1.0-exp(-one_cm*2.*CS%WaveNum_Cen(b))) / &
                  (one_cm*2.*CS%WaveNum_Cen(b))
      elseif (CS%PartitionMode==1) then
         ! In frequency we are not averaging over level and taking top
        CMN_FAC = 1.0
      endif
      ! X-direction
      do II = G%iscB,G%iecB ; do jj = G%jsc,G%jec
         CS%US0_x(II,jj) = CS%US0_x(II,jj) + CS%STKx0(II,jj,b)*CMN_FAC
      enddo ; enddo
      ! Y-direction
      do ii = G%isc,G%iec ; do JJ = G%jscB,G%jecB
        CS%US0_y(ii,JJ) = CS%US0_y(ii,JJ) + CS%STKy0(ii,JJ,b)*CMN_FAC
      enddo ; enddo
    enddo

    ! 2. Second compute the level averaged Stokes drift
    !
    ! X-direction
    do II = G%iscB,G%iecB ; do jj = G%jsc,G%jec
      bottom = 0.0
      do kk = 1,G%ke
        Top = Bottom
        MidPoint = Bottom - GV%H_to_Z*0.25*(h(II,jj,kk)+h(II-1,jj,kk))
        Bottom = Bottom - GV%H_to_Z*0.5*(h(II,jj,kk)+h(II-1,jj,kk))
        do b = 1,CS%NumBands
          if (CS%PartitionMode==0) then
            ! In wavenumber we are averaging over level
            CMN_FAC = (exp(Top*2.*CS%WaveNum_Cen(b))-exp(Bottom*2.*CS%WaveNum_Cen(b)))&
                      / ((Top-Bottom)*(2.*CS%WaveNum_Cen(b)))
          elseif (CS%PartitionMode==1) then
            if (CS%StkLevelMode==0) then
              ! Take the value at the midpoint
              CMN_FAC = exp(MidPoint*2.*(2.*CS%PI*CS%Freq_Cen(b)*US%T_to_s)**2/(US%L_to_Z**2*GV%g_Earth))
            elseif (CS%StkLevelMode==1) then
              ! Use a numerical integration and then
              ! divide by layer thickness
              WN = (2.*CS%PI*CS%Freq_Cen(b)*US%T_to_s)**2 / (US%L_to_Z**2*GV%g_Earth) !bgr bug-fix missing g
              CMN_FAC = (exp(2.*WN*Top)-exp(2.*WN*Bottom)) / (2.*WN*(Top-Bottom))
            endif
          endif
          CS%US_x(II,jj,kk) = CS%US_x(II,jj,kk) + CS%STKx0(II,jj,b)*CMN_FAC
        enddo
      enddo
    enddo ; enddo
    ! Y direction
    do ii = G%isc,G%iec ;do JJ = G%jscB,G%jecB
      ! Compute the level averages.
      bottom = 0.0
      do kk = 1,G%ke
        Top = Bottom
        MidPoint = Bottom - GV%H_to_Z*0.25*(h(ii,JJ,kk)+h(ii,JJ-1,kk))
        Bottom = Bottom - GV%H_to_Z*0.5*(h(ii,JJ,kk)+h(ii,JJ-1,kk))
        do b = 1,CS%NumBands
          if (CS%PartitionMode==0) then
            ! In wavenumber we are averaging over level
            CMN_FAC = (exp(Top*2.*CS%WaveNum_Cen(b)) - &
                       exp(Bottom*2.*CS%WaveNum_Cen(b))) / &
                      ((Top-Bottom)*(2.*CS%WaveNum_Cen(b)))
          elseif (CS%PartitionMode==1) then
            if (CS%StkLevelMode==0) then
              ! Take the value at the midpoint
              CMN_FAC = exp(MidPoint*2.*(2.*CS%PI*CS%Freq_Cen(b)*US%T_to_s)**2/(US%L_to_Z**2*GV%g_Earth))
            elseif (CS%StkLevelMode==1) then
              ! Use a numerical integration and then
              ! divide by layer thickness
              WN = (2.*CS%PI*CS%Freq_Cen(b)*US%T_to_s)**2 / (US%L_to_Z**2*GV%g_Earth)
              CMN_FAC = (exp(2.*WN*Top)-exp(2.*WN*Bottom)) / (2.*WN*(Top-Bottom))
            endif
          endif
          CS%US_y(ii,JJ,kk) = CS%US_y(ii,JJ,kk) + CS%STKy0(ii,JJ,b)*CMN_FAC
        enddo
      enddo
    enddo ; enddo
  elseif (CS%WaveMethod==DHH85) then
    if (.not.(CS%StaticWaves .and. CS%DHH85_is_set)) then
      do II = G%iscB,G%iecB ; do jj = G%jsc,G%jec
        bottom = 0.0
        do kk = 1,G%ke
          Top = Bottom
          MidPoint = Bottom - GV%H_to_Z*0.25*(h(II,jj,kk)+h(II-1,jj,kk))
          Bottom = Bottom - GV%H_to_Z*0.5*(h(II,jj,kk)+h(II-1,jj,kk))
          ! 4-point avg from V to U points
          tauy = (G%mask2dCv(II+1,jj)*forces%tauy(II+1,jj) +&
                 G%mask2dCv(II+1,jj-1)*forces%tauy(II+1,jj-1) +&
                 G%mask2dCv(II,jj)*forces%tauy(II,jj) +&
                 G%mask2dCv(II,jj-1)*forces%tauy(II,jj-1)) /&
                 (G%mask2dCv(II+1,jj)+G%mask2dCv(II+1,jj-1)+&
                 G%mask2dCv(II,jj)+G%mask2dCv(II,jj-1))
          ustar = sqrt(sqrt(forces%taux(II,jj)**2+tauy**2)/GV%rho0)
          ustarX = ustar*forces%taux(II,jj)/sqrt(forces%taux(II,jj)**2+tauy**2)
          if (abs(ustarX).gt.0.0) then
            call ust_2_u10_coare3p5(US%Z_to_m*US%s_to_T*ustarX*sqrt(US%R_to_kg_m3*GV%Rho0/1.225),&
                 u10, GV, US)
          else
            u10 = 0.0
          endif
          !u10 = CS%WaveWind -> Now computed from stress
          call DHH85_mid(CS, GV, US, u10, CS%WaveAge, MidPoint, UStokes)
          CS%US_x(II,jj,kk) = UStokes
        enddo
      enddo; enddo
      do ii = G%isc,G%iec ; do JJ = G%jscB,G%jecB
        bottom = 0.0
        do kk = 1,G%ke
          Top = Bottom
          MidPoint = Bottom - GV%H_to_Z*0.25*(h(II,jj,kk)+h(II-1,jj,kk))
          Bottom = Bottom - GV%H_to_Z*0.5*(h(II,jj,kk)+h(II-1,jj,kk))
          ! 4-point avg from U to V points
          taux = (G%mask2dCu(ii,JJ+1)*forces%taux(ii,JJ+1) +&
                 G%mask2dCu(ii-1,JJ+1)*forces%taux(ii-1,JJ+1) +&
                 G%mask2dCu(ii,JJ)*forces%taux(ii,JJ) +&
                 G%mask2dCu(ii-1,JJ)*forces%taux(ii-1,JJ)) /&
                 (G%mask2dCu(ii,JJ+1)+G%mask2dCu(ii-1,JJ+1)+&
                 G%mask2dCu(ii,JJ)+G%mask2dCu(ii-1,JJ))
          ustar = sqrt(sqrt(forces%tauy(ii,JJ)**2*G%mask2dCv(ii,JJ)+taux**2)/GV%rho0)
          ustarY = ustar*forces%tauy(ii,JJ)/sqrt(forces%tauy(ii,JJ)**2+taux**2)
          if (abs(ustarY).gt.0.0) then
            call ust_2_u10_coare3p5(US%Z_to_m*US%s_to_T*ustarY*sqrt(US%R_to_kg_m3*GV%Rho0/1.225),&
                                    u10, GV, US)
          else
            u10 = 0.0
          endif
          !u10 = CS%WaveWind -> Now computed from stress
          call DHH85_mid(CS, GV, US, u10, CS%WaveAge, MidPoint, UStokes)
          CS%US_y(ii,JJ,kk) = UStokes
        enddo
      enddo; enddo
      CS%us0_y(:,:) = CS%us_y(:,:,1)
      CS%us0_x(:,:) = CS%us_x(:,:,1)
      CS%DHH85_is_set = .true.
    endif
  elseif (CS%WaveMethod==WIND) then
    do II = G%iscB,G%iecB ; do jj = G%jsc,G%jec
      ! 4-point avg from V to U points
      tauy = (G%mask2dCv(II+1,jj)*forces%tauy(II+1,jj) +&
             G%mask2dCv(II+1,jj-1)*forces%tauy(II+1,jj-1) +&
             G%mask2dCv(II,jj)*forces%tauy(II,jj) +&
             G%mask2dCv(II,jj-1)*forces%tauy(II,jj-1)) /&
             (G%mask2dCv(II+1,jj)+G%mask2dCv(II+1,jj-1)+&
             G%mask2dCv(II,jj)+G%mask2dCv(II,jj-1))
      ustar = sqrt(sqrt(forces%taux(II,jj)**2+tauy**2)/GV%rho0)
      ustarX = ustar*forces%taux(II,jj)/sqrt(forces%taux(II,jj)**2+tauy**2)
      if (abs(ustarX).gt.0.0) then
        call ust_2_u10_coare3p5(US%Z_to_m*US%s_to_T*abs(ustarX)*sqrt(US%R_to_kg_m3*GV%Rho0/1.225),&
             u10, GV, US)
      else
        u10 = 0.0
      endif
      u10 = sign(u10,ustarX)
      bottom = 0.0
      CS%us0_x(II,jj) = u10*0.03
      do kk = 1,G%ke
        Top = Bottom
        MidPoint = Bottom - GV%H_to_Z*0.25*(h(II,jj,kk)+h(II-1,jj,kk))
        Bottom = Bottom - GV%H_to_Z*0.5*(h(II,jj,kk)+h(II-1,jj,kk))
        CS%US_x(II,jj,kk) = u10*0.03*exp(MidPoint*2.*3.14/50.)
      enddo
    enddo; enddo
    do ii = G%isc,G%iec ; do JJ = G%jscB,G%jecB
      ! 4-point avg from U to V points
      taux = (G%mask2dCu(ii,JJ+1)*forces%taux(ii,JJ+1) +&
              G%mask2dCu(ii-1,JJ+1)*forces%taux(ii-1,JJ+1) +&
              G%mask2dCu(ii,JJ)*forces%taux(ii,JJ) +&
              G%mask2dCu(ii-1,JJ)*forces%taux(ii-1,JJ)) /&
              (G%mask2dCu(ii,JJ+1)+G%mask2dCu(ii-1,JJ+1)+&
              G%mask2dCu(ii,JJ)+G%mask2dCu(ii-1,JJ))
      ustar = sqrt(sqrt(forces%tauy(ii,JJ)**2*G%mask2dCv(ii,JJ)+taux**2)/GV%rho0)
      ustarY = ustar*forces%tauy(ii,JJ)/sqrt(forces%tauy(ii,JJ)**2+taux**2)
      if (abs(ustary).gt.0.0) then
        call ust_2_u10_coare3p5(US%Z_to_m*US%s_to_T*abs(ustarY)*sqrt(US%R_to_kg_m3*GV%Rho0/1.225),&
             u10, GV, US)
      else
        u10 = 0.0
      endif
      u10 = sign(u10,ustarY)
      !u10 = CS%WaveWind -> Now computed from stress
      bottom = 0.0
      CS%us0_y(II,jj) = u10*0.03
      do kk = 1,G%ke
        Top = Bottom
        MidPoint = Bottom - GV%H_to_Z*0.25*(h(II,jj,kk)+h(II-1,jj,kk))
        Bottom = Bottom - GV%H_to_Z*0.5*(h(II,jj,kk)+h(II-1,jj,kk))
        CS%US_y(ii,JJ,kk) = u10*0.03*exp(MidPoint*2.*3.14/50.)
      enddo
    enddo; enddo
  else! Keep this else, fallback to 0 Stokes drift
    do kk= 1,G%ke
      do II = G%iscB,G%iecB ; do jj = G%jsc,G%jec
        CS%Us_x(II,jj,kk) = 0.
      enddo ; enddo
      do ii = G%isc,G%iec ; do JJ = G%jscB,G%jecB
        CS%Us_y(ii,JJ,kk) = 0.
      enddo ; enddo
    enddo
  endif

  call pass_vector(CS%US_x,CS%Us_y, G%Domain, To_ALL)
  call pass_vector(CS%US0_x,CS%Us0_y, G%Domain, To_ALL)

  ! Turbulent Langmuir number is computed here and available to use anywhere.
  ! SL Langmuir number requires mixing layer depth, and therefore is computed
  ! in the routine it is needed by (e.g. KPP or ePBL).
  do ii = G%isc,G%iec ; do jj = G%jsc, G%jec
    Top = h(ii,jj,1)*GV%H_to_Z
    call get_Langmuir_Number( La, G, GV, US, Top, forces%ustar(ii,jj), ii, jj, &
                              H(ii,jj,:),Override_MA=.false.,CS=CS)
    CS%La_turb(ii,jj) = La
  enddo ; enddo

  CS%ddt_us_x(:,:,:) = (CS%US_x(:,:,:) - CS%ddt_us_x(:,:,:)) * idt
  CS%ddt_us_y(:,:,:) = (CS%US_y(:,:,:) - CS%ddt_us_y(:,:,:)) * idt
  
  ! Output any desired quantities
  if (CS%id_surfacestokes_y>0) &
    call post_data(CS%id_surfacestokes_y, CS%us0_y, CS%diag)
  if (CS%id_surfacestokes_x>0) &
    call post_data(CS%id_surfacestokes_x, CS%us0_x, CS%diag)
  if (CS%id_3dstokes_y>0) &
    call post_data(CS%id_3dstokes_y, CS%us_y, CS%diag)
  if (CS%id_ddt_3dstokes_x>0) &
       call post_data(CS%id_ddt_3dstokes_x, CS%ddt_us_x, CS%diag)
  if (CS%id_ddt_3dstokes_y>0) &
    call post_data(CS%id_ddt_3dstokes_y, CS%ddt_us_y, CS%diag)
  if (CS%id_3dstokes_x>0) &
    call post_data(CS%id_3dstokes_x, CS%us_x, CS%diag)
  if (CS%id_La_turb>0) &
    call post_data(CS%id_La_turb, CS%La_turb, CS%diag)

end subroutine Update_Stokes_Drift

!> A subroutine to fill the Stokes drift from a NetCDF file
!! using the data_override procedures.
subroutine Surface_Bands_by_data_override(day_center, G, GV, US, CS)
  use NETCDF
  type(time_type),          intent(in) :: day_center !< Center of timestep
  type(wave_parameters_CS), pointer    :: CS         !< Wave structure
  type(ocean_grid_type), intent(inout) :: G          !< Grid structure
  type(verticalGrid_type),  intent(in) :: GV         !< Vertical grid structure
  type(unit_scale_type),    intent(in) :: US         !< A dimensional unit scaling type
  ! Local variables
  real    :: temp_x(SZI_(G),SZJ_(G)) ! Pseudo-zonal Stokes drift of band at h-points [m s-1]
  real    :: temp_y(SZI_(G),SZJ_(G)) ! Psuedo-meridional Stokes drift of band at h-points [m s-1]
  real    :: Top, MidPoint
  integer :: b
  integer :: i, j
  integer, dimension(4) :: start, counter, dims, dim_id
  character(len=12)  :: dim_name(4)
  character(20) :: varname, varread1, varread2
  integer :: rcode_fr, rcode_wn, ncid, varid_fr, varid_wn, id, ndims

  if (.not.CS%dataOverrideIsInitialized) then
    call data_override_init(Ocean_domain_in=G%Domain%mpp_domain)
    CS%dataOverrideIsInitialized = .true.

    ! Read in number of wavenumber bands in file to set number to be read in
    ! Hardcoded filename/variables
    varread1 = 'wavenumber' !Old method gives wavenumber
    varread2 = 'frequency'  !New method gives frequency
    rcode_wn = NF90_OPEN(trim(CS%SurfBandFileName), NF90_NOWRITE, ncid)
    if (rcode_wn /= 0) then
      call MOM_error(FATAL,"error opening file "//trim(CS%SurfBandFileName)//&
            " in MOM_wave_interface.")
    endif

    ! Check if rcode_wn or rcode_fr is 0 (checks if input has wavenumber or frequency)
    rcode_wn = NF90_INQ_VARID(ncid, varread1, varid_wn)
    rcode_fr = NF90_INQ_VARID(ncid, varread2, varid_fr)

    if (rcode_wn /= 0 .and. rcode_fr /= 0) then
      call MOM_error(FATAL,"error finding variable "//trim(varread1)//&
         " or "//trim(varread2)//" in file "//trim(CS%SurfBandFileName)//" in MOM_wave_interface.")

    elseif (rcode_wn == 0) then
      ! wavenumbers found:
      CS%PartitionMode = 0
      rcode_wn = NF90_INQUIRE_VARIABLE(ncid, varid_wn, ndims=ndims, &
           dimids=dims)
      if (rcode_wn /= 0) then
        call MOM_error(FATAL, &
             'error inquiring dimensions MOM_wave_interface.')
      endif
      rcode_wn = NF90_INQUIRE_DIMENSION(ncid, dims(1), dim_name(1), len=id)
      if (rcode_wn /= 0) then
        call MOM_error(FATAL,"error reading dimension 1 data for "// &
             trim(varread1)//" in file "// trim(CS%SurfBandFileName)//          &
             " in MOM_wave_interface.")
      endif
      rcode_wn = NF90_INQ_VARID(ncid, dim_name(1), dim_id(1))
      if (rcode_wn /= 0) then
        call MOM_error(FATAL,"error finding variable "//trim(dim_name(1))//&
          " in file "//trim(CS%SurfBandFileName)//" in MOM_wave_interace.")
      endif
      ! Allocating size of wavenumber bins
      allocate( CS%WaveNum_Cen(1:id) )
      CS%WaveNum_Cen(:) = 0.0
    elseif (rcode_fr == 0) then
      ! frequencies found:
      CS%PartitionMode = 1
      rcode_fr = NF90_INQUIRE_VARIABLE(ncid, varid_fr, ndims=ndims, &
           dimids=dims)
      if (rcode_fr /= 0) then
        call MOM_error(FATAL,&
             'error inquiring dimensions MOM_wave_interface.')
      endif
      rcode_fr = NF90_INQUIRE_DIMENSION(ncid, dims(1), dim_name(1), len=id)
      if (rcode_fr /= 0) then
        call MOM_error(FATAL,"error reading dimension 1 data for "// &
             trim(varread2)//" in file "// trim(CS%SurfBandFileName)// &
             " in MOM_wave_interface.")
      endif
      rcode_fr = NF90_INQ_VARID(ncid, dim_name(1), dim_id(1))
      if (rcode_fr /= 0) then
        call MOM_error(FATAL,"error finding variable "//trim(dim_name(1))//&
             " in file "//trim(CS%SurfBandFileName)//" in MOM_wave_interace.")
      endif
      ! Allocating size of frequency bins
      allocate( CS%Freq_Cen(1:id) )
      CS%Freq_Cen(:) = 0.0
      ! Allocating size of wavenumber bins
      allocate( CS%WaveNum_Cen(1:id) )
      CS%WaveNum_Cen(:) = 0.0
      allocate( CS%STKx0(G%isdB:G%iedB,G%jsd:G%jed,1:id))
      CS%STKx0(:,:,:) = 0.0
      allocate( CS%STKy0(G%isd:G%ied,G%jsdB:G%jedB,1:id))
      CS%STKy0(:,:,:) = 0.0
    endif

    ! Reading wavenumber bins/Frequencies
    start(:) = 1 ! Set all start to 1
    counter(:) = 1 ! Set all counter to 1
    counter(1) = id ! Set counter(1) to id (number of frequency bins)
    if (CS%PartitionMode==0) then
      rcode_wn = NF90_GET_VAR(ncid, dim_id(1), CS%WaveNum_Cen, start, counter)
      if (rcode_wn /= 0) then
        call MOM_error(FATAL,&
             "error reading dimension 1 values for var_name "// &
             trim(varread1)//",dim_name "//trim(dim_name(1))//  &
             " in file "// trim(CS%SurfBandFileName)//" in MOM_wave_interface")
      endif
      CS%NumBands = ID
      do B = 1,CS%NumBands ; CS%WaveNum_Cen(b) = US%Z_to_m*CS%WaveNum_Cen(b) ; enddo
    elseif (CS%PartitionMode==1) then
      rcode_fr = NF90_GET_VAR(ncid, dim_id(1), CS%Freq_Cen, start, counter)
      if (rcode_fr /= 0) then
        call MOM_error(FATAL,&
             "error reading dimension 1 values for var_name "// &
             trim(varread2)//",dim_name "//trim(dim_name(1))//  &
             " in file "// trim(CS%SurfBandFileName)//" in MOM_wave_interface")
      endif
      CS%NumBands = ID
      do B = 1,CS%NumBands
        CS%WaveNum_Cen(b) = (2.*CS%PI*CS%Freq_Cen(b)*US%T_to_s)**2 / (US%L_to_Z**2*GV%g_Earth)
      enddo
    endif

  endif

  do b = 1,CS%NumBands
    temp_x(:,:) = 0.0
    temp_y(:,:) = 0.0
    varname = '                    '
    write(varname,"(A3,I0)")'Usx',b
    call data_override('OCN',trim(varname), temp_x, day_center)
    varname = '                    '
    write(varname,'(A3,I0)')'Usy',b
    call data_override('OCN',trim(varname), temp_y, day_center)
    ! Disperse into halo on h-grid
    call pass_vector(temp_x, temp_y, G%Domain, To_All, AGRID)
    !Filter land mask values.  Needed?
    do j = G%jsd,G%jed ; do i = G%Isd,G%Ied
      ! Assume anything with magnitude more than 10 M/S is mask contiminated?
      if (abs(temp_x(i,j)) > 10. .or. abs(temp_y(i,j)) > 10.) then
        ! Assume land-mask and zero out
        temp_x(i,j) = 0.0
        temp_y(i,j) = 0.0
      endif
    enddo ; enddo

    ! Interpolate to u/v grids
    do j = G%jsc,G%jec
      do I = G%IscB,G%IecB
        CS%STKx0(I,j,b) = 0.5 * (temp_x(i,j) + temp_x(i+1,j))
      enddo
    enddo
    do J = G%JscB,G%JecB
      do i = G%isc,G%iec
        CS%STKy0(i,J,b) = 0.5 * (temp_y(i,j) + temp_y(i,j+1))
      enddo
    enddo

  enddo !Closes b-loop
  ! Disperse into halo on u/v grids (moved from 2d to 3d call)
  call pass_vector(CS%STKx0(:,:,:),CS%STKy0(:,:,:), G%Domain, To_ALL)

  return
end subroutine Surface_Bands_by_data_override

!> Interface to get Langmuir number based on options stored in wave structure
subroutine get_Langmuir_Number( LA, G, GV, US, HBL, ustar, i, j, &
                                H, U_H, V_H, Override_MA, CS )
  type(ocean_grid_type), intent(in) :: G  !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV !< Ocean vertical grid structure
  type(unit_scale_type), intent(in) :: US !< A dimensional unit scaling type
  integer, intent(in) :: i !< Meridional index of h-point
  integer, intent(in) :: j !< Zonal index of h-point
  real, intent(in) :: ustar !< Friction velocity [Z T-1 ~> m s-1].
  real, intent(in) :: HBL !< (Positive) thickness of boundary layer [Z ~> m].
  logical, optional, intent(in) :: Override_MA !< Override to use misalignment in LA
                                               !! calculation. This can be used if diagnostic
                                               !! LA outputs are desired that are different than
                                               !! those used by the dynamical model.
  real, dimension(SZK_(GV)), optional,&
       intent(in) :: H !< Grid layer thickness [H ~> m or kg m-2]
  real, dimension(SZK_(GV)), optional,&
       intent(in) :: U_H !< Zonal velocity at H point [m s-1]
  real, dimension(SZK_(GV)), optional,&
       intent(in) :: V_H !< Meridional velocity at H point [m s-1]
  type(Wave_parameters_CS),&
       pointer :: CS !< Surface wave control structure.
  real, intent(out) :: LA !< Langmuir number

!Local Variables
  real :: Top, bottom, midpoint
  real :: Dpt_LASL, ShearDirection, WaveDirection
  real :: LA_STKx, LA_STKy, LA_STK ! Stokes velocities in [m s-1]
  logical :: ContinueLoop, USE_MA
  real, dimension(SZK_(G)) :: US_H, VS_H
  real, dimension(CS%NumBands) :: StkBand_X, StkBand_Y
  integer :: KK, BB

 ! Compute averaging depth for Stokes drift (negative), bounded at 0.1 m minimum
  Dpt_LASL = min(-0.1*US%m_to_Z, -CS%LA_FracHBL*HBL)

  USE_MA = CS%LA_Misalignment
  if (present(Override_MA)) USE_MA = Override_MA

  ! If requesting to use misalignment in the Langmuir number compute the Shear Direction
  if (USE_MA) then
    if (.not.(present(H).and.present(U_H).and.present(V_H))) then
      call MOM_error(Fatal,'Get_LA_waves requested to consider misalignment.')
    endif
    ContinueLoop = .true.
    bottom = 0.0
    do kk = 1,G%ke
      Top = Bottom
      MidPoint = Bottom + GV%H_to_Z*0.5*h(kk)
      Bottom = Bottom + GV%H_to_Z*h(kk)
      if (MidPoint > Dpt_LASL .and. kk > 1 .and. ContinueLoop) then
        ShearDirection = atan2(V_H(1)-V_H(kk),U_H(1)-U_H(kk))
        ContinueLoop = .false.
      endif
    enddo
  endif

  if (CS%WaveMethod==TESTPROF) then
    do kk = 1,G%ke
      US_H(kk) = 0.5*(CS%US_X(I,j,kk)+CS%US_X(I-1,j,kk))
      VS_H(kk) = 0.5*(CS%US_Y(i,J,kk)+CS%US_Y(i,J-1,kk))
    enddo
    call Get_SL_Average_Prof( GV, Dpt_LASL, H, US_H, LA_STKx)
    call Get_SL_Average_Prof( GV, Dpt_LASL, H, VS_H, LA_STKy)
    LA_STK = sqrt(LA_STKX*LA_STKX+LA_STKY*LA_STKY)
  elseif (CS%WaveMethod==SURFBANDS) then
    do bb = 1,CS%NumBands
      StkBand_X(bb) = 0.5*(CS%STKx0(I,j,bb)+CS%STKx0(I-1,j,bb))
      StkBand_Y(bb) = 0.5*(CS%STKy0(i,J,bb)+CS%STKy0(i,J-1,bb))
    enddo
    call Get_SL_Average_Band(GV, Dpt_LASL, CS%NumBands, CS%WaveNum_Cen, StkBand_X, LA_STKx )
    call Get_SL_Average_Band(GV, Dpt_LASL, CS%NumBands, CS%WaveNum_Cen, StkBand_Y, LA_STKy )
    LA_STK = sqrt(LA_STKX**2 + LA_STKY**2)
  elseif (CS%WaveMethod==DHH85 .or. CS%WaveMethod==WIND) then
    ! Temporarily integrating profile rather than spectrum for simplicity
    do kk = 1,GV%ke
      US_H(kk) = 0.5*(CS%US_X(I,j,kk)+CS%US_X(I-1,j,kk))
      VS_H(kk) = 0.5*(CS%US_Y(i,J,kk)+CS%US_Y(i,J-1,kk))
    enddo
    call Get_SL_Average_Prof( GV, Dpt_LASL, H, US_H, LA_STKx)
    call Get_SL_Average_Prof( GV, Dpt_LASL, H, VS_H, LA_STKy)
    LA_STK = sqrt(LA_STKX**2 + LA_STKY**2)
  elseif (CS%WaveMethod==LF17) then
    call get_StokesSL_LiFoxKemper(CS, ustar, hbl*CS%LA_FracHBL, GV, US, LA_STK, LA)
  elseif (CS%WaveMethod==Null_WaveMethod) then
    call MOM_error(FATAL, "Get_Langmuir_number called without defining a WaveMethod. "//&
                          "Suggest to make sure USE_LT is set/overridden to False or "//&
                          "choose a wave method (or set USE_LA_LI2016 to use statistical "//&
                          "waves.")
  endif

  if (.not.(CS%WaveMethod==LF17)) then
    ! This is an arbitrary lower bound on Langmuir number.
    ! We shouldn't expect values lower than this, but
    ! there is also no good reason to cap it here other then
    ! to prevent large enhancements in unconstrained parts of
    ! the curve fit parameterizations.
    ! Note the dimensional constant background Stokes velocity of 10^-10 m s-1.
    LA = max(CS%La_min, sqrt(US%Z_to_m*US%s_to_T*ustar / (LA_STK+1.e-10)))
  endif

  if (Use_MA) then
    WaveDirection = atan2(LA_STKy, LA_STKx)
    LA = LA / sqrt(max(1.e-8, cos( WaveDirection - ShearDirection)))
  endif

  return
end subroutine get_Langmuir_Number

!> Get SL averaged Stokes drift from Li/FK 17 method
!!
!! Original description:
!! - This function returns the enhancement factor, given the 10-meter
!!   wind [m s-1], friction velocity [m s-1] and the boundary layer depth [m].
!!
!! Update (Jan/25):
!! - Converted from function to subroutine, now returns Langmuir number.
!! - Computs 10m wind internally, so only ustar and hbl need passed to
!!   subroutine.
!!
!! Qing Li, 160606
!! - BGR port from CVMix to MOM6 Jan/25/2017
!! - BGR change output to LA from Efactor
!! - BGR remove u10 input
!! - BGR note: fixed parameter values should be changed to "get_params"
subroutine get_StokesSL_LiFoxKemper(CS, ustar, hbl, GV, US, UStokes_SL, LA)
  type(Wave_parameters_CS), pointer :: CS !< Surface wave related control structure.
  real, intent(in) :: ustar !< water-side surface friction velocity [Z T-1 ~> m s-1].
  real, intent(in) :: hbl !< boundary layer depth [Z ~> m].
  type(verticalGrid_type), intent(in) :: GV !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in) :: US !< A dimensional unit scaling type
  real, intent(out) :: UStokes_SL !< Surface layer averaged Stokes drift [m s-1]
  real, intent(out) :: LA !< Langmuir number
  ! Local variables
  ! parameters
  real, parameter :: &
       ! ratio of U19.5 to U10 (Holthuijsen, 2007)
       u19p5_to_u10 = 1.075, &
       ! ratio of mean frequency to peak frequency for
       ! Pierson-Moskowitz spectrum (Webb, 2011)
       fm_into_fp = 1.296, &
       ! ratio of surface Stokes drift to U10
       us_to_u10 = 0.0162, &
       ! loss ratio of Stokes transport
       r_loss = 0.667
  real :: UStokes, hm0, fm, fp, vstokes, kphil, kstar
  real :: z0, z0i, r1, r2, r3, r4, tmp, lasl_sqr_i
  real :: u10

  if (ustar > 0.0) then
    ! Computing u10 based on u_star and COARE 3.5 relationships
    call ust_2_u10_coare3p5(US%Z_to_m*US%s_to_T*ustar*sqrt(US%R_to_kg_m3*GV%Rho0/1.225), u10, GV, US)
    ! surface Stokes drift
    UStokes = us_to_u10*u10
    !
    ! significant wave height from Pierson-Moskowitz
    ! spectrum (Bouws, 1998)
    hm0 = 0.0246 *u10**2
    !
    ! peak frequency (PM, Bouws, 1998)
    tmp = 2.0 * CS%PI * u19p5_to_u10 * u10
    fp = 0.877 * US%L_T_to_m_s**2*US%m_to_Z * GV%g_Earth / tmp
    !
    ! mean frequency
    fm = fm_into_fp * fp
    !
    ! total Stokes transport (a factor r_loss is applied to account
    !  for the effect of directional spreading, multidirectional waves
    !  and the use of PM peak frequency and PM significant wave height
    !  on estimating the Stokes transport)
    vstokes = 0.125 * CS%PI * r_loss * fm * hm0**2
    !
    ! the general peak wavenumber for Phillips' spectrum
    ! (Breivik et al., 2016) with correction of directional spreading
    kphil = 0.176 * UStokes / vstokes
    !
    ! surface layer averaged Stokes dirft with Stokes drift profile
    ! estimated from Phillips' spectrum (Breivik et al., 2016)
    ! the directional spreading effect from Webb and Fox-Kemper, 2015
    ! is also included
    kstar = kphil * 2.56
    ! surface layer
    z0 = abs(US%Z_to_m*hbl)
    z0i = 1.0 / z0
    ! term 1 to 4
    r1 = ( 0.151 / kphil * z0i -0.84 ) * &
         ( 1.0 - exp(-2.0 * kphil * z0) )
    r2 = -( 0.84 + 0.0591 / kphil * z0i ) * &
         sqrt( 2.0 * CS%PI * kphil * z0 ) * &
         erfc( sqrt( 2.0 * kphil * z0 ) )
    r3 = ( 0.0632 / kstar * z0i + 0.125 ) * &
         (1.0 - exp(-2.0 * kstar * z0) )
    r4 = ( 0.125 + 0.0946 / kstar * z0i ) * &
         sqrt( 2.0 * CS%PI *kstar * z0) * &
         erfc( sqrt( 2.0 * kstar * z0 ) )
    UStokes_sl = UStokes * (0.715 + r1 + r2 + r3 + r4)
    LA = sqrt(US%Z_to_m*US%s_to_T*ustar / UStokes_sl)
  else
    UStokes_sl = 0.0
    LA=1.e8
  endif

end subroutine Get_StokesSL_LiFoxKemper

!> Get SL Averaged Stokes drift from a Stokes drift Profile
subroutine Get_SL_Average_Prof( GV, AvgDepth, H, Profile, Average )
  type(verticalGrid_type),  &
       intent(in)   :: GV       !< Ocean vertical grid structure
  real, intent(in)  :: AvgDepth !< Depth to average over (negative) [Z ~> m].
  real, dimension(SZK_(GV)), &
       intent(in)   :: H        !< Grid thickness [H ~> m or kg m-2]
  real, dimension(SZK_(GV)), &
       intent(in)   :: Profile  !< Profile of quantity to be averaged [arbitrary]
                                !! (used here for Stokes drift)
  real, intent(out) :: Average  !< Output quantity averaged over depth AvgDepth [arbitrary]
                                !! (used here for Stokes drift)
  !Local variables
  real :: top, midpoint, bottom ! Depths, negative downward [Z ~> m].
  real :: Sum
  integer :: kk

  ! Initializing sum
  Sum = 0.0

  ! Integrate
  bottom = 0.0
  do kk = 1, GV%ke
    Top = Bottom
    MidPoint = Bottom - GV%H_to_Z * 0.5*h(kk)
    Bottom = Bottom - GV%H_to_Z * h(kk)
    if (AvgDepth < Bottom) then ! The whole cell is within H_LA
      Sum = Sum + Profile(kk) * (GV%H_to_Z * H(kk))
    elseif (AvgDepth < Top) then ! A partial cell is within H_LA
      Sum = Sum + Profile(kk) * (Top-AvgDepth)
      exit
    else
      exit
    endif
  enddo

  ! Divide by AvgDepth or the depth in the column, whichever is smaller.
  if (abs(AvgDepth) <= abs(Bottom)) then
    Average = Sum / abs(AvgDepth)
  elseif (abs(Bottom) > 0.0) then
    Average = Sum / abs(Bottom)
  else
    Average = 0.0
  endif

end subroutine Get_SL_Average_Prof

!> Get SL averaged Stokes drift from the banded Spectrum method
subroutine Get_SL_Average_Band( GV, AvgDepth, NB, WaveNumbers, SurfStokes, Average )
  type(verticalGrid_type),  &
       intent(in)     :: GV          !< Ocean vertical grid
  real, intent(in)    :: AvgDepth    !< Depth to average over [Z ~> m].
  integer, intent(in) :: NB          !< Number of bands used
  real, dimension(NB), &
       intent(in)     :: WaveNumbers !< Wavenumber corresponding to each band [Z-1 ~> m-1]
  real, dimension(NB), &
       intent(in)     :: SurfStokes  !< Surface Stokes drift for each band [m s-1]
  real, intent(out)   :: Average     !< Output average Stokes drift over depth AvgDepth [m s-1]

  ! Local variables
  integer :: bb

  ! Loop over bands
  Average = 0.0
  do bb = 1, NB
    ! Factor includes analytical integration of e(2kz)
    !  - divided by (-H_LA) to get average from integral.
    Average = Average + SurfStokes(BB) * &
              (1.-EXP(-abs(AvgDepth * 2.0 * WaveNumbers(BB)))) / &
              abs(AvgDepth * 2.0 * WaveNumbers(BB))
  enddo

  return
end subroutine Get_SL_Average_Band

!> Compute the Stokes drift at a given depth
!!
!! Taken from Qing Li (Brown)
!! use for comparing MOM6 simulation to his LES
!! computed at z mid point (I think) and not depth averaged.
!! Should be fine to integrate in frequency from 0.1 to sqrt(-0.2*grav*2pi/dz
subroutine DHH85_mid(CS, GV, US, WaveAge, U10, zpt, UStokes)
  type(wave_parameters_CS), pointer :: CS !< Wave parameter control structure
  type(verticalGrid_type), intent(in) :: GV !< Ocean vertical grid
  type(unit_scale_type), intent(in) :: US !< A dimensional unit scaling type
  real, intent(in) :: zpt !< Depth to get Stokes drift [Z ~> m].
  real, intent(in) :: WaveAge  !< Wave-age (non-dim)
  real, intent(in) :: U10   !< Wind [m s-1]
  real, intent(out) :: UStokes !< Stokes drift [m s-1]
  !
  real :: ann, Bnn, Snn, Cnn, Dnn
  real :: omega_peak, omega, domega
  real :: omega_min, omega_max, wavespec, Stokes
  real :: g_Earth ! Gravitational acceleration [m s-2]
  integer :: Nomega, OI

  g_Earth = US%L_T_to_m_s**2*US%m_to_Z * GV%g_Earth

  !/
  omega_min = 0.1 ! Hz
  ! Cut off at 30cm for now...
  omega_max = 10. ! ~sqrt(0.2*g_Earth*2*pi/0.3)
  NOmega = 1000
  domega = (omega_max-omega_min)/real(NOmega)

  !
  if (CS%WaveAgePeakFreq) then
    omega_peak = g_Earth / (WaveAge * u10)
  else
    omega_peak = 2. * CS%pi * 0.13 * g_Earth / U10
  endif
  !/
  Ann = 0.006 * WaveAge**(-0.55)
  Bnn = 1.0
  Snn = 0.08 * (1.0 + 4.0 * WaveAge**3)
  Cnn = 1.7
  if (WaveAge < 1.) then
    Cnn = Cnn - 6.0*log10(WaveAge)
  endif
  !/
  UStokes = 0.0
  omega = omega_min + 0.5*domega
  do oi = 1,nomega-1
    Dnn = exp ( -0.5 * (omega-omega_peak)**2 / (Snn**2 * omega_peak**2) )
    ! wavespec units = m2s
    wavespec = (Ann * g_Earth**2 / (omega_peak*omega**4 ) ) * &
               exp(-bnn*(omega_peak/omega)**4)*Cnn**Dnn
    ! Stokes units m  (multiply by frequency range for units of m/s)
    Stokes = 2.0 * wavespec * omega**3 * &
         exp( 2.0 * omega**2 * US%Z_to_m*zpt / g_Earth) / g_Earth
    UStokes = UStokes + Stokes*domega
    omega = omega + domega
  enddo

  return
end subroutine DHH85_mid

!> Explicit solver for Stokes mixing.
!! Still in development do not use.
subroutine StokesMixing(G, GV, dt, h, u, v, Waves )
  type(ocean_grid_type), &
       intent(in)    :: G     !< Ocean grid
  type(verticalGrid_type), &
       intent(in)    :: GV    !< Ocean vertical grid
  real, intent(in)   :: dt    !< Time step of MOM6 [T ~> s] for explicit solver
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),&
       intent(in)    :: h     !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), &
       intent(inout) :: u     !< Velocity i-component [m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
       intent(inout) :: v     !< Velocity j-component [m s-1]
  type(Wave_parameters_CS), &
       pointer       :: Waves !< Surface wave related control structure.
  ! Local variables
  real :: dTauUp, dTauDn ! Vertical momentum fluxes [Z T-1 m s-1]
  real :: h_Lay  ! The layer thickness at a velocity point [Z ~> m].
  integer :: i,j,k

! This is a template to think about down-Stokes mixing.
! This is not ready for use...

  do k = 1, G%ke
    do j = G%jsc, G%jec
      do I = G%iscB, G%iecB
        h_lay = GV%H_to_Z*0.5*(h(i,j,k)+h(i+1,j,k))
        dTauUp = 0.0
        if (k > 1) &
          dTauUp = 0.5*(waves%Kvs(i,j,k)+waves%Kvs(i+1,j,k)) * &
               (waves%us_x(i,j,k-1)-waves%us_x(i,j,k)) / &
               (0.5*(h_lay + GV%H_to_Z*0.5*(h(i,j,k-1)+h(i+1,j,k-1)) ))
        dTauDn = 0.0
        if (k < G%ke-1) &
          dTauDn = 0.5*(waves%Kvs(i,j,k+1)+waves%Kvs(i+1,j,k+1)) * &
               (waves%us_x(i,j,k)-waves%us_x(i,j,k+1)) / &
               (0.5*(h_lay + GV%H_to_Z*0.5*(h(i,j,k+1)+h(i+1,j,k+1)) ))
        u(i,j,k) = u(i,j,k) + dt * (dTauUp-dTauDn) / h_Lay
      enddo
    enddo
  enddo

  do k = 1, G%ke
    do J = G%jscB, G%jecB
      do i = G%isc, G%iec
        h_Lay = GV%H_to_Z*0.5*(h(i,j,k)+h(i,j+1,k))
        dTauUp = 0.
        if (k > 1) &
          dTauUp = 0.5*(waves%Kvs(i,j,k)+waves%Kvs(i,j+1,k)) * &
               (waves%us_y(i,j,k-1)-waves%us_y(i,j,k)) / &
               (0.5*(h_lay + GV%H_to_Z*0.5*(h(i,j,k-1)+h(i,j+1,k-1)) ))
        dTauDn = 0.0
        if (k < G%ke-1) &
          dTauDn =0.5*(waves%Kvs(i,j,k+1)+waves%Kvs(i,j+1,k+1)) * &
               (waves%us_y(i,j,k)-waves%us_y(i,j,k+1)) / &
               (0.5*(h_lay + GV%H_to_Z*0.5*(h(i,j,k+1)+h(i,j+1,k+1)) ))
        v(i,J,k) = v(i,J,k) + dt * (dTauUp-dTauDn) / h_Lay
      enddo
    enddo
  enddo

end subroutine StokesMixing

!> Computes tendency due to Stokes pressure gradient force
subroutine Stokes_PGF(G, GV, h, u, v, PFu_Stokes, PFv_Stokes, CS )
  type(ocean_grid_type), &
       intent(in)    :: G     !< Ocean grid
  type(verticalGrid_type), &
       intent(in)    :: GV    !< Ocean vertical grid
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),&
       intent(in)    :: h       !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), &
       intent(in) :: u          !< Velocity i-component [m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
       intent(in) :: v          !< Velocity j-component [m s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), &
       intent(out) :: PFu_Stokes !< PGF Stokes-shear i-component [L T-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
       intent(out) :: PFv_Stokes !< PGF Stokes-shear j-component [m s-1]
  type(Wave_parameters_CS), &
       pointer       :: CS !< Surface wave related control structure.

  ! Local variables
  real :: P_Stokes_l, P_Stokes_r ! Contribution of Stokes shear to pressure (left/right index) [L2 T-2 ~> m2 s-2]
  real :: u_l, u_r, v_l, v_r ! Velocity components
  real :: dUs_dz_l, dUs_dz_r ! Vertical derivative of zonal Stokes drift (left/right index) [T-1 ~> s-1]
  real :: dVs_dz_l, dVs_dz_r ! Vertical derivative of meridional Stokes drift (left/right index) [T-1 ~> s-1]
  real :: z_top_l, z_top_r   ! The height of the top of the cell (left/right index) [Z ~> m].
  real :: z_mid_l, z_mid_r   ! The height of the middle of the cell (left/right index) [Z ~> m].
  real :: h_l, h_r   ! The thickness of the cell (left/right index) [Z ~> m].
  real :: wavenum,TwoKexpL, TwoKexpR !TMP DELETE THIS
  integer :: i,j,k

  ! Comput the Stokes contribution to the pressure gradient force

  PFu_Stokes(:,:,:) = 0.0
  PFv_Stokes(:,:,:) = 0.0

  wavenum = (2.*3.14)/50.
  do j = G%jsc, G%jec ; do I = G%iscB, G%iecB
    if (G%mask2dCu(I,j)>0.5) then
      z_top_l = 0.0
      z_top_r = 0.0
      P_Stokes_l = 0.0
      P_Stokes_r = 0.0
      do k = 1, G%ke
        h_l = h(i,j,k)
        h_r = h(i+1,j,k)
        z_mid_l = z_top_l - 0.5*h_l
        z_mid_r = z_top_r - 0.5*h_r
        TwoKexpL = (2.*wavenum)*exp(2*wavenum*z_mid_l)
        TwoKexpR = (2.*wavenum)*exp(2*wavenum*z_mid_r)
        !UL -> I-1 & I, j
        !UR -> I & I+1, j 
        !VL -> i, J-1 & J
        !VR -> i+1, J-1 & J
        dUs_dz_l = TwoKexpL*0.5 * &
                   (CS%Us0_x(I-1,j)*G%mask2dCu(I-1,j) + &
                   CS%Us0_x(I,j)*G%mask2dCu(I,j))
        dUs_dz_r = TwoKexpR*0.5 * &
                   (CS%Us0_x(I,j)*G%mask2dCu(I,j) + &
                   CS%Us0_x(I+1,j)*G%mask2dCu(I+1,j))
        dVs_dz_l = TwoKexpL*0.5 * &
                   (CS%Us0_y(i,J-1)*G%mask2dCv(i,J-1) + &
                   CS%Us0_y(i,J)*G%mask2dCv(i,J))
        dVs_dz_r = TwoKexpR*0.5 * &
                   (CS%Us0_y(i+1,J-1)*G%mask2dCv(i+1,J-1) + &
                    CS%Us0_y(i+1,J)*G%mask2dCv(i+1,J))
        u_l = 0.5*(u(I-1,j,k)*G%mask2dCu(I-1,j) + &
                   u(I,j,k)*G%mask2dCu(I,j))
        u_r = 0.5*(u(I,j,k)*G%mask2dCu(I,j) + &
                   u(I+1,j,k)*G%mask2dCu(I+1,j))
        v_l = 0.5*(v(i,J-1,k)*G%mask2dCv(i,J-1) + &
                   v(i,J,k)*G%mask2dCv(i,J))
        v_r = 0.5*(v(i+1,J-1,k)*G%mask2dCv(i+1,J-1) + &
                   v(i+1,J,k)*G%mask2dCv(i+1,J))
        if (G%mask2dT(i,j)>0.5) &
             P_Stokes_l = P_Stokes_l + h_l*(dUs_dz_l*u_l+dVs_dz_l*v_l)
        if (G%mask2dT(i+1,j)>0.5) &
             P_Stokes_r = P_Stokes_r + h_r*(dUs_dz_r*u_r+dVs_dz_r*v_r)
        PFu_Stokes(I,j,k) = (P_Stokes_r - P_Stokes_l)*G%IdxCu(I,j)
        z_top_l = z_top_l - h_l
        z_top_r = z_top_r - h_r
      enddo
    endif
  enddo ;  enddo
  do J = G%jscB, G%jecB ; do i = G%isc, G%iec
    if (G%mask2dCv(i,J)>0.5) then
      z_top_l = 0.0
      z_top_r = 0.0
      P_Stokes_l = 0.0
      P_Stokes_r = 0.0
      do k = 1, G%ke
        h_l = h(i,j,k)
        h_r = h(i,j+1,k)
        z_mid_l = z_top_l - 0.5*h_l
        z_mid_r = z_top_r - 0.5*h_r
        TwoKexpL = (2.*wavenum)*exp(2*wavenum*z_mid_l)
        TwoKexpR = (2.*wavenum)*exp(2*wavenum*z_mid_r)
        !UL -> I-1 & I, j
        !UR -> I-1 & I, j+1
        !VL -> i, J & J-1
        !VR -> i, J & J+1
        dUs_dz_l = TwoKexpL*0.5 * &
                   (CS%Us0_x(I-1,j)*G%mask2dCu(I-1,j) + &
                   CS%Us0_x(I,j)*G%mask2dCu(I,j))
        dUs_dz_r = TwoKexpR*0.5 * &
                   (CS%Us0_x(I-1,j+1)*G%mask2dCu(I-1,j+1) + &
                   CS%Us0_x(I,j+1)*G%mask2dCu(I,j+1))
        dVs_dz_l = TwoKexpL*0.5 * &
                   (CS%Us0_y(i,J-1)*G%mask2dCv(i,J-1) + &
                   CS%Us0_y(i,J)*G%mask2dCv(i,J))
        dVs_dz_r = TwoKexpR*0.5 * &
                   (CS%Us0_y(i,J)*G%mask2dCv(i,J) + &
                    CS%Us0_y(i,J+1)*G%mask2dCv(i,J+1))
        u_l = 0.5*(u(I-1,j,k)*G%mask2dCu(I-1,j) + &
                   u(I,j,k)*G%mask2dCu(I,j))
        u_r = 0.5*(u(I-1,j+1,k)*G%mask2dCu(I-1,j+1) + &
                   u(I,j+1,k)*G%mask2dCu(I,j+1))
        v_l = 0.5*(v(i,J-1,k)*G%mask2dCv(i,J-1) + &
                   v(i,J,k)*G%mask2dCv(i,J))
        v_r = 0.5*(v(i,J,k)*G%mask2dCv(i,J) + &
                   v(i,J+1,k)*G%mask2dCv(i,J+1))
        if (G%mask2dT(i,j)>0.5) &
             P_Stokes_l = P_Stokes_l + h_l*(dUs_dz_l*u_l+dVs_dz_l*v_l)
        if (G%mask2dT(i,j+1)>0.5) &
             P_Stokes_r = P_Stokes_r + h_r*(dUs_dz_r*u_r+dVs_dz_r*v_r)
        PFv_Stokes(i,J,k) = (P_Stokes_r - P_Stokes_l)*G%IdyCv(i,J)
        z_top_l = z_top_l - h_l
        z_top_r = z_top_r - h_r
      enddo
    endif
  enddo ; enddo

  if (CS%id_PFv_Stokes>0) &
    call post_data(CS%id_PFv_Stokes, PFv_Stokes, CS%diag)
  if (CS%id_PFu_Stokes>0) &
    call post_data(CS%id_PFu_Stokes, PFu_Stokes, CS%diag)
  
end subroutine Stokes_PGF

!> Solver to add Coriolis-Stokes to model
!! Still in development and not meant for general use.
!! Can be activated (with code intervention) for LES comparison
!! CHECK THAT RIGHT TIMESTEP IS PASSED IF YOU USE THIS**
!!
!! Not accessed in the standard code.
subroutine CoriolisStokes(G, GV, DT, h, u, v, CS, US)
  type(ocean_grid_type), &
       intent(in)    :: G     !< Ocean grid
  type(verticalGrid_type), &
       intent(in)   :: GV     !< Ocean vertical grid
  real, intent(in)  :: Dt     !< Time step of MOM6 [s] CHECK IF PASSING RIGHT TIMESTEP
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  &
       intent(in)    :: h     !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), &
       intent(inout) :: u     !< Velocity i-component [m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
       intent(inout) :: v     !< Velocity j-component [m s-1]
  type(Wave_parameters_CS), &
       pointer       :: CS !< Surface wave related control structure.
  type(unit_scale_type),   intent(in) :: US     !< A dimensional unit scaling type
  ! Local variables
  real :: DVel ! A rescaled velocity change [m s-1 T-1 ~> m s-2]
  integer :: i,j,k

  do k = 1, G%ke
    do j = G%jsc, G%jec
      do I = G%iscB, G%iecB
        DVel = 0.25*(CS%us_y(i,j+1,k)+CS%us_y(i-1,j+1,k))*G%CoriolisBu(i,j+1) + &
               0.25*(CS%us_y(i,j,k)+CS%us_y(i-1,j,k))*G%CoriolisBu(i,j)
        u(I,j,k) = u(I,j,k) + DVEL*US%s_to_T*DT
      enddo
    enddo
  enddo

  do k = 1, G%ke
    do J = G%jscB, G%jecB
      do i = G%isc, G%iec
        DVel = 0.25*(CS%us_x(i+1,j,k)+CS%us_x(i+1,j-1,k))*G%CoriolisBu(i+1,j) + &
               0.25*(CS%us_x(i,j,k)+CS%us_x(i,j-1,k))*G%CoriolisBu(i,j)
        v(i,J,k) = v(i,j,k) - DVEL*US%s_to_T*DT
      enddo
    enddo
  enddo
end subroutine CoriolisStokes

!> Computes wind speed from ustar_air based on COARE 3.5 Cd relationship
!! Probably doesn't belong in this module, but it is used here to estimate
!! wind speed for wind-wave relationships.  Should be a fine way to estimate
!! the neutral wind-speed as written here.
subroutine ust_2_u10_coare3p5(USTair, U10, GV, US)
  real, intent(in)                    :: USTair !< Wind friction velocity [m s-1]
  real, intent(out)                   :: U10    !< 10-m neutral wind speed [m s-1]
  type(verticalGrid_type), intent(in) :: GV     !< vertical grid type
  type(unit_scale_type),   intent(in) :: US     !< A dimensional unit scaling type

  ! Local variables
  real, parameter :: vonkar = 0.4 ! Should access a get_param von karman
  real, parameter :: nu=1e-6 ! Should access a get_param air-viscosity
  real, parameter :: z0sm_coef=0.11
  real :: z0sm, z0, z0rough, u10a, alpha, CD

  integer :: CT

  ! Uses empirical formula for z0 to convert ustar_air to u10 based on the
  !  COARE 3.5 paper (Edson et al., 2013)
  ! alpha=m*U10+b
  ! Note in Edson et al. 2013, eq. 13 m is given as 0.017.  However,
  ! m=0.0017 reproduces the curve in their figure 6.

  z0sm = z0sm_coef * nu * US%m_to_Z / USTair !Compute z0smooth from ustar guess
  u10 = USTair/sqrt(0.001)  !Guess for u10, 0.001 is non-dimensional drag
  u10a = 1000

  CT=0
  do while (abs(u10a/u10-1.) > 0.001)
    CT=CT+1
    u10a = u10
    alpha = min(0.028, 0.0017 * u10 - 0.005)
    z0rough = alpha * (US%m_s_to_L_T*USTair)**2 / GV%g_Earth ! Compute z0rough from ustar guess
    z0 = z0sm + z0rough
    CD = ( vonkar / log(10.*US%m_to_Z / z0) )**2 ! Compute CD from derived roughness
    u10 = USTair/sqrt(CD)  ! Compute new u10 from derived CD, while loop
                           ! ends and checks for convergence...CT counter
                           ! makes sure loop doesn't run away if function
                           ! doesn't converge.  This code was produced offline
                           ! and converged rapidly (e.g. 2 cycles)
                           ! for ustar=0.0001:0.0001:10.
    if (CT>20) then
      u10 = USTair/sqrt(0.0015) ! I don't expect to get here, but just
                              !  in case it will output a reasonable value.
      exit
    endif
  enddo
  return
end subroutine ust_2_u10_coare3p5

!> Clear pointers, deallocate memory
subroutine Waves_end(CS)
  type(wave_parameters_CS), pointer :: CS !< Control structure

  if (allocated(CS%WaveNum_Cen)) deallocate( CS%WaveNum_Cen )
  if (allocated(CS%Freq_Cen))    deallocate( CS%Freq_Cen )
  if (allocated(CS%Us_x))        deallocate( CS%Us_x )
  if (allocated(CS%Us_y))        deallocate( CS%Us_y )
  if (allocated(CS%La_SL))       deallocate( CS%La_SL )
  if (allocated(CS%La_turb))     deallocate( CS%La_turb )
  if (allocated(CS%STKx0))       deallocate( CS%STKx0 )
  if (allocated(CS%STKy0))       deallocate( CS%STKy0 )
  if (allocated(CS%KvS))         deallocate( CS%KvS )
  if (allocated(CS%Us0_y))       deallocate( CS%Us0_y )
  if (allocated(CS%Us0_x))       deallocate( CS%Us0_x )

  deallocate( CS )

  return
end subroutine Waves_end

!> \namespace  mom_wave_interface
!!
!! \author Brandon Reichl, 2018.
!!
!! This module should be moved as wave coupling progresses and
!! likely will should mirror the iceberg or sea-ice model set-up.
!!
!! This module is meant to contain the routines to read in and
!! interpret surface wave data for MOM6. In its original form, the
!! capabilities include setting the Stokes drift in the model (from a
!! variety of sources including prescribed, empirical, and input
!! files).  In short order, the plan is to also ammend the subroutine
!! to accept Stokes drift information from an external coupler.
!! Eventually, it will be necessary to break this file apart so that
!! general wave information may be stored in the control structure
!! and the Stokes drift effect can be isolated from processes such as
!! sea-state dependent momentum fluxes, gas fluxes, and other wave
!! related air-sea interaction and boundary layer phenomenon.
!!
!! The Stokes drift are stored on the C-grid with the conventional
!! protocol to interpolate to the h-grid to compute Langmuir number,
!! the primary quantity needed for Langmuir turbulence
!! parameterizations in both the ePBL and KPP approach.  This module
!! also computes full 3d Stokes drift profiles, which will be useful
!! if second-order type boundary layer parameterizations are
!! implemented (perhaps via GOTM, work in progress).

end module MOM_wave_interface
