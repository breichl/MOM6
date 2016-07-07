!> Initial conditions and forcing for the single column model (SCM) CVmix
!! test set.
module MOM_wave_interface

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_forcing_type, only : forcing, allocate_forcing_type
use MOM_grid, only : ocean_grid_type
use MOM_verticalgrid, only: verticalGrid_type
use MOM_safe_alloc, only : safe_alloc_ptr
use MOM_time_manager, only : time_type, operator(+), operator(/), get_time,&
                             time_type_to_real
use MOM_variables, only : thermo_var_ptrs, surface
use data_override_mod, only : data_override_init, data_override
implicit none ; private

#include <MOM_memory.h>

public MOM_wave_interface_init
public Import_Stokes_Drift


!> Container for wave related parameters
type, public:: wave_parameters_CS ;
private
  logical :: UseWaves    !< True to Compute Wave parameters
  logical, public :: LagrangianMixing !If Stokes drift is present and viscous mixing
                              ! should be applied to Lagrangian current
  logical, public :: LangmuirEnhanceW   !If turbulent velocity scales should be
                                        ! modified due to presence of Langmuir
                                        ! mixing.
  logical, public :: LangmuirEnhanceVt2 !If unresolved turbulent velocity scale
                                        ! should be modified due to presence
                                        ! of Langmuir mixing.
  logical, public :: LangmuirEnhanceK   !If diffusivity/viscosity should be
                                        ! modified due to presence of Langmuir
                                        ! mixing
  logical, public :: StokesShearInRIb   !If Stokes drift should be included
                                        ! in current shear calculation for
                                        ! bulk Richardson number.
  logical, public ::SurfaceStokesInRIb
  integer :: WaveMethod  !< Options for various wave methods
  integer :: SpecMethod  !< Options for various wave spectra
  integer :: NumBands    !< Number of wavenumber bands to recieve
  real ALLOCABLE_, dimension(:) :: WaveNum_Cen !Wavenumber bands for read/coupled
  real ALLOCABLE_, dimension( NIMEMB_, NJMEM_,NKMEM_), public :: &
       Us_x ! Stokes drift (zonal) 
  real ALLOCABLE_, dimension( NIMEM_, NJMEMB_,NKMEM_), public :: &
       Us_y ! Stokes drift (meridional) 
  real ALLOCABLE_, dimension( NIMEM_, NJMEM_) ::&
       LangNum !Langmuir number
  real ALLOCABLE_, dimension( NIMEM_, NJMEM_),public ::    &
       OBLdepth, LangmuirEF_W, LangmuirEF_Vt2, LangmuirEF_K, &
       US0_x, US0_y, US10pct_x, US10pct_y
  real ALLOCABLE_, dimension( NIMEMB_, NJMEM_,NKMEM_), public :: &
       STKx0
  real ALLOCABLE_, dimension( NIMEM_, NJMEMB_,NKMEM_), public :: &
       STKy0  

  logical :: dataoverrideisinitialized
end type

! This include declares and sets the variable "version".
#include "version_variable.h"

character(len=40)  :: mod = "MOM_wave_interface" ! This module's name.

! Switches needed in import_stokes_drift
integer, parameter :: NO_SCHEME = 0, FROMFILE = 1, DATAOVERRIDE =2,&
                      PARAMETRIC = 3, TESTPROF = 99;
integer, parameter :: ELFOUHAILY = 1

! For Test Prof
Real :: TP_STKX0, TP_STKY0, TP_WVL

!-------------------------------------------------------------
CONTAINS
!
!> Initializes parameters related to MOM_wave_interface
subroutine MOM_wave_interface_init(G,GV,param_file, CS)
  type(ocean_grid_type),                  intent(in)  :: G !< Grid structure
  type(verticalGrid_type),                intent(in)  :: GV!< Vertical grid structure
  type(param_file_type),                  intent(in)  :: param_file !< Input parameter structure
  type(wave_parameters_CS),              pointer     :: CS
  ! Local variables

  ! I/O
  character*(13) :: TMPSTRING1,TMPSTRING2
  character*(10), parameter :: NULL_STRING = "EMPTYEMPTY"
  character*(10), parameter :: PARAMETRIC_STRING = "PARAMETRIC" 
  character*(10), parameter :: FROMFILE_STRING = "FROM_FILE"
  character*(13), parameter :: DATAOVERRIDE_STRING = "DATA_OVERRIDE"
  character*(10), parameter :: TESTPROF_STRING = "TEST_PROF"
  character*(10), parameter :: ELF97_STRING = "ELF97"

  if (associated(CS)) then
     call MOM_error(WARNING, "wave_interface_init called with an associated"//&
                             "control structure.")
     return
  endif
  
  allocate(CS)

  ! Add any initializations needed here
  CS%dataOverrideIsInitialized = .false.

  call log_version(param_file, mod, version)
  call get_param(param_file,mod,"USE_WAVES",CS%UseWaves, &
                 'Main switch to use wave input', units='',default=.false.)
  call get_param(param_file, mod, "LAGRANGIAN_MIXING", CS%LagrangianMixing, &
       "Flag to use Lagrangian Mixing", units="", &
       Default=.false.)
  call get_param(param_file, mod, "LANGMUIR_ENHANCE_W", CS%LangmuirEnhanceW, &
       'Flag for Langmuir turbulence enhancement of turbulent'//&
       'velocity scale.', units="", Default=.false.) 
  call get_param(param_file, mod, "LANGMUIR_ENHANCE_VT2", CS%LangmuirEnhanceVt2, &
       'Flag for Langmuir turbulence enhancement of turbulent'//&
       'velocity scale.', units="", Default=.false.) 
  call get_param(param_file, mod, "LANGMUIR_ENHANCE_K", CS%LangmuirEnhanceK, &
       'Flag for Langmuir turbulence enhancement of turbulent'//&
       'velocity scale.', units="", Default=.false.) 
 call get_param(param_file, mod, "STOKES_IN_RIB", CS%StokesShearInRIb, &
       'Flag for using Stokes drift in RIb calculation.'&
       , units="", Default=.false.) 
 call get_param(param_file, mod, "SURFACE_STOKES_IN_RIB", CS%SurfaceStokesInRIb, &
       'Flag for using surface Stokes drift in RIb calculation.'&
       , units="", Default=.false.) 

  if ( (CS%LagrangianMixing.or.CS%LangmuirEnhanceW) .and. (.not.CS%UseWaves)) then
     call MOM_error(FATAL,"MOM_vert_friction(visc): "// &
          "LagrangianMixing and WaveEnhancedDiff cannot"//&
          "be called without USE_WAVES = .true.")
  endif

  if (CS%UseWaves) then 
     ! 1. Get Wave Method and write to integer WaveMethod
     call get_param(param_file,mod,"WAVE_METHOD",TMPSTRING1, &
                    'Choice of wave method, valid options...',units='',&
                    default=NULL_STRING)
     select case (TRIM(TMPSTRING1))
       case (NULL_STRING)! No Waves
          CS%WaveMethod = NO_SCHEME
          Print*,'You did not specify a wave method, so no waves are used.'
       case (FROMFILE_STRING)! From File
          CS%WaveMethod = FROMFILE
       case (DATAOVERRIDE_STRING)! From File
          CS%WaveMethod = DATAOVERRIDE
      case (PARAMETRIC_STRING)! Parameteric Spc
          CS%WaveMethod = PARAMETRIC
          call get_param(param_file,mod,"SPECTRUM_CHOICE",TMPSTRING2, &
                         'Choice of empirical wave spectrum, valid...',&
                         units='',default=NULL_STRING)
          select case(TRIM(TMPSTRING2))
            case (NULL_STRING)
               CS%WaveMethod = NO_SCHEME
               Print*,' '
               Print*,'You did not specify an empirical wave spectrum,'
               Print*,' so no waves are used.'
            case (ELF97_STRING)
               CS%SpecMethod = ELFOUHAILY
          endselect
        case (TESTPROF_STRING)
           CS%WaveMethod = TESTPROF
           call get_param(param_file,mod,"TP_STKX_SURF",TP_STKX0,&
                          'Surface Stokes (x) for test profile',&
                          units='m/s',default=0.1)
           call get_param(param_file,mod,"TP_STKY_SURF",TP_STKY0,&
                          'Surface Stokes (y) for test profile',&
                          units='m/s',default=0.0)
           call get_param(param_file,mod,"TP_WVL",TP_WVL,&
                          units='m',default=50.0)
     endselect
     ! 2. Allocate and initialize
     !    Stokes drift
     ALLOC_ (CS%Us_x(G%isdB:G%IedB,G%jsd:G%jed,G%ke)) ; CS%Us_x(:,:,:) = 0.0
     ALLOC_ (CS%Us_y(G%isd:G%Ied,G%jsdB:G%jedB,G%ke)) ; CS%Us_y(:,:,:) = 0.0
     !    Langmuir number
     ALLOC_ (CS%LangNum(G%isc:G%iec,G%jsc:G%jec)) ; CS%LangNum(:,:) = 1e10
     ALLOC_ (CS%LangmuirEF_W(G%isc:G%iec,G%jsc:G%jec)) ; CS%LangmuirEF_W(:,:) = 1.
     ALLOC_ (CS%LangmuirEF_Vt2(G%isc:G%iec,G%jsc:G%jec)) ; CS%LangmuirEF_Vt2(:,:) = 1.
     ALLOC_ (CS%LangmuirEF_K(G%isc:G%iec,G%jsc:G%jec)) ; CS%LangmuirEF_K(:,:) = 1.
     ALLOC_ (CS%OBLdepth(G%isc:G%iec,G%jsc:G%jec)) ; CS%OBLdepth(:,:) = 0.
     ALLOC_ (CS%US0_x(G%isdB:G%iedB,G%jsd:G%jed)) ; CS%US0_x(:,:) = 0.
     ALLOC_ (CS%US0_y(G%isd:G%ied,G%jsdB:G%jedB)) ; CS%US0_y(:,:) = 0.  
     ALLOC_ (CS%US10pct_x(G%isdB:G%iedB,G%jsd:G%jed)) ; CS%US10pct_x(:,:) = 0.
     ALLOC_ (CS%US10pct_y(G%isd:G%ied,G%jsdB:G%jedB)) ; CS%US10pct_y(:,:) = 0. 
  endif
  !/BGRTEMP{
  print*,' '
  print*,'-----------------------------------------------'
  print*,'You chose this wave method: ',CS%WaveMethod
  if(CS%WaveMethod==PARAMETRIC) then
     print*,'You chose this specrum: ',CS%SpecMethod
  endif
  if (CS%WaveMethod==TESTPROF) then
     print*,' '
     print*,'You chose the following for the test profile'
     print*,'--------------------------------------------'
     print*,'Surface Stk X [m/s]: ',TP_STKX0
     print*,'Surface Stk Y [m/s]: ',TP_STKY0
     print*,'Mean Wavelength [m]: ',TP_WVL
  endif
  print*,'-----------------------------------------------'
  print*,' '
  !\BGRTEMP}

end subroutine MOM_wave_interface_init
!/
!/
!/
! Constructs the Stokes Drift profile on the model grid based on 
! desired coupling options
subroutine Import_Stokes_Drift(G,GV,Day,DT,CS,h,FLUXES)
  type(wave_parameters_CS),              pointer        :: CS
  type(ocean_grid_type),                  intent(in)    :: G !< Grid structure
  type(verticalGrid_type),                intent(in)    :: GV!< Vertical grid structure
  type(time_type), intent(in)                           :: Day
  type(time_type), intent(in)                           :: DT
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in) :: h
  type(forcing), intent(in)                             :: FLUXES
  ! local variables
  real    :: USy20pct, USx20pct, H20pct, H10pct
  real    :: Top, MidPoint, Bottom
  real    :: DecayScale
  type(time_type) :: Day_Center
  integer :: ii, jj, kk, b, iim1, jjm1

  Day_Center = Day + DT/2

  !print*,' '
  !print*,'Into Import_Stokes_Drift'

  if (CS%WaveMethod==TESTPROF) then
     DecayScale = 12.5663706/TP_WVL !4pi
     !print*,'Test Profile Construction'
     Bottom = 0.0
     MidPoint = 0.0
     !print*,'DecayScale:',TP_WVL,TP_STKX0,TP_STKY0
     do kk=1, G%ke
        Top = Bottom
        !****************************************************
        !NOTE THIS H WILL NOT BE CORRECT FOR NON-UNIFORM GRID
        MidPoint = Bottom - GV%H_to_m * h(1,1,kk)/2.
        Bottom = Bottom - GV%H_to_m * h(1,1,kk)
        do ii=G%isdB,G%iedB
           do jj=G%jsd,G%jed
              CS%Us_x(ii,jj,kk) = TP_STKX0 * EXP(MIDPOINT*DecayScale)
           enddo
        enddo
        do ii=G%isd,G%ied
           do jj=G%jsdB,G%jedB
              CS%Us_y(ii,jj,kk) = TP_STKY0 * EXP(MIDPOINT*DecayScale)
           enddo
        enddo
        !print*,MIDPOINT,CS%US_x(1,1,kk),CS%Us_y(1,1,kk)
     enddo     
  elseif (CS%WaveMethod==DATAOVERRIDE) then
       call Stokes_Drift_by_data_override(day_center,G,GV,CS)
       CS%Us_x(:,:,:) = 0.0
       CS%Us_y(:,:,:) = 0.0
       CS%Us0_x(:,:) = 0.0
       CS%Us0_y(:,:) = 0.0
       CS%Us10pct_x(:,:) = 0.0
       CS%Us10pct_y(:,:) = 0.0
     ! ---------------------------------------------------------|
     ! This computes the average Stokes drift based on the      |
     !  analytical integral over the layer divided by the layer |
     !  thickness.                                              |
     ! ---------------------------------------------------------|
     do ii=G%isdB,G%iedB
        do jj=G%jsd,G%jed
           do b=1,CS%NumBands
              CS%US0_x(ii,jj)=CS%US0_x(ii,jj) + CS%STKx0(ii,jj,b) *&
                   (1.0 - EXP(-0.01*2*CS%WaveNum_Cen(b))) / (0.01)/&
                   (2*CS%WaveNum_Cen(b))
           enddo
           bottom = 0.0
           do kk=1, G%ke
              Top = Bottom
              iim1 = max(ii-1,1)
              MidPoint = Bottom - GV%H_to_m * (h(ii,jj,kk)+h(iim1,jj,kk))/4.
              Bottom = Bottom - GV%H_to_m *  (h(ii,jj,kk)+h(iim1,jj,kk))/2.
              do b=1,CS%NumBands
                 CS%US_x(ii,jj,kk)=CS%US_x(ii,jj,kk) + CS%STKx0(ii,jj,b) *&
                      (EXP(TOP*2*CS%WaveNum_Cen(b))- &
                      EXP(BOTTOM*2*CS%WaveNum_Cen(b))) / (Top-Bottom)/&
                       (2*CS%WaveNum_Cen(b))
              enddo
           enddo
        enddo
     enddo
     do ii=G%iscB,G%iecB
        do jj=G%jsc,G%jec
           H10pct=min(-0.1,-CS%OBLdepth(ii,jj)*0.2);
           bottom = 0.0
           do b=1,CS%NumBands
              CS%US10pct_x(ii,jj)=CS%US10pct_x(ii,jj) + &
                   0.5*(CS%STKx0(ii,jj,b)+CS%STKx0(ii-1,jj,b)) *&
                   (1.0 - EXP(H10pct*CS%WaveNum_Cen(b))) /H10pct/&
                   (2*CS%WaveNum_Cen(b))
           enddo
        enddo
     enddo

     do ii=G%isd,G%ied
        do jj=G%jsdB,G%jedB
           do b=1,CS%NumBands
              CS%US0_y(ii,jj)=CS%US0_y(ii,jj) + CS%STKy0(ii,jj,b) *&
                   (1.0 - EXP(-0.01*2*CS%WaveNum_Cen(b))) / (0.01)/&
                   (2*CS%WaveNum_Cen(b))
           enddo
           bottom = 0.0
           do kk=1, G%ke
              Top = Bottom
              jjm1 = max(jj-1,1)
              MidPoint = Bottom - GV%H_to_m * (h(ii,jj,kk)+h(ii,jjm1,kk))/4.
              Bottom = Bottom - GV%H_to_m *  (h(ii,jj,kk)+h(ii,jjm1,kk))/2.
              do b=1,CS%NumBands
                 CS%US_y(ii,jj,kk)=CS%US_y(ii,jj,kk) + CS%STKy0(ii,jj,b) *&
                      (EXP(TOP*2*CS%WaveNum_Cen(b))- &
                      EXP(BOTTOM*2*CS%WaveNum_Cen(b))) / (Top-Bottom)/&
                       (2*CS%WaveNum_Cen(b))
              enddo
           enddo
        enddo
     enddo
     do ii=G%isc,G%iec
        do jj=G%jscB,G%jecB
           H10pct=min(-0.1,-CS%OBLdepth(ii,jj)*0.2);
           bottom = 0.0
           do b=1,CS%NumBands
              CS%US10pct_y(ii,jj)=CS%US10pct_y(ii,jj) + &
                   0.5*(CS%STKy0(ii,jj,b)+CS%STKy0(ii,jj-1,b)) *&
                   (1.0 - EXP(H10pct*CS%WaveNum_Cen(b))) /H10pct/&
                   (2*CS%WaveNum_Cen(b))
           enddo
        enddo
     enddo


     !At h points for Langmuir number
     do ii=G%isc,G%iec
        do jj=G%jsc,G%jec
           USy20pct = 0.0;USx20pct = 0.0; 
           H20pct=min(-0.1,-CS%OBLdepth(ii,jj)*0.2);
           do b=1,CS%NumBands
              USy20pct=USy20pct + &
                   0.5*(CS%STKy0(ii,jj,b)+CS%STKy0(ii,jj-1,b)) *&
                   (1.0 - EXP(H20pct*2*CS%WaveNum_Cen(b))) &
                   / (0.0-H20pct) / (2*CS%WaveNum_Cen(b))
              USx20pct=USx20pct + &
                   0.5*(CS%STKx0(ii,jj,b)+CS%STKx0(ii-1,jj,b)) *&
                   (1.0 - EXP(H20pct*2*CS%WaveNum_Cen(b))) &
                   / (0.0-H20pct) / (2*CS%WaveNum_Cen(b))
           enddo
           CS%LangNum(ii,jj) = sqrt(FLUXES%ustar(ii,jj) / &
                max(1.e-10,sqrt(USx20pct**2 + USy20pct**2)))
           if (CS%LangmuirEnhanceW) then
              !McWilliams et al., 2000
              !CS%LangmuirEF_W(ii,jj) = sqrt(1+0.08/CS%LangNum(ii,jj)**4)
              !VanRoekel et a. 2012
              CS%LangmuirEF_W(ii,jj) = sqrt(1.+(1.5*CS%LangNum(ii,jj))**(-2) + &
                                        (5.4*CS%LangNum(ii,jj))**(-4))
           endif
           if (CS%LangmuirEnhanceVt2) then
             CS%LangmuirEF_Vt2(ii,jj) = min(50.,1. + 2.3/sqrt(CS%LangNum(ii,jj)))
           endif
           if (CS%LangmuirEnhanceK) then
             CS%LangmuirEF_K(ii,jj) = min(2.25, 1. + 1./CS%LangNum(ii,jj))
           endif
        enddo
     enddo
  else!Keep this else, fallback to 0 Stokes drift
     do ii=G%isdB,G%iedB
           do jj=G%jsd,G%jed
              CS%Us_x(ii,jj,kk) = 0
           enddo
        enddo
        do ii=G%isd,G%ied
           do jj=G%jsdB,G%jedB
              CS%Us_y(ii,jj,kk) = 0
           enddo
        enddo
  endif

end subroutine Import_Stokes_Drift
!
subroutine Stokes_Drift_by_data_override(day_center,G,GV,CS)
  use NETCDF
  type(time_type),             intent(in)  :: day_center
  type(wave_parameters_CS),    pointer     :: CS
  type(ocean_grid_type),       intent(in)  :: G !< Grid structure
  type(verticalGrid_type),     intent(in)  :: GV!< Vertical grid structure
  ! local variables
  real    :: Top, MidPoint, Bottom
  real    :: DecayScale
  integer :: b
  integer :: i, j

  integer, dimension(4) :: start, count, dims, dim_id 
  character(len=12)  :: dim_name(4)
  character(20) :: varname, filename, varread
  integer :: rcode, ncid, varid, id, ndims

  if (.not.CS%dataOverrideIsInitialized) then
    print*,'into init'
    call data_override_init(Ocean_domain_in=G%Domain%mpp_domain)
    print*,'out of init'
    CS%dataOverrideIsInitialized = .true.
    
    ! Read in number of wavenumber bands in file to set number to be read in
    filename = 'StkSpec.nc'
    varread = 'wavenumber'

    rcode = NF90_OPEN(trim(filename), NF90_NOWRITE, ncid)
    if (rcode .ne. 0) call MOM_error(FATAL,"error opening file "//trim(filename)//&
         " in MOM_wave_interface.")

    rcode = NF90_INQ_VARID(ncid, varread, varid)
    if (rcode .ne. 0) call MOM_error(FATAL,"error finding variable "//trim(varread)//&
         " in file "//trim(filename)//" in MOM_wave_interface.")

    rcode = NF90_INQUIRE_VARIABLE(ncid, varid, ndims=ndims, dimids=dims)
    if (rcode .ne. 0) call MOM_error(FATAL,'error inquiring dimensions MOM_wave_interface.')

    rcode = NF90_INQUIRE_DIMENSION(ncid, dims(1), dim_name(1), len=id)
    if (rcode .ne. 0) call MOM_error(FATAL,"error reading dimension 1 data for "// &
         trim(varread)//" in file "// trim(filename)//" in MOM_wave_interface.")

    rcode = NF90_INQ_VARID(ncid, dim_name(1), dim_id(1))
    if (rcode .ne. 0) call MOM_error(FATAL,"error finding variable "//trim(dim_name(1))//&
         " in file "//trim(filename)//" in MOM_wave_interace.")

    ! Allocating size of wavenumber bins
    ALLOC_ ( CS%WaveNum_Cen(1:id) ) ; CS%WaveNum_Cen(:)=0.0
    ALLOC_ ( CS%STKx0(G%isdB:G%iedB,G%jsd:G%jed,1:id)) ; CS%STKx0(:,:,:) = 0.0
    ALLOC_ ( CS%STKy0(G%isd:G%ied,G%jsdB:G%jedB,1:id)) ; CS%STKy0(:,:,:) = 0.0

    ! Reading wavenumber bins
    start = 1; count = 1; count(1) = id
    rcode = NF90_GET_VAR(ncid, dim_id(1), CS%WaveNum_Cen, start, count)
    if (rcode .ne. 0) call MOM_error(FATAL,"error reading dimension 1 values for var_name "// &
         trim(varread)//",dim_name "//trim(dim_name(1))//" in file "// trim(filename)//" in MOM_wave_interface")

    CS%NUMBANDS = ID
    print*,CS%WaveNum_Cen

    print*,'******************'
    print*,'End of NetCDF Read'
  endif
  
  do b=1,CS%NumBands
     varname = '                    '
     write(varname,"(A3,I0)")'Usx',b
     call data_override('OCN',trim(varname), CS%STKx0(:,:,b), day_center)
     varname = '                    '
     write(varname,'(A3,I0)')'Usy',b
     call data_override('OCN',trim(varname), CS%STKy0(:,:,b), day_center)     
  enddo
  
  !Brandon: Hacking to update HALO until properly resolved
  do i=G%isdB,G%iedB
     do j=G%isd,G%ied
        if (i.lt.G%iscB) then
           CS%STKX0(i,j,:)=CS%STKX0(G%iscB,j,:)
        elseif (i.gt.G%iecB) then
           CS%STKX0(i,j,:)=CS%STKX0(G%iecB,j,:)
        endif
     enddo
  enddo
  do i=G%isd,G%ied
     do j=G%isdB,G%iedB
        if (j.lt.G%jscB) then
           CS%STKY0(i,j,:)=CS%STKY0(i,G%iscB,:)
        elseif (j.gt.G%jecB) then
           CS%STKY0(i,j,:)=CS%STKY0(i,G%jecB,:)
        endif
     enddo
  enddo
  
  !print*,' '
  !print*,'-------------------------------------'
  !print*,'End of Stokes Drift By Data Override.'
  !print*,'-------------------------------------'
  !print*,' '


end subroutine Stokes_Drift_by_Data_Override
end module MOM_wave_interface
