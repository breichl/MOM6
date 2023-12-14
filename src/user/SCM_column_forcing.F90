!> Adds routines to provide external columnwise forcing, useful for SCM experiments.
module SCM_column_forcing

! This file is part of MOM6. See LICENSE.md for the license.

! History
!--------
! December 2022: Origination.
!

use time_interp_external_mod, only: init_external_field, time_interp_external
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_alloc
use MOM_diag_mediator, only : diag_ctrl
use MOM_error_handler, only : MOM_error, FATAL
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_io,            only : file_exists, get_var_sizes, read_variable
use MOM_safe_alloc, only : safe_alloc_ptr
use MOM_time_manager, only : time_type, operator(+), operator(/), real_to_time
use MOM_unit_scaling,  only : unit_scale_type
use MOM_verticalGrid, only : verticalGrid_type
use MOM_variables,             only : thermo_var_ptrs
use MOM_domains,              only : pass_var, pass_vector

implicit none ; private

#include <MOM_memory.h>

public SCM_column_forcing_init
public SCM_column_forcing_calculate
public SCM_column_forcing_apply_thermo
public SCM_column_forcing_apply_dynamics

!> Container for parameters describing idealized wind structure
type, public :: SCM_column_forcing_CS ; private

  real, allocatable, dimension(:,:,:) :: dT_dt_input, dT_dt
  real, allocatable, dimension(:,:,:) :: dS_dt_input, dS_dt
  real, allocatable, dimension(:,:,:) :: dU_dt_input, dU_dt
  real, allocatable, dimension(:,:,:) :: dV_dt_input, dV_dt

  character(len=40) :: z_name_temp
  character(len=40) :: varname_temp
  real, allocatable, dimension(:) :: z_input_temp
  integer :: NZ_input_temp

  character(len=40) :: z_name_salt
  character(len=40) :: varname_salt
  real, allocatable, dimension(:) :: z_input_salt
  integer :: NZ_input_salt

  character(len=40) :: z_name_ucur
  character(len=40) :: varname_ucur
  real, allocatable, dimension(:) :: z_input_ucur
  integer :: NZ_input_ucur

  character(len=40) :: z_name_vcur
  character(len=40) :: varname_vcur
  real, allocatable, dimension(:) :: z_input_vcur
  integer :: NZ_input_vcur

  logical :: apply_tendency_temp
  logical :: apply_tendency_salt
  logical :: apply_tendency_ucur
  logical :: apply_tendency_vcur

  character(len=40)  :: filename_temp !< Temp filename if using data_override
  integer            :: file_ID_temp
  character(len=40)  :: filename_salt !< Salt filename if using data_override
  integer            :: file_ID_salt
  character(len=40)  :: filename_ucur !< Current filename if using data_override
  integer            :: file_ID_ucur
  character(len=40)  :: filename_vcur !< Current filename if using data_override
  integer            :: file_ID_vcur

  logical :: override_init !< Flag for data override initialization

  type(time_type), pointer :: Time !< A pointer to the ocean model's clock.
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                                   !! timing of diagnostic output.

  !>@{ Diagnostic handles
  integer :: id_dU_dt = -1, id_dV_dt = -1, id_dT_dt = -1, id_dS_dt = -1
  !>@}
  end type SCM_column_forcing_CS

! This include declares and sets the variable "version".
#include "version_variable.h"

character(len=40)  :: mdl = "SCM_column_forcing" !< This module's name.

contains

!> Initializes SCM_forcing parameters
subroutine SCM_column_forcing_init(Time, G, GV, US, param_file, CS, diag )
  type(time_type), target,    intent(in) :: Time   !< Model time
  type(ocean_grid_type),      intent(in) :: G      !< Grid structure
  type(verticalGrid_type), intent(in)    :: GV         !< Vertical grid structure
  type(unit_scale_type),      intent(in) :: US     !< A dimensional unit scaling type
  type(param_file_type),      intent(in) :: param_file !< Input parameter structure
  type(SCM_column_forcing_CS),  pointer  :: CS     !< Parameter container for this module
  type(diag_ctrl), target, intent(inout) :: diag       !< Diagnostic Pointer
  integer :: ndims
  integer, dimension(4) :: sizes
  character(len=48) :: dim_name(4)

  ! This include declares and sets the variable "version".
# include "version_variable.h"

  if (associated(CS)) then
    call MOM_error(FATAL, "SCM_forcing_init called with an associated control structure.")
    return
  endif

  allocate(CS)

  CS%diag => diag
  CS%Time => Time

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  ! Checking for external tendencies and files.

  ! 1. Temperature
  call get_param(param_file, mdl, "EXTERNAL_SCM_TEMP_TENDENCY", CS%apply_tendency_temp, &
                 "Logical to toggle the external SCM model temperature tendency term.",&
                 units='nondim', default=.false.)
  if (CS%apply_tendency_temp) then
    call get_param(param_file, mdl, "EXTERNAL_SCM_TEMP_FILENAME", CS%filename_temp, &
                 "Filename for the external SCM model temperature tendency term.",&
                 default='')
    call get_param(param_file, mdl, "SCM_TEMP_DEPTH_COORDINATE", CS%z_name_temp, &
                 "Variable name for the external SCM model temperature tendency term vertical coordinate.",&
                 default='depth')
    ! Check for file
    if (.not.file_exists(CS%filename_temp)) &
      call MOM_error(FATAL, "SCM_column_forcing is unable to find file "//trim(CS%filename_temp))
    call get_param(param_file, mdl, "EXTERNAL_SCM_TEMP_VARNAME", CS%varname_temp, &
                 "Variable name for the external SCM model temperature tendency term.",&
                 default='dTdt')
    ! Check for variable
    call get_var_sizes(CS%filename_temp, trim(CS%varname_temp), &
         ndims, sizes, dim_names=dim_name)
    if (ndims < 0) &
      call MOM_error(FATAL,"No variable "//trim(CS%varname_temp)//" in "//&
                           trim(CS%filename_temp))
    call get_param(param_file, mdl, "EXTERNAL_SCM_TEMP_ZNAME", CS%z_name_temp, &
                 "Depth coordinate name for the external SCM model temperature tendency term.",&
                 default='depth')
    ! Check for depth coordinate
    call get_var_sizes(CS%filename_temp, trim(CS%z_name_temp), ndims, sizes, dim_names=dim_name)
    if (ndims < 0) &
      call MOM_error(FATAL,"No depth variable "//trim(CS%z_name_temp)//" in "//&
                           trim(CS%filename_temp))
    CS%NZ_input_temp = sizes(1)
    allocate( CS%z_input_temp(CS%NZ_input_temp), source=0.0 )
    call read_variable(CS%filename_temp, dim_name(1), CS%z_input_temp, scale=US%Z_to_m)
    allocate( CS%dT_dt_input(G%isc:G%iec,G%jsc:G%jec,CS%NZ_input_temp), source=0.0 )
    allocate( CS%dT_dt(G%isd:G%ied,G%jsd:G%jed,GV%ke), source=0.0 )

    CS%file_ID_temp = init_external_field(CS%filename_temp, CS%varname_temp)

  endif
  ! 2. Salinity
  call get_param(param_file, mdl, "EXTERNAL_SCM_SALT_TENDENCY", CS%apply_tendency_salt, &
                 "Logical to toggle the external SCM model salinity tendency term.",&
                 units='nondim', default=.false.)
  if (CS%apply_tendency_salt) then
    call get_param(param_file, mdl, "EXTERNAL_SCM_SALT_FILENAME", CS%filename_salt, &
                 "Filename for the external SCM model salinity tendency term.",&
                 default='')
    ! Check for file
    if (.not.file_exists(CS%filename_salt)) &
      call MOM_error(FATAL, "SCM_column_forcing is unable to find file "//trim(CS%filename_salt))
    call get_param(param_file, mdl, "EXTERNAL_SCM_SALT_VARNAME", CS%varname_salt, &
                 "Variable name for the external SCM model salinity tendency term.",&
                 default='dSdt')
    ! Check for variable
    call get_var_sizes(CS%filename_salt, trim(CS%varname_salt), &
         ndims, sizes, dim_names=dim_name)
    if (ndims < 0) &
      call MOM_error(FATAL,"No variable "//trim(CS%varname_salt)//" in "//&
                           trim(CS%filename_salt))
    call get_param(param_file, mdl, "EXTERNAL_SCM_SALT_ZNAME", CS%z_name_salt, &
                 "Depth coordinate name for the external SCM model salinity tendency term.",&
                 default='depth')
    ! Check for depth coordinate
    call get_var_sizes(CS%filename_salt, trim(CS%z_name_salt), ndims, sizes, dim_names=dim_name)
    if (ndims < 0) &
     call MOM_error(FATAL,"No depth variable "//trim(CS%z_name_salt)//" in "//&
                           trim(CS%filename_salt))
    CS%NZ_input_salt = sizes(1)
    allocate( CS%z_input_salt(CS%NZ_input_salt), source=0.0 )
    call read_variable(CS%filename_salt, dim_name(1), CS%z_input_salt, scale=US%Z_to_m)
    allocate( CS%dS_dt_input(G%isc:G%iec,G%jsc:G%jec,CS%NZ_input_salt), source=0.0 )
    allocate( CS%dS_dt(G%isd:G%ied,G%jsd:G%jed,GV%ke), source=0.0 )

    CS%file_ID_salt = init_external_field(CS%filename_salt, CS%varname_salt)

  endif

  ! 3. currents (done separately)
  call get_param(param_file, mdl, "EXTERNAL_SCM_XCURRENT_TENDENCY", CS%apply_tendency_ucur, &
                 "Logical to toggle the external SCM model u-current tendency term.",&
                 units='nondim', default=.false.)
  if (CS%apply_tendency_ucur) then
    call get_param(param_file, mdl, "EXTERNAL_SCM_XCURRENT_FILENAME", CS%filename_ucur, &
                 "Filename for the external SCM model x-current tendency term.",&
                 default='')
    ! Check for file
    if (.not.file_exists(CS%filename_ucur)) &
      call MOM_error(FATAL, "SCM_column_forcing is unable to find file "//trim(CS%filename_ucur))
    call get_param(param_file, mdl, "EXTERNAL_SCM_UCURRENT_VARNAME", CS%varname_ucur, &
                 "Variable name for the external SCM model x-current tendency term.",&
                 default='dUdt')
    ! Check for variable
    call get_var_sizes(CS%filename_ucur, trim(CS%varname_ucur), &
         ndims, sizes, dim_names=dim_name)
    if (ndims < 0) &
      call MOM_error(FATAL,"No variable "//trim(CS%varname_ucur)//" in "//&
                           trim(CS%filename_ucur))
    call get_param(param_file, mdl, "EXTERNAL_SCM_UCURRENT_ZNAME", CS%z_name_ucur, &
                 "Depth coordinate name for the external SCM model x-current tendency term.",&
                 default='depth')
    ! Check for depth coordinate
    call get_var_sizes(CS%filename_ucur, trim(CS%z_name_ucur), ndims, sizes, dim_names=dim_name)
    if (ndims < 0) &
      call MOM_error(FATAL,"No depth variable "//trim(CS%z_name_ucur)//" in "//&
                           trim(CS%filename_ucur))
    CS%NZ_input_ucur = sizes(1)
    allocate( CS%z_input_ucur(CS%NZ_input_ucur), source=0.0 )
    call read_variable(CS%filename_ucur, dim_name(1), CS%z_input_ucur, scale=US%Z_to_m)
    allocate( CS%dU_dt_input(G%isc:G%iec,G%jsc:G%jec,CS%NZ_input_ucur), source=0.0 )
    allocate( CS%dU_dt(G%isdB:G%iedB,G%jsd:G%jed,GV%ke), source=0.0 )

    CS%file_ID_ucur = init_external_field(CS%filename_ucur, CS%varname_ucur)

  endif

  call get_param(param_file, mdl, "EXTERNAL_SCM_YCURRENT_TENDENCY", CS%apply_tendency_vcur, &
                 "Logical to toggle the external SCM model v-current tendency term.",&
                 units='nondim', default=.false.)
  if (CS%apply_tendency_vcur) then
call get_param(param_file, mdl, "EXTERNAL_SCM_YCURRENT_FILENAME", CS%filename_vcur, &
                 "Filename for the external SCM model y-current tendency term.",&
                 default='dVdt')
    ! Check for file
    if (.not.file_exists(CS%filename_vcur)) &
      call MOM_error(FATAL, "SCM_column_forcing is unable to find file "//trim(CS%filename_vcur))
    call get_param(param_file, mdl, "EXTERNAL_SCM_VCURRENT_VARNAME", CS%varname_vcur, &
                 "Variable name for the external SCM model y-current tendency term.",&
                 default='dVdt')
    ! Check for variable
    call get_var_sizes(CS%filename_vcur, trim(CS%varname_vcur), &
         ndims, sizes, dim_names=dim_name)
    if (ndims < 0) &
      call MOM_error(FATAL,"No variable "//trim(CS%varname_vcur)//" in "//&
                           trim(CS%filename_vcur))
    call get_param(param_file, mdl, "EXTERNAL_SCM_VCURRENT_ZNAME", CS%z_name_vcur, &
                 "Depth coordinate name for the external SCM model y-current tendency term.",&
                 default='depth')
    ! Check for depth coordinate
    call get_var_sizes(CS%filename_vcur, trim(CS%z_name_vcur), ndims, sizes, dim_names=dim_name)
    if (ndims < 0) &
      call MOM_error(FATAL,"No depth variable "//trim(CS%z_name_vcur)//" in "//&
                           trim(CS%filename_vcur))
    CS%NZ_input_vcur = sizes(1)
    allocate( CS%z_input_vcur(CS%NZ_input_vcur), source=0.0 )
    call read_variable(CS%filename_vcur, dim_name(1), CS%z_input_vcur, scale=US%Z_to_m)
    allocate( CS%dV_dt_input(G%isc:G%iec,G%jsc:G%jec,CS%NZ_input_vcur), source=0.0 )
    allocate( CS%dV_dt(G%isd:G%ied,G%jsdB:G%jedB,GV%ke), source=0.0 )

    CS%file_ID_vcur = init_external_field(CS%filename_vcur, CS%varname_vcur)

  endif

  if (CS%apply_tendency_ucur) then
    CS%id_dU_dt = register_diag_field('ocean_model','dU_dt_external', &
         CS%diag%axesCuL,Time,'dU_dt from external file', 'm s-2', conversion=US%L_T2_to_m_s2)
  endif
  if (CS%apply_tendency_vcur) then
    CS%id_dV_dt = register_diag_field('ocean_model','dV_dt_external', &
         CS%diag%axesCvL,Time,'dV_dt from external file', 'm s-2', conversion=US%L_T2_to_m_s2)
  endif
  if (CS%apply_tendency_temp) then
    CS%id_dT_dt = register_diag_field('ocean_model','dT_dt_external', &
         CS%diag%axesTL,Time,'dT_dt from external file', 'm s-2')
  endif
  if (CS%apply_tendency_salt) then
    CS%id_dS_dt = register_diag_field('ocean_model','dS_dt_external', &
         CS%diag%axesTL,Time,'dS_dt from external file', 'm s-2')
  endif

end subroutine SCM_column_forcing_init

!> Initializes SCM_forcing parameters
subroutine SCM_column_forcing_calculate(Time, G, GV, US, CS, h, dt)
 type(time_type),               intent(in) :: Time   !< Model time
  type(ocean_grid_type),         intent(in) :: G      !< Grid structure
  type(verticalGrid_type), intent(in)    :: GV         !< Vertical grid structure
  type(unit_scale_type),         intent(in) :: US     !< A dimensional unit scaling type
  type(SCM_column_forcing_CS),  pointer    :: CS     !< Parameter container for this module
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
       intent(in)    :: h     !< Thickness [H ~> m or kg m-2]
  real,         intent(in)  :: dt  !< Timestep as a time-type

  type(time_type) :: Time_center
  real :: zu, zl, zc, z_input, zl_input, zu_input
  integer :: i, j, k, k_input
  integer, dimension(GV%ke) :: kup, klo
  real, dimension(GV%ke) :: kup_weight, klo_weight

  Time_center = Time + real_to_time(US%T_to_s*dt/2.)


  ! Read in fields, interpolated in time but on the native file depth.

  if (CS%apply_tendency_temp) call time_interp_external(CS%file_ID_temp,Time_center,CS%dT_dt_input)
  if (CS%apply_tendency_salt) call time_interp_external(CS%file_ID_salt,Time_center,CS%dS_dt_input)
  if (CS%apply_tendency_ucur) call time_interp_external(CS%file_ID_ucur,Time_center,CS%dU_dt_input)
  if (CS%apply_tendency_vcur) call time_interp_external(CS%file_ID_vcur,Time_center,CS%dV_dt_input)

  ! Map tendencies to model grid
    ! Uses a simple interpolation (could be done fancier?)
  ! Need dimensional rescaling.

  ! Temperature
  if (CS%apply_tendency_temp) then
    do j=G%jsc,G%jec
      do i=G%isc,G%iec
        zu = 0.0
        k_input = 1
        z_input = abs(CS%z_input_temp(1))
        do k=1,GV%ke
          zc = zu + 0.5*h(i,j,k)
          zl = zu + h(i,j,k)
          ! search through input depth grid until it is deeper than zc
          do while(z_input<zc .and. k_input<CS%NZ_input_temp)
            k_input = k_input+1
            z_input = abs(CS%z_input_temp(k_input))
          enddo
          if (k_input>1) then
            kup(k) = k_input-1
            zu_input = abs(CS%z_input_temp(k_input-1))
            klo(k) = k_input
            zl_input = abs(CS%z_input_temp(k_input))
            klo_weight(k) = (zc-zu_input)/(zl_input-zu_input)
            kup_weight(k) = (zl_input-zc)/(zl_input-zu_input)
          else
             kup(k) = 1 !This could be anything, but this is valid.
            kup_weight(k) = 0.0
            klo(k) = 1
            klo_weight(k) = 1.0
          endif
          zu = zl
        enddo
        do k=1,GV%ke
          CS%dT_dt(i,j,k) = CS%dT_dt_input(i,j,kup(k))*kup_weight(k) + &
                            CS%dT_dt_input(i,j,klo(k))*klo_weight(k)
        enddo
      enddo
    enddo
  endif

  ! Salt
  if (CS%apply_tendency_salt) then
    do j=G%jsc,G%jec
      do i=G%isc,G%iec
        zu = 0.0
        k_input = 1
        z_input = abs(CS%z_input_salt(1))
        do k=1,GV%ke
          zc = zu + 0.5*h(i,j,k)
          zl = zu + h(i,j,k)
          ! search through input depth grid until it is deeper than zc
          do while(abs(z_input)<zc .and. k_input<CS%NZ_input_salt)
            k_input = k_input+1
            z_input = abs(CS%z_input_salt(k_input))
          enddo
          if (k_input>1) then
            kup(k) = k_input-1
            zu_input = abs(CS%z_input_salt(k_input-1))
            klo(k) = k_input
            zl_input = abs(CS%z_input_salt(k_input))
            kup_weight(k) = (zc-zu_input)/(zl_input-zu_input)
            klo_weight(k) = (zl_input-zc)/(zl_input-zu_input)
          else
            kup(k) = k_input !This could be anything, but this is valid.
            kup_weight(k) = 0.0
            klo(k) = k_input
            klo_weight(k) = 1.0
          endif
          zu = zl
        enddo
        do k=1,GV%ke
          CS%dS_dt(i,j,k) = CS%dS_dt_input(G%isc,G%jsc,kup(k))*kup_weight(k) + &
                            CS%dS_dt_input(G%isc,G%jsc,klo(k))*klo_weight(k)
        enddo
      enddo
    enddo
  endif

  ! U-current
  if (CS%apply_tendency_ucur) then
    do j=G%jsc,G%jec
      do I=G%iscB,G%iecB
        zu = 0.0
        k_input = 1
        z_input = abs(CS%z_input_ucur(1))
        do k=1,GV%ke
          zc = zu + 0.25*(h(I,j,k)+h(I+1,j,k))
          zl = zu + 0.5*(h(I,j,k)+h(I+1,j,k))
          ! search through input depth grid until it is deeper than zc
          do while(abs(z_input)<zc .and. k_input<CS%NZ_input_ucur)
            k_input = k_input+1
            z_input = abs(CS%z_input_ucur(k_input))
          enddo
          if (k_input>1) then
            kup(k) = k_input-1
            zu_input = abs(CS%z_input_ucur(k_input-1))
            klo(k) = k_input
            zl_input = abs(CS%z_input_ucur(k_input))
            kup_weight(k) = (zc-zu_input)/(zl_input-zu_input)
            klo_weight(k) = (zl_input-zc)/(zl_input-zu_input)
          else
            kup(k) = k_input !This could be anything, but this is valid.
            kup_weight(k) = 0.0
            klo(k) = k_input
            klo_weight(k) = 1.0
          endif
          zu = zl
        enddo
        do k=1,GV%ke
          CS%dU_dt(I,j,k) = CS%dU_dt_input(G%isc,G%jsc,kup(k))*kup_weight(k) + &
                            CS%dU_dt_input(G%isc,G%jsc,klo(k))*klo_weight(k)
        enddo
      enddo
    enddo
  endif

  ! V-current
  if (CS%apply_tendency_vcur) then
    do J=G%jscB,G%jecB
      do i=G%isc,G%iec
        zu = 0.0
        k_input = 1
        z_input = abs(CS%z_input_vcur(1))
        do k=1,GV%ke
          zc = zu + 0.25*(h(i,J,k)+h(i,J+1,k))
          zl = zu + 0.5*(h(i,J,k)+h(i,J+1,k))
          ! search through input depth grid until it is deeper than zc
          do while(abs(z_input)<zc .and. k_input<CS%NZ_input_vcur)
            k_input = k_input+1
            z_input = abs(CS%z_input_vcur(k_input))
          enddo
          if (k_input>1) then
             kup(k) = k_input-1
             zu_input = abs(CS%z_input_vcur(k_input-1))
            klo(k) = k_input
            zl_input = abs(CS%z_input_vcur(k_input))
            kup_weight(k) = (zc-zu_input)/(zl_input-zu_input)
            klo_weight(k) = (zl_input-zc)/(zl_input-zu_input)
          else
            kup(k) = k_input !This could be anything, but this is valid.
            kup_weight(k) = 0.0
            klo(k) = k_input
            klo_weight(k) = 1.0
          endif
          zu = zl
        enddo
        do k=1,GV%ke
          CS%dV_dt(i,J,k) = CS%dV_dt_input(G%isc,G%jsc,kup(k))*kup_weight(k) + &
                            CS%dV_dt_input(G%isc,G%jsc,klo(k))*klo_weight(k)
        enddo
      enddo
    enddo
  endif

  !Write out some diagnostics
  if (CS%id_dU_dt>0 .and. CS%apply_tendency_ucur) &
    call post_data(CS%id_dU_dt, CS%dU_dt, CS%diag)
  if (CS%id_dV_dt>0 .and. CS%apply_tendency_vcur) &
    call post_data(CS%id_dV_dt, CS%dV_dt, CS%diag)
  if (CS%id_dT_dt>0 .and. CS%apply_tendency_temp) &
    call post_data(CS%id_dT_dt, CS%dT_dt, CS%diag)
  if (CS%id_dS_dt>0 .and. CS%apply_tendency_salt) &
    call post_data(CS%id_dS_dt, CS%dS_dt, CS%diag)


end subroutine SCM_column_forcing_calculate

!> This subroutine applys a column "forcing" time tendency on temperature and salinity, if activated.
subroutine SCM_column_forcing_apply_thermo(G, GV, tv, CS, dt)
  type(ocean_grid_type),         intent(in) :: G      !< Grid structure
  type(verticalGrid_type), intent(in)    :: GV         !< Vertical grid structure
  type(thermo_var_ptrs),    intent(inout) :: tv     !< A structure pointing to various thermodynamic variables
  type(SCM_column_forcing_CS),  pointer    :: CS     !< Parameter container for this module
  real,                     intent(in)    :: dt      !< The time interval over which to advance [T ~> s]

  integer :: i, j, k
  real :: inc, maxinc

  if (CS%apply_tendency_temp) then
    maxinc = 0.0
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      do k=1,GV%ke
        inc = CS%dT_dt(i,j,k)*dt
        maxinc = max(abs(inc),maxinc)
        tv%T(i,j,k) = tv%T(i,j,k) + inc
      enddo
    enddo ; enddo
  endif

  if (CS%apply_tendency_salt) then
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      tv%S(i,j,:) = tv%S(i,j,:) + CS%dS_dt(i,j,:)*dt
    enddo; enddo
  endif

  call pass_var(tv%S,G%domain)
  call pass_var(tv%T,G%domain)

end subroutine SCM_column_forcing_apply_thermo

!> This subroutine applys a column "forcing" acceleration on u and v currents, if activated.
subroutine SCM_column_forcing_apply_dynamics(G, GV, u, v, CS, dt)
  type(ocean_grid_type),         intent(in) :: G      !< Grid structure
  type(verticalGrid_type), intent(in)    :: GV         !< Vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                            intent(inout) :: u      !< zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                            intent(inout) :: v      !< meridional velocity [L T-1 ~> m s-1]
  type(SCM_column_forcing_CS),  pointer    :: CS     !< Parameter container for this module
  real,                     intent(in)    :: dt      !< The time interval over which to advance [T ~> s]

  real :: maxinc, inc
  integer :: i, j, k


  if (CS%apply_tendency_ucur) then
    maxinc=0.0
    do j=G%jsc,G%jec ; do i=G%iscB,G%iecB
      do k=1,GV%ke
        inc = CS%dU_dt(I,j,k)*dt
        maxinc = max(abs(inc),maxinc)
        u(I,j,k) = u(I,j,k) + inc
      enddo
    enddo; enddo
  endif

  if (CS%apply_tendency_vcur) then
    maxinc=0.0
    do j=G%jscB,G%jecB ; do i=G%isc,G%iec
      do k=1,GV%ke
        inc = CS%dV_dt(i,J,k)*dt
        maxinc = max(abs(inc),maxinc)
        v(i,J,k) = v(i,J,k) + inc
      enddo
    enddo; enddo
  endif

  call pass_vector(u, v, G%Domain)

end subroutine SCM_column_forcing_apply_dynamics

end module SCM_column_forcing
