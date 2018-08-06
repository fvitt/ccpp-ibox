!> \file ibox_type_defs.f90
!!  Contains type definitions for ibox variables and physics-related variables

module ibox_vardefs

use machine, only: kind_phys

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! The following definition sets up the variables for use within ibox
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

! Filter with CPP for PGI compiler
#ifndef __PGI
!> \section arg_table_ibox_vardefs
!! | local_name             | standard_name               | long_name               | units       | rank | type      |    kind   | intent | optional |
!! |------------------------|-----------------------------|-------------------------|-------------|------|-----------|-----------|--------|----------|
!! | my_co(:)               | my_volume_mixing_ratio_co   | volume mixing ratio co  | mole mole-1 |    1 | real      | kind_phys | none   | F        |
!! | my_o3(:)               | my_volume_mixing_ratio_o3   | volume mixing ratio o3  | mole mole-1 |    1 | real      | kind_phys | none   | F        |
!! | state_host%Temperature | air_temperature             | temperature             | K           |    0 | real      | kind_phys | none   | F        | 
!! | errmsg                 | error_message               | CCPP error message      | none        |    0 | character | len=*     | out    | F        | 
!! | errflg                 | error_flag                  | CCPP error flag         | flag        |    0 | integer   |           | out    | F        |
!! | tune_factor            | tuning factor for solver    | tuning factor for solver| mole mole-1 |    0 | real      | kind_phys | none   | F        |
!! | k_rateConst            | k_rate_constants            | k Rate Constants        | none        |    1 | real      | kind_phys | none   | F        |
!!

#endif

  real (kind=kind_phys), pointer :: my_co(:)
  real (kind=kind_phys), pointer :: my_o3(:)
  character(len=512)     :: errmsg
  integer :: errflg

  type my_state
    real(kind_phys) :: Temperature
  end type my_state
  type(my_state), target :: state_host

  real(kind_phys)  :: tune_factor
  real, pointer    :: k_rateConst(:)

end module ibox_vardefs
