module ibox_main

use machine, only: kind_phys
use input_file, only: input_file_type
use input_slice, only : slice_type
use json_loader, only : create_cnst_info_array
use const_props_mod, only: const_props_type

implicit none

contains

subroutine ibox_main_sub()

  ! Add the CCPP specific types and functions
  use :: ccpp_api,                           &
         only: ccpp_t,                       &
               ccpp_init,                    &
               ccpp_finalize,                &
               ccpp_physics_init,            &
               ccpp_physics_run,             &
               ccpp_physics_finalize,        &
               ccpp_field_add

! NOTE -- The variables managed by the CCPP are included in the the ccpp_modules.inc file in the "use" statements
#include "ccpp_modules.inc"

  implicit none


  integer                           :: i, j

  ! Create the CCPP required cdata structure
  type(ccpp_t), allocatable, target                      :: cdata(:)

  integer                                                :: ierr
  integer ,parameter :: ncols=1
  integer ,parameter :: nlevs=8
  integer :: ntimes = 0

  character(len=*), parameter :: inputfile = '/terminator-data1/fvitt/micm_inputs/camchem4.cam.h0.0000-01-01-00000.nc'
  type(input_file_type) :: infile
  type(slice_type) :: slice
  real, pointer :: tdata(:,:,:,:)=>null()

  type(const_props_type), pointer :: cnst_info(:)
  integer :: ncnst
  
  character(len=*), parameter :: jsonfile = '/terminator-data1/home/fvitt/ccpp/inputs/tagfileoutput.195.json'
!  character(len=*), parameter :: jsonfile = '/terminator-data1/home/fvitt/ccpp/inputs/tagfileoutput.161.json'
  
  cnst_info => create_cnst_info_array( jsonfile )

  ncnst = size(cnst_info)
  do i = 1,ncnst
     call cnst_info(i)%print()
  enddo
  
  call infile%open( inputfile )

  ntimes = infile%get_ntimes()
  
  ! Allocate the host variables
  allocate(k_rateConst(3))
  allocate(my_co(nlevs))
  allocate(my_o3(nlevs))
  allocate(cdata(ncols))

  ! Initialize the host variables
  state_host%Temperature = 200.
  my_co(:) = 100_kind_phys
  my_o3(:) = 1e-6_kind_phys

  do i = 1, ncols

      ! Use the suite information to setup the run
      call ccpp_init( '../suites/suite_ibox_test_simple1.xml', cdata(i), ierr)
      if (ierr/=0) then
          write(*,'(a,i0,a)') 'An error occurred in ccpp_init for column ', i, '. Exiting...'
          stop
      end if

! use ccpp_fields.inc to call ccpp_field_add for all variables to be exposed to CCPP (this is auto-generated from /src/ccpp/scripts/ccpp_prebuild.py - the script parses tables in the ibox_var_defs.f90)
#include "ccpp_fields.inc"

      ! initialize each column's physics
      call ccpp_physics_init(cdata(i), ierr=ierr)
      if (ierr/=0) then
          write(*,'(a,i0,a)') 'An error occurred in ccpp_physics_init for column ', i, '. Exiting...'
          stop
      end if

  end do

  write(6,*) ' '
  write(6,*) 'After initialization, my_co(1)=',my_co(1)
  write(6,*) 'After initialization, my_o3(1)=',my_o3(1)
  write(6,*) ' '

  slice = infile%set_slice( beglat=0.,endlat=0., beglon=180.,endlon=180., beglev=1000.,endlev=1000.)
  slice%ntimes = 1

  ! loop over all time steps
  do j = 1, ntimes
    slice%begtime = j
    !call slice%print()
    tdata => infile%extract('T',slice)
    state_host%Temperature = tdata(1,1,1,1)
    deallocate(tdata)
    nullify(tdata)
    write(6,*) 'At time step', j, 'in host model state_host%Temperature =', state_host%Temperature
    do i = 1, ncols
       call ccpp_physics_run(cdata(i), ierr=ierr)
       if (ierr/=0) then
           write(*,*) errmsg
           write(*,'(a,i0,a)') 'An error occurred in ccpp_physics_run for column ', i, '. Exiting...'
           stop
       end if

       write(6,*) ' At time j=',j,' my_co(1)=',my_co(1)
       write(6,*) ' At time j=',j,' my_o3(1)=',my_o3(1)
     end do
  end do


  do i=1, ncols
      call ccpp_finalize(cdata(i), ierr)
      if (ierr/=0) then
          write(*,'(a,i0,a)') 'An error occurred in ccpp_finalize for column ', i, '. Exiting...'
          stop
      end if
  end do

end subroutine ibox_main_sub

end module ibox_main

!> \brief Main SCM program that calls the main SCM subroutine
!!
!! The Doxygen documentation system cannot handle in-body comments in Fortran main programs, so the "main" program was put in the
!! subroutine \ref ibox_main_sub above.
program ibox
  use ibox_main
  call ibox_main_sub()
end program ibox
