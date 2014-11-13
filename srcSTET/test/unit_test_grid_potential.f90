program unit_test_grid_potential
  use ModSeGrid
  use ModSeMpi

  use ModMPI
  use CON_planet, ONLY: init_planet_const, set_planet_defaults

  integer :: iError
  
  !-----------------------------------------------------------------------------

  !****************************************************************************
  ! Initiallize MPI and get number of processors and rank of given processor
  !****************************************************************************

  write(*,*) 'Initiallizing MPI'

  !---------------------------------------------------------------------------
  call MPI_INIT(iError)
  iComm = MPI_COMM_WORLD

  call MPI_COMM_RANK(iComm,iProc,iError)
  call MPI_COMM_SIZE(iComm,nProc,iError)

  !\
  ! Initialize the planetary constant library and set Earth
  ! as the default planet.
  !/
  write(*,*) 'Initiallizing Planet'

  call init_planet_const
  call set_planet_defaults

  write(*,*) 'starting se_grid_test'

  call grid_potential_test

end program unit_test_grid_potential

!============================================================================
! UNIT test for SE update states
subroutine grid_potential_test
  use ModSeGrid, only:create_se_test_grid,calc_potential,set_energy_bounds,&
       locate_reference_alt_for_mu0, nMu0RefAlt_II&
       nLine,nPoint,nIono,nPoint,nTop,Efield_IC,FieldLineGrid_IC
  use ModSePlot, only: plot_along_field

  integer :: iLine=1, flag=1
  logical :: DoSavePreviousAndReset = .true.
  integer :: nStep
  logical :: IsOpen
  real,parameter :: cCmToM = 1.0e-2
  !--------------------------------------------------------------------------
  
  ! First set up the grid that we will update the state in (this is the same 
  ! as the unit test for the grid).
  write(*,*) 'creating grid'
  call create_se_test_grid

  ! input electric field (assume constant in plasmapshere and zero in iono)
  ! choose values to have a 5 V drop from top of iono to equator
  Efield_IC(iLine,:)=0.0
  Efield_IC(iLine,nIono+1:nTop-1)=5.0&
       /((FieldLineGrid_IC(iLine,nTop)-FieldLineGrid_IC(iLine,nIono))*cCmToM)
  Efield_IC(iLine,nTop+1:nPoint-nIono)=-5.0&
       /((FieldLineGrid_IC(iLine,nTop)-FieldLineGrid_IC(iLine,nIono))*cCmToM)

  ! calculate the potential energy
  call calc_potential(iLine)

  ! set the bounds on the energy calculation as a function of altitude
  call set_energy_bounds(iLine)
     
  ! plot the potential, and min total energy considered
  call plot_along_field

  !get mu0 reference altitude
  call locate_reference_alt_for_mu0
  write(*,*) nMu0RefAlt_II(iLine,:)
end subroutine grid_potential_test



!============================================================================
! The following subroutines are here so that we can use SWMF library routines
! Also some features available in SWMF mode only require empty subroutines
! for compilation of the stand alone code.
!============================================================================
subroutine CON_stop(StringError)
  use ModSeMpi, ONLY : iProc,iComm
  use ModMpi
  implicit none
  character (len=*), intent(in) :: StringError

  ! Local variables:
  integer :: iError,nError
  !----------------------------------------------------------------------------

  write(*,*)'Stopping execution! me=',iProc,&
       ' with msg:'
  write(*,*)StringError
  call MPI_abort(iComm, nError, iError)
  stop

end subroutine CON_stop

subroutine CON_set_do_test(String,DoTest,DoTestMe)
  implicit none
  character (len=*), intent(in)  :: String
  logical          , intent(out) :: DoTest, DoTestMe

  DoTest = .false.; DoTestMe = .false.

end subroutine CON_set_do_test

subroutine CON_io_unit_new(iUnit)

  use ModIoUnit, ONLY: io_unit_new
  implicit none
  integer, intent(out) :: iUnit

  iUnit = io_unit_new()

end subroutine CON_io_unit_new

