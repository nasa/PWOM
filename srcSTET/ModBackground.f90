Module ModSeBackground
  implicit none

  private !except

  real, allocatable :: eThermalDensity_IC(:,:),eThermalTemp_IC(:,:)
  
contains
  !subroutines to fill in the neutral atmosphere and thermal plasma
  
  !==========================================================================
  subroutine allocate_background_arrays(nLine,nPoint)
    integer, intent(in) :: nLine, nPoint
    
    if(.not.allocated(eThermalDensity_IC)) &
         allocate(eThermalDensity_IC(nLine,nPoint))
    if(.not.allocated(eThermalTemp_IC)) &
         allocate(eThermalTemp_IC(nLine,nPoint))
  end subroutine allocate_grid_arrays

  !============================================================================
end Module ModSeBackground
