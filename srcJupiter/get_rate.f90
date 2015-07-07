!******************************************************************************
! get_rate takes Alt as input and returns the rate coef for the 
! H+ + H2 --> H2+ + H reaction. Data from Moses and Bass, 2000
! Created by Alex Glocer, 07/06
!******************************************************************************

Subroutine get_rate( Temp,Rate)

  real, intent(in)  :: Temp
  real, intent(out) :: Rate
  
  real :: k0 = 2.0e-9
  integer :: i
  
  
  real,dimension(15)  :: MeasuredAlt, MeasuredRate
  
  do i=4,8
     Rate = Rate + exp(-real(i)*6100.0/Temp)
  enddo
  
  Rate = k0*Rate
  
end Subroutine get_rate
