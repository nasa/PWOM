Module ModSeBackground
  implicit none

  private !except

  real, public,allocatable :: eThermalDensity_IC(:,:),eThermalTemp_IC(:,:)
  
  public :: allocate_background_arrays
contains
  !subroutines to fill in the neutral atmosphere and thermal plasma
  
  !=============================================================================
  subroutine fill_thermal_plasma_emprical(AP,F107,F107A,
    **  Set up IRI thermal plasma
      IF (facn.GE.0. .OR. t.EQ.0) THEN
	CALL GetIRI(sym)
c  The next few lines are for thinning the topside ionosphere densities
      IF (ABS(facn).GT.0.) THEN
	do i=Inum1+Inum2+1,Iono
	  j=IMAX-i+1
	  factor=(s(i)-s(Inum1+Inum2))/(s(Iono)-s(Inum1+Inum2))
	  factor=10**(ALOG10(ABS(facn))*factor)
	  ZE(i)=ZE(i)/factor
	  ZE(j)=ZE(j)/factor
	end do
      END IF
**  Alter the topside ionospheric thermal temperatures
      IF (fact.NE.0.) THEN
	do i=Inum1+Inum2+1,Iono
	  j=IMAX-i+1
	  factor=(s(i)-s(Inum1+Inum2))/(s(Iono)-s(Inum1+Inum2))
	  TE(i)=TE(i)+factor*(fact-TE(i))
	  TE(j)=TE(j)+factor*(fact-TE(j))
	end do
      END IF
!      write(*,*) 'test2'
**  Fill in the plasmaspheric thermal densities
      j=Iono			! Simplifying notation, not en. index
      k=Iono+Ip+1		! Simplifying notation, not angle index
      IF (ZEP.EQ.1.) THEN	! Here because it's the usual case
	do i=j+1,k-1
	  factor=(s(i)-s(j))/(s(k)-s(j))
	  ZE(i)=((1.-factor)*ZE(j)+factor*ZE(k))*Bfield(i)/Bfield(j)
	end do
      ELSE
	do i=j+1,k-1
	  factor=(s(i)-s(j))/(s(k)-s(j))
	  ZE(i)=((1.-factor)*ZE(j)+factor*ZE(k))
     &		*(Bfield(i)/Bfield(j))**ZEP
	end do
      END IF
    
**  Fill in the plasmaspheric thermal temperatures
      do i=j+1,k-1
	factor=(s(i)-s(j))/(s(k)-s(j))
	TE(i)=((1.-factor)*TE(j)+factor*TE(k))
      end do
      END IF		! End IRI setup

  end subroutine fill_thermal_plasma_emprical
  !=============================================================================
  subroutine allocate_background_arrays(nLine,nPoint)
    integer, intent(in) :: nLine, nPoint
    
    if(.not.allocated(eThermalDensity_IC)) &
         allocate(eThermalDensity_IC(nLine,nPoint))
    if(.not.allocated(eThermalTemp_IC)) &
         allocate(eThermalTemp_IC(nLine,nPoint))
  end subroutine allocate_background_arrays

  !=============================================================================
end Module ModSeBackground
