Module ModSeBackground
  implicit none

  private !except

  real, public,allocatable :: eThermalDensity_IC(:,:),eThermalTemp_IC(:,:)
  
  public :: allocate_background_arrays
contains
  !subroutines to fill in the neutral atmosphere and thermal plasma
  
  !=============================================================================
  subroutine fill_thermal_plasma_empirical(iLine,AP,F107,F107A,
    
    integer :: iIono,iPlas,nTopIono1+1,nTopIono2

    !  Use IRI to fill the thermal plasma
    IF (facn.GE.0. .OR. t.EQ.0) THEN
       CALL GetIRI(sym)
       
       !  The next few lines are for thinning the topside ionosphere densities
       IF (ABS(facn).GT.0.) THEN
          do iIono=nIono1+nIono2+1,nIono
             iIono2=nPoint-iIono+1
             factor= (FieldLineGrid_IC(iLine,iIono) &
                  - FieldLineGrid_IC(iLine,nIono1+nIono2))&
                  /(FieldLineGrid_IC(iLine,nIono) &
                  - FieldLineGrid_IC(iLine,nIono1+nIono2))
             factor=10**(ALOG10(ABS(facn))*factor)
             eThermalDensity_IC(iLine,iIono)=eThermalDensity_IC(iLine,iIono) &
                  /factor
             eThermalDensity_IC(iLine,iIono2)=eThermalDensity_IC(iLine,iIono2)/factor
          end do
       END IF

       ! Alter the topside ionospheric thermal temperatures
       IF (fact.NE.0.) THEN
          do iIono=nIono1+nIono2+1,nIono
             iIono2=nPoint-iIono+1
             factor= (FieldLineGrid_IC(iLine,iIono) &
                  - FieldLineGrid_IC(iLine,nIono1+nIono2))&
                  /(FieldLineGrid_IC(iLine,nIono) &
                  - FieldLineGrid_IC(iLine,nIono1+nIono2))

             eThermalTemp_IC(iLine,iIono) = eThermalTemp_IC(iLine,iIono) &
                  + factor*(fact-eThermalTemp_IC(iLine,iIono))
             eThermalTemp_IC(iLine,iIono2) = eThermalTemp_IC(iLine,iIono2) &
                  + factor*(fact-eThermalTemp_IC(iLine,iIono2))
          end do
       END IF
       

       ! Fill in the plasmaspheric thermal densities
       nTopIono1 = nIono		! Simplifying notation, not en. index
       nTopIono2 = nIono+nPlas+1	! Simplifying notation, not angle index

       do iPlas=nTopIono1+1,nTopIono2-1
          factor=&
               ( FieldLineGrid_IC(iLine,iPlas) &
               - FieldLineGrid_IC(iLine,nTopIono1) ) &
               / ( FieldLineGrid_IC(iLine,nTopIono2) &
               - FieldLineGrid_IC(iLine,nTopIono1) )
          eThermalDensity_IC(iLine,iPlas) = &
               ((1.-factor)*eThermalDensity_IC(iLine,nTopIono1) + &
               factor*eThermalDensity_IC(iLine,nTopIono2))        &
               * (Bfield_IC(iLine,iPlas)/Bfield(nTopIono1))**ZEP
       end do

    
       ! Fill in the plasmaspheric thermal temperatures
       do iPlas=nTopIono1+1,nTopIono2-1
          factor=&
               ( FieldLineGrid_IC(iLine,iPlas) &
               - FieldLineGrid_IC(iLine,nTopIono1) ) &
               / ( FieldLineGrid_IC(iLine,nTopIono2)&
               - FieldLineGrid_IC(iLine,nTopIono1) )
          eThermalTemp_IC(iLine,iPlas)= & 
               ((1.-factor)*eThermalTemp_IC(iLine,nTopIono1) &
               + factor*eThermalTemp_IC(iLine,nTopIono2))
       end do
    END IF		! End IRI setup

  end subroutine fill_thermal_plasma_empirical
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
