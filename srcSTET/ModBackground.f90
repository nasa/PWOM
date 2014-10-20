Module ModSeBackground
  implicit none

  private !except

  real, public,allocatable :: eThermalDensity_IC(:,:),eThermalTemp_IC(:,:)

  real, public,allocatable :: mLat_I(:), mLon_I(:), gLat1_I(:), gLon1_I(:), &
                              gLat2_I(:), gLon2_I(:)


  real, public :: UT = 0.0
  integer      :: Idate=97046 !day in the form of YYDDD

  ! when the dipole and rotation axis are aligned
  logical      :: DoAlignDipoleRot = .false.

  ! variables for thinning the topside ionosophere 
  real :: facn, fact

  ! exponent for extending solution above IRI or PWOM solution 
  real :: ZEP

  public :: allocate_background_arrays
  public :: fill_thermal_plasma_empirical
  public :: background_test
contains
  !subroutines to fill in the neutral atmosphere and thermal plasma
  
  !=============================================================================
  subroutine fill_thermal_plasma_empirical(iLine,AP,F107,F107A,t)
    use ModSeGrid, only: nIono1,nIono2,nIono,nPlas, nPoint, &
                         FieldLineGrid_IC,Bfield_IC
    
    integer, intent(in) :: iLine
    real   , intent(in) :: AP, F107, F107A, t
    integer :: iIono,iIono2,iPlas,nTopIono1,nTopIono2
    real    :: factor

    !  Use IRI to fill the thermal plasma
    IF (facn.GE.0. .OR. t.EQ.0) THEN
       write(*,*) 'calling get_iri'
       CALL get_iri(iLine,F107A)
       write(*,*) 'finish get_iri'
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
             eThermalDensity_IC(iLine,iIono2)=&
                  eThermalDensity_IC(iLine,iIono2)/factor
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
               * (Bfield_IC(iLine,iPlas)/Bfield_IC(iLine,nTopIono1))**ZEP
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
  !*  Subroutine get_iri calls IRI-90 for each ionosphere.
  SUBROUTINE get_iri(iLine,F107A)
    use ModSeGrid, only: nIono,nPoint,FieldLineGrid_IC
    
    real,    intent(in) :: F107A
    integer, intent(in) :: iLine

    real    :: RZ12 ! parameter in MSIS equal to negative of F107A
    integer :: MMDD ! -day of year
    integer :: i,j  !Generic loop indices
    real    :: STL1, STL2 ! Solar Local Time in iono 1 or 2
    
    ! IRI has output for a number of variables for a provided alt array
    real,allocatable   :: IriOutput_VC(:,:)
    integer, parameter :: nIriOutputs = 11
    real               :: OARR(30) !additional IRI outputs
    integer :: JMAG
    LOGICAL JF(12)
 
    ! unit conversion parameters
    real, parameter :: KtoeV=8.6149E-5, PerM3toPerCm3=1.E-6
    !---------------------------------------------------------------------------

    ! Set IRI output array
    if(.not.allocated(IriOutput_VC)) allocate(IriOutput_VC(nIriOutputs,nIono))

    !  Set up IRI inputs
    do i=1,12
       JF(i)=.TRUE.
    end do
    JF(4)=.FALSE.
    JF(5)=.FALSE.
    JMAG=0
    RZ12=-F107A
    MMDD=-(Idate-(Idate/1000)*1000)
    
    !\
    ! Work on first ionosphere
    !/
    !  Find the geographic coordinates for our geomagnetic coordinates
    write(*,*) 'calling geomag'
    CALL GEOMAG(1,gLon1_I(iLine),gLat1_I(iLine),mLon_I(iLine),mLat_I(iLine))
    write(*,*) 'finish geomag'
    IF (DoAlignDipoleRot) THEN
       gLat1_I(iLine)=mLat_I(iLine)      ! Use these two lines if you want the
       gLon1_I(iLine)=mLon_I(iLine)      ! magnetic and geographic poles aligned
    END IF

    !  Calculate the local solar time
    STL1=(UT/240+gLon1_I(iLine))/15
    IF (STL1.LT.0.) STL1=STL1+24.
    IF (STL1.GT.24.) STL1=STL1-24.
    write(*,*) 'calling iri'
    CALL IRI90(JF,JMAG,gLat1_I(iLine),gLon1_I(iLine),RZ12,MMDD,STL1, &
         FieldLineGrid_IC(iLine,i)/1e5,nIono,' ',IriOutput_VC,OARR)
    write(*,*) 'finish iri'
    do i=nIono,1,-1
       eThermalDensity_IC(iLine,i)=IriOutput_VC(1,i)*PerM3toPerCm3 
       IF (IRIOUTPUT_VC(4,i).LT.0.) IriOutput_VC(4,i)=IriOutput_VC(4,i+1)
       eThermalTemp_IC(iLine,i)=IriOutput_VC(4,i)*KtoeV
    end do
    
    !\
    ! Work on second ionosphere
    !/
    CALL GEOMAG(1,gLon2_I(iLine),gLat2_I(iLine),mLon_I(iLine),-mLat_I(iLine))
    IF (DoAlignDipoleRot) THEN
       gLat2_I(iLine)=mLat_I(iLine)      ! Use these two lines if you want t
       gLon2_I(iLine)=mLon_I(iLine)      ! magnetic and geographic poles aligne
    END IF
    
    !  Calculate the local solar time
    STL2=(UT/240+gLon2_I(iLine))/15
    IF (STL2.LT.0.) STL2=STL2+24.
    IF (STL2.GT.24.) STL2=STL2-24.
    
    !  Call IRI for the second ionosphere
    CALL IRI90(JF,JMAG,gLat2_I(iLine),gLon2_I(iLine),RZ12,MMDD,STL1, &
         FieldLineGrid_IC(iLine,i)/1e5,nIono,' ',IriOutput_VC,OARR)
    
    do i=nIono,1,-1
       j=nPoint-i+1
       eThermalDensity_IC(iLine,j)=IriOutput_VC(1,i)*PerM3toPerCm3
       IF (IriOutput_VC(4,i).LT.0.) IriOutput_VC(4,i)=IriOutput_VC(4,i+1)
       eThermalTemp_IC(iLine,j)=IriOutput_VC(4,i)*KtoeV
    end do
    
  end SUBROUTINE get_iri
  
  !============================================================================
  
  subroutine plot_background
    use ModSeGrid,     ONLY: FieldLineGrid_IC, nLine, nPoint

    use ModIoUnit,     ONLY: UnitTmp_
    use ModPlotFile,   ONLY: save_plot_file
    use ModNumConst,   ONLY: cRadToDeg
    real, allocatable   :: Coord_I(:), PlotState_IV(:,:)
    integer, parameter :: nDim =1, nVar=2, eDens_=1,eTemp_=2
    character(len=100),parameter :: NamePlotVar='S ne te g r'
    character(len=100) :: NamePlot
    character(len=*),parameter :: NameHeader='background output'
    character(len=5) :: TypePlot='ascii'
    integer :: iLine,iPoint
    logical,save :: IsFirstCall =.true.
    !--------------------------------------------------------------------------
    allocate(Coord_I(nPoint), PlotState_IV(nPoint,nVar))
    
    do iLine=1,nLine
       PlotState_IV = 0.0
       Coord_I     = 0.0
       
       !Set Coordinates along field line and PA
       do iPoint=1,nPoint
          Coord_I(iPoint) = FieldLineGrid_IC(iLine,iPoint)/6375.0e5
          PlotState_IV(iPoint,eDens_) = eThermalDensity_IC(iLine,iPoint)
          PlotState_IV(iPoint,eTemp_) = eThermalTemp_IC(iLine,iPoint)
       enddo
       
       ! set name for plotfile
       write(NamePlot,"(a,i4.4,a)") 'background_iLine',iLine,'.out'
       
       !Plot grid for given line. Overwrite old results on firstcall
       if(IsFirstCall) then
          call save_plot_file(NamePlot, TypePositionIn='rewind', &
               TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
               NameVarIn = NamePlotVar, nStepIn= 1,TimeIn=1.0,     &
               nDimIn=nDim,CoordIn_I=Coord_I,                &
               VarIn_IV = PlotState_IV, ParamIn_I = (/1.6, 1.0/))
          IsFirstCall = .false.
       else
          call save_plot_file(NamePlot, TypePositionIn='append', &
               TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
               NameVarIn = NamePlotVar, nStepIn= 1,TimeIn=1.0,     &
               nDimIn=nDim,CoordIn_I=Coord_I,                &
               VarIn_IV = PlotState_IV, ParamIn_I = (/1.6, 1.0/))
       end if
    end do
    
    deallocate(Coord_I, PlotState_IV)
  end subroutine plot_background


  !============================================================================

  !============================================================================
  ! UNIT test for SE update states
  subroutine background_test
    use ModSeGrid, only:create_se_test_grid,nLine,nPoint,nIono,nPlas
    
    integer :: iLine=1, flag=1
    logical :: DoSavePreviousAndReset = .true.
    
    real :: Ap =0.0, F107=80, F107A=80, t=0
    !--------------------------------------------------------------------------

    ! First set up the grid that we will update the state in (this is the same 
    ! as the unit test for the grid).
    write(*,*) 'creating grid'
    call create_se_test_grid

    ! Allocate the background right
    write(*,*) 'allocating background arrays'
    call allocate_background_arrays(nLine,nPoint)
    
    ! set location of field line
    mLat_I(iLine)=60.0 
    mLon_I(iLine)=0.0
    
    ! Fill the background arrays
    write(*,*) 'filling background arrays'
    call fill_thermal_plasma_empirical(iLine,AP,F107,F107A,t)

    ! plot initial state
    call plot_background

  end subroutine background_test
  !=============================================================================
  subroutine allocate_background_arrays(nLine,nPoint)
    integer, intent(in) :: nLine, nPoint
    
    if(.not.allocated(eThermalDensity_IC)) &
         allocate(eThermalDensity_IC(nLine,nPoint))
    if(.not.allocated(eThermalTemp_IC)) &
         allocate(eThermalTemp_IC(nLine,nPoint))
    if(.not.allocated(mLat_I)) &
         allocate(mLat_I(nLine))
    if(.not.allocated(mLon_I)) &
         allocate(mLon_I(nLine))
    if(.not.allocated(gLat1_I)) &
         allocate(gLat1_I(nLine))
    if(.not.allocated(gLon1_I)) &
         allocate(gLon1_I(nLine))
    if(.not.allocated(gLat2_I)) &
         allocate(gLat2_I(nLine))
    if(.not.allocated(gLon2_I)) &
         allocate(gLon2_I(nLine))
    
  end subroutine allocate_background_arrays
  
  !=============================================================================
end Module ModSeBackground
