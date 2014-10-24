Module ModSeBackground
  implicit none

  private !except

  real, public,allocatable :: eThermalDensity_IC(:,:),eThermalTemp_IC(:,:)

  real, public,allocatable :: mLat_I(:), mLon_I(:), gLat1_I(:), gLon1_I(:), &
                              gLat2_I(:), gLon2_I(:)


  real, public :: UT = 43200.0 ! default at noon
  integer      :: Idate=97046 !day in the form of YYDDD

  ! when the dipole and rotation axis are aligned
  logical      :: DoAlignDipoleRot = .false.

  ! variables for thinning the topside ionosophere 
  real :: facn, fact

  ! exponent for extending solution above IRI or PWOM solution 
  real :: ZEP

  ! Neutral Atmosphere arrays and variables
  integer, parameter :: nNeutralSpecies = 3
  real, allocatable  :: NeutralDens1_IIC(:,:,:),NeutralDens2_IIC(:,:,:)
  real, allocatable  :: NeutralTemp1_IC(:,:),NeutralTemp2_IC(:,:)
  integer,parameter  :: O_=1, O2_=2, N2_=3

  ! Logical variables for if we are calculating photo e spectrum in iono 1 or 2
  logical :: DoCalcPeIono1=.true., DoCalcPeIono2=.true.

  ! Arrays that hold the photo electron production spectrum in iono 1 or 2
  real, allocatable :: ePhotoProdSpec1_IIC(:,:,:),ePhotoProdSpec2_IIC(:,:,:)

  public :: allocate_background_arrays
  public :: fill_thermal_plasma_empirical
  public :: background_test
contains
  !subroutines to fill in the neutral atmosphere and thermal plasma
  !=============================================================================
  subroutine set_footpoint_locations(iLine)
    integer, intent(in) :: iLine
    !---------------------------------------------------------------------------

    !  Find the geographic coordinates for our geomagnetic coordinates iono1
    CALL GEOMAG(1,gLon1_I(iLine),gLat1_I(iLine),mLon_I(iLine),mLat_I(iLine))
    write(*,*) 'finish geomag'
    IF (DoAlignDipoleRot) THEN
       gLat1_I(iLine)=mLat_I(iLine)      ! Use these two lines if you want the
       gLon1_I(iLine)=mLon_I(iLine)      ! magnetic and geographic poles aligned
    END IF
    
    !  Find the geographic coordinates for our geomagnetic coordinates iono2
    CALL GEOMAG(1,gLon2_I(iLine),gLat2_I(iLine),mLon_I(iLine),-mLat_I(iLine))
    IF (DoAlignDipoleRot) THEN
       gLat2_I(iLine)=mLat_I(iLine)      ! Use these two lines if you want t
       gLon2_I(iLine)=mLon_I(iLine)      ! magnetic and geographic poles aligne
    END IF
  end subroutine set_footpoint_locations
  !=============================================================================
  subroutine fill_thermal_plasma_empirical(iLine,F107,F107A,t)
    use ModSeGrid, only: nIono1,nIono2,nIono,nPlas, nPoint, &
                         FieldLineGrid_IC,Bfield_IC
    
    integer, intent(in) :: iLine
    real   , intent(in) :: F107, F107A, t
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
 
    !  Calculate the local solar time
    STL1=(UT/240+gLon1_I(iLine))/15
    IF (STL1.LT.0.) STL1=STL1+24.
    IF (STL1.GT.24.) STL1=STL1-24.
    write(*,*) 'calling iri'
    CALL IRI90(JF,JMAG,gLat1_I(iLine),gLon1_I(iLine),RZ12,MMDD,STL1, &
         FieldLineGrid_IC(iLine,1:nIono)/1e5,nIono,' ',IriOutput_VC,OARR)
    write(*,*) 'finish iri'
    do i=nIono,1,-1
       eThermalDensity_IC(iLine,i)=IriOutput_VC(1,i)*PerM3toPerCm3 
       IF (IRIOUTPUT_VC(4,i).LT.0.) IriOutput_VC(4,i)=IriOutput_VC(4,i+1)
       eThermalTemp_IC(iLine,i)=IriOutput_VC(4,i)*KtoeV
    end do
    
    !\
    ! Work on second ionosphere
    !/
    
    !  Calculate the local solar time
    STL2=(UT/240+gLon2_I(iLine))/15
    IF (STL2.LT.0.) STL2=STL2+24.
    IF (STL2.GT.24.) STL2=STL2-24.
    
    !  Call IRI for the second ionosphere
    CALL IRI90(JF,JMAG,gLat2_I(iLine),gLon2_I(iLine),RZ12,MMDD,STL1, &
         FieldLineGrid_IC(iLine,1:nIono)/1e5,nIono,' ',IriOutput_VC,OARR)
    
    do i=nIono,1,-1
       j=nPoint-i+1
       eThermalDensity_IC(iLine,j)=IriOutput_VC(1,i)*PerM3toPerCm3
       IF (IriOutput_VC(4,i).LT.0.) IriOutput_VC(4,i)=IriOutput_VC(4,i+1)
       eThermalTemp_IC(iLine,j)=IriOutput_VC(4,i)*KtoeV
    end do
    
  end SUBROUTINE get_iri
  !============================================================================
  ! subroutine that fills the neutral atmosphere and PE production spectrum
  subroutine get_neutrals_and_pe_spectrum(iLine,F107,F107A,AP)
    use ModSeGrid,      only: nIono,nEnergy,nPoint,FieldLineGrid_IC
    use ModSeProduction,only: RCOLUM,ESPEC,SOLZEN,SSFLUX
    use ModSeCross,     only: cross
    use EUA_ModMsis90,  only: GTD6,TSELEC
    use ModNumConst,    only: cDegToRad

    integer, intent(in) :: iLine
    real   , intent(in) :: F107, F107A,AP(7)
    
    integer :: iIono, iEnergy ! loop variables
    
    real    :: STL1, STL2 ! Solar Local Time in iono 1 or 2
    ! production variables
    real    :: SZA1,SZA2
    real,allocatable :: ColumnDens_IC(:,:)
   
    !MSIS variables
    integer,parameter :: msisO_=2, msisO2_=4, msisN2_=3 
    real :: SW(25),DN(8),TN(2)
    DATA sw/8*1.,-1.,16*1./
    !--------------------------------------------------------------------------
 
    !set the solar flux
    CALL SSFLUX(0,F107,F107A,0.,0.,0.,0.,1.)

    !\
    ! Work on ionosphere 1
    !/
    !  Calculate the local solar time
    STL1=(UT/240+gLon1_I(iLine))/15
    IF (STL1.LT.0.) STL1=STL1+24.
    IF (STL1.GT.24.) STL1=STL1-24.
    
    !  Call MSIS to get the neutral densities and temperature
    CALL TSELEC(SW)
    do  iIono=1,nIono
       CALL GTD6(Idate,UT,FieldLineGrid_IC(iLine,iIono)/1e5, &
            gLat1_I(iLine),gLon1_I(iLine),STL1,F107A,F107,AP,48,DN,TN)
       NeutralDens1_IIC(iLine,O_ ,iIono)=DN(msisO_)
       NeutralDens1_IIC(iLine,O2_,iIono)=DN(msisO2_)
       NeutralDens1_IIC(iLine,N2_,iIono)=DN(msisN2_)
       NeutralTemp1_IC (iLine,iIono)=TN(2)
    end do
    
    ! Calculate the solar zenith angle
    CALL SOLZEN(Idate,UT,gLat1_I(iLine),gLon1_I(iLine),SZA1)
    SZA1=SZA1*cDegToRad

    !  Set the slant path column densities for O,O2 and N2
    CALL RCOLUM(SZA1,FieldLineGrid_IC(iLine,1:nIono), &
         NeutralDens1_IIC(iLine,:,:),NeutralTemp1_IC(iLine,:),nIono)    
    
    write(*,*) 'SZA1',SZA1
    !  Calculate the photoelectron production spectrum
    IF ((SZA1.LT.2.).AND.(DoCalcPeIono1)) THEN
       CALL ESPEC(NeutralDens1_IIC(iLine,:,:),ePhotoProdSpec1_IIC(iLine,:,:),&
            nIono,0)
    ELSE
       do iIono=1,nIono
          do iEnergy=1,nEnergy
             ePhotoProdSpec1_IIC(iLine,iEnergy,iIono)=0.
	  end do
       end do
    END IF

    !\
    ! Work on ionosphere 2
    !/
    !  Calculate the local solar time
    STL2=(UT/240+gLon2_I(iLine))/15
    IF (STL2.LT.0.) STL2=STL2+24.
    IF (STL2.GT.24.) STL2=STL2-24.
    
    !  Call MSIS to get the neutral densities and temperature
    CALL TSELEC(SW)
    do  iIono=1,nIono
       CALL GTD6(Idate,UT,FieldLineGrid_IC(iLine,iIono)/1e5, &
            gLat2_I(iLine),gLon2_I(iLine),STL2,F107A,F107,AP,48,DN,TN)
       NeutralDens2_IIC(iLine,O_ ,iIono)=DN(msisO_)
       NeutralDens2_IIC(iLine,O2_,iIono)=DN(msisO2_)
       NeutralDens2_IIC(iLine,N2_,iIono)=DN(msisN2_)
       NeutralTemp2_IC (iLine,iIono)=TN(2)
    end do
    
    ! Calculate the solar zenith angle
    CALL SOLZEN(Idate,UT,gLat2_I(iLine),gLon2_I(iLine),SZA2)
    SZA2=SZA2*cDegToRad

    !  Set the slant path column densities for O,O2 and N2
    CALL RCOLUM(SZA2,FieldLineGrid_IC(iLine,1:nIono), &
         NeutralDens2_IIC(iLine,:,:),NeutralTemp2_IC(iLine,:),nIono)    

    write(*,*) 'SZA1',SZA2    
    !  Calculate the photoelectron production spectrum
    IF ((SZA2.LT.2.).AND.(DoCalcPeIono2)) THEN
       CALL ESPEC(NeutralDens2_IIC(iLine,:,:),ePhotoProdSpec2_IIC(iLine,:,:),&
            nIono,nPoint)
    ELSE
       do iIono=1,nIono
          do iEnergy=1,nEnergy
             ePhotoProdSpec2_IIC(iLine,iEnergy,iIono)=0.
	  end do
       end do
    END IF

    ! set the cross sections (perhaps this should only be called once?)
    call cross
  end subroutine get_neutrals_and_pe_spectrum
  !============================================================================
  
  subroutine plot_background(iLine,nStep,time)
    use ModSeGrid,     ONLY: FieldLineGrid_IC, nLine, nPoint,nIono
    use ModIoUnit,     ONLY: UnitTmp_
    use ModPlotFile,   ONLY: save_plot_file
    use ModNumConst,   ONLY: cRadToDeg

    integer, intent(in) :: iLine,nStep
    real,    intent(in) :: time

    real, allocatable   :: Coord_I(:), PlotState_IV(:,:)
    integer, parameter :: nDim =1, nVar=5, eDens_=1,eTemp_=2, &
                          VarO_=3,VarO2_=4,VarN2_=5
    character(len=100),parameter :: NamePlotVar='S ne te nO nO2 nN2 g r'
    character(len=100) :: NamePlot
    character(len=*),parameter :: NameHeader='background output'
    character(len=5) :: TypePlot='ascii'
    integer :: iPoint
    logical,save :: IsFirstCall =.true.
    !--------------------------------------------------------------------------
    allocate(Coord_I(nPoint), PlotState_IV(nPoint,nVar))
    
    PlotState_IV = 0.0
    Coord_I     = 0.0
    
    !Set Coordinates along field line and PA
    do iPoint=1,nPoint
       Coord_I(iPoint) = FieldLineGrid_IC(iLine,iPoint)/6375.0e5
       PlotState_IV(iPoint,eDens_) = eThermalDensity_IC(iLine,iPoint)
       PlotState_IV(iPoint,eTemp_) = eThermalTemp_IC(iLine,iPoint)
       if (iPoint <= nIono) then
          PlotState_IV(iPoint,VarO_)  = NeutralDens1_IIC(iLine,O_,iPoint)
          PlotState_IV(iPoint,VarO2_) = NeutralDens1_IIC(iLine,O2_,iPoint)
          PlotState_IV(iPoint,VarN2_) = NeutralDens1_IIC(iLine,N2_,iPoint)
       elseif(iPoint>nPoint-nIono) then
          PlotState_IV(iPoint,VarO_)  = &
               NeutralDens2_IIC(iLine,O_,nPoint-iPoint+1)
          PlotState_IV(iPoint,VarO2_) = &
               NeutralDens2_IIC(iLine,O2_,nPoint-iPoint+1)
          PlotState_IV(iPoint,VarN2_) = &
               NeutralDens2_IIC(iLine,N2_,nPoint-iPoint+1)
       else
          !neutral atmosphere not considered in plasmaphere
          PlotState_IV(iPoint,VarO_)  = 0.0
          PlotState_IV(iPoint,VarO2_) = 0.0
          PlotState_IV(iPoint,VarN2_) = 0.0
       endif
    enddo
    
    ! set name for plotfile
    write(NamePlot,"(a,i4.4,a)") 'background_iLine',iLine,'.out'
    
    !Plot grid for given line. Overwrite old results on firstcall
    if(IsFirstCall) then
       call save_plot_file(NamePlot, TypePositionIn='rewind', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn= nStep,TimeIn=time,     &
            nDimIn=nDim,CoordIn_I=Coord_I,                &
            VarIn_IV = PlotState_IV, ParamIn_I = (/1.6, 1.0/))
       IsFirstCall = .false.
    else
       call save_plot_file(NamePlot, TypePositionIn='append', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn= nStep,TimeIn=time,     &
            nDimIn=nDim,CoordIn_I=Coord_I,                &
            VarIn_IV = PlotState_IV, ParamIn_I = (/1.6, 1.0/))
    end if
     
    deallocate(Coord_I, PlotState_IV)
  end subroutine plot_background


  !============================================================================
  !============================================================================
  ! save state plot for verification
  subroutine plot_ephoto_prod(iLine,nStep,time)
    use ModSeGrid,     ONLY: FieldLineGrid_IC,nIono,nEnergy, nPoint, &
         DeltaE_I,EnergyGrid_I
    use ModIoUnit,     ONLY: UnitTmp_
    use ModPlotFile,   ONLY: save_plot_file
    use ModNumConst,   ONLY: cRadToDeg,cPi

    integer, intent(in) :: iLine, nStep
    real,    intent(in) :: time

    real, allocatable   :: Coord_DII(:,:,:), PlotState_IIV(:,:,:)
    real, parameter     :: rEarthCM = 6375.0e5
    !grid parameters
    integer, parameter :: nDim =2, nVar=2, E_=1, S_=2
    integer, parameter :: spec1_=1, spec2_=2

    character(len=100),parameter :: NamePlotVar='E[eV] Alt[km] eProd1[cm-3eV-1s-1sr-1] eProd2[cm-3eV-1s-1sr-1] g r'
    character(len=*),parameter :: NameHeader='ePhoto Production Spectrum output'
    character(len=5) :: TypePlot='ascii'
    integer :: iEnergy,iIono
    character(len=100) :: NamePlot
    logical,save :: IsFirstCall =.true.
    !--------------------------------------------------------------------------
    allocate(Coord_DII(nDim,nEnergy,nIono),PlotState_IIV(nEnergy,nIono,nVar))


       PlotState_IIV = 0.0
       Coord_DII     = 0.0
       
       !Set Coordinates along field line and PA
       do iEnergy=1,nEnergy
          do iIono=1,nIono
             Coord_DII(E_,iEnergy,iIono) = EnergyGrid_I(iEnergy)             
             Coord_DII(S_,iEnergy,iIono) = FieldLineGrid_IC(iLine,iIono)/1e5
             PlotState_IIV(iEnergy,iIono,spec1_)  = &
                  ePhotoProdSpec1_IIC(iLine,iEnergy,iIono)&
                  /4.0/cPi/DeltaE_I(iEnergy)
             PlotState_IIV(iEnergy,iIono,spec2_)  = &
                  ePhotoProdSpec2_IIC(iLine,iEnergy,iIono)&
                  /4.0/cPi/DeltaE_I(iEnergy)
          enddo
       enddo

       ! set name for plotfile
       write(NamePlot,"(a,i4.4,a)") 'ephotoprod_iLine',iLine,'.out'
       
       !Plot grid for given line
       if(IsFirstCall) then
          call save_plot_file(NamePlot, TypePositionIn='rewind', &
               TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
               NameVarIn = NamePlotVar, nStepIn=nStep,TimeIn=time,     &
               nDimIn=nDim,CoordIn_DII=Coord_DII,                &
               VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/))
          IsFirstCall = .false.
       else
          call save_plot_file(NamePlot, TypePositionIn='append', &
               TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
               NameVarIn = NamePlotVar, nStepIn=nStep,TimeIn=time,     &
               nDimIn=nDim,CoordIn_DII=Coord_DII,                &
               VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/))
       endif
    
    deallocate(Coord_DII, PlotState_IIV)
  end subroutine plot_ephoto_prod
  
  !============================================================================
  ! UNIT test for SE update states
  subroutine background_test
    use ModSeGrid, only:create_se_test_grid,nLine,nPoint,nIono,nPlas
    
    integer :: iLine=1, flag=1, nStep=0
    real    :: time=0
    logical :: DoSavePreviousAndReset = .true.
    
    real :: Ap(7), F107=80, F107A=80, t=0
    !--------------------------------------------------------------------------

    ! First set up the grid that we will update the state in (this is the same 
    ! as the unit test for the grid).
    write(*,*) 'creating grid'
    call create_se_test_grid

    ! Allocate the background right
    write(*,*) 'allocating background arrays'
    call allocate_background_arrays
    
    !set zep and thining parameters
    ZEP = 1
    facn =-1
    fact = 0
    
    !align dipole and rotation
    DoAlignDipoleRot = .true.
    
    ! set location of field line
    mLat_I(iLine)=60.0 
    mLon_I(iLine)=0.0
    
    !set glat and glon coords
    call set_footpoint_locations(iLine)
    
    ! Fill the background arrays
    write(*,*) 'filling background arrays'
    call fill_thermal_plasma_empirical(iLine,F107,F107A,t)
    
    ! Get the neutral atmosphere and photo e production spectrum
    AP(:)=4.0
    call get_neutrals_and_pe_spectrum(iLine,F107,F107A,AP)

    ! plot initial state
    call plot_background(iLine,nStep,time)
    
    ! plot ephoto production
    call plot_ephoto_prod(iLine,nStep,time)

  end subroutine background_test
  !=============================================================================
  subroutine allocate_background_arrays
    use ModSeGrid,     ONLY: nLine, nPoint, nIono, nEnergy
        
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

    if(.not.allocated(NeutralDens1_IIC)) &
         allocate(NeutralDens1_IIC(nLine,nNeutralSpecies,nIono))
    if(.not.allocated(NeutralDens2_IIC)) &
         allocate(NeutralDens2_IIC(nLine,nNeutralSpecies,nIono))

    if(.not.allocated(NeutralTemp1_IC)) &
         allocate(NeutralTemp1_IC(nLine,nIono))
    if(.not.allocated(NeutralTemp2_IC)) &
         allocate(NeutralTemp2_IC(nLine,nIono))

    if(.not.allocated(ePhotoProdSpec1_IIC)) &
         allocate(ePhotoProdSpec1_IIC(nLine,nEnergy,nIono))
    if(.not.allocated(ePhotoProdSpec2_IIC)) &
         allocate(ePhotoProdSpec2_IIC(nLine,nEnergy,nIono))


  end subroutine allocate_background_arrays
  



  !=============================================================================
end Module ModSeBackground
