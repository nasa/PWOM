Module ModParticle
  Use ModRandomNumber, ONLY: random_real
  implicit none
  
  private
  
  type particle
     integer :: iSpecies
     real    :: vpar, vperp !velocity in [cm/s]
     real    :: Alt ! position in configurational space [cm]
     logical :: IsOpen   ! defines if particle is in domain or index is avail.
  end type particle

  !array that holds the particles to push
  type(particle), allocatable :: Particles_I(:)

  ! grid variables
  integer :: nAlt    ! points on grid
  integer :: nCells  ! number of cells including ghost cells
  real,allocatable  :: Alt_G(:)  ! position of cell center with ghost [cm]
  real,allocatable  :: dAlt_G(:) ! width of cell with ghost [cm]
  real,allocatable  :: AltBot_F(:)  ! position of cell face [cm]
  real,allocatable  :: AltTop_F(:)  ! position of cell face [cm]
  real,allocatable  :: Volume_G(:)  ! Cell volume

  ! array relating cell index to particle index
  !integer :: iCell_II(:,:)

  ! Electric field variables [statV/cm]
  real,allocatable  :: Efield_G(:) 
  
  ! Time step for moving particles [s]
  real :: DtMove=1.0

  !particle variables
  real,allocatable :: Mass_I(:)
  integer, parameter :: O_=1, H_=2, He_=3
  integer :: nParticle
  
  integer :: nSpecies

  integer :: iSeed = 1
  
  integer :: nNumPerParticle = 1.0e6

  real, parameter :: cBoltzmannCGS = 1.3807E-16

  real :: Time=0.0

  public :: test_sample
contains

  !============================================================================
  subroutine init_particle(nAltIn,AltMin,AltMax,TypeGrid)
    use ModPlanetConst, ONLY: Planet_, NamePlanet_I, rPlanet_I
    integer,intent(in) :: nAltIn
    real,   intent(in) :: AltMin, AltMax
    character(len=100),intent(in):: TypeGrid
    real, parameter :: cGramsPerAMU=1.66054e-24
    real :: rPlanetCM, alpha, Area, AreaBot,AreaTop, dAlt
    integer :: iAlt
    real, parameter :: cMtoCm=1e2
    !--------------------------------------------------------------------------
    !set number of species and mass
    select case(NamePlanet_I(Planet_))
    case('EARTH')
       nSpecies=2
       allocate(Mass_I(nSpecies))
       Mass_I(O_) = cGramsPerAMU*16.0
       Mass_I(H_) = cGramsPerAMU
       !Mass_I(He_) = cGramsPerAMU*4.0
       
    case DEFAULT
       call con_stop('particles not for planet')
    end select
       
    ! Initialize the grid, set size and allocate arrays
    nAlt = nAltIn
    allocate(Alt_G(0:nAlt+1))
    allocate(dAlt_G(0:nAlt+1))
    allocate(AltBot_F(0:nAlt+1))
    allocate(AltTop_F(0:nAlt+1))
    allocate(Volume_G(0:nAlt+1))
    
    ! based on grid type set the grid
    select case(TypeGrid)
    case('Uniform')
       dAlt=(AltMin-AltMax)/nAlt
       dAlt_G=dAlt

       ! set area coef if A=alpha r^3 assuming crossection of 1 cm2 at base
       rPlanetCM=rPlanet_I(Planet_)*cMtoCm

       alpha= 1.0/(rPlanet_I(Planet_)*cMtoCm+AltMin)**3

       do iAlt=0,nAlt+1
          ! set location of cell center and top and bottom faces
          Alt_G(iAlt)=AltMin-dAlt+dAlt*iAlt
          AltBot_F(iAlt)=Alt_G(iAlt)-0.5*dAlt
          AltTop_F(iAlt)=Alt_G(iAlt)+0.5*dAlt
          
          ! calculate the cell volume assuming calculated by assuming 
          !an area function in the form A=alpha r^3
          !and then assuming each cell is a truncated cone.
          AreaBot = alpha*(rPlanetCM+AltBot_F(iAlt))**3
          AreaTop = alpha*(rPlanetCM+AltTop_F(iAlt))**3
          
          Volume_G(iAlt) = 1.0/3.0 * dAlt *&
               ( AreaBot + AreaTop + (AreaBot*AreaTop)**0.5 )
       enddo
    case DEFAULT
       call con_stop('Gridtype not recognized')
    end select
  end subroutine init_particle
  
  !============================================================================
  ! push guiding center of particles
  subroutine push_guiding_center
    use ModPlanetConst, ONLY: Planet_, NamePlanet_I,mPlanet_I,rPlanet_I    
    use ModConst,ONLY:cGravitation
    use ModInterpolate, only: linear
    integer :: iParticle
    real :: rCoord, Efield,acceleration,gravity,AltStart
    real, parameter :: cMtoCm=1e2
    real, parameter :: cElecChargeCGS= 4.80320425e-10 !statcoulombs
    !--------------------------------------------------------------------------
    do iParticle=1,nParticle
       AltStart=Particles_I(iParticle)%Alt
       !set the radial distance and the gravitational acceleration 
       ! at particle location. note need rCoord in cm but modplanetconst in SI
       rCoord=rPlanet_I(Planet_)*cMtoCm+AltStart
       gravity=cMtoCm**3*cGravitation*mPlanet_I(Planet_)/rCoord**2
       
       
       ! interpolate electric field to particle
       Efield= linear(Efield_G(:),0,nAlt,Particles_I(iParticle)%Alt,Alt_G)
       
       !determine acceleration (mirror force-gravity-eField)
       acceleration=1.5*Particles_I(iParticle)%vperp**2/rCoord - gravity&
            +eField*(cElecChargeCGS/Mass_I(Particles_I(iParticle)%iSpecies))
       
       !from initial velocity and acceleration update state
       Particles_I(iParticle)%vpar=Particles_I(iParticle)%vpar&
            +acceleration*DtMove
       Particles_I(iParticle)%Alt=AltStart+Particles_I(iParticle)%vpar*DtMove
       Particles_I(iParticle)%vperp=&
            Particles_I(iParticle)%vperp&
            *(AltStart/Particles_I(iParticle)%Alt)**1.5
    
       !check is particle leaves computational domain
       if (Particles_I(iParticle)%Alt<AltBot_F(1) .or. &
            Particles_I(iParticle)%Alt<AltTop_F(nAlt)) then
          Particles_I(iParticle)%IsOpen = .true.
       endif
          
    end do
  end subroutine push_guiding_center

  !=============================================================================
  ! sample maxwellian in cell
  subroutine sample_maxwellian_cell(iCell,iSpecies,Density,uBulk,Temperature)
    use ModNumConst, ONLY: cPi
    ! index of cell to sample
    integer,intent(in):: iCell

    ! Species to Sample
    integer,intent(in):: iSpecies
    
    ! parameters of maxwellian for sampling
    real, intent(in) :: density ![cm-3]
    real, intent(in) :: uBulk ![cm/s]
    real, intent(in) :: Temperature ![k]
    
    real :: uTherm,uMax,uMin,uRange, uRand, uMostProb,RandNum, PitchAngle
    real :: uMagRel, uPar,uPerp
    real :: fmax, ftemp
    !variable to hold new particle info before inserting it to particles array
    type(particle),allocatable :: NewParticle_I(:)
    integer :: nNew,iParticle
    logical :: DoTest=.true.
    !--------------------------------------------------------------------------
    
    ! Find number of particles to create by taking Ntrue=density*volume and 
    !dividing by nNumPerParticle
    nNew=Density * Volume_G(iCell)/nNumPerParticle
    write(*,*) nNew

    ! allocate array to hold new particles
    allocate(NewParticle_I(nNew))

    ! set bounds on velocity magnitude range to sample from based on thermal 
    !velocity. note that there are very few particles beyond 5*vtherm 
    uTherm=sqrt(8.0*cBoltzmannCGS*Temperature/Mass_I(iSpecies)/cPi)
    
!    uMin = max(uBulk-50.0*uTherm,0.0)
    uMin = uBulk-5.0*uTherm
    uMax = uBulk+5.0*uTherm
    uRange = uMax-uMin
    
    ! Find peak probability by evaluating maxwellian at most probable velocity
    !found when df/dv=0 
    uMostProb = sqrt(2.0*cBoltzmannCGS*Temperature/Mass_I(iSpecies))
    !fmax=maxwellian(Mass_I(iSpecies),uMostProb,Temperature)
    fmax=maxwellian(Mass_I(iSpecies),0.0,Temperature)
    write(*,*) 'fmax',fmax

    !use accept-reject algorithm for distribution sampling
    iParticle=1
    sample_loop: do while (iParticle<nNew)
       
       ! randomly choose velocity in range
       RandNum=random_real(iSeed)
       uRand = RandNum*uRange+uMin
       !now randomly choose pitchangle and azimuth to set
       RandNum=random_real(iSeed)
       PitchAngle=RandNum*cPi
       
       !set par and perp vel
       uPar=uRand*cos(PitchAngle)
       uPerp=abs(uRand*sin(PitchAngle))
       
       !set velocity magnitude relative to bulk
       uMagRel=sqrt((uPar-uBulk)**2+uPerp**2)
       
       !evaluate Maxwellian with this random velocity
       ftemp=maxwellian(Mass_I(iSpecies),uMagRel,Temperature)
       
       !determine probabalisitcally if this is a good choice
       RandNum=random_real(iSeed)
       if (ftemp>RandNum*fmax)then
          !choice found
          
          NewParticle_I(iParticle)%vpar =uPar
          NewParticle_I(iParticle)%vperp=uPerp
          
          !randomly place particle in cell
          RandNum=random_real(iSeed)
          NewParticle_I(iParticle)%Alt=RandNum*dAlt_G(iCell)+AltBot_F(iCell)
          
          !assign remaining particle properties
          NewParticle_I(iParticle)%iSpecies=iSpecies
          NewParticle_I(iParticle)%IsOpen=.true.
          
          !increment particle counter
          iParticle=iParticle+1
       else
          !try again
          cycle sample_loop
       endif
    end do sample_loop
    write(*,*) 'nNew',nNew
    if(DoTest) call plot_distribution_cell(iSpecies,iCell,nNew,NewParticle_I)

    !Now newly sampled particles need to be put into main particle array

    !deallocate to save memory
    deallocate(NewParticle_I)
  end subroutine sample_maxwellian_cell

  !=============================================================================
  ! sample maxwellian in cell
  subroutine sample_maxwellian_cell_boxmuller(iCell,iSpecies,Density,uBulk,Temperature)
    use ModNumConst, ONLY: cPi,cTwoPi
    ! index of cell to sample
    integer,intent(in):: iCell

    ! Species to Sample
    integer,intent(in):: iSpecies
    
    ! parameters of maxwellian for sampling
    real, intent(in) :: density ![cm-3]
    real, intent(in) :: uBulk ![cm/s]
    real, intent(in) :: Temperature ![k]
    
    real :: RandNum1,RandNum2
    real :: uTherm,uMax,uMin,uRange, uRand, uMostProb,RandNum, PitchAngle
    real :: uMagRel, uPar,uPerp
    real :: fmax, ftemp
    !variable to hold new particle info before inserting it to particles array
    type(particle),allocatable :: NewParticle_I(:)
    integer :: nNew,iParticle
    logical :: DoTest=.true.
    !--------------------------------------------------------------------------
    
    ! Find number of particles to create by taking Ntrue=density*volume and 
    !dividing by nNumPerParticle
    nNew=Density * Volume_G(iCell)/nNumPerParticle
    write(*,*) nNew

    ! allocate array to hold new particles
    allocate(NewParticle_I(nNew))
   
    !use accept-reject algorithm for distribution sampling
    iParticle=1
    sample_loop: do while (iParticle<nNew)
       
       ! randomly choose velocity in range
       RandNum1=random_real(iSeed)
       RandNum2=random_real(iSeed)

       uPar=uBulk+sqrt(-2.0*(cBoltzmannCGS/Mass_I(iSpecies))*Temperature*log(RandNum1))*cos(cTwoPi*RandNum2)
       uPerp=sqrt(-2.0*(cBoltzmannCGS/Mass_I(iSpecies))*Temperature*log(RandNum1))*sin(cTwoPi*RandNum2)
          
       NewParticle_I(iParticle)%vpar =uPar
       NewParticle_I(iParticle)%vperp=abs(uPerp)
       
       !randomly place particle in cell
       RandNum=random_real(iSeed)
       NewParticle_I(iParticle)%Alt=RandNum*dAlt_G(iCell)+AltBot_F(iCell)
       
       !assign remaining particle properties
       NewParticle_I(iParticle)%iSpecies=iSpecies
       NewParticle_I(iParticle)%IsOpen=.true.
       
       !increment particle counter
       iParticle=iParticle+1
    end do sample_loop
    write(*,*) 'nNew',nNew
    if(DoTest) call plot_distribution_cell(iSpecies,iCell,nNew,NewParticle_I)

    !Now newly sampled particles need to be put into main particle array

    !deallocate to save memory
    deallocate(NewParticle_I)
  end subroutine sample_maxwellian_cell_boxmuller

  !=============================================================================
  real function maxwellian(mass,vel,Temp)
    use ModNumConst, ONLY: cTwoPi  
    real,intent(in) :: mass,vel,Temp
    
    !maxwellian = sqrt((mass/(cTwoPi*cBoltzmannCGS*Temp))**3)*2.0*cTwoPi*vel**2&
    !     *exp(-mass*vel**2/(2.0*cBoltzmannCGS*Temp))

    maxwellian = (mass/(cTwoPi*cBoltzmannCGS*Temp))**1.5&
         *exp(-mass*vel**2/(2.0*cBoltzmannCGS*Temp))
  end function maxwellian
  !=============================================================================
  
  
  subroutine plot_distribution_cell(iSpecies,iCell,nParticleInCell,CellParticle_I)
    use ModNumConst, ONLY: cPi,cTwoPi
    use ModPlotFile,   ONLY: save_plot_file
    integer,intent(in) :: iSpecies,iCell,nParticleInCell
    type(particle), intent(in) :: CellParticle_I(nParticleInCell)
    
    real :: density,uBulkPar,uBulkPerp,Pressure,Temp, uTherm
    real :: dVel, vParMin, vParMax, vPerpMin, vPerpMax
    integer, parameter :: nVel = 100
    real :: vPar_C(nVel), vPerp_C(nVel)
    integer :: iVel, iParticle, iVpar,iVperp
    real, allocatable   :: Coord_DII(:,:,:), PlotState_IIV(:,:,:)
    !grid parameters
    integer, parameter :: nDim =2, Vperp_=1,Vpar_=2,nVar=1, PSD_=1
    !plot variables
    character(len=100) :: NamePlot
    character(len=100),parameter :: NamePlotVar='Vperp[cm/s] Vpar[cm/s] Particles g r'
    character(len=*),parameter :: NameHeader='distribution function'
    character(len=5) :: TypePlot='ascii'
    logical, save :: IsFirstCall=.true.
    !---------------------------------------------------------------------------

    allocate(Coord_DII(nDim,nVel,nVel),PlotState_IIV(nVel,nVel,nVar))

    ! get moments in cell so we can calculate thermal velocity
    call calc_moments_cell(iSpecies,iCell,nParticleInCell,&
         CellParticle_I,density,uBulkPar,uBulkPerp,Pressure,Temp)
    
    write(*,*) 'density,uBulkPar,uBulkPerp,Pressure,Temp'&
         ,density,uBulkPar,uBulkPerp,Pressure,Temp
    ! calculate thermal velocity

    !kludge
    Temp=1000.

    uTherm=sqrt(8.0*cBoltzmannCGS*Temp/Mass_I(iSpecies)/cPi)

    ! discretize velocity space, center around bulk velocity to 
    !5 uTherm in every direction with grid size .1 uTherm
    vParMin=uBulkPar-5.0*uTherm
    vParMax=uBulkPar+5.0*uTherm
    vPerpMin=uBulkPerp-5.0*uTherm
    vPerpMax=uBulkPerp+5.0*uTherm
    dVel=0.1*uTherm
    do iVel=1,nVel
       vPar_C(iVel)=vParMin+iVel*dVel
       vPerp_C(iVel)=vPerpMin+iVel*dVel
    enddo

    do iVpar=1,nVel
       do iVperp=1,nVel
          Coord_DII(Vperp_,iVperp,iVpar)=VPerp_C(iVperp)
          Coord_DII(Vpar_,iVperp,iVpar)=VPar_C(iVpar)
       enddo
    enddo
    
    PlotState_IIV=0.0
    !Sort particles into bins
    do iParticle=1,nParticleInCell
       iVpar =floor((CellParticle_I(iParticle)%vpar -vParMin )/dVel)
       iVperp=floor((CellParticle_I(iParticle)%vperp-vPerpMin)/dVel)
       
       if (iVpar>0 .and.iVpar<=nVel .and.iVperp>0 .and.iVperp<=nVel) then
          PlotState_IIV(iVperp,iVpar,PSD_)=&
               PlotState_IIV(iVperp,iVpar,PSD_)+nNumPerParticle

       endif
    enddo
    
    ! divide the particles in each velocity bin by the density and the 
    !velocity space cell volume. d3v = dvpar*dvperp*2pi*vperp assuming 
    !rings in velocity space around the vperp axis.
    PlotState_IIV(:,:,PSD_) = PlotState_IIV(:,:,PSD_)&
         /nParticleInCell/nNumPerParticle!/(cTwoPi*Coord_DII(Vperp_,:,:)**2*dVel**2)

    !Plot 
    write(NamePlot,"(a,i4.4,a)") 'DistFuncCell_',iCell,'.out'
    if(IsFirstCall) then
       call save_plot_file(NamePlot, TypePositionIn='rewind', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn=0,TimeIn=time,     &
            nDimIn=nDim,CoordIn_DII=Coord_DII,                &
            VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/))
       IsFirstCall = .false.
    else
       call save_plot_file(NamePlot, TypePositionIn='append', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn=0,TimeIn=time,     &
            nDimIn=nDim,CoordIn_DII=Coord_DII,                &
            VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/))
    endif
    
    deallocate(Coord_DII, PlotState_IIV)
  end subroutine plot_distribution_cell

  !=============================================================================
  subroutine calc_moments_cell(iSpecies,iCell,nParticleInCell,CellParticle_I,&
       density,uBulkPar,uBulkPerp,Pressure,Temp)
    integer,intent(in) :: iSpecies,iCell, nParticleInCell
    type(particle), intent(in) :: CellParticle_I(nParticleInCell)
    real, intent(out) :: density,uBulkPar,uBulkPerp,Pressure,Temp
    !--------------------------------------------------------------------------
    density=nNumPerParticle/Volume_G(iCell)*nParticleInCell
    uBulkPar  = sum(CellParticle_I(:)%vpar) /nParticleInCell
!    uBulkPerp = sum(CellParticle_I(:)%vperp)/nParticleInCell
    uBulkPerp = 0.0
    
    !from eq6.12 Gombosi Gas Kinetic Theory book 
    Pressure = Mass_I(iSpecies)*density&
         *sum(((CellParticle_I(:)%vpar-uBulkPar)**2&
         +2.0*(CellParticle_I(:)%vperp-uBulkPerp)**2))&
         /(3.0*nParticleInCell)

    !get temperature from ideal gas law P=nkT
    Temp = Pressure/(density*cBoltzmannCGS)

  end subroutine calc_moments_cell

  !============================================================================
  ! unit test subroutine for sampling
  subroutine test_sample
    integer :: nAltIn, iCell, iSpecies
    real :: AltMin, AltMax, Density, uBulk, Temperature
    character(len=100):: TypeGrid
    !--------------------------------------------------------------------------
    nAltIn=2
    AltMin=1000.0e5
    AltMin=1020.0e5
    TypeGrid='Uniform'
    iCell = 1
    iSpecies=1
    Density=1e5
    uBulk=0.0
    Temperature=1000.0
    write(*,*) 'init_particle'
    call init_particle(nAltIn,AltMin,AltMax,TypeGrid)

    write(*,*) 'test maxwellian function'
    write(*,*) maxwellian(Mass_I(iSpecies),0.0,Temperature)
    write(*,*) maxwellian(Mass_I(iSpecies),1.0e5,Temperature)
    

    write(*,*) 'sample_maxwellian_cell'
    !call sample_maxwellian_cell(iCell,iSpecies,Density,uBulk,Temperature)
    call sample_maxwellian_cell_boxmuller(iCell,iSpecies,Density,uBulk,Temperature)
    
    
  end subroutine test_sample


end Module ModParticle
