Module ModParticle
  Use ModRandomNumber, ONLY: random_real
  implicit none
  
  private
  
  !basic particle type
  type particle
     integer :: iSpecies
     integer :: iCell ! cell index for particle
     real    :: vpar, vperp !velocity in [cm/s]
     real    :: Alt ! position in configurational space [cm]
     logical :: IsOpen   ! defines if particle is in domain or index is avail.
  end type particle

  !type for holding pointer for particles of a specific species in a cell
  type particleCellSpecies
     type(particle),pointer :: Particle 
  end type particleCellSpecies

  !array that holds the particles to push
  type(particle), target,allocatable :: Particles_I(:)
  
  !pointer to reference the particles of a particular species in a particular 
  !cell.
  type(particleCellSpecies), allocatable:: SortParticles_III(:,:,:)
  !number of particles of particluar type in cell for a given species
  integer,allocatable :: nSortedParticle_II(:,:)

  !total number of particles in the simulation
  integer :: nParticle
  
  !maximum number of particles before throwing an error (not set yet)
  integer :: MaxParticles
  
  
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

  
  integer :: nSpecies

  character(len=2),allocatable :: NameSpecies_I(:) 
  
  integer :: iSeed = 1
  
  integer :: nNumPerParticle = 1.0e6

  real, parameter :: cBoltzmannCGS = 1.3807E-16

  real :: Time=0.0

  public :: test_sample
  public :: test_pusher
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
       allocate(NameSpecies_I(nSpecies))
       NameSpecies_I(O_)='O_'
       NameSpecies_I(H_)='H_'
    case DEFAULT
       call con_stop('particles not for planet')
    end select
       

    ! Initialize the grid, set size and allocate arrays. 
    ! note -1 and nAlt+2 ghost cells only needed for E field interpolation.
    ! 0 and nAlt ghost cells needed for sampling.
    nAlt = nAltIn
    nCells=nAlt+4
    
    allocate(Alt_G(-1:nAlt+2))
    allocate(dAlt_G(-1:nAlt+2))
    allocate(AltBot_F(-1:nAlt+2))
    allocate(AltTop_F(-1:nAlt+2))
    allocate(Volume_G(-1:nAlt+2))
    
    allocate(Efield_G(-1:nAlt+2))
    
    ! allocate array to tell number of particles of a given type in each cell
    write(*,*) nAlt+1
    allocate(nSortedParticle_II(nSpecies,0:nAlt+1))


    ! based on grid type set the grid
    select case(TypeGrid)
    case('Uniform')
       dAlt=(AltMax-AltMin)/nAlt
       dAlt_G=dAlt

       ! set area coef if A=alpha r^3 assuming crossection of 1 cm2 at base
       rPlanetCM=rPlanet_I(Planet_)*cMtoCm

       alpha= 1.0/(rPlanet_I(Planet_)*cMtoCm+AltMin)**3

       do iAlt=-1,nAlt+2
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
    integer :: iParticle,iAlt
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
       Efield= linear(Efield_G(:),-1,nAlt+1,Particles_I(iParticle)%Alt,Alt_G)
             
       !determine acceleration (mirror force-gravity-eField)
       acceleration=1.5*Particles_I(iParticle)%vperp**2/rCoord - gravity&
            +Efield*(cElecChargeCGS/Mass_I(Particles_I(iParticle)%iSpecies))
       
       
       !from initial velocity and acceleration update state
       Particles_I(iParticle)%vpar=Particles_I(iParticle)%vpar&
            +acceleration*DtMove
       Particles_I(iParticle)%Alt=AltStart+Particles_I(iParticle)%vpar*DtMove
       Particles_I(iParticle)%vperp=&
            Particles_I(iParticle)%vperp&
            *(AltStart/Particles_I(iParticle)%Alt)**1.5


       !check is particle leaves computational domain
       if (Particles_I(iParticle)%Alt<AltBot_F(1) .or. &
            Particles_I(iParticle)%Alt>AltTop_F(nAlt)) then
          Particles_I(iParticle)%IsOpen = .true.
       endif
       
       !Assign cell index to particle
       CELL_ASSIGN: do iAlt=0,nAlt+1
          if(Particles_I(iParticle)%Alt>AltBot_F(iAlt) &
               .and. Particles_I(iParticle)%Alt<AltTop_F(iAlt)) then
             Particles_I(iParticle)%iCell=iAlt
             exit CELL_ASSIGN
          endif
       end do CELL_ASSIGN
    
    end do
  end subroutine push_guiding_center

!  !=============================================================================
!  ! sample maxwellian in cell
!  subroutine sample_maxwellian_cell(iCell,iSpecies,Density,uBulk,Temperature)
!    use ModNumConst, ONLY: cPi
!    ! index of cell to sample
!    integer,intent(in):: iCell
!
!    ! Species to Sample
!    integer,intent(in):: iSpecies
!    
!    ! parameters of maxwellian for sampling
!    real, intent(in) :: density ![cm-3]
!    real, intent(in) :: uBulk ![cm/s]
!    real, intent(in) :: Temperature ![k]
!    
!    real :: uTherm,uMax,uMin,uRange, uRand, uMostProb,RandNum, PitchAngle
!    real :: uMagRel, uPar,uPerp
!    real :: fmax, ftemp
!    !variable to hold new particle info before inserting it to particles array
!    type(particle),allocatable :: NewParticle_I(:)
!    integer :: nNew,iParticle
!    logical :: DoTest=.true.
!    !--------------------------------------------------------------------------
!    
!    ! Find number of particles to create by taking Ntrue=density*volume and 
!    !dividing by nNumPerParticle
!    nNew=Density * Volume_G(iCell)/nNumPerParticle
!    write(*,*) nNew
!
!    ! allocate array to hold new particles
!    allocate(NewParticle_I(nNew))
!
!    ! set bounds on velocity magnitude range to sample from based on thermal 
!    !velocity. note that there are very few particles beyond 5*vtherm 
!    uTherm=sqrt(8.0*cBoltzmannCGS*Temperature/Mass_I(iSpecies)/cPi)
!    
!!    uMin = max(uBulk-50.0*uTherm,0.0)
!    uMin = uBulk-5.0*uTherm
!    uMax = uBulk+5.0*uTherm
!    uRange = uMax-uMin
!    
!    ! Find peak probability by evaluating maxwellian at most probable velocity
!    !found when df/dv=0 
!    uMostProb = sqrt(2.0*cBoltzmannCGS*Temperature/Mass_I(iSpecies))
!    !fmax=maxwellian(Mass_I(iSpecies),uMostProb,Temperature)
!    fmax=maxwellian(Mass_I(iSpecies),0.0,Temperature)
!    write(*,*) 'fmax',fmax
!
!    !use accept-reject algorithm for distribution sampling
!    iParticle=1
!    sample_loop: do while (iParticle<nNew)
!       
!       ! randomly choose velocity in range
!       RandNum=random_real(iSeed)
!       uRand = RandNum*uRange+uMin
!       !now randomly choose pitchangle and azimuth to set
!       RandNum=random_real(iSeed)
!       PitchAngle=RandNum*cPi
!       
!       !set par and perp vel
!       uPar=uRand*cos(PitchAngle)
!       uPerp=abs(uRand*sin(PitchAngle))
!       
!       !set velocity magnitude relative to bulk
!       uMagRel=sqrt((uPar-uBulk)**2+uPerp**2)
!       
!       !evaluate Maxwellian with this random velocity
!       ftemp=maxwellian(Mass_I(iSpecies),uMagRel,Temperature)
!       
!       !determine probabalisitcally if this is a good choice
!       RandNum=random_real(iSeed)
!       if (ftemp>RandNum*fmax)then
!          !choice found
!          
!          NewParticle_I(iParticle)%vpar =uPar
!          NewParticle_I(iParticle)%vperp=uPerp
!          
!          !randomly place particle in cell
!          RandNum=random_real(iSeed)
!          NewParticle_I(iParticle)%Alt=RandNum*dAlt_G(iCell)+AltBot_F(iCell)
!          
!          !assign remaining particle properties
!          NewParticle_I(iParticle)%iSpecies=iSpecies
!          NewParticle_I(iParticle)%IsOpen=.true.
!          
!          !increment particle counter
!          iParticle=iParticle+1
!       else
!          !try again
!          cycle sample_loop
!       endif
!    end do sample_loop
!    write(*,*) 'nNew',nNew
!    if(DoTest) call plot_distribution_cell(iSpecies,iCell,nNew,NewParticle_I)
!
!    !Now newly sampled particles need to be put into main particle array
!
!    !deallocate to save memory
!    deallocate(NewParticle_I)
!  end subroutine sample_maxwellian_cell

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
    integer :: nNew,iParticle, iParticleOld
    logical :: DoTest=.true.
    
    !variables for assignment of particles
    type(particle),allocatable ::ParticlesOld_I(:)
    integer :: nParticleOld,nAvail,i,j
    integer,allocatable :: IndexAvail_I(:)
    !--------------------------------------------------------------------------
    
    ! Find number of particles to create by taking Ntrue=density*volume and 
    !dividing by nNumPerParticle
    nNew=Density * Volume_G(iCell)/nNumPerParticle
    !write(*,*) nNew

    ! allocate array to hold new particles
    allocate(NewParticle_I(nNew))
   
    !use accept-reject algorithm for distribution sampling
    iParticle=1
    sample_loop: do while (iParticle<=nNew)
       
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
       NewParticle_I(iParticle)%iCell=iCell
       NewParticle_I(iParticle)%IsOpen=.false.
       
       !increment particle counter
       iParticle=iParticle+1
    end do sample_loop
!    write(*,*) 'nNew',nNew
!    if(DoTest) then 
!       call plot_distribution_cell(iSpecies,iCell,nNew,NewParticle_I)
!       return
!    endif

    ! Now newly sampled particles need to be put into main particle array
    !on the first call there is no Particles array allocated so allocate
    if (.not.allocated(Particles_I))then
       nParticle=nNew
       allocate(Particles_I(nParticle))
       Particles_I=NewParticle_I
    else
       ! when Particles_I already allocated (usual case) then calculate the
       !number of open spots (nAvail) then save Particles_I and allocate a new 
       !Particles_I array to save the good particles and the newly created ones
       allocate(IndexAvail_I(nParticle))
       where(Particles_I%IsOpen)
          IndexAvail_I=1
       elsewhere
          IndexAvail_I=0
       end where
       nAvail=sum(IndexAvail_I)
       deallocate(IndexAvail_I)
       
       !save old particle array information
       nParticleOld=nParticle
       allocate(ParticlesOld_I(nParticleOld))
       ParticlesOld_I=Particles_I
       deallocate(Particles_I)
       
       !allocate new Particles_I array
       nParticle=nParticleOld-nAvail+nNew
       allocate(Particles_I(nParticle))

       !now fill new Particle_I array with old particles that are not open and 
       !the newly created particles
       iParticle=1
       do iParticleOld=1,nParticleOld
          if(.not.ParticlesOld_I(iParticleOld)%IsOpen) then
             Particles_I(iParticle)=ParticlesOld_I(iParticleOld)
             iParticle=iParticle+1
          endif
       end do
       Particles_I(nParticleOld-nAvail+1:nParticle)=NewParticle_I

       deallocate(ParticlesOld_I)
    endif

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
  
  
  subroutine plot_distribution_cell(iSpecies,iCell)
    use ModNumConst, ONLY: cPi,cTwoPi
    use ModPlotFile,   ONLY: save_plot_file
    integer,intent(in) :: iSpecies,iCell
    integer :: nParticleInCell
    
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

    nParticleInCell=nSortedParticle_II(iSpecies,iCell)
    if(nParticleInCell==0)return
    allocate(Coord_DII(nDim,nVel,nVel),PlotState_IIV(nVel,nVel,nVar))

    ! get moments in cell so we can calculate thermal velocity
    call calc_moments_cell(iSpecies,iCell,&
         density,uBulkPar,uBulkPerp,Pressure,Temp)
    
    write(*,*) 'density,uBulkPar,uBulkPerp,Pressure,Temp'&
         ,density,uBulkPar,uBulkPerp,Pressure,Temp
    ! calculate thermal velocity

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
       iVpar =floor((SortParticles_III(iSpecies,iCell,iParticle)%Particle%vpar &
            -vParMin )/dVel)
       iVperp=floor((SortParticles_III(iSpecies,iCell,iParticle)%Particle%vperp&
            -vPerpMin)/dVel)
       
       if (iVpar>0 .and.iVpar<=nVel .and.iVperp>0 .and.iVperp<=nVel) then
          PlotState_IIV(iVperp,iVpar,PSD_)=&
               PlotState_IIV(iVperp,iVpar,PSD_)+nNumPerParticle
          
       endif
    enddo
    
    ! divide the particles in each velocity bin by the density and the 
    !velocity space cell volume. d3v = dvpar*dvperp*2pi*vperp assuming 
    !rings in velocity space around the vperp axis.
    if (nParticleInCell>0) then
       PlotState_IIV(:,:,PSD_) = PlotState_IIV(:,:,PSD_)&
            /nParticleInCell/nNumPerParticle!/(cTwoPi*Coord_DII(Vperp_,:,:)**2*dVel**2)
    else
       PlotState_IIV(:,:,PSD_) = 0.0
    endif

    !Plot 
    write(NamePlot,"(a,a,a,i5.5,a)") 'DistFunc_',NameSpecies_I(iSpecies),'Alt',floor(Alt_G(iCell)*1e-5),'km.out'
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
  subroutine plot_distribution_cell_orig(iSpecies,iCell,nParticleInCell,CellParticle_I)
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
    call calc_moments_cell_orig(iSpecies,iCell,nParticleInCell,&
         CellParticle_I,density,uBulkPar,uBulkPerp,Pressure,Temp)
    
    write(*,*) 'density,uBulkPar,uBulkPerp,Pressure,Temp'&
         ,density,uBulkPar,uBulkPerp,Pressure,Temp
    ! calculate thermal velocity

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
  end subroutine plot_distribution_cell_orig

  !=============================================================================
  subroutine calc_moments_cell(iSpecies,iCell,&
       density,uBulkPar,uBulkPerp,Pressure,Temp)
    integer,intent(in) :: iSpecies,iCell
    integer :: nParticleInCell,iParticle
    real :: uParTmp,uPerpTmp
    real, intent(out) :: density,uBulkPar,uBulkPerp,Pressure,Temp
    !--------------------------------------------------------------------------
    
    nParticleInCell=nSortedParticle_II(iSpecies,iCell)
    
    !density is particles in cell over volume
    density=nNumPerParticle/Volume_G(iCell)*nParticleInCell

    !bulk velocity is average of velocity
    uBulkPar=0.0
    do iParticle=1,nParticleInCell
       uParTmp =SortParticles_III(iSpecies,iCell,iParticle)%Particle%vpar
       uBulkPar  = uBulkPar + uParTmp        
    end do
    uBulkPar=uBulkPar/nParticleInCell
!    uBulkPerp = sum(CellParticle_I(:)%vperp)/nParticleInCell
    uBulkPerp = 0.0
    
    !from eq6.12 Gombosi Gas Kinetic Theory book
    Pressure=0.0
    do iParticle=1,nParticleInCell
       uParTmp =SortParticles_III(iSpecies,iCell,iParticle)%Particle%vpar
       uPerpTmp=SortParticles_III(iSpecies,iCell,iParticle)%Particle%vperp
       Pressure = Pressure+&
            (uParTmp-uBulkPar)**2&
            +2.0*(uPerpTmp-uBulkPerp)**2
    enddo
    Pressure=Pressure*Mass_I(iSpecies)*density/(3.0*nParticleInCell)

    !get temperature from ideal gas law P=nkT
    Temp = Pressure/(density*cBoltzmannCGS)

  end subroutine calc_moments_cell
  !============================================================================
  subroutine calc_moments_cell_orig(iSpecies,iCell,nParticleInCell,CellParticle_I,&
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

  end subroutine calc_moments_cell_orig

  !============================================================================
  ! sort particles so they can be referenced by species, cell, and particle 
  ! index in cell. Note that our sorting uses a derive type holding pointers 
  ! so all sorted data is referenced ultimately to the Particles_I array.
  ! so modification to the data in the sort should change the data in the 
  ! main particle array.
  subroutine sort_particles
    !integer and cell indices for which to assign to pointer
    integer :: iCell, iSpecies
    
    integer,allocatable :: Index_I(:),iParticleCell_II(:,:)
    integer :: iParticle,MaxSortedParticle
    !--------------------------------------------------------------------------
    
    ! Find the number of particles of a given species in a cell and build 
    ! an array of indicies connecting the local index and the particle index
    allocate(Index_I(nParticle))
    do iSpecies=1,nSpecies
       do iCell=0,nAlt+1
          where(Particles_I%iCell==iCell .and. Particles_I%iSpecies==iSpecies)
             Index_I=1
          elsewhere
             Index_I=0
          end where
          nSortedParticle_II(iSpecies,iCell)=sum(Index_I)
       enddo
    enddo
    deallocate(Index_I)
    
    !find maximum number of particles in a given cell
    MaxSortedParticle=maxval(nSortedParticle_II)
    
    !deallocate previous pointer and reallocate
    if (allocated(SortParticles_III)) then
       deallocate(SortParticles_III)
    endif
    allocate(SortParticles_III(nSpecies,0:nAlt+1,maxval(nSortedParticle_II)))

    ! loop through particles and assign pointer index to target particle from 
    !global list
    allocate(iParticleCell_II(nSpecies,0:nAlt+1))
    iParticleCell_II=1
    PARTICLE_LOOP: do iParticle=1,nParticle
       iCell=Particles_I(iParticle)%iCell
       if(iCell<0 .or.iCell>nAlt+1) cycle PARTICLE_LOOP
       iSpecies=Particles_I(iParticle)%iSpecies
       if(iSpecies==0)write(*,*)iPArticle,nParticle
       SortParticles_III(iSpecies,iCell,&
            iParticleCell_II(iSpecies,iCell))%Particle& 
            =>Particles_I(iParticle)
       !          write(*,*)'test'
       !          write(*,*) SortParticles_III(iParticleCell)%Particle%vpar,Particles_I(iParticle)%vpar
       iParticleCell_II(iSpecies,iCell)=iParticleCell_II(iSpecies,iCell)+1
       
    enddo PARTICLE_LOOP
    deallocate(iParticleCell_II)
  end subroutine sort_particles
  !============================================================================
  ! unit test subroutine for sampling
  subroutine test_sample
    integer :: nAltIn, iCell, iSpecies
    real :: AltMin, AltMax, Density, uBulk, Temperature
    character(len=100):: TypeGrid
    !--------------------------------------------------------------------------
    nAltIn=2
    AltMin=1000.0e5
    AltMax=1020.0e5
    TypeGrid='Uniform'
    iCell = 1
    iSpecies=1
    Density=1e5
    uBulk=100000.0
    Temperature=10000.0
    write(*,*) 'init_particle'
    call init_particle(nAltIn,AltMin,AltMax,TypeGrid)

    write(*,*) 'sample_maxwellian_cell'
    !call sample_maxwellian_cell(iCell,iSpecies,Density,uBulk,Temperature)
    call sample_maxwellian_cell_boxmuller(iCell,iSpecies,Density,uBulk,Temperature)
!    call plot_distribution_cell_orig(iSpecies,iCell,nParticle,Particles_I)
    call sample_maxwellian_cell_boxmuller(iCell+1,iSpecies,Density,uBulk,Temperature)
    write(*,*) 'test pointer extraction'
!    Particles_I%vperp=Particles_I%vperp*2.0
    call sort_particles
    !now plot from sorted
    call plot_distribution_cell(iSpecies,iCell)
    call plot_distribution_cell(iSpecies,iCell+1)
!    call plot_distribution_cell(iSpecies,iCell,nSortedParticle_II,&
!         SortParticles_III(iSpecies)%CellParticle)
  end subroutine test_sample

  !============================================================================
  ! unit test subroutine for sampling
  subroutine test_pusher
    integer :: nAltIn, iCell, iSpecies,nTime,iTime
    real :: AltMin, AltMax, Density, uBulk, Temperature
    character(len=100):: TypeGrid
    !--------------------------------------------------------------------------
    nTime=100
    nAltIn=10
    AltMin=1000.0e5
    AltMax=1200.0e5
    TypeGrid='Uniform'
    iCell = 0
    iSpecies=1
    Density=1e5
    uBulk=100000.0
    Temperature=10000.0
    write(*,*) 'init_particle'
    call init_particle(nAltIn,AltMin,AltMax,TypeGrid)

    write(*,*) 'initialize the ghost cell'
    call sample_maxwellian_cell_boxmuller(iCell,iSpecies,Density,uBulk,Temperature)
    !push the guiding center 100 times and reinitialize ghost cell each time
    do iTime=1,nTime
       call push_guiding_center
       call sample_maxwellian_cell_boxmuller(iCell,iSpecies,Density,uBulk,Temperature)    
    enddo
!    do iCell=1,nAlt
!       call sample_maxwellian_cell_boxmuller(iCell,iSpecies,Density,uBulk,Temperature*2*(iCell))    
!    enddo
    
    write(*,*) 'Sort the particles'
    call sort_particles

    !now plot distribution function at each alt
    do iCell=0,nAlt
       write(*,*)nSortedParticle_II(iSpecies,iCell)
       call plot_distribution_cell(iSpecies,iCell)
    enddo

  end subroutine test_pusher



end Module ModParticle
