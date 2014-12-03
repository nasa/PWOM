! Initializes the stet code, 

subroutine stet_init(nLineIn,Coord_ID,Ap_I,F107,F107A, TimeIn)
  use ModSeGrid, only:init_se_grid,set_grid_pot,allocate_grid_arrays,&
       nLine,nPoint,nIono,nPlas,Efield_IC,DoIncludePotential,&
       DrIono1,nIono1,DrIono2,nIono2,DrIono3,nIono3,DrIono4,nIono4,TypeGridE,&
       nEnergy,EnergyMax,DeltaE, nTheta_II, Lshell_I, nAngle
  use ModSeBackground,only: allocate_background_arrays,mLat_I,mLon_I, &
       set_footpoint_locations,fill_thermal_plasma_empirical,plot_background,&
       plot_ephoto_prod,DoAlignDipoleRot,get_neutrals_and_pe_spectrum
  use ModNumConst,    ONLY: cDegToRad
  use ModSeState,     ONLY: Time
  implicit none 
  
  integer, intent(in):: nLineIn ! number of lines to be worked on
  real, intent(in) :: Coord_ID(nLineIn,2) ! Lat and Lon in degrees for each line
  
  !Thermospheric inputs
  real, intent(in) :: Ap_I(7),F107,F107A
  
  !Time input
  real, intent(in) :: TimeIn
  
  integer,parameter :: Lat_=1 ,Lon_=2 !named parameters for Coord_ID
  
  integer :: iLine !loop variable for line
  !-----------------------------------------------------------------------------
  
  ! Set nLine to nLineIn
  nLine=nLineIn
  
  ! Set the Time to TimeIn
  Time=TimeIn
  
  !\
  ! Set up the grid
  !/
  ! For now this is the same as the unit test grid.
  write(*,*) 'creating grid'
  !  call create_se_test_grid
  
  DrIono1 = 1e6
  nIono1  = 12
  DrIono2 = 2e6
  nIono2  = 10
  DrIono3 = 3e6
  nIono3  = 10
  DrIono4 = 5e6
  nIono4  = 2    
  
  nIono=nIono1+nIono2+nIono3+nIono4
  ! set the energy parameters for the energy grid
  TypeGridE = 'ConstDE'
  nEnergy=100
  !    nEnergy=94
  EnergyMax=100.5
  DeltaE = 1.0
  
  ! Allocated the grid arrays and populate the bfield, sgrid, and PA grid  
  write(*,*) 'allocating arrays'
  call allocate_grid_arrays
  
  ! default values for theta grid. 
  nTheta_II(:,1)=5
  nTheta_II(:,2)=20
  nTheta_II(:,3)=90
  nTheta_II(:,4)=20
  nAngle = 135

  ! Set the Lshell for each line based on the input latitude
  do iLine=1,nLine
     Lshell_I(iLine) = (cos(Coord_ID(iLine,Lat_)*cDegToRad))**-2.0
  end do

  
  write(*,*) 'calling init_se_grid'
  call init_se_grid

  ! when including a potential a new grid is needed
  if (DoIncludePotential) then
     do iLine=1,nLine
        call set_grid_pot(iLine)
     enddo
  endif

  !\
  ! Set the background arrays, sources, and locations
  !/
  ! Allocate the background right
  write(*,*) 'allocating background arrays'
  call allocate_background_arrays
  
  !align dipole and rotation
  DoAlignDipoleRot = .true.
  
  do iLine=1,nLine
     ! set location of field line from inputs
     mLat_I(iLine)=Coord_ID(iLine,Lat_) 
     mLon_I(iLine)=Coord_ID(iLine,Lon_) 
     
     !set glat and glon coords
     call set_footpoint_locations(iLine)
     write(*,*) 'finished setting footpoints'
     ! Get the neutral atmosphere and photo e production spectrum
     call get_neutrals_and_pe_spectrum(iLine,F107,F107A,AP_I)
          write(*,*) 'finished getting neutrals and pe spec'
     ! Fill the background arrays
     write(*,*) 'filling background arrays'
     call fill_thermal_plasma_empirical(iLine,F107,F107A,Time)
  end do
  

end subroutine stet_init
