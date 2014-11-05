!============================================================================
! run STET to get a new steady state or to advance some amount of time
subroutine stet_run
  use ModSeGrid, only:create_se_test_grid,nLine,nPoint,nIono,nPlas
  use ModSeBackground,only: allocate_background_arrays,mLat_I,mLon_I, &
       set_footpoint_locations,fill_thermal_plasma_empirical,plot_background,&
       plot_ephoto_prod,DoAlignDipoleRot,get_neutrals_and_pe_spectrum, &
       eThermalDensity_IC,eThermalTemp_IC,nNeutralSpecies,&
       NeutralDens1_IIC,ePhotoProdSpec1_IIC,NeutralDens2_IIC,ePhotoProdSpec2_IIC
  use ModSeState,only: allocate_state_arrays, iphiup,iphidn,phiup,phidn,&
       initiono,check_time,delt,epsilon,update_se_state_iono,update_se_state,&
       liphiup,liphidn,lphiup,lphidn,specup, specdn, initplas
  use ModSeCross,only: SIGS,SIGI,SIGA    
  use ModSePlot, only: plot_state,plot_omni_iono  
  integer :: iLine=1, flag=1, nStep=0
  real    :: time=0
  logical :: DoSavePreviousAndReset = .true.
  logical :: IsOpen = .false.

  real :: Ap(7), F107=80, F107A=80, t=0
  
  logical,parameter :: IsIono1=.true.
  !--------------------------------------------------------------------------
  
  !\
  ! Set up the grid
  !/
  ! For now this is the same as the unit test grid.
  write(*,*) 'creating grid'
  call create_se_test_grid
  
  !\
  ! Set the background arrays, sources, and locations
  !/
  ! Allocate the background right
  write(*,*) 'allocating background arrays'
  call allocate_background_arrays
  
  !align dipole and rotation
  DoAlignDipoleRot = .true.
  
  ! set location of field line
  mLat_I(iLine)=60.0 
  mLon_I(iLine)=0.0
  
  !set glat and glon coords
  call set_footpoint_locations(iLine)
    
  ! Get the neutral atmosphere and photo e production spectrum
  AP(:)=4.0
  call get_neutrals_and_pe_spectrum(iLine,F107,F107A,AP)
  
  ! Fill the background arrays
  write(*,*) 'filling background arrays'
  call fill_thermal_plasma_empirical(iLine,F107,F107A,t)

  ! plot initial state
  call plot_background(iLine,nStep,time)
  
  ! plot ephoto production
  call plot_ephoto_prod(iLine,nStep,time)
  
  !\
  ! The main update
  !/
  call allocate_state_arrays

  ! Set the timestep and convergence criteria
  delt=1.0e5
  epsilon = 0.4

  ! Define the initial state
  iphiup(iLine,:,:,:)=0.0
  iphidn(iLine,:,:,:)=0.0
  
  liphiup(iLine,:,:,:)=0.0
  liphidn(iLine,:,:,:)=0.0
  
  phiup(iLine,:,:,:)=0.0
  phidn(iLine,:,:,:)=0.0
  
  lphiup(iLine,:,:,:)=0.0
  lphidn(iLine,:,:,:)=0.0
  
  nStep = 0

  !Start Timeloop
  TIME_LOOP: do while (flag == 1)
     ! Initialize the ionosphere
     write(*,*) 'Initializing iono'
     call initiono(iLine,DoSavePreviousAndReset)
     
     ! update the SE state for iono1
     write(*,*) 'update se state for iono1'
     call update_se_state_iono(iLine,IsIono1,eThermalDensity_IC(iLine,:),&
          eThermalTemp_IC(iLine,:),nNeutralSpecies,&
          NeutralDens1_IIC(iLine,:,:),SIGS,SIGI,SIGA,&
          ePhotoProdSpec1_IIC(iLine,:,:))
     
     ! update the SE state for iono2
     write(*,*) 'update se state for iono2'
     call update_se_state_iono(iLine,.not.IsIono1,eThermalDensity_IC(iLine,:),&
          eThermalTemp_IC(iLine,:),nNeutralSpecies,&
          NeutralDens2_IIC(iLine,:,:),SIGS,SIGI,SIGA,&
          ePhotoProdSpec2_IIC(iLine,:,:))
     
     ! Initialize the plasmasphere
     write(*,*) 'Initializing plasmasphere'
     call initplas(1,DoSavePreviousAndReset)
     
     ! update the SE state
     write(*,*) 'update se state'
     call update_se_state(1, eThermalDensity_IC(iLine,:),&
          eThermalTemp_IC(iLine,:),IsOpen)

     ! check convergence
     write(*,*) 'check for convergence'
     call check_time(1,flag)
     
     ! increment step
     nStep=nStep+1
     
  end do TIME_LOOP

  ! plot output
  call plot_state(1,nStep,time,iphiup,iphidn,phiup,phidn)
  
  call plot_omni_iono(iLine,nStep,time,specup,specdn,.true.)
  call plot_omni_iono(iLine,nStep,time,specup,specdn,.false.)
  


end subroutine stet_run
