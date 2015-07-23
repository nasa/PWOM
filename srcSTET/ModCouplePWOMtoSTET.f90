Module ModCouplePWOMtoSTET
  implicit none
  
  private !except
  
  integer :: nAltPw
  real,allocatable :: AltPw_C(:)
  real :: dAltPw
  
  !public methods
  public :: init_pwom_stet_coupling
  public :: get_stet_for_pwom
  
contains
  !=============================================================================
  ! set up coupling between pwom and stet for all lines. 
  ! thermal e density and temperature for each line and interpolate to STET 
  ! grid 
  subroutine init_pwom_stet_coupling(nAltPwIn,nLinePw,iLineGlobalPw_I,&
       AltPwIn_C,PrecipEminPwIn,PrecipEmaxPwIn,PrecipEmeanPwIn,PrecipEfluxPwIn)
    use ModSeGrid, only: nLine,iLineGlobal_I,DoIncludePotential
    use ModSeBackground,only: allocate_background_arrays,DoAlignDipoleRot,ZEP,&
                              DoUsePWOM
    use ModSeState,only: allocate_state_arrays,PrecipEmin, PrecipEmax, &
         PrecipEmean,PrecipEflux,UsePrecipitation
    integer, intent(in) :: nAltPwIn, nLinePw,iLineGlobalPw_I(nLinePw)
    real,    intent(in) :: AltPwIn_C(nAltPwIn)
    real, optional, intent(in)::PrecipEminPwIn,PrecipEmaxPwIn, &
                                PrecipEmeanPwIn,PrecipEfluxPwIn
    !---------------------------------------------------------------------------
    DoUsePWOM = .true.    
    
    ! Set lines on proc to match lines in PWOM
    nLine = nLinePw
    
    ! Assume thermosphere drops with B2 for PW case
    ZEP=2.0
    
    ! Set up the grid parameters, actual grid will be set later by each line
    
    write(*,*) 'setting grid dimensions'
    call set_grid_dimensions_default
    
    ! Set the global line number
    iLineGlobal_I=iLineGlobalPw_I


    ! Allocate the background arrays
    write(*,*) 'allocating background arrays'
    call allocate_background_arrays
    
    ! align dipole and rotation
    DoAlignDipoleRot = .false.

    ! Include the parallel electric field when using PWOM 
    DoIncludePotential=.true.

    ! store PWOM altitude grid in module for use by interpolation routines
    ! for passing info between STET and PWOM
    nAltPw=nAltPwIn
    if (.not.allocated(AltPw_C)) allocate(AltPw_C(nAltPw))
    AltPw_C=AltPwIn_C
    
    ! save the grid spacing in PWOM
    dAltPw = AltPw_C(2)-AltPw_C(1)

    !allocate the SE state arrays
    call allocate_state_arrays

    ! Set the incomming precipitation
    if (present(PrecipEminPwIn).and.present(PrecipEmaxPwIn) &
         .and.present(PrecipEmeanPwIn).and.present(PrecipEfluxPwIn)) then
       UsePrecipitation = .true.
       PrecipEmax=PrecipEmaxPwIn
       PrecipEmin=PrecipEminPwIn
       PrecipEmean=PrecipEmeanPwIn
       PrecipEflux=PrecipEfluxPwIn
    elseif(present(PrecipEminPwIn).or.present(PrecipEmaxPwIn) &
         .or.present(PrecipEmeanPwIn).or.present(PrecipEfluxPwIn)) then
       call con_stop('PW_error: STET precip update values incomplete')
    endif
       
  end subroutine init_pwom_stet_coupling
  
  !=============================================================================
  ! input the pwom grid, thermal e density, and Efield. run stet for iLine
  subroutine get_stet_for_pwom(TimePw,UtPw,iLine,Coord_D,CoordG_D,CoordG2_D,&
       eDensPW_C,eTempPW_C,EfieldPW_C,Ap_I,F107,F107A,IYD,&
       SeDensPW_C, SeFluxPW_C, SeHeatPW_C, IonRatePW_C)
    use ModSeGrid, only: Lshell_I,update_grid,Efield_IC,iLineGlobal_I
    use ModSeBackground,only: mLat_I,mLon_I, gLat1_I,gLat2_I,gLon1_I,gLon2_I,&
         Idate, UT,set_footpoint_locations,fill_thermal_plasma_empirical,&
         plot_background,plot_ephoto_prod,get_neutrals_and_pe_spectrum
    use ModSeState,only: Time
    implicit none
    ! Incomming time from PWOM
    real, intent(in) :: TimePw
    ! Universal Time for background
    real, intent(in) :: UtPW
    ! index of line we are working on
    integer, intent(in) :: iLine
    ! magnetic and geo Lat and Lon coord for base of line in degrees
    real,  intent(in) :: Coord_D(2),CoordG_D(2),CoordG2_D(2)
    ! thermal e dens [/cc], temp [k] and E|| [V/m] from PWOM
    real,  intent(in) :: eDensPW_C(nAltPw), eTempPW_C(nAltPw),EfieldPW_C(nAltPw)
    ! Ap and F107 values from PWOM to set thermosphere in STET
    real,  intent(in) :: Ap_I(7), F107, F107A
    integer,intent(in):: IYD !same as idate, but set in PWOM
    ! SE dens, flux, and heat from STET (interpolated to PWOM grid)
    real,  intent(out)::SeDensPW_C(nAltPw),SeFluxPW_C(nAltPw),SeHeatPW_C(nAltPw)
    ! Ionization rate from STET(interpolated to PWOM grid)
    real,optional,  intent(out):: IonRatePW_C(nAltPw)
    ! named parameters for coordinates
    integer,parameter :: Lat_=1 ,Lon_=2 !named parameters for Coord_ID
    ! Is line open, for now always assume yes, but this could be passed
    logical :: IsOpen=.true.
    ! Conversion between statV/cm in PWOM to Volts/m needed for STET
    real, parameter :: cSTATVperCMtoVperM=30000.0
    !---------------------------------------------------------------------------
    write(*,*) '!!!!!!!!!!!!!!!!!!!'
    write(*,*) 'CALLING get_stet_for_pwom at time=',TimePw
    write(*,*) 'Working on iLine (local,global)=',iLine,iLineGlobal_I(iLine)
    write(*,*) '!!!!!!!!!!!!!!!!!!!'

    ! Set the Time in STET to match the time in PWOM
    Time=TimePw
    
    ! set idate needed for MSIS call
    Idate=IYD
    
    ! set UT in ModBackground to match UT from PWOM
    UT=UtPw

    ! Set the Lat and Lon coordinates and Lshell for field line
    mLat_I(iLine)=Coord_D(Lat_) 
    mLon_I(iLine)=Coord_D(Lon_) 
    gLat1_I(iLine)=CoordG_D(Lat_) 
    gLon1_I(iLine)=CoordG_D(Lon_) 
    gLat2_I(iLine)=CoordG2_D(Lat_) 
    gLon2_I(iLine)=CoordG2_D(Lon_) 
    ! Set Efield to zero everywhere. So only PWOM efield is used.
    Efield_IC(iLine,:)=0.0
    
    !Lshell_I(iLine) = cos(Coord_D(Lat_)*cDegToRad)


    ! From footpoint locations  from PWOM we can set the spatial part 
    ! of the grid which we need to set the background and interpolate the 
    ! PWOM values. The PA ranges will be wrong and the region of exisitance.
    ! this will be corrected after the Efield from PWOM is interpolated 
    ! to the STET grid and update_grid is called again below.
    call update_grid(iLine,Coord_D,DoOnlySpatial=.true.)

    ! Get the neutral atmosphere and photo e production spectrum
    call get_neutrals_and_pe_spectrum(iLine,F107,F107A,AP_I)
    
    ! Fill the background arrays
    call fill_thermal_plasma_empirical(iLine,F107,F107A,Time)
    
    ! plot background before. interp
    call plot_background(iLine,1,time)
    
    ! Interpolate the PW values onto SE. This overwrites the thermal plasma 
    ! in the overlap region of the grids
    call interpolate_pwom_to_stet(iLine,eDensPW_C,eTempPW_C,&
         EfieldPW_C*cSTATVperCMtoVperM)

    write(*,*) 'maxval(Efield_IC(iLine,:))',maxval(Efield_IC(iLine,:))
    write(*,*) 'maxval(EfieldPW_C(:))',maxval(EfieldPW_C(:)*cSTATVperCMtoVperM)
    ! plot background after interp
!    call plot_background(iLine,1,time)
!    call con_stop('')     

    ! Now that we know the footpoint locations we can set or update the SE grid
    call update_grid(iLine,Coord_D)
    
    ! Get a new steady state solution for an open line 
    call stet_run(iLine,IsOpen,.true.)
    
    ! Interpolate the output back to PWOM grid
    if (present(IonRatePW_C)) then
       call interpolate_stet_to_PWOM(iLine,SeDensPW_C,SeFluxPW_C,SeHeatPW_C,&
            IonRatePW_C=IonRatePW_C)
    else
       call interpolate_stet_to_PWOM(iLine,SeDensPW_C,SeFluxPW_C,SeHeatPW_C)
    endif
    
    write(*,*) '!!!!!!!!!!!!!!!!!!!'
    write(*,*) 'FINISH CALLING get_stet_for_pwom at time=',TimePw
    write(*,*) '!!!!!!!!!!!!!!!!!!!'

!    call con_stop('')     
  end subroutine get_stet_for_pwom
  
  
  !=============================================================================
  ! interpolate PWOM thermal electron density and temperature onto STET grid

  subroutine interpolate_stet_to_PWOM(iLine,SeDensPW_C,SeFluxPW_C,SeHeatPW_C,&
       IonRatePW_C)
    use ModSeGrid,      only: FieldLineGrid_IC,nPoint
    use ModSeState,     only: NumberDens_IC,NumberFlux_IC,HeatingRate_IC, &
                              TotalIonizationRate_IC
    use ModInterpolate, only: linear
    implicit none
    ! index of line we are working on
    integer, intent(in) :: iLine
    ! thermal e dens [/cc], temp [k] and E|| [V/m] from PWOM
    real,  intent(out)::SeDensPW_C(nAltPw),SeFluxPW_C(nAltPw),SeHeatPW_C(nAltPw)
    real,  optional, intent(out)::IonRatePW_C(nAltPw)
    
    real :: Coord ! coordinate for interpolation
    integer :: iAlt
    real, parameter :: cEVtoErg=1.60217657e-12 
    !--------------------------------------------------------------------------
    
    !loop over each point on PWOM grid STET outputs to PWOM grid
    ALONG_LINE: do iAlt=1,nAltPw
       Coord= AltPw_C(iAlt)
       SeDensPW_C(iAlt) = &
            linear(NumberDens_IC(iLine,:),1,nPoint,Coord,&
            FieldLineGrid_IC(iLine,:))
       SeFluxPW_C(iAlt) = &
            linear(NumberFlux_IC(iLine,:),1,nPoint,Coord,&
            FieldLineGrid_IC(iLine,:))
       SeHeatPW_C(iAlt) = &
            linear(HeatingRate_IC(iLine,:),1,nPoint,Coord,&
            FieldLineGrid_IC(iLine,:))*cEVtoErg
       if (present(IonRatePW_C)) then
          IonRatePW_C(iAlt) = &
               linear(TotalIonizationRate_IC(iLine,:),1,nPoint,Coord,&
               FieldLineGrid_IC(iLine,:))
       endif
       !write(*,*) AltPw_C(iAlt)/1e5,SeDensPW_C(iAlt),SeFluxPW_C(iAlt),SeHeatPW_C(iAlt)
    enddo ALONG_LINE
  end subroutine interpolate_stet_to_PWOM
  !=============================================================================
  ! interpolate STET density, number flux, and volume heating rate onto PWOMgrid
  ! below PWOM lower boundary use IRI conditions that are already there on 
  ! STET grids. Above PWOM boundary assume density falls with B^2. For now 
  ! use IRI in second hemisphere
  subroutine interpolate_pwom_to_stet(iLine,eDensPW_C,eTempPW_C,EfieldPW_C)
    use ModSeGrid,      only: FieldLineGrid_IC,Efield_IC,nIono,nPoint,nTop
    use ModSeBackground,only:eThermalDensity_IC,eThermalTemp_IC
    use ModInterpolate, only: linear
    implicit none
    ! index of line we are working on
    integer, intent(in) :: iLine
    ! thermal e dens [/cc], temp [k] and E|| [V/m] from PWOM
    real,  intent(in) :: eDensPW_C(nAltPw), eTempPW_C(nAltPw),EfieldPW_C(nAltPw)
    
    real :: NormCoord ! normalized coordinate for interpolation
    integer :: iPoint
    real, parameter :: cKelvinToEV=8.621738e-5
    !--------------------------------------------------------------------------
    
    !assume Efield is zero outside of overlap region
    Efield_IC(iLine,:)=0.0        

    !loop over each point on line and interpolate PWOM values to STET
    ALONG_LINE: do iPoint=1,nTop
       if (FieldLineGrid_IC(iLine,iPoint) < AltPw_C(nAltPw)) then
          if (FieldLineGrid_IC(iLine,iPoint) > AltPw_C(1))then
             NormCoord=&
                  (FieldLineGrid_IC(iLine,iPoint)-AltPw_C(1))/dAltPw
             eThermalDensity_IC(iLine,iPoint)=&
                  linear(eDensPW_C,1,nAltPw,NormCoord)
             eThermalTemp_IC(iLine,iPoint)=&
                  linear(eTempPW_C,1,nAltPw,NormCoord)*cKelvinToEV
             ! special treatment for Efield since we do not set field in iono
             if(iPoint>nIono) then
                Efield_IC(iLine,iPoint)=&
                     linear(EfieldPW_C,1,nAltPw,NormCoord)
             endif
             !set Efield in conjugate hemisphere to match 
             !(ensures good grid setup)
             Efield_IC(iLine,nPoint-iPoint)=-Efield_IC(iLine,iPoint)      
          endif
       else
          !once we are above the PW grid stop the do loop
          exit ALONG_LINE
       endif
    enddo ALONG_LINE
  end subroutine interpolate_pwom_to_stet
  
end Module ModCouplePWOMtoSTET
