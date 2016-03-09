! module to hold the photoelectron (well total SE) solution for a given line 
!and associated parameters
module ModPhotoElectron
  real, allocatable :: SeDens_C(:), SeFlux_C(:), SeHeat_C(:)

  !couple time to call update the SE flux
  real :: DtGetSe=120.0
  
  !minimum thermal density of electrons [/cc]
  real :: eThermalDensMin=1.0

  logical :: DoCoupleSTET = .true., UseFeedbackFromSTET=.true.

  !Fixed precipitation to pass to STET
  logical :: UseFixedPrecip = .false.
  real :: PrecipEnergyMin, PrecipEnergyMax, PrecipEnergyMean, PrecipEnergyFlux

  !PolarRain precipitation to pass to STET
  logical :: UsePolarRain = .false.
  real :: PolarRainEMin, PolarRainEMax, PolarRainEMean, PolarRainEFlux

  !Should STET be verbose with output?
  logical :: IsVerboseSTET=.false.

end module ModPhotoElectron
