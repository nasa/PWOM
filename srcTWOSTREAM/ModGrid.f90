Module ModSeGrid
  implicit none
  
  private !except
  ! should the code write in a verbose manner?
  logical, public :: IsVerbose=.false.

  real,public :: Time
 
  !E Grid
  integer,parameter,public :: nEnergy=190
  real, public :: DeltaE_I(nEnergy),EnergyGrid_I(nEnergy)
  real,   public :: EnergyMin, EnergyMax, DeltaE
!  real, allocatable,public  :: KineticEnergy_IIC(:,:,:)

  !Alt Grid
  integer,parameter,public :: nAlt=120
  real,public :: Alt_C(nAlt) !altitude grid in cm

  ! information on how global line number (as opposed to line number on 
  ! given proc)
  integer, public :: iLineGlobal
  
  ! public methods
!  public :: allocate_grid_arrays
  public :: set_egrid
  public :: set_altgrid
  public :: BINNUM

  real, public :: rPlanetCM


contains
  !=============================================================================
  ! Subroutine Based on EGRID from GLOW sets up electron energy grid
  !
  ! Stan Solomon, 1/92
  !
  ! This software is part of the GLOW model.  Use is governed by the Open Source
  ! Academic Research License Agreement contained in the file glowlicense.txt.
  ! For more information see the file glow.txt.
  !
  SUBROUTINE set_egrid
    implicit none

    ! Local:
    Integer iEnergy
    
    DO iEnergy=1,nEnergy
       IF (iEnergy .LE. 21) THEN
          EnergyGrid_I(iEnergy) = 0.5 * REAL(iEnergy)
       ELSE
          EnergyGrid_I(iEnergy) = EXP (0.05 * REAL(iEnergy+26))
       ENDIF
    end do
    
    DeltaE_I(1) = 0.5
    DO iEnergy=2,nEnergy
       DeltaE_I(iEnergy) = EnergyGrid_I(iEnergy)-EnergyGrid_I(iEnergy-1)
    End Do
    
    DO iEnergy=1,nEnergy
       EnergyGrid_I(iEnergy) = EnergyGrid_I(iEnergy) - DeltaE_I(iEnergy)/2.0
    End Do
    
    EnergyMin=EnergyGrid_I(1)
    EnergyMax=EnergyGrid_I(nEnergy)
    
    
  End Subroutine set_egrid

  !=============================================================================
  ! Set Altitude Grid

  SUBROUTINE set_altgrid
    use ModPlanetConst, ONLY: Planet_, rPlanet_I
    implicit none

    ! Local:
    real,parameter :: cKmToCm=1.0e5, cMtoCM = 1.0e2
    Integer ::iAlt
    real :: AltKM_C(nAlt)
    DATA AltKM_C/     80., 81., 82., 83., 84., 85., 86., 87., 88., 89., &
         90., 91., 92., 93., 94., 95., 96., 97., 98., 99., &
         100.,101.,102.,103.,104.,105.,106.,107.,108.,109., &
         110.,111.5,113.,114.5,116.,118.,120.,122.,124.,126., &
         128.,130.,132.,134.,136.,138.,140.,142.,144.,146., &
         148.,150.,153.,156.,159.,162.,165.,168.,172.,176., &
         180.,185.,190.,195.,200.,205.,211.,217.,223.,230., &
         237.,244.,252.,260.,268.,276.,284.,292.,300.,309., &
         318.,327.,336.,345.,355.,365.,375.,385.,395.,406., &
         417.,428.,440.,453.,467.,482.,498.,515.,533.,551., &
         570.,590.,610.,630.,650.,670.,690.,710.,730.,750., &
         770.,790.,810.,830.,850.,870.,890.,910.,930.,950./
    
    Alt_C(:) = AltKM_C * cKmToCm
    
    rPlanetCM=rPlanet_I(Planet_)*cMtoCM

  End Subroutine set_altgrid
  


  !=============================================================================
  ! ------------------------------------------------------------------ ** 
  ! This function finds the energy grid number of the input energy E.
  !  
  INTEGER FUNCTION BINNUM(E)
    
!    use ModSeGrid,      ONLY: FieldLineGrid_IC,nIono,nEnergy, nPoint, &
!         DeltaE_I,EnergyGrid_I, EnergyMin
    
    real, intent(in) :: E
    integer :: N, i
    logical :: IsBinFound
    !---------------------------------------------------------------------------
    IsBinFound=.false.
    i=1
    DO WHILE ((.not.IsBinFound).AND.(i.LE.nEnergy))
       IF (E.LE.EnergyGrid_I(i)+DeltaE_I(i)/2) IsBinFound=.true.
       i=i+1
    END DO
    IF (IsBinFound) THEN
       BINNUM=i-1
    ELSE
       BINNUM=nEnergy+1
    END IF
    IF (E.LE.EnergyMin) BINNUM=0
    RETURN
  END FUNCTION BINNUM


!  subroutine allocate_grid_arrays
!    
!    ! Allocate PA related grid
!
!    if(.not.allocated(DeltaE_I))        allocate(DeltaE_I(nEnergy+1))
!    if(.not.allocated(EnergyGrid_I))    allocate(EnergyGrid_I(nEnergy))
!
!!    if (.not.allocated(iLineGlobal_I))  allocate(iLineGlobal_I(nLine))
!  end subroutine allocate_grid_arrays

 
end Module ModSeGrid
