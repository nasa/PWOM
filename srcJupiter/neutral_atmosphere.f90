!                                                                      C
!**********************************************************************C
!                                                                      C
!  Neutral atmosphere using JGITM 1-D output
!                                                                      C
!**********************************************************************C
!

SUBROUTINE JupiterAtmos (nAlt,Alt,nH2,nH,nH2O,nCH4,Temp)
      
  use ModInterpolate, ONLY: linear
  use ModIoUnit, ONLY : UnitTmp_
  integer, intent(in) :: nAlt
  real,    intent(in) :: Alt(nAlt)

  real, intent(out) :: nH2(nAlt),nH(nAlt),nH2O(nAlt),nCH4(nAlt)
  real, intent(out) :: Temp(nAlt)

  character (len=50), parameter :: DatafileName='JGITM-1D-atmos.dat'
  integer,            parameter :: nSpecies = 4   ! 4 species in file
  integer,            parameter :: nFileAlt = 10000 ! file altitude grid

  character (len=189) :: line1
  real :: AtmosArray(5+nSpecies,nFileAlt)
  real :: nHe(nAlt)
  integer iAlt
    
  open(UnitTmp_,FILE=DatafileName,STATUS='OLD')
  
  read(UnitTmp_,'(a)') line1
  read(UnitTmp_,*) AtmosArray
  
  close(UnitTmp_)
  
  do iAlt=1,nAlt
     nH2(iAlt) = linear(AtmosArray(5,:), &        ! H2
          1,nFileAlt,Alt(iAlt)/1e5,AtmosArray(1,:))
     nHe(iAlt) = linear(AtmosArray(6,:), &        ! He
          1,nFileAlt,Alt(iAlt)/1e5,AtmosArray(1,:))
     nH(iAlt) = linear(AtmosArray(7,:), &        ! H
          1,nFileAlt,Alt(iAlt)/1e5,AtmosArray(1,:))
     nCH4(iAlt) = linear(AtmosArray(8,:), &        ! CH4
          1,nFileAlt,Alt(iAlt)/1e5,AtmosArray(1,:))
     Temp(iAlt) = linear(AtmosArray(2,:), &        ! Temp
          1,nFileAlt,Alt(iAlt)/1e5,AtmosArray(1,:))
     nH2O(iAlt) = 0.0
  end do

END SUBROUTINE JupiterAtmos
