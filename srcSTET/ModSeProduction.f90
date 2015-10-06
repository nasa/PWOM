Module ModSeProduction
  ! This module contains all of the converted awfulness 
  save
  private !except
  public :: RCOLUM
  public :: espec
  public :: SOLZEN
  public :: SSFLUX
  
  !number of wavelengthincrements
  integer, parameter :: LMAX=59
  !number of excited states
  integer, parameter :: NEI=10
  
  !number of major neutral species
  integer,save :: nNeutral

  ! WAVE1   Array of minimum boundaries for radiation intervals; Ang
  ! WAVE2   Array of maximum boundaries for radiation intervals: Ang
  ! SFLUX   Array of solar flux intensities in the intervals defined 
  real  :: WAVE1(LMAX),WAVE2(LMAX),SFLUX(LMAX)

  ! store the slant column density
  real, allocatable ::ZCOL(:,:)

  ! Store the crosseections
  real :: SIGION(nNeutral,LMAX),SIGABS(nNeutral,LMAX)

  ! Branching ratios
  real :: PROB(nStatesMax,nNeutral,LMAX)
  
  !atomic mass of neutrals
  real,allocatable :: NeutralMassAMU_I(:)

  !Auger e- energies and thresholds
  real,allocatable :: AugEnergy_I(:), AugThreshold_I(:)

  !surface gravity in cm/2^2
  real :: gSurface = 0.0
  
  integer :: nStatesPerSpecies_I(:), nStatesMax 

  ! Store the ionization potentials
  real,allocatable :: TPOT(:,:)
  
contains
  
  !=============================================================================
  ! RCOLUM and required functions and subroutines
  !
  !=============================================================================
  ! Subroutine RCOLUM
  !
  ! Stan Solomon, 1988, 1991
  !
  ! Calculates the column density ZCOL for each species ZMAJ above height
  ! ZZ at zenith angle CHI.  Uses a Chapman function fit [Smith and Smith,
  ! JGR 77, 3592, 1972].  If CHI is less than 90 degrees, column
  ! densities are calculated directly; if CHI is greater than 90 degrees
  ! the column density at grazing height for 90 degrees is calculated and
  ! doubled and the column density above ZZ(J) is subtracted; if CHI is
  ! such that the grazing height is less than the radius of the earth the
  ! column densities are set to 'infinity', i.e., 1.0E30.  Densities
  ! supplied in array ZMAJ are used in the calculation except where
  ! grazing height is below the lowest level specified, in this case
  ! values are interpolated logarithmically from a US Standard Atmosphere
  ! at sea level, the tropopause, the stratopause, and the mesopause.
  ! >> CHI should be in radians
  !
  SUBROUTINE RCOLUM (CHI, ZZ, ZMAJ, TN, IONO)
    ! *** planet specific
!    INCLUDE 'numbers.h'
    PARAMETER (NM=3) 
    PARAMETER (NU=4) !heights in us standard model
    !
    DIMENSION ZZ(IONO), TN(IONO), &
          ZCG(NM), ZUS(NU), TNUS(NU), ZCUS(NM,NU)
    real, allocatable :: ZMAJ(:,:),ZVCD(:,:)
    !
    DATA PI/3.1415926535/, RE/6.37E8/
    ! These are the heights, temperatures and densities for US standard model
    DATA ZUS/0., 1.5E6, 5.E6, 9.E6/, TNUS/288., 217., 271., 187./
    DATA ZCUS/8.00E17, 4.54E24, 1.69E25, &
         8.00E17, 5.46E23, 2.03E24, &
         8.00E17, 3.63E21, 1.35E22, &
         7.80E17, 8.48E18, 3.16E19/
    !
    !
    if (.not.allocated(ZMAJ)) allocate(ZMAJ(nNeutral,IONO))
    if (.not.allocated(ZVCD)) allocate(ZVCD(nNeutral,IONO))

    if (.not.allocated(ZCOL)) allocate(ZCOL(nNeutral,IONO))
    CALL VCD (ZZ, ZMAJ, ZVCD, IONO, nNeutral)
    !
    IF (CHI .GE. 2.) THEN
       do I=1,nNeutral
          do J=1,IONO
             ZCOL(I,J) = 1.0E30
          enddo
       enddo
       RETURN
    ENDIF
    !
    IF (CHI .LE. PI/2.) THEN
       DO I=1,nNeutral
          DO J=1,IONO
             ZCOL(I,J) = ZVCD(I,J) * CHAP(CHI,ZZ(J),TN(J),I)
          enddo
       enddo
    ELSE
       do J=1,IONO
          GHRG=(RE+ZZ(J))*SIN(CHI)
          GHZ=GHRG-RE
          IF (GHZ .LE. 0.) THEN
             do I=1,nNeutral
                ZCOL(I,J) = 1.0E30
             end do
             cycle
          ENDIF
          IF (GHZ .GE. ZZ(1)) THEN
             DO JG=1,J-1
                IF (ZZ(JG) .LE. GHZ .AND. ZZ(JG+1) .GT. GHZ) GOTO 120
             enddo
120          TNG = TN(JG)+(TN(JG+1)-TN(JG))*(GHZ-ZZ(JG))/(ZZ(JG+1)-ZZ(JG))
             do I=1,nNeutral
                ZCG(I) = ZVCD(I,JG) * (ZVCD(I,JG+1) / ZVCD(I,JG)) ** &
                     ((GHZ-ZZ(JG)) / (ZZ(JG+1)-ZZ(JG)))
             enddo
          ELSE
             !Here the grazing altitude is less than the bottom of the model 
             !therefore the values are interpolated from US standard atmosphere 
             !at sea level. Only good for Earth, Not suitable for Jupiter!!!
             ! for Jupiter we have densities down to 0 altitude so just use 
             ! GITM values. For jupiter zz should go to 0 alt.
             do JG=1,3
                IF (ZUS(JG) .LT. GHZ .AND. ZUS(JG+1) .GT. GHZ) GOTO 180
             enddo
180          TNG = TNUS(JG) &
                  + (TNUS(JG+1)-TNUS(JG))*(GHZ-ZUS(JG))/(ZUS(JG+1)-ZUS(JG))
             do I=1,nNeutral
                ZCG(I) = ZCUS(I,JG) * (ZCUS(I,JG+1) / ZCUS(I,JG)) ** &
                     ((GHZ-ZUS(JG)) / (ZUS(JG+1)-ZUS(JG)))
             end do
          ENDIF
          do I=1,nNeutral
             ZCOL(I,J) = 2. * ZCG(I) * CHAP(PI/2.,GHZ,TNG,I) &
                  - ZVCD(I,J) * CHAP(CHI,ZZ(J),TN(J),I)
          end do
       end do
    ENDIF
    !
    RETURN
  END SUBROUTINE RCOLUM
  
  !
  !
  !
  !
  FUNCTION CHAP (CHI, Z, T, I)
    use ModSeGrid,      ONLY: rPlanetCM
    use ModNumConst,    ONLY: cPi
    GR=gSurface*(rPlanetCM/(rPlanetCM+Z))**2
    HN=1.38E-16*T/(NeutralMassAMU_I(I)*1.662E-24*GR)
    HG=(rPlanetCM+Z)/HN
    HF=0.5*HG*(COS(CHI)**2)
    SQHF=SQRT(HF)
    CHAP=SQRT(0.5*cPi*HG)*SPERFC(SQHF)
    RETURN
  END FUNCTION CHAP
  !
  !
  !
  !
  FUNCTION SPERFC(DUMMY)
! error function erfc(y) = 1 - erf(y) 
! from Smith & Smith [1972], eq 12
    IF (DUMMY .LE. 8.) THEN
       SPERFC = (1.0606963+0.55643831*DUMMY) / &
            (1.0619896+1.7245609*DUMMY+DUMMY*DUMMY)
    ELSE
       SPERFC=0.56498823/(0.06651874+DUMMY)
    ENDIF
    RETURN
  END FUNCTION SPERFC
  !
  !
  !
  !
  SUBROUTINE VCD(ZZ,ZMAJ,ZVCD,IONO,nNeutralIn)
    DIMENSION ZZ(IONO), ZMAJ(nNeutralIn,IONO), ZVCD(nNeutralIn,IONO)
    !
    DO I=1,nNeutralIn
       ZVCD(I,IONO) =   ZMAJ(I,IONO) &
            * (ZZ(IONO)-ZZ(IONO-1)) &
            / ALOG(ZMAJ(I,IONO-1)/ZMAJ(I,IONO))
       DO J=IONO-1,1,-1
          RAT = ZMAJ(I,J+1) / ZMAJ(I,J)
          ZVCD(I,J) =   ZVCD(I,J+1) &
               + ZMAJ(I,J) * (ZZ(J)-ZZ(J+1)) / ALOG(RAT) * (1.-RAT)
       enddo
    end DO
    RETURN
  END SUBROUTINE VCD

  !*  Subroutine EPHOTO
!*  Actually, ESPEC, because this allows for any array to go into ZMAJ.
!*  Thus, allowing the use of two conjugate atmospheres stored in
!*  ZMAJ1 and ZMAJ2. Last change: 1/27/94.
!
! Adapted from Banks & Nagy 2-stream input code by Stan Solomon, 6/88
! Modified to handle Auger electrons, Stan Solomon, 7/90
!
! This subroutine calculates photoionization, rates, certain
! photodissociative excitation rates, and the photoelectron production
! spectrum as a function of altitude.  Uses continuously variable energy
! grid.
! For Earth - 3 major neutral species : O, O2, N2
! NO is no longer tracked
! For Jupiter/Saturn - 4 major neutral species : H2, H, CH4, He
!
! Supplied by calling routine:
! WAVE1   wavelength array, upper bound; Angstroms
! WAVE2   wavelength array, lower bound; Angstroms
! SFLUX   solar flux array; photons cm-2 sec-1
! ZZ      altitude array; cm above planet surface
! ZMAJ    density array per species, altitude; cm-3
! ZCOL    slant column density per species, altitude; cm-2
! ENER    energy grid for photoelectrons; eV
! DEL     array of energy grid increments; eV
!
! Calculated by subroutine:
! PESPEC  photoelectron production spectrum for each altitude; cm-3 s-1
! PHOTOI  photoionization rates for state, species, altitude; cm-3 s-1
! PHOTOD  photodissoc./exc. rates for state, species, alt.; cm-3 s-1
!
! Other definitions:
! DSPECT  ionization rate in particular wavelength bin; cm-3 s-1
! TAU     optical depth, dimensionless
! FLUX    solar flux at altitude; cm-2 s-1
! SIGABS  photoabsorption cross sections per species; cm2
! SIGION  photoionization cross sections per species; cm2
! cross section data arrays now moved to init_production
!   Earth -> SIGAO, SIGAO2, SIGAN2, SIGIO, SIGIO2, SIGIN2
!   J/S   -> SIGAH2, SIGAH, SIGACH4, SIGAHe
!            SIGIH2, SIGIH, SIGICH4, SIGIHe
! nStatesPerSpecies     number of states for each species
! TPOT    ionization potentials for each species, state; eV
! PROB    branching ratios for each state, species, and wavelength bin:
!         EARTH:
!           O+ states: 4S, 2Do, 2Po, 4Pe, 2Pe
!           O2+ states: X, a+A, b, dissoc.
!           N2+ states: X, A, B, C, F, dissoc.
!         Jupiter:
!           H2+ states: 
!           H+  states: H+, 
!           CH4+states: H+ + CH3, CH+ +H2 +H, CH2+ + H2,CH3+ +H, CH4+
!           He+ states: 

! PROBO, PROBO2, PROBN2; branching ratio data arrays
! BSO2    yield of O(1S) from dissociation of O2
! EPSIL1  energy loss lower bound for state, species, wavelength; eV
! EPSIL2  energy loss upper bound for state, species, wavelength; eV
! AUGE    Mean energy of Auger electrons for each species; eV
! AUGL    Wavelength threshold for Auger electrons; Angstroms
!
! Array dimensions:
! Elen    number of energetic electron energy bins
! LMAX    number of wavelength intervals for solar flux
! nNeutral    number of major species
! nStatesMax     number of states produced by photoionization/dissociation

  !following not used!
! JMAX    number of altitude levels (actually I am now just using IONO)
! ZNO     density of NO at each altitude; cm-3 
! PHONO   photoionization/dissoc./exc. rates for NO; cm-3 s-1
! SIGNO   NO photoionization xsect at Ly-alpha
! NEX     number of ionized/excited species
! NW      number of airglow emission wavelengths
! NC      number of component production terms for each emission
! NEI     number of states produced by electron impact
! NF      number of available types of auroral fluxes
!
!
  SUBROUTINE ESPEC(ZMAJ,PESPEC,Iono,Ioff)
    !
    !
    
    !INCLUDE 'numbers.h'
    use ModSeGrid,only:nEnergy,del=>DeltaE_I,ener=>EnergyGrid_I,Emin=>EnergyMin
!    PARAMETER (NEX=20)
!    PARAMETER (NW=20)
!    PARAMETER (NC=10)
!    PARAMETER (NF=4)
    !
!    COMMON /CGLOW/ &
!         ZCOL(nNeutral,IONO),WAVE1(LMAX),WAVE2(LMAX),SFLUX(LMAX)
    !
    !      COMMON /CENERGY/ ener(Elen),del(Elen),Emin,Jo
    !
    DIMENSION FLUX(LMAX,IONO), &
         PESPEC(NEnergy,IONO), &
         ZMAJ(nNeutral,IONO), &
         EPSIL1(nStatesMax,nNeutral,LMAX), EPSIL2(nStatesMax,nNeutral,LMAX), &
         PHOTOI(nStatesMax,nNeutral,IONO), PHOTOD(nStatesMax,nNeutral,IONO), &
         BSO2(LMAX), TAU(LMAX), LAUG(nNeutral), &
         EPA(nStatesMax,nStatesMax,nNeutral,LMAX),EPB1(nNeutral,LMAX),EPB2(nNeutral,LMAX)
    !
    SAVE SIGION, SIGABS, PROB, EPSIL1, EPSIL2, EPA, EPB1, EPB2
    !
    DATA  IFIRST/1/, LIMIN/16/ ! No PROBs below L=16
    !
    !

    !
    DATA BSO2/12*0.,.01,.03,.10,.09,.10,.09,.07,.07,.03,.01,37*0./
    !
    !
    DATA C1/12397.7/               ! Converting wavelengths to energie
    !
    !
    !
    !
    ! First time only:  pack photoabsorption, photoioniation cross sections
    ! and convert to cm2; pack branching ratios; calculate energy losses:
    !
    

    IF (IFIRST .EQ. 1) THEN
       IFIRST = 0
          !
       DO  L=1,LMAX
          DO  I=1,nNeutral
    
             IF (WAVE1(L).LE.AugThreshold_I(I)) THEN
                EPB1(I,L)=C1/WAVE1(L)-AugEnergy_I(I)
                EPB2(I,L)=C1/WAVE2(L)-AugEnergy_I(I)
                DO  K1=1,nStatesPerSpecies_I(I)
                   DO  K2=1,nStatesPerSpecies_I(I)
                      EPA(K1,K2,I,L)=AugEnergy_I(I)-TPOT(K1,I)-TPOT(K2,I)
                   end do
                end do
             ELSE
                DO  K=1,nStatesPerSpecies_I(I)
                   EPSIL1(K,I,L)=C1/WAVE1(L)-TPOT(K,I)
                   EPSIL2(K,I,L)=C1/WAVE2(L)-TPOT(K,I)
                enddo
             END IF
          end do
       end do
       write(*,*) 'test2'

!
    ENDIF

    !
    !
    ! Zero arrays:
    !
    DO J=1,IONO
       DO I=1,nNeutral
          DO K=1,nStatesMax
             PHOTOI(K,I,J) = 0.
             PHOTOD(K,I,J) = 0.
          end do
       end do
       DO  M=1,nEnergy
          PESPEC(M,J) = 0.
       end do
    end do


    !
    !
    ! Calculate attenuated solar flux at all altitudes and wavelengths:
    !
    DO  L=1,LMAX
       DO  J=1,Iono
          TAU(L)=0.
          DO  I=1,nNeutral
             TAU(L)=TAU(L)+SIGABS(I,L)*ZCOL(I,J)
          end do
          IF (TAU(L) .LT. 20.) THEN
             FLUX(L,J)=SFLUX(L)*EXP(-TAU(L))
          ELSE
             FLUX(L,J) = 0.0
          ENDIF
          !
          !
          ! Calculate SRC photodissociation of O2, dissociative excitation of
          ! O(1S), photodissociation of N2, and photoionization of NO by solar
          ! Ly-alpha:
          ! *** planet specific, but doesn't seem to be used currently
          IF (WAVE1(L) .LT. 1751 .AND. WAVE2(L) .GT. 1349.) &
               PHOTOD(1,2,J) = PHOTOD(1,2,J)+ZMAJ(2,J)*SIGABS(2,L)*FLUX(L,J)
          PHOTOD(2,2,J) = PHOTOD(2,2,J) + ZMAJ(2,J)*SIGABS(2,L)*FLUX(L,J) &
               * BSO2(L)
          PHOTOD(1,3,J) = PHOTOD(1,3,J) + &
               ZMAJ(3,J)*(SIGABS(3,L)-SIGION(3,L))*FLUX(L,J)
       end do
    end do
    !
    !
    Emax=ENER(nEnergy)+DEL(nEnergy)/2.
    !
    ! Calculate ionization rates and photoelectron production:
    !
    !
    ! Loop over ionizing wavelengths:
    !
    DO  L=LIMIN,LMAX
       !
       !
       ! Loop over species:
       !
       DO  I=1,nNeutral
          !
          !
          ! Choose between ionization possibilities
          !
          IF (WAVE1(L).GT.AugThreshold_I(I)) THEN      ! No Auger electron production
             !
             ! Loop over altitude:
             !
             DO  J=1,Iono
                J1=ABS(Ioff-J)
                !
                !
                ! Loop over states:
                !
                DO  K=1,nStatesPerSpecies_I(I)
                   !
                   E1= EPSIL1(K,I,L)
                   E2= EPSIL2(K,I,L)
                   Y=E2-E1
                   IF (E2 .LT. Emin) cycle
                   IF (E1 .GT. Emax) cycle
                   !
                   !
                   ! Calculate ionization rates:
                   !
                   DSPECT = ZMAJ(I,J)*SIGION(I,L)*FLUX(L,J)*PROB(K,I,L)
                   PHOTOI(K,I,J) = PHOTOI(K,I,J) + DSPECT
                   !
                   !
                   ! Find box numbers M1, M2 corresponding to energies E1, E2:
                   !
                   CALL BOXNUM (E1,E2,M1,M2,R1,R2,Emax)
                   !
                   !
                   ! Fill the boxes from M1 to M2:
                   !
                   DO  N=M1,M2
                      IF (M1 .EQ. M2) THEN
                         FAC = 1.
                         IF ((M1.EQ.1).AND.(E1.LT.Emin)) FAC = R2 / Y
                         IF ((M1.EQ.nEnergy).AND.(E2.GT.Emax)) FAC = R1 / Y
                      ELSE
                         IF (N .EQ. M1) THEN
                            FAC = R1 / Y
                         ELSE
                            IF (N .EQ. M2) THEN
                               FAC = R2 / Y
                            ELSE
                               FAC = DEL(N) / Y
                            ENDIF
                         ENDIF
                      ENDIF
                      PESPEC(N,J) = PESPEC(N,J) + DSPECT * FAC
                   enddo
                   !
                enddo   ! End of ion state loop
                !
             enddo            ! End of altitude loop
             !
             !
          ELSE            ! Begin second branch (double ionization)
             !
             ! Loop over altitude:
             !
             DO  J=1,Iono
                J1=ABS(Ioff-J)
                !
                !
                ! Calculate the electron with energy AugEnergy_I-TPOT(K1)-TPOT(K2)
                !
                DO  K1=1,nStatesPerSpecies_I(I)! Excited state of ion for initial e-
                   !
                   DO  K2=1,nStatesPerSpecies_I(I) ! Excited state of final ion
                      !
                      E1= EPA(K1,K2,I,L)
                      E2= E1
                      IF (E1.LT.Emin .OR. E1.GT.Emax) cycle
                      DSPECT = &
                           ZMAJ(I,J)*SIGION(I,L)*FLUX(L,J)&
                           *PROB(K1,I,L)*PROB(K2,I,LAUG(I))
                      PHOTOI(K1,I,J) = PHOTOI(K1,I,J) + DSPECT      ! Technically, it's
                      PHOTOI(K2,I,J) = PHOTOI(K2,I,J) + DSPECT      ! double ionization
                      CALL BOXNUM (E1,E2,M1,M2,R1,R2,Emax)       ! not two single ions
                      PESPEC(M1,J) = PESPEC(M1,J) + DSPECT
                      !
                   enddo            ! End of ion states loops
                enddo
                !
                ! Calculate the electron with energy C1/WAVE-AugEnergy_I
                !
                E1= EPB1(I,L)
                E2= EPB2(I,L)
                Y=E2-E1
                IF (.not.(E2 .LT. Emin .or. E1 .GT. Emax)) then
                   
                   !
                   !
                   ! Calculate ionization rates:
                   !
                   DSPECT = ZMAJ(I,J)*SIGION(I,L)*FLUX(L,J)      ! Already added into
                   !                                    ! PHOTOI above
                   !
                   ! Find box numbers M1, M2 corresponding to energies E1, E2:
                   !
                   CALL BOXNUM (E1,E2,M1,M2,R1,R2,Emax)
                   
                   !
                   ! Fill the boxes from M1 to M2:
                   !
                   DO N=M1,M2
                      IF (M1 .EQ. M2) THEN
                         FAC = 1.
                         IF ((M1.EQ.1).AND.(E1.LT.Emin)) FAC = R2 / Y
                         IF ((M1.EQ.nEnergy).AND.(E2.GT.Emax)) FAC = R1 / Y
                      ELSE
                         IF (N .EQ. M1) THEN
                            FAC = R1 / Y
                         ELSE
                            IF (N .EQ. M2) THEN
                               FAC = R2 / Y
                            ELSE
                               FAC = DEL(N) / Y
                            ENDIF
                         ENDIF
                      ENDIF
                      PESPEC(N,J) = PESPEC(N,J) + DSPECT * FAC
                   end do
                   
                   
                endif            ! "Out of our E range" skip
                !
             enddo            ! End of altitude loop
             !
          END IF            ! End of second branch
          !
       end DO            ! End of species loop
!
    end DO            ! End of wavelength loop
!
!
    RETURN
    !
  END SUBROUTINE ESPEC
  !
  !
  !
  !
  SUBROUTINE BOXNUM (E1,E2,M1,M2,R1,R2,Emax)
    !
    ! This subroutine finds the box numbers corresponding to
    ! energies E1 and E2, and calls them M1 and M2
    !
    ! R1 is the upper edge of the lower box, R2 is the lower edge of the
    ! upper box.
    !
    use ModSeGrid,only:nEnergy,del=>DeltaE_I,ener=>EnergyGrid_I,Emin=>EnergyMin
    !INCLUDE 'numbers.h'
    !COMMON /CENERGY/ ENER(nEnergy),DEL(nEnergy),Emin,nEnergy
    !
    I=0
    M1=0
    DO  WHILE (M1.EQ.0)
       I=I+1
       IF (E1.LT.ENER(I)+DEL(I)/2.) M1=I
    end DO
    R1=ENER(M1)+DEL(M1)/2.-E1
    IF (E1.LT.Emin) R1=DEL(1)
    !
    I=0
    M2=0
    DO  WHILE (M2.EQ.0)
       I=I+1
       IF (E2.LT.ENER(I)+DEL(I)/2.) M2=I
       IF (I.EQ.nEnergy) M2=I
    end do
    R2=E2-ENER(M2)+DEL(M2)/2.
    IF (E2.GT.Emax) R2=DEL(nEnergy)
    
    RETURN
  END SUBROUTINE BOXNUM

  !=============================================================================
  ! Subroutine SOLZEN
  !
  ! Stan Solomon, 1988
  !
  ! Returns Solar Zenith Angle SZA in degrees for specified date in form
  ! yyddd, universal time in seconds, geographic latitude and longitude
  ! in degrees.
  !
  ! note for jupiter or another planet will need planet specific calculation
  SUBROUTINE SOLZEN (IDATE, UT, GLAT, GLONG, SZA)
    !
    DATA PI/3.1415926536/
    !
    RLAT = GLAT * PI/180.
    RLONG = GLONG * PI/180.
    CALL SUNCOR (IDATE, UT, SDEC, SRASN, GST)
    RH = SRASN - (GST+RLONG)
    COSSZA = SIN(SDEC)*SIN(RLAT) + COS(SDEC)*COS(RLAT)*COS(RH)
    SZA = ACOS(COSSZA) * 180./PI
    RETURN
  END SUBROUTINE SOLZEN
  !
  !
  !
  !
  ! Subroutine SUNCOR returns the declination SDEC and right ascension
  ! SRASN of the sun in GEI coordinates, radians, for a given date IDATE
  ! in yyddd format and universal time UT in seconds.  Greenwich Sidereal
  ! Time GST in radians is also returned.  Reference:  C.T. Russell,
  ! Geophysical Coordinate Transforms.
  !
  SUBROUTINE SUNCOR (IDATE, UT, SDEC, SRASN, GST)
    DATA PI/3.1415926536/
    !
    FDAY=UT/86400.
    IYR=IDATE/1000
    IDAY=IDATE-IYR*1000
    DJ=365*IYR+(IYR-1)/4+IDAY+FDAY-0.5
    T=DJ/36525.
    VL=AMOD(279.696678+.9856473354*DJ,360.)
    GST=AMOD(279.696678+.9856473354*DJ+360.*FDAY+180.,360.) * PI/180.
    G=AMOD(358.475845+.985600267*DJ,360.) * PI/180.
    SLONG=VL+(1.91946-.004789*T)*SIN(G)+.020094*SIN(2.*G)
    OBLIQ=(23.45229-0.0130125*T) *PI/180.
    SLP=(SLONG-.005686) * PI/180.
    SIND=SIN(OBLIQ)*SIN(SLP)
    COSD=SQRT(1.-SIND**2)
    SDEC=ATAN(SIND/COSD)
    SRASN=3.14159-ATAN2(1./TAN(OBLIQ)*SIND/COSD,-COS(SLP)/COSD)
    RETURN
  END SUBROUTINE SUNCOR

  ! Subroutine SSFLUX
  !
  ! Stan Solomon, 2/92
  !
  ! Subroutine SSFLUX calculates the solar EUV and FUV flux in the range
  ! 1 to 1750 Angstroms for a specified level of solar activity.  The
  ! calling routine supplies a scaling switch ISCALE, the daily 10.7 cm
  ! flux F107, its 81-day centered average F107A, and, optionally,
  ! the H Lyman-Beta 1026A ratio to its solar minimum value HLYBR, the
  ! Fe XVI 335A ratio to its solar minimum value FEXVIR, the H Lyman-alpha
  ! 1216A flux HLYA, the He I 10830A equivalent width HEIEW, and an XUV
  ! enhancement factor XUVFAC.  Any optional calling parameters not used
  ! should be set to zero.  The subroutine returns the longwave boundary
  ! WAVE1 and shortwave boundary WAVE2 of the wavelenth bins, and the
  ! solar flux in each bin SFLUX.
  !    One of five methods is used, depending on the value of ISCALE,
  ! as follows:
  !   If ISCALE=0 the flux is scaled using parameterization methods based
  ! on F107 and F107A.  For ionizing EUV, Hinteregger's contrast ratio
  ! method (Hinteregger et al., GRL, 8, 1147, 1981) is used, based on the
  ! Torr and Torr (JGR 90, 6675, 1985) bin structure for reference
  ! spectrum SC#21REFW.  If the H Lyman-Beta (1026A) or Fe XVI (335A)
  ! enhancement ratios are provided (>0) as calling arguments, they are
  ! used to scale the EUV spectrum.  Otherwise, enhancement ratios for
  ! H Ly-B and Fe XVI are calculated from F107 and F107A using
  ! Hinteregger's formula, employing coefficients which reduce to the
  ! reference values at F107=67.6, F107A=71.5.  The 'best fit'
  ! coefficients are not used as they produce some negative values at low
  ! solar activity, but remain in a 'commented out' data statement for
  ! reference.  The EUV spectrum is then scaled from these modeled ratios.
  ! Scaling factors were calculated from contrast ratios in the SC#21REFW
  ! data file. For FUV from 1050A-1200A, 50A interval averaging from
  ! SC#21REFW is also used.  Note that O2 cross sections used by EPHOTO
  ! in the 1050-1350A region are approximate since 50A bins are not really
  ! adequate for O2 absorption and dissociation calculations in this
  ! region.  For H Lyman-alpha at 1216A, which is treated separately as
  ! an individual line, the correlation relationship with F10.7 derived
  ! from SME data by Barth et al. (GRL, 17, 571, 1990) for the ascending
  ! phase of solar cycle 22 is used, unless a value for HLYA (>0) is
  ! provided by the calling program, in which case this value is
  ! subsituted into the spectrum.  For 1200-1350A and for the
  ! Schumann-Runge continuum (1350-1750A), SME and LASP rocket flight
  ! data  (Rottman, JGR 86, 6697, 1981; Adv. Space Res. 8, (7)53, 1988;
  ! Mount and Rottman, JGR 88, 5403, 1983; JGR 90, 13031, 1985) are scaled
  ! by linear interpolation using F10.7.  The solar minimum spectrum and
  ! scaling factors for 1200-1750A should be considered preliminary at
  ! this time.  50A bins are also used, following Torr et al. (JGR 80,
  ! 6063, 1980).
  !   If ISCALE=1 linear interpolation between high and low activity
  ! spectra is performed, based on F107 alone, and assuming that the low
  ! activity spectrum corresponds to F107=68 and the high activitiy
  ! spectrum to F107=243.  The Hinteregger SC#21REFW and F79050N spectra
  ! as binned by Torr and Torr are used for ionizing EUV.  The SC#21REFW
  ! and F79050N spectra are also binned into 50A intervals from
  ! 1050-1200A.  From 1200-1750A, linear interpolation on F10.7 of SME
  ! data as above is used regardless of the value of ISCALE.
  !   If ISCALE=2, the EUV flux (16-1050A) is scaled using the EUV91 model
  ! (Tobiska, "EUV91: Solar EUV Flux Model", J. Atm. Terr. Phys., 31, in
  ! press, 1991, as updated by Tobiska et al., AGU, fall 1991) This model
  ! uses the H Lyman-alpha flux HLYA and HeI 10830A equivalent width HEIEW
  ! to scale chromospheric fluxes, and F107 and F107A to scale coronal
  ! fluxes.  If only one of HLYA or HEIEW is provided (>0), the subroutine
  ! calculates the other from a linear correlation relationship.  If
  ! neither is supplied, the subroutine calculates HLYA from F107 using
  ! the SME-derived formula as above, and uses that in turn to calculate
  ! HEIEW.  Bins not covered by the EUV91 model are scaled using linear
  ! interpolation on F10.7 as above.  Since the EUV91 model uses a single
  ! bin from 18-31A instead of 16-23A and 23-32A, these two are adjusted
  ! accordingly.
  !    If ISCALE=3, the Woods & Rottman rocket spectrum measured on
  ! 10 November, 1988 is employed from 300-1100A, with the EUV91 model
  ! spectrum for that day from 16-300A and linear interpolation from
  ! 1-16A.
  !    If ISCALE=4, the Woods & Rottman rocket spectrum measured on
  ! 20 June, 1989 is employed from 400-1100A, with the EUV91 model
  ! spectrum for that day from 16-400A and linear interpolation from
  ! 1-16A.
  !    "XUV" fluxes between 16A and 250A may be adjusted by a constant
  ! factor (e.g., by a factor of 2.0, cf. Richards and Torr, 1984).  The
  ! XUVFAC parameter specifies this quantity.
  !    None of the above models extends shortward of 18A, so from 1-18 A
  ! an amalgam of sources are used to derive an estimated flux, e.g.,
  ! DeJager, in Astronomical Observations from Space Vehicles, Steinberg,
  ! ed., 1964; Smith & Gottlieb, SSR 16, 771, 1974; Manson, in The Solar
  ! Output and its Variation, White, ed., 1977; Kreplin et al, ibid;
  ! Horan & Kreplin, Solar Physics 74, 265, 1981; Wagner, Adv. Space Res.
  ! 8, (7)67, 1988.
  !
  ! Modification history:
  !   Stan Solomon, 12/88  Basic Hinteregger EUV, approx. SME FUV
  !   C.S. Gaskill,  7/89  Added early Tobiska model
  !   Stan Solomon,  8/89  Corrections to above
  !   Stan Solomon,  1/90  Tobiska SERF2; added W & R spectra
  !   Stan Solomon,  6/91  Tobiska EUV 91; Hntggr Ly-B, Fe XVI scaling
  !   Stan Solomon,  2/92  Updated Tobiska EUV91; corrected SME FUV
  !
  ! Arguments:
  ! ISCALE   =0 for Hinteregger contrast ratio method
  !          =1 for Hinteregger linear interpolation
  !          =2 for Tobiska EUV91 model
  !          =3 for Woods & Rottman 10 Nov. 1988 measurement
  !          =4 for Woods & Rottman 20 Jun. 1989 measurement
  ! F107     daily 10.7 cm flux (1.E-22 W m-2 Hz-1)
  ! F107A    81-day centered average 10.7 cm flux
  ! HLYBR    ratio of H Ly-b 1026A flux to solar minimum value (optional)
  ! FEXVIR   ratio of Fe XVI 335A flux to solar minimum value (optional)
  ! HLYA     H Lyman-alpha flux (photons cm-2 s-1) (optional)
  ! HEIEW    He I 10830A equivalent width, (milliAngstroms) (optional)
  ! XUVFAC   factor for scaling flux 16-250A (optional)
  ! WAVE1    longwave bound of spectral intervals (Angstroms)
  ! WAVE2    shortwave bound of intervals (= WAVE1 for indiv. lines)
  ! SFLUX    scaled solar flux returned by subroutine (photons cm-2 s-1)
  !
  ! Other definitions:
  ! LMAX     dimension of flux and scaling arrays, currently = 59
  ! WAVEL    = WAVE1
  ! WAVES    = WAVE2
  ! RFLUX    low solar activity reference flux, from SC#21REFW
  ! XFLUX    high solar activity flux, from F79050
  ! SCALE1   scaling factors for H LyB-keyed chromospheric emissions
  ! SCALE2   scaling factors for FeXVI-keyed coronal emissions
  ! TCHR0    intercept for chromospheric fluxes, EUV91 model
  ! TCHR1    HLYA slope for    "
  ! TCHR2    HEIEW slope for    "
  ! TCOR0    intercept for coronal fluxes, EUV91 model
  ! TCOR1    F107 slope for  "
  ! TCOR2    F107A slope for "
  ! WAR1     Woods & Rottman 10 Nov. 1988 spectrum
  ! WAR2     Woods & Rottman 20 Jun. 1989 spectrum
  ! B1       fit coefficients for H LyB
  ! B2       fit coefficients for FeXVI
  ! R1       enhancement ratio for H LyB
  ! R2       enhancement ratio for FeXVI
  ! HLYMOD   model H Ly-a, = HLYA if provided > 0
  ! HEIMOD   model He I 10830 equivalent width converted to an H Ly-a flux
  !
  !
  SUBROUTINE SSFLUX (ISCALE, F107, F107A, HLYBR, FEXVIR, HLYA, &
       HEIEW, XUVFAC)
    !
!    PARAMETER (LMAX=59)
    !
    DIMENSION WAVEL(LMAX), WAVES(LMAX), RFLUX(LMAX), XFLUX(LMAX), &
         SCALE1(LMAX), SCALE2(LMAX), &
         TCHR0(LMAX), TCHR1(LMAX), TCHR2(LMAX), &
         TCOR0(LMAX), TCOR1(LMAX), TCOR2(LMAX), &
         WAR1(LMAX), WAR2(LMAX),FISM(LMAX), &
         B1(3), B2(3)
    !
    ! regression coefficients which reduce to solar min. spectrum:
    DATA B1/1.0, 0.0138, 0.005/, B2/1.0, 0.59425, 0.3811/
    !
    ! 'best fit' regression coefficients, commented out, for reference:
    !     DATA B1/1.31, 0.01106, 0.00492/, B2/-6.618, 0.66159, 0.38319/
    !
    !
    DATA WAVEL/ 1750.00, 1700.00, 1650.00, 1600.00, 1550.00, 1500.00, &
         1450.00, 1400.00, 1350.00, 1300.00, 1250.00, 1215.67, &
         1200.00, 1150.00, 1100.00, 1050.00, 1031.91, 1025.72, &
         1000.00,  977.02,  950.00,  900.00,  850.00,  800.00, &
         789.36,  770.41,  765.15,  750.00,  703.31,  700.00, &
         650.00,  629.73,  609.76,  600.00,  584.33,  554.37, &
         550.00,  500.00,  465.22,  450.00,  400.00,  368.07, &
         350.00,  303.78,  303.31,  300.00,  284.15,  256.30, &
         250.00,  200.00,  150.00,  100.00,   50.00,   32.00, &
         23.00,   16.00,    8.00,    4.00,    2.00/
    DATA WAVES/ 1700.00, 1650.00, 1600.00, 1550.00, 1500.00, 1450.00, &
         1400.00, 1350.00, 1300.00, 1250.00, 1200.00, 1215.67, &
         1150.00, 1100.00, 1050.00, 1000.00, 1031.91, 1025.72, &
         950.00,  977.02,  900.00,  850.00,  800.00,  750.00, &
         789.36,  770.41,  765.15,  700.00,  703.31,  650.00, &
         600.00,  629.73,  609.76,  550.00,  584.33,  554.37, &
         500.00,  450.00,  465.22,  400.00,  350.00,  368.07, &
         300.00,  303.78,  303.31,  250.00,  284.15,  256.30, &
         200.00,  150.00,  100.00,   50.00,   32.00,   23.00, &
         16.00,    8.00,    4.00,    2.00,    1.00/
    DATA RFLUX/  322.00,  168.00,   95.00,   62.00,   44.00,   25.00, &
         16.90,   11.80,   19.50,    4.10,   11.10,  249.00, &
         2.78,    0.70,    3.07,    3.64,    3.18,    4.38, &
         1.78,    5.96,    4.22,    4.43,    1.93,    0.87, &
         0.79,    0.24,    0.20,    0.17,    0.39,    0.22, &
         0.17,    1.50,    0.45,    0.48,    1.58,    0.80, &
         0.51,    0.31,    0.18,    0.39,    0.21,    0.74, &
         0.87,    6.00,    0.24,    0.84,    0.10,    0.27, &
         0.92,    1.84,    0.13,    0.38,  0.0215,  0.0067, &
         1.E-3,   2.E-3,   1.E-5,   5.E-8,   1.E-10/
    DATA XFLUX/  354.00,  191.00,  110.00,   76.00,   55.00,   28.00, &
         19.60,   14.30,   25.30,    5.00,   17.20,  401.00, &
         6.26,    1.51,    6.11,    8.66,    9.04,   13.12, &
         4.42,   13.18,   12.03,   13.29,    5.01,    2.18, &
         1.59,    0.67,    0.43,    0.43,    0.72,    0.46, &
         0.48,    3.02,    1.46,    1.02,    4.86,    1.59, &
         1.57,    1.67,    0.36,    0.99,    2.20,    1.39, &
         5.63,   11.28,    2.50,    4.14,    3.16,    0.59, &
         3.70,    4.85,    0.34,    1.15,    0.18,    0.08, &
         2.5E-2,   5.E-2,   8.E-4,   3.E-5,   5.E-7/
    DATA SCALE1/    0.0,     0.0,     0.0,     0.0,     0.0,     0.0, &
         0.0,     0.0,     0.0,     0.0,     0.0,     0.0, &
         1692.09,  405.95, 1516.20, 2731.70, 3314.57, 4375.00, &
         1316.91, 3621.91, 3908.56, 4432.54, 1541.21,  531.73, &
         364.83,    0.00,  116.00,  129.41,  162.48,   94.07, &
         41.29,  709.50,    0.00,  268.47, 1561.05,  367.64, &
         290.06,  184.36,    0.00,   86.15,    7.50,    0.00, &
         0.00, 2220.00,    0.00,   61.00,    0.00,   86.95, &
         206.00,  135.89,   60.35,  157.12,    7.06,    0.75, &
         0.00,    0.00,    0.00,    0.00,    0.00/
    DATA SCALE2/   0.00,    0.00,    0.00,    0.00,    0.00,    0.00, &
         0.00,    0.00,    0.00,    0.00,    0.00,    0.00, &
         0.00,    0.00,    0.00,    0.00,    0.00,    0.00, &
         0.00,    0.00,    0.00,    0.00,    0.00,    0.00, &
         0.00,    5.34,    0.00,    0.00,    0.00,    0.54, &
         3.30,    0.00,   12.60,    0.00,    0.00,    0.00, &
         5.34,   11.63,    2.28,    5.56,   24.93,    8.16, &
         60.69,    0.00,   28.20,   45.90,   40.80,    1.27, &
         35.47,   42.80,    1.12,    6.19,    1.26,    0.69, &
         0.23,    0.46,  7.6E-3,  2.9E-4,  4.8E-6/
    DATA TCHR0/ &
         0.0,       0.0,       0.0,       0.0,       0.0,       0.0, &
         0.0,       0.0,       0.0,       0.0,       0.0,       0.0, &
         0.0,       0.0,       0.0,-4.290E+00,-5.709E+00,-8.493E+00, &
         -1.161E+00,-3.429E+00,-5.464E+00,-6.502E+00,-1.912E+00,-4.034E-01, &
         -1.448E-01, 0.000E+00,-9.702E-02,-6.591E-02,-2.338E-02,-1.273E-01, &
         -2.406E-01,-3.351E-01, 0.000E+00,-1.465E+00,-2.405E+00,-7.975E-02, &
         -4.197E-01,-1.971E-01, 0.000E+00,-5.895E-02,-5.815E-03, 0.000E+00, &
         0.000E+00, 2.138E-01, 0.000E+00,-7.713E-02, 0.000E+00,-3.035E-02, &
         -2.039E-01,-1.749E-01,-1.041E-01,-2.638E-01,-1.094E-02, 0.000E+00, &
         0.0,       0.0,       0.0,       0.0,       0.0/
    DATA TCHR1/ &
         0.0,       0.0,       0.0,       0.0,       0.0,       0.0, &
         0.0,       0.0,       0.0,       0.0,       0.0,       0.0, &
         0.0,       0.0,       0.0,-3.023E-13,-3.745E-13,-5.385E-13, &
         -1.211E-13,-3.868E-13,-3.646E-13,-4.125E-13,-1.527E-13,-4.753E-14, &
         -3.411E-14, 0.000E+00,-1.190E-14,-1.034E-14,-1.343E-14,-1.539E-14, &
         -5.174E-14,-6.934E-14, 0.000E+00,-1.215E-13,-1.537E-13,-2.024E-14, &
         -4.596E-14,-1.562E-14, 0.000E+00,-1.221E-14,-1.123E-15, 0.000E+00, &
         0.000E+00,-2.263E-13, 0.000E+00,-1.508E-14, 0.000E+00,-1.744E-14, &
         -2.100E-14,-1.805E-14,-8.224E-15,-1.919E-14,-7.944E-16, 0.000E+00, &
         0.0,       0.0,       0.0,       0.0,       0.0/
    DATA TCHR2/ &
         0.0,       0.0,       0.0,       0.0,       0.0,       0.0, &
         0.0,       0.0,       0.0,       0.0,       0.0,       0.0, &
         0.0,       0.0,       0.0, 3.275E-11, 4.057E-11, 6.160E-11, &
         1.312E-11, 4.189E-11, 4.167E-11, 4.716E-11, 1.654E-11, 5.150E-12, &
         3.901E-12, 0.000E+00, 1.289E-12, 1.120E-12, 1.455E-12, 1.667E-12, &
         5.604E-12, 7.931E-12, 0.000E+00, 1.317E-11, 1.757E-11, 2.194E-12, &
         4.978E-12, 1.693E-12, 0.000E+00, 1.324E-12, 1.285E-13, 0.000E+00, &
         0.000E+00, 2.586E-11, 0.000E+00, 1.724E-12, 0.000E+00, 1.889E-12, &
         2.400E-12, 2.063E-12, 8.911E-13, 2.193E-12, 9.090E-14, 0.000E+00, &
         0.0,       0.0,       0.0,       0.0,       0.0/
    DATA TCOR0/ &
         0.0,       0.0,       0.0,       0.0,       0.0,       0.0, &
         0.0,       0.0,       0.0,       0.0,       0.0,       0.0, &
         0.0,       0.0,       0.0, 0.000E+00, 0.000E+00, 0.000E+00, &
         0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,-6.060E-02, &
         0.000E+00,-3.399E-02, 0.000E+00, 0.000E+00, 0.000E+00, 4.866E-02, &
         -1.762E-01, 0.000E+00,-2.412E-01, 0.000E+00, 0.000E+00, 0.000E+00, &
         -4.743E-01,-9.713E-01, 5.891E-02,-1.263E-01,-1.246E+00, 2.870E-01, &
         -4.659E+00, 0.000E+00,-1.058E+00,-3.821E+00,-1.874E+00, 0.000E+00, &
         -1.896E+00,-8.505E-01,-2.101E-04,-2.012E-01,-6.097E-02,-2.925E-02, &
         -4.875E-03,       0.0,       0.0,       0.0,       0.0/
    DATA TCOR1/ &
         0.0,       0.0,       0.0,       0.0,       0.0,       0.0, &
         0.0,       0.0,       0.0,       0.0,       0.0,       0.0, &
         0.0,       0.0,       0.0, 0.000E+00, 0.000E+00, 0.000E+00, &
         0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 2.877E-03, &
         0.000E+00, 1.760E-03, 0.000E+00, 0.000E+00, 0.000E+00, 3.313E-04, &
         3.643E-03, 0.000E+00, 5.225E-03, 0.000E+00, 0.000E+00, 0.000E+00, &
         4.085E-03, 1.088E-02, 8.447E-04, 3.237E-03, 1.907E-02, 2.796E-03, &
         4.460E-02, 0.000E+00, 1.007E-02, 3.481E-02, 1.604E-02, 0.000E+00, &
         2.029E-02, 2.160E-02, 6.342E-04, 3.594E-03, 5.503E-04, 2.687E-04, &
         4.479E-05,       0.0,       0.0,       0.0,       0.0/
    DATA TCOR2/ &
         0.0,       0.0,       0.0,       0.0,       0.0,       0.0, &
         0.0,       0.0,       0.0,       0.0,       0.0,       0.0, &
         0.0,       0.0,       0.0, 0.000E+00, 0.000E+00, 0.000E+00, &
         0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.846E-03, &
         0.000E+00, 1.127E-03, 0.000E+00, 0.000E+00, 0.000E+00, 1.891E-04, &
         2.326E-03, 0.000E+00, 2.801E-03, 0.000E+00, 0.000E+00, 0.000E+00, &
         2.446E-03, 7.121E-03, 5.204E-04, 1.983E-03, 1.204E-02, 1.721E-03, &
         2.911E-02, 0.000E+00, 7.177E-03, 2.272E-02, 9.436E-03, 0.000E+00, &
         1.316E-02, 1.398E-02, 4.098E-04, 2.328E-03, 3.574E-04, 1.745E-04, &
         2.909E-05,       0.0,       0.0,       0.0,       0.0/
    DATA WAR1/     0.00,    0.00,    0.00,    0.00,    0.00,    0.00, &
         0.00,    0.00,    0.00,    0.00,    0.00,    0.00, &
         0.00,    0.00,    3.80,    6.25,    4.93,    6.06, &
         2.70,    7.07,    8.62,    9.60,    4.54,    2.37, &
         0.82,    0.33,    0.24,    0.67,    0.28,    0.55, &
         1.56,    1.11,    0.77,    1.32,    1.71,    0.44, &
         1.11,    0.95,    0.39,    0.81,    2.00,    1.49, &
         6.81,    5.07,    1.63,    5.62,    2.08,    0.59, &
         3.89,    5.19,    0.35,    1.18,    0.099,   0.04, &
         0.007,   0.00,    0.00,    0.00,    0.00/
    DATA WAR2/     0.00,    0.00,    0.00,    0.00,    0.00,    0.00, &
         0.00,    0.00,    0.00,    0.00,    0.00,    0.00, &
         0.00,    0.00,   20.80,   17.90,    9.30,   14.30, &
         6.90,   12.00,   15.60,   18.60,   10.10,    4.30, &
         12.40,    8.00,    3.60,    1.80,    0.50,    1.40, &
         3.90,    2.60,    1.60,    3.40,    4.10,    0.70, &
         4.30,    4.30,    3.80,    2.60,    6.08,    1.35, &
         12.60,    9.78,    2.96,   10.20,    4.11,    6.68, &
         6.62,    8.07,    0.47,    1.73,    0.17,    0.075, &
         0.012,   0.00,    0.00,    0.00,    0.00/
    
    ! FISM spectrum for 46th day of 1997
    DATA FISM/    305.263, 176.664, 93.7853, 64.2060,  43.7263, 24.415, &
         16.4750, 12.9340,24.5531,  7.09041, 37.1989,336.710, &
         11.5175,  3.33318,2.55812, 3.60368,  1.51270, 3.032, &
         2.32707, 4.82513,5.20694, 6.12267,  3.70282, 2.498, &
         0.482684,0.273679,0.0777309,0.792697,0.275682,0.54, &
         0.492563,1.92352, 0.601734,0.684815,1.40116,0.6644, &
         1.00781, 0.925885,0.287596,0.606941,0.663962,0.472, &
         1.33970, 5.87943, 0.308720,0.941772,0.0145345,0.40, &
         1.98800, 3.55109, 0.354092,0.399132,0.0591637,0.00, &
         0.0188158,0.00341563,3.91238e-05,6.08654e-08, 2.52/
    
    
    !
    ! Linear Interpolation between SC#21REFW and F79050:
    !
    FRAT = (F107-68.) / (243.-68.)
    DO  L=1,LMAX
       SFLUX(L) = RFLUX(L) + (XFLUX(L)-RFLUX(L)) * FRAT
    end do
    !
    ! Hinteregger contrast ratio method:
    !
    IF (ISCALE .EQ. 0) THEN
       IF (HLYBR .GT. 0.001) THEN
          R1 = HLYBR
       ELSE
          R1 =  B1(1) + B1(2)*(F107A-71.5) + B1(3)*(F107-F107A+3.9)
       ENDIF
       IF (FEXVIR .GT. 0.001) THEN
          R2 = FEXVIR
       ELSE
          R2 =  B2(1) + B2(2)*(F107A-71.5) + B2(3)*(F107-F107A+3.9)
       ENDIF
       DO  L=13,LMAX
          SFLUX(L) = (RFLUX(L) + ((R1-1.)*SCALE1(L) &
               +  (R2-1.)*SCALE2(L)) / 1000.)
       end do
    ENDIF
    !
    ! Tobiska EUV91 Method:
    !
    IF (ISCALE .EQ. 2) THEN
       IF (HLYA .GT. 0.001) THEN
          HLYMOD = HLYA
       ELSE
          IF (HEIEW .GT. 0.001) THEN
             HLYMOD = HEIEW * 3.77847E9 + 8.40317E10
          ELSE
             HLYMOD = 8.70E8 * F107 + 1.90E11
          ENDIF
       ENDIF
       IF (HEIEW .GT. 0.001) THEN
          HEIMOD = HEIEW * 3.77847E9 + 8.40317E10
       ELSE
          HEIMOD = HLYMOD
       ENDIF
       DO  L=16,55
          SFLUX(L) = TCHR0(L) + TCHR1(L)*HLYMOD + TCHR2(L)*HEIMOD &
               + TCOR0(L) + TCOR1(L)*F107 + TCOR2(L)*F107A
       enddo
    ENDIF
    !
    ! Woods and Rottman (10 Nov. 1988) spectrum:
    !
    IF (ISCALE .EQ. 3) THEN
       DO  L=15,55
          SFLUX(L) = WAR1(L)
       end do
    ENDIF
    !
    ! Woods and Rottman (20 June 1989) spectrum:
    !
    IF (ISCALE .EQ. 4) THEN
       DO  L=15,55
          SFLUX(L) = WAR2(L)
       end do
    ENDIF
    !
    
    ! ! FISM spectrum for 46th day of 1997
    !
    IF (ISCALE .EQ. 5) THEN
       DO L=15,55
          SFLUX(L) = FISM(L)
       enddo
    ENDIF
    !
    
    ! Substitute in H Lyman-alpha and XUVFAC if provided:
    !
    IF (HLYA .GT. 0.001) SFLUX(12) = HLYA / 1.E9
    IF (XUVFAC .GT. 0.001) THEN
       XUVF = XUVFAC
    ELSE
       XUVF = 1.0
    ENDIF
    !
    ! Convert from gigaphotons to photons, etc.:
    !
    DO L=1,LMAX
       WAVE1(L) = WAVEL(L)
       WAVE2(L) = WAVES(L)
       IF (SFLUX(L) .LT. 0.0) SFLUX(L) = 0.0
       IF (WAVEL(L).LT.251.0 .AND. WAVES(L).GT.15.0) &
            SFLUX(L)=SFLUX(L)*XUVF
       SFLUX(L) = SFLUX(L) * 1.E9
    end do
    !


    RETURN
    
  END SUBROUTINE SSFLUX

  !============================================================================
  subroutine init_production
    use ModPlanetConst, ONLY: Planet_, rPlanet_I

    !\
    ! Earth
    !/
    integer, parameter :: nStatesMaxEarth=6 !only for setting data arrays
    ! branching ratios (O,O2,N2)
    real :: PROBO(nStatesMaxEarth,LMAX), PROBO2(nStatesMaxEarth,LMAX) 
    real :: PROBN2(nStatesMaxEarth,LMAX)
    ! absorption crossections (O,O2,N2)
    real :: SIGAO(LMAX), SIGAO2(LMAX), SIGAN2(LMAX)
    ! ionization crossections (O,O2,N2)
    real :: SIGIO(LMAX), SIGIO2(LMAX), SIGIN2(LMAX)
    ! named parameters for neutral species
    integer, parameter :: O_=1, O2_=2, N2_=3
    
    !\
    ! Jupiter/Saturn
    !/
    ! branching ratios (H2,H,CH4,He)
    real,allocatable :: PROBH2(:,:), PROBH(:,:),PROBCH4(:,:),PROBHe(:,:)
    ! absorption crossections (O,O2,N2)
    real :: SIGAH2(LMAX), SIGAH(LMAX), SIGACH4(LMAX),SIGAHe(LMAX)
    ! ionization crossections (O,O2,N2)
    real :: SIGIH2(LMAX), SIGIH(LMAX), SIGICH4(LMAX), SIGIHe(LMAX)
    ! named parameters for neutral species
    integer, parameter :: H2_=1, H_=2, CH4_=3, He_=4

  DATA ((PROBO(K,L),K=1,6),L=1,38) &
         / 120 * 0.00, &
         1.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
         1.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
         1.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
         1.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
         1.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
         1.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
         1.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
         0.52, 0.48, 0.00, 0.00, 0.00, 0.00, &
         0.43, 0.57, 0.00, 0.00, 0.00, 0.00, &
         0.40, 0.55, 0.05, 0.00, 0.00, 0.00, &
         0.30, 0.45, 0.25, 0.00, 0.00, 0.00, &
         0.41, 0.45, 0.24, 0.00, 0.00, 0.00, &
         0.30, 0.45, 0.25, 0.00, 0.00, 0.00, &
         0.29, 0.45, 0.26, 0.00, 0.00, 0.00, &
         0.29, 0.45, 0.25, 0.00, 0.00, 0.00, &
         0.29, 0.45, 0.26, 0.00, 0.00, 0.00, &
         0.28, 0.45, 0.26, 0.00, 0.00, 0.00, &
         0.28, 0.45, 0.27, 0.00, 0.00, 0.00/
    DATA ((PROBO(K,L),K=1,6),L=39,52) &
         / 0.28, 0.45, 0.27, 0.00, 0.00, 0.00, &
         0.27, 0.42, 0.26, 0.05, 0.00, 0.00, &
         0.26, 0.40, 0.25, 0.08, 0.00, 0.00, &
         0.26, 0.40, 0.25, 0.08, 0.00, 0.00, &
         0.26, 0.40, 0.25, 0.09, 0.00, 0.00, &
         0.25, 0.37, 0.24, 0.09, 0.04, 0.00, &
         0.25, 0.37, 0.24, 0.09, 0.04, 0.00, &
         0.25, 0.36, 0.23, 0.10, 0.06, 0.00, &
         0.25, 0.37, 0.23, 0.10, 0.05, 0.00, &
         0.25, 0.36, 0.23, 0.10, 0.06, 0.00, &
         0.30, 0.31, 0.20, 0.11, 0.07, 0.00, &
         0.37, 0.26, 0.17, 0.13, 0.06, 0.00, &
         0.29, 0.32, 0.21, 0.10, 0.08, 0.00, &
         0.30, 0.32, 0.21, 0.09, 0.08, 0.00/
    DATA ((PROBO(K,L),K=1,6),L=53,59) &
         / 0.30, 0.32, 0.21, 0.09, 0.08, 0.00, &
         0.30, 0.32, 0.21, 0.09, 0.08, 0.00, &
         0.30, 0.32, 0.21, 0.09, 0.08, 0.00, &
         0.30, 0.32, 0.21, 0.09, 0.08, 0.00, &
         0.30, 0.32, 0.21, 0.09, 0.08, 0.00, &
         0.30, 0.32, 0.21, 0.09, 0.08, 0.00, &
         0.30, 0.32, 0.21, 0.09, 0.08, 0.00/
    !
    DATA ((PROBO2(K,L),K=1,6),L=1,33) &
         /  90 * 0.00, &
         1.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
         1.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
         1.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
         1.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
         1.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
         1.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
         1.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
         1.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
         1.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
         1.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
         0.96, 0.04, 0.00, 0.00, 0.00, 0.00, &
         0.90, 0.10, 0.00, 0.00, 0.00, 0.00, &
         0.56, 0.44, 0.00, 0.00, 0.00, 0.00, &
         0.56, 0.40, 0.00, 0.00, 0.00, 0.00, &
         0.40, 0.46, 0.14, 0.00, 0.00, 0.00, &
         0.24, 0.36, 0.35, 0.05, 0.00, 0.00, &
         0.24, 0.36, 0.35, 0.05, 0.00, 0.00, &
         0.23, 0.38, 0.31, 0.07, 0.00, 0.00/
    DATA ((PROBO2(K,L),K=1,6),L=34,52) &
         / 0.35, 0.28, 0.21, 0.16, 0.00, 0.00, &
         0.30, 0.33, 0.21, 0.16, 0.00, 0.00, &
         0.36, 0.24, 0.22, 0.18, 0.00, 0.00, &
         0.36, 0.28, 0.13, 0.23, 0.00, 0.00, &
         0.42, 0.25, 0.12, 0.21, 0.00, 0.00, &
         0.42, 0.25, 0.12, 0.21, 0.00, 0.00, &
         0.42, 0.24, 0.12, 0.22, 0.00, 0.00, &
         0.40, 0.22, 0.12, 0.26, 0.00, 0.00, &
         0.37, 0.21, 0.12, 0.30, 0.00, 0.00, &
         0.36, 0.20, 0.12, 0.32, 0.00, 0.00, &
         0.35, 0.19, 0.11, 0.35, 0.00, 0.00, &
         0.35, 0.19, 0.11, 0.35, 0.00, 0.00, &
         0.34, 0.18, 0.11, 0.37, 0.00, 0.00, &
         0.34, 0.18, 0.11, 0.37, 0.00, 0.00, &
         0.34, 0.18, 0.11, 0.37, 0.00, 0.00, &
         0.33, 0.18, 0.10, 0.39, 0.00, 0.00, &
         0.30, 0.16, 0.09, 0.45, 0.00, 0.00, &
         0.20, 0.11, 0.07, 0.62, 0.00, 0.00, &
         0.10, 0.06, 0.04, 0.80, 0.00, 0.00/
    DATA ((PROBO2(K,L),K=1,6),L=53,59) &
         / 0.10, 0.06, 0.04, 0.80, 0.00, 0.00, &
         0.10, 0.06, 0.04, 0.80, 0.00, 0.00, &
         0.00, 0.00, 0.00, 1.00, 0.00, 0.00, &
         0.00, 0.00, 0.00, 1.00, 0.00, 0.00, &
         0.00, 0.00, 0.00, 1.00, 0.00, 0.00, &
         0.00, 0.00, 0.00, 1.00, 0.00, 0.00, &
         0.00, 0.00, 0.00, 1.00, 0.00, 0.00/
    !
    DATA ((PROBN2(K,L),K=1,6),L=1,38) &
         / 120 * 0.00, &
         0.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
         0.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
         0.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
         1.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
         1.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
         1.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
         1.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
         0.53, 0.47, 0.00, 0.00, 0.00, 0.00, &
         0.64, 0.36, 0.00, 0.00, 0.00, 0.00, &
         0.34, 0.66, 0.00, 0.00, 0.00, 0.00, &
         0.31, 0.59, 0.10, 0.00, 0.00, 0.00, &
         0.31, 0.59, 0.10, 0.00, 0.00, 0.00, &
         0.31, 0.59, 0.10, 0.00, 0.00, 0.00, &
         0.33, 0.57, 0.10, 0.00, 0.00, 0.00, &
         0.32, 0.58, 0.10, 0.00, 0.00, 0.00, &
         0.35, 0.55, 0.10, 0.00, 0.00, 0.00, &
         0.38, 0.53, 0.09, 0.00, 0.00, 0.00, &
         0.40, 0.47, 0.09, 0.00, 0.00, 0.04/
    DATA ((PROBN2(K,L),K=1,6),L=39,52) &
         / 0.42, 0.46, 0.08, 0.00, 0.00, 0.04, &
         0.44, 0.44, 0.08, 0.00, 0.00, 0.04, &
         0.35, 0.46, 0.10, 0.03, 0.00, 0.06, &
         0.33, 0.46, 0.10, 0.03, 0.00, 0.07, &
         0.25, 0.43, 0.10, 0.05, 0.01, 0.16, &
         0.22, 0.38, 0.09, 0.06, 0.05, 0.20, &
         0.22, 0.38, 0.09, 0.06, 0.05, 0.20, &
         0.20, 0.33, 0.07, 0.03, 0.10, 0.26, &
         0.20, 0.36, 0.07, 0.04, 0.09, 0.24, &
         0.19, 0.27, 0.07, 0.04, 0.12, 0.31, &
         0.18, 0.21, 0.07, 0.04, 0.15, 0.35, &
         0.17, 0.18, 0.07, 0.04, 0.17, 0.36, &
         0.17, 0.18, 0.07, 0.04, 0.17, 0.36, &
         0.17, 0.18, 0.07, 0.04, 0.17, 0.36/
    DATA ((PROBN2(K,L),K=1,6),L=53,59) &
         / 0.17, 0.18, 0.07, 0.04, 0.17, 0.36, &
         0.02, 0.02, 0.00, 0.00, 0.00, 0.96, &
         0.02, 0.02, 0.00, 0.00, 0.00, 0.96, &
         0.02, 0.02, 0.00, 0.00, 0.00, 0.96, &
         0.02, 0.02, 0.00, 0.00, 0.00, 0.96, &
         0.02, 0.02, 0.00, 0.00, 0.00, 0.96, &
         0.02, 0.02, 0.00, 0.00, 0.00, 0.96/
    !
    ! NB - absorption and ionization cross sections are multiplied by 1.E-18
    ! on first call.
    !
    DATA SIGAO /  18 * 0.00, &
         0.00, 0.00, 1.66, 3.85, 4.06, 4.08, &
         4.08, 4.08, 4.08, 7.06, 8.52, 8.98, &
         13.10,13.19,13.30,12.88,13.20,12.44, &
         12.23,12.00,11.18,11.04, 9.64, 9.79, &
         8.68, 7.69, 7.68, 6.63, 7.13, 6.04, &
         5.22, 2.95, 1.73, 0.61, 0.16, 0.05, &
         0.51, 0.07, .012, .002, .0002/
    !
    DATA SIGAO2/ 0.50, 1.50, 3.40, 6.00,10.00,13.00, &
         15.00,12.00, 2.20, 0.40,13.00, 0.01, &
         1.40, 0.40, 1.00, 1.23, 1.15, 1.63, &
         22.15, 4.00,12.12, 8.54,16.63,24.32, &
         26.66,18.91,20.82,28.55,27.48,21.49, &
         25.97,27.33,25.19,26.64,22.81,25.95, &
         24.56,22.84,21.79,20.13,18.19,18.40, &
         17.35,16.64,16.61,14.74,15.69,13.54, &
         11.04, 7.11, 3.76, 1.21, 0.32, 0.10, &
         1.02, 0.14, .024, .004, .0004/
    !
    DATA SIGAN2/  18 * 0.00, &
         38.40, 0.70,19.43,34.88,15.06,16.91, &
         16.50,21.19,35.46,24.26,21.82,26.42, &
         23.36,23.37,22.80,22.78,22.40,24.13, &
         24.63,23.47,23.17,21.64,16.44,16.91, &
         13.79,11.70,11.67,10.57,10.90,10.21, &
         8.52, 4.80, 2.29, 0.72, 0.24, 1.16, &
         0.48, 0.09, .015, .003, .0003/
    !
    DATA SIGIO /  18 * 0.00, &
         0.00, 0.00, 1.66, 3.85, 4.06, 4.08, &
         4.08, 4.08, 4.08, 7.06, 8.52, 8.98, &
         13.10,13.19,13.30,12.88,13.20,12.44, &
         12.23,12.00,11.18,11.04, 9.64, 9.79, &
         8.68, 7.69, 7.68, 6.63, 7.13, 6.04, &
         5.22, 2.95, 1.73, 0.61, 0.16, 0.05, &
         0.51, 0.07, .012, .002, .0002/
    !
    DATA SIGIO2/ 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
         0.00, 0.00, 0.00, 0.00, 0.00, 0.00, &
         0.00, 0.00, 0.00, 0.19, 0.00, 1.00, &
         16.36, 2.50, 7.94, 5.17, 6.24,12.02, &
         12.86,10.87,11.21,23.58,23.95,20.80, &
         25.95,27.33,25.19,26.64,22.81,25.95, &
         24.56,22.84,21.79,20.13,18.19,18.40, &
         17.35,16.64,16.61,14.74,15.69,13.54, &
         11.04, 7.11, 3.76, 1.21, 0.32, 0.10, &
         1.02, 0.14, .024, .004, .0004/
    !
    DATA SIGIN2/  18 * 0.00, &
         0.00, 0.00, 0.00, 0.00, 0.00, 9.90, &
         8.67,11.89,23.77,21.03,21.02,25.06, &
         23.36,23.37,22.80,22.78,22.40,24.13, &
         24.63,23.47,23.17,21.64,16.44,16.91, &
         13.79,11.70,11.67,10.57,10.90,10.21, &
         8.52, 4.80, 2.29, 0.72, 0.24, 1.16, &
         0.48, 0.09, .015, .003, .0003/
   


    select case(NamePlanet_I(Planet_))
    case('EARTH')
       ! Earth
       ! set the surface gravity
       gSurface = 978.1

       ! set the neutral parameters
       nNeutral=3
       
       call allocate_arrays

       ! set the neutral mass per species
       NeutralMassAMU_I = (/16.0,32.0,28.0/)

       ! set auger e- mean energies and wavelength thresholds
       AugEnergy_I = (/533., 533., 402./)
       AugThreshold_I = (/23., 23., 32./)
! ***something else Auger e- related that seems to belong here
!    index of wavelength bin?
       LAUG = (/54,54,53/)

       ! set number of states per species
       nStatesPerSpecies_I=(/5,4,6/)

       nStatesMax = max(nStatesPerSpecies)

       ! allocate & set ionization potentials per state, per species
       if (.not.allocated(TPOT)) &
            allocate(TPOT(nStatesMax,nNeutral))
       TPOT(:,O_) = (/13.61, 16.93, 18.63, 28.50, 40.00,  0.00/)
       TPOT(:,O2_)= (/12.07, 16.10, 18.20, 20.00,  0.00,  0.00/) 
       TPOT(:,N2_)= (/15.60, 16.70, 18.80, 30.00, 34.80, 25.00/)

       !set branching ratios and crossections
       DO  L=1,LMAX
          DO  K=1,nStatesMax
             PROB(K,1,L) = PROBO(K,L)
             PROB(K,2,L) = PROBO2(K,L)
             PROB(K,3,L) = PROBN2(K,L)
          end do
       end DO
       
       DO  L=1,LMAX
          SIGABS(1,L) = SIGAO(L)  * 1.E-18
          SIGABS(2,L) = SIGAO2(L) * 1.E-18
          SIGABS(3,L) = SIGAN2(L) * 1.E-18
          SIGION(1,L) = SIGIO(L)  * 1.E-18
          SIGION(2,L) = SIGIO2(L) * 1.E-18
          SIGION(3,L) = SIGIN2(L) * 1.E-18
       end DO
       
       
    case('JUPITER')
       ! Jupiter
       ! set the surface gravity
       gSurface = 2479.0

       ! set the neutral parameters
       nNeutral=3

       call allocate_arrays

       ! set the neutral mass per species
       NeutralMassAMU_I = (/2.0,1.0,16.0,4.0/)

       ! set auger e- mean energies and wavelength thresholds
       ! both set to zero for no Auger production
       AugEnergy_I = (/0.0, 0.0, 0.0, 0.0/)
       AugThreshold_I = (/0.0, 0.0, 0.0, 0.0/)
       ! *** set LAUG here?

       ! set number of states per species
       nStatesPerSpecies_I=(/2,1,5,1/)

       nStatesMax = max(nStatesPerSpecies)

       ! allocate & set ionization potentials per state, per species
       if (.not.allocated(TPOT)) &
            allocate(TPOT(nStatesMax,nNeutral))
       TPOT(:,H2_) = (/0.00, 0.00, 0.00, 0.00, 0.00/)
       TPOT(:,H_)  = (/13.5, 0.00, 0.00, 0.00, 0.00/) 
       TPOT(:,CH4_)= (/18.0, 20.0, 15.0, 19.0, 12.1/)
       TPOT(:,He_) = (/0.00, 0.00, 0.00, 0.00, 0.00/) 

       !set branching ratios and crossections
       DO  L=1,LMAX
          DO  K=1,nStatesMax
             PROB(K,1,L) = PROBH2(K,L)
             PROB(K,2,L) = PROBH(K,L)
             PROB(K,3,L) = PROBCH4(K,L)
             PROB(K,4,L) = PROBHe(K,L)
          end do
       end DO
       
       DO  L=1,LMAX
          SIGABS(1,L) = SIGAH2(L)  * 1.E-18
          SIGABS(2,L) = SIGAH(L) * 1.E-18
          SIGABS(3,L) = SIGACH4(L) * 1.E-18
          SIGABS(4,L) = SIGAHe(L) * 1.E-18
          SIGION(1,L) = SIGIH2(L)  * 1.E-18
          SIGION(2,L) = SIGIH(L) * 1.E-18
          SIGION(3,L) = SIGICH4(L) * 1.E-18
          SIGION(4,L) = SIGIHe(L) * 1.E-18
       end DO

    end select


  end subroutine init_production


  subroutine allocate_arrays
    ! array to hold neutral mass
    if (.not.allocated(NeutralMassAMU_I)) &
         allocate(NeutralMassAMU_I(nNeutral))
    ! states per species
    if (.not.allocated(nStatesPerSpecies_I)) &
         allocate(nStatesPerSpecies_I(nNeutral))
    ! auger e- mean energy per species
    if (.not.allocated(AugEnergy_I)) &
         allocate(AugEnergy_I(nNeutral))
    ! auger e- wavelength threshold per species
    if (.not.allocated(AugThreshold_I)) &
         allocate(AugThreshold_I(nNeutral))
    ! auger something
    if (.not.allocated(LAUG)) &
         allocate(LAUG(nNeutral))

    
  end subroutine allocate_arrays

end Module ModSeProduction
