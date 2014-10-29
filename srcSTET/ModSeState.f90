Module ModSeState
  implicit none
  
  private !except
  
  real, public :: Time
  !\
  ! FLUX VARIABLES; cm-2 s-1 eV-1 sr-1
  !/
  
  ! Flux in the plasmasphere from first to second (up) ionosphere 
  real, public,allocatable :: phiup(:,:,:,:),phidn(:,:,:,:)
  
  ! Flux in the ionosphere from first to second (up) ionosphere 
  real, public,allocatable :: iphiup(:,:,:,:),iphidn(:,:,:,:)
  
  ! Flux from previous time step
  real, public,allocatable :: lphiup(:,:,:,:),lphidn(:,:,:,:)
  real, public,allocatable :: liphiup(:,:,:,:),liphidn(:,:,:,:)
  
  !Omnidirectional flux for upward/downward-directed hemisphere
  real, public,allocatable :: specup(:,:,:),specdn(:,:,:)
  
  real, public :: delt,epsilon
  
  ! Ring current parameters for coulomb collisions (note that ring.dat would 
  ! need to be read in the future to use this feature.
  integer,parameter :: IRC=0 
  real :: NRC(IRC),TRC(IRC),MRC(IRC)          !,Iterm1,Iterm2

  real,allocatable :: SRC(:)
  integer,parameter :: NRing=6 ! ring current Maxwellian fits

  ! ilocal defines the number of points at the bottom of each ionosphere where 
  ! transport is explicitly turned off. default is shown
  integer :: ilocal = 12

  ! Arrays for electron production in ionosphere. These are only used to store 
  ! these values for output. They are already added into Qstar for the purposes 
  ! of calculation.
  real,allocatable :: Qestar_ICI(:,:,:) ! secondary production in ionosphere 
  real,allocatable :: Qpstar_ICI(:,:,:) ! Beam-induced production 
                                        ! (non-degredating beam)
  
  
  !public methods
  public :: update_se_state
  public :: se_update_state_test
contains
  SUBROUTINE update_se_state(iLine,IsOpen)
    !*  This subroutine used to be the main program, until other operators
    !*  added (3/23/95). Now it is one of several operators called by the
    !*  driver program, 'stet1.f' (Super_Thermal_electron_Transport_Model).
    
    USE ModSeGrid, only: nEnergy, EnergyGrid_I, DeltaE_I, EnergyMin, EnergyMax,&
         EqAngleGrid_IG,nThetaAlt_II,dThetaEnd_II,mu_III,&
         FieldLineGrid_IC,nTop,Bfield_IC, BFieldEq_I, &
         nIono, nPlas, nPoint,nAngle
    
    USE ModSeBackground, only: eThermalDensity_IC,eThermalTemp_IC

    use ModMath, only: midpnt_int
    
    use ModNumConst,    ONLY: cPi

    IMPLICIT NONE
!    INCLUDE 'numbers.h'
    integer, intent(in) :: iLine
    logical, intent(in) :: IsOpen
    REAL h,p,Fcheck, sigmaO,muO,velt,coef, lbeta(nPoint),alpha(nAngle), &
         sigma(nAngle),Flastj,kk, del1,del2, newphi,beta(nPoint),sigO1, &
         Fsum, delE,Flasti, Flastt,Qstar, Theta,s1,s2
    INTEGER i,j,k,iPlas,jj,Ist, &
         flag,count,warning,SPick,space,ii

    !set maximum iterations
    integer, parameter :: countmax=600

    real :: cascade, lossum
    
    !Altitude range variables
    integer :: nAltMin, nAltMax
    
    real,parameter :: cEVtoCMperS = 5.88e7 ! convert energy to velocity

    !---------------------------------------------------------------------------
    
    !initialize lbeta to 0
    lbeta(:)=0.0
    
    !      initialize warning to 0
    warning =0
    alpha(1)=1.
    sigma(1)=0.
    kk=1.
    !  Start the energy loop
    !      PRINT *, 'MainPlas, t=',t
    ENERGY: DO j=nEnergy,1,-1
       delE=.5*(DeltaE_I(j+1)+DeltaE_I(j))
       IF (j.EQ.nEnergy) delE=DeltaE_I(j)
       velt=cEVtoCMperS*SQRT(EnergyGrid_I(j))*delt
       flag=1
       count=0
       !  Start the iteration loop
       ITERATION: DO WHILE ((flag.EQ.1).AND.(count.LT.countmax))
          flag=0
          count=count+1
          IF (count.GT.countmax-2) WRITE (10,*) 'Count: ',count
          !*  Upward Region #2: the plasmasphere
          do k=0,nThetaAlt_II(iLine,nIono)
             phiup(iLine,k,0,j)=iphiup(iLine,k,nIono,j)
          end do
          h=FieldLineGrid_IC(iLine,nIono+1)-FieldLineGrid_IC(iLine,nIono)
          
          !check if field line is open and set bounds for up and down loop
          if (IsOpen) then
             ! open line so detach hemispheres
             nAltMin = nIono+1
             nAltMax = nTop
          else
             ! closed line so keep hemispheres attached
             nAltMin = nIono+1
             nAltMax = nIono+nPlas
          endif
          
          FIELDLINE_UPWARD: DO i=nAltMin,nAltMax
!             write(*,*) 'Start FieldLine up'
             CALL CoulVar(i,beta(i),sigO1,eThermalTemp_IC(iLine,i),eThermalDensity_IC(iLine,i),EnergyGrid_I(j))

             !KLUDGE kill pitchangle scattering
             !sigO1=0.0
             !END KLUDGE
             iPlas=i-nIono
             coef=BFieldEq_I(iLine)/Bfield_IC(iLine,i)
             CALL IonoVar3(Qstar, cascade, lossum)
             do k=1,nThetaAlt_II(iLine,i)-1
                CALL ThetaVar1(nAngle,Theta,del1,del2,k,EqAngleGrid_IG(iLine,0:nAngle))
                IF (k.EQ.nThetaAlt_II(iLine,i)-1) del2=dThetaEnd_II(iLine,i)
                muO=cos(Theta)
                sigmaO=sigO1
                s1=sigmaO
                s2=1.
                Flasti=phiup(iLine,k,iPlas-1,j)
                Flastj=phiup(iLine,k,iPlas,j+1)
                Flastt=lphiup(iLine,k,iPlas,j)
                CALL NumCalcVar(alpha(k+1),alpha(k),sigma(k+1),sigma(k), &
                     mu_III(iLine,k,i),coef,muO,del1,del2,beta(i),lbeta(i), &
                     Flasti,Flastj,Flastt,velt,h,delE,Qstar,kk,s1,s2,&
                     cascade,lossum)
                !            CALL CheckWarn(sigma(k+1),warning,1,*9999)
                !            CALL CheckWarn(alpha(k+1),warning,1,*9999)
             end do
             IF (i.EQ.nTop) THEN
                p=-2.*s1/dThetaEnd_II(iLine,i)
             ELSE
                p=-s1/SQRT(ABS(1./coef-1.))
             END IF
             Flastt=lphiup(iLine,nThetaAlt_II(iLine,i),iPlas,j)
             Fsum=phiup(iLine,nThetaAlt_II(iLine,i)-1,iPlas,j)+phidn(iLine,nThetaAlt_II(iLine,i)-1,iPlas,j)
             Fcheck=phiup(iLine,nThetaAlt_II(iLine,i),iPlas,j)
             Flastj=phiup(iLine,nThetaAlt_II(iLine,i),iPlas,j+1)
             newphi=(Flastt/velt-p*Fsum/(2.*dThetaEnd_II(iLine,i)) &
                  +s2*(lbeta(i)* Flastj/delE+Qstar+cascade))/ &
                  (1/velt-p/dThetaEnd_II(iLine,i)+s2*(beta(i)/delE+lossum))
             !            CALL CheckWarn(newphi,warning,2,*9999)
             
             CALL CheckConv(newphi,Fcheck,epsilon,flag)
             phiup(iLine,nThetaAlt_II(iLine,i),iPlas,j)=newphi
             do k=nThetaAlt_II(iLine,i),1,-1
                newphi=alpha(k)*phiup(iLine,k,iPlas,j)+sigma(k)
                Fcheck=phiup(iLine,k-1,iPlas,j)
                !            CALL CheckWarn(newphi,warning,3,*9999)
                CALL CheckConv(newphi,Fcheck,epsilon,flag)
                phiup(iLine,k-1,iPlas,j)=newphi
             end do
             CALL midpnt_int(specup(iLine,j,i),phiup(iLine,0,iPlas,j),&
                  mu_III(iLine,0,i),1, nThetaAlt_II(iLine,i)+1,nAngle+1,1)
             specup(iLine,j,i)=-.5*specup(iLine,j,i)
             !  Write to a file (if the solution isn't converging)
             IF (count.GT.countmax-2) THEN
                WRITE (10,*) time,j,i,iPlas,nThetaAlt_II(iLine,i)/4,&
                     flag,(phiup(iLine,k,iPlas,j),k=0,nThetaAlt_II(iLine,i))
                WRITE (10,*) time,j,i,iPlas,nThetaAlt_II(iLine,i),flag,&
                     (phiup(iLine,k,iPlas,j),k=nThetaAlt_II(iLine,i)/2+1, &
                     nThetaAlt_II(iLine,i),2)
             END IF
             h=FieldLineGrid_IC(iLine,i+1)-FieldLineGrid_IC(iLine,i)            
             
! End Upward Region #2 loop
          end DO FIELDLINE_UPWARD
          !*  Downward Region #2: the plasmasphere

          !check if field line is open and set bounds for up and down loop
          if (IsOpen) then
             ! open line so detach hemispheres
             nAltMin = nIono+1
             nAltMax = nTop
             !set step
             h = &
                  FieldLineGrid_IC(iLine,nTop) &
                  -FieldLineGrid_IC(iLine,nTop+1)
             !fill bc from iono
             do k=0,nThetaAlt_II(iLine,nTop)
                phidn(iLine,k,nTop+1,j) = 0.0
             end do
          else
             ! closed line so keep hemispheres attached
             nAltMin = nIono+1
             nAltMax = nIono+nPlas
             !set step
             h = &
                  FieldLineGrid_IC(iLine,nIono+nPlas) &
                  -FieldLineGrid_IC(iLine,nIono+nPlas+1)
             !fill bc from iono
             do k=0,nThetaAlt_II(iLine,nIono)
                phidn(iLine,k,nPlas+1,j)=iphidn(iLine,k,nIono+1,j)
             end do
          endif
          
          FIELDLINE_DOWN: DO i=nAltMax,nAltMin,-1
!             write(*,*) 'Start FieldLine Down'
             CALL CoulVar(i,beta(i),sigO1,eThermalTemp_IC(iLine,i),eThermalDensity_IC(iLine,i),EnergyGrid_I(j))
             !KLUDGE kill pitchangle scattering
             !sigO1=0.0
             !END KLUDGE
             iPlas=i-nIono
             coef=BFieldEq_I(iLine)/Bfield_IC(iLine,i)
             CALL IonoVar3(Qstar, cascade, lossum)
             DO k=1,nThetaAlt_II(iLine,i)-1
                CALL ThetaVar1(nAngle,Theta,del1,del2,k,EqAngleGrid_IG(iLine,0:nAngle))
                Theta=cPi-Theta
                del1=-del1
                del2=-del2
                IF (k.EQ.nThetaAlt_II(iLine,i)-1) del2=-dThetaEnd_II(iLine,i)
                muO=cos(Theta)
                sigmaO=sigO1
                s1=sigmaO
                s2=1.
                Flasti=phidn(iLine,k,iPlas+1,j)
                Flastj=phidn(iLine,k,iPlas,j+1)
                Flastt=lphidn(iLine,k,iPlas,j)
                CALL NumCalcVar(alpha(k+1),alpha(k),sigma(k+1),sigma(k), &
                     -mu_III(iLine,k,i),coef,muO,del1,del2,beta(i),lbeta(i), &
                     Flasti,Flastj,Flastt,velt,h,delE,Qstar,kk,s1,s2, &
                     cascade,lossum)
                !            CALL CheckWarn(sigma(k+1),warning,4,*9999)
                !            CALL CheckWarn(alpha(k+1),warning,4,*9999)
             enddo
             IF (i.EQ.nTop) THEN
                p=2.*s1/dThetaEnd_II(iLine,i)
             ELSE
                p=s1/SQRT(ABS(1./coef-1.))
             END IF
             Flastt=phidn(iLine,nThetaAlt_II(iLine,i),iPlas,j)
             Fsum=phidn(iLine,nThetaAlt_II(iLine,i)-1,iPlas,j)+phiup(iLine,nThetaAlt_II(iLine,i)-1,iPlas,j)
             Fcheck=phidn(iLine,nThetaAlt_II(iLine,i),iPlas,j)
             Flastj=phidn(iLine,nThetaAlt_II(iLine,i),iPlas,j+1)
             newphi=(Flastt/velt+p*Fsum/(2.*dThetaEnd_II(iLine,i))+s2*(lbeta(i)* &
                  Flastj/delE+Qstar+cascade))&
                  /(1/velt+p/dThetaEnd_II(iLine,i)+s2*(beta(i)/delE+lossum))
             !            CALL CheckWarn(newphi,warning,5,*9999)
             CALL CheckConv(newphi,Fcheck,epsilon,flag)
             phidn(iLine,nThetaAlt_II(iLine,i),iPlas,j)=newphi
             DO k=nThetaAlt_II(iLine,i),1,-1
                newphi=alpha(k)*phidn(iLine,k,iPlas,j)+sigma(k)
                Fcheck=phidn(iLine,k-1,iPlas,j)
                !            CALL CheckWarn(newphi,warning,6,*9999)
                CALL CheckConv(newphi,Fcheck,epsilon,flag)
                phidn(iLine,k-1,iPlas,j)=newphi
             end DO
             CALL midpnt_int(specdn(iLine,j,i),phidn(iLine,0,iPlas,j),&
                  mu_III(iLine,0,i),1, nThetaAlt_II(iLine,i)+1,nAngle+1,1)
             specdn(iLine,j,i)=.5*specdn(iLine,j,i)
             !  Write to a file (if the solution isn't converging)
             IF (count.GT.countmax-2) THEN
                WRITE (10,*) time,j,i,iPlas,flag,&
                     (phidn(iLine,k,iPlas,j),k=0,nThetaAlt_II(iLine,i)/2,2)
                WRITE (10,*) time,j,i,iPlas,flag,&
                     (phidn(iLine,k,iPlas,j),k=nThetaAlt_II(iLine,i)/2+1, &
                     nThetaAlt_II(iLine,i),2)
             END IF
             h=FieldLineGrid_IC(iLine,i-1)-FieldLineGrid_IC(iLine,i)            ! End Downward Region #2 loop
          end DO FIELDLINE_DOWN
       END DO ITERATION                  ! End of iteration loop
       CALL CheckFlag(flag,1,warning,-1,*9999)
       !  Save the calculated values for the next time and energy steps
       do i=nIono+1,nPoint-nIono
          lbeta(i)=beta(i)
       end do
       !      PRINT *, 'Energy convergence:',j,EnergyGrid_I(j),count
    end DO ENERGY
    
9999 IF (warning.GT.0) THEN
       PRINT *, 'Negative Densities Occurred:',newphi,warning
       PRINT *, time,j,i,k,Bfield_IC(iLine,i),BFieldEq_I(iLine),&
            BFieldEq_I(iLine),FieldLineGrid_IC(iLine,i),h,Theta,muO,&
            mu_III(iLine,k,i),del1, dThetaEnd_II(iLine,i),beta(i),sigmaO,p,count,&
            Qstar
       PRINT *, phiup(iLine,k,iPlas,j),phidn(iLine,k,iPlas,j)
       STOP
    ELSE IF (warning.LT.0) THEN
       PRINT *, 'No convergence: ',count,Time,warning,j
       STOP
    END IF
    !  Format for the nonconvergent output
!!!27  FORMAT (100E10.3,4I4,1P,100E10.3)
    
    RETURN
  END SUBROUTINE update_se_state
  
  !============================================================================
  SUBROUTINE update_se_state_iono(iLine,IsIono1,nNeutral,NeutralDens_IC,&
       SIGS,SIGI,SIGA,ePhotoProdSpec_IC)
    !*  Same as update_se_state but for ionosphere, includes sources
    
    USE ModSeGrid, only: nEnergy, EnergyGrid_I, DeltaE_I, EnergyMin, EnergyMax,&
         EqAngleGrid_IG,nThetaAlt_II,dThetaEnd_II,mu_III,&
         FieldLineGrid_IC,nTop,Bfield_IC, BFieldEq_I, &
         nIono, nPlas, nPoint,nAngle
    
    USE ModSeBackground, only: eThermalDensity_IC,eThermalTemp_IC

    use ModMath, only: midpnt_int
    
    use ModNumConst,    ONLY: cPi

    IMPLICIT NONE
!    INCLUDE 'numbers.h'
    integer, intent(in) :: iLine
    logical, intent(in) :: IsIono1
    
    !neutral atmosphere inputs
    integer, intent(in) :: nNeutral
    real   , intent(in) :: NeutralDens_IC(nNeutral,nIono)

    !incomming crossections
    real   , intent(in) :: SIGS(nNeutral,nEnergy)
    real   , intent(in) :: SIGI(nNeutral,nEnergy,nEnergy)
    real   , intent(in) :: SIGA(nNeutral,nEnergy,nEnergy)
    
    ! incomming photo electron production spectrum
    real   , intent(in) :: ePhotoProdSpec_IC(nEnergy,nIono)
    
    REAL h,p,Fcheck, sigmaO,muO,velt,coef, lbeta(nPoint),alpha(nAngle), &
         sigma(nAngle),Flastj,kk, del1,del2, newphi,beta(nPoint),sigO1, &
         Fsum, delE,Flasti, Flastt,Qstar, Theta,s1,s2
    INTEGER i,j,k,jj,Ist, &
         flag,count,warning,SPick,space,ii
    !inidicies that store full ionosphere grid(2*nIono points) and half 
    ! ionosphere grid (nIono points). These convert position along the field 
    ! line into these sub grids.
    integer :: iIono, iIonoHalf

    !set maximum iterations
    integer, parameter :: countmax=600

    real :: cascade, lossum
    
    !Altitude range variables
    integer :: nAltMin, nAltMax
    
    ! array to hold temporarily the values being overwritten by BCs
    real :: fhold(0:nAngle)
    real,parameter :: cEVtoCMperS = 5.88e7 ! convert energy to velocity

    !---------------------------------------------------------------------------
    
    !initialize lbeta to 0
    lbeta(:)=0.0
    
    !      initialize warning to 0
    warning =0
    alpha(1)=1.
    sigma(1)=0.
    kk=1.
    !  Start the energy loop
    !      PRINT *, 'MainPlas, t=',t
    ENERGY: DO j=nEnergy,1,-1
       delE=.5*(DeltaE_I(j+1)+DeltaE_I(j))
       IF (j.EQ.nEnergy) delE=DeltaE_I(j)
       velt=cEVtoCMperS*SQRT(EnergyGrid_I(j))*delt
       flag=1
       count=0
       !  Start the iteration loop
       ITERATION: DO WHILE ((flag.EQ.1).AND.(count.LT.countmax))
          flag=0
          count=count+1
          IF (count.GT.countmax-2) WRITE (10,*) 'Count: ',count

          ! Upward Region

          ! set BC when in iono2. Put BC in iphiup at top of iono1 and hold 
          ! the current value to put it back after the calculation
          if (.not.IsIono1)then
             do k=0,nThetaAlt_II(iLine,nIono)
                Fhold(k)=iphiup(iLine,k,nIono,j)
             end do

             ! fill BC from the the plasmasphere
             do k=0,nThetaAlt_II(iLine,nIono)
                iphiup(iLine,k,nIono,j)=phiup(iLine,i,nPlas,j)
             end do

             ! Add precip info here
          endif
             
          !set bounds for ionosphere loop depending if you are in iono 1 or 2
          if (IsIono1) then
             ! set upward bounds for iono1
             nAltMin = 1
             nAltMax = nIono
             h=FieldLineGrid_IC(iLine,nAltMin+1)-FieldLineGrid_IC(iLine,nAltMin)
          else
             ! set upward bounds for iono2
             nAltMin = nPoint-nIono+1
             nAltMax = nPoint
             h=FieldLineGrid_IC(iLine,nAltMin)-FieldLineGrid_IC(iLine,nAltMin-1)
          endif


          
          FIELDLINE_UPWARD: DO i=nAltMin,nAltMax
!             write(*,*) 'Start FieldLine up'
             CALL CoulVar(i,beta(i),sigO1,eThermalTemp_IC(iLine,i),eThermalDensity_IC(iLine,i),EnergyGrid_I(j))

             ! set iIono and iIonoHalf indices based on 1st or 2nd ionosphere
             ! remember, iIono is index spanning both ionospheres, 
             ! iIonoHalf is index spanning just one ionosphere
             if(IsIono1) then
                iIono=i
                iIonoHalf=i
             else
                iIono=i-nPlas
                iIonoHalf = nPoint-i+1
             end if
             
             !set the b-field ratio coeficient
             coef=BFieldEq_I(iLine)/Bfield_IC(iLine,i)
             
             ! get sigma0 and electron production (primary+secondary)
             CALL get_sigma0_and_eprod(j,nNeutral,sigO1,&
                  ePhotoProdSpec_IC(j,iIonoHalf), &
                  NeutralDens_IC(:,iIonoHalf),SIGS,SIGI,SIGA,Qstar, &
                  specup(iLine,:,i),specdn(iLine,:,i),&
                  Qestar_ICI(iLine,iIono,j),Qpstar_ICI(iLine,iIono,j))

             do k=1,nThetaAlt_II(iLine,i)-1
                call get_cascade_and_lossum(iLine,nNeutral,SIGA,&
                     NeutralDens_IC(:,iIonoHalf),j,iIono,k,cascade,lossum,1)
                CALL ThetaVar1(nAngle,Theta,del1,del2,k,&
                     EqAngleGrid_IG(iLine,0:nAngle))
                IF (k.EQ.nThetaAlt_II(iLine,i)-1) del2=dThetaEnd_II(iLine,i)
                muO=cos(Theta)
                sigmaO=sigO1
                s1=sigmaO
                s2=1.
                Flasti=iphiup(iLine, k,iIono-1,j)
                Flastj=iphiup(iLine, k,iIono,j+1)
                Flastt=liphiup(iLine,k,iIono,j)
                
                !When we are in the ilocal region kill the transport terms in 
                ! coeficients by setting kk to 0
                kk=1
                if (iIonoHalf <= iLocal) kk=0
                CALL NumCalcVar(alpha(k+1),alpha(k),sigma(k+1),sigma(k), &
                     mu_III(iLine,k,i),coef,muO,del1,del2,beta(i),lbeta(i), &
                     Flasti,Flastj,Flastt,velt,h,delE,Qstar,kk,s1,s2,&
                     cascade,lossum)
                !            CALL CheckWarn(sigma(k+1),warning,1,*9999)
                !            CALL CheckWarn(alpha(k+1),warning,1,*9999)
             end do
             
             p=-s1/SQRT(ABS(1./coef-1.))

             Flastt=liphiup(iLine,nThetaAlt_II(iLine,i),iIono,j)
             Fsum=iphiup(iLine,nThetaAlt_II(iLine,i)-1,iIono,j)&
                  +iphidn(iLine,nThetaAlt_II(iLine,i)-1,iIono,j)
             Fcheck=iphiup(iLine,nThetaAlt_II(iLine,i),iIono,j)
             Flastj=iphiup(iLine,nThetaAlt_II(iLine,i),iIono,j+1)
             call get_cascade_and_lossum(iLine,nNeutral,SIGA,&
                  NeutralDens_IC(:,iIonoHalf),j,iIono,nThetaAlt_II(iLine,i),&
                  cascade,lossum,1)
             newphi=(Flastt/velt-p*Fsum/(2.*dThetaEnd_II(iLine,i)) &
                  +s2*(lbeta(i)* Flastj/delE+Qstar+cascade))/ &
                  (1/velt-p/dThetaEnd_II(iLine,i)+s2*(beta(i)/delE+lossum))
             !            CALL CheckWarn(newphi,warning,2,*9999)
             
             CALL CheckConv(newphi,Fcheck,epsilon,flag)
             iphiup(iLine,nThetaAlt_II(iLine,i),iIono,j)=newphi
             do k=nThetaAlt_II(iLine,i),1,-1
                newphi=alpha(k)*iphiup(iLine,k,iIono,j)+sigma(k)
                Fcheck=iphiup(iLine,k-1,iIono,j)
                !            CALL CheckWarn(newphi,warning,3,*9999)
                CALL CheckConv(newphi,Fcheck,epsilon,flag)
                iphiup(iLine,k-1,iIono,j)=newphi
             end do
             CALL midpnt_int(specup(iLine,j,i),iphiup(iLine,:,iIono,j),&
                  mu_III(iLine,:,i),1, nThetaAlt_II(iLine,i)+1,nAngle+1,1)
             specup(iLine,j,i)=-.5*specup(iLine,j,i)
             !  Write to a file (if the solution isn't converging)
             IF (count.GT.countmax-2) THEN
                WRITE (10,*) time,j,i,iIono,nThetaAlt_II(iLine,i),&
                     flag,(iphiup(iLine,k,iIono,j),k=0,nThetaAlt_II(iLine,i))
             END IF
             h=FieldLineGrid_IC(iLine,i+1)-FieldLineGrid_IC(iLine,i)            
             
             ! End Upward Region of ionosphere
          end DO FIELDLINE_UPWARD
          !put back the iphiup at top of iono1 (which was storing plasmasph bc)
          do k=0,nThetaAlt_II(iLine,nIono)
             iphiup(iLine,k,nIono,j)=Fhold(k)
          end do

          
          !*  Downward Region of the ionosphere

          ! set BC when in iono1. Put BC in iphidn at top of iono2 and hold 
          ! the current value to put it back after the calculation
          if (IsIono1)then
             do k=0,nThetaAlt_II(iLine,nIono)
                Fhold(k)=iphidn(iLine,k,nIono+1,j)
             end do

             ! fill BC from the the plasmasphere
             do k=0,nThetaAlt_II(iLine,nIono)
                iphidn(iLine,k,nIono+1,j)=phidn(iLine,i,1,j)
             end do

             ! Add precip info here
          endif

          !set bounds for ionosphere loop depending if you are in iono 1 or 2
          if (IsIono1) then
             ! set upward bounds for iono1
             nAltMin = 1
             nAltMax = nIono
             h=FieldLineGrid_IC(iLine,nAltMax)-FieldLineGrid_IC(iLine,nAltMax+1)
          else
             ! set upward bounds for iono2
             nAltMin = nPoint-nIono+1
             nAltMax = nPoint
             h=FieldLineGrid_IC(iLine,nAltMax-1)-FieldLineGrid_IC(iLine,nAltMax)
          endif
          

          FIELDLINE_DOWN: DO i=nAltMax,nAltMin,-1
!             write(*,*) 'Start FieldLine Down'
             CALL CoulVar(i,beta(i),sigO1,eThermalTemp_IC(iLine,i),eThermalDensity_IC(iLine,i),EnergyGrid_I(j))

             ! set iIono and iIonoHalf indices based on 1st or 2nd ionosphere
             ! remember, iIono is index spanning both ionospheres, 
             ! iIonoHalf is index spanning just one ionosphere
             if(IsIono1) then
                iIono=i
                iIonoHalf=i
             else
                iIono=i-nPlas
                iIonoHalf = nPoint-i+1
             end if

             !set the b-field ratio coeficient
             coef=BFieldEq_I(iLine)/Bfield_IC(iLine,i)

             ! get sigma0 and electron production (primary+secondary)
             CALL get_sigma0_and_eprod(j,nNeutral,sigO1,&
                  ePhotoProdSpec_IC(j,iIonoHalf), &
                  NeutralDens_IC(:,iIonoHalf),SIGS,SIGI,SIGA,Qstar, &
                  specup(iLine,:,i),specdn(iLine,:,i),&
                  Qestar_ICI(iLine,iIono,j),Qpstar_ICI(iLine,iIono,j))

             DO k=1,nThetaAlt_II(iLine,i)-1
                call get_cascade_and_lossum(iLine,nNeutral,SIGA,&
                     NeutralDens_IC(:,iIonoHalf),j,iIono,k,cascade,lossum,2)
                CALL ThetaVar1(nAngle,Theta,del1,del2,k, &
                     EqAngleGrid_IG(iLine,0:nAngle))
                Theta=cPi-Theta
                del1=-del1
                del2=-del2
                IF (k.EQ.nThetaAlt_II(iLine,i)-1) del2=-dThetaEnd_II(iLine,i)
                muO=cos(Theta)
                sigmaO=sigO1
                s1=sigmaO
                s2=1.
                Flasti=iphidn(iLine,k,iIono+1,j)
                Flastj=iphidn(iLine,k,iIono,j+1)
                Flastt=liphidn(iLine,k,iIono,j)
                
                !When we are in the ilocal region kill the transport terms in 
                ! coeficients by setting kk to 0
                kk=1
                if (iIonoHalf <= iLocal) kk=0
                CALL NumCalcVar(alpha(k+1),alpha(k),sigma(k+1),sigma(k), &
                     -mu_III(iLine,k,i),coef,muO,del1,del2,beta(i),lbeta(i), &
                     Flasti,Flastj,Flastt,velt,h,delE,Qstar,kk,s1,s2, &
                     cascade,lossum)
                !            CALL CheckWarn(sigma(k+1),warning,4,*9999)
                !            CALL CheckWarn(alpha(k+1),warning,4,*9999)
             enddo

             p=s1/SQRT(ABS(1./coef-1.))

             Flastt=liphidn(iLine,nThetaAlt_II(iLine,i),iIono,j)
             Fsum=iphidn(iLine,nThetaAlt_II(iLine,i)-1,iIono,j)&
                  +iphiup(iLine,nThetaAlt_II(iLine,i)-1,iIono,j)
             Fcheck=iphidn(iLine,nThetaAlt_II(iLine,i),iIono,j)
             Flastj=iphidn(iLine,nThetaAlt_II(iLine,i),iIono,j+1)
             call get_cascade_and_lossum(iLine,nNeutral,SIGA,&
                  NeutralDens_IC(:,iIonoHalf),j,iIono,nThetaAlt_II(iLine,i), &
                  cascade,lossum,2)
             newphi=(Flastt/velt+p*Fsum/(2.*dThetaEnd_II(iLine,i))+s2*(lbeta(i)* &
                  Flastj/delE+Qstar+cascade))&
                  /(1/velt+p/dThetaEnd_II(iLine,i)+s2*(beta(i)/delE+lossum))
             !            CALL CheckWarn(newphi,warning,5,*9999)
             CALL CheckConv(newphi,Fcheck,epsilon,flag)
             iphidn(iLine,nThetaAlt_II(iLine,i),iIono,j)=newphi
             DO k=nThetaAlt_II(iLine,i),1,-1
                newphi=alpha(k)*iphidn(iLine,k,iIono,j)+sigma(k)
                Fcheck=iphidn(iLine,k-1,iIono,j)
                !            CALL CheckWarn(newphi,warning,6,*9999)
                CALL CheckConv(newphi,Fcheck,epsilon,flag)
                iphidn(iLine,k-1,iIono,j)=newphi
             end DO
             CALL midpnt_int(specdn(iLine,j,i),iphidn(iLine,:,iIono,j),&
                  mu_III(iLine,:,i),1, nThetaAlt_II(iLine,i)+1,nAngle+1,1)
             specdn(iLine,j,i)=.5*specdn(iLine,j,i)
             !  Write to a file (if the solution isn't converging)
             IF (count.GT.countmax-2) THEN
                WRITE (10,*) time,j,i,iIono,flag,&
                     (iphidn(iLine,k,iIono,j),k=0,nThetaAlt_II(iLine,i))
             END IF
             if (i==1) then
                h=FieldLineGrid_IC(iLine,i)-FieldLineGrid_IC(iLine,i+1)        
             else
                h=FieldLineGrid_IC(iLine,i-1)-FieldLineGrid_IC(iLine,i)         
             endif
             ! End Downward Region of ionosphere
          end DO FIELDLINE_DOWN
          !put back the iphidn at top of iono2 (which was storing plasmasph bc)
          do k=0,nThetaAlt_II(iLine,nIono)
             iphidn(iLine,k,nIono+1,j)=Fhold(k)
          end do

       END DO ITERATION                  ! End of iteration loop
       CALL CheckFlag(flag,1,warning,-1,*9999)
       !  Save the calculated values for the next time and energy steps
       if (IsIono1) then
          lbeta(1:nIono)=beta(1:nIono)
       else
          lbeta(nPoint-nIono+1:nPoint)=beta(nPoint-nIono+1:nPoint)
       endif
    end DO ENERGY
    
9999 IF (warning.GT.0) THEN
       PRINT *, 'Negative Densities Occurred:',newphi,warning
       PRINT *, time,j,i,k,Bfield_IC(iLine,i),BFieldEq_I(iLine),&
            BFieldEq_I(iLine),FieldLineGrid_IC(iLine,i),h,Theta,muO,&
            mu_III(iLine,k,i),del1, dThetaEnd_II(iLine,i),beta(i),sigmaO,p,count,&
            Qstar
       PRINT *, iphiup(iLine,k,iIono,j),iphidn(iLine,k,iIono,j)
       STOP
    ELSE IF (warning.LT.0) THEN
       PRINT *, 'No convergence: ',count,Time,warning,j
       STOP
    END IF
    !  Format for the nonconvergent output
!!!27  FORMAT (100E10.3,4I4,1P,100E10.3)
    
    RETURN
  END SUBROUTINE update_se_state_iono
  
  !=============================================================================

  !=============================================================================
  !* ------------------------------------------------------------------ **
!  Subroutine NumCalcVar finds alpha(k) and sigma(k).
!*  VARIABLE DESCRIPTIONS
!*      alpha   Flux attenuator for this pitch angle
!*      alp      Flux attenuator for the previous pitch angle
!*      sigma   Flux offset for this angle: [flux]
!*      sig      Flux offset for the previous angle; [flux]
!*      del1    Pitch angle step between this and the previous step; rad
!*      del2    Pitch angle step between this and the next step; rad
!*      muO     Cosine of the equatorial pitch angle
!*      mu      Cosine of the actual pitch angle
!*      BB      Magnetic field strength for this step; G
!*      BFieldEq_I(iLine)      Magnetic field strength at the equator; G
!*      beta    Energy loss parameters; ????
!*      lbeta   Previous energy step's beta; ????
!*      Fj      Flux for previous E for same space, time, and angle step
!*      Fi      Flux for previous space for same angle, time, and E step
!*      Ft      Flux for previous time for same E, space, and angle step
!*      vt      Velocity of an electron at this energy * delt; cm
!*      h       Size of this step; cm
!*      dE      Energy step size; eV
!*      Qstar   Electron production rate(primary plus secondary;
!*                cm-3 s-1 rad-1 ????
!*      p
!*      q
!*      a,b,c,d
!*
  SUBROUTINE NumCalcVar(alpha,alp,sigma,sig,mu,coef,muO,del1, &
       del2,beta,lbeta,Fi,Fj,Ft,vt,h,dE,Qstar,kk,s1,s2,cascade,lossum)
    REAL alpha,alp,sigma,sig,mu,coef,muO,del1,del2,beta,lbeta,s1,s2, &
         Fi,Fj,Ft,vt,h,dE,a,b,c,d,p,q,Qstar,lossum,cascade,kk
    
    !write(*,*) 'muO, BB, BFieldEq_I(iLine)',muO, BB, BFieldEq_I(iLine)
    p=kk*s1*(muO**4+coef-1)/((muO**3)*(1-muO**2)**.5)
    q=kk*s1*ABS((coef-1+muO**2)/muO**2)
    a=(2*q-p*del2)/(del1*(del1+del2))
    b=(2*q+p*del1)/(del2*(del1+del2))
    c=(2*q+(del1-del2)*p)/(del1*del2)+s2*(kk*mu/h+beta/dE+lossum)+1/vt
    d=s2*(lbeta*Fj/dE+kk*mu*Fi/h+Qstar+cascade)+Ft/vt

    IF (b.LT.0.) b=0.
    alpha=b/(c-a*alp)
    sigma=(a*sig+d)/(c-a*alp)
    IF ((d.LT.0).OR.(d-d.NE.0)) THEN
       PRINT *,'NEG D ',p,q,a,b,c,d,sigma,sig,alpha,alp
       PRINT *,mu,h,dE,vt,Fi,Fj,Ft,Qstar,cascade
    ELSE IF ((sigma.LT.0).OR.(sigma-sigma.NE.0)) THEN
       PRINT *,'NEG SIGMA ',p,q,a,b,c,d,sigma,sig,alpha,alp
       PRINT *,mu,h,dE,vt,Fi,Fj,Ft,Qstar
    ELSE IF ((alpha.LT.0).OR.(alpha-alpha.NE.0)) THEN
       PRINT *,'NEG ALPHA ',p,q,a,b,c,d,sigma,sig,alpha,alp
       PRINT *,mu,h,dE,vt,Fi,Fj,Ft,Qstar
    END IF
    !      PRINT 10, BB,BFieldEq_I(iLine)/BB,muO,del1,del2,p,q,a,b,alpha
10  FORMAT (10(1PG11.4,1X))
    RETURN
  END SUBROUTINE NumCalcVar
  !=============================================================================

  !* ------------------------------------------------------------------ **
  !  Subroutine get_sigmaO_and_eprod (was IonoVar1) 
  !*  calculates beta, sigmaO, and Qstar for the ionosphere.
  !*  VARIABLE DESCRIPTIONS
  !*      EnIon   Energies for ionization states for each species; eV
  !*      SIGS    Elastic collision cross sections, species, energy; cm2
  !*      SIGI    Differential cross sections, state,species,sec,pri; cm2
  !*      ionsum  Sum of energy losses due to ionization; ???? [beta]
  !*      elassum Sum of elastic scattering cross sections; ???? [beta]
  !*      exsum   Sum of energy losses due to excitation; ???? [beta]
  !*      prodsum Sum of secondary electron production; cm-2 s-1 eV-1 sr-1
  !*      BINNUM  Function that finds the array position for a given E
  !*      jj,m,n  Energy, state, and species loop counters
  !*      PE      Photoelectron production rate in the ionos.; cm-3 s-1
  !*      pi      3.1415926536
  !*      j      Current energy increment
  !*      dens    Thermal electron density for this step; cm-3
  !*      Qstar   Electron production rate(primary plus secondary;
  !*                cm-3 s-1 rad-1 ????
  !*  FLUX VARIABLES; cm-2 s-1 eV-1 sr-1
  !*      phiup  Flux in the ionospheres (moving from first to second)
  !*      phidn  Flux in the ionospheres (moving from second to first)
  !*
  SUBROUTINE get_sigma0_and_eprod(iEnergyIn,nNeutral,&
       sigO1,PE,NeutralDens_I,SIGS,SIGI,SIGA,Qstar,&
       OmniDirFluxUp_I,OmniDirFluxDn_I,Qe,Qp)
    use ModSeGrid,      ONLY: FieldLineGrid_IC,nIono,nEnergy, nPoint, &
         DeltaE_I,EnergyGrid_I,BINNUM
    use ModMath,        ONLY: midpnt_int
    use ModNumConst,    ONLY: cPi

    IMPLICIT NONE
    integer, intent(in) :: iEnergyIn,nNeutral
    real   , intent(inout) :: sigO1
    real   , intent(in) :: PE,NeutralDens_I(nNeutral)
    !incomming crossections
    real   , intent(in) :: SIGS(nNeutral,nEnergy)
    real   , intent(in) :: SIGI(nNeutral,nEnergy,nEnergy)
    real   , intent(in) :: SIGA(nNeutral,nEnergy,nEnergy)
    !outgoing electron production rate
    real   , intent(out):: Qstar
    
    real   , intent(in) :: OmniDirFluxUp_I(nEnergy),OmniDirFluxDn_I(nEnergy)
    !outgoing electron production rate from secondary production
    real   , intent(out):: Qe
    !incomming electron production rate from precip
    real   , intent(in) :: Qp

    real   :: NetFlux, SecProd(nEnergy)
    REAL   :: elassum
    
    ! Minimum ionization threshold
    real,parameter :: Eplus = 12.0 
    INTEGER j,m,n,LL,jj
    !---------------------------------------------------------------------------
        
    ! this part gives the sum elastic scattering cross sections
    elassum=0.
    Qe=0.
    SecProd(:)=0.

    do n=1,nNeutral
       elassum=elassum+NeutralDens_I(n)*SIGS(n,j)
    end do
    sigO1=sigO1+.5*elassum
    
    ! this calculates the electron total production
    DO  n=1,nNeutral
       IF (2*EnergyGrid_I(iEnergyIn)+Eplus &
            < EnergyGrid_I(nEnergy)+.5*DeltaE_I(nEnergy)) THEN
          LL=BINNUM(2*EnergyGrid_I(iEnergyIn)+Eplus)
          DO jj=LL,nEnergy
             NetFlux=OmniDirFluxUp_I(jj)-OmniDirFluxDn_I(jj)
             SecProd(jj)=SecProd(jj)+NeutralDens_I(n)*SIGI(n,iEnergyIn,jj)*NetFlux
          enddo
       END IF
    enddo
    CALL midpnt_int(Qe,SecProd,DeltaE_I,1,nEnergy,nEnergy,2)
    Qstar=PE/(4*cPi*DeltaE_I(iEnergyIn))+Qe+Qp
    RETURN
  end SUBROUTINE get_sigma0_and_eprod


  !=============================================================================
  !  Subroutine get_cascade_and_lossum (was IonoVar2) finds cascade and lossum.
  !     m       Specifies the use of iphiup or iphidn
  SUBROUTINE get_cascade_and_lossum(iLine,nNeutral,SIGA,NeutralDens_I, &
       iEnergyIn,iIono,iAngle,cascade,lossum,m)

    use ModSeGrid,      ONLY: FieldLineGrid_IC,nIono,nEnergy, nPoint, &
         DeltaE_I,EnergyGrid_I, EnergyMin

    IMPLICIT NONE    
    integer, intent(in) :: iLine, nNeutral
    real   , intent(in) :: SIGA(nNeutral,nEnergy,nEnergy)
    real   , intent(in) :: NeutralDens_I(nNeutral)
    integer, intent(in) :: iEnergyIn,iIono,iAngle
    real   , intent(out):: cascade,lossum
    integer, intent(in) :: m
    
    real :: flux
    integer :: jj,n,LL
    !---------------------------------------------------------------------------
    lossum=0.
    cascade=0.
    
    ! New method (with Swartz fix and cross section splitting)           
    do jj=1,iEnergyIn
       do n=1,nNeutral
          lossum=lossum+NeutralDens_I(n)*SIGA(n,jj,iEnergyIn)
       end do
    end do
    IF (iEnergyIn.LT.nEnergy) THEN
       do jj=iEnergyIn+1,nEnergy
          LL=jj-iEnergyIn
          ! IF (m.EQ.1) flux=AngIonJ(i,j,k,i,jj,iphiup(0,i1,jj))
          ! IF (m.EQ.2) flux=AngIonJ(i,j,k,i,jj,iphidn(0,i1,jj))
          IF (m.EQ.1) flux=iphiup(iLine,iAngle,iIono,jj)
          IF (m.EQ.2) flux=iphidn(iLine,iAngle,iIono,jj)
          do n=1,nNeutral
             cascade=cascade+NeutralDens_I(n)*SIGA(n,LL,jj)*flux
          end do
       end do
    END IF
    
    RETURN
  end SUBROUTINE get_cascade_and_lossum


  !=============================================================================
  !* ------------------------------------------------------------------ **
  !*  Subroutine IonoVar3 sets collisional energy loss terms to zero,
  !*  the plasmasperic case.
  SUBROUTINE IonoVar3(Qstar,lossum,cascade)
    real, intent(out) :: Qstar,lossum,cascade
    Qstar=0.
    lossum=0.
    cascade=0.
  END SUBROUTINE IonoVar3

  !=============================================================================
  !* ------------------------------------------------------------------ **
  !*  Finds beta and sigO, the coulomb collision coefficients.
  !NOTE: CHECK HOW SRC is defined
  SUBROUTINE CoulVar(i,beta,sigO1,TE,ZE,KE)
    use ModMath, only: G,erf
    REAL beta,sigO1,TE,ZE,KE,A,MC,X
    INTEGER n,i                              !,Iterm1,Iterm2
    
    DATA A/2.6E-12/, MC/1837./
    !write(*,*) 'KE, TE',KE, TE
    X=SQRT(KE/TE)
    beta=ZE/TE*G(X)                  ! ion infl. --> zero
    sigO1=ZE*(1.+erf(X)-G(X))            ! ion infl. = ZE (thus the 1.
    do n=1,IRC                  ! Ring current influence
       X=SQRT(KE*MC*MRC(n)/TRC(n))
       beta=beta+SRC(i)*NRC(n)/TRC(n)*G(X)
       sigO1=sigO1+SRC(i)*NRC(n)*(erf(X)-G(X))
    end do
    beta=2.*A*beta
    sigO1=.25*A*sigO1/(KE*KE)
    RETURN
  END SUBROUTINE CoulVar

  !=============================================================================
  !* ------------------------------------------------------------------ **
  !  Subroutine ThetaVar1 finds Theta, del1, and del2.
  SUBROUTINE ThetaVar1(nAngle,Theta,del1,del2,k,grid)
    INTEGER, intent(in) :: nAngle
    REAL Theta,del1,del2,grid(0:nAngle)
    INTEGER k,kk
    !*
    kk=k
    Theta=grid(kk)
    del1=grid(kk)-grid(kk-1)
    del2=grid(kk+1)-grid(kk)
    RETURN
  END SUBROUTINE ThetaVar1
  !=============================================================================
  
  !* ------------------------------------------------------------------ **
  !  Subroutine CheckConv sees if a flux has converged or not.
  !*  VARIABLE DESCRIPTIONS
  !*      epsil   Convergence parameter in the iteration loop
  !*      flag    Indicates whether the fluxes have converged
  !*  FLUX VARIABLES; cm-2 s-1 eV-1 sr-1
  !*        flux      Current flux value to be checked
  !*      oldflux      Previous value to check against
  !*
  SUBROUTINE CheckConv(flux,oldflux,epsil,flag, DoReportError, ReturnError)
    REAL   ,  intent(in) :: flux,oldflux,epsil
    INTEGER, intent(out) :: flag
    LOGICAL,optional,intent(in) :: DoReportError
    real,optional,intent(out) :: ReturnError

    real :: error
    !--------------------------------------------------------------------------
    
    IF (flux.GT.1E-20) THEN
       error=ABS(flux-oldflux)/flux
       IF (error.GT.epsil) flag=1
    END IF
    if (present(DoReportError) .and. DoReportError) then
       write(*,*) 'Error is:',error,'Epsilon is:',epsil
    endif
    
    if (present(ReturnError)) ReturnError=error
    RETURN
  END SUBROUTINE CheckConv
  !=============================================================================
  
  !* ------------------------------------------------------------------ **
  !  Subroutine CheckFlag sees if the loop  ended due to convergence or
  !  not.
  !*  VARIABLE DESCRIPTIONS
  !*      flag    Indicates whether the fluxes have converged
  !*      check      The number to check flag against
  !*      warn    Indicates (nonzero) negative densities have occurred
  !*      mark      Value warn is set to if check fails (for diagnostics)
  !*
  SUBROUTINE CheckFlag(flag,check,warn,mark,*)
    INTEGER flag,check,warn,mark
    IF (flag.EQ.check) THEN
       warn=mark
       PRINT *, 'Flag=',flag
       RETURN 1
    END IF
    RETURN
  END SUBROUTINE CheckFlag

  !=============================================================================
  !* ------------------------------------------------------------------ **
  !*  Pass fluxes to "last time step fluxes" and set flux arrays to zero.
  SUBROUTINE InitPlas(iLine,DoSavePreviousAndReset)
    USE ModSeGrid, only: nEnergy, nPlas,nIono,nThetaAlt_II
    IMPLICIT NONE

    integer, intent(in) :: iLine
    logical, intent(in) :: DoSavePreviousAndReset
    real :: c
    INTEGER iPlas,j,iAngle,iAlt

    ! Set upper energy boundary condition on the flux. 
    ! Set fluxes at nEnergy+1 to fluxes at nEnergy times a constant <= 1 
    do iPlas=1,nPlas
       iAlt=iPlas+nIono
       do iAngle=0,nThetaAlt_II(iLine,iAlt)
          
          ! set bc for fluxes up
          c=0.
          IF (phiup(iLine,iAngle,iPlas,nEnergy-1).GT.0.) &
               c=min(1., .75*phiup(iLine,iAngle,iPlas,nEnergy)&
               /phiup(iLine,iAngle,iPlas,nEnergy-1))
         
          phiup(iLine,iAngle,iPlas,nEnergy+1)=&
               phiup(iLine,iAngle,iPlas,nEnergy)*c

          ! set bc for fluxes down
          c=0.
          IF (phidn(iLine,iAngle,iPlas,nEnergy-1).GT.0.) &
               c=min(1., .75*phidn(iLine,iAngle,iPlas,nEnergy)&
               /phidn(iLine,iAngle,iPlas,nEnergy-1))
         
          phidn(iLine,iAngle,iPlas,nEnergy+1)=&
               phidn(iLine,iAngle,iPlas,nEnergy)*c
       end do
    end do
    
    ! Only update nEnergy+1 values?
    if (.not.DoSavePreviousAndReset) RETURN

!*  Move current time fluxes to previous time fluxes
    lphiup(iLine,:,1:nPlas,1:nEnergy) = phiup(iLine,:,1:nPlas,1:nEnergy)
    lphidn(iLine,:,1:nPlas,1:nEnergy) = phidn(iLine,:,1:nPlas,1:nEnergy)

!*  Reset current time fluxes
    phiup(iLine,:,:,:)=0.0
    phidn(iLine,:,:,:)=0.0

    RETURN
  END SUBROUTINE InitPlas
  !=============================================================================
  !* ------------------------------------------------------------------ **
  !*  Pass fluxes to "last time step fluxes" and set flux arrays to zero.
  SUBROUTINE InitIono(iLine,DoSavePreviousAndReset)
    USE ModSeGrid, only: nEnergy, nPlas,nIono,nThetaAlt_II
    IMPLICIT NONE

    integer, intent(in) :: iLine
    logical, intent(in) :: DoSavePreviousAndReset
    real :: c
    INTEGER iIono,j,iAngle,iAlt

    ! Set upper energy boundary condition on the flux. 
    ! Set fluxes at nEnergy+1 to fluxes at nEnergy times a constant <= 1 
    do iIono=1,nIono
       iAlt=iIono
       do iAngle=0,nThetaAlt_II(iLine,iAlt)
          
          ! set bc for fluxes up
          c=0.
          IF (iphiup(iLine,iAngle,iIono,nEnergy-1).GT.0.) c=min(1., &
               .75*iphiup(iLine,iAngle,iIono,nEnergy) &
               /iphiup(iLine,iAngle,iIono,nEnergy-1))
          iphiup(iLine,iAngle,iIono,nEnergy+1)=&
               iphiup(iLine,iAngle,iIono,nEnergy)*c

          ! set bc for fluxes down
          c=0.
          IF (iphidn(iLine,iAngle,iIono,nEnergy-1).GT.0.) c=min(1., &
               .75*iphidn(iLine,iAngle,iIono,nEnergy) &
               /iphidn(iLine,iAngle,iIono,nEnergy-1))
          iphidn(iLine,iAngle,iIono,nEnergy+1)=&
               iphidn(iLine,iAngle,iIono,nEnergy)*c

       end do
    end do
    
    ! Only update nEnergy+1 values?
    if (.not.DoSavePreviousAndReset) RETURN

!*  Move current time fluxes to previous time fluxes
    liphiup(iLine,:,1:2*nIono,1:nEnergy) = iphiup(iLine,:,1:2*nIono,1:nEnergy)
    liphidn(iLine,:,1:2*nIono,1:nEnergy) = iphidn(iLine,:,1:2*nIono,1:nEnergy)

!*  Reset current time fluxes
    iphiup(iLine,:,:,:)=0.0
    iphidn(iLine,:,:,:)=0.0

    RETURN
  END SUBROUTINE InitIono

  !============================================================================
  !  This subroutine checks for "time" convergence for a particular line
  SUBROUTINE check_time(iLine,flag)
    use ModSeGrid, only: nAngle, nPlas, nIono, nEnergy, nPoint, nLine,&
         nThetaAlt_II

    integer,intent(in) :: iLine
    integer,intent(out):: flag
    integer :: iEnergy, iAngle, iAlt, iIono,iPlas
    integer ibad
    real    :: flux, oldflux
    
    real ::   error, ErrorMax, FluxAtMax, OldFluxAtMax
    integer:: iAltAtMax,iAngleAtMAx,iEnergyAtMax
    !--------------------------------------------------------------------------

    Ibad=0
    flag=0
    errormax=0.0
    ! loop over all energy, altitudes, and angles to evaluate convergence
    do iEnergy=nEnergy,1,-1
       do iAlt=1,nPoint
          do iAngle=0,nThetaAlt_II(iLine,iAlt)
             !choose what region to set the flux from
             IF (iAlt <= nIono) THEN
                !ionosphere 1
                flux=iphiup(iLine,iAngle,iAlt,iEnergy)
                oldflux=liphiup(iLine,iAngle,iAlt,iEnergy)
             elseif(iAlt >nIono .and. iAlt <=nPoint-nIono) then
                !plasmasphere
                iPlas=iAlt-nIono
                flux=phiup(iLine,iAngle,iPlas,iEnergy)
                oldflux=lphiup(iLine,iAngle,iPlas,iEnergy)
             else
                !ionosphere 2
                iIono=iAlt-nPlas
                flux=iphiup(iLine,iAngle,iIono,iEnergy)
                oldflux=liphiup(iLine,iAngle,iIono,iEnergy)
             endif
             
             ! check the convergence
             call CheckConv(flux,oldflux,epsilon,flag,DoReportError=.false.,&
                  ReturnError=error)
             
             !find maximum error and associated values to report at end of check
             ErrorMax = max(error,errormax)
             ! If a new errormax update the associated iAltAtMax and fluxes
             if (ErrorMax==error) then
                iAltAtMax   = iAlt
                iAngleAtMax = iAngle
                iEnergyAtMax= iEnergy
                FluxAtMax   = flux
                OldFluxAtMax= oldflux
             endif

             if (flag ==1 .and. oldflux > 0.0) Ibad=Ibad+1
             
          end do
       end do
    end do
    !report results of check_time
    write(*,*), 'FINISHED check_time, Report:'
    write(*,*), 'check_time: Ibad=',Ibad
    write(*,*), 'check_time: ErrorMax=',ErrorMax
    write(*,*), 'check_time: iAltAtMax=',iAltAtMax
    write(*,*), 'check_time: iAngleAtMax=',iAngleAtMax
    write(*,*), 'check_time: nThetaAtMax=',nThetaAlt_II(iLine,iAltAtMax)
    write(*,*), 'check_time: iEnergyAtMax=',iEnergyAtMax
    write(*,*), 'check_time: FluxAtMax=',FluxAtMax
    write(*,*), 'check_time: OldFluxMax=',OldFluxAtMax
    RETURN
  END SUBROUTINE check_time
  
  !============================================================================
  subroutine allocate_state_arrays
    use ModSeGrid, only: nAngle, nPlas, nIono, nEnergy, nPoint, nLine
    
    
    if(.not.allocated(phiup)) allocate(phiup(nLine,0:nAngle,0:nPlas+1,nEnergy+1))
    if(.not.allocated(phidn)) allocate(phidn(nLine,0:nAngle,0:nPlas+1,nEnergy+1))

    if(.not.allocated(iphiup)) &
         allocate(iphiup(nLine,0:nAngle,0:2*nIono,nEnergy+1))
    if(.not.allocated(iphidn)) &
         allocate(iphidn(nLine,0:nAngle,0:2*nIono,nEnergy+1))

    if(.not.allocated(specup)) allocate(specup(nLine,nEnergy,nPoint))
    if(.not.allocated(specdn)) allocate(specdn(nLine,nEnergy,nPoint))

    if(.not.allocated(lphiup)) allocate(lphiup(nLine,0:nAngle,nPlas,nEnergy))
    if(.not.allocated(lphidn)) allocate(lphidn(nLine,0:nAngle,nPlas,nEnergy))

    if(.not.allocated(liphiup))allocate(liphiup(nLine,0:nAngle,2*nIono,nEnergy))
    if(.not.allocated(liphidn))allocate(liphidn(nLine,0:nAngle,2*nIono,nEnergy))

    if(.not.allocated(SRC))    allocate(SRC(nPoint))

    if(.not.allocated(Qestar_ICI))allocate(Qestar_ICI(nLine,2*nIono,nEnergy))
    if(.not.allocated(Qpstar_ICI))allocate(Qpstar_ICI(nLine,2*nIono,nEnergy))
  end subroutine allocate_state_arrays

  !============================================================================
  ! UNIT test for SE update states
  subroutine se_update_state_test
    use ModSeBackground
    use ModSeGrid, only:create_se_test_grid,nLine,nPoint,nIono,nPlas
    use ModSePlot, only: plot_state
    
    integer :: iLine=1, flag=1
    logical :: DoSavePreviousAndReset = .true.
    integer :: nStep
    logical :: IsOpen
    !--------------------------------------------------------------------------

    ! First set up the grid that we will update the state in (this is the same 
    ! as the unit test for the grid).
    write(*,*) 'creating grid'
    call create_se_test_grid

    ! Allocate the state arrays
    write(*,*) 'allocating state arrays'
    call allocate_state_arrays
    
    ! Allocate the background right
    write(*,*) 'allocating background arrays'
    call allocate_background_arrays
    
    ! Fill the background arrays
!    eThermalDensity_IC(:,:) = 0.00001
    eThermalDensity_IC(:,:) = 1.0e2
    eThermalTemp_IC(:,:) = 1.0
    
    ! Define the initial state in the ionosphere
    do iLine=1,nLine
       iphiup(iLine,:,:,:)=1.0e5
       iphidn(iLine,:,:,:)=1.0e5

       liphiup(iLine,:,:,:)=1.0e5
       liphidn(iLine,:,:,:)=1.0e5

       phiup(iLine,:,:,:)=0.00001
       phidn(iLine,:,:,:)=0.00001
    end do
    
    ! Set the timestep and convergence criteria
    delt=1.0e5
    epsilon = 0.4

    write(*,*) 'Starting Time loop'
    nStep = 0

    ! plot initial state
    call plot_state(1,nStep,time,iphiup,iphidn,phiup,phidn)

    IsOpen = .true.
    TIME_LOOP: do while (flag == 1)
       ! Initialize the plasmasphere
       write(*,*) 'Initializing plasmasphere'
       call initplas(1,DoSavePreviousAndReset)
       
       ! update the SE state
       write(*,*) 'update se state'
       call update_se_state(1, IsOpen)
       
       ! check convergence
       write(*,*) 'check for convergence'
       call check_time(1,flag)
       
       ! increment step
       nStep=nStep+1

       ! plot output
       call plot_state(1,nStep,time,iphiup,iphidn,phiup,phidn)
       

    end do TIME_LOOP
    
  end subroutine se_update_state_test
end Module ModSeState
