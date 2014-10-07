Module ModSeState
  implicit none
  
  private !except
  
  real, public :: Time
  !\
  ! FLUX VARIABLES; cm-2 s-1 eV-1 sr-1
  !/
  
  ! Flux in the plasmasphere from first to second (up) ionosphere 
  real, public,allocatable :: phiup(:,:,:),phidn(:,:,:)
  
  ! Flux in the ionosphere from first to second (up) ionosphere 
  real, public,allocatable :: iphiup(:,:,:),iphidn(:,:,:)
  
  ! Flux from previous time step
  real, public,allocatable :: lphiup(:,:,:),lphidn(:,:,:)
  real, public,allocatable :: liphiup(:,:,:),liphidn(:,:,:)
  
  !Omnidirectional flux for upward/downward-directed hemisphere
  real, public,allocatable :: specup(:,:),specdn(:,:)
  
  real, public :: delt,epsilon

  ! Ring current parameters for coulomb collisions
  real :: NRC,TRC,MRC,IRC           !,Iterm1,Iterm2
  real,allocatable :: SRC(:)
  integer,parameter :: NRing=6 ! ring current Maxwellian fits
contains
  SUBROUTINE update_se_state(iLine)
    !*  This subroutine used to be the main program, until other operators
    !*  added (3/23/95). Now it is one of several operators called by the
    !*  driver program, 'stet1.f' (Super_Thermal_electron_Transport_Model).
    
    USE ModSeGrid, only: nEnergy, EnergyGrid_I, DeltaE_I, EnergyMin, EnergyMax,&
         EqAngleGrid_IG,nThetaAlt_II,dTheta_II,mu_II,&
         FieldLineGrid_IC,nTop,Bfield_IC, BFieldIono_I, &
         nIono, nPlas, nPoint,nAngle
    
    USE ModSeBackground, only: eThermalDensity_IC,eThermalTemp_IC

    use ModMath, only: midpnt_int
    
    use ModNumConst,    ONLY: cPi

    IMPLICIT NONE
!    INCLUDE 'numbers.h'
    integer, intent(in) :: iLine
    REAL h,p,Fcheck, sigmaO,muO,velt,coef, lbeta(nPoint),alpha(nAngle), &
         sigma(nAngle),Flastj,kk, del1,del2, newphi,beta(nPoint),sigO1, &
         Fsum, delE,Flasti, Flastt,Qstar, Theta,s1,s2
    INTEGER i,j,k,i1,jj,Ist, &
         flag,count,warning,SPick,space,ii

    !set maximum iterations
    integer, parameter :: countmax=300

    real :: cascade, lossum
    
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
             phiup(k,0,j)=iphiup(k,nIono,j)
          end do
          h=FieldLineGrid_IC(iLine,nIono+1)-FieldLineGrid_IC(iLine,nIono)
          FIELDLINE_UPWARD: DO i=nIono+1,nIono+nPlas
             CALL CoulVar(i,beta(i),sigO1,eThermalTemp_IC(iLine,i),eThermalDensity_IC(iLine,i),EnergyGrid_I(j))
             !KLUDGE kill pitchangle scattering
             !sigO1=0.0
             !END KLUDGE
             i1=i-nIono
             coef=BFieldIono_I(iLine)/Bfield_IC(iLine,i)
             CALL IonoVar3(Qstar, cascade, lossum)
             do k=1,nThetaAlt_II(iLine,i)-1
                CALL ThetaVar1(nAngle,Theta,del1,del2,k,EqAngleGrid_IG(iLine,0:nAngle))
                !write(*,*) 'k,i,ko(i),angles,Theta*180.0/3.14159', k,i,
                IF (k.EQ.nThetaAlt_II(iLine,i)-1) del2=dTheta_II(iLine,i)
                muO=cos(Theta)
                sigmaO=sigO1
                s1=sigmaO
                s2=1.
                Flasti=phiup(k,i1-1,j)
                Flastj=phiup(k,i1,j+1)
                Flastt=lphiup(k,i1,j)
                !write(*,*)' i,j,k', i,j,k
                CALL NumCalcVar(alpha(k+1),alpha(k),sigma(k+1),sigma(k), &
                     mu_II(k,i),coef,muO,del1,del2,beta(i),lbeta(i), &
                     Flasti,Flastj,Flastt,velt,h,delE,Qstar,kk,s1,s2,&
                     cascade,lossum)
                !            CALL CheckWarn(sigma(k+1),warning,1,*9999)
                !            CALL CheckWarn(alpha(k+1),warning,1,*9999)
             end do
             IF (i.EQ.nTop) THEN
                p=-2.*s1/dTheta_II(iLine,i)
             ELSE
                p=-s1/SQRT(ABS(1./coef-1.))
             END IF
             Flastt=lphiup(nThetaAlt_II(iLine,i),i1,j)
             Fsum=phiup(nThetaAlt_II(iLine,i)-1,i1,j)+phidn(nThetaAlt_II(iLine,i)-1,i1,j)
             Fcheck=phiup(nThetaAlt_II(iLine,i),i1,j)
             Flastj=phiup(nThetaAlt_II(iLine,i),i1,j+1)
             newphi=(Flastt/velt-p*Fsum/(2.*dTheta_II(iLine,i))+s2*(lbeta(i)* &
                  Flastj/delE+Qstar+cascade))/ &
                  (1/velt-p/dTheta_II(iLine,i)+s2*(beta(i)/delE+lossum))
             !            CALL CheckWarn(newphi,warning,2,*9999)
             CALL CheckConv(newphi,Fcheck,epsilon,flag)
             phiup(nThetaAlt_II(iLine,i),i1,j)=newphi
             do k=nThetaAlt_II(iLine,i),1,-1
                newphi=alpha(k)*phiup(k,i1,j)+sigma(k)
                Fcheck=phiup(k-1,i1,j)
                !            CALL CheckWarn(newphi,warning,3,*9999)
                CALL CheckConv(newphi,Fcheck,epsilon,flag)
                phiup(k-1,i1,j)=newphi
             end do
             CALL midpnt_int(specup(j,i),phiup(0,i1,j),mu_II(0,i),1, &
                  nThetaAlt_II(iLine,i)+1,nAngle+1,1)
             specup(j,i)=-.5*specup(j,i)
             !  Write to a file (if the solution isn't converging)
             IF (count.GT.countmax-2) THEN
                WRITE (10,27) time,j,i,i1,nThetaAlt_II(iLine,i)/4,&
                     flag,(phiup(k,i1,j),k=0,nThetaAlt_II(iLine,i))
                WRITE (10,27) time,j,i,i1,nThetaAlt_II(iLine,i),flag,&
                     (phiup(k,i1,j),k=nThetaAlt_II(iLine,i)/2+1, &
                     nThetaAlt_II(iLine,i),2)
             END IF
             h=FieldLineGrid_IC(iLine,i+1)-FieldLineGrid_IC(iLine,i)            ! End Upward Region #2 loop
          end DO FIELDLINE_UPWARD
          !*  Downward Region #2: the plasmasphere
          do k=0,nThetaAlt_II(iLine,nIono)
             phidn(k,nPlas+1,j)=iphidn(k,nIono+1,j)
          end do
          h=FieldLineGrid_IC(iLine,nIono+nPlas)-FieldLineGrid_IC(iLine,nIono+nPlas+1)
          FIELDLINE_DOWN: DO i=nIono+nPlas,nIono+1,-1
             CALL CoulVar(i,beta(i),sigO1,eThermalTemp_IC(iLine,i),eThermalDensity_IC(iLine,i),EnergyGrid_I(j))
             !KLUDGE kill pitchangle scattering
             !sigO1=0.0
             !END KLUDGE
             i1=i-nIono
             coef=BFieldIono_I(iLine)/Bfield_IC(iLine,i)
             CALL IonoVar3(Qstar, cascade, lossum)
             DO k=1,nThetaAlt_II(iLine,i)-1
                CALL ThetaVar1(nAngle,Theta,del1,del2,k,EqAngleGrid_IG(iLine,0:nAngle))
                Theta=cPi-Theta
                del1=-del1
                del2=-del2
                IF (k.EQ.nThetaAlt_II(iLine,i)-1) del2=-dTheta_II(iLine,i)
                muO=cos(Theta)
                sigmaO=sigO1
                s1=sigmaO
                s2=1.
                Flasti=phidn(k,i1+1,j)
                Flastj=phidn(k,i1,j+1)
                Flastt=lphidn(k,i1,j)
                CALL NumCalcVar(alpha(k+1),alpha(k),sigma(k+1),sigma(k), &
                     -mu_II(k,i),coef,muO,del1,del2,beta(i),lbeta(i), &
                     Flasti,Flastj,Flastt,velt,h,delE,Qstar,kk,s1,s2, &
                     cascade,lossum)
                !            CALL CheckWarn(sigma(k+1),warning,4,*9999)
                !            CALL CheckWarn(alpha(k+1),warning,4,*9999)
             enddo
             IF (i.EQ.nTop) THEN
                p=2.*s1/dTheta_II(iLine,i)
             ELSE
                p=s1/SQRT(ABS(1./coef-1.))
             END IF
             Flastt=phidn(nThetaAlt_II(iLine,i),i1,j)
             Fsum=phidn(nThetaAlt_II(iLine,i)-1,i1,j)+phiup(nThetaAlt_II(iLine,i)-1,i1,j)
             Fcheck=phidn(nThetaAlt_II(iLine,i),i1,j)
             Flastj=phidn(nThetaAlt_II(iLine,i),i1,j+1)
             newphi=(Flastt/velt+p*Fsum/(2.*dTheta_II(iLine,i))+s2*(lbeta(i)* &
                  Flastj/delE+Qstar+cascade))/ &
                  (1/velt+p/dTheta_II(iLine,i)+s2*(beta(i)/delE+lossum))
             !            CALL CheckWarn(newphi,warning,5,*9999)
             CALL CheckConv(newphi,Fcheck,epsilon,flag)
             phidn(nThetaAlt_II(iLine,i),i1,j)=newphi
             DO k=nThetaAlt_II(iLine,i),1,-1
                newphi=alpha(k)*phidn(k,i1,j)+sigma(k)
                Fcheck=phidn(k-1,i1,j)
                !            CALL CheckWarn(newphi,warning,6,*9999)
                CALL CheckConv(newphi,Fcheck,epsilon,flag)
                phidn(k-1,i1,j)=newphi
             end DO
             CALL midpnt_int(specdn(j,i),phidn(0,i1,j),mu_II(0,i),1, &
                  nThetaAlt_II(iLine,i)+1,nAngle+1,1)
             specdn(j,i)=.5*specdn(j,i)
             !  Write to a file (if the solution isn't converging)
             IF (count.GT.countmax-2) THEN
                WRITE (10,27) time,j,i,i1,flag,&
                     (phidn(k,i1,j),k=0,nThetaAlt_II(iLine,i)/2,2)
                WRITE (10,27) time,j,i,i1,flag,&
                     (phidn(k,i1,j),k=nThetaAlt_II(iLine,i)/2+1, &
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
       PRINT *, time,j,i,k,Bfield_IC(iLine,i),BFieldIono_I(iLine),&
            BFieldIono_I(iLine),FieldLineGrid_IC(iLine,i),h,Theta,muO,mu_II(k,i),&
            del1, dTheta_II(iLine,i),beta(i),sigmaO,p,count,Qstar
       PRINT *, phiup(k,i1,j),phidn(k,i1,j)
       STOP
    ELSE IF (warning.LT.0) THEN
       PRINT *, 'No convergence: ',count,Time,warning,j
       STOP
    END IF
    !  Format for the nonconvergent output
27  FORMAT (5I4,1P,100E10.3)
    
    RETURN
  END SUBROUTINE update_se_state
  
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
!*      BFieldIono_I(iLine)      Magnetic field strength at the equator; G
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
    
    !write(*,*) 'muO, BB, BFieldIono_I(iLine)',muO, BB, BFieldIono_I(iLine)
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
    !      PRINT 10, BB,BFieldIono_I(iLine)/BB,muO,del1,del2,p,q,a,b,alpha
10  FORMAT (10(1PG11.4,1X))
    RETURN
  END SUBROUTINE NumCalcVar


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
    USE ModSeGrid, only: nPoint
    EXTERNAL G,erf
    REAL beta,sigO1,TE,ZE,KE,G,MRC(NRing),NRC(NRing),TRC(NRing),A,MC, &
         X,erf
    INTEGER n,IRC,i                              !,Iterm1,Iterm2
    
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
    kk=k+1            ! Because grid should start at 0, not 1
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
  SUBROUTINE CheckConv(flux,oldflux,epsil,flag)
    REAL   ,  intent(in) :: flux,oldflux,epsil
    INTEGER, intent(out) :: flag
    
    real :: error
    !--------------------------------------------------------------------------
    
    IF (flux.GT.1E-20) THEN
       error=ABS(flux-oldflux)/flux
       IF (error.GT.epsil) flag=1
    END IF
    RETURN
  END SUBROUTINE CheckConv
  !=============================================================================
  
  subroutine allocate_state_arrays
    use ModSeGrid, only: nAngle, nPlas, nIono, nEnergy, nPoint
    
    
    if(.not.allocated(phiup)) allocate(phiup(0:nAngle,0:nPlas,nEnergy))
    if(.not.allocated(phidn)) allocate(phidn(0:nAngle,0:nPlas,nEnergy))

    if(.not.allocated(iphiup)) allocate(iphiup(0:nAngle,0:2*nIono,nEnergy))
    if(.not.allocated(iphidn)) allocate(iphidn(0:nAngle,0:2*nIono,nEnergy))

    if(.not.allocated(specup)) allocate(specup(nEnergy,nPoint))
    if(.not.allocated(specdn)) allocate(specdn(nEnergy,nPoint))

    if(.not.allocated(lphiup)) allocate(lphiup(0:nAngle,nPlas,nEnergy))
    if(.not.allocated(lphidn)) allocate(lphidn(0:nAngle,nPlas,nEnergy))    

    if(.not.allocated(liphiup)) allocate(liphiup(0:nAngle,2*nIono,nEnergy))
    if(.not.allocated(liphidn)) allocate(liphidn(0:nAngle,2*nIono,nEnergy))    

    if(.not.allocated(SRC))    allocate(SRC(nPoint))
  end subroutine allocate_state_arrays

  !============================================================================
  ! UNIT test for SE update states
  subroutine se_update_state_test
    use ModSeBackground
    Use ModSeGrid, only:create_se_test_grid,nLine,nPoint

    ! First set up the grid that we will update the state in (this is the same 
    ! as the unit test for the grid).
    call create_se_test_grid

    ! Allocate the state arrays
    call allocate_state_arrays
    
    ! Allocate the background right
    call allocate_background_arrays(nLine,nPoint)
    
    ! Fill the background arrays
    eThermalDensity_IC(:,:) = 0.00001
    eThermalTemp_IC(:,:) = 1.0
    
    ! Define the initial state in the ionosphere
    iphiup(:,:,:)=1.0
    iphidn(:,:,:)=1.0

    ! Set the timestep and convergence criteria
    delt=1.0e5
    epsilon = 0.01

    ITERATION_LOOP: do

       ! update the SE state
       call update_se_state(1)

    end do ITERATION_LOOP

  end subroutine se_update_state_test
end Module ModSeState
