Module ModSePlot
  implicit none
  
  private !except
  
  public :: plot_state

contains
  !============================================================================
  ! save state plot for verification
  subroutine plot_state(iLine,nStep,time,iphiup,iphidn,phiup,phidn)
    use ModSeGrid,     ONLY: FieldLineGrid_IC,EqAngleGrid_IG,Bfield_IC,&
                             nIono,nPlas, nAngle, nEnergy, nLine, nPoint,&
                             nThetaAlt_II
    use ModIoUnit,     ONLY: UnitTmp_
    use ModPlotFile,   ONLY: save_plot_file
    use ModNumConst,   ONLY: cRadToDeg,cPi

    integer, intent(in) :: iLine, nStep
    real,    intent(in) :: time
    real,    intent(in) :: iphiup(nLine,0:nAngle,0:2*nIono,nEnergy+1),&
                           iphidn(nLine,0:nAngle,0:2*nIono,nEnergy+1), &
                           phiup(nLine,0:nAngle,0:nPlas+1,nEnergy+1), &
                           phidn(nLine,0:nAngle,0:nPlas+1,nEnergy+1)
    real, allocatable   :: Coord_DII(:,:,:), PlotState_IIV(:,:,:)
    real, parameter     :: rEarthCM = 6375.0e5
    !grid parameters
    integer, parameter :: nDim =2, nVar=7, S_=1, PA_=2,B_=1
    integer, parameter :: E1_=3, E2_=4, E3_=5, E4_=6, E5_=7
    !set the corresponding energy channels for E1-E5
    integer, parameter :: EE1_=5, EE2_=10, EE3_=20, EE4_=40, EE5_=80
    character(len=100),parameter :: NamePlotVar='S PA B PA E1 E2 E3 E4 E5 g r'
    character(len=13) :: NamePlot='StatePlot.out'
    character(len=*),parameter :: NameHeader='SE output'
    character(len=5) :: TypePlot='ascii'
    integer :: iAngle,iAngleDn,iPoint,iIono,iPlas
    !--------------------------------------------------------------------------
    allocate(Coord_DII(nDim,nPoint,2*nAngle),PlotState_IIV(nPoint,2*nAngle,nVar))

!    do iLine=1,nLine
       PlotState_IIV = 0.0
       Coord_DII     = 0.0
       
       !Set Coordinates along field line and PA
       do iPoint=1,nPoint
          do iAngle=1,2*nAngle
             !set coord based on up or down region
             if(iAngle<nAngle) then
                Coord_DII(S_,iPoint,iAngle) = FieldLineGrid_IC(iLine,iPoint)&
                     /rEarthCM
                Coord_DII(PA_,iPoint,iAngle)= EqAngleGrid_IG(iLine,iAngle)
             else
                iAngleDn = 2*nAngle-iAngle
                Coord_DII(S_,iPoint,iAngle) = FieldLineGrid_IC(iLine,iPoint)&
                     /rEarthCM
                Coord_DII(PA_,iPoint,iAngle)= &
                     cPi-EqAngleGrid_IG(iLine,iAngleDn)
             endif
             
             !set plotstate based on up or down region  
             if (iAngle <= nThetaAlt_II(iLine,iPoint))then
                PlotState_IIV(iPoint,iAngle,B_)  = Bfield_IC(iLine,iPoint)
                PlotState_IIV(iPoint,iAngle,PA_) = EqAngleGrid_IG(iLine,iAngle)
                !choose phiup or iphiup by spatial region
                if (iPoint <= nIono)then
                   iIono=iPoint
                   PlotState_IIV(iPoint,iAngle,E1_) = &
                        iphiup(iLine,iAngle,iIono,EE1_)
                   PlotState_IIV(iPoint,iAngle,E2_) = &
                        iphiup(iLine,iAngle,iIono,EE2_)
                   PlotState_IIV(iPoint,iAngle,E3_) = &
                        iphiup(iLine,iAngle,iIono,EE3_)
                   PlotState_IIV(iPoint,iAngle,E4_) = &
                        iphiup(iLine,iAngle,iIono,EE4_)
                   PlotState_IIV(iPoint,iAngle,E5_) = &
                        iphiup(iLine,iAngle,iIono,EE5_)
                elseif(iPoint >nIono .and. iPoint <=nPoint-nIono) then
                   iPlas=iPoint-nIono
                   PlotState_IIV(iPoint,iAngle,E1_) = &
                        phiup(iLine,iAngle,iPlas,EE1_)
                   PlotState_IIV(iPoint,iAngle,E2_) = &
                        phiup(iLine,iAngle,iPlas,EE2_)
                   PlotState_IIV(iPoint,iAngle,E3_) = &
                        phiup(iLine,iAngle,iPlas,EE3_)
                   PlotState_IIV(iPoint,iAngle,E4_) = &
                        phiup(iLine,iAngle,iPlas,EE4_)
                   PlotState_IIV(iPoint,iAngle,E5_) = &
                        phiup(iLine,iAngle,iPlas,EE5_)
                else
                   iIono=iPoint-nPlas
                   PlotState_IIV(iPoint,iAngle,E1_) = &
                        iphiup(iLine,iAngle,iIono,EE1_)
                   PlotState_IIV(iPoint,iAngle,E2_) = &
                        iphiup(iLine,iAngle,iIono,EE2_)
                   PlotState_IIV(iPoint,iAngle,E3_) = &
                        iphiup(iLine,iAngle,iIono,EE3_)
                   PlotState_IIV(iPoint,iAngle,E4_) = &
                        iphiup(iLine,iAngle,iIono,EE4_)
                   PlotState_IIV(iPoint,iAngle,E5_) = &
                        iphiup(iLine,iAngle,iIono,EE5_)
                endif
             elseif(iAngle >= 2*nAngle-nThetaAlt_II(iLine,iPoint)) then
                iAngleDn = 2*nAngle-iAngle
                PlotState_IIV(iPoint,iAngle,B_)  = Bfield_IC(iLine,iPoint)
                PlotState_IIV(iPoint,iAngle,PA_) = &
                     cPi-EqAngleGrid_IG(iLine,iAngleDn)
                !choose phidn or iphidn by spatial region
                if (iPoint <= nIono)then
                   iIono=iPoint
                   PlotState_IIV(iPoint,iAngle,E1_) = &
                        iphidn(iLine,iAngleDn,iIono,EE1_)
                   PlotState_IIV(iPoint,iAngle,E2_) = &
                        iphidn(iLine,iAngleDn,iIono,EE2_)
                   PlotState_IIV(iPoint,iAngle,E3_) = &
                        iphidn(iLine,iAngleDn,iIono,EE3_)
                   PlotState_IIV(iPoint,iAngle,E4_) = &
                        iphidn(iLine,iAngleDn,iIono,EE4_)
                   PlotState_IIV(iPoint,iAngle,E5_) = &
                        iphidn(iLine,iAngleDn,iIono,EE5_)
                elseif(iPoint >nIono .and. iPoint <=nPoint-nIono) then
                   iPlas=iPoint-nIono
                   PlotState_IIV(iPoint,iAngle,E1_) = &
                        phidn(iLine,iAngleDn,iPlas,EE1_)
                   PlotState_IIV(iPoint,iAngle,E2_) = &
                        phidn(iLine,iAngleDn,iPlas,EE2_)
                   PlotState_IIV(iPoint,iAngle,E3_) = &
                        phidn(iLine,iAngleDn,iPlas,EE3_)
                   PlotState_IIV(iPoint,iAngle,E4_) = &
                        phidn(iLine,iAngleDn,iPlas,EE4_)
                   PlotState_IIV(iPoint,iAngle,E5_) = &
                        phidn(iLine,iAngleDn,iPlas,EE5_)
                else
                   iIono=iPoint-nPlas
                   PlotState_IIV(iPoint,iAngle,E1_) = &
                        iphidn(iLine,iAngleDn,iIono,EE1_)
                   PlotState_IIV(iPoint,iAngle,E2_) = &
                        iphidn(iLine,iAngleDn,iIono,EE2_)
                   PlotState_IIV(iPoint,iAngle,E3_) = &
                        iphidn(iLine,iAngleDn,iIono,EE3_)
                   PlotState_IIV(iPoint,iAngle,E4_) = &
                        iphidn(iLine,iAngleDn,iIono,EE4_)
                   PlotState_IIV(iPoint,iAngle,E5_) = &
                        iphidn(iLine,iAngleDn,iIono,EE5_)
                endif
             else
                PlotState_IIV(iPoint,iAngle,:)=0.0
             endif
          enddo
       enddo
       
       !Plot grid for given line
       call save_plot_file(NamePlot, TypePositionIn='append', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn=nStep,TimeIn=time,     &
            nDimIn=nDim,CoordIn_DII=Coord_DII,                &
            VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/))
 !   end do
    
    deallocate(Coord_DII, PlotState_IIV)
  end subroutine plot_state
  
  !============================================================================
  ! 1D output plots of integrated quantities along the field (for now just pot)
  subroutine plot_along_field
    use ModSeGrid,     ONLY: FieldLineGrid_IC, DeltaPotential_C, nLine, nPoint

    use ModIoUnit,     ONLY: UnitTmp_
    use ModPlotFile,   ONLY: save_plot_file
    use ModNumConst,   ONLY: cRadToDeg
    real, allocatable   :: Coord_DII(:,:,:), PlotState_IV(:,:)
    integer, parameter :: nDim =1, nVar=1, S_=1, Pot_=1
    character(len=100),parameter :: NamePlotVar='S Pot g r'
    character(len=100) :: NamePlot
    character(len=*),parameter :: NameHeader='Pot output'
    character(len=5) :: TypePlot='ascii'
    integer :: iLine,iPoint
    logical,save :: IsFirst
    !--------------------------------------------------------------------------
    allocate(Coord_DI(nDim,nPoint), PlotState_IV(nPoint,nVar))
    
    do iLine=1,nLine
       PlotState_IIV = 0.0
       Coord_DII     = 0.0
       
       !Set Coordinates along field line and PA
       do iPoint=1,nPoint
          Coord_DI(S_,iPoint) = FieldLineGrid_IC(iLine,iPoint)/6375.0e5
          PlotState_IV(iPointPot_) = DeltaPotential_C(iLine,iPoint)
       enddo
       
       ! set name for plotfile
       write(NamePlot,"(a,i4.4,a)") 'STET_1D_iLine',iLine,'.out'
       
       !Plot grid for given line. Overwrite old results on firstcall
       if(IsFirst) then
          call save_plot_file(NamePlot, TypePositionIn='rewind', &
               TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
               NameVarIn = NamePlotVar, nStepIn= 1,TimeIn=1.0,     &
               nDimIn=nDim,CoordIn_DI=Coord_DI,                &
               VarIn_IV = PlotState_IV, ParamIn_I = (/1.6, 1.0/))
          IsFirstCall = .false.
       else
          call save_plot_file(NamePlot, TypePositionIn='append', &
               TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
               NameVarIn = NamePlotVar, nStepIn= 1,TimeIn=1.0,     &
               nDimIn=nDim,CoordIn_DI=Coord_DI,                &
               VarIn_IV = PlotState_IV, ParamIn_I = (/1.6, 1.0/))
       end if
    end do
    
    deallocate(Coord_DI, PlotState_IV)
  end subroutine plot_along_field


  !============================================================================

end Module ModSePlot
