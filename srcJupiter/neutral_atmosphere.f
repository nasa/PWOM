
                                                                 
      SUBROUTINE MODATM (ALT,XNH2,XNH,XNH2O,XNCH4,TEMP)

      REAL ALT1,ALT,XNH2,XNH,XNH2O,XNCH4,TEMP,Scaleheight_H2,
     ; Scaleheight_H,Scaleheight_CH4
!!
      
      integer, parameter ::nH2=7,nH=9,nCH4=6,nT=7,H2_=1,H_=2,CH4_=3
      integer, parameter ::H2O_=4,T_=5
      integer :: i
      real :: Alt_VC(5,10), LogDens_VC(4,10),T_C(10) 
      real :: MixingRatio = 0.0015


      
      Alt_VC(H2_,1)     =2.10566E2
      LogDens_VC(H2_,1) =1.48980E1
      Alt_VC(H2_,2)     =3.64169E2
      LogDens_VC(H2_,2) =1.31590E1
      Alt_VC(H2_,3)     =5.70484E2
      LogDens_VC(H2_,3) =1.13664E1
      Alt_VC(H2_,4)     =1.11750E3
      LogDens_VC(H2_,4) =9.20574E0
      Alt_VC(H2_,5)     =1.93269E3
      LogDens_VC(H2_,5) =7.46675E0
      Alt_VC(H2_,6)     =2.87973E3
      LogDens_VC(H2_,6) =5.99072E0
      Alt_VC(H2_,7)     =3.98366E3
      LogDens_VC(H2_,7) =4.51408E0


      Alt_VC(H_,1)     =2.17871E2
      LogDens_VC(H_,1) =9.61521E0
      Alt_VC(H_,2)     =2.78797E2
      LogDens_VC(H_,2) =1.01983E1
      Alt_VC(H_,3)     =3.94418E2
      LogDens_VC(H_,3) =1.00426E1
      Alt_VC(H_,4)     =4.69632E2
      LogDens_VC(H_,4) =9.35725E0
      Alt_VC(H_,5)     =7.34254E2
      LogDens_VC(H_,5) =8.67405E0
      Alt_VC(H_,6)     =1.17742E3
      LogDens_VC(H_,6) =7.67408E0
      Alt_VC(H_,7)     =2.03327E3
      LogDens_VC(H_,7) =6.78037E0
      Alt_VC(H_,8)     =3.17925E3
      LogDens_VC(H_,8) =6.04434E0
      Alt_VC(H_,9)     =3.97697E3
      LogDens_VC(H_,9) =5.46501E0

      
      Alt_VC(CH4_,1)     =2.11770E2
      LogDens_VC(CH4_,1) =1.16754E1
      Alt_VC(CH4_,2)     =2.71851E2
      LogDens_VC(CH4_,2) =1.04623E1
      Alt_VC(CH4_,3)     =3.40794E2
      LogDens_VC(CH4_,3) =8.66787E0
      Alt_VC(CH4_,4)     =3.96821E2
      LogDens_VC(CH4_,4) =6.60872E0
      Alt_VC(CH4_,5)     =4.62188E2
      LogDens_VC(CH4_,5) =4.39108E0
      Alt_VC(CH4_,6)     =5.11884E2
      LogDens_VC(CH4_,6) =2.64851E0

      Alt_VC(H2O_,:) = Alt_VC(H2_,:)
      LogDens_VC(H2O_,:) = LogDens_VC(H2_,:) + alog10(MixingRatio)

      T_C(1)      = 1.64E2
      Alt_VC(T_,1)= 1.89E2
      T_C(2)      = 1.93E2
      Alt_VC(T_,2)= 3.11E2
      T_C(3)      = 3.29E2
      Alt_VC(T_,3)= 4.10E2
      T_C(4)      = 6.27E2
      Alt_VC(T_,4)= 5.57E2
      T_C(5)      = 8.22E2
      Alt_VC(T_,5)= 7.35E2
      T_C(6)      = 9.11E2
      Alt_VC(T_,6)= 1.097E3
      T_C(7)      = 9.41E2
      Alt_VC(T_,7)= 2.250E3

      
      ALT1 = ALT/1.E+05

!\
! Get T
!/
      if (Alt1 >= Alt_VC(T_,nT)) then
         TEMP = T_C(nT)
      elseif (Alt1 < Alt_VC(T_,1)) then
         TEMP = T_C(1) 
      else
         do i = 1,nT-1
            if (Alt1 >= Alt_VC(T_,i) .and.  Alt1 < Alt_VC(T_,i+1)) then
               Temp = 
     &              (T_C(i+1)-T_C(i))
     &              /(Alt_VC(T_,i+1)-Alt_VC(T_,i))
     &              *(Alt1-Alt_VC(T_,i)) + T_C(i)
               exit
            else
               cycle
            endif
         enddo
      endif

!\
! Get H2 neutral Density
!/
      if (Alt1 >= Alt_VC(H2_,nH2)) then
         Scaleheight_H2 =1.380658e-26*T_C(nT)/(3.3452462e-27*22.88)
         XNH2  = (10.** LogDens_VC(H2_,nH2)) 
     &        *exp(-(ALT1-Alt_VC(H2_,nH2))/Scaleheight_H2)     
      elseif (Alt1 < Alt_VC(H2_,1)) then
         XNH2  = (10.** LogDens_VC(H2_,1))
      else
         do i = 1,nH2-1
            if (Alt1 >= Alt_VC(H2_,i) .and.  Alt1 < Alt_VC(H2_,i+1)) then
               XNH2 = 
     &              (10.0**LogDens_VC(H2_,i+1)-10.0**LogDens_VC(H2_,i))
     &              /(Alt_VC(H2_,i+1)-Alt_VC(H2_,i))
     &              *(Alt1-Alt_VC(H2_,i)) + 10.0**LogDens_VC(H2_,i)
               exit
            else
               cycle
            endif
         enddo
      endif


!\
! Get H neutral Density
!/
      if (Alt1 >= Alt_VC(H_,nH)) then
         Scaleheight_H =1.380658e-26*T_C(nT)/(3.3452462e-27*22.88)
         XNH  = (10.** LogDens_VC(H_,nH)) 
     &        *exp(-(ALT1-Alt_VC(H_,nH))/Scaleheight_H)     
      elseif (Alt1 < Alt_VC(H_,1)) then
         XNH  = (10.** LogDens_VC(H_,1))
      else
         do i = 1,nH-1
            if (Alt1 >= Alt_VC(H_,i) .and.  Alt1 < Alt_VC(H_,i+1)) then
               XNH = 
     &              (10.0** LogDens_VC(H_,i+1)-10.0** LogDens_VC(H_,i))
     &              /(Alt_VC(H_,i+1)-Alt_VC(H_,i))
     &              *(Alt1-Alt_VC(H_,i)) + 10.0**LogDens_VC(H_,i)
               exit
            else
               cycle
            endif
         enddo
      endif

!\
! Get CH4 neutral Density
!/
      if (Alt1 >= Alt_VC(CH4_,nCH4)) then
         Scaleheight_CH4 =1.380658e-26*T_C(4)/(3.3452462e-27*22.88)
         XNCH4  = (10.** LogDens_VC(CH4_,nCH4)) 
     &        *exp(-(ALT1-Alt_VC(CH4_,nCH4))/Scaleheight_CH4)     
      elseif (Alt1 < Alt_VC(CH4_,1))then
         XNCH4  = (10.** LogDens_VC(CH4_,1))
      else
         do i = 1,nCH4-1
            if (Alt1 >= Alt_VC(CH4_,i) .and.  Alt1 < Alt_VC(CH4_,i+1)) then
               XNCH4 = 
     &              (10.0** LogDens_VC(CH4_,i+1)-10.0** LogDens_VC(CH4_,i))
     &              /(Alt_VC(CH4_,i+1)-Alt_VC(CH4_,i))
     &              *(Alt1-Alt_VC(CH4_,i)) + 10.0**LogDens_VC(CH4_,i)
               exit
            else
               cycle
            endif
         enddo
      endif

! Set H2O 
      XNH2O = MixingRatio * XNH2 



      END
      !========================================================================
      subroutine test_neutral_atmosphere
      
      integer,parameter :: nAlt = 5,nNeutral=4,H2_=1,H_=2,CH4_=3,H2O_=4
      integer :: i 
      real :: Alt_C(nAlt), NeutralDens_VC(nNeutral,nAlt),T_C(nAlt)
      real :: AltMin = 500.0, dAlt = 200.0
      !------------------------------------------------------------------------
      
      do i = 1, nAlt
         Alt_C(i) = AltMin + real(i-1)*dAlt
      enddo
      
      do i = 1,nAlt
         call MODATM(
     &        ALT_C(i)* 1.0e5,NeutralDens_VC(H2_,i),NeutralDens_VC(H_,i),
     &        NeutralDens_VC(H2O_,i),NeutralDens_VC(CH4_,i),T_C(i))
      enddo
      
      NeutralDens_VC = alog10(NeutralDens_VC)

      write(*,*) 'Alt   log[H2]   log[H]   log[H2O]   log[CH4]   T'
      do i=1,nAlt
         write(*,*) Alt_C(i),NeutralDens_VC(H2_,i),NeutralDens_VC(H_,i),
     &        NeutralDens_VC(H2O_,i),NeutralDens_VC(CH4_,i),T_C(i)
      enddo
      end
