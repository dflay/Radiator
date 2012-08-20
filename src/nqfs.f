
C
C  CALCULATE QUASIELASTIC, TWO-NUCLEON, RESONANCE, AND X-SCALING
C  ELECTRON SCATTERING CROSS SECTION FOR ARBITRARY NUCLEI
C
C  WRITTEN BY J.W. LIGHTBODY JR. AND J.S. O'CONNELL
C  NATIONAL BUREAU OF STANDARDS, GAITHERSBURG, MD 20899
C  OCTOBER 1987

c  09/16/02 K. SlIFer, Temple U.
c  ModIFied subroutine "RADIATE" to include external radiation using
c  formalism from S. Stein et al. Phys. Rev. D 12 7
c  In the resonance region, the new code gives identical results
c  to original code IF external radiation is turned off.
c  At larger nu, the new code predicts a much larger tail than the original.
c
c  Included option 'Qdep' to introduce a Q^2 dependent correction to the
c  Q.E. peak, based on a comparison with Saclay and SLAC carbon data.
c  The correction is only to the peak height.  A more sophistacated correction
c  would adjust the gaussian width as well.
c  
c  The most recent F.F. are provided as an option, but they are commented 
c  out at present. See FUNCTIONs GEP,GEN for example  
c 
c  If the radiated cross section exhibits glitches adjust the parameters
c  "PREC" or "DEL" in subroutine RADIATE
c
c  08/08/05 P. Solvignon, Temple University
c  Made several changes to improve Nitrogen xs model:
c       - change the dipole parameter AD1 in the delta region from 
c         700 to 685 MeV.
c       - multiplie delta cross section by 0.1+0.75*Ep
c       - multiplie dip cross section by 0.65
c       - multiplie SIG_1500 by 4.0/E
c       - multiplie SIG_1700 by Q2/100
c       - divide SIG_DIS by 2.7/(4.0*E*sin(THR/2.0))
c  Used radiation length proper to E01-012 experiment.  
c
c  05/29/11 D. Flay, Temple University
c  Made the following changes to improve the Nitrogen xs model:
c       - multiplied the delta cross section by 0.1-2.65*Ep
c       - multiplied the dip cross section by 20*(1+W/E)
c           - provides good agreement at 4- and 5-pass

      PROGRAM QFS
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL Qdep,INTERACTIVE,RATES,RAD,extrad
      CHARACTER TAG*4,LINE*80,nucleus*3,units*2
      COMMON/PAR/E,TH,W,Z,A,EPS,EPSD,PF,SPENCE
      COMMON/ADDONS/SCALE,Tb,Ta,extrad
      COMMON/QDEPENDENCE/XPAR0,XPAR1,Qdep
      REAL*8 MYPAR1,MYPAR2,MYPAR3,MYPAR4,MYPAR5,MYPAR6,
     &       MYPAR7,MYPAR8,MYPAR9,MYPAR10,MYPAR11,MYPAR12,
     &       ANGPARX,ESPAR,QEANGPAR,QEPAR,DELANGPAR,DELPAR, 
     &       DISANGPAR,DISPAR,R1PAR,R2PAR
      REAL*8 FPAR1,FPAR2,FPAR3,FPAR4,FPAR5,FPAR6
      REAL*8 SF1,SF2,SF3,SF4,SF5,SF6 
      COMMON/MYPAR/MYPAR1,MYPAR2,MYPAR3,MYPAR4,MYPAR5,MYPAR6,
     &             MYPAR7,MYPAR8,MYPAR9,MYPAR10,MYPAR11,MYPAR12,
     &             ANGPARX,ESPAR,QEANGPAR,QEPAR,DELANGPAR,DELPAR, 
     &             DISANGPAR,DISPAR,R1PAR,R2PAR

C      REAL*8 MYX,MYX2,MYG,MYS,MYE,MYA
      REAL*8 current, density,zlength,domega,deltap,P0(30),P0N(30)
      REAL*8 sigmatot(1000), xnu(1000),sigave(30)
      DATA sigmatot/1000*0.0/,xnu/1000*0.0/,sigave/30*0.0/
C     Variables for Patricia's kinematics
      PARAMETER(NPTS=200)
      REAL MyW(NPTS),MyEp(NPTS)
      CHARACTER PATH*26

      PM=939.D0
      DM=1232.D0
      ALPH=1./137.03604D0
      HBARC=197.32858D0
      PI=ACOS(-1.)

      INTERACTIVE=.FALSE. ! Prompt user for input? Otherwise READ from input.dat
      RATES      =.FALSE. ! Calculat the rates? (THIS FEATURE NOT QUITE DEBUGGED YET)
      RAD        =.FALSE. ! turn on R.C. (internal only)
      extrad     =.FALSE. ! turns on external R.C. 
      units      ='pb'    ! units/MeV-sr
      Qdep       =.FALSE. ! Correct Q.E. peak to agree with world data.

c      OPEN(20,FILE='input_N2.dat',STATUS='old')
       OPEN(20,FILE='./input/input.dat',STATUS='old')
C      Read in xs parameters
       OPEN(27,FILE='./input/input-params.dat',STATUS='old') 
       READ(27,'(A)',ERR=999) LINE 
       READ(LINE,*) MYPAR1,MYPAR2,MYPAR3,MYPAR4,MYPAR5,MYPAR6  
       OPEN(98,FILE='./input/input-params-f.dat',STATUS='old') 
       READ(98,'(A)',ERR=999) LINE 
       READ(LINE,*) FPAR1,FPAR2,FPAR3,FPAR4,FPAR5,FPAR6  
C      Read in scale factor (systematic error study)  
       OPEN(97,FILE='./input/scale-factor.dat',STATUS='old') 
       READ(97,'(A)',ERR=999) LINE 
       READ(LINE,*) SF1,SF2,SF3,SF4,SF5,SF6   

C       p1--6 are from the 10/27/11 talk 
C       MYPAR1=1.25709
C       MYPAR2=1.39379
C       MYPAR3=0.04198
C       MYPAR4=0.99429
C       MYPAR5=0.98722
C       MYPAR6=0.75728
C       MYPAR7=0.2
C       MYPAR8=-1.12
C       MYPAR9=1.10
C       MYPAR10=0.27
C       MYPAR11=0.80

C       WRITE(*,*) '[nqfs.f]: QFS Parameters: '
C       WRITE(*,*) 'Free: '
C       WRITE(*,'(A,1F7.5)') 'p0  =   ',MYPAR1
C       WRITE(*,'(A,1F7.5)') 'p1  =   ',MYPAR2
C       WRITE(*,'(A,1F7.5)') 'p2  =   ',MYPAR3
C       WRITE(*,'(A,1F7.5)') 'p3  =   ',MYPAR4
C       WRITE(*,'(A,1F7.5)') 'p4  =   ',MYPAR5
C       WRITE(*,'(A,1F7.5)') 'p5  =   ',MYPAR6
C      Scale each parameter by a scale factor  
C       MYPAR1 = MYPAR1*SF1 
C       MYPAR2 = MYPAR2*SF2
C       MYPAR3 = MYPAR3*SF3
C       MYPAR4 = MYPAR4*SF4
C       MYPAR5 = MYPAR5*SF5
C       MYPAR6 = MYPAR6*SF6
C       WRITE(*,*) 'Scale Factors: '
C       WRITE(*,'(A,1F7.5)') 'SF1  =   ',SF1
C       WRITE(*,'(A,1F7.5)') 'SF2  =   ',SF2
C       WRITE(*,'(A,1F7.5)') 'SF3  =   ',SF3
C       WRITE(*,'(A,1F7.5)') 'SF4  =   ',SF4
C       WRITE(*,'(A,1F7.5)') 'SF5  =   ',SF5
C       WRITE(*,'(A,1F7.5)') 'SF6  =   ',SF6
           
C       WRITE(*,'(A,1F7.5)') 'p7  =   ',MYPAR7
C       WRITE(*,'(A,1F7.5)') 'p8  =   ',MYPAR8
C       WRITE(*,'(A,1F7.5)') 'p9  =   ',MYPAR9
C       WRITE(*,'(A,1F7.5)') 'p10 =   ',MYPAR10
C       WRITE(*,'(A,1F7.5)') 'p11 =   ',MYPAR11
C       WRITE(*,'(A,1F7.5)') 'p11 =   ',MYPAR12

990   IF(INTERACTIVE) THEN
        CALL PROMPT(Z,A,E,TH,PF,EPS,EPSD,Tb,Ta)
        WRITE(6,'(A)') 'ENTER WMIN,WMAX,DELTAW: '
        READ(*,*) WMIN,WMAX,DELTAW
c         WIN = 0.0
c         WMAX = 3500.0
c         DELTAW = 10.0

         TAG     = '4_NI'
         nucleus = 'Nit'
         GOTO 991
      ENDIF

      READ(20,'(A)',ERR=999) LINE
      IF (LINE(1:1).EQ.'#') GOTO 990
      READ(LINE,*) TAG,nucleus,TH,PF,EPS,EPSD

991   CONTINUE

      WRITE(*,*) TAG,' ',nucleus 

C      Read in data for Patricia's kinematics [ (W,Ep) for Nu calculations]  
      IF(TAG.EQ.'4_DU') THEN
        MyMax = 109
      ELSEIF(TAG.EQ.'6_DU') THEN
        MyMax = 95
      ENDIF

      IF((TAG.EQ.'4_DU').OR.(TAG.EQ.'6_DU')) THEN
        PATH='./infiles/kin/'//TAG//'_kin.dat'
        OPEN(28,FILE=PATH,STATUS='old')
        DO J=1,MyMax
           READ(28,*) MyW(J),MyEp(J)
        ENDDO
      ENDIF

      CALL TARGETCELL(TAG,INTERACTIVE,E,WMIN,WMAX,DELTAW,Tb,Ta,
     &                 density,zlength)
      CALL GETZA(nucleus,INTERACTIVE,Z,A)
      CALL OPENFILES(TAG,INTERACTIVE,RATES,nucleus)
      CALL GETSCALE(units,SCALE,UNSCALE)

      CALL GETQDEP(A,Qdep,XPAR0,XPAR1)
      ! WRITE(8,'(2A,x,2F4.0,3A)')'#',nucleus,Z,A,':  ',units,'/MeV-sr'
      ! WRITE(8,'(A,2F7.1)')      '#E0 (MeV),TH : ',E,TH
      ! WRITE(8,'(A,3F7.1)')      '#PF,EPS,EPSD : ',PF,EPS   ,EPSD
      ! WRITE(8,'(A,L7,2F7.4)')   '#extrad,Tb,Ta: ',extrad,Tb,Ta
      ! WRITE(8,'(A,L7,2F7.4)')   '#Q Dependence: ',Qdep,XPAR0,XPAR
C------

C      DELTAW = 1.0

      Pmax=E-WMIN
      Pmin=E-WMAX
      
      IF (RATES) THEN
        CALL GETRATES(deltaP,current,domega)
        CALL INITIALIZE(Pmax,deltaP,xnu,sigmatot,sigave,P0,P0N) 
      ENDIF

      W=WMIN

      EQE    = 4.2386 
      EQE2   = 4.018 
 
      THR    = TH*PI/180.
      TH0R   = 15.5*PI/180.
      TH1R   = 25.0*PI/180.
      TH2R   = 32.0*PI/180.
      TH3R   = 45.0*PI/180.

C------------ 'Derived' Parameters for Each Sub-process ------

      QEPAR  = -5.6690 + 2.2371*(E*1E-3) - 2.1351E-1*(E*1E-3)**2
      QEPAR  = 1.0 - QEPAR 

      DELPAR = QEPAR 

      R1PAR  = 1.0771 - 1.2988E-1*TH + 4.7559E-3*TH**2  
     &          - 5.0074E-5*TH**3
      R1PAR  = 1 - R1PAR 

      R2PAR  = R1PAR

C      DISPAR = -3.3643E-1 + 1.7805E-2*TH - 8.1388E-5*TH**2 
      DISPAR = -2.7726E-1 + 1.4641E-2*TH - 4.0299E-5*TH**2 
      DISPAR = 1.0 - DISPAR 

      ! NOTE: DIP REGION DOES NOT HAVE AN F PARAMETER. WOULD BE INDEX 2!  
      QEPAR  = QEPAR *(1.0 + FPAR1)  
      DELPAR = DELPAR*(1.0 + FPAR3)
      R1PAR  = R1PAR *(1.0 + FPAR4) 
      R2PAR  = R2PAR *(1.0 + FPAR5)
      DISPAR = DISPAR*(1.0 + FPAR6) 
  
C-------------- Test parameterization for Karl's data 

C      QEPAR  = 8.3531E-1 - 4.6343E-1*(E*1E-3) + 6.6650E-2*(E*1E-3)**2
C      QEPAR  = 1.0 - QEPAR 
C      DELPAR = QEPAR 
C
C      R1PAR  = 1.0771 - 1.2988E-1*TH + 4.7559E-3*TH**2  
C     &          - 5.0074E-5*TH**3
C      R1PAR  = 1 - R1PAR 
C
C      R2PAR  = R1PAR
C
C      DISPAR = -3.3643E-1 + 1.7805E-2*TH - 8.1388E-5*TH**2 
C      DISPAR = 1.0 - DISPAR 

C      Reset parmeters (for comparison of original to optimized model) 
C      QEPAR  = 1.0
C      DELPAR = 1.0
C      R1PAR  = 1.0
C      R2PAR  = 1.0
C      DISPAR = 1.0
C
C      MYPAR1 = 1.0
C      MYPAR2 = 1.0
C      MYPAR3 = 1.0
C      MYPAR4 = 1.0
C      MYPAR5 = 1.0
C      MYPAR6 = 1.0

    
      ii=0

      WRITE(*,'(A,F10.1)') 'E0:    ',E
      WRITE(*,'(A,F10.1)') 'TH:    ',TH
C      WRITE(*,*) 'F PARAMETER PERCENT OFFSETS: '
C      WRITE(*,'(A,1F7.5)') 'QE:  ',FPAR1 
C      WRITE(*,'(A,1F7.5)') 'DEL: ',FPAR3
C      WRITE(*,'(A,1F7.5)') 'R1:  ',FPAR4
C      WRITE(*,'(A,1F7.5)') 'R2:  ',FPAR5
C      WRITE(*,'(A,1F7.5)') 'DIS: ',FPAR6
C      WRITE(*,*) 'DERIVED PARAMETERS: '
C      WRITE(*,'(A,1F7.5)') 'pTHX:  ',ANGPARX
C      WRITE(*,'(A,1F7.5)') 'pEs:   ',ESPAR
C      WRITE(*,'(A,1F7.5)') 'fQE:  ',QEPAR
C      WRITE(*,'(A,1F7.5)') 'fDEL: ',DELPAR
C      WRITE(*,'(A,1F7.5)') 'fR1:  ',R1PAR
C      WRITE(*,'(A,1F7.5)') 'fR2:  ',R2PAR
C      WRITE(*,'(A,1F7.5)') 'fDIS: ',DISPAR 

10    CONTINUE 

      ii=ii+1
      IF((TAG.EQ.'4_DU').OR.(TAG.EQ.'6_DU')) THEN
        W = (1./(2.*PM))*( MyW(ii)**2 + 4.0*E*MyEp(ii)*
     &           SIN(THR/2.)**2 - PM**2 )
      ELSE
        W = W + DELTAW  ! E-E'
      ENDIF

      IF((TAG.EQ.'4_DU').OR.(TAG.EQ.'6_DU')) THEN
         IF(ii.GE.MyMax) GO TO 30
      ENDIF

      IF(W.GE.WMAX)GO TO 30
      IF(W.LE.1.)GO TO 30

      SIGQFZA=0
      SIGDA  =0
      SIG2NA =0 
      SIGR1A =0
      SIGR2A =0
      SIGXA  =0 

      QMS=4.*E*(E-W)*SIN(THR/2.)**2
      SIGQFZA=SIGQFS(E,TH,W,Z,A,EPS,PF)*SCALE ! Q.E. Peak
      SIGDA=SIGDEL(E,TH,W,A,EPSD,PF)*SCALE    ! Delta
      SIG2NA=SIG2N(E,TH,W,Z,A,PF)*SCALE       ! Dip region
      SIGR1A=SIGR1(E,TH,W,A,PF)*SCALE         ! Resonance 1500MeV
      SIGR2A=SIGR2(E,TH,W,A,PF)*SCALE         ! Resonance 1700MeV
      SIGXA=SIGX(E,TH,W,A)*SCALE              ! DIS
      SIG=SIGQFZA+SIGDA+SIGR1A+SIGXA+SIGR2A+SIG2NA

c      WRITE (8,'(F7.1,8E10.3)')xW,SIGQFZA,SIG2NA,SIGDA,SIGR1A,SIGR2A,
c     &                         SIGXA,SIG,SIGQFZA+SIGDA+SIGR1A

C     Photon scattering
c----------------------
      xW = PM**2+2.*PM*W-QMS                  ! missing mass squared
      xW = sqrt(abs(xW)) 
      QVS=QMS+W**2
      EPD=E-(DM-PM)*(DM+PM)/2./PM
      EPD=EPD/(1.+2.*E*SIN(THR/2.)**2/PM)
      EPD=EPD-EPSD
      EPQ=4.*E**2*SIN(THR/2.)**2/2./PM
      EPQ=EPQ/(1.+2.*E*SIN(THR/2.)**2/PM)+EPS
      EPQ=E-EPQ
      EKAPPA=W-QMS/2./PM
      IF(EKAPPA.GT.-PM/2.)THEN
        CMTOT=SQRT(PM**2+2.*PM*EKAPPA)
      ELSE
        CMTOT=PM
      ENDIF
      FLUX=(ALPH/2./PI**2)*((E-W)/E)*((2.*PM*W-QMS)/2./PM/QMS)
      POLARI=1./(1.+2.*QVS*TAN(THR/2.)**2/QMS)
      FLUX=FLUX/(1.-POLARI)
      IF(EKAPPA.LT.0.)PHOTSIG=0.
      PHOTSIG=SIG/FLUX        ! PHOTSIG never used

C------

      IF(.NOT.RAD) THEN
        WRITE (8,'(3F7.1,2E10.3)') E,W,xW,SIG,0.0 
        GO TO 10
      ENDIF

      CALL RADIATE(E,TH,W,SIG,SIGRAD)
c-- OUTPUT
c      WRITE (8,'(2F7.1,2E10.3)') W,xW,SIGRAD,SIG   ! nu,W,exp. & born XS in "units"
      ! E,nu,W, born XS and 0.0 in "units"  for radcor.f   
      WRITE (8,'(3F7.1,A,1E10.3,A,1E10.3)') E,W,xW,'  ',SIGRAD,'  ',0.0 

      IF(TAG.EQ.'3REF') THEN
         WRITE (33,'(F7.1,E10.3,4F7.1)') xW,SIGRAD,0.0,W,0.0,0.0   ! W,exp., 0.0, nu,0.0,0.0 for smooth.f
      ELSEIF(TAG.EQ.'4REF') THEN 
         WRITE (34,'(F7.1,E10.3,4F7.1)') xW,SIGRAD,0.0,W,0.0,0.0   ! W,exp., 0.0, nu,0.0,0.0 for smooth.f
      ELSEIF(TAG.EQ.'5REF') THEN 
         WRITE (35,'(F7.1,E10.3,4F7.1)') xW,SIGRAD,0.0,W,0.0,0.0   ! W,exp., 0.0, nu,0.0,0.0 for smooth.f
      ELSEIF(TAG.EQ.'6REF') THEN 
         WRITE (36,'(F7.1,E10.3,4F7.1)') xW,SIGRAD,0.0,W,0.0,0.0   ! W,exp., 0.0, nu,0.0,0.0 for smooth.f
      ENDIF
c-----------
      IF(RATES) THEN
        sigmatot(ii)=SIGRAD*UNSCALE                  ! convert back to cm^2 for rates.
        xnu(ii)=W
      ENDIF

      GO TO 10
  30  CONTINUE
      NBIN=ii
      close(8) 

C-----COMPUTE THE RATES 
      IF (RATES) THEN
      IF (Z.EQ.2) THEN ! Helium-3
        XLum=(current/1.6E-19)*density*zlength ! Luminosity

        do i=1,NBIN
           do j=1,30
             Pup= P0(j)*(1.0+deltap)
             Plow=P0(j)*(1.0-deltap)
             escat=E-xnu(i)
             IF (escat.LE.Pup.AND.ESCAT.GE.Plow) THEN
               sigave(j)=sigave(j)+sigmatot(i)
               P0N(j)=P0N(j)+1
             ENDIF
           enddo
        enddo
        do j=1,30 
          WRITE(6,*) "DEBUG3",Pmin,P0N(j),P0(j)
          IF ( (P0(j).GT.Pmin).and.(P0N(j).gt.0) ) THEN
            sigave(j)=sigave(j)/P0N(j)
            deltaE=P0(j)*(2*deltaP)
            rate=Xlum*sigave(j)*domega*deltaE
            WRITE(10,*) P0(j), rate
          ENDIF
        enddo
      ENDIF
      ENDIF
C------

      IF((TAG.EQ.'4_DU').OR.(TAG.EQ.'6_DU')) THEN
        CALL ZERO(MyW,MyMax)
        CALL ZERO(MyEp,MyMax)
      ENDIF

      GOTO 990

999   CONTINUE
      STOP
      END
C-------------------------------------------------------------------
      SUBROUTINE ZERO(X,N)

      DIMENSION X(N)

      DO I = 1, N
      X(I) = 0.
      ENDDO

      RETURN
      END
C-------------------------------------------------------------------
*FD
      REAL*8 FUNCTION FD(QMS,A)
      IMPLICIT REAL*8 (A-H,O-Z)
      FD=1./(1.+QMS/A**2)**2
      RETURN
      END
C-------------------------------------------------------------------
*FM
      REAL*8 FUNCTION FM(QMS,A)
      IMPLICIT REAL*8 (A-H,O-Z)
      FM=1./(1.+QMS/A**2)
      RETURN
      END
C-------------------------------------------------------------------
*FPHENOM
      REAL*8 FUNCTION FPHENOM(QMS)
      IMPLICIT REAL*8 (A-H,O-Z)
      A1=.55
      A2=20./1.E6
      B1=.45
      B2=.45/1.E6
      C1=0.03
      C2=0.2/1.E12
      FPHENOM=A1*EXP(-A2*QMS)+B1*EXP(-B2*QMS)
      FPHENOM=FPHENOM+C1*EXP(-C2*(QMS-4.5E6)**2)
      FPHENOM=SQRT(FPHENOM)
      RETURN
      END
C-------------------------------------------------------------------
*FYUKAWA
      REAL*8 FUNCTION FYUKAWA(QMS,A)
      IMPLICIT REAL*8 (A-H,O-Z)
      IF(QMS.LT.1.E-5.OR.A.LT.1.E-5)THEN
      FYUKAWA=0.
      ELSE
      ARG=SQRT(QMS/2.)/A
      FYUKAWA=ATAN(ARG)/ARG
      ENDIF
      RETURN
      END
C-------------------------------------------------------------------
*SIGMOT
      REAL*8 FUNCTION SIGMOT(E,THR)
      IMPLICIT REAL*8 (A-H,O-Z)
      ALPH=1./137.03604
      HBARC=197.3286
      SIGMOT=(ALPH*HBARC*COS(THR/2.)/2./E/SIN(THR/2.)**2)**2
C  FM**2/SR
      RETURN
      END
C-------------------------------------------------------------------
*RECOIL
      REAL*8 FUNCTION RECOIL(E,THR,TM)
      IMPLICIT REAL*8 (A-H,O-Z)
      RECOIL=1./(1.+2.*E*SIN(THR/2.)**2/TM)
      RETURN
      END
C-------------------------------------------------------------------
*SIGX = DIS cross section
      REAL*8 FUNCTION SIGX(E,TH,W,A)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 MYPAR1,MYPAR2,MYPAR3,MYPAR4,MYPAR5,MYPAR6,
     &       MYPAR7,MYPAR8,MYPAR9,MYPAR10,MYPAR11,MYPAR12,
     &       ANGPARX,ESPAR,QEANGPAR,QEPAR,DELANGPAR,DELPAR, 
     &       DISANGPAR,DISPAR,R1PAR,R2PAR
      COMMON/MYPAR/MYPAR1,MYPAR2,MYPAR3,MYPAR4,MYPAR5,MYPAR6,
     &             MYPAR7,MYPAR8,MYPAR9,MYPAR10,MYPAR11,MYPAR12,
     &             ANGPARX,ESPAR,QEANGPAR,QEPAR,DELANGPAR,DELPAR, 
     &             DISANGPAR,DISPAR,R1PAR,R2PAR

      ALPH=1./137.03604
      PI=ACOS(-1.)
C     SIG0=111.*1.E-30
      SIG0=100.D-4
C     SIG1=60.*1.E-27
      SIG1=54.*1.D-1
      PIMASS=140.
      PM=939.
C     GAM0=550.
      GAM0=650.
C     R=0.10
      AQ=250.    ! NEVER USED  
      THR=TH*PI/180.
      IF(W.LT.1.E-5)GO TO 4
      QMS=4.*E*(E-W)*SIN(THR/2.)**2
      ARG0=W-QMS/2./PM-PIMASS-PIMASS**2/2./PM
      ARG1=ARG0/GAM0
      ARG=ARG1**2/2.
      IF(ARG1.GT.8.)THEN
      SHAPE=1.+SIG1/SIG0/ARG0
      ELSEIF(ARG1.LT.1.E-5)THEN
      SHAPE=0.
      ELSEIF(ARG1.LT.0.1)THEN
      SHAPE=SIG1*ARG0/2./GAM0**2/SIG0
      ELSE
      SHAPE=(1.-EXP(-ARG))*(1.+SIG1/SIG0/ARG0)
      ENDIF
      EKAPPA=W-QMS/2./PM
      SIGGAM=SIG0*SHAPE
      QS=QMS+W**2
      EPS=1./(1.+2.*QS*TAN(THR/2.)**2/QMS)
      FLUX=ALPH*EKAPPA*(E-W)/2./PI**2/QMS/E/(1.-EPS)
      IF(FLUX.LT.1.E-20)FLUX=0.
      SIGEE=FLUX*SIGGAM*FPHENOM(QMS)**2
C     SIGEE=FLUX*SIGGAM
      R=0.56*1.E6/(QMS+PM**2)
      FACTOR1=1.+EPS*R
      SIGEE=SIGEE*FACTOR1
 4    IF (A.eq.14.0) THEN
C        SIGX=A*SIGEE*(2.7*(E-W)*sin(THR/2.0)/(QMS*1E-3))
C        SIGX=A*SIGEE*(1.0/sin(THR/2.0))*(E/(E-W))*(1.0/16.5)*(1.25)
        SIGX=A*SIGEE*(1.0/sin(THR/2.0))*(E/(E-W))*1.0/16.5
      ELSE IF (A.eq.3.0) THEN
C        SIGX=A*SIGEE*(4.0*(E-W)*sin(THR/2.0)/(QMS*1E-3))*0.8*(E/4730.0)
C        SIGX=A*SIGEE*(4.0/(SIN(THR/2.0)))*0.86*(1E+2/(E-W))*
C     &       (5009/E)*(SIN(32.0*PI/180.0)/SIN(THR))
C        SIGX=A*SIGEE*(4.0/(QMS*1E-6))*(1E-3*E)*(1E+2/(E-W))
C        IF (TH .EQ. 15.5) THEN
C          SIGX=A*SIGEE*0.85*(SIN(15.5*(PI/180.0))/SIN(THR))
C        ELSE IF (TH . EQ. 25 ) THEN 
C          SIGX=A*SIGEE*0.85*(SIN(15.5*(PI/180.0))/SIN(THR))
C     &         *EXP(0.1*W*1E-3)*(W*1E-3)*0.65
C        ELSE IF (TH . EQ. 32 ) THEN
C          SIGX=A*SIGEE*0.85*(SIN(15.5*(PI/180.0))/SIN(THR))
C     &         *EXP(0.1*W*1E-3)*(W*1E-3)*0.35
C        ELSE IF (TH . EQ. 45 ) THEN
C          SIGX=A*SIGEE*0.85*(SIN(15.5*(PI/180.0))/SIN(THR))
C     &         *EXP(0.1*W*1E-3)*(1E+1/(E-W))
C        ENDIF
C        SIGX=A*SIGEE*0.85*(SIN(15.5*(PI/180.0))/SIN(THR)) 
        SIGX=A*SIGEE*MYPAR6*DISPAR
      ELSE 
        SIGX=A*SIGEE
      ENDIF
      RETURN
      END
C-------------------------------------------------------------------
*SIGR1
      REAL*8 FUNCTION SIGR1(E,TH,W,A,PF)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 MYPAR1,MYPAR2,MYPAR3,MYPAR4,MYPAR5,MYPAR6,
     &       MYPAR7,MYPAR8,MYPAR9,MYPAR10,MYPAR11,MYPAR12,
     &       ANGPARX,ESPAR,QEANGPAR,QEPAR,DELANGPAR,DELPAR, 
     &       DISANGPAR,DISPAR,R1PAR,R2PAR
      COMMON/MYPAR/MYPAR1,MYPAR2,MYPAR3,MYPAR4,MYPAR5,MYPAR6,
     &             MYPAR7,MYPAR8,MYPAR9,MYPAR10,MYPAR11,MYPAR12,
     &             ANGPARX,ESPAR,QEANGPAR,QEPAR,DELANGPAR,DELPAR, 
     &             DISANGPAR,DISPAR,R1PAR,R2PAR

      PI=ACOS(-1.)
      PM=939.
      PIMASS=140.
      THR=TH*PI/180.
      PFR=230.
      RM=1500.
      EPSR=0.
      AR0=1000.
      AR1=1000.
      GAMQFR=120.
      GAMSPRD=140.
      GAMR=110.
      GAMPI=5.
      QFRP=1.20D-7
      QMSQFR=4.*730.*(730.-115.)*SIN(37.1*PI/180./2.)**2
      QVSQFR=QMSQFR+115.**2
      QMSRR=4.*10000.*(10000.-1240.)*SIN(6.*PI/180./2.)**2
      QVSRR=QMSRR+1240.**2
      SIGREF=FD(QMSRR,AR0)**2*QVSRR
      SIGREF=SIGREF*(QMSRR/2./QVSRR+TAN(6.*PI/180./2.)**2)
      SIGREF=SIGREF*SIGMOT(10000.D0,6.*PI/180.)
      NA=INT(A)
      IF(NA.EQ.1)THEN
      QFR=QFRP
      GSPRDA=0.
      AR=AR0
      ELSEIF(NA.LT.4)THEN
      QFR=QFRP
      GSPRDA=(A-1.)*GAMSPRD/3.
      AR=AR0+(A-1.)*(AR1-AR0)/3.
      ELSE
      AR=AR1
      GSPRDA=GAMSPRD
      QFR=QFRP
      ENDIF
      QMS=4.*E*(E-W)*SIN(THR/2.)**2
      QVS=QMS+W**2
      IF(NA.GT.1)THEN
      GAMQ=GAMQFR*PF*SQRT(QVS)/PFR/SQRT(QVSQFR)
      ELSE
      GAMQ=0.
      ENDIF
      CMTOT2=PM**2+2.*PM*W-QMS
      WTHRESH=4.*E**2*SIN(THR/2.)**2+PIMASS**2+2.*PIMASS*PM
      WTHRESH=WTHRESH/2./PM
      THRESHD=1.+PF/PM+PF**2/2./PM**2+2.*E*SIN(THR/2.)**2/PM
      WTHRESH=WTHRESH/THRESHD
c
      WTHRESH=0.0
c
      IF(W.GT.WTHRESH)THEN
      THRESH=1.-EXP(-(W-WTHRESH)/GAMPI)
      ELSE
      THRESH=0.
      ENDIF
      EPR=E-(RM-PM)*(RM+PM)/2./PM
      EPR=EPR/(1.+2.*E*SIN(THR/2.)**2/PM)
      EPR=EPR-EPSR
      WR=E-EPR
      GAM=SQRT(GAMR**2+GAMQ**2+GSPRDA**2)
      SIGR=QFR*(GAMR/GAM)/SIGREF
      SIGR=SIGR*CMTOT2*GAM**2
      SIGR=SIGR/((CMTOT2-(RM+EPSR)**2)**2+CMTOT2*GAM**2)
      SIGR=SIGR*QVS*FD(QMS,AR)**2
      SIGR=SIGR*(QMS/2./QVS+TAN(THR/2.)**2)
      SIGR=SIGR*SIGMOT(E,THR)
      SIGR1=A*THRESH*SIGR
      IF (A.eq.14.0) THEN
        SIGR1=SIGR1
C        SIGR1=SIGR1
      ELSE IF (A.eq.3.0) THEN
C        SIGR1=SIGR1*(4.0/(E*1E-3))
C        SIGR1=SIGR1*(1.0/(SIN(THR)))
C        SIGR1=SIGR1*(3.0/(QMS*1E-6))
C        SIGR1=SIGR1*(SIN((32.0)*(PI/180.0))/SIN(THR))*0.8
        SIGR1=SIGR1*MYPAR4*R1PAR
      ELSE
        SIGR1=SIGR1
      ENDIF
      RETURN
      END
C-------------------------------------------------------------------
*SIGR2
      REAL*8 FUNCTION SIGR2(E,TH,W,A,PF)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 MYPAR1,MYPAR2,MYPAR3,MYPAR4,MYPAR5,MYPAR6,
     &       MYPAR7,MYPAR8,MYPAR9,MYPAR10,MYPAR11,MYPAR12,
     &       ANGPARX,ESPAR,QEANGPAR,QEPAR,DELANGPAR,DELPAR, 
     &       DISANGPAR,DISPAR,R1PAR,R2PAR
      COMMON/MYPAR/MYPAR1,MYPAR2,MYPAR3,MYPAR4,MYPAR5,MYPAR6,
     &             MYPAR7,MYPAR8,MYPAR9,MYPAR10,MYPAR11,MYPAR12,
     &             ANGPARX,ESPAR,QEANGPAR,QEPAR,DELANGPAR,DELPAR,
     &             DISANGPAR,DISPAR,R1PAR,R2PAR

      PI=ACOS(-1.)
      PM=939.
      PIMASS=140.
      THR=TH*PI/180.
      PFR=230.
      RM=1700.
      EPSR=0.
      AR0=1200.
      AR1=1200.
      GAMQFR=120.
      GAMSPRD=140.
      GAMR=110.
      GAMPI=5.
      QFRP=0.68D-7
      QMSQFR=4.*730.*(730.-115.)*SIN(37.1*PI/180./2.)**2
      QVSQFR=QMSQFR+115.**2
      QMSRR=4.*10000.*(10000.-1520.)*SIN(6.*PI/180./2.)**2
      QVSRR=QMSRR+1520.**2
      SIGREF=FD(QMSRR,AR0)**2*QVSRR
      SIGREF=SIGREF*(QMSRR/2./QVSRR+TAN(6.*PI/180./2.)**2)
      SIGREF=SIGREF*SIGMOT(10000.D0,6.*PI/180.)
      NA=INT(A)
      IF(NA.EQ.1)THEN
      QFR=QFRP
      GSPRDA=0.
      AR=AR0
      ELSEIF(NA.LT.4)THEN
      QFR=QFRP
      GSPRDA=(A-1.)*GAMSPRD/3.
      AR=AR0+(A-1.)*(AR1-AR0)/3.
      ELSE
      AR=AR1
      GSPRDA=GAMSPRD
      QFR=QFRP
      ENDIF
      QMS=4.*E*(E-W)*SIN(THR/2.)**2
      QVS=QMS+W**2
      IF(NA.GT.1)THEN
      GAMQ=GAMQFR*PF*SQRT(QVS)/PFR/SQRT(QVSQFR)
      ELSE
      GAMQ=0.
      ENDIF
      CMTOT2=PM**2+2.*PM*W-QMS
      WTHRESH=4.*E**2*SIN(THR/2.)**2+PIMASS**2+2.*PIMASS*PM
      WTHRESH=WTHRESH/2./PM
      THRESHD=1.+PF/PM+PF**2/2./PM**2+2.*E*SIN(THR/2.)**2/PM
      WTHRESH=WTHRESH/THRESHD
      IF(W.GT.WTHRESH)THEN
      THRESH=1.-EXP(-(W-WTHRESH)/GAMPI)
      ELSE
      THRESH=0.
      ENDIF
      EPR=E-(RM-PM)*(RM+PM)/2./PM
      EPR=EPR/(1.+2.*E*SIN(THR/2.)**2/PM)
      EPR=EPR-EPSR
      WR=E-EPR
      GAM=SQRT(GAMR**2+GAMQ**2+GSPRDA**2)
      SIGR=QFR*(GAMR/GAM)/SIGREF
      SIGR=SIGR*CMTOT2*GAM**2
      SIGR=SIGR/((CMTOT2-(RM+EPSR)**2)**2+CMTOT2*GAM**2)
      SIGR=SIGR*QVS*FD(QMS,AR)**2
      SIGR=SIGR*(QMS/2./QVS+TAN(THR/2.)**2)
      SIGR=SIGR*SIGMOT(E,THR)
      SIGR2=A*THRESH*SIGR
      IF (A.eq.14.0) THEN
        SIGR2=SIGR2*(0.01*QMS*1E-6)*(E/W)*5
        SIGR2=SIGR2
      ELSE IF (A.eq.3.0) THEN
C        SIGR2=SIGR2*(0.01*QMS*1E-6)
C        SIGR2=SIGR2*(0.5*QMS*1E-6)
C        SIGR2=SIGR2*1.1*(SIN(THR)/SIN(15.5*(PI/180.0))) 
        SIGR2=SIGR2*MYPAR5*R2PAR
      ELSE
        SIGR2=SIGR2
      ENDIF
      RETURN
      END
C-------------------------------------------------------------------
*SIG2N
      REAL*8 FUNCTION SIG2N(E,TH,W,Z,A,PF)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 MYPAR1,MYPAR2,MYPAR3,MYPAR4,MYPAR5,MYPAR6,
     &       MYPAR7,MYPAR8,MYPAR9,MYPAR10,MYPAR11,MYPAR12,
     &       ANGPARX,ESPAR,QEANGPAR,QEPAR,DELANGPAR,DELPAR, 
     &       DISANGPAR,DISPAR,R1PAR,R2PAR
      COMMON/MYPAR/MYPAR1,MYPAR2,MYPAR3,MYPAR4,MYPAR5,MYPAR6,
     &             MYPAR7,MYPAR8,MYPAR9,MYPAR10,MYPAR11,MYPAR12,
     &             ANGPARX,ESPAR,QEANGPAR,QEPAR,DELANGPAR,DELPAR, 
     &             DISANGPAR,DISPAR,R1PAR,R2PAR

      PI=ACOS(-1.)
      THR=TH*PI/180.
      DM=1232.
      PIMASS=140.
      PM=940.
      A2=550.
      PFR=60.
      GAM2N=20.
      GAMQFR=40.
      GAMREF=300.
      GAMR=GAMREF
      SIGREF=0.20D-7
      QMSR=4.*596.8*(596.8-380.)*SIN(60.*PI/180./2.)**2
      QVSR=QMSR+380.**2
      SIGKIN=0.5*SIGMOT(596.8D0,60.*PI/180.)
      SIGKIN=SIGKIN*(QMSR/2./QVSR+TAN(60.*PI/180./2.)**2)
      SIGKIN=SIGKIN*QVSR*FD(QMSR,A2)**2
      SIGKIN=SIGKIN*GAMR/GAMREF
      SIGCON=SIGREF/SIGKIN
      QMS=4.*E*(E-W)*SIN(THR/2.)**2
      QVS=QMS+W**2
      GAMQF=GAMQFR*(PF/PFR)*(SQRT(QVS)/SQRT(QVSR))
      EFFMASS=(PM+DM)/2.
      SIG=(Z*(A-Z)/A)*SIGMOT(E,THR)
      SIG=SIG*(QMS/2./QVS+TAN(THR/2.)**2)
      SIG=SIG*QVS*FD(QMS,A2)**2
      EKAPPA=W-QMS/2./PM
      CMTOT2=PM**2+2.*PM*EKAPPA
C     GAM=SQRT(GAMR**2+GAMQF**2)
      GAM=GAMR
      SIG=SIG*CMTOT2*GAM**2
      SIG=SIG/((CMTOT2-EFFMASS**2)**2+CMTOT2*GAM**2)
      SIG=SIG*(GAMR/GAM)*SIGCON
      SIG2N=SIG
      WTHRESH=QMS/4./PM
c      
      WTHRESH=0.0
c
      IF(W.GT.WTHRESH)THEN
      THRESH=1.-EXP(-(W-WTHRESH)/GAM2N)
      ELSE
      THRESH=0.
      ENDIF
      IF (A.eq.14.0) THEN
c        SIG2N=SIG2N*THRESH*0.65
C        SIG2N=SIG2N*THRESH*12.0*(1.0+E/W)
C        SIG2N=SIG2N*THRESH*(68.0*E/1000)*(1.0/(E-W))*1E+2*(1.0+E/W)
        SIG2N=SIG2N*THRESH*2.0*(E-W)*1E-2*(1.0+E/W)
      ELSE IF (A.eq.3.0) THEN
C        SIG2N=SIG2N*THRESH*0.9*(E-W)*1E-2*(1.0+E/W)
C        SIG2N=SIG2N*THRESH*(E/5009)*(W*1E-3)
C     &        *(SIN((32.0/2.0)*(PI/180.0)))/SIN(THR/2.0)
C        SIG2N=SIG2N*THRESH*QMS*1E-6*(4018/E)
C        SIG2N=SIG2N*THRESH
C        MYPAR3 = QMS*1E-6
C        MYPAR5 = E 
        SIG2N=SIG2N*THRESH*MYPAR2
      ELSE
        SIG2N=SIG2N*THRESH
      ENDIF
      RETURN
      END
C-------------------------------------------------------------------
*SIGDEL
      REAL*8 FUNCTION SIGDEL(E,TH,W,A,EPSD,PF)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 MYPAR1,MYPAR2,MYPAR3,MYPAR4,MYPAR5,MYPAR6,
     &       MYPAR7,MYPAR8,MYPAR9,MYPAR10,MYPAR11,MYPAR12,
     &       ANGPARX,ESPAR,QEANGPAR,QEPAR,DELANGPAR,DELPAR, 
     &       DISANGPAR,DISPAR,R1PAR,R2PAR
      COMMON/MYPAR/MYPAR1,MYPAR2,MYPAR3,MYPAR4,MYPAR5,MYPAR6,
     &             MYPAR7,MYPAR8,MYPAR9,MYPAR10,MYPAR11,MYPAR12,
     &             ANGPARX,ESPAR,QEANGPAR,QEPAR,DELANGPAR,DELPAR,
     &             DISANGPAR,DISPAR,R1PAR,R2PAR

      PM     = 939.
      PIMASS = 140.
c      DM     = 1219.
      DM     = 1232.
      AD1    = 685.
      AD0    = 774.
      PI     = ACOS(-1.)
      ALPH   = 1./137.03604
      HBARC  = 197.32858
      GAMDP  = 110.
      GAMSPRD= 140.
      GAMR   = 120.
      GAMPI  = 5.
      QFDP   = 1.02D-7
      PFR    = 230.
      QMSR   = 4.*730.*(730.-390.)*SIN(37.1*PI/180./2.)**2
      QVSR   = QMSR+390.**2
      QMSRQ  = 4.*730.*(730.-115.)*SIN(37.1*PI/180./2.)**2
      QVSRQ  = QMSRQ+115.**2

      NA=INT(A)
      IF(NA.EQ.1)THEN
        QFD=QFDP
        GSPRDA=0.
        AD=AD0
      ELSEIF(NA.LT.4)THEN
        QFD=QFDP
        GSPRDA=(A-1.)*GAMSPRD/3.
        AD=AD0+(A-1.)*(AD1-AD0)/3.
      ELSE
        AD=AD1
        GSPRDA=GAMSPRD
        QFD=QFDP
      ENDIF
      THR = TH*PI/180.
      QMS = 4.*E*(E-W)*SIN(THR/2.)**2
      QVS = QMS+W**2
      EKAPPA = W-QMS/2./PM
      CMTOT2 = PM**2+2.*PM*EKAPPA
C  BEGIN DELTA CALCULATION
      IF(NA.GT.1)THEN
        GAMQ=GAMR*PF*SQRT(QVS)/PFR/SQRT(QVSRQ)
      ELSE
        GAMQ=0.
      ENDIF

      EPD = E-(DM-PM)*(DM+PM)/2./PM
      EPD = EPD/(1.+2.*E*SIN(THR/2.)**2/PM)
      EPD = EPD-EPSD
      WD  = E-EPD
      QMSPK = 4.*E*EPD*SIN(THR/2.)**2
      QVSPK = QMSPK+WD**2
C
C NOTE WIDTH INCLUDES E-DEPENDENCE,FERMI BROADENING,& SPREADING
C
      WTHRESH = 4.*E**2*SIN(THR/2.)**2+PIMASS**2+2.*PIMASS*PM
      WTHRESH = WTHRESH/2./PM
      THRESHD = 1.+PF/PM+PF**2/2./PM**2+2.*E*SIN(THR/2.)**2/PM
      WTHRESH = WTHRESH/THRESHD
c
      WTHRESH = 0.0
c
      IF(W.GT.WTHRESH)THEN
        THRESH=1.-EXP(-(W-WTHRESH)/GAMPI)
      ELSE
        THRESH=0.
      ENDIF
      GAMD = GAMDP
      GAM  = SQRT(GAMD**2+GAMQ**2+GSPRDA**2)
      SIGD = QFDP*(GAMDP/GAM)
      SIGD = SIGD*CMTOT2*GAM**2
      SIGD = SIGD/((CMTOT2-(DM+EPSD)**2)**2+CMTOT2*GAM**2)
      SIGD = SIGD*FD(QMS,AD)**2/FD(QMSR,AD)**2
      TEST = QVS/QVSR
      SIGD = SIGD*TEST
      SIGD = SIGD*(QMS/2./QVS+TAN(THR/2.)**2)
      SIGD = SIGD/(QMSR/2./QVSR+TAN(37.1*PI/180./2.)**2)
      SIGD = SIGD*SIGMOT(E,THR)/SIGMOT(730.D0,37.1*PI/180.)
      SIGD = SIGD*A
      SIGD = SIGD*THRESH
      IF (A.eq.14.0) THEN
C        SIGDEL = SIGD*( 2.5*E/(4000.0)+ (E-W) )/W
C        SIGDEL = SIGD*(0.1-2.65*(E-W)*1E-3)
        SIGDEL = SIGD*(2.0+0.75*(E-W)*1E-3)
      ELSE IF (A.eq.3.0) THEN
C        SIGDEL = SIGD*( 4.*E/(4000.0)+ (E-W) )/W
C        SIGDEL = SIGD*0.4
        SIGDEL = SIGD*MYPAR3*DELPAR
      ELSE
        SIGDEL = SIGD
      ENDIF
      RETURN
      END
C-------------------------------------------------------------------
*SIGQFS
      REAL*8 FUNCTION SIGQFS(E,TH,W,Z,A,EPS,PF)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 MYPAR1,MYPAR2,MYPAR3,MYPAR4,MYPAR5,MYPAR6,
     &       MYPAR7,MYPAR8,MYPAR9,MYPAR10,MYPAR11,MYPAR12,
     &       ANGPARX,ESPAR,QEANGPAR,QEPAR,DELANGPAR,DELPAR, 
     &       DISANGPAR,DISPAR,R1PAR,R2PAR
      COMMON/MYPAR/MYPAR1,MYPAR2,MYPAR3,MYPAR4,MYPAR5,MYPAR6,
     &             MYPAR7,MYPAR8,MYPAR9,MYPAR10,MYPAR11,MYPAR12,
     &             ANGPARX,ESPAR,QEANGPAR,QEPAR,DELANGPAR,DELPAR, 
     &             DISANGPAR,DISPAR,R1PAR,R2PAR
      LOGICAL Qdep
      COMMON/QDEPENDENCE/XPAR0,XPAR1,Qdep
      PM   = 939.
      UP   = 2.7928456
      UN   = -1.91304184
      AP0  = 840.
      AP1  = 750.
      ALPH = 1./137.03604
      HBARC= 197.32858
      PI   = ACOS(-1.)
      GAMR = 120.
      PFR  = 230.
      QMSRQ= 4.*730.*(730.-115.)*SIN(37.1*PI/180./2.)**2
      QVSRQ= QMSRQ+115.**2
      NA   = INT(A)

      IF(NA.EQ.1)THEN
        AP=AP0
      ELSEIF(NA.LT.4)THEN
        AP=AP0+(A-1.)*(AP1-AP0)/3.
      ELSE
        AP=AP1
      ENDIF
      THR = TH*PI/180.
      QMS = 4.*E*(E-W)*SIN(THR/2.)**2
      QVS = QMS+W**2

C  START QFS SECTION
      SIGNS  = SIGMOT(E,THR)*RECOIL(E,THR,PM)

      SIGEP = GEP(QMS,AP)**2 + TAU(QMS) * GMP(QMS,AP)**2
      SIGEP = SIGEP/(1.0+TAU(QMS) )
      SIGEP = SIGEP+2.0*TAU(QMS)*GMP(QMS,AP)**2 * (TAN(THR/2.))**2
      SIGEP = SIGNS*SIGEP

      SIGEN = GEN(QMS,AP)**2 + TAU(QMS) * GMN(QMS,AP)**2
      SIGEN = SIGEN/(1.0+TAU(QMS) )
      SIGEN = SIGEN+2.0*TAU(QMS)*GMN(QMS,AP)**2 * (TAN(THR/2.))**2
      SIGEN = SIGNS*SIGEN

      EPQ    = 4.*E**2*SIN(THR/2.)**2/2./PM
      EPQ    = EPQ/(1.+2.*E*SIN(THR/2.)**2/PM)+EPS
      EPQ    = E-EPQ

CDEBUG   QSEPQ=4.*E*EPQ*SIN(THR/2.)**2  ! Q^2 at Q.E. Peak

      IF(INT(A).EQ.1)THEN
        ARG = (E-W-EPQ)/SQRT(2.)/1.
        DEN = 2.51
      ELSE
        GAMQ = GAMR*PF*SQRT(QVS)/PFR/SQRT(QVSRQ)
        ARG  = (E-W-EPQ)/1.20/(GAMQ/2.)
        DEN  = 2.13*(GAMQ/2.)
      ENDIF
      NQ = INT(ARG)
      IF(ABS(NQ).GT.10)THEN
        SIGQ = 0.
      ELSE
        SIGQ = (Z*SIGEP+(A-Z)*SIGEN)*EXP(-ARG**2)/DEN
      ENDIF
      IF (A.eq.14.0) THEN
        SIGQFS=SIGQ*E*(E/4000.0)*(1/1178.0)
      ELSE IF (A.eq.3.0) THEN
C        SIGQFS=SIGQ*(E/W)*(E/4000.0)
C        SIGQFS=SIGQ*(SIN((32.0/2.0)*(PI/180.0)))/SIN(THR/2.0)
C        SIGQFS=SIGQ*1.2
         SIGQFS=SIGQ*MYPAR1*QEPAR 
      ELSE
        SIGQFS=SIGQ
      ENDIF
CDEBUG Q2 DEPENDENCE 
      IF (Qdep) THEN
        QQ    = QMS/1.E6
        XCOR  = XPAR1/QQ + XPAR0   
        SIGQFS= SIGQFS/XCOR
      ENDIF

      RETURN
      END
C-------------------------------------------------------------------
*     RADIATE
      SUBROUTINE RADIATE(E,TH,W,SIGNR,SIGRAD)
C     DOES NOT INCLUDE CONTRIBUTION FROM ELASTIC SCATTERING
C
C-----K. SlIFer. 09/16/02
C
C     Rewrote subroutine to include external bremsstrahlung using
C     formalism from S. Stein et al. Phys. Rev. D 12 7. Equation (A82)
C     Where possible the equation number is given with each expression.
C----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL extrad
      INTEGER NUM 
      COMMON/PAR/E0,TH0,W0,Z,A,EPS,EPSD,PF,SPENCE
      COMMON/ADDONS/SCALE,Tb,Ta,extrad
      EXTERNAL FES,FEP

C      REAL ESMIN 

      DEL   = 10
      PREC  = .0001  ! 0.001 gets rid of glitches

C     For Simpson integration
      NUM=80
      STEP=1.0
 
      IF (.NOT.extrad) THEN  ! don't apply external rad. cor.
        Tb=0.0
        Ta=0.0
      ENDIF
      xb    = 4./3.  ! A45
      XM    = 931.49 ! Mass of the nucleon
      XMT   = A*XM   ! Mass of the target
      ALPH  = 1./137.03604
      EMASS = 0.511
      PI    = ACOS(-1.)
      THR   = TH*PI/180.
      ARG   = COS(THR/2.)**2

      SPENCE= PI**2/6.-LOG(ARG)*LOG(1.-ARG)
      DO 10 NSP=1,50
 10     SPENCE = SPENCE-ARG**NSP/FLOAT(NSP)**2

      QMS   = 4.*E*(E-W)*SIN(THR/2.)**2

      D1=(2.*ALPH/PI)*(LOG(QMS/EMASS**2)-1.)      ! =2b*tr (A57)
      tr=D1/2./xb

      D2 = 13.*(LOG(QMS/EMASS**2)-1.)/12.-17./36. ! this term dominates D2
      D2 = D2 +0.5*(PI**2/6.-SPENCE)
      D2 = D2 -1./4.*( LOG( E/(E-W) ) )**2        ! Correct. to peak. appr.
      D2 = D2*(2.*ALPH/PI)                 
      D2 = D2+0.5772*xb*(Tb+Ta)                   ! Here D2= F-1
      xF = (1.+D2)                                ! (A44)

      Tpb = tr + Tb
      Tpa = tr + Ta  
   
      R   = ( XMT+E*(1-COS(THR)) )/( XMT-(E-W)*(1-COS(THR)) ) ! (A83)
      eta = LOG(1440.*Z**(-2./3.) )/LOG(183.*Z**(-1./3.) )    ! (A46)
      xi  = (PI*EMASS/2./ALPH)*(Ta+Tb)
      xi  = xi/( (Z+eta)*LOG(183.*Z**(-1./3.)) )              ! (A52)

      SIGRAD = SIGNR * xF
      SIGRAD = SIGRAD*( (R*DEL/E  )**(xb*Tpb) )
      SIGRAD = SIGRAD*( (DEL/(E-W))**(xb*Tpa) )
      SIGRAD = SIGRAD*(1. - xi/DEL/( 1.-xb*(Tpb+Tpa)) )

      TERM1=(R*DEL/E  )**(xb*Tpb)
      TERM2=(DEL/(E-W))**(xb*Tpa)
      TERM3=(1. - xi/DEL/( 1.-xb*(Tpb+Tpa)) )
      TERM4=xF

C
C-----Stein's 1st integral wrt dEs' (A82)
C
C     limits of 0 to W-DEL give almost same results
C
      X1LO   = (E-W)*( XMT/( XMT-2.*(E-W)*(SIN(THR/2.))**2) -1.0 )
      X1HI   = W-R*DEL

      ANS_Es = 0.
      IF (X1HI.GT.X1LO) THEN
C        CALL ROM(X1LO,X1HI,PREC,ANS_Es,KF,1)
        CALL SIM(X1LO,X1HI,STEP,NUM,FES,ANS_Es)
      ELSE
        WRITE(6,*) "Integral dEs. SKIPPING:nu,lower,higher ",W,X1LO,X1HI
      ENDIF
      ANS_Es = ANS_Es * SCALE
C
C-----Stein's 2nd integral wrt dEp' (A82)
C
C     limits of 0 to W-DEL give almost same results
C

      X2LO   = E*( 1.0-XMT/( XMT+2.*E*(SIN(THR/2.))**2) )
      X2HI   = W-DEL 

      ANS_Ep = 0.
      IF (X2HI.GT.X2LO) THEN
C        CALL ROM(X2LO,X2HI,PREC,ANS_Ep,KF,2)
        CALL SIM(X2LO,X2HI,STEP,NUM,FEP,ANS_Ep)
      ELSE
         WRITE(6,*) "Integral dEp. SKIPPING:nu,lower,higher ",W,X2LO,
     &              X2HI
      ENDIF
      ANS_Ep = ANS_Ep * SCALE

CCDEBUG
C      WRITE(6,'(F7.1,3F8.2)') W,(SIGRAD-SIGNR)/SIGNR*100., 
C     &                          ANS_Es/SIGNR*100.,
C     &                          ANS_Ep/SIGNR*100.
C

      SIGRAD = SIGRAD + ANS_Es + ANS_Ep

C      IF( (E-W).EQ.600.0) THEN 
C        WRITE(*,*) E-W,' ',X1LO+E-W,' ',ESMIN,' ',X1HI+E-W,' ',E-X2LO,
C     &             ' ',E-X2HI 
C      ENDIF

C      WRITE(*,*) E-W,'  ',TERM1*TERM2*TERM3*TERM4,'  ',
C     &           ANS_Es,'  ',ANS_Ep

CDEBUG
C      WRITE(6,'(F7.1,6F8.2)') 
C     &         W,TERM1,TERM2,TERM3,TERM4,
C     &         ANS_Es/SIGRAD*100.0,ANS_Ep/SIGRAD*100.0

   
      RETURN
      END
C-------------------------------------------------------------------
*VALY
      SUBROUTINE VALY(X,F,IFUNC)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL extrad
C     REAL*8 FUNCTION F(X)
      COMMON/PAR/E,TH,W,Z,A,EPS,EPSD,PF,SPENCE
      COMMON/ADDONS/SCALE,Tb,Ta,extrad

      IF (.NOT.extrad) THEN  ! apply external rad. cor.
        Tb=0.0
        Ta=0.0
      ENDIF
      ALPH  = 1./137.03604
      EMASS = 0.511
      PI    = ACOS(-1.)
      THR   = TH*PI/180.
      xb    = 4./3.
      XM    = 931.49 ! Mass of the nucleon
      XMT   = A*XM   ! Mass of the target

      eta = LOG(1440.*Z**(-2./3.) )/LOG(183.*Z**(-1./3.) )
      xi  = (PI*EMASS/2./ALPH)*(Ta+Tb)
      xi  = xi/( (Z+eta)*LOG(183.*Z**(-1./3.)) )
      R   = ( XMT+E*(1-COS(THR)) )/( XMT-(E-W)*(1-COS(THR)) )
C
C-----Stein's 2nd integral dEp'
C
      QMS2  = 4.*E*(E-X)*SIN(THR/2.)**2  ! 1/15/03
      tr2   = 1./xb*(ALPH/PI)*(LOG(QMS2/EMASS**2)-1.)
      Tpb   = tr2 + Tb
      Tpa   = tr2 + Ta

      D2    = 13.*(LOG(QMS2/EMASS**2)-1.)/12.-17./36.
      D2    = D2 - 1./4.*( LOG( E/(E-W) ) )**2 !KS. Correction to peak. approx.
      D2    = D2 + 0.5*(PI**2/6.-SPENCE)
      D2    = D2 * (2.*ALPH/PI)
      D2    = D2 + 0.5772*xb*(Tb+Ta)

      SIG2  = SIGQFS(E,TH,X,Z,A,EPS,PF)
      SIG2  = SIG2 + SIGDEL(E,TH,X,A,EPSD,PF)
      SIG2  = SIG2 + SIGX(E,TH,X,A)
      SIG2  = SIG2 + SIGR1(E,TH,X,A,PF)
      SIG2  = SIG2 + SIGR2(E,TH,X,A,PF)
      SIG2  = SIG2 + SIG2N(E,TH,X,Z,A,PF)
      F2    = ( xb*Tpa/(W-X) ) *phi((W-X)/(E-X))
      F2    = F2 + xi/(2.*(W-X)**2)
      F2    = F2 * SIG2*(1.+D2)
      F2    = F2 * ( (W-X)/(E-X) )**(xb*Tpa)
      F2    = F2 * ( (W-X)*R/(E) )**(xb*Tpb)
C
C-----Stein's 1st integral dEs'
C
      QMS1  = 4.*(E-W+X)*(E-W)*SIN(THR/2.)**2    !    1/15/03
      tr1   = 1./xb*(ALPH/PI)*(LOG(QMS1/EMASS**2)-1.)
      Tpb   = tr1 + Tb
      Tpa   = tr1 + Ta

      D2    = 13.*(LOG(QMS1/EMASS**2)-1.)/12.-17./36.
      D2    = D2 - 1./4.*( LOG( E/(E-W) ) )**2 !Corr. to peak. approx.
      D2    = D2 + 0.5*(PI**2/6.-SPENCE)
      D2    = D2 * (2.*ALPH/PI)
      D2    = D2 + 0.5772*xb*(Tb+Ta) ! 1/14/02

      SIG1  = SIGQFS(E-W+X,TH,X,Z,A,EPS,PF)
      SIG1  = SIG1 + SIGDEL(E-W+X,TH,X,A,EPSD,PF)
      SIG1  = SIG1 +   SIGX(E-W+X,TH,X,A)
      SIG1  = SIG1 +  SIGR1(E-W+X,TH,X,A,PF)
      SIG1  = SIG1 +  SIGR2(E-W+X,TH,X,A,PF)
      SIG1  = SIG1 +  SIG2N(E-W+X,TH,X,Z,A,PF)

      F1    = ( xb*Tpb/(W-X) ) *phi((W-X)/(E))   ! 
      F1    = F1 + xi/(2.*(W-X)**2)
      F1    = F1 * SIG1*(1.+D2)
      F1    = F1 * ( (W-X)/((E-W)*R) )**(xb*Tpa)
      F1    = F1 * ( (W-X)/ (E)      )**(xb*Tpb) ! 

      IF(IFUNC.EQ.2) THEN      ! dEp'
        F=F2
      ELSEIF (IFUNC.EQ.1) THEN ! dEs'
        F=F1
      ELSE
        WRITE(6,*) "PROBLEM. "
        STOP
      ENDIF
      RETURN
      END
C-------------------------------------------------------------------
*  ROM 
      SUBROUTINE ROM(A,B,EPS,ANS,K,IFUNC)
      IMPLICIT REAL*8 (A-H,O-Z)
C  ROMBERG METHOD OF INTEGRATION
      DIMENSION W(50,50)
      H=B-A
      K=0
      CALL VALY(A,FA,IFUNC)
      CALL VALY(B,FB,IFUNC)
      W(1,1)=(FA+FB)*H/2.
    4 K=K+1
      IF(K.GE.49)GO TO 5
      H=H/2.
      SIG=0.
      M=2**(K-1)
      DO 1 J=1,M
      J1=2*J-1
      X=A+FLOAT(J1)*H
      CALL VALY(X,F,IFUNC)
C      WRITE(6,*) "DEBUG: ",k,IFUNC,f
    1 SIG=SIG+F
      W(K+1,1)=W(K,1)/2.+H*SIG
      DO 2 L=1,K
      IU=K+1-L
      IV=L+1
    2 W(IU,IV)=(4.**(IV-1)*W(IU+1,IV-1)-W(IU,IV-1))/(4.**(IV-1)-1.)
      E=(W(IU,IV)-W(IU,IV-1))/W(IU,IV)
      IF(ABS(E).LT.EPS) GO TO 3
      IF(ABS(E).EQ.EPS) GO TO 3
      IF(ABS(E).GT.EPS) GO TO 4
    3 ANS=W(1,IV)
      RETURN
    5 PRINT 100
  100 FORMAT(' K OVERFLOW')
      CALL EXIT(0)
      END
C-------------------------------------------------------------------
      REAL*8 FUNCTION phi(x)
      IMPLICIT REAL*8 (A-H,O-Z)
      phi=1.0-x+3./4.*x**2
      RETURN
      END
C-------------------------------------------------------------------
      SUBROUTINE OPENFILES(TAG,INTERACTIVE,RATES,nucleus)
      CHARACTER TAG*4,nucleus*3
      LOGICAL INTERACTIVE,RATES
 
      IF (nucleus.EQ.'He3') THEN ! He-3
         OPEN(8,FILE=TAG//'_He_sig_qfs.out',STATUS='UNKNOWN')
      ELSEIF (RATES) THEN
         OPEN(10,FILE=TAG//'_He_rates.out',STATUS='UNKNOWN')
      ELSEIF (nucleus.EQ.'Nit') THEN ! Nitrogen
         OPEN(8,FILE=TAG//'_N2_sig_qfs.out',STATUS='UNKNOWN')
      ELSEIF (nucleus.EQ.'Sil') THEN ! Silicon
         OPEN(8,FILE=TAG//'_Si_sig_qfs.out',STATUS='UNKNOWN')
      ELSEIF (nucleus.EQ.'Irn') THEN ! Iron
         OPEN(8,FILE=TAG//'_Fe_sig_qfs.out',STATUS='UNKNOWN')
      ELSEIF (INTERACTIVE) THEN ! will prompt user for input
         OPEN(8,FILE='sig_qfs.out',STATUS='UNKNOWN')
C        OPEN(18,FILE='sig_qfs.fak',STATUS='UNKNOWN')
      ELSE
        WRITE(6,*) 'Unknown target '
        STOP
      ENDIF
      RETURN
      END
C-------------------------------------------------------------------
      SUBROUTINE GETSCALE(units,SCALE,UNSCALE)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER units*2

C-----DETERMINE DESIRED UNITS OF CROSS SECTION
      IF(units.EQ.'cm') THEN     ! cm^2/(MeV sr)
        SCALE=1.D-26
        UNSCALE=1.0
      ELSEIF(units.EQ.'mb')THEN  ! mb/(MeV sr)
        SCALE=1.D+01
        UNSCALE=1.D-27
      ELSEIF(units.EQ.'ub')THEN  ! ub/(MeV sr)
        SCALE=1.D+04
        UNSCALE=1.D-30
      ELSEIF(units.EQ.'nb')THEN  ! nb/(MeV sr)
        SCALE=1.D+07
        UNSCALE=1.D-33
      ELSEIF(units.EQ.'pb') THEN ! pb/(MeV sr)
        SCALE=1.D+10
        UNSCALE=1.D-36
      ELSE
        WRITE(6,*) 'Unknown scale'
        STOP
      ENDIF
C      WRITE(8,'(A,A)') '#',units
      RETURN
      END
C-------------------------------------------------------------------
      SUBROUTINE GETZA(nucleus,INTERACTIVE,Z,A)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER nucleus*3 
      LOGICAL INTERACTIVE
C-----DETERMINE nucleus charge, mass
      IF (nucleus.EQ.'He3') THEN ! He-3
        Z=2
        A=3
      ELSEIF (nucleus.EQ.'Nit') THEN ! Nitrogen
        Z=7
        A=14
      ELSEIF (nucleus.EQ.'Sil') THEN ! Silicon
        Z=14
        A=28
      ELSEIF (nucleus.EQ.'Irn') THEN ! Iron
        Z=26
        A=56
      ELSEIF (INTERACTIVE) THEN
        CONTINUE
      ELSE
        WRITE(6,*) 'Problem'
        stop
      ENDIF
C------
      RETURN
      END
C-------------------------------------------------------------------
      SUBROUTINE GETQDEP(A,Qdep,XPAR0,XPAR1)
C     Introduce Q-dependent correction to Q.E. Peak height
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL Qdep
     
      REAL*8 DUM

      IF(Qdep) DUM = 1.*A 
        
      XPAR1 = 0.0453828 ! GeV^2  ! 12/13/02 All available CARBON(only) FIT
      XPAR0 = 0.686704  ! GeV^2

      ! XPAR1 = 0.0683579 ! GeV^2  ! He3 Fit. Bates data, T.S. Ueng Thesis
      ! XPAR0 = 0.757269  ! GeV^2

      ! XPAR1 = 0.0476549 ! GeV^2 ! He3 Fit. E94010 unfolded xs.
      ! XPAR0 = 0.801282  ! GeV^2

      RETURN
      END
C-------------------------------------------------------------------
      SUBROUTINE TARGETCELL(TAG,INTERACTIVE,E,WMIN,WMAX,DELTAW,Tb,Ta,
     &                      density,zlength)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER TAG*4
      LOGICAL INTERACTIVE

C-----DETERMINE TARGET CELL CHARACTERISTICS
C----- CSR Targets     -----------------------------------------------
      IF (TAG.EQ.'CSR6') THEN 
        E      = 646.0 
        WMIN   = 6.0
        WMAX   = 496.0 
        DELTAW = 10.0 
        Tb     = 1.000E-02 
        Ta     = 1.000E-02 
        density= 1.0
        zlength= 39.0 
      ELSEIF (TAG.EQ.'CSR7') THEN 
        E      = 739.0 
        WMIN   = 19.0
        WMAX   = 589.0 
        DELTAW = 10.0 
        Tb     = 1.000E-02 
        Ta     = 1.000E-02 
        density= 1.0
        zlength= 39.0
       ELSEIF (TAG.EQ.'CSR8') THEN 
        E      = 844.0 
        WMIN   = 54.0
        WMAX   = 694.0 
        DELTAW = 10.0 
        Tb     = 1.000E-02 
        Ta     = 1.000E-02 
        density= 1.0
        zlength= 39.0
       ELSEIF (TAG.EQ.'CSR9') THEN 
        E      = 957.0 
        WMIN   = 107.0
        WMAX   = 807.0 
        DELTAW = 10.0 
        Tb     = 1.000E-02 
        Ta     = 1.000E-02 
        density= 1.0
        zlength= 39.0 
C----- E01-012 Targets -----------------------------------------------
      ELSEIF (TAG.EQ.'4_DU') THEN   ! Target is "DUKE" 
        E      = 4017.9         ! Incident Energy in MeV
        WMIN   = 1105.69            
        WMAX   = 2267.69
        DELTAW = 1.0
        Tb     = 2.540E-3       ! Thickness before scatt., including 1/2 target
        Ta     = 5.435E-2       ! Thickness after scatt., including targ contribution
        density= 11.6           ! in Amagats
        zlength = 39.4          ! in cm
      ELSEIF (TAG.EQ.'5_DU') THEN ! DUKE   
        E      = 5009.0          
        WMIN   = 1580.
        WMAX   = 2930.
        DELTAW = 1.0
        Tb     = 2.540E-3      
        Ta     = 5.480E-2
        density=  11.6        
        zlength = 39.4       
      ELSEIF(TAG.EQ.'6_DU') THEN ! DUKE 
        E      = 5009.0
        WMIN   = 2268.50
        WMAX   = 3073.50 
        DELTAW = 1.0
        Tb     = 2.540E-3
        Ta     = 4.480E-2
        density= 11.6  
        zlength = 39.4 
      ELSEIF(TAG.EQ.'1_DU') THEN ! DUKE
        E      = 1046.0
        WMIN   = 0.
        WMAX   = 1050.0
        DELTAW = 5.0
        Tb     = 2.540E-3
        Ta     = 7.935E-2
        density= 11.6    
        zlength = 39.4   
      ELSEIF(TAG.EQ.'3_EX') THEN ! EXODUS
        E      = 3028.1
        WMIN   = 600.
        WMAX   = 2250.
        DELTAW = 10.0
        Tb     = 2.460E-3
        Ta     = 4.660E-2
        density= 12.0          
        zlength = 39.6
      ELSEIF(TAG.EQ.'6_EX') THEN ! EXODUS
        E      = 5009.0
        WMIN   = 2170.
        WMAX   = 3280.0
        DELTAW = 10.0
        Tb     = 2.460E-3
        Ta     = 3.870E-2
        density= 12.0
        zlength = 39.6        
C---------- E06-014 Targets ----------------------------------
      ELSEIF(TAG.EQ.'1_SA') THEN  ! SAMANTHA (d2n)
        E      = 1000.0               ! Incident energy (MeV) 
        WMIN   = 100.0
        WMAX   = 900.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 
        Ta     = 3.623E-02            ! Thickness in #X0 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]   
      ELSEIF(TAG.EQ.'1_S0') THEN  ! SAMANTHA (d2n)
        E      = 1500.0               ! Incident energy (MeV) 
        WMIN   = 530.0
        WMAX   = 1200.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 
        Ta     = 3.623E-02            ! Thickness in #X0 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]  
      ELSEIF(TAG.EQ.'1_T0') THEN  ! SAMANTHA (d2n)
        E      = 1575.0               ! Incident energy (MeV) 
        WMIN   = 230.0
C        WMIN   = 530.0
        WMAX   = 1400.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 
        Ta     = 3.623E-02            ! Thickness in #X0 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]  
      ELSEIF(TAG.EQ.'1_T1') THEN  ! SAMANTHA (d2n)
        E      = 1700.0               ! Incident energy (MeV) 
        WMIN   = 230.0
C        WMIN   = 530.0
        WMAX   = 1600.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 
        Ta     = 3.623E-02            ! Thickness in #X0 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]  
      ELSEIF(TAG.EQ.'1_S1') THEN  ! SAMANTHA (d2n)
        E      = 1625.0               ! Incident energy (MeV) 
        WMIN   = 250.0
        WMAX   = 1300.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 
        Ta     = 3.623E-02            ! Thickness in #X0 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]   
      ELSEIF(TAG.EQ.'1_S2') THEN  ! SAMANTHA (d2n)
        E      = 1750.0               ! Incident energy (MeV) 
        WMIN   = 290.0
        WMAX   = 1450.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 
        Ta     = 3.623E-02            ! Thickness in #X0 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]   
      ELSEIF(TAG.EQ.'1_S3') THEN  ! SAMANTHA (d2n)
        E      = 1875.0               ! Incident energy (MeV) 
        WMIN   = 400.0
        WMAX   = 1550.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 
        Ta     = 3.623E-02            ! Thickness in #X0 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]   
      ELSEIF(TAG.EQ.'2_SA') THEN  ! SAMANTHA (d2n)
        E      = 2000.0               ! Incident energy (MeV) 
        WMIN   = 800.0
        WMAX   = 1820.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 
        Ta     = 3.623E-02            ! Thickness in #X0 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]  
      ELSEIF(TAG.EQ.'2_S0') THEN  ! SAMANTHA (d2n)
        E      = 2125.0               ! Incident energy (MeV) 
        WMIN   = 650.0
        WMAX   = 1820.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 
        Ta     = 3.623E-02            ! Thickness in #X0 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]  
      ELSEIF(TAG.EQ.'2_S1') THEN  ! SAMANTHA (d2n)
        E      = 2250.0               ! Incident energy (MeV) 
        WMIN   = 740.0
        WMAX   = 1870.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 
        Ta     = 3.623E-02            ! Thickness in #X0 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]  
      ELSEIF(TAG.EQ.'2_S2') THEN  ! SAMANTHA (d2n)
        E      = 2375.0               ! Incident energy (MeV) 
        WMIN   = 820.0
        WMAX   = 1870.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 
        Ta     = 3.623E-02            ! Thickness in #X0 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]  
      ELSEIF(TAG.EQ.'2_S3') THEN  ! SAMANTHA (d2n)
        E      = 2500.0               ! Incident energy (MeV) 
        WMIN   = 920.0
        WMAX   = 2220.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 
        Ta     = 3.623E-02            ! Thickness in #X0 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]  
      ELSEIF(TAG.EQ.'2_S4') THEN  ! SAMANTHA (d2n)
        E      = 2625.0               ! Incident energy (MeV) 
        WMIN   = 1000.0
        WMAX   = 2320.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 
        Ta     = 3.623E-02            ! Thickness in #X0 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]  
      ELSEIF(TAG.EQ.'2_S5') THEN  ! SAMANTHA (d2n)
        E      = 2750.0               ! Incident energy (MeV) 
        WMIN   = 1070.0
        WMAX   = 2470.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 
        Ta     = 3.623E-02            ! Thickness in #X0 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]  
      ELSEIF(TAG.EQ.'2_S6') THEN  ! SAMANTHA (d2n)
        E      = 2875.0               ! Incident energy (MeV) 
        WMIN   = 1190.0
        WMAX   = 2470.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 
        Ta     = 3.623E-02            ! Thickness in #X0 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]  
       ELSEIF(TAG.EQ.'3_SA') THEN  ! SAMANTHA (d2n)
        E      = 3000.0               ! Incident energy (MeV) 
        WMIN   = 1480.0
        WMAX   = 2730.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 
        Ta     = 3.623E-02            ! Thickness in #X0 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]   
       ELSEIF(TAG.EQ.'3_S0') THEN  ! SAMANTHA (d2n)
        E      = 3125.0               ! Incident energy (MeV) 
        WMIN   = 1380.0
        WMAX   = 2830.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 
        Ta     = 3.623E-02            ! Thickness in #X0 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]   
       ELSEIF(TAG.EQ.'3_S1') THEN  ! SAMANTHA (d2n)
        E      = 3250.0               ! Incident energy (MeV) 
        WMIN   = 1470.0
        WMAX   = 3090.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 
        Ta     = 3.623E-02            ! Thickness in #X0 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]   
       ELSEIF(TAG.EQ.'3_S2') THEN  ! SAMANTHA (d2n)
        E      = 3375.0               ! Incident energy (MeV) 
        WMIN   = 1570.0
        WMAX   = 3190.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 
        Ta     = 3.623E-02            ! Thickness in #X0 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]   
       ELSEIF(TAG.EQ.'3_S3') THEN  ! SAMANTHA (d2n)
        E      = 3500.0               ! Incident energy (MeV) 
        WMIN   = 1670.0
        WMAX   = 3200.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 
        Ta     = 3.623E-02            ! Thickness in #X0 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]   
       ELSEIF(TAG.EQ.'3_S4') THEN  ! SAMANTHA (d2n)
        E      = 3625.0               ! Incident energy (MeV) 
        WMIN   = 1770.0
        WMAX   = 3200.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 
        Ta     = 3.623E-02            ! Thickness in #X0 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]   
       ELSEIF(TAG.EQ.'3_S5') THEN  ! SAMANTHA (d2n)
        E      = 3750.0               ! Incident energy (MeV) 
        WMIN   = 1870.0
        WMAX   = 3470.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 
        Ta     = 3.623E-02            ! Thickness in #X0 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]   
       ELSEIF(TAG.EQ.'3_S6') THEN  ! SAMANTHA (d2n)
        E      = 3875.0               ! Incident energy (MeV) 
        WMIN   = 1970.0
        WMAX   = 3470.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 
        Ta     = 3.623E-02            ! Thickness in #X0 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]   
      ELSEIF(TAG.EQ.'4_SA') THEN  ! SAMANTHA (d2n)
        E      = 4730.0               ! Incident energy (MeV) 
        WMIN   = 2840.0
        WMAX   = 4150.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 
        Ta     = 3.623E-02            ! Thickness in #X0 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]   
      ELSEIF(TAG.EQ.'4_S0') THEN      ! SAMANTHA [3He] (d2n)
        E      = 4000.0               ! Incident energy (MeV) 
        WMIN   = 2344.0
        WMAX   = 3421.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 [all materials before interaction]
        Ta     = 3.623E-02            ! Thickness in #X0 [all materials after interaction ]
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]       
      ELSEIF(TAG.EQ.'4_S1') THEN      ! SAMANTHA [3He] (d2n)
        E      = 4125.0               ! Incident energy (MeV) 
        WMIN   = 2344.0
        WMAX   = 3421.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 [all materials before interaction]
        Ta     = 3.623E-02            ! Thickness in #X0 [all materials after interaction ]
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]       
      ELSEIF(TAG.EQ.'4_S2') THEN      ! SAMANTHA [3He] (d2n)
        E      = 4250.0               ! Incident energy (MeV) 
        WMIN   = 2464.0
        WMAX   = 3721.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 [all materials before interaction]
        Ta     = 3.623E-02            ! Thickness in #X0 [all materials after interaction ]
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]       
      ELSEIF(TAG.EQ.'4_S3') THEN      ! SAMANTHA [3He] (d2n)
        E      = 4375.0               ! Incident energy (MeV) 
        WMIN   = 2464.0
        WMAX   = 3721.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 [all materials before interaction]
        Ta     = 3.623E-02            ! Thickness in #X0 [all materials after interaction ]
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]       
      ELSEIF(TAG.EQ.'4_S4') THEN      ! SAMANTHA [3He] (d2n)
        E      = 4500.0               ! Incident energy (MeV) 
        WMIN   = 2664.0
        WMAX   = 3921.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 [all materials before interaction]
        Ta     = 3.623E-02            ! Thickness in #X0 [all materials after interaction ]
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]          
      ELSEIF(TAG.EQ.'4_S5') THEN      ! SAMANTHA [3He] (d2n)
        E      = 4625.0               ! Incident energy (MeV) 
        WMIN   = 2664.0
        WMAX   = 3921.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 [all materials before interaction]
        Ta     = 3.623E-02            ! Thickness in #X0 [all materials after interaction ]
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]          
      ELSEIF(TAG.EQ.'4_S6') THEN      ! SAMANTHA [3He] (d2n)
        E      = 4875.0               ! Incident energy (MeV) 
        WMIN   = 2974.0
        WMAX   = 4600.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 [all materials before interaction]
        Ta     = 3.623E-02            ! Thickness in #X0 [all materials after interaction ]
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]          
      ELSEIF(TAG.EQ.'5_SA') THEN  ! SAMANTHA (d2n)
        E      = 5890.0               ! Incident energy (MeV) 
        WMIN   = 3830.0
        WMAX   = 5311.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 
        Ta     = 3.623E-02            ! Thickness in #X0 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]
      ELSEIF(TAG.EQ.'5_S0') THEN      ! SAMANTHA [3He] (d2n)
        E      = 5000.0               ! Incident energy (MeV) 
        WMIN   = 3060.0
        WMAX   = 4610.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 [all materials before interaction]
        Ta     = 3.623E-02            ! Thickness in #X0 [all materials after interaction ] 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]       
      ELSEIF(TAG.EQ.'TEST') THEN      ! SAMANTHA [3He] (d2n)
        E      = 4243.145             ! Incident energy (MeV) 
        WMIN   = 2543.145
        WMAX   = 3743.145
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 [all materials before interaction]
        Ta     = 3.623E-02            ! Thickness in #X0 [all materials after interaction ] 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]       
      ELSEIF(TAG.EQ.'5_S1') THEN      ! SAMANTHA [3He] (d2n)
        E      = 5125.0               ! Incident energy (MeV) 
        WMIN   = 3060.0
        WMAX   = 4610.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 [all materials before interaction]
        Ta     = 3.623E-02            ! Thickness in #X0 [all materials after interaction ] 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]       
      ELSEIF(TAG.EQ.'5_S2') THEN      ! SAMANTHA [3He] (d2n)
        E      = 5250.0               ! Incident energy (MeV) 
        WMIN   = 3320.0
        WMAX   = 4910.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 [all materials before interaction]
        Ta     = 3.623E-02            ! Thickness in #X0 [all materials after interaction ] 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]       
      ELSEIF(TAG.EQ.'5_S3') THEN      ! SAMANTHA [3He] (d2n)
        E      = 5375.0               ! Incident energy (MeV) 
        WMIN   = 3320.0
        WMAX   = 4910.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 [all materials before interaction]
        Ta     = 3.623E-02            ! Thickness in #X0 [all materials after interaction ] 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]       
      ELSEIF(TAG.EQ.'5_S4') THEN      ! SAMANTHA [3He] (d2n)
        E      = 5500.0               ! Incident energy (MeV) 
        WMIN   = 3580.0
        WMAX   = 5210.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 [all materials before interaction]
        Ta     = 3.623E-02            ! Thickness in #X0 [all materials after interaction ] 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]   
      ELSEIF(TAG.EQ.'5_S5') THEN      ! SAMANTHA [3He] (d2n)
        E      = 5625.0               ! Incident energy (MeV) 
        WMIN   = 3580.0
        WMAX   = 5210.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 [all materials before interaction]
        Ta     = 3.623E-02            ! Thickness in #X0 [all materials after interaction ] 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]   
      ELSEIF(TAG.EQ.'5_S6') THEN      ! SAMANTHA [3He] (d2n)
        E      = 5750.0               ! Incident energy (MeV) 
        WMIN   = 3730.0
        WMAX   = 5210.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 [all materials before interaction]
        Ta     = 3.623E-02            ! Thickness in #X0 [all materials after interaction ] 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]   
      ELSEIF(TAG.EQ.'5_T4') THEN      ! SAMANTHA [3He] (d2n)
        E      = 5840.0               ! Incident energy (MeV) 
        WMIN   = 3730.0
        WMAX   = 5210.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 [all materials before interaction]
        Ta     = 3.623E-02            ! Thickness in #X0 [all materials after interaction ] 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]   
      ELSEIF(TAG.EQ.'5_T3') THEN      ! SAMANTHA [3He] (d2n)
        E      = 5790.0               ! Incident energy (MeV) 
        WMIN   = 3730.0
        WMAX   = 5210.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 [all materials before interaction]
        Ta     = 3.623E-02            ! Thickness in #X0 [all materials after interaction ] 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]   
      ELSEIF(TAG.EQ.'5_T2') THEN      ! SAMANTHA [3He] (d2n)
        E      = 5740.0               ! Incident energy (MeV) 
        WMIN   = 3730.0
        WMAX   = 5210.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 [all materials before interaction]
        Ta     = 3.623E-02            ! Thickness in #X0 [all materials after interaction ] 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]   
      ELSEIF(TAG.EQ.'5_T1') THEN      ! SAMANTHA [3He] (d2n)
        E      = 5690.0               ! Incident energy (MeV) 
        WMIN   = 3630.0
        WMAX   = 5210.0
        DELTAW = 10.0
        Tb     = 2.928E-03            ! Thickness in #X0 [all materials before interaction]
        Ta     = 3.623E-02            ! Thickness in #X0 [all materials after interaction ] 
        density= 10.8                 ! (APPROX) Density [in Amg]
        zlength = 38.1                ! Target length [cm]   
      ELSEIF(TAG.EQ.'4GMA') THEN      ! GMA [Nitrogen] (d2n) 
        E      = 4730.0               ! Incident energy (MeV) 
        WMIN   = 2970.0
        WMAX   = 4150.0
        DELTAW = 10.0
        Tb     = 3.209E-03            ! Thickness in #X0 
        Ta     = 3.724E-02            ! Thickness in #X0 
        density= 7.71                 ! (APPROX) Density [in Amg]
        zlength = 39.65               ! Target length [cm]       
      ELSEIF(TAG.EQ.'5GMA') THEN      ! GMA [Nitrogen] (d2n) 
        E      = 5890.0               ! Incident energy (MeV) 
        WMIN   = 4130.0
        WMAX   = 5311.0
        DELTAW = 10.0
        Tb     = 3.209E-03            ! Thickness in #X0 
        Ta     = 3.724E-02            ! Thickness in #X0 
        density= 7.71                 ! (APPROX) Density [in Amg]
        zlength = 39.65               ! Target length [cm]          
C---------- XiaoChao Targets ----------------------------------
       ELSEIF(TAG.EQ.'GORE') THEN     ! XiaoChao's target
c        E      = 1230.0              ! Incident energy (MeV) 
c        WMIN   = 0.0
c        WMAX   = 1090.0
c        DELTAW = 10.0
        Tb     = 2.098E-03           ! Thickness in #X0 
        Ta     = 5.785E-02           ! Thickness in #X0 
        density= 9.10                ! (APPROX) Density [in Amg]
        zlength = 25.0               ! Target length [cm] 
        ELSEIF(TAG.EQ.'TIL') THEN    ! XiaoChao's target (2)
c        E      = 1230.0             ! Incident energy (MeV) 
c        WMIN   = 0.0
c        WMAX   = 1090.0
c        DELTAW = 10.0
        Tb     = 2.277E-03           ! Thickness in #X0 
        Ta     = 3.439E-02           ! Thickness in #X0 
        density= 8.28                ! (APPROX) Density [in Amg]
        zlength = 25.0               ! Target length [cm] 
C---------- E94-010 Targets ---------------------------------
      ELSEIF (TAG.EQ.'8_SY') THEN    ! Target is "Sysiphos" 
        E      = 862.0               ! Incident Energy in MeV
        WMIN   = 15.
        WMAX   = 600.
        DELTAW = 5.0
        Tb     = 7.184E-3            ! Thickness before scatt., including 1/2 target
        Ta     = 4.166E-2            ! Thickness after scatt., including targ contibution
        density= 9.939               ! in Amagats
        zlength = 40.0               ! in cm
      ELSEIF (TAG.EQ.'8hSY') THEN
        E      = 862.0
        WMIN   = 10.5
        WMAX   = 700.5
        DELTAW = 1.0
        Tb     = 7.184E-3
        Ta     = 4.273E-2
        density= 9.939
        zlength = 40.0
      ELSEIF(TAG.EQ.'1_DW') THEN ! Dont worry
        E      = 1717.9
        WMIN   = 45.
        WMAX   = 975.
        DELTAW = 5.0
        Tb     = 7.338E-3
        Ta     = 4.189E-2
        density= 12.25
        zlength = 40.0
      ELSEIF(TAG.EQ.'1hDW') THEN ! Dont worry
        E      = 1717.9
        WMIN   = 50.5
        WMAX   = 1500.5
        DELTAW = 1.0
        Tb     = 7.338E-3
        Ta     = 4.337E-2
        density= 12.25
        zlength = 40.0
      ELSEIF(TAG.EQ.'1_BH') THEN ! Be Happy
        E      = 1717.9
        WMIN   = 50.5
        WMAX   = 1500.5
        DELTAW = 5.0
        Tb     = 7.271E-3
        Ta     = 4.128E-2
        density= 11.53
        zlength = 40.0
      ELSEIF(TAG.EQ.'1hBH') THEN ! Be Happy
        E      = 1717.9
        WMIN   = 50.5
        WMAX   = 1500.5
        DELTAW = 1.0
        Tb     = 7.271E-3
        Ta     = 4.128E-2
        density= 11.53
        zlength = 40.0
      ELSEIF(TAG.EQ.'1REF') THEN ! Reference Cell
        E      = 1720.0  ! to compare with Pibero's XS
        WMIN   = 10.5
        WMAX   = 1500.5
        DELTAW = 5.0
        Tb     = 4.931E-3  ! Educated Guess 
        Ta     = 5.041E-2
        density=  0.0
        zlength = 0.0
      ELSEIF(TAG.EQ.'4REF') THEN ! ref cell.  
        E      = 4238.6
        WMIN   = 405.0
        WMAX   = 4015.0
        DELTAW = 10.0
        Tb     = 4.931E-3  ! Educated Guess
        Ta     = 5.041E-2
        density= 0.0
        zlength= 0.0
      ELSEIF(TAG.EQ.'1_AR') THEN ! Armegeddon
        E      = 1716.9
        WMIN   = 45.0
        WMAX   = 970.
C        DELTAW = 1.0
        DELTAW=5.
        Tb     = 7.922E-3
        Ta     = 7.492E-2
        density= 12.32
        zlength = 40.0
      ELSEIF(TAG.EQ.'1hAR') THEN ! Armegeddon
        E      = 1716.9
        WMIN   = 10.0
        WMAX   = 1500.5
        DELTAW = 1.0
        Tb     = 7.922E-3
C       Ta     = 8.505E-2
        Ta     = 9.92E-2    ! max possible win thickness.
        density= 12.32
        zlength = 40.0
      ELSEIF(TAG.EQ.'2_JN') THEN ! JIN
        E      = 2580.5
        WMIN   = 85.0
        WMAX   = 1905.
        DELTAW = 10.0
        Tb     = 7.269E-3
        Ta     = 3.739E-2
        density= 10.23
        zlength = 40.0
      ELSEIF(TAG.EQ.'2hJN') THEN ! JIN
        E      = 2580.5
        WMIN   = 50.0
        WMAX   = 2395.0
        DELTAW = 10.0
        Tb     = 7.269E-3
        Ta     = 3.902E-2
        density= 10.23
        zlength = 40.0
      ELSEIF(TAG.EQ.'3_AR') THEN ! Armageddon
        E      = 3381.8
        WMIN   = 240.0
        WMAX   = 2440.0
        DELTAW = 20.0
        Tb     = 7.922E-3
        Ta     = 7.492E-2
        density= 12.32
        zlength = 40.0
      ELSEIF(TAG.EQ.'3min') THEN ! Armageddon
        E      = 3381.8
        WMIN   = 255.0
        WMAX   = 3135.0
        DELTAW = 10.0
        Tb     = 7.922E-3
        Ta     = 6.968E-2
        density= 12.32
        zlength = 40.0
      ELSEIF(TAG.EQ.'3hAR') THEN ! Armageddon
        E      = 3381.8
        WMIN   = 255.0
        WMAX   = 3135.0
        DELTAW = 10.0
        Tb     = 7.922E-3
C        Ta     = 8.505 ! alex's best estimate 
        Ta     = 9.92E-2    ! max possible win thickness. (NOt true)
        density= 12.32
        zlength = 40.0
      ELSEIF(TAG.EQ.'4_NE') THEN ! Nephali
        E      = 4238.6
        WMIN   = 465.0
        WMAX   = 2925.0
        DELTAW = 30.0
        Tb     = 7.381E-3
        Ta     = 4.116E-2
        density= 13.77
        zlength = 40.0
      ELSEIF(TAG.EQ.'4hNE') THEN ! Nephali
        E      = 4238.6
        WMIN   = 405.0
        WMAX   = 4015.0
        DELTAW = 10.0
        Tb     = 7.381E-3
        Ta     = 4.116E-2
        density= 13.77
        zlength = 40.0
      ELSEIF(TAG.EQ.'4_SY') THEN ! Sysiphos
        E      = 4238.6
        WMIN   = 465.0
        WMAX   = 2925.0
        DELTAW = 30.0
        Tb     = 7.184E-3
        Ta     = 4.166E-2
        density= 9.939
        zlength = 40.0
      ELSEIF(TAG.EQ.'4hSY') THEN ! Sysiphos
        E      = 4238.6
        WMIN   = 405.0
        WMAX   = 4015.0
        DELTAW = 10.0
        Tb     = 7.184E-3
        Ta     = 4.273E-2
        density= 9.939
        zlength = 40.0
      ELSEIF(TAG.EQ.'5_JN') THEN ! JIN
        E      = 5058.2
        WMIN   = 400.0
        WMAX   = 4000.0
        DELTAW = 50.0
        Tb     = 7.269E-3
        Ta     = 3.739E-2
        density= 10.23
        zlength = 40.0
      ELSEIF(TAG.EQ.'5hJN') THEN ! JIN
        E      = 5058.2
        WMIN   = 405.0
        WMAX   = 4015.0
        DELTAW = 10.0
        Tb     = 7.269E-3
        Ta     = 3.902E-2
        density= 10.23
        zlength = 40.0
C---------- Miscellaneous  ----------------------------------
c     FIXME: GMA Cell -- Tb, Ta and density are NOT correct! 
      ELSEIF(TAG.EQ.'GMA') THEN ! Reference Cell (d2n -- H2 cell)
        E      = 1230.0  
        WMIN   = 0.
        WMAX   = 1300.0
        DELTAW = 5.0
        Tb     = 8.43E-3
        Ta     = 0.0832
        density= 9.61 * 2.0          
        zlength = 39.4               ! measured
c   END d2n Targets
C-- REFERENCE CELL
      ELSEIF(TAG.EQ.'1REF') THEN ! Reference Cell (Duke's ref. cell)
        E      = 1046.0  
        WMIN   = 0.
        WMAX   = 1050.0
        DELTAW = 5.0
        Tb     = 8.43E-3
        Ta     = 0.0832
        density= 9.61 * 2.0  ! Pavg = 142 psig, cell leaked
        zlength = 39.2           ! measured
      ELSEIF(TAG.EQ.'3REF') THEN ! Reference Cell (Exodus's ref. cell)
        E      = 3028.1
        WMIN   = 600.
        WMAX   = 2250.
        DELTAW = 20.0
        Tb = 8.76E-3
        Ta = 0.0606
        density= 9.94 * 2.0  ! P = 147 psig
        zlength= 39.2            ! approx. 
      ELSEIF(TAG.EQ.'4REF') THEN ! Reference Cell (Duke's ref. cell) 
        E      = 4017.9
        WMIN   = 1060.
        WMAX   = 2500.
        DELTAW = 20.0
        Tb     = 8.36E-3
        Ta     = 0.0569
        density= 9.50 * 2.0  ! Pavg = 139 psig, cell leaked
        zlength=  39.2           ! measured
      ELSEIF(TAG.EQ.'5REF') THEN ! Reference Cell (Duke's ref. cell) 
        E      = 5009.0
        WMIN   = 1580.
        WMAX   = 2930.
        DELTAW = 20.0
        Tb     = 8.34E-3 
        Ta     = 0.0595
        density= 9.48 * 2.0  ! Pavg = 140.5 psig, cell leaked
        zlength= 39.2            ! measured
      ELSEIF(TAG.EQ.'6REF') THEN ! Reference Cell (Duke's ref. cell)
        E      = 5009.0
        WMIN   = 2170.
        WMAX   = 3280.
        DELTAW = 20.0
        Tb     = 8.40E-3
        Ta     = 0.0469
        density= 9.57 * 2.0  ! Pavg = 141.3 psig, cell leaked
        zlength= 39.2            ! measured
      ELSEIF(TAG.EQ.'7REF') THEN ! Reference Cell (Exodus's ref. cell)
        E      = 5009.0
        WMIN   = 2170.
        WMAX   = 3280.
        DELTAW = 20.0
        Tb     = 8.75E-3
        Ta     = 0.0498
        density= 9.92 * 2.0  ! P = 147 psig
        zlength= 39.2            ! approx.
      ELSEIF(TAG.EQ.'2_JP') THEN ! JP's thesis
        E      = 2699.7
        Tb     = 0.05910/2.+0.000285
        Ta     = 0.05910/2.+(0.00346+0.00259+0.00124)
        density= 7.83           ! g/cm**3
        zlength = 0.02908       ! inches
      ELSEIF(INTERACTIVE) THEN ! will prompt user for input
        density = 0.0           ! g/cm**3
        zlength = 0.0           ! inches
        CONTINUE
      ELSE
        WRITE(6,*) TAG,':  Unknown Energy'
        STOP
      ENDIF
C      density = 7.4
      density=density*2.6868E19 ! Convert from Amagat to (cm)^-3
C      WRITE(8,'(A,F6.1,x,2(F9.5,x),E10.5,x,F5.1)') 
C     &         '#',E,    Tb,Ta,    density,zlength 

C     For Huan's code! These values have only the target cell + 1/2 target material  
C      Tb = 2.165E-3
C      Ta = 3.393E-2
C     For systematic error studies (d2n ONLY) 
C     Tb = Tb*(1 + 0.1) increase by 10% 
C      Tb = 3.221E-3 
C     Tb = Tb*(1-0.1) decrease by 10% 
C      Tb = 2.899E-3 
C     Ta = Ta*(1+0.1) increase by 10% 
C      Ta = 3.985E-2 
C     Ta = Ta*(1-0.1) decrease by 10% 
C      Ta = 3.261E-2 

      RETURN 
      END
C-------------------------------------------------------------------------- 
      SUBROUTINE PROMPT(Z,A,E,TH,PF,EPS,EPSD,Tb,Ta)
      IMPLICIT REAL*8 (A-H,O-Z)
C     Prompt user for input
      WRITE(6,'(A)')  'Enter Z,A: '
      READ(*,*) Z,A
c      Z = 7.0
c      A = 14.0
      WRITE(6,'(A)' ) 'Incident energy(MeV) and scatt angle(deg.): '
      READ(*,*) E,TH
c      WRITE(6,'(A)' ) 'P_fermi,Nucleon and Delta separation (MeV): '
c      READ(*,*) PF,EPS,EPSD
      PF   = 221.0
      EPS  = 25.0
      EPSD = 15.0
      WRITE(6,'(A)') 'Target thickness Tb, Ta (in radiation lengths):'
      READ(*,*) Tb,Ta
c      Tb = 0.0
c      Ta = 0.0      
      RETURN
      END
C-------------------------------------------------------------------
      SUBROUTINE GETRATES(deltaP,current,domega)
      IMPLICIT REAL*8 (A-H,O-Z)
C-----INPUTS NEEDED FOR RATES CALCULATIONS...
C   --E94010 RUNPLAN VALUES--
      deltaP  = .04
      current = 10.0
      domega  = 3.9
c      density = 7.4
C   --MORE REALISTIC VALUES--
C      deltaP = .045           ! Spectrometer Momentum acceptance. DeltaP/P
C      current= 10.0           ! Amps
C      domega = 6.7            ! Acceptance in msr

      current= current*1.0E-6 ! Convert from microAmps to Amps
      domega = domega*1.0E-3  ! Convert from msr to sr
      RETURN
      END
C-------------------------------------------------------------------
      SUBROUTINE INITIALIZE(Pmax,deltaP,xnu,sigmatot,sigave,P0,P0N)
      IMPLICIT REAL*8 (A-H,O-Z)
      real*8 P0(30),P0N(30)
      real*8 sigmatot(1000), xnu(1000),sigave(30)
C----Initialize some arrays used for rates calculations
      do i =1,1000
        sigmatot(i)=0
        xnu(i)     =0
      enddo
      do n=1,30
        P0(n)    =Pmax*(1.0-deltap)**(2*n-1) ! find central momentums
        P0N(n)   =0                          ! Number of bins covered
        sigave(n)=0.0
      enddo
C------
      RETURN
      END
C-------------------------------------------------------------------
CDEBUG      These are the most recent form factors available. 
C           Do not use option 'Qdep' with these FF.
CDEBUG      REAL*8 FUNCTION GEP(QMS,AP)
CDEBUG      IMPLICIT REAL*8 (A-H,O-Z)
CDEBUG      GEP=1./(1.+QMS/AP**2)**2
CDEBUG      RETURN
CDEBUG      END 
CDEBUG      
CDEBUG      REAL*8 FUNCTION GEN(QMS,AP)
CDEBUG      IMPLICIT REAL*8 (A-H,O-Z)
CDEBUGC     H.ZHU et al. PRL 87, Number 8, 081801-1 (2001)
CDEBUG
CDEBUG      PM   = 939.
CDEBUG      UN   = -1.91304184
CDEBUG      xp   = 4.36   ! +- 1.11
CDEBUG      xa   = 0.895  ! +- 0.039 
CDEBUG
CDEBUG      GEN = -UN * xa * TAU(QMS)/( 1.0+xp*TAU(QMS) )
CDEBUG      GEN = GEN * GEP(QMS,AP)
CDEBUG      RETURN
CDEBUG      END
CDEBUG
CDEBUG      REAL*8 FUNCTION GMP(QMS,AP)
CDEBUG      IMPLICIT REAL*8 (A-H,O-Z)
CDEBUGC     Gayou et al. PRL 88,number 9, 092301-1 (2002)
CDEBUG      UP   =  2.7928456
CDEBUG      GMP  =  UP * GEP(QMS,AP)
CDEBUG      GMP  =  GMP/(1.0 - 0.13*(QMS-0.04) )
CDEBUG      RETURN
CDEBUG      END
CDEBUG
CDEBUG      REAL*8 FUNCTION GMN(QMS,AP)
CDEBUG      IMPLICIT REAL*8 (A-H,O-Z)
CDEBUG      UN   = -1.91304184
CDEBUG      GMN  = UN * GEP(QMS,AP)
CDEBUG      RETURN
CDEBUG      END
CDEBUG
C-------------------------------------------------------------------
      REAL*8 FUNCTION TAU(QMS)
      IMPLICIT REAL*8 (A-H,O-Z)
      PM   = 939.
      TAU  = QMS/4.0/PM**2
      RETURN
      END
C-------------------------------------------------------------------
      REAL*8 FUNCTION GEP(QMS,AP)
      IMPLICIT REAL*8 (A-H,O-Z)
      GEP=1./(1.+QMS/AP**2)**2
      RETURN
      END
C-------------------------------------------------------------------
      REAL*8 FUNCTION GEN(QMS,AP)
      IMPLICIT REAL*8 (A-H,O-Z)
      PM   = 939.
      UN   = -1.91304184

      GEN = -UN
      GEN = GEN * TAU(QMS)/( 1.0+5.6*TAU(QMS) )
      GEN = GEN * GEP(QMS,AP)
      RETURN
      END
C-------------------------------------------------------------------
      REAL*8 FUNCTION GMP(QMS,AP)
      IMPLICIT REAL*8 (A-H,O-Z)
      UP   =  2.7928456
      GMP  =  UP * GEP(QMS,AP)
      RETURN
      END
C-------------------------------------------------------------------
      REAL*8 FUNCTION GMN(QMS,AP)
      IMPLICIT REAL*8 (A-H,O-Z)
      UN   = -1.91304184
      GMN  = UN * GEP(QMS,AP)
      RETURN
      END
C-------------------------------------------------------------------
      SUBROUTINE SIM(FL,FU,DEL,N,F,ANS)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER N,NX,II

C      EXTERNAL F       
C      WRITE(*,*) 'FL = ',FL
C      WRITE(*,*) 'FU = ',FU
C      WRITE(*,*) 'ANS = ',ANS 
C      WRITE(*,*) 'DEL = ',DEL
C      WRITE(*,*) 'N = ',N

C      FL  = Lower bound
C      FU  = Upper bound 
C      DEL = Step size
C      N   = Minimum number of steps
C      F   = Function 
C      ANS = Placeholder to store the answer 

      NX = ((FU - FL)/DEL/2.)
      NX = 2*NX

      IF(NX.LT.N) NX=(N/2)*2

      DT=(FU-FL)/FLOAT(NX)
      SI = F(FL)
      SI = SI + F(FU)

      XX = FL + DT
      SI = SI + F(XX) *4.

      AA=FL
      II=0
5     II=II+2

      IF (II.GE.NX)  GO TO 7

      AA =FLOAT (II)*DT
      AA=AA+FL
      BB=AA+DT

      SI = SI + 2.*F(AA) + 4.*F(BB)
      GO TO 5

7      ANS = DT * SI/3.0
C 8      FORMAT (' FL,FU,DEL,DT,N,NX=', 6(G11.4,1X), ' ANS ', G11.4)

      RETURN
      END
C----------------------------------------------------------------
*FES 
      REAL*8 FUNCTION FES(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL extrad
      COMMON/PAR/E,TH,W,Z,A,EPS,EPSD,PF,SPENCE
      COMMON/ADDONS/SCALE,Tb,Ta,extrad
      REAL MYPHI

C      REAL X,VAR,EF ! use this to convert X (which is really W')  

      EF  = E-W
      VAR = E-W+X ! = Es-(Es-Ep)+(Es'-Ep) = Es' 

      IF (.NOT.extrad) THEN  ! apply external rad. cor.
        Tb=0.0
        Ta=0.0
      ENDIF

      ALPH  = 1./137.03604
      EMASS = 0.511
      PI    = ACOS(-1.)
      THR   = TH*PI/180.
      xb    = 4./3.
      XM    = 931.49 ! Mass of the nucleon
      XMT   = A*XM   ! Mass of the target

      eta = LOG(1440.*Z**(-2./3.) )/LOG(183.*Z**(-1./3.) )
      xi  = (PI*EMASS/2./ALPH)*(Ta+Tb)
      xi  = xi/( (Z+eta)*LOG(183.*Z**(-1./3.)) )
      R   = ( XMT+E*(1-COS(THR)) )/( XMT-(EF)*(1-COS(THR)) )
C
C-----Stein's 1st integral dEs'
C
      QMS1  = 4.*(VAR)*(EF)*SIN(THR/2.)**2    !    1/15/03
      tr1   = 1./xb*(ALPH/PI)*(LOG(QMS1/EMASS**2)-1.)
      Tpb   = tr1 + Tb
      Tpa   = tr1 + Ta

      D2    = 13.*(LOG(QMS1/EMASS**2)-1.)/12.-17./36.
      D2    = D2 - 1./4.*( LOG( E/(EF) ) )**2 !Corr. to peak. approx.
      D2    = D2 + 0.5*(PI**2/6.-SPENCE)
      D2    = D2 * (2.*ALPH/PI)
      D2    = D2 + 0.5772*xb*(Tb+Ta) ! 1/14/02

C      SIG1  = SIGQFS(E-W+X,TH,X,Z,A,EPS,PF)
C      SIG1  = SIG1 + SIGDEL(E-W+X,TH,X,A,EPSD,PF)
C      SIG1  = SIG1 +   SIGX(E-W+X,TH,X,A)
C      SIG1  = SIG1 +  SIGR1(E-W+X,TH,X,A,PF)
C      SIG1  = SIG1 +  SIGR2(E-W+X,TH,X,A,PF)
C      SIG1  = SIG1 +  SIG2N(E-W+X,TH,X,Z,A,PF)

C      F1   = ( xb*Tpb/(W-X) ) *phi((W-X)/(E))   ! 
C      F1   = F1 + xi/(2.*(W-X)**2)
C      F1   = F1 * SIG1*(1.+D2)
C      F1   = F1 * ( (W-X)/(EF*R) )**(xb*Tpa)
C      F1   = F1 * ( (W-X)/ (E)      )**(xb*Tpb) ! 

C     NOTE: SIG = SIG(Es',Ep)... So, using VAR and X is 
C     indeed correct here -- VAR = Es', X = Es'-Ep.
C     => SIG = SIG(Es',Ep). 
      SIG1=0 
      SIG1  = SIGQFS(VAR,TH,X,Z,A,EPS,PF)
      SIG1  = SIG1 + SIGDEL(VAR,TH,X,A,EPSD,PF)
      SIG1  = SIG1 +  SIG2N(VAR,TH,X,Z,A,PF)
      SIG1  = SIG1 +  SIGR1(VAR,TH,X,A,PF)
      SIG1  = SIG1 +  SIGR2(VAR,TH,X,A,PF)
      SIG1  = SIG1 +   SIGX(VAR,TH,X,A)

      F1   = ( xb*Tpb/(E-VAR) ) *phi((E-VAR)/(E))   ! 
      F1   = F1 + xi/(2.*(E-VAR)**2)
      F1   = F1 * SIG1*(1.+D2)
      F1   = F1 * ( (E-VAR)/(EF*R) )**(xb*Tpa)
      F1   = F1 * ( (E-VAR)/(E)    )**(xb*Tpb) ! 

      FES = F1

      MYPHI = phi( (E-VAR)/E )


C      IF(EF.EQ.600 .AND. VAR.EQ.3000.0) THEN
C        WRITE(*,*) SIG1*SCALE 
C      ENDIF 

C     ORDER: Es', b, Tr+Tb, Tr+Ta, phi, xi, sig  

C      IF(EF.GT.1858.0 .AND. EF.LT.1859.0 ) THEN
C        WRITE(*,*) VAR,  '  ',xb,'  '
C     &          ,Tpb,'  ',Tpa,'  ',MYPHI,'  ',xi,'  ',SIG1*SCALE
C      ENDIF


      RETURN
      END
C--------------------------------------------------------------
*FEP
      REAL*8 FUNCTION FEP(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL extrad
      COMMON/PAR/E,TH,W,Z,A,EPS,EPSD,PF,SPENCE
      COMMON/ADDONS/SCALE,Tb,Ta,extrad

      REAL*8 VAR,EF
      REAL MYPHI

      EF  = E-W
      VAR = E-X ! = Es-(Es-Ep') = Ep'   

      if (.NOT.extrad) THEN  ! apply external rad. cor.
        Tb=0.0
        Ta=0.0
      endif
      ALPH  = 1./137.03604
      EMASS = 0.511
      PI    = ACOS(-1.)
      THR   = TH*PI/180.
      xb    = 4./3.
      XM    = 931.49 ! Mass of the nucleon
      XMT   = A*XM   ! Mass of the target

      eta = LOG(1440.*Z**(-2./3.) )/LOG(183.*Z**(-1./3.) )
      xi  = (PI*EMASS/2./ALPH)*(Ta+Tb)
      xi  = xi/( (Z+eta)*LOG(183.*Z**(-1./3.)) )
      R   = ( XMT+E*(1-COS(THR)) )/( XMT-(EF)*(1-COS(THR)) )
C
C-----Stein's 2nd integral dEp'
C
      QMS2  = 4.*E*(VAR)*SIN(THR/2.)**2  ! 1/15/03
      tr2   = 1./xb*(ALPH/PI)*(LOG(QMS2/EMASS**2)-1.)
      Tpb   = tr2 + Tb
      Tpa   = tr2 + Ta

      D2    = 13.*(LOG(QMS2/EMASS**2)-1.)/12.-17./36.
      D2    = D2 - 1./4.*( LOG( E/(EF) ) )**2 !KS. Correction to peak. approx.
      D2    = D2 + 0.5*(PI**2/6.-SPENCE)
      D2    = D2 * (2.*ALPH/PI)
      D2    = D2 + 0.5772*xb*(Tb+Ta)

C      SIG2  = SIGQFS(E,TH,X,Z,A,EPS,PF)
C      SIG2  = SIG2 + SIGDEL(E,TH,X,A,EPSD,PF)
C      SIG2  = SIG2 + SIGX(E,TH,X,A)
C      SIG2  = SIG2 + SIGR1(E,TH,X,A,PF)
C      SIG2  = SIG2 + SIGR2(E,TH,X,A,PF)
C      SIG2  = SIG2 + SIG2N(E,TH,X,Z,A,PF)
C
C      F2   = ( xb*Tpa/(W-X) ) *phi((W-X)/(E-X))
C      F2   = F2 + xi/(2.*(W-X)**2)
C      F2   = F2 * SIG2*(1.+D2)
C      F2   = F2 * ( (W-X)/(E-X) )**(xb*Tpa)
C      F2   = F2 * ( (W-X)*R/(E) )**(xb*Tpb)

C     NOTE: SIG = SIG(Es,Ep').  Using E,X is correct 
C     here, since X = Es-Ep'.
C     => SIG = SIG(Es,Ep'). 
      SIG2=0
      SIG2  = SIGQFS(E,TH,X,Z,A,EPS,PF)
      SIG2  = SIG2 + SIGDEL(E,TH,X,A,EPSD,PF)
      SIG2  = SIG2 + SIG2N(E,TH,X,Z,A,PF)
      SIG2  = SIG2 + SIGR1(E,TH,X,A,PF)
      SIG2  = SIG2 + SIGR2(E,TH,X,A,PF)
      SIG2  = SIG2 + SIGX(E,TH,X,A)

      F2   = ( xb*Tpa/(VAR-EF) ) *phi((VAR-EF)/(VAR))
      F2   = F2 + xi/(2.*(VAR-EF)**2)
      F2   = F2 * SIG2*(1.+D2)
      F2   = F2 * ( (VAR-EF)/(VAR) )**(xb*Tpa)
      F2   = F2 * ( (VAR-EF)*R/(E) )**(xb*Tpb)

      FEP = F2

      MYPHI = phi( (VAR-EF)/VAR )

C      IF(E.EQ.5890.0 .AND. EF.EQ.600.0 ) THEN
C        WRITE(*,*) VAR,  '  ',xb,'  '
C     &          ,Tpb,'  ',Tpa,'  ',MYPHI,'  ',xi,'  ',SIG2*SCALE
C      ENDIF


      RETURN
      END
C--------------------------------------------------------------


