C ********    LES37ST2.F   ********************************************
C *  Large Eddy Simulation of                                         *
C *                       Turbulent Flow in a Two-Dimensional Channel *
C * ----------------------------------------------------------------- *
C *     SGS Model     ....  Smagorinsky Eddy Viscosity Model          *
C *     Spacial FDM   ....  2th order, Central, Staggered mesh        *
C *     Time Marching ....  Nonlinear & Viscous : 2nd Adams-Bashforth *
C *                         Pressure & Continuity: Backward-Euler     *
C * ----------------------------------------------------------------- *
C *   KAJISHIMA, T.                                                   *
C *      Dept. Mechanical Engineering,   Osaka University             *
C *              Yamadaoka, Suita, Osaka, 565 Japan                   *
C *              Phone: +81 6 879 7249   Fax: +81 6 879 7250          *
C *              E-mail: kajisima@mech.eng.osaka-u.ac.jp              *
C *********************************************************************

      PROGRAM LES37ST2
      PARAMETER (NX=32,NY=64,NZ=32)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*6  FILE1,FILE2
      COMMON /FILE1/FILE1 /FILE2/FILE2 
      COMMON /RMS/URMSX,URMSC,VRMSC,WRMSC
      COMMON /RET/RET /DT/DT /STEP/ISTEP /TSTEP/TSTEP /START/ISTART
      COMMON /DIV/DIVX /COU/COMX /ENE/ENE /UME/UME,UMX /CSGS/CSGS
      COMMON /ERRP/POIERR /ITRP/ITRP /IRMS/IRMS 
         FILE1='cs3701'
         FILE2='cs3702'
         DT=  2.5D-4
         RET=360.D0
         CSGS=0.1D0
         ICOUNT=0
         IRMS=0
      CALL SBRMSH
      CALL SBRDTR
C      WRITE(6,*) 'How many time step counts to advance ?'
C      READ(5,*) ISTOP
      ISTOP=1000

      WRITE( 6,1000)
        CALL SBRUMR
        CALL SBRST1
        CALL SBRCHK
        WRITE( 6,1100) ISTART,TSTEP,ITRP,POIERR,DIVX,COMX
     &   ,UME,UMX,INT(RET*UME),INT(RET*UMX),ENE,URMSX,URMSC,VRMSC,WRMSC
 1000 FORMAT(3X,4HSTEP,7X,1HT,1X,4HITRP,4X,6HPOIerr,4X,6HDIVmax
     &,2X,6HCOUmax,3X,5HUmean,4X,4HUmax,4X,3HRem,4X,3HRex
     &,4X,6HEnergy,3X,5HU'max,5X,3HU'c,5X,3HV'c,5X,3HW'c,2X,6HCPU(s))
 1100 FORMAT(I7,F8.4,I5,2E10.2,3F8.3,2I7,F10.5,5F8.4)

   10 CONTINUE
C +++ CPU TIME : CALL CLOCK(DT0)
      CALL SBRCON
      CALL SBRVGT
      CALL SBRSGS
      CALL SBRFLX
      CALL SBRVIS
      IF(ICOUNT.GE.1.AND.IRMS.EQ.0) THEN
        CALL SBRPRE(1.5D0,-0.5D0)
        ELSE
        CALL SBRPRE(1.0D0, 0.0D0)
        ENDIF
      CALL SBRRHP
      CALL SBRSOR(20)
      CALL SBRCOR

      ICOUNT=ICOUNT+1
      ISTEP=ISTART+ICOUNT
      TSTEP=TSTEP+DT
C +++ CPU TIME : CALL CLOCK(DT1)

      IF(MOD(ISTEP,50).eq.0.OR.ICOUNT.LE.50) THEN
        CALL SBRUMR
        CALL SBRST1
        CALL SBRCHK
        WRITE( 6,1100) ISTEP,TSTEP,ITRP,POIERR,DIVX,COMX
     &   ,UME,UMX,INT(RET*UME),INT(RET*UMX),ENE,URMSX,URMSC,VRMSC,WRMSC
C +++ &   ,DT1-DT0
      ENDIF

      IF(MOD(ISTEP,5000).EQ.0) CALL SBRDTS
      IF(ISTEP.GE.ISTOP) GO TO 20

      IF(ISTEP.GE.50.AND.DIVX.GE.1.D+2) GO TO 30
      GO TO 10

   20 CONTINUE
      CALL SBRDTS
      STOP

   30 CONTINUE
      WRITE(6,*) 'Stopped  because of the numerical instability'
      STOP
      END

C *********************************************************************
C *     SBR. MSH  :  MESH PARAMETERS                                  *
C *********************************************************************

      SUBROUTINE SBRMSH
      PARAMETER (NX=32,NY=64,NZ=32,NY1=NY+1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /H/HX,HZ /Y/YV(0:NY),YP(0:NY1)
      COMMON /D/DX,DZ /DY/DY(NY) /RET/RET /DT/DT
      CALL MESHXZ
      CALL MESHY

      WRITE(6,1000) NX,HX,DX,RET*HX,RET*DX,NZ,HZ,DZ,RET*HZ,RET*DZ
     &,NY,1.,DY(1),DY(NY/2),RET*DY(1),RET*DY(NY/2),YP(1),RET*YP(1)
     &,INT(RET),DT
 1000 FORMAT(/10(1H*),"  IN SBR.MSH  ",10(1H*)/
     &/5X,"NX=",I4,5X,"HX=",F7.3,3X,"DX=",F7.5
     &,5X,"HX+=",F7.1,3X,"DX+=",F7.3
     &/5X,"NZ=",I4,5X,"HZ=",F7.3,3X,"DZ=",F7.5
     &,5X,"HZ+=",F7.1,3X,"DZ+=",F7.3
     &/5X,"NY=",I4,5X,"HY=",F7.3,3X,"DYmin=",F7.5,3X,"DYmax=",F7.5
     &/30X,"DYmin+=",F6.3,3X,"DYmax+=",F6.2
     &/5X,"YP(1)=",F7.5,3X,"YP(1)+=",F6.2
     &/5X,"RET=",I6/5X,"DT=",F9.5)

      CALL MESH1V
      CALL MESH1P
      CALL MESH2V
      CALL MESH2P
      CALL MESHIN
      RETURN
      END
C---------------------------------------------------------------------
C-----  Mesh generation .....  for X, Z directions  ------------------
      SUBROUTINE MESHXZ
      PARAMETER (NX=32,NY=64,NZ=32)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /H/HX,HZ /D/DX,DZ /RET/RET
      COMMON /IP/IP(NX) /IM/IM(NX) /KP/KP(NZ) /KM/KM(NZ)
      COMMON /TXZA/TXA,TZA /BXZA/BXA,BZA 
      COMMON /CMPX/CM1X,C00X,CP1X /CMPZ/CM1Z,C00Z,CP1Z
      DX=36.D0/RET
      DZ=18.D0/RET
      HX=NX*DX
      HZ=NZ*DZ
      BXA=1.D0/DX
      BZA=1.D0/DZ
      TXA=BXA/2.D0
      TZA=BZA/2.D0

C  ---  List vectors to identify the neighbouring point
      DO 50 I=1,NX
      IP(I)=I+1
      IM(I)=I-1
        IF(IP(I).GT.NX) IP(I)=IP(I)-NX
        IF(IM(I).LT. 1) IM(I)=IM(I)+NX
   50 CONTINUE
      DO 60 K=1,NZ
      KP(K)=K+1
      KM(K)=K-1
        IF(KP(K).GT.NZ) KP(K)=KP(K)-NZ
        IF(KM(K).LT. 1) KM(K)=KM(K)+NZ
   60 CONTINUE

C  ---  Parameters for finite-difference of Poisson equation
      CM1X= 1.D0/DX**2
      C00X=-2.D0/DX**2
      CP1X= 1.D0/DX**2
      CM1Z= 1.D0/DZ**2
      C00Z=-2.D0/DZ**2
      CP1Z= 1.D0/DZ**2
      RETURN
      END
C---------------------------------------------------------------------
C-----  Mesh generation .....  for Y direction  ----------------------
      SUBROUTINE MESHY
      PARAMETER (NX=32,NY=64,NZ=32,NY1=NY+1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /Y/YV(0:NY),YP(0:NY1) /D/DX,DZ /DY/DY(NY) /DYC/DYC(NY-1)
      COMMON /YPLS/YPLSV(0:NY),YPLSP(NY) /CAK/CAK(NY)
      COMMON /RET/RET /CSGS/CSGS
      ALG=0.95D0
      AT=DLOG((1.D0+ALG)/(1.D0-ALG))/2.D0
      YV(0)=0.
      DO 10 J=1,NY-1
      ETA=AT*(-1.D0+2.D0*DBLE(J)/DBLE(NY))
      YV(J)=(DTANH(ETA)/ALG+1.D0)/2.D0
   10 CONTINUE
      YV(NY)=1.D0
      DO 15 J=1,NY
      ETA=AT*(-1.D0+2.D0*(DBLE(J)-0.5D0)/DBLE(NY))
      YP(J)=(DTANH(ETA)/ALG+1.D0)/2.D0
   15 CONTINUE
C ... Outer points (half mesh)
        YP(0)=2.D0*YV(0)-YP(1)
        YP(NY1)=2.D0*YV(NY)-YP(NY)

      DO 20 J=1,NY
        DY(J)=-YV(J-1)+YV(J)
C     for Smagorinsky model
        YPLSP(J)=RET*DMIN1(YP(J),(1.D0-YP(J)))
        FS=1.D0-DEXP(-YPLSP(J)/25.D0)
        DS=(DX*DY(J)*DZ)**(1.D0/3.D0)
        CAK(J)=(CSGS*DS*FS)**2
   20 CONTINUE
      DO 30 J=1,NY-1
        DYC(J)=-YP(J)+YP(J+1)
   30 CONTINUE
      DO 40 J=0,NY
        YPLSV(J)=RET*DMIN1(YV(J),(1.D0-YV(J)))
   40 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
C---- D1VNM,D1VNP:  FDM operators at V points using data 
C                   on P points and INSIDE walls 
C                   with Neumann B.C.
C---- D1VM,D1VP:    FDM operators at V points using data 
C                   on P points and AT walls 
C---- D0VM,D0VP:    Interpolation operators at V points using data on P points

      SUBROUTINE MESH1V
      PARAMETER (NY=64,NY1=NY+1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /Y/YV(0:NY),YP(0:NY1) /DYC/DYC(NY-1)
      COMMON /D1VN/D1VNM(0:NY),D1VNP(0:NY) 
      COMMON /D0V/D0VM(NY-1),D0VP(NY-1) /D1V/D1VM(0:NY),D1VP(0:NY)
      DO 10 J=1,NY-1
      D1VNM(J)=-1.D0/DYC(J)
      D1VNP(J)= 1.D0/DYC(J)

      D0VM(J)=(YP(J+1)-YV(J))/DYC(J)
      D0VP(J)=(YV(J)-YP(J))/DYC(J)
      D1VM(J)=-1.D0/DYC(J)
      D1VP(J)= 1.D0/DYC(J)
   10 CONTINUE

C     For the Neumann B.C. at the wall (DP/DY=0) for Pressure
      D1VNM(0)= 0.D0
      D1VNP(0)= 0.D0
      D1VNM(NY)= 0.D0
      D1VNP(NY)= 0.D0

Caution!  Stencils are shifted at the wall point
C         (3-points FDM to keep the order of accuracy)
C         using  Uw (=0) at y=0,  U(1) at YP(1),  U(2) at YP(2)
C         DU/DY= -(YP(1)+YP(2))/YP(1)/YP(2)*Uw
C                +D1VM(0)*U(1)+D1VP(0)*U(2)
      DYP0=YP(2)-YP(1)
      D1VM(0)= YP(2)/YP(1)/DYP0
      D1VP(0)=-YP(1)/YP(2)/DYP0
      DYPN=YP(NY)-YP(NY-1)
      D1VM(NY)= (1.-YP(NY))/(1.-YP(NY-1))/DYPN
      D1VP(NY)=-(1.-YP(NY-1))/(1.-YP(NY))/DYPN
      RETURN
      END

C-----------------------------------------------------------------------
C---  FDM operators for velocity at P points using data on V points ----
      SUBROUTINE MESH1P
      PARAMETER (NY=64,NY1=NY+1,NYH1=NY-1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /Y/YV(0:NY),YP(0:NY1) /DY/DY(NY)
      COMMON /D0P/D0PM(NY),D0PP(NY) /D1P/D1PM(NY),D1PP(NY)
      DO 10 J=1,NY
      D0PM(J)=(YV(J)-YP(J))/DY(J)
      D0PP(J)=(YP(J)-YV(J-1))/DY(J)
      D1PM(J)=-1.D0/DY(J)
      D1PP(J)= 1.D0/DY(J)
   10 CONTINUE
      RETURN
      END
C---  FDM operators for diffusion terms   ------------------------------
      SUBROUTINE MESH2V
      PARAMETER (NY=64,NY1=NY+1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /D2V/D2VM(NY-1),D2V0(NY-1),D2VP(NY-1) 
      COMMON /D2P/D2PM(2:NY),D2P0(NY),D2PP(NY-1) 
      COMMON /D1V/D1VM(0:NY),D1VP(0:NY) /D1P/D1PM(NY),D1PP(NY)
C - at V points  -------------------
      DO 10 J=1,NY-1
      D2VM(J)=D1VM(J)*D1PM(J)
      D2V0(J)=D1VM(J)*D1PP(J)+D1VP(J)*D1PM(J+1)
      D2VP(J)=D1VP(J)*D1PP(J+1)
   10 CONTINUE
C - at P points  -------------------
      J=1
      D2P0(J)=D1PM(J)*D1VM(J-1)+D1PP(J)*D1VM(J)
      D2PP(J)=D1PM(J)*D1VP(J-1)+D1PP(J)*D1VP(J)
      DO 20 J=2,NY-1
      D2PM(J)=D1PM(J)*D1VM(J-1)
      D2P0(J)=D1PM(J)*D1VP(J-1)+D1PP(J)*D1VM(J)
      D2PP(J)=                  D1PP(J)*D1VP(J)
   20 CONTINUE
      J=NY
      D2PM(J)=D1PM(J)*D1VM(J-1)+D1PP(J)*D1VM(J)
      D2P0(J)=D1PM(J)*D1VP(J-1)+D1PP(J)*D1VP(J)
      RETURN
      END
C-----------------------------------------------------------------------
C---  FDM operators for pressure's Poisson equation  -----------------
      SUBROUTINE MESH2P
      PARAMETER (NY=64,NY1=NY+1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /DPP/DPPM(2:NY),DPP0(NY),DPPP(NY-1)
      COMMON /D1VN/D1VNM(0:NY),D1VNP(0:NY) /D1P/D1PM(NY),D1PP(NY)
      DO 10 J=1,NY
      IF(J.GT. 1) DPPM(J)=D1PM(J)*D1VNM(J-1)
      DPP0(J)=D1PM(J)*D1VNP(J-1)+D1PP(J)*D1VNM(J)
      IF(J.LT.NY) DPPP(J)=D1PP(J)*D1VNP(J)
   10 CONTINUE
      RETURN
      END
C---------------------------------------------------------------------
C-----  Mesh generation .....  for Y direction  ----------------------
      SUBROUTINE MESHIN
      PARAMETER (NX=32,NY=64,NZ=32,NY1=NY+1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /DY/DY(NY) /DYC/DYC(NY-1)
      COMMON /D1INP/D1IMP(NY),D1I0P(NY),D1IPP(NY)
      COMMON /D1INV/D1IMV(NY-1),D1I0V(NY-1),D1IPV(NY-1)
      DO 10 J=1,NY
      IF(J.EQ.1) THEN
         HM=DY(1)/2.D0
         ELSE
         HM=DYC(J-1)
         ENDIF
      IF(J.EQ.NY) THEN
         HP=DY(NY)/2.D0
         ELSE
         HP=DYC(J)
         ENDIF
      D1IMP(J)=-HP/HM/(HM+HP)
      D1I0P(J)=(-HM+HP)/HM/HP
      D1IPP(J)= HM/HP/(HM+HP)
   10 CONTINUE

      DO 20 J=1,NY-1
         HM=DY(J)
         HP=DY(J+1)
      D1IMV(J)=-HP/HM/(HM+HP)
      D1I0V(J)=(-HM+HP)/HM/HP
      D1IPV(J)= HM/HP/(HM+HP)
   20 CONTINUE
      RETURN
      END

C *********************************************************************
C *     SBR. DTR  :  DATA READ                                        *
C *********************************************************************

      SUBROUTINE SBRDTR
      PARAMETER (NX=32,NY=64,NZ=32,NY1=NY+1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*6 FILE1
      COMMON /FILE1/FILE1 /START/ISTART /TSTEP/TSTEP
      COMMON /U/U(NZ,NX,NY) /V/V(NZ,NX,0:NY) /W/W(NZ,NX,NY)
      COMMON /P/P(NZ,NX,NY) 
      OPEN(10,file=FILE1//'u.d',status='old',form='formatted')
        READ(10,1000) ISTART,TSTEP
        READ(10,2000) U
        CLOSE(10)
      OPEN(10,file=FILE1//'v.d',status='old',form='formatted')
        READ(10,1000) ISTART,TSTEP
        READ(10,2000) V
        CLOSE(10)
      OPEN(10,file=FILE1//'w.d',status='old',form='formatted')
        READ(10,1000) ISTART,TSTEP
        READ(10,2000) W
        CLOSE(10)
      OPEN(10,file=FILE1//'p.d',status='old',form='formatted')
        READ(10,1000) ISTART,TSTEP
        READ(10,3000) P
        CLOSE(10)
 1000 FORMAT(I10,F12.6)
 2000 FORMAT(8F10.6)
 3000 FORMAT(8F10.5)
      WRITE(6,1200) ISTART,TSTEP
 1200 FORMAT(/'Data READ',10X,'STEP=',I7,10X,'T=',F10.3)
      RETURN
      END

C *********************************************************************
C *     SBR. CON  :  NONLINEAR TERM  --- CONSISTENT FORM              *
C *********************************************************************

      SUBROUTINE SBRCON
        CALL SBRNLU
        CALL SBRNLV
        CALL SBRNLW
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE SBRNLU
      PARAMETER (NX=32,NY=64,NZ=32)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /UF/UF(NZ,NX,NY)
      COMMON /U/U(NZ,NX,NY) /V/V(NZ,NX,0:NY) /W/W(NZ,NX,NY)
      COMMON /UP/UP(NZ,NX,NY) /UV/UV(NZ,NX,0:NY) /WU/WU(NZ,NX,NY)
      COMMON /TXZA/TXA,TZA /BXZA/BXA,BZA 
      COMMON /D1V/D1VM(0:NY),D1VP(0:NY) /D0P/D0PM(NY),D0PP(NY)
      COMMON /IP/IP(NX) /IM/IM(NX) /KP/KP(NZ) /KM/KM(NZ)

C ... UP:-U^XU_X  at P,  WU:-W^XU_Z  at P
      DO 10 J=1,NY
      DO 10 I=1,NX
      DO 10 K=1,NZ
      UP(K,I,J)=-TXA*(U(K,IM(I),J)+U(K,I,J))*(-U(K,IM(I),J)+U(K,I,J))
      WU(K,I,J)=-TZA*(W(K,I,J)+W(K,IP(I),J))*(-U(K,I,J)+U(KP(K),I,J))
   10 CONTINUE

C ... UV:-V^XU_Y  at UV
      DO 31 J=1,NY-1
      DO 31 I=1,NX
      DO 31 K=1,NZ
      UV(K,I,J)=-(D1VM(J)*U(K,I,J)+D1VP(J)*U(K,I,J+1))
     &          *(V(K,I,J)+V(K,IP(I),J))/2.D0
   31 CONTINUE
      DO 32 J=0,NY,NY
      DO 32 I=1,NX
      DO 32 K=1,NZ
      UV(K,I,J)=0.D0
   32 CONTINUE

C ... UF:-(U^XU_X)^X-(V^XU_Y)^Y-(W^XU_Z)^Z  at U
      DO 50 J=1,NY
      DO 50 I=1,NX
      DO 50 K=1,NZ
      UF(K,I,J)=(UP(K,I,J)+UP(K,IP(I),J)+WU(KM(K),I,J)+WU(K,I,J))/2.D0
     &         +D0PM(J)*UV(K,I,J-1)+D0PP(J)*UV(K,I,J)
   50 CONTINUE
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE SBRNLV
      PARAMETER (NX=32,NY=64,NZ=32)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /VF/VF(NZ,NX,NY-1)
      COMMON /U/U(NZ,NX,NY) /V/V(NZ,NX,0:NY) /W/W(NZ,NX,NY)
      COMMON /VP/VP(NZ,NX,NY) /UV/UV(NZ,NX,0:NY) /VW/VW(NZ,NX,0:NY)
      COMMON /TXZA/TXA,TZA /BXZA/BXA,BZA 
      COMMON /D0V/D0VM(NY-1),D0VP(NY-1) 
      COMMON /D0P/D0PM(NY),D0PP(NY) /D1P/D1PM(NY),D1PP(NY)
      COMMON /IP/IP(NX) /IM/IM(NX) /KP/KP(NZ) /KM/KM(NZ)

C ... VP:-(V^Y*V_Y)^Y at P
      DO 22 J=1,NY
      DO 22 I=1,NX
      DO 22 K=1,NZ
      VP(K,I,J)=-(D0PM(J)*V(K,I,J-1)+D0PP(J)*V(K,I,J))
     &          *(D1PM(J)*V(K,I,J-1)+D1PP(J)*V(K,I,J))
   22 CONTINUE

C ... UV:-(U^Y*V_X) at UV, VW:-(W^Y*V_Z)  at VW
      DO 30 J=1,NY-1
      DO 30 I=1,NX
      DO 30 K=1,NZ
      UV(K,I,J)=-(D0VM(J)*U(K,I,J)+D0VP(J)*U(K,I,J+1))
     &           *BXA*(-V(K,I,J)+V(K,IP(I),J))
      VW(K,I,J)=-(D0VM(J)*W(K,I,J)+D0VP(J)*W(K,I,J+1))
     &           *BZA*(-V(K,I,J)+V(KP(K),I,J))
   30 CONTINUE

C ... VF:-(U^Y*V_X)^X-(V^Y*V_Y)^Y-(W^Y*V_Z)^Z  at V
      DO 50 J=1,NY-1
      DO 50 I=1,NX
      DO 50 K=1,NZ
      VF(K,I,J)=(UV(K,IM(I),J)+UV(K,I,J)+VW(KM(K),I,J)+VW(K,I,J))/2.D0
     &         +D0VM(J)*VP(K,I,J)+D0VP(J)*VP(K,I,J+1)
   50 CONTINUE
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE SBRNLW
      PARAMETER (NX=32,NY=64,NZ=32)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /WF/WF(NZ,NX,NY)
      COMMON /U/U(NZ,NX,NY) /V/V(NZ,NX,0:NY) /W/W(NZ,NX,NY)
      COMMON /WU/WU(NZ,NX,NY) /VW/VW(NZ,NX,0:NY) /WP/WP(NZ,NX,NY)
      COMMON /TXZA/TXA,TZA /BXZA/BXA,BZA 
      COMMON /D1V/D1VM(0:NY),D1VP(0:NY) /D0P/D0PM(NY),D0PP(NY)
      COMMON /IP/IP(NX) /IM/IM(NX) /KP/KP(NZ) /KM/KM(NZ)

C ... WU:-U^ZW_X  at UW,  WP:-W^ZW_Z  at P
      DO 10 J=1,NY
      DO 10 I=1,NX
      DO 10 K=1,NZ
      WU(K,I,J)=-TXA*(U(K,I,J)+U(KP(K),I,J))*(-W(K,I,J)+W(K,IP(I),J))
      WP(K,I,J)=-TZA*(W(KM(K),I,J)+W(K,I,J))*(-W(KM(K),I,J)+W(K,I,J))
   10 CONTINUE

C ... VW:-V^ZW_Y  at VW
      DO 31 J=1,NY-1
      DO 31 I=1,NX
      DO 31 K=1,NZ
      VW(K,I,J)=-(D1VM(J)*W(K,I,J)+D1VP(J)*W(K,I,J+1))
     &         *(V(K,I,J)+V(KP(K),I,J))/2.D0
   31 CONTINUE
      DO 32 J=0,NY,NY
      DO 32 I=1,NX
      DO 32 K=1,NZ
      VW(K,I,J)=0.D0
   32 CONTINUE

C ... WF:-(U^ZW_X)^X-(V^ZW_Y)^Y-(W^ZW_Z)^Z  at W
      DO 50 J=1,NY
      DO 50 I=1,NX
      DO 50 K=1,NZ
      WF(K,I,J)=(WU(K,IM(I),J)+WU(K,I,J)+WP(K,I,J)+WP(KP(K),I,J))/2.D0
     &         +D0PM(J)*VW(K,I,J-1)+D0PP(J)*VW(K,I,J)
   50 CONTINUE
      RETURN
      END

C *********************************************************************
C *     SBR. VGT  :  VELOCITY GRADIENT TENSOR                         *
C *********************************************************************
C
      SUBROUTINE SBRVGT
      PARAMETER (NX=32,NY=64,NZ=32)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /U/U(NZ,NX,NY) /V/V(NZ,NX,0:NY) /W/W(NZ,NX,NY)
      COMMON /UP/UP(NZ,NX,NY) /VP/VP(NZ,NX,NY) /WP/WP(NZ,NX,NY)
     &/UV/UV(NZ,NX,0:NY) /VW/VW(NZ,NX,0:NY) /WU/WU(NZ,NX,NY)
      COMMON /D/DX,DZ /DY/DY(NY)
      COMMON /D1V/D1VM(0:NY),D1VP(0:NY) /DYC/DYC(NY-1)
      COMMON /IP/IP(NX) /IM/IM(NX) /KP/KP(NZ) /KM/KM(NZ)
C
      DO 10 J=1,NY
      DO 10 I=1,NX
      DO 10 K=1,NZ
      UP(K,I,J)=(-U(K,IM(I),J)+U(K,I,J))/DX
      VP(K,I,J)=(-V(K,I,J-1)+V(K,I,J))/DY(J)
      WP(K,I,J)=(-W(KM(K),I,J)+W(K,I,J))/DZ
      WU(K,I,J)=
     & (-W(K,I,J)+W(K,IP(I),J))/DX+(-U(K,I,J)+U(KP(K),I,J))/DZ
   10 CONTINUE
C
      DO 20 J=1,NY-1
      DO 20 I=1,NX
      DO 20 K=1,NZ
      UV(K,I,J)=
     & (-U(K,I,J)+U(K,I,J+1))/DYC(J)+(-V(K,I,J)+V(K,IP(I),J))/DX
      VW(K,I,J)=
     & (-V(K,I,J)+V(KP(K),I,J))/DZ+(-W(K,I,J)+W(K,I,J+1))/DYC(J)
   20 CONTINUE
      DO 30 I=1,NX
      DO 30 K=1,NZ
      UV(K,I, 0)=D1VM(0)*U(K,I,1)+D1VP(0)*U(K,I,2)
      UV(K,I,NY)=D1VM(NY)*U(K,I,NY-1)+D1VP(NY)*U(K,I,NY)
      VW(K,I, 0)=D1VM(0)*W(K,I,1)+D1VP(0)*W(K,I,2)
      VW(K,I,NY)=D1VM(NY)*W(K,I,NY-1)+D1VP(NY)*W(K,I,NY)
   30 CONTINUE
      RETURN
      END

C *********************************************************************
C *     SBR. SGS  :  SGS EDDY VISCOSITY COEFFICIENT                   *
C *********************************************************************
C
      SUBROUTINE SBRSGS
      PARAMETER (NX=32,NY=64,NZ=32,NY1=NY+1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /RET/RET /K/AK(NZ,NX,NY)
      COMMON /Y/YV(0:NY),YP(0:NY1) /D/DX,DZ /DY/DY(NY) /CAK/CAK(NY)
      COMMON /IP/IP(NX) /IM/IM(NX) /KP/KP(NZ) /KM/KM(NZ)
      COMMON /UP/UP(NZ,NX,NY) /VP/VP(NZ,NX,NY) /WP/WP(NZ,NX,NY)
     &/UV/UV(NZ,NX,0:NY) /VW/VW(NZ,NX,0:NY) /WU/WU(NZ,NX,NY)
C
      DO 10 J=1,NY
      DO 10 I=1,NX
      DO 10 K=1,NZ
      DUYDVX=(UV(K,IM(I),J-1)+UV(K,I,J-1)+UV(K,IM(I),J)+UV(K,I,J))/4.D0
      DVZDWY=(VW(KM(K),I,J-1)+VW(K,I,J-1)+VW(KM(K),I,J)+VW(K,I,J))/4.D0
      DWXDUZ=(WU(KM(K),IM(I),J)+WU(K,IM(I),J)+WU(KM(K),I,J)+WU(K,I,J))
     &/4.D0
      AHD=2.D0*(UP(K,I,J)**2+VP(K,I,J)**2+WP(K,I,J)**2)
     &         +DUYDVX**2+DVZDWY**2+DWXDUZ**2
      AK(K,I,J)=CAK(J)*DSQRT(AHD)
   10 CONTINUE
      RETURN
      END

C *********************************************************************
C *     SBR. FLX  :  MOMENTUM FLUX EXEPT FOR PRESSURE                 *
C *********************************************************************

      SUBROUTINE SBRFLX
      PARAMETER (NX=32,NY=64,NZ=32)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /U/U(NZ,NX,NY) /V/V(NZ,NX,0:NY) /W/W(NZ,NX,NY)
      COMMON /UP/UP(NZ,NX,NY) /VP/VP(NZ,NX,NY) /WP/WP(NZ,NX,NY)
     &/UV/UV(NZ,NX,0:NY) /VW/VW(NZ,NX,0:NY) /WU/WU(NZ,NX,NY)
      COMMON /RET/RET /K/AK(NZ,NX,NY)
      COMMON /D/DX,DZ /DY/DY(NY)
      COMMON /D1V/D1VM(0:NY),D1VP(0:NY) /DYC/DYC(NY-1)
      COMMON /IP/IP(NX) /IM/IM(NX) /KP/KP(NZ) /KM/KM(NZ)

      DO 10 J=1,NY
      DO 10 I=1,NX
      DO 10 K=1,NZ
      AS=2.D0*(1.D0/RET+AK(K,I,J))
C *** Conservative form
C *     US=(U(K,IM(I),J)+U(K,I,J))/2.
C *     VS=(V(K,I,J-1)+V(K,I,J))/2.
C *     WS=(W(KM(K),I,J)+W(K,I,J))/2.
C *     UP(K,I,J)=-US**2+AS*UP(K,I,J)
C *     VP(K,I,J)=-VS**2+AS*VP(K,I,J)
C *     WP(K,I,J)=-WS**2+AS*WP(K,I,J)
      UP(K,I,J)=AS*UP(K,I,J)
      VP(K,I,J)=AS*VP(K,I,J)
      WP(K,I,J)=AS*WP(K,I,J)
      AS=1.D0/RET
     &+(AK(K,I,J)+AK(K,IP(I),J)+AK(KP(K),I,J)+AK(KP(K),IP(I),J))/4.D0
C *     WU(K,I,J)=-(W(K,I,J)+W(K,IP(I),J))*(U(K,I,J)+U(KP(K),I,J))/4.D0
C *    &+AS*WU(K,I,J)
      WU(K,I,J)=AS*WU(K,I,J)
   10 CONTINUE

      DO 20 J=1,NY-1
      DO 20 I=1,NX
      DO 20 K=1,NZ
      AS=1.D0/RET
     &+(AK(K,I,J)+AK(K,IP(I),J)+AK(K,I,J+1)+AK(K,IP(I),J+1))/4.D0
C *     UV(K,I,J)=-(U(K,I,J)+U(K,I,J+1))*(V(K,I,J)+V(K,IP(I),J))/4.D0
C *    &+AS*UV(K,I,J)
      UV(K,I,J)=AS*UV(K,I,J)

      AS=1.D0/RET
     &+(AK(K,I,J)+AK(K,I,J+1)+AK(KP(K),I,J)+AK(KP(K),I,J+1))/4.D0
C *     VW(K,I,J)=-(V(K,I,J)+V(KP(K),I,J))*(W(K,I,J)+W(K,I,J+1))/4.D0
C *    &+AS*VW(K,I,J)
      VW(K,I,J)=AS*VW(K,I,J)
   20 CONTINUE

      DO 30 J=0,NY,NY
      DO 30 I=1,NX
      DO 30 K=1,NZ
      UV(K,I,J)=UV(K,I,J)/RET
      VW(K,I,J)=VW(K,I,J)/RET
   30 CONTINUE
      RETURN
      END

C *********************************************************************
C *     SBR. VIS  :  VISCOUS TERM                                     *
C *********************************************************************

      SUBROUTINE SBRVIS
      PARAMETER (NX=32,NY=64,NZ=32,NY1=NY+1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /UF/UF(NZ,NX,NY) /VF/VF(NZ,NX,NY-1) /WF/WF(NZ,NX,NY)
      COMMON /IP/IP(NX) /IM/IM(NX) /KP/KP(NZ) /KM/KM(NZ) 
      COMMON /D/DX,DZ /DY/DY(NY) /DYC/DYC(NY-1)
      COMMON /UP/UP(NZ,NX,NY) /VP/VP(NZ,NX,NY) /WP/WP(NZ,NX,NY)
     &/UV/UV(NZ,NX,0:NY) /VW/VW(NZ,NX,0:NY) /WU/WU(NZ,NX,NY)
      DO 10 J=1,NY
      DO 10 I=1,NX
      DO 10 K=1,NZ
      UF(K,I,J)=UF(K,I,J)+(-UP(K,I,J)+UP(K,IP(I),J))/DX
     &+(-UV(K,I,J-1)+UV(K,I,J))/DY(J)+(-WU(KM(K),I,J)+WU(K,I,J))/DZ
      WF(K,I,J)=WF(K,I,J)+(-WU(K,IM(I),J)+WU(K,I,J))/DX
     &+(-VW(K,I,J-1)+VW(K,I,J))/DY(J)+(-WP(K,I,J)+WP(KP(K),I,J))/DZ
   10 CONTINUE

      DO 20 J=1,NY-1
      DO 20 I=1,NX
      DO 20 K=1,NZ
      VF(K,I,J)=VF(K,I,J)+(-UV(K,IM(I),J)+UV(K,I,J))/DX
     &+(-VP(K,I,J)+VP(K,I,J+1))/DYC(J)+(-VW(KM(K),I,J)+VW(K,I,J))/DZ
   20 CONTINUE
      RETURN
      END

C *********************************************************************
C *     SBR. PRE  :  PREDICTION STEP                                  *
C *********************************************************************

      SUBROUTINE SBRPRE(AB,BB)
      PARAMETER (NX=32,NY=64,NZ=32)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /DT/DT /U/U(NZ,NX,NY) /V/V(NZ,NX,0:NY) /W/W(NZ,NX,NY)
      COMMON /UF/UF(NZ,NX,NY) /VF/VF(NZ,NX,NY-1) /WF/WF(NZ,NX,NY)
      COMMON /UB/UB(NZ,NX,NY) /VB/VB(NZ,NX,NY-1) /WB/WB(NZ,NX,NY)
      DAB=AB*DT
      DBB=BB*DT
      DO 10 J=1,NY
      DO 10 I=1,NX
      DO 10 K=1,NZ
      UF(K,I,J)=UF(K,I,J)+2.D0
      U(K,I,J)=U(K,I,J)+DAB*UF(K,I,J)+DBB*UB(K,I,J)
      W(K,I,J)=W(K,I,J)+DAB*WF(K,I,J)+DBB*WB(K,I,J)
      UB(K,I,J)=UF(K,I,J)
      WB(K,I,J)=WF(K,I,J)
   10 CONTINUE

      DO 20 J=1,NY-1
      DO 20 I=1,NX
      DO 20 K=1,NZ
      V(K,I,J)=V(K,I,J)+DAB*VF(K,I,J)+DBB*VB(K,I,J)
      VB(K,I,J)=VF(K,I,J)
   20 CONTINUE
      RETURN
      END

C *********************************************************************
C *     SBR. RHP  :  R.H.S. OF POISSON EQ.                            *
C *********************************************************************

C --- R.H.S. FOR SCALER POTANTIAL -------------------------------------
      SUBROUTINE SBRRHP
      PARAMETER (NX=32,NY=64,NZ=32)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /U/U(NZ,NX,NY) /V/V(NZ,NX,0:NY) /W/W(NZ,NX,NY)
      COMMON /W3/Q(NZ,NX,NY) /DT/DT /D/DX,DZ
      COMMON /D1P/D1PM(NY),D1PP(NY) /IM/IM(NX) /KM/KM(NZ)
      DO 10 J=1,NY
      DO 10 I=1,NX
      DO 10 K=1,NZ
      Q(K,I,J)=((-U(K,IM(I),J)+U(K,I,J))/DX
     &         +D1PM(J)*V(K,I,J-1)+D1PP(J)*V(K,I,J)
     &         +(-W(KM(K),I,J)+W(K,I,J))/DZ )/DT
   10 CONTINUE
      RETURN
      END

C *********************************************************************
C *     SBR. SOR  :  S.O.R. SCHEME FOR POISSON EQ.                    *
C *********************************************************************

      SUBROUTINE SBRSOR(ISOR)
      PARAMETER (NX=32,NY=64,NZ=32,NY1=NY+1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /DPP/DPPM(2:NY),DPP0(NY),DPPP(NY-1)
      COMMON /P/P(NZ,NX,NY) /W3/Q(NZ,NX,NY) /ERRP/POIERR /ITRP/ITRP
      COMMON /IP/IP(NX) /IM/IM(NX) /CMPX/CM1X,C00X,CP1X
      COMMON /KP/KP(NZ) /KM/KM(NZ) /CMPZ/CM1Z,C00Z,CP1Z
      ITRP=0
      SUMS=0.D0
      DO 10 J=1,NY
      DO 10 I=1,NX
      DO 10 K=1,NZ
      SUMS=SUMS+Q(K,I,J)**2
   10 CONTINUE
  100 ITRP=ITRP+1

      SUMR=0.
      DO 20 J=1,NY
      C00=C00X+C00Z+DPP0(J)
      DPC=-1.5D0/C00
      IF(J.EQ.1) THEN
      DO 31 I=1,NX
C      DO 31 L=1,2
C *** CHECKER BORD S.O.R. METHOD
C    *VDIR NODEP(P)
C      DO 31 K=L,NZ,2
      DO 31 K=1,NZ
      RESI=CM1X*P(K,IM(I),J)+CM1Z*P(KM(K),I,J)+C00*P(K,I,J)
     &    +CP1Z*P(KP(K),I,J)+CP1X*P(K,IP(I),J)+DPPP(J)*P(K,I,J+1)
     &    -Q(K,I,J)
      P(K,I,J)=P(K,I,J)+RESI*DPC
      SUMR=SUMR+RESI**2
   31 CONTINUE

      ELSE
      IF(J.LT.NY) THEN
      DO 30 I=1,NX
C      DO 30 L=1,2
C *** CHECKER BORD S.O.R. METHOD
C     *VDIR NODEP(P)
C      DO 30 K=L,NZ,2
      DO 30 K=1,NZ
      RESI=DPPM(J)*P(K,I,J-1)+CM1X*P(K,IM(I),J)+CM1Z*P(KM(K),I,J)
     &    +C00*P(K,I,J)
     &    +CP1Z*P(KP(K),I,J)+CP1X*P(K,IP(I),J)+DPPP(J)*P(K,I,J+1)
     &    -Q(K,I,J)
      P(K,I,J)=P(K,I,J)+RESI*DPC
      SUMR=SUMR+RESI**2
   30 CONTINUE

      ELSE
      DO 32 I=1,NX
C      DO 32 L=1,2
C *** CHECKER BORD S.O.R. METHOD
C    *VDIR NODEP(P)
C      DO 32 K=L,NZ,2
      DO 32 K=1,NZ
      RESI=DPPM(J)*P(K,I,J-1)+CM1X*P(K,IM(I),J)+CM1Z*P(KM(K),I,J)
     &    +C00*P(K,I,J)+CP1Z*P(KP(K),I,J)+CP1X*P(K,IP(I),J)
     &    -Q(K,I,J)
      P(K,I,J)=P(K,I,J)+RESI*DPC
      SUMR=SUMR+RESI**2
   32 CONTINUE
      ENDIF
      ENDIF
   20 CONTINUE

      POIERR=DSQRT(SUMR/SUMS)
      IF(ITRP.GE.ISOR.OR.POIERR.LT.1.D-5) GO TO 200
      GO TO 100
  200 CONTINUE
      RETURN
      END

C *********************************************************************
C *     SBR. COR  :  CORRECTION STEP                                  *
C *********************************************************************

      SUBROUTINE SBRCOR
      PARAMETER (NX=32,NY=64,NZ=32,NY1=NY+1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /U/U(NZ,NX,NY) /V/V(NZ,NX,0:NY) /W/W(NZ,NX,NY)
      COMMON /P/P(NZ,NX,NY) /DT/DT /D/DX,DZ /DY/DY(NY)
      COMMON /IP/IP(NX) /KP/KP(NZ) /D1VN/D1VNM(0:NY),D1VNP(0:NY)
      TX=DT/DX
      TZ=DT/DZ
      DO 10 J=1,NY
      DO 10 I=1,NX
      DO 10 K=1,NZ
      U(K,I,J)=U(K,I,J)-TX*(-P(K,I,J)+P(K,IP(I),J))
      W(K,I,J)=W(K,I,J)-TZ*(-P(K,I,J)+P(KP(K),I,J))
   10 CONTINUE

      DO 20 J=1,NY-1
      DO 20 I=1,NX
      DO 20 K=1,NZ
      V(K,I,J)=V(K,I,J)-DT*(D1VNM(J)*P(K,I,J)+D1VNP(J)*P(K,I,J+1))
   20 CONTINUE
      RETURN
      END

C *********************************************************************
C *     SBR. UMR  :  MEAN & RMS VALUES                                *
C *********************************************************************

      SUBROUTINE SBRUMR
      PARAMETER (NX=32,NY=64,NZ=32,NY1=NY+1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /U/U(NZ,NX,NY) /V/V(NZ,NX,0:NY) /W/W(NZ,NX,NY)
      COMMON /K/AK(NZ,NX,NY) /D/DX,DZ
      COMMON /UM/UM(NY),VM(0:NY),WM(NY)
      COMMON /S12/S12L(0:NY),S12S(0:NY),S12V(0:NY) /AKM/AKM(NY)
      COMMON /D0V/D0VM(NY-1),D0VP(NY-1) /DYC/DYC(NY-1)
      COMMON /IP/IP(NX) /IM/IM(NX) /KP/KP(NZ) /KM/KM(NZ)
      COMMON /BVU/BVU2(NY),BVU3(NY),BVU4(NY)
      COMMON /BVV/BVV2(0:NY),BVV3(0:NY),BVV4(0:NY)
      COMMON /BVW/BVW2(NY),BVW3(NY),BVW4(NY)
      ARXZ=1./DBLE(NX*NZ)
      DO 20 J=1,NY
      SUMU1=0.D0
      SUMU2=0.D0
      SUMW1=0.D0
      SUMW2=0.D0
      SUMA1=0.D0
      DO 10 I=1,NX
      DO 10 K=1,NZ
      SUMU1=SUMU1+U(K,I,J)
      SUMU2=SUMU2+U(K,I,J)**2
      SUMW1=SUMW1+W(K,I,J)
      SUMW2=SUMW2+W(K,I,J)**2
      SUMA1=SUMA1+AK(K,I,J)
   10 CONTINUE
      UM(J)=ARXZ*SUMU1
      BVU2(J)=ARXZ*SUMU2
      WM(J)=ARXZ*SUMW1
      BVW2(J)=ARXZ*SUMW2
      AKM(J)=ARXZ*SUMA1
   20 CONTINUE
      DO 50 J=1,NY-1
      SUMV1=0.D0
      SUMV2=0.D0
      SUMUV=0.D0
      SUMRS=0.D0
      DO 40 I=1,NX
      DO 40 K=1,NZ
      SUMV1=SUMV1+V(K,I,J)
      SUMV2=SUMV2+V(K,I,J)**2
      SUMUV=SUMUV+(D0VM(J)*U(K,I,J)+D0VP(J)*U(K,I,J+1))
     &           *(V(K,I,J)+V(K,IP(I),J))
      SUMRS=SUMRS
     &- (AK(K,I,J)+AK(K,IP(I),J)+AK(K,I,J+1)+AK(K,IP(I),J+1))/4.D0
     &*((-U(K,I,J)+U(K,I,J+1))/DYC(J)+(-V(K,I,J)+V(K,IP(I),J))/DX)
   40 CONTINUE
      VM(J)=ARXZ*SUMV1
      BVV2(J)=ARXZ*SUMV2
      S12L(J)=ARXZ*SUMUV/2.D0
      S12S(J)=ARXZ*SUMRS
   50 CONTINUE
      RETURN
      END

C *********************************************************************
C *     SBR. CRR  :  AUTO CORRELATION                                 *
C *********************************************************************
C      SUBROUTINE SBRCRR
C *********************************************************************
C *     SBR. SFP  :  SKEWNESS & FLATNESS,   & STATISTICS FOR PRESSURE *
C *********************************************************************
C      SUBROUTINE SBRSFP
C *********************************************************************
C *     SBR. DSS  :  ENERGY DISSIPATION RATE                          *
C *********************************************************************
C      SUBROUTINE SBRDSS
C *********************************************************************
C *     SBR. DIF  :  ENERGY DIFFUSION FLUX                            *
C *********************************************************************
C      SUBROUTINE SBRDIF

C *********************************************************************
C *     SBR. STT  :  TURBULENCE STATISTICS                            *
C *********************************************************************
C --- Mean Velocity and Turbulence Intensity----------------------------
      SUBROUTINE SBRST1
      PARAMETER (NX=32,NY=64,NZ=32,NY1=NY+1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /UM/UM(NY),VM(0:NY),WM(NY) /PM/PM(NY)
      COMMON /UR/UR(NY),VR(0:NY),WR(NY) /PR/PR(NY)
      COMMON /ENE/ENE /UME/UME,UMX /RMS/URMSX,URMSC,VRMSC,WRMSC
      COMMON /DY/DY(NY) /DYC/DYC(NY-1)
      COMMON /BVU/BVU2(NY),BVU3(NY),BVU4(NY)
      COMMON /BVV/BVV2(0:NY),BVV3(0:NY),BVV4(0:NY)
      COMMON /BVW/BVW2(NY),BVW3(NY),BVW4(NY)
      COMMON /BVP/BVP2(NY)
        ENE=0.D0
        UME=0.D0
        UMX=0.D0
        URMSX=0.D0
      DO 10 J=1,NY
      UR(J)=DSQRT(BVU2(J)-UM(J)**2)
      WR(J)=DSQRT(BVW2(J)-WM(J)**2)
      PR(J)=DSQRT(BVP2(J)-PM(J)**2)
      ENE=ENE+DY(J)*(UR(J)**2+WR(J)**2)
      UME=UME+DY(J)*UM(J)
      UMX=DMAX1(UMX,UM(J))
      URMSX=DMAX1(URMSX,UR(J))
   10 CONTINUE
      DO 20 J=1,NY-1
      VR(J)=DSQRT(BVV2(J)-VM(J)**2)
      ENE=ENE+DYC(J)*VR(J)**2
   20 CONTINUE
      URMSC=(UR(NY/2)+UR(NY/2+1))/2.D0
      VRMSC=VR(NY/2)
      WRMSC=(WR(NY/2)+WR(NY/2+1))/2.D0
      ENE=ENE/2.D0
      RETURN
      END

C --- Skewness, Flatness Factors and Production Terms ------------------
C      SUBROUTINE SBRST2
C --- Dissipation Terms ------------------------------------------------
C      SUBROUTINE SBRST3
C --- Diffusion and Pressure-Related Terms -----------------------------
C      SUBROUTINE SBRST4
C --- Auto-Correlation coefficient  ------------------------------------
C      SUBROUTINE SBRST5
C --- Averagr in symmetric plane ---------------------------------------
C      SUBROUTINE TAISHO

C *********************************************************************
C *     SBR. CHK  :  CHECK of DIVERGENCE and COURANT-NUMBER           *
C *********************************************************************

      SUBROUTINE SBRCHK
      PARAMETER (NX=32,NY=64,NZ=32)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /U/U(NZ,NX,NY) /V/V(NZ,NX,0:NY) /W/W(NZ,NX,NY)
      COMMON /D/DX,DZ /DY/DY(NY) /DT/DT /DIV/DIVX /COU/COMX
      COMMON /D0P/D0PM(NY),D0PP(NY) /D1P/D1PM(NY),D1PP(NY)
      COMMON /UR/UR(NY),VR(0:NY),WR(NY)
      COMMON /IM/IM(NX) /KM/KM(NZ)
      DIVX=0.D0
      COMX=0.D0
      DO 40 J=1,NY
         DNOMAL=(DX*DY(J)*DZ)**(1.D0/3.D0)/UR(J)
      DO 40 I=1,NX
      DO 40 K=1,NZ
      DIV=(-U(K,IM(I),J)+U(K,I,J))/DX
     &   +D1PM(J)*V(K,I,J-1)+D1PP(J)*V(K,I,J)
     &   +(-W(KM(K),I,J)+W(K,I,J))/DZ
      UCP=(U(K,IM(I),J)+U(K,I,J))/2.D0
      WCP=(W(KM(K),I,J)+W(K,I,J))/2.D0
      VCP=D0PM(J)*V(K,I,J-1)+D0PP(J)*V(K,I,J)
      DIVX=DMAX1(DIVX,DNOMAL*DIV)
      COU=DT*(DABS(UCP)/DX+DABS(VCP)/DY(J)+DABS(WCP)/DZ)
      COMX=DMAX1(COMX,COU)
   40 CONTINUE
      RETURN
      END

C *********************************************************************
C *     SBR. DTW  :  DATA  PRINT OUT                                  *
C *********************************************************************
C      SUBROUTINE SBRDTW

C *********************************************************************
C *     SBR. DTS  :  DATA SAVE                                        *
C *********************************************************************

      SUBROUTINE SBRDTS
      PARAMETER (NX=32,NY=64,NZ=32,NY1=NY+1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*6  FILE2
      COMMON /FILE2/FILE2 /STEP/ISTEP /TSTEP/TSTEP
      COMMON /U/U(NZ,NX,NY) /V/V(NZ,NX,0:NY) /W/W(NZ,NX,NY)
      COMMON /P/P(NZ,NX,NY) 
      OPEN(10,file=FILE2//'u.d',status='unknown',form='formatted')
        WRITE(10,1000) ISTEP,TSTEP
        WRITE(10,2000) U
        CLOSE(10)
      OPEN(10,file=FILE2//'v.d',status='unknown',form='formatted')
        WRITE(10,1000) ISTEP,TSTEP
        WRITE(10,2000) V
        CLOSE(10)
      OPEN(10,file=FILE2//'w.d',status='unknown',form='formatted')
        WRITE(10,1000) ISTEP,TSTEP
        WRITE(10,2000) W
        CLOSE(10)
      OPEN(10,file=FILE2//'p.d',status='unknown',form='formatted')
        WRITE(10,1000) ISTEP,TSTEP
        WRITE(10,3000) P
        CLOSE(10)
 1000 FORMAT(I10,F12.6)
 2000 FORMAT(8F10.6)
 3000 FORMAT(8F10.5)
      WRITE(6,1200) ISTEP,TSTEP
 1200 FORMAT(/'Data SAVE',10X,'STEP=',I7,10X,'T=',F10.3/)
      RETURN
      END

C *********************************************************************
C *     SBR. TAV  :  TIME AVERAGING                                   *
C *********************************************************************
C      SUBROUTINE TAVCLR
C      SUBROUTINE TAVDTR
C      SUBROUTINE TAVRAG
