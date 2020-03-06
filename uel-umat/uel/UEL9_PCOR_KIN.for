C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C     GENERAL FINITE ELEMENT LIBRARY TO BE USED WITHIN ABAQUS          C
C     UNIVERSIDAD EAFIT                                                C
C     LABORATORIO DE MECANICA APLICADA                                 C
C     BLOQUE 14-PISO 2                                                 C
C     MEDELLIN, COLOMBIA                                               C
C                                                                      C
C     LAST UPDATED APRIL 16/2015                                       C
C                                                                      C
C     CURRENTLY:                                                       C
C                                                                      C
C     8 NODED ISOPARAMETRIC REDUCED INTEGRATION PENALTHY METHOD        C
C     REDUCED COSSERAT COUPLE STRESS SOLID ELEMENT                     C
C----------------------------------------------------------------------C
C                                                                      C
C----------------------------------------------------------------------C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C       --U S E R   E L E M E N T    S U B R O U T I N E S---          C
C                                                                      C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C  UEL9_PCOR.for  8-Noded Plasticity Couple Stress theory              C
C                 Return Mapping Algorithm Linear Isotropic/Kinematic  C
C                 Hardening                                       .    C
C                                                                      C
C  8 NODED ISOPARAMETRIC PENALTHY METHOD  COSSERAT SOLID ELEMENT       C
C  3X3 GAUSS INTEGRATION FOR STRESS AND COUPLE STRESS TERMS            C
C  2X2 GAUSS INTEGRATION FOR THE ASSYMETRIC STRESS COMPONENT           C
C                                                                      C
C  CREATED BY JUAN GOMEZ                                               C
C  STATE VARIABLES DEFINITIONS/GAUSS  POINT                            C
C                                                                      C
C    FULL INTEGRATION PART                                             C
C                                                                      C
C     1-6: STRESS VECTOR                       (6)                     C
C    7-12: TOTAL STRAIN VECTOR                 (6)                     C
C   13-18: ELASTIC STRAIN VECTOR               (6)                     C
C   19-24: PLASTIC STRAIN VECTOR               (6)                     C
C   25-30: BACK STRESS VECTOR                  (6)                     C
C      31: EQUIVALENT PLASTIC STRAIN           (1)                     C
C      32: EQUVALENT STRESS                    (1)                     C
C      33: PLASTIC LENGTH SCALE                (1)                     C
C                                   TOTAL     297 SVARS                C
C                                                                      C
C    REDUCED INTEGRATION PART(4 GAUSS POINTS)                          C
C                                                                      C
C      34: ASYMMETRIC STRESS                   (1)                     C
C      35: RELATIVE ROTATION                   (1)                     C
C                                   TOTAL       8 SVARS                C
C                                   TOTAL     305 SVARS                C
C                                                                      C
C    NODAL VALUES OF ONE STRESS AND STRAIN COMPONENT                   C
C                                                                      C
C     306: STRESS AT GAUSS POINT SPECIFIED                             C
C     307  STRAINT GAUSS POINT SPECIFIED                               C
C                                                                      C
C                                   TOTAL     307 SVARS                C
C                                                                      C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4PERIOD)
C     
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
C  
      PARAMETER (ZERO=0.D0,HALF=0.5D0,ONE=1.D0,NTENS=6,TWO=2.D0,
     1           THREE=3.D0)
C
C     Parameters required by UMAT.f
C
      PARAMETER (NDI=3,NSTATV=20,SSE=0.D0,SCD=0.D0,
     1           RPL=0.D0,DRPLDT=0.D0,TEMP=0.D0,DTEMP=0.D0,NSHR=1,
     2           CELENT=2.D0,LAYER=1,KSPT=1)
C
C     Parameter arrays from UEL.f
C
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),
     1SVARS(NSVARS),ENERGY(8),PROPS(*),COORDS(MCRD,NNODE),
     2U(NDOFEL),DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),
     3PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     4DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
     5JPROPS(*)
C
C     User defined arrays
C
      DIMENSION B(NTENS,NDOFEL),BT(NDOFEL,NTENS),BLM(1,NDOFEL),
     1BLMT(NDOFEL,1),FRST1(NDOFEL,NDOFEL),FRST2(NDOFEL,1),
     2XX(2,NNODE),XWF(9),XPF(2,9),AUX1(NTENS,NDOFEL),
     3AUX2(1,NDOFEL),AUX3(NDOFEL,NDOFEL),AUX4(1,NDOFEL),
     4AUX5(NDOFEL,NDOFEL),STRESS(NTENS),STRAN(NTENS),
     5DSTRAN(NTENS),RELROT(1),DRELROT(1),DDBAR(1,1),STRASYM(1),
     6FRST3(NDOFEL,NDOFEL),FRST4(NDOFEL,1),XWR(4),XPR(2,4),
     7EELAS(NTENS),EPLAS(NTENS),XBACK(NTENS)
C
C     Arrays to be used in UMAT.f
C
      DIMENSION STATEV(NSTATV),DDSDDT(NTENS),DRPLDE(NTENS),
     1UPREDEF(1), DPRED(1),UCOORDS(3),DROT(3, 3), DFGRD0(3, 3),
     2DFGRD1(3, 3),DDSDDE(NTENS, NTENS)
C
C     Initializes parameters required by UMAT.f
C
      CALL KCLEARV(STRESS,NTENS)
      CALL KCLEARV(STRAN,NTENS)
      CALL KCLEARV(RELROT,1)
      CALL KCLEARV(DRELROT,1)
      CALL KCLEARV(DSTRAN,NTENS)
      CALL KCLEARV(DDSDDT,NTENS)
      CALL KCLEARV(DRPLDE,NTENS)
      CALL KCLEAR(DROT,3,3)
      CALL KCLEAR(DFGRD0,3,3)
      CALL KCLEAR(DFGRD1,3,3)
      CALL KCLEAR(DDSDDE,NTENS,NTENS)
      CALL KCLEAR(DDBAR,1,1)
      UPREDEF(1)=ZERO
      DPRED(1)=ZERO
C
C     Clears RHS vector and Stiffness matrix
C
      CALL KCLEARV(RHS,NDOFEL)
      CALL KCLEAR(AMATRX,NDOFEL,NDOFEL)
      CALL KCLEAR(AUX1,NTENS,NDOFEL)
      CALL KCLEAR(AUX2,1,NDOFEL)
      CALL KCLEAR(AUX3,NDOFEL,NDOFEL)
      CALL KCLEAR(AUX4,2,NDOFEL)
      CALL KCLEAR(AUX5,NDOFEL,NDOFEL)
C
      NSVARS_N=2
      NSVARS_R=8
      NSVARS_F=NSVARS-NSVARS_R-NSVARS_N
C
      CMNAME='COSSERAT'
      PLS=PROPS(3)
      PENNUM=PROPS(4)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C               Starts full integration part.                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      CALL KGPOINTS3X3(XPF,XWF)
      NGPT_F=9
C
C     Loops around all Gauss points
C
      DO NN=1,NGPT_F
C
        RII=XPF(1,NN)
        SII=XPF(2,NN)
        ALF=XWF(NN)
C
C       Compute State variable index corresponding
C       to current Gauss point and load stress,total strain
C       elastic strain,plastic strain and equivalent plastic
C       strain from state variables as defined in USER ELEMENT. 
C       Different variables are required for different constitutive models.
C
        ISVINT_F=1+(NN-1)*NSVARS_F/NGPT_F
        JJ=1
        DO II=ISVINT_F,ISVINT_F+5
          STRESS(JJ)=SVARS(II)
          STRAN(JJ )=SVARS(II+6)
          EELAS(JJ )=SVARS(II+12)
          EPLAS(JJ )=SVARS(II+18)
          XBACK(JJ )=SVARS(II+24)
          JJ=JJ+1
        END DO
        EQPLAS=SVARS(ISVINT_F+30)
        SMISES=SVARS(ISVINT_F+31)
C
C       Plastic Length scale as a function of equivalent plastic strain.
C
        IF(ABS(EQPLAS).LE.0.00001) THEN
           PLSC=SVARS(ISVINT_F+32)+0.0000001
        ELSE
           PLSC=SVARS(ISVINT_F+32)
        END IF
C
C       Starts STATEV array definition as required by UMAT.f
C       Different variables for different constitutive models.
C
        DO II=1,NTENS
          STATEV(II)=EELAS(II)
          STATEV(II+NTENS)=EPLAS(II)
          STATEV(II+2*NTENS)=XBACK(II)
        END DO
        STATEV(3*NTENS+1)=EQPLAS
        STATEV(3*NTENS+2)=SMISES
C
C       Ends STATEV array definition as required by UMAT.f
C
C       Assembles strain-displacement and relative rotation-displacement
C       matrices B and BLM=(L-M)
C
        CALL KSTDM(JELEM,NNODE,NDOFEL,NTENS,COORDS,B,BLM,DDET,RII,
     1  SII,XBAR)
C
C       Captures length scale and penalty parameter.
C
        DO K1=1,NDOFEL
          B(5,K1)=B(5,K1)*PLS
          B(6,K1)=B(6,K1)*PLS
        END DO
C
C       Computes strain increment and rel.rotation increment
C       and updates strain and rel.rotation vectors.
C
        CALL KMAVEC(B,NTENS,NDOFEL,DU,DSTRAN)
        CALL KUPDVEC(STRAN,NTENS,DSTRAN)
C
C       Assembles material matrix and updates state variables
C
        CALL UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,
     1  DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP,
     2  DTEMP, UPREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV,
     3  PROPS, NPROPS, UCOORDS, DROT, PNEWDT, CELENT, DFGRD0,
     4  DFGRD1, JELEM, NN, LAYER, KSPT, KSTEP, KINC)
C
C       Assembles STIFFNESS matrix
c                                        T
C       Generalized stress component.===B M B (stress and couple stress)
C
        CALL KMMULT(DDSDDE,NTENS,NTENS,B,NTENS,NDOFEL,AUX1)
        CALL KMTRAN(B,NTENS,NDOFEL,BT)
        CALL KMMULT(BT,NDOFEL,NTENS,AUX1,NTENS,NDOFEL,FRST1)
C
C       Assembles RHS vector contribution
C
        CALL KMAVEC(BT,NDOFEL,NTENS,STRESS,FRST2)
C
C       Considers Gauss weight and Jacobian determinant representing
C       volume differential.
C
        CALL KSMULT(FRST1,NDOFEL,NDOFEL,ALF*DDET*XBAR)
        CALL KSMULT(FRST2,NDOFEL,1,-ALF*DDET*XBAR)
C
C       Updates Stiffness matrix and RHS vector
C       
        CALL KUPDMAT(AMATRX,NDOFEL,NDOFEL,FRST1)
        CALL KUPDVEC(RHS,NDOFEL,FRST2)
C
C       Clears material Jacobian and temporary stiffness matrix
C       array for new Gauss point 
C
        CALL KCLEAR(DDSDDE,NTENS,NTENS)
        CALL KCLEAR(FRST1,NDOFEL,NDOFEL)
        CALL KCLEAR(FRST2,NDOFEL,1)
C
C       Starts updating of state variables with updated values from UMAT.f
C
        DO II=1,NTENS
          EELAS(II)=STATEV(II)
          EPLAS(II)=STATEV(II+NTENS)
          XBACK(II)=STATEV(II+2*NTENS)
        END DO
        EQPLAS=STATEV(3*NTENS+1)
        SMISES=STATEV(3*NTENS+2)
C
        JJ=1
        DO II=ISVINT_F,ISVINT_F+5
          SVARS(II)=STRESS(JJ)
          SVARS(II+6  )=STRAN(JJ)
          SVARS(II+12 )=EELAS(JJ)
          SVARS(II+18 )=EPLAS(JJ)
          SVARS(II+24 )=XBACK(JJ)
          JJ=JJ+1
        END DO
        SVARS(ISVINT_F+30)=EQPLAS
        SVARS(ISVINT_F+31)=SMISES
C
C       Computes plastic length scale as a function of equivalent plastic strain
C
        SVARS(ISVINT_F+32)=PLS*(1.0-1.0*EXP(-1000.0*EQPLAS))
C
C       Ends updating of state variables with updated values from UMAT.f
C
      END DO
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C               Starts reduced integration part.                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      CALL KGPOINTS2X2(XPR,XWR)
      NGPT_R=4
C
C     Loops around all Gauss points
C
      DO NN=1,NGPT_R
C
        RII=XPR(1,NN)
        SII=XPR(2,NN)
        ALF=XWR(NN)
C
C       Compute State variable index
C       Different variables are required for different constitutive models.
C
        ISVINT_R=(ISVINT_F+32)+1+(NN-1)*NSVARS_R/NGPT_R
        JJ=1
        II=ISVINT_R
        STRASYM(JJ)=SVARS(II)
        RELROT(JJ) =SVARS(II+1)
C
C       Assembles strain-displacement and relative rotation-displacement
C       matrices B and BLM=(L-M)
C
        CALL KSTDM(JELEM,NNODE,NDOFEL,NTENS,COORDS,B,BLM,DDET,RII,
     1  SII,XBAR)
C
C       Computes relative rotation increment
C       and updates rel.rotation vector.
C
        CALL KMAVEC(BLM,1,NDOFEL,DU,DRELROT)
        CALL KUPDVEC(RELROT,1,DRELROT)
C
C       Updates asymmetric stress component.
C
        STRASYM(1)=STRASYM(1)+PENNUM*DRELROT(1)
C
C       Assembles STIFFNESS matrix
C                                                                  T_
C       Relative rotation component.(i.e constrained part) === (L-M)D(L-M)
C
        DDBAR(1,1)=PENNUM
        CALL KMMULT(DDBAR,1,1,BLM,1,NDOFEL,AUX4)
        CALL KMTRAN(BLM,1,NDOFEL,BLMT)
        CALL KMMULT(BLMT,NDOFEL,1,AUX4,1,NDOFEL,FRST3)
C
C       Assembles RHS vector contribution
C
        DO KK=1,NDOFEL
           FRST4(KK,1)=BLMT(KK,1)*STRASYM(1)
        END DO
C
C       Considers Gauss weight and Jacobian determinant representing
C       volume differential.
C
        CALL KSMULT(FRST3,NDOFEL,NDOFEL,ALF*DDET*XBAR)
        CALL KSMULT(FRST4,NDOFEL,1,-ALF*DDET*XBAR)
C
C       Updates Stiffness matrix and RHS vector
C       
        CALL KUPDMAT(AMATRX,NDOFEL,NDOFEL,FRST3)
        CALL KUPDVEC(RHS,NDOFEL,FRST4)
C
C       Clears material Jacobian and temporary stiffness matrix
C       array for new Gauss point 
C
        CALL KCLEAR(FRST3,NDOFEL,NDOFEL)
        CALL KCLEAR(FRST4,NDOFEL,1)
C
C       Starts updating of state variables with updated values from UMAT.f
C
        II=ISVINT_R
        SVARS(II)=STRASYM(JJ)
        SVARS(II+1)=RELROT(JJ)
C
      END DO
C
C     EXTRAPOLATE STRAIN TO THE NODES
C
      RII=-1.0
      SII= 1.0
      CALL KSTDM(JELEM,NNODE,NDOFEL,NTENS,COORDS,B,BLM,DDET,RII,
     1           SII,XBAR)
      CALL KMAVEC(B,NTENS,NDOFEL,U,STRAN)
      SVARS(NSVARS-1)=STRESS(1)
      SVARS(NSVARS  )=STRAN(1)
C
      RETURN
C      
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C     UMAT_PCOR_KIN.for    Plasticity Couple Stress Theory             C
C                          Return Mapping Algorithm                    C
C                                                                      C
C     SUBROUTINE UMAT                                                  C
C     PLANE STRAIN ANALYSIS                                            C
C     RATE INDEPENDENT LINAR ISOTROPIC/KINEMATIC HARDENING PLASTICITY  C
C                                                                      C
C     CREATED BY JUAN GOMEZ                                            C
C     June 24/2004                                                     C
C     ELECTRONIC PACKAGING LAB                                         C
C     102 KETTER HALL                                                  C
C     UNIVERSITY AT BUFFALO                                            C
C     BUFFALO, NY 14260                                                C
C                                                                      C
C     PLANE STRAIN                                                     C
C     NTENS: LENGTH OF STRESS VECTOR                                   C
C     NDI:   NUMBER OF NORMAL STRESS COMPONENTS                        C
C     LAST UPDATED JUNE 24/2004                                        C
C     USES RETURN MAPPING INTEGRATION SCHEME                           C
C                                                                      C
C     STATE VARIABLES DEFINITION                                       C
C                                                                      C
C     1- 6 : ELASTIC STRAINS                                           C
C     7-12 : PLASTIC STRAINS                                           C
C    13-18 : BACK STRESS                                               C
C       19 : EQUIVALENT PLASTIC STRAIN                                 C
C       20 : EQUIVALENT MISES STRESS                                   C
C                                                                      C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,
     1DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP,
     2DTEMP, PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV,
     3PROPS, NPROPS, COORDS, DROT, PNEWDT, CELENT, DFGRD0,
     4DFGRD1, NOEL, NPT, LAYER, KSPT, KSTEP, KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
C 
      DIMENSION STRESS(NTENS), STATEV(NSTATV),
     1DDSDDE(NTENS, NTENS),DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),
     2DSTRAN(NTENS),TIME(2),PREDEF(1), DPRED(1), PROPS(NPROPS),
     3COORDS(3),DROT(3, 3), DFGRD0(3, 3), DFGRD1(3, 3)
C
      DIMENSION AUX1(NTENS,NTENS),AUX2(NTENS,NTENS),AUX3(NTENS,NTENS),
     1AUX4(NTENS,NTENS),P(NTENS,NTENS),SINVAR(1,1),SDEV(NTENS,1),
     2XBACK(NTENS),STSREL(NTENS),DSTRTHER(NTENS),EELAS(NTENS),
     3EPLAS(NTENS),FLOW(NTENS)
C
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, SIX=6.0D0,
     1ENUMAX=0.4999D0, NEWTON=10, TOLER=1.0D-7,MAXITER=10)
C
C**********************************************************************
C     I S O T R O P I C  M I S E S  E L A S T O P L A S T I C I T Y
C     P L A N E  S T R A I N  A N A L Y S I S
C     C O U P L E  S T R E S S  T H E O R Y
C**********************************************************************
C
C     Recover equivalent plastic strain,elastic strains plastic
C     strains.
C
      EQPLAS=STATEV(3*NTENS+1)
      SMISES=STATEV(3*NTENS+2)
      DO K1=1, NTENS
        EELAS(K1)=STATEV(K1)
        EPLAS(K1)=STATEV(K1+NTENS)
        XBACK(K1)=STATEV(K1+2*NTENS)
        FLOW(K1)=ZERO
      END DO
C
C     Elastic properties
C
      EMOD=PROPS(1)
      ENU=MIN(PROPS(2),ENUMAX)
      SIG0=PROPS(5)
      SIGSAT=PROPS(6)
      HRDRATE=PROPS(7)
      HMOD=PROPS(8)
C
      EBULK3=EMOD/(ONE-TWO*ENU)
      EG2=EMOD/(ONE+ENU)
      EG=EG2/TWO
      EG3=THREE*EG
      ELAM=(EBULK3-EG2)/THREE
C
C     Elastic stiffness
C
      DO K1=1, NDI
        DO K2=1, NDI
          DDSDDE(K2, K1)=ELAM
        END DO
        DDSDDE(K1, K1)=EG2+ELAM
      END DO
      DO K1=NDI+1,NTENS-2
        DDSDDE(K1, K1)=EG
      END DO
      DDSDDE(5,5)=EG2
      DDSDDE(6,6)=EG2
C
C     Calculate predictor stress and elastic strains
C
      DO K1=1, NTENS
        DO K2=1, NTENS
          STRESS(K2)=STRESS(K2)+DDSDDE(K2, K1)*DSTRAN(K1)
        END DO
        EELAS(K1)=EELAS(K1)+DSTRAN(K1)
      END DO
C
      SHYDRO=(STRESS(1)+STRESS(2)+STRESS(3))/THREE
C
C     Calculate relative trial stress
C
      CALL KPROYECTOR(P)
      CALL KMMULT(P,NTENS,NTENS,STRESS,NTENS,1,SDEV)
C
      DO K1=1,NTENS
        STSREL(K1)=SDEV(K1,1)-XBACK(K1)
      END DO
C
C     Calculate equivalent Von Mises stress
C
      SMISES=(STSREL(1)-STSREL(2))**2
     1+(STSREL(2)-STSREL(3))**2+(STSREL(3)-STSREL(1))**2
      DO K1=NDI+1, NTENS-2
        SMISES=SMISES+SIX*STSREL(K1)**2
      END DO
      SMISES=SMISES+THREE*STSREL(5)**2+THREE*STSREL(6)**2
      SMISES=DSQRT(SMISES/TWO)
      FBAR=DSQRT(TWO/THREE)*SMISES
C
C     Get the yield stress from the specifid hardening function.
C
      CALL KUHARDNLIN(SYIEL0,SYIELDK,HMOD,SIG0,SIGSAT,HRDRATE,
     1                      EHARDI,EHARDK,EQPLAS)
C
C     Determine if actively yielding
C
      SYIELD=SYIEL0
      SYIELDK0=SYIELDK
      IF(FBAR.GT.(ONE+TOLER)*SYIEL0) THEN
C
C       Actively yielding-Perform local Newton iterations
C       to find consistncy parameter and equivalent plastic
C       strain
C
        DO K1=1, NTENS
          FLOW(K1)=STSREL(K1)/FBAR
        END DO
C
C       Starts iterations
C
        ITER=1
        GAM_PAR=ZERO
        IFLAG=0
        DO
          FGAM=-SYIELD+FBAR-(EG2*GAM_PAR
     &    +DSQRT(TWO/THREE)*(SYIELDK-SYIELDK0))
          TETA2=ONE+((EHARDI+EHARDK)/(EG3))
          FJAC=-EG2*TETA2
C
C         Updates
C
          GAM_PAR=GAM_PAR-FGAM/FJAC
          EQPLAS1=EQPLAS+DSQRT(TWO/THREE)*GAM_PAR
          CALL KUHARDNLIN(SYIELD,SYIELDK,HMOD,SIG0,SIGSAT,HRDRATE,
     1                    EHARDI,EHARDK,EQPLAS1)
C
          IF(ABS(FGAM/FJAC).LT.TOLER) THEN
            IFLAG=0
            GOTO 801
          ELSE
            IF(ITER.GT.MAXITER) THEN
              IFLAG=1
              GOTO 802
            END IF
          END IF
C
          ITER=ITER+1
        END DO
C
  801   CONTINUE
C
C       Local Newton algorithm converged
C       Update back stresses, plastic strains, equivalent plastic
C       strains, stresses
C
C       Normal components
C
        DO K1=1,NDI
          XBACK(K1)=XBACK(K1)+DSQRT(TWO/THREE)*(SYIELDK-SYIELDK0)
     &    *FLOW(K1)
          EPLAS(K1)=EPLAS(K1)+GAM_PAR*FLOW(K1)
          EELAS(K1)=EELAS(K1)-EPLAS(K1)
          STRESS(K1)=FLOW(K1)*SYIELD+XBACK(K1)+SHYDRO
        END DO
C
C       Shear components
C
        DO K1=NDI+1,NTENS-2
          XBACK(K1)=XBACK(K1)+DSQRT(TWO/THREE)*(SYIELDK-SYIELDK0)
     &    *FLOW(K1)
          EPLAS(K1)=EPLAS(K1)+TWO*GAM_PAR*FLOW(K1)
          EELAS(K1)=EELAS(K1)-EPLAS(K1)
          STRESS(K1)=FLOW(K1)*SYIELD+XBACK(K1)
        END DO
C
C       Couple Stress components
C
        DO K1=NTENS-1,NTENS
          XBACK(K1)=XBACK(K1)+DSQRT(TWO/THREE)*(SYIELDK-SYIELDK0)
     &    *FLOW(K1)
          EPLAS(K1)=EPLAS(K1)+TWO*GAM_PAR*FLOW(K1)
          EELAS(K1)=EELAS(K1)-EPLAS(K1)
          STRESS(K1)=FLOW(K1)*SYIELD+XBACK(K1)
        END DO
C
        EQPLAS=EQPLAS+DSQRT(TWO/THREE)*GAM_PAR
C
C       Formulate the Jacobian (Material tangent)
C
        CALL KCLEAR(DDSDDE,NTENS,NTENS)
C
C       Calculate effective properties 
C
        BETA1=ONE-EG2*GAM_PAR/FBAR
        C1=ONE+((EHARDI+EHARDK)/EG3)
        BETABAR=(ONE/C1)-(EG2*GAM_PAR/FBAR)
        EFFG=EG*BETA1
        EFFG2=TWO*EFFG
        EFFLAM=ELAM+(EG2/THREE)*(ONE-BETA1)
        EFFHRD=-EG2*BETABAR
C
C      Cauchy stress terms (4x4 Submtrix)
C
        DO K1=1, NDI
          DO K2=1, NDI
            DDSDDE(K2, K1)=EFFLAM
          END DO
          DDSDDE(K1, K1)=(EFFG2+EFFLAM)
        END DO
        DO K1=NDI+1, NTENS-2
          DDSDDE(K1, K1)=EFFG
        END DO
        DO K1=1, NTENS-2
          DO K2=1, NTENS-2
            DDSDDE(K2, K1)=DDSDDE(K2, K1)+EFFHRD
     &                     *FLOW(K2)*FLOW(K1)
          END DO
        END DO
C
C      Cauchy-Couple stress terms (4x2 Submtrix)
C
        DO K1=1, NTENS-2
          DO K2=NTENS-1,NTENS
            DDSDDE(K1, K2)=DDSDDE(K1,K2)-(EG2/C1)*FLOW(K1)*FLOW(K2)
          END DO
        END DO
C
C      Couple-Cauchy stress terms (2x4 Submtrix)
C
        DO K1=NTENS-1,NTENS
          DO K2=1,NTENS-2
            DDSDDE(K1, K2)=DDSDDE(K1,K2)-(EG2/C1)*FLOW(K1)*FLOW(K2)
          END DO
        END DO
C
C      Couple stress terms (2x2 Submtrix)
C
        DO K1=NTENS-1,NTENS
          DO K2=NTENS-1,NTENS
            DDSDDE(K2, K1)=DDSDDE(K2,K1)+EFFHRD
     &                     *FLOW(K2)*FLOW(K1)
          END DO
          DDSDDE(K1, K1)=DDSDDE(K1, K1)+EFFG2
        END DO
      END IF
C
C     Store updated state variables
C
      DO K1=1, NTENS
        STATEV(        K1)=EELAS(K1)
        STATEV(NTENS+K1  )=EPLAS(K1)
        STATEV(2*NTENS+K1)=XBACK(K1)
      END DO
      STATEV(3*NTENS+1)=EQPLAS
      STATEV(3*NTENS+2)=SMISES
C
  802 IF (IFLAG.EQ.1) THEN
         WRITE(*,*)
         WRITE(*,*) 'LOCAL PLASTICITY ALGORITHM DID NOT CONVREGED'
         WRITE(*,*) 'AT GAUSS POINT=',NPT, 'ELEMENT=',NOEL
         WRITE(*,*) 'AFTER=',ITER,' ITERATIONS'
         WRITE(*,*) 'LAST CORRECTION=',-fgam/fjac
         CALL XIT
      END IF
C
      RETURN
C
      END
C
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C     SUBROUTINE KUHARDNLIN                                            C
C                                                                      C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE KUHARDNLIN(SYIELDI,SYIELDK,HMOD,SIG0,SIGSAT,HRDRATE,
     1                      EHARDI,EHARDK,EQPLAS)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0)
C
C     Compute yield stress and hardening modulus corresponding to EQPLAS
C     Isotropic hardenining==> Nonlinear Chaboche model
C     Kinematic hardening  ==> Linear
C
      SYIELDI=DSQRT(TWO/THREE)*(SIG0+SIGSAT*(ONE-EXP(-HRDRATE*EQPLAS)))
      EHARDI=SIGSAT*HRDRATE*(EXP(-HRDRATE*EQPLAS))
      SYIELDK=HMOD*EQPLAS
      EHARDK=HMOD
C
      RETURN
C
      END
C
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C     SUBROUTINE PROYECTOR                                             C
C     PROJECTS THE STRESS TENSOR INTO THE DEVIATORIC SPACE             C
C     FOR A PLANE STRAIN PROBLEM                                       C
C                                                                      C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE KPROYECTOR(P)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0)
C
      DIMENSION P(6,6)
C
      CALL KCLEAR(P,6,6)
      P(1,1)=TWO/THREE
      P(1,2)=-ONE/THREE
      P(1,3)=-ONE/THREE
      P(2,1)=-ONE/THREE
      P(2,2)=TWO/THREE
      P(2,3)=-ONE/THREE
      P(3,1)=-ONE/THREE
      P(3,2)=-ONE/THREE
      P(3,3)=TWO/THREE
      P(4,4)=ONE
      P(5,5)=ONE
      P(6,6)=ONE
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C         --G E N E R A L  F E M  S U B R O U T I N E S--              C
C                                                                      C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE KSTDM:GENERATES THE STRAIN-DISPLACEMENT MATRIX B      C
C     AND JACOBIAN DETERMINANT DDET AT THE POINT r ,s                  C
C                                                 i  j                 C
C     FOR AN 9-NODED 2D ELEMENT-PLANE STRAIN                           C
C     B    =STRAIN-DISPLACEMENT MATRIX                                 C
C     BLM  =RELATIVE STRAIN DISPLACEMENT MATRIX                        C
C     DDET =JACOBIAN DETERMINANT                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE KSTDM(IDEL,NNE,NDOFEL,NTENS,XX,B,BLM,DDET,R,S,XBAR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ONE=1.D0,TWO=2.D0,THREE=3.D0,HALF=0.5D0,ZERO=0.D0)
C
      DIMENSION XX(2,NNE),B(NTENS,NDOFEL),BLM(1,NDOFEL),
     1P(2,NNE),AUX1(2,NNE),XJ(2,2),XJI(2,2),SP(NNE)
C
C     Initialize arrays
C
      XBAR=ONE
      CALL KCLEAR(B,NTENS,NDOFEL)
      CALL KCLEAR(XJ,2,2)
      CALL KCLEAR(XJI,2,2)
      CALL KCLEAR(P,2,NNE)
      CALL KCLEAR(SP,NNE,1)
C
C     Shape functions derivatives w.r.t natural coordinates
C
      CALL KSFDER(IELT,NDOFEL,NNE,R,S,P)
      CALL KSHAPEFUN(NNE,R,S,SP)
C
C     Computes the Jacobian operator and its determinant at point (r,s)
C
      CALL KJACOPER(NNE,XJ,XX,P)
      DDET=XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1)
C
C     Computes the inverse of the Jacobiam operator
C
      CALL KJACINVE(XJ,DDET,XJI)
C
C     Jacobian Inverse times Shape Functions derivatives w.r.t natural coordinates
C     to produce shape function derivatives w.r.t x,y coordinates.
C          
      CALL KMMULT(XJI,2,2,P,2,NNE,AUX1)
C
C     Assembles B matrix for a
C     Cosserat element.
C
      DO I=1,NNE
        II=3*(I-1)+1
        B(1,II)=AUX1(1,I)
        B(2,II+1)=AUX1(2,I)
        B(4,II)=AUX1(2,I)
        B(4,II+1)=AUX1(1,I)
        B(5,II+2)=AUX1(2,I)
        B(6,II+2)=AUX1(1,I)
      END DO
C
C     Assembles (L-M) matrix to compute
C     relative strain vector.
C
      DO I=1,NNE
        II=3*(I-1)+1
        BLM(1,II)=HALF*AUX1(2,I)
        BLM(1,II+1)=-HALF*AUX1(1,I)
        BLM(1,II+2)=HALF*SP(I)
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE KSFDER:GENERATES THE SHAPE FUNCTION DERIVATIVES       C
C     ACCORDING TO ELEMENT TYPE AT THE POINT r ,s                      C
C                                             i  j                     C
C     B    =STRAIN-DISPLACEMENT MATRIX                                 C
C     DDET =JACOBIAN DETERMINANT                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KSFDER(IELT,NDOFEL,NNE,R,S,P)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(ONE=1.D0,TWO=2.D0, HALF=0.5D0,QUARTER=0.25D0,FOUR=4.D0)
C
      DIMENSION P(2,NNE)
C
      RP= ONE+R
      SP= ONE+S
      RM= ONE-R
      SM= ONE-S
      RMS=ONE-R**TWO
      SMS=ONE-S**TWO
C
C     8-NODED ELEMENT
C
C     Derivatives w.r.t the natural coordinates
C     w.r.t.r
C
      P(1,8)=-HALF*SMS
      P(1,7)=-R*SP
      P(1,6)= HALF*SMS
      P(1,5)=-R*SM
      P(1,4)=-QUARTER*SP-HALF*P(1,7)-HALF*P(1,8)
      P(1,3)= QUARTER*SP-HALF*P(1,7)-HALF*P(1,6)
      P(1,2)= QUARTER*SM-HALF*P(1,5)-HALF*P(1,6)
      P(1,1)=-QUARTER*SM-HALF*P(1,8)-HALF*P(1,5)
C
C        w.r.t.s
C
      P(2,8)=-S*RM
      P(2,7)= HALF*RMS
      P(2,6)=-S*RP
      P(2,5)=-HALF*RMS
      P(2,4)= QUARTER*RM-HALF*P(2,7)-HALF*P(2,8)
      P(2,3)= QUARTER*RP-HALF*P(2,7)-HALF*P(2,6)
      P(2,2)=-QUARTER*RP-HALF*P(2,5)-HALF*P(2,6)
      P(2,1)=-QUARTER*RM-HALF*P(2,8)-HALF*P(2,5)
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE KSHAPEFUN:GENERATES THE SHAPE FUNCTION DERIVATIVES    C
C     ACCORDING TO ELEMENT TYPE AT THE POINT r ,s                      C
C                                             i  j                     C
C     B    =STRAIN-DISPLACEMENT MATRIX                                 C
C     DDET =JACOBIAN DETERMINANT                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KSHAPEFUN(NNE,R,S,SF)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(ONE=1.D0,TWO=2.D0, HALF=0.5D0,QUARTER=0.25D0,FOUR=4.D0)
C
      DIMENSION SF(NNE)
C
      RP= ONE+R
      SP= ONE+S
      RM= ONE-R
      SM= ONE-S
      RMS=ONE-R**TWO
      SMS=ONE-S**TWO
C
C     8-NODED ELEMENT
C
      SF(8)=HALF*SMS*RM
      SF(7)=HALF*RMS*SP
      SF(6)=HALF*SMS*RP
      SF(5)=HALF*RMS*SM
      SF(4)=QUARTER*RM*SP-HALF*SF(7)-HALF*SF(8)
      SF(3)=QUARTER*RP*SP-HALF*SF(7)-HALF*SF(6)
      SF(2)=QUARTER*RP*SP-HALF*SF(5)-HALF*SF(6)
      SF(1)=QUARTER*RM*SP-HALF*SF(8)-HALF*SF(5)
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C 3X3 GAUSS POINTS GENERATION                                          C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE KGPOINTS3X3(XP,XW)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ZERO=0.D0,ONE=1.D0)
C
      DIMENSION XP(2,9),XW(9),RLGP(2,9)
C
      RLGP(1,1)=-ONE
      RLGP(1,2)=ZERO
      RLGP(1,3)= ONE
      RLGP(1,4)=-ONE
      RLGP(1,5)=ZERO
      RLGP(1,6)= ONE
      RLGP(1,7)=-ONE
      RLGP(1,8)=ZERO
      RLGP(1,9)= ONE
C
      RLGP(2,1)=-ONE
      RLGP(2,2)=-ONE
      RLGP(2,3)=-ONE
      RLGP(2,4)=ZERO
      RLGP(2,5)=ZERO
      RLGP(2,6)=ZERO
      RLGP(2,7)= ONE
      RLGP(2,8)= ONE
      RLGP(2,9)= ONE
C
      XW(1)=0.555555555555556D0**2
      XW(2)=0.555555555555556D0*0.888888888888889D0
      XW(3)=0.555555555555556D0**2
C
      XW(4)=0.555555555555556D0*0.888888888888889D0
      XW(5)=0.888888888888889D0**2
      XW(6)=0.555555555555556D0*0.888888888888889D0
C
      XW(7)=0.555555555555556D0**2
      XW(8)=0.555555555555556D0*0.888888888888889D0
      XW(9)=0.555555555555556D0**2
C
      G=DSQRT(0.60D0)
      DO I=1,2
        DO J=1,9 
          XP(I,J)=G*RLGP(I,J)
        END DO
      END DO
C     
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   SUBROUTINE KGPOINTS2X2                                             C
C                                                                      C
C   2X2 GAUSS POINTS GENERATION                                        C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KGPOINTS2X2(XP,XW)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ZERO=0.D0,ONE=1.D0)
C
      DIMENSION XP(2,4),XW(4)
C
      CALL KCLEARV(XW,4)
      CALL KCLEAR(XP,2,4)
C
      DO K1=1,4
        XW(K1)=ONE
      END DO
C
      XP(1,1)=-0.577350269189626
      XP(2,1)=-0.577350269189626
      XP(1,2)= 0.577350269189626
      XP(2,2)=-0.577350269189626
      XP(1,3)=-0.577350269189626
      XP(2,3)= 0.577350269189626
      XP(1,4)= 0.577350269189626
      XP(2,4)= 0.577350269189626
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   SUBROUTINE KJACOPER                                                C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KJACOPER(NNE,XJA,XCORD,RDER)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(ZERO=0.0D0)
C
      DIMENSION XJA(2,2),XCORD(2,NNE),RDER(2,NNE)
C
      CALL KCLEAR(XJA,2,2)
C
      DUM=ZERO
      DO K=1,2
        DO J=1,2
          DO I=1,NNE
            DUM=DUM+RDER(K,I)*XCORD(J,I)
          END DO
          XJA(K,J)=DUM
          DUM=ZERO
        END DO
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      2X2 Jacobian inverse                                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KJACINVE(XJA,DD,XJAI)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION XJA(2,2),XJAI(2,2),XADJ(2,2)
C
      XADJ(1,1)=XJA(2,2)
      XADJ(1,2)=XJA(1,2)
      XADJ(2,1)=XJA(2,1)
      XADJ(2,2)=XJA(1,1)
      DO J=1,2
        DO K=1,2
          COFA=((-1)**(J+K))*XADJ(J,K)
          XJAI(J,K)=COFA/DD
        END DO
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C             M A T R I X   H A N D L I N G                            C
C-------------U T I L I T I E S   B L O C K--------------              C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE KCLEAR(A,N,M)                                        C
C      Clear a real matrix                                             C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE KCLEAR(A,N,M)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(ZERO=0.0D0)
      DIMENSION A(N,M)
C
      DO I=1,N
        DO J=1,M
          A(I,J)=ZERO
        END DO
      END DO
C
      RETURN
C
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE KMMULT(A,NRA,NCA,B,NRB,NCB,C)                        C
C      Real matrix product                                             C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE KMMULT(A,NRA,NCA,B,NRB,NCB,C)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(ZERO=0.D0)
      DIMENSION A(NRA,NCA),B(NRB,NCB),C(NRA,NCB)
C
      CALL KCLEAR(C,NRA,NCB)
      DUM=ZERO
      DO I=1,NRA
        DO J=1,NCB
         DO K=1,NCA
           DUM=DUM+A(I,K)*B(K,J)
          END DO
          C(I,J)=DUM
          DUM=ZERO
        END DO
      END DO
C
      RETURN
C
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE KSMULT(A,NR,NC,S)                                    C
C      Matrix times a scalar.                                          C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE KSMULT(A,NR,NC,S)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(NR,NC)
C
      DO I=1,NR
        DO J=1,NC
          DUM=A(I,J)
          A(I,J)=S*DUM
          DUM=0.D0
        END DO  
      END DO
C
      RETURN
C
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE KUPDMAT(A,NR,NC,B)                                   C
C      Updates an existing matrix with an incremental matrix.          C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE KUPDMAT(A,NR,NC,B)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(ZERO=0.D0)
C
      DIMENSION A(NR,NC),B(NR,NC)
C
      DO I=1,NR
        DO J=1,NC
          DUM=A(I,J)
          A(I,J)=ZERO
          A(I,J)=DUM+B(I,J)
          DUM=ZERO
        END DO
      END DO
C
      RETURN
C
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE KMTRAN(A,NRA,NCA,B)                                  C      
C      Matrix transpose                                                C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE KMTRAN(A,NRA,NCA,B)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(NRA,NCA),B(NCA,NRA)
C
      CALL KCLEAR(B,NCA,NRA)
      DO I=1,NRA
       DO J=1,NCA
         B(J,I)=A(I,J)
        END DO
      END DO
C
      RETURN
C
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE KMAVEC(A,NRA,NCA,B,C)                                C
C      Real matrix times vector                                        C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE KMAVEC(A,NRA,NCA,B,C)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(ZERO=0.D0)
      DIMENSION A(NRA,NCA),B(NCA),C(NRA)
C
      CALL KCLEARV(C,NRA)
C
      DO K1=1,NRA
        DO K2=1,NCA
          C(K1)=C(K1)+A(K1,K2)*B(K2)	    
        END DO
      END DO     
C
      RETURN
C
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE KCLEARV(A,N)                                         C
C      Clear a real vector                                             C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE KCLEARV(A,N)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(ZERO=0.0D0)
C
      DIMENSION A(N)
C
      DO I=1,N
        A(I)=ZERO
      END DO
C
      RETURN
C
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE KUPDVEC(A,NR,B)                                      C
C      Updates an existing vector with an incremental vector.          C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE KUPDVEC(A,NR,B)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(ZERO=0.D0)
C
      DIMENSION A(NR),B(NR)
C
      DO I=1,NR
        DUM=A(I)
        A(I)=ZERO
        A(I)=DUM+B(I)
        DUM=ZERO
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE KVECSUB(A,NRA,B,NRB,C)                               C
C      Substracts one column vector from another column vector         C
C      IFLAG=0 for substraction                                        C
C      IFLAG=1 for addition                                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE KVECSUB(A,NRA,B,NRB,C,IFLAG)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (ONE=1.0D0, ONENEG=-1.0D0)
C
      DIMENSION A(NRA,1),B(NRB,1),C(NRB,1)
C
      SCALAR=ONENEG
C
      IF (IFLAG.EQ.1) SCALAR=ONE
C
      DO I=1,NRA
        C(I,1)=A(I,1)+B(I,1)*SCALAR
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE KMATSUB(A,NRA,NCA,B,C,IFLAG)                         C
C      Substracts one rectangular matrix from another rectangular      C
C      matrix                                                          C
C      IFLAG=0 for substraction                                        C
C      IFLAG=1 for addition                                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE KMATSUB(A,NRA,NCA,B,C,IFLAG)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (ONE=1.0D0, ONENEG=-1.0D0)
C
      DIMENSION A(NRA,NCA),B(NRA,NCA),C(NRA,NCA)
C
      CALL KCLEAR(C,NRA,NCA)
C
      SCALAR=ONENEG
C
      IF (IFLAG.EQ.1) SCALAR=ONE
C
      DO I=1,NRA
        DO J=1,NCA
          C(I,J)=A(I,J)+B(I,J)*SCALAR
        END DO
      END DO
C
      RETURN
C
      END
C
C
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C     SUBROUTINE IDENTITY                                              C
C     CREATES AN IDENTITY MATRIX OF DIMENSIONS NDIM,NDIM               C
C                                                                      C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE KIDENTITY(DEL,NDIM)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(ONE=1.D0)
C
      DIMENSION DEL(NDIM,NDIM)
C
      CALL KCLEAR(DEL,NDIM,NDIM)
C
      DO K1=1,NDIM
        DEL(K1,K1)=ONE
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   SUBROUTINE KINVERSE                                                 C
C                                                                      C
C   IVEERSE OF A MATRIX USING LU DECOMPOSITION                         C
C   TAKEN FROM NUMERICAL RECIPES By Press et al                        C
C                                                                      C
C   A   Matrix to be inverted.                                         C
C   Y   Inverse of A                                                   C
C   N   Dimension                                                      C
C                                                                      C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KINVERSE(A,Y,NP,N)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (ZERO=0.D0,ONE=1.D0)
C
      DIMENSION A(NP,NP),Y(NP,NP),INDX(NP),AUX(NP,NP)
C
      CALL KCLEAR(AUX,NP,NP)
      CALL KCOPYMAT(A,AUX,N)
C
      DO I=1,N
        DO J=1,N
          Y(I,J)=ZERO
        END DO
        Y(I,I)=ONE
      END DO
      CALL KLUDCMP(AUX,N,NP,INDX,D)
      DO J=1,N
        CALL KLUBKSB(AUX,N,NP,INDX,Y(1,J))
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   SUBROUTINE KLUDCMP                                                 C
C                                                                      C
C   LU MATRIX DECOMPOSITION                                            C
C   TAKEN FROM NUMERICAL RECIPES By Press et al                        C
C                                                                      C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KLUDCMP(A,N,NP,INDX,D)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(NMAX=500,TINY=1.0E-20,ZERO=0.D0,ONE=1.D0)
C
      DIMENSION INDX(N),A(NP,NP),VV(NMAX)
C
      D=ONE
      DO I=1,N
        AAMAX=ZERO
        DO J=1,N
          IF(ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
        END DO
        IF(AAMAX.EQ.0.) PAUSE 'SINGULAR MATRIX IN LUDCMP'
        VV(I)=ONE/AAMAX
      END DO
C
      DO J=1,N
        DO I=1,J-1
          SUM=A(I,J)
          DO K=1,I-1
            SUM=SUM-A(I,K)*A(K,J)
          END DO
          A(I,J)=SUM
        END DO
        AAMAX=ZERO
        DO I=J,N
          SUM=A(I,J)
          DO K=1,J-1
            SUM=SUM-A(I,K)*A(K,J)
          END DO
          A(I,J)=SUM
          DUM=VV(I)*ABS(SUM)
          IF(DUM.GE.AAMAX) THEN
             IMAX=I
             AAMAX=DUM
          END IF
        END DO
        IF(J.NE.IMAX) THEN
           DO K=1,N
             DUM=A(IMAX,K)
             A(IMAX,K)=A(J,K)
             A(J,K)=DUM
           END DO
           D=-D
           VV(IMAX)=-VV(J)
        END IF
        INDX(J)=IMAX
        IF(A(J,J).EQ.0.) A(J,J)=TINY
        IF(J.NE.N) THEN
           DUM=ONE/A(J,J)
           DO I=J+1,N
             A(I,J)=A(I,J)*DUM
           END DO
        END IF
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   SUBROUTINE KLUBKSB                                                 C
C                                                                      C
C   FORWARD SUBSTITUTION                                               C
C   TAKEN FROM NUMERICAL RECIPES By Press et al                        C
C                                                                      C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE KLUBKSB(A,N,NP,INDX,B)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (ZERO=0.D0)
C
      DIMENSION INDX(N),A(NP,NP),B(NP)
C
      II=0
      DO I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF(II.NE.0) THEN
           DO J=II,I-1
             SUM=SUM-A(I,J)*B(J)
           END DO
        ELSE IF(SUM.NE.ZERO) THEN
           II=I
        END IF
        B(I)=SUM
      END DO
C
      DO I=N,1,-1
        SUM=B(I)
        DO J=I+1,N
          SUM=SUM-A(I,J)*B(J)
        END DO
        B(I)=SUM/A(I,I)
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C     SUBROUTINE KCOPYMAT                                              C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE KCOPYMAT(A,B,N)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION A(N,N),B(N,N)
C
      CALL KCLEAR(B,N,N)
C
      DO K1=1,N
        DO K2=1,N
          B(K1,K2)=A(K1,K2)
        END DO
      END DO
C
      RETURN
C
      END