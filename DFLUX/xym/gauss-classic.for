c Version6.0 
      SUBROUTINE DFLUX(FLUX,SOL,JSTEP,JINC,TIME,NOEL,NPT,COORDS,JLTYP,
     1                 TEMP,PRESS,SNAME)
C
      INCLUDE 'ABA_PARAM.INC'


      DIMENSION COORDS(3),FLUX(2),TIME(2)
      CHARACTER*80 SNAME

      real a,Qs,rs,pi,Qv,H,b,rv,Q,aa,x0,y0,z0
c     Qs为面热源功率, a为面热源能量集中系数，rs为面热源作用范围
c     Qv为体热源功率，H为体热源深度，b为体热源能量衰减系数，rv为体热源有效作用半径
c     Q为热源功率，aa为热源有效吸收系数,(x0,y0,z0)为当前热源中心位置
	Q=1200000;	aa=0.99;	Qs=Q*aa*0.2;	Qv=Q-Qs
	a=0.3;	rs=2
	H=2.8;	b=0.15;	rv=1
	x0=0;	y0=5*TIME(1);	z0=0
	pi=3.14

C     JLTYP＝1，表示为体热源

      IF(JLTYP.EQ.1) THEN
          FLUX(1)=(6*Qv*(H-b*(z0-COORDS(3)))/(pi*rv*rv*H*H*(2-b)))
     $		*exp((-3)*sqrt((COORDS(1)-x0)**2+
     $		(COORDS(2)-y0)**2)/(rv*rv))
      END IF
C     JLTYP＝0，表示为面热源      
      IF(JLTYP.EQ.0) THEN
          FLUX(1)=(a*Qs/(pi*rs*rs))*exp(-1*a*((COORDS(1)-x0)**2+
     $		(COORDS(2)-y0)**2)/(rs**2))
      ENDIF
      
      RETURN
      END  
