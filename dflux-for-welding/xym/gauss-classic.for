c Version6.0 
      SUBROUTINE DFLUX(FLUX,SOL,JSTEP,JINC,TIME,NOEL,NPT,COORDS,JLTYP,
     1                 TEMP,PRESS,SNAME)
C
      INCLUDE 'ABA_PARAM.INC'


      DIMENSION COORDS(3),FLUX(2),TIME(2)
      CHARACTER*80 SNAME

      real a,Qs,rs,pi,Qv,H,b,rv,Q,aa,x0,y0,z0
c     QsΪ����Դ����, aΪ����Դ��������ϵ����rsΪ����Դ���÷�Χ
c     QvΪ����Դ���ʣ�HΪ����Դ��ȣ�bΪ����Դ����˥��ϵ����rvΪ����Դ��Ч���ð뾶
c     QΪ��Դ���ʣ�aaΪ��Դ��Ч����ϵ��,(x0,y0,z0)Ϊ��ǰ��Դ����λ��
	Q=1200000;	aa=0.99;	Qs=Q*aa*0.2;	Qv=Q-Qs
	a=0.3;	rs=2
	H=2.8;	b=0.15;	rv=1
	x0=0;	y0=5*TIME(1);	z0=0
	pi=3.14

C     JLTYP��1����ʾΪ����Դ

      IF(JLTYP.EQ.1) THEN
          FLUX(1)=(6*Qv*(H-b*(z0-COORDS(3)))/(pi*rv*rv*H*H*(2-b)))
     $		*exp((-3)*sqrt((COORDS(1)-x0)**2+
     $		(COORDS(2)-y0)**2)/(rv*rv))
      END IF
C     JLTYP��0����ʾΪ����Դ      
      IF(JLTYP.EQ.0) THEN
          FLUX(1)=(a*Qs/(pi*rs*rs))*exp(-1*a*((COORDS(1)-x0)**2+
     $		(COORDS(2)-y0)**2)/(rs**2))
      ENDIF
      
      RETURN
      END  
