C �ƶ���Դ�ӳ���DFLUX 
C ȫ�����ù��ʵ�λ�Ƶ�λ
C ���� ���� 300; ת�� 400
      SUBROUTINE DFLUX(FLUX,SOL,KSTEP,KINC,TIME,NOEL,NPT,COORDS,
     & JLTYP,TEMP,PRESS,SNAME)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION FLUX(2), TIME(2), COORDS(3)
      CHARACTER*80 SNAME
	  real x1,y1,z1,x2,t1,t2,v,r1,r2,h,p,r,m,e,tr,k,xn,rn,x,y,z
      REAL PI; PARAMETER ( PI= 3.14159 )
C ----------------------------------------------------------
      x=COORDS(1)
      y=COORDS(2)
      z=COORDS(3)
C ���ӳ�����x������Ϊλ�Ʒ���	  
C ���ӿ�ʼ����ͷλ�� m
      x1=0.005
      y1=0.0 
	  z1=0.003            !�۳ر������꣨����ԴΪԲ���Σ�
C ���ӽ�������ͷλ�� m
      x2=0.025
C ���ӿ�ʼ����ͷͣ��ʱ�� s
      t1=3.0     
C ���ӽ�������ͷͣ��ʱ�� s
      t2=1.0 
C �����ٶ� m/s
      v=0.005
C ���뾶 m�� (����ʵ������޸ģ��˴�������ͷ��ΪԲ��������׶�ȣ�
      r1=0.005 
C ������뾶 m
      r2=0.00125 
C �����볤�� m
      h=0.003 
C �������ѹǿ Pa-��ѹ��Ϊ12KN
      p=152866242.0 
C ����ͷת�� rad/s
      r= 600*2*PI/60 
C Ħ��ϵ��
	  m=0.37 
C ��Ч��
	  e=0.98 
C �����������������֮��:k 
      k=0.33
C ��ǰ����ͷ����Xλ��:xn   
C ��ǰ����ͷ���ĵ����ֵ�ľ���:rn 
C ����Ħ���� pa  :Tr

C ȷ��Ħ������solΪ�ڵ��¶�ֵ����������Ӧ�����¶�����ͼȷ��
	IF (SOL .LT. 390 ) THEN
	Tr=m*p
	ELSE IF(SOL .LT. 450) THEN
	Tr=0.577*(98.1-(98.1-47)*(SOL-390)/(450-390))*1e6
	ELSE IF(SOL .LT. 550) THEN
	Tr=0.577*(47-(47-5)*(SOL-450)/(550-450))*1e6
	ELSE
	Tr=0.577*1*1e6
	END IF
    
	
C ���Ƚ�Flux��1������
	FLUX(1)=0.0
	
C ȷ������ͷ����Xλ��
      IF ( TIME(2).LE.t1 ) THEN
        xn=x1                  !��ǰͣ��
      ELSE IF ( TIME(2).GT.(t1+((x2-x1)/v)) ) THEN
        xn=x2                  !����ͣ��
      ELSE
        xn=x1+((TIME(2)-t1)*v)
      END IF

C ȷ������ͷ���ĵ���ǰ���ֵ�ľ���
      
	  rn=SQRT((x-xn)**2+(y-y1)**2)

C �жϻ��ֵ�����������
      IF ((JLTYP.EQ.0).AND.(rn.LE.r1))THEN
		       FLUX(1)=Tr*r*rn*e
	ELSE IF (JLTYP.EQ.1 
     &		.AND. COORDS(3).GE.(z1-h)
     &		.AND. rn.LE.r2 ) THEN
		FLUX(1)=(2*k*Tr*r*e*(r1**3-r2**3))/(3*h*r2**2)
		
	  ELSE
	       FLUX(1)=0.0
	  END IF
	
	  RETURN 
      END     
C �ƶ���Դ�ӳ������