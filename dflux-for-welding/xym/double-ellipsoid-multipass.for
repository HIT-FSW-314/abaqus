      SUBROUTINE DFLUX(FLUX,SOL,JSTEP,JINC,TIME,NOEL,NPT,COORDS,JLTYP,
     1                 TEMP,PRESS,SNAME)
C
      INCLUDE 'ABA_PARAM.INC'


      DIMENSION COORDS(3),FLUX(2),TIME(2)
      CHARACTER*80 SNAME
      x=COORDS(1)
      y=COORDS(2)
      z=COORDS(3)
	  PI=3.1415926	 
	  
C     第一条焊缝焊接参数	
  
      wu1=110.0
      wi1=8.0
      effi1=0.8
      v1=0.005
      q1=wu1*wi1*effi1
	  t1=0                 !焊接起始时间（相对总时间）
      d1=v1*(TIME(2)-t1)

      x1=0.0               !焊接起始时间X坐标
      y1=0.0               !焊接起始时间X坐标
      z1=0.002             !焊接起始时间X坐标
	  
      af1=0.003            !熔池前半轴长度
	  b1=0.003             !熔池宽度
	  c1=0.003             !熔池深度
	  ar1=0.006            !熔池后半轴长度
	  ff1=0.7              !熔池前半区热流量分配比例

C     第一二条焊缝焊接层间冷却时间	   
	  t12=10

C     第二条焊缝焊接参数	
  
      wu2=110.0
      wi2=8.0
      effi2=0.8
      v2=0.005
      q2=wu2*wi2*effi2
	  t2=17                 !焊接起始时间（相对总时间）
      d2=v2*(TIME(2)-t2)

      x2=0.0               !焊接起始时间X坐标
      y2=0.0               !焊接起始时间X坐标
      z2=0.003             !焊接起始时间X坐标
	  
      af2=0.003            !熔池前半轴长度
	  b2=0.003             !熔池宽度
	  c2=0.003             !熔池深度
	  ar2=0.006            !熔池后半轴长度
	  ff2=0.7              !熔池前半区热流量分配比例	  
	
C     第一道焊缝焊接
      if(TIME(2).lt.(t2-t12)) then   
        q=q1
		d=d1
		x0=x1
		y0=y1
		z0=z1
		af=af1
		b=b1
		c=c1
		ar=ar1
		ff=ff1

C     第一二人道焊缝层间冷区
		
      else if (TIME(2).lt.t2.and.TIME(2).ge.(t2-t12)) then
        q=0        !冷却时热输入为0
		d=d1
		x0=x1
		y0=y1
		z0=z1
		af=af1
		b=b1
		c=c1
		ar=ar1
		ff=ff1	

C     第二道焊缝焊接
		
      else if (TIME(2).ge.t2) then
        q=q2
		d=d2
		x0=x2
		y0=y2
		z0=z2
		af=af2
		b=b2
		c=c2
		ar=ar2
		ff=ff2

        end if		

C     JLTYP＝1，表示为体热源
      JLTYP=1
	  
	  heat1=6.0*sqrt(3.0)*q/(af*b*c*PI*sqrt(PI))*ff
	  heat2=6.0*sqrt(3.0)*q/(ar*b*c*PI*sqrt(PI))*(2.0-ff)
	  
	  shape1=exp(-3.0*(x-x0-d)**2/af**2-3.0*(y-y0)**2/b**2
     $	-3.0*(z-z0)**2/c**2)
      shape2=exp(-3.0*(x-x0-d)**2/ar**2-3.0*(y-y0)**2/b**2
     $	-3.0*(z-z0)**2/c**2)
	 
	  IF(x .GE.(x0+d)) THEN
	  FLUX(1)=heat1*shape1
	  ELSE
	  FLUX(1)=heat2*shape2
	  ENDIF
      RETURN
      END