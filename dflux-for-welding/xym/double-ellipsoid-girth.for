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
	  
	  	  
      !焊接工艺参数，使用同一热源
      wu=110.0
      wi=10.0
      effi=0.8
      v=0.010    !圆周方向速度
	  v1=0.005   !纵向方向速度
      q=wu*wi*effi
	  
      !焊接热源形状参数(球热源)
      a=0.003
	  b=0.003
	  c=0.003


      !起点坐标
      x1=0
      y1=-0.015
      z1=0.06
     
      

      !焊缝几何形状参数
      r=0.015

      !热源分布
	  heat=6.0*sqrt(3.0)*q/(a*b*c*PI*sqrt(PI))
        
          

          d=v*TIME(2)   !圆弧扫过的圆心角
		  d1=v1*TIME(2)      !纵向方向位移
          !焊缝圆弧扫过的圆心角
          theta=d/r
          !热源中心在总体坐标系中的坐标
          xd=x1+r*sin(theta)
		  yd=y1+r*(1-cos(theta))
          zd=z1-d1
          
          xx=x-xd
          yy=y-yd
          zz=z-zd
          
          shape=exp(-3.0*xx**2/a**2-3.0*yy**2/b**2-3.0*zz**2/c**2)
       
            
C     JLTYP＝1，表示为体热源
      JLTYP=1                                                                                                                
      FLUX(1)=heat*shape
      RETURN
      END