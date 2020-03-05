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
      v=0.005
      q=wu*wi*effi
	  
      !焊接热源形状参数(球热源)
      a=0.003
	  b=0.003
	  c=0.003


      !各个关键点坐标及焊枪经过时间(相对总时间累计shijian)
      x1=0
      y1=-0.003
      z1=0.003
	  t1=0.0
      
      x2=0.012
      y2=-0.003
      z2=0.003
	  t2=2.4
      
      x3=0.015
      y3=0
      z3=0.003
	  t3=3.3425
      
      x4=0.018
      y4=0.003
      z4=0.003
	  t4=4.285
      

      !焊缝几何形状参数
      r=0.003

      !热源分布
	  heat=6.0*sqrt(3.0)*q/(a*b*c*PI*sqrt(PI))
        
        !第一段直线
        if(TIME(2).lt.t2) then      
          d=v*TIME(2)
        !热源中心在总体坐标系的坐标
          xd=x1+d
          yd=y1
          zd=z1
        !热源局部坐标
          xx=x-xd
          yy=y-yd
          zz=z-zd
          
	   shape=exp(-3.0*xx**2/a**2-3.0*yy**2/b**2-3.0*zz**2/c**2)
	    

        !第一段圆弧
        else if (TIME(2). lt. t3. and. TIME(2). ge. t2)then
          d=v*(TIME(2)-t2) 
          !焊缝圆弧扫过的圆心角
          theta=d/r
          !热源中心在总体坐标系中的坐标
          xd=x2+r*sin(theta)
		  yd=y2+r*(1-cos(theta))
          zd=z2
          
          xx=x-xd
          yy=y-yd
          zz=z-zd
          
          shape=exp(-3.0*xx**2/a**2-3.0*yy**2/b**2-3.0*zz**2/c**2)
        
         !第二段圆弧
        else if (TIME(2). lt. t4. and. TIME(2). ge. t3)then
          d=v*(TIME(2)-t3) 
          !焊缝圆弧扫过的圆心角
          theta=d/r
          !热源中心在总体坐标系中的坐标
          xd=x3+r*(1-cos(theta))
		  yd=y3+r*sin(theta)
          zd=z3
          
          xx=x-xd
          yy=y-yd
          zz=z-zd
          
          shape=exp(-3.0*xx**2/a**2-3.0*yy**2/b**2-3.0*zz**2/c**2)
                 
          
        else if (TIME(2).ge.t4) then
          d=v*(TIME(2)-t4)
          xd=x4+d
          yd=y4
          zd=z4
          
          xx=x-xd
          yy=y-yd
          zz=z-zd
          
          shape=exp(-3.0*xx**2/a**2-3.0*yy**2/b**2-3.0*zz**2/c**2)
          
       
        end if
            
C     JLTYP＝1，表示为体热源
      JLTYP=1                                                                                                                
      FLUX(1)=heat*shape
      RETURN
      END