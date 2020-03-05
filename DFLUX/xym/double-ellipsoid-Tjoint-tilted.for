C     文件头，定义函数参数  
      SUBROUTINE DFLUX(FLUX,SOL,JSTEP,JINC,TIME,NOEL,NPT,COORDS,JLTYP,
     1                 TEMP,PRESS,SNAME)
C
      INCLUDE 'ABA_PARAM.INC'

C     定义基本数组，坐标，热通量值，时间
      DIMENSION COORDS(3),FLUX(2),TIME(2)
      CHARACTER*80 SNAME
C     定义热源参数，注意单位统一
      u=380.0    !电压
      i=2.0      !电流
      effi=0.8   !焊接效率
      v=0.005    !焊接速度
      q=u*i*effi     !热输入大小
      d=v*TIME(2)    !热源中心移动距离
C     定义坐标
      x=COORDS(1)
      y=COORDS(2)
      z=COORDS(3)
C     定义焊接起点
      x0=-0.0035
      y0=0.001
      z0=0.03
C       定义热源参数
        a=0.002
	  b=0.002
	  c=0.002
	  aa=0.004     !双椭球各方向半径 
	  f1=1.0       !前半球份额
	  PI=3.1415926
	 
      beta=45.0    !旋转角度，本例为绕z轴旋转,逆时针为正，z轴负方向为前进方向
	  betahudu=beta/180.0*PI  !将角度转换为弧度制
	  
	  heat1=6.0*sqrt(3.0)*q/(a*b*c*PI*sqrt(PI))*f1         !前半部分热输入
	  heat2=6.0*sqrt(3.0)*q/(aa*b*c*PI*sqrt(PI))*(2.0-f1)  !后半部分热输入
C       T型接头的高斯坐标xy平面需要沿z轴逆时针旋转45度，需要进行相应的坐标旋转变换处理
C	    此时x'=x*cos(betahudu)+y*sin(betahudu) ,
C       y'=Y*cos(betahudu)-X*sin(betahudu).将上述两式带入形状方程 
	  shape1=exp(-3.0*((x-x0)*cos(betahudu)+(y-y0)*sin(betahudu))**2/a**2-3.0*((y-y0)*cos(betahudu)-
     $     (x-x0)*sin(betahudu))**2/b**2-3.0*((z-z0+d))**2/c**2)
      shape2=exp(-3.0*((x-x0)*cos(betahudu)+(y-y0)*sin(betahudu))**2/aa**2-3.0*((y-y0)*cos(betahudu)-
     $     (x-x0)*sin(betahudu))**2/b**2-3.0*((z-z0+d))**2/c**2)

C     JLTYP＝1，表示为体热源
      JLTYP=1
	  IF(z .GE.(z0-d)) THEN      !热源移动方向为z轴负方向
	  FLUX(1)=heat1*shape1       !前半部分热通量
	  ELSE
	  FLUX(1)=heat2*shape2       !后半部分热通量
	  ENDIF
      RETURN
      END