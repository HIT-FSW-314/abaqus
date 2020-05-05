C 移动热源子程序DFLUX 
C 全部采用国际单位制单位
C 参数 焊速 300; 转速 400
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
C 本子程序以x正方向为位移方向	  
C 焊接开始搅拌头位置 m
      x1=0.005
      y1=0.0 
	  z1=0.003            !熔池表面坐标（体热源为圆柱形）
C 焊接结束搅拌头位置 m
      x2=0.025
C 焊接开始搅拌头停留时间 s
      t1=3.0     
C 焊接结束搅拌头停留时间 s
      t2=1.0 
C 焊接速度 m/s
      v=0.005
C 轴肩半径 m， (依据实际情况修改，此处讲搅拌头简化为圆柱，忽略锥度）
      r1=0.005 
C 搅拌针半径 m
      r2=0.00125 
C 搅拌针长度 m
      h=0.003 
C 轴肩向下压强 Pa-下压力为12KN
      p=152866242.0 
C 搅拌头转速 rad/s
      r= 600*2*PI/60 
C 摩擦系数
	  m=0.37 
C 热效率
	  e=0.98 
C 搅拌针产热与轴肩产热之比:k 
      k=0.33
C 当前搅拌头中心X位置:xn   
C 当前搅拌头中心到积分点的距离:rn 
C 定义摩擦力 pa  :Tr

C 确定摩擦力，sol为节点温度值。根据流变应力与温度曲线图确定
	IF (SOL .LT. 390 ) THEN
	Tr=m*p
	ELSE IF(SOL .LT. 450) THEN
	Tr=0.577*(98.1-(98.1-47)*(SOL-390)/(450-390))*1e6
	ELSE IF(SOL .LT. 550) THEN
	Tr=0.577*(47-(47-5)*(SOL-450)/(550-450))*1e6
	ELSE
	Tr=0.577*1*1e6
	END IF
    
	
C 首先将Flux（1）清零
	FLUX(1)=0.0
	
C 确定搅拌头中心X位置
      IF ( TIME(2).LE.t1 ) THEN
        xn=x1                  !焊前停留
      ELSE IF ( TIME(2).GT.(t1+((x2-x1)/v)) ) THEN
        xn=x2                  !焊后停留
      ELSE
        xn=x1+((TIME(2)-t1)*v)
      END IF

C 确定搅拌头中心到当前积分点的距离
      
	  rn=SQRT((x-xn)**2+(y-y1)**2)

C 判断积分点所处的区域
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
C 移动热源子程序结束