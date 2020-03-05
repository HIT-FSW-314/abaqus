      SUBROUTINE DFLUX(FLUX,SOL,JSTEP,JINC,TIME,NOEL,NPT,COORDS,JLTYP,
     1                 TEMP,PRESS,SNAME)
C
      INCLUDE 'ABA_PARAM.INC'


      DIMENSION COORDS(3),FLUX(2),TIME(2)
      CHARACTER*80 SNAME

      wu=110.0
      wi=10.0
      effi=0.8
      v=0.005
      q=wu*wi*effi
      d=v*TIME(2)

      x=COORDS(1)
      y=COORDS(2)
      z=COORDS(3)

      x0=0.0
      y0=0.0
      z0=0.003
	  
      a=0.003
	b=0.003
	c=0.003
	aa=0.006
	f1=0.7
	PI=3.1415926
	beta=0.0
	aacos=a/cos(beta/180.0*PI)
	cccos=c*cos(beta/180.0*PI)
	aaacos=aa/cos(beta/180.0*PI)
	cccos=c*cos(beta/180.0*PI)
	 
	heat1=6.0*sqrt(3.0)*q/(a*b*c*PI*sqrt(PI))*f1
	heat2=6.0*sqrt(3.0)*q/(aa*b*c*PI*sqrt(PI))*(2.0-f1)
	  
	shape1=exp(-3.0*(x-x0-d)**2/aacos**2-3.0*(y-y0)**2/b**2
     $	-3.0*(z-z0)**2/cccos**2)
      shape2=exp(-3.0*(x-x0-d)**2/aaacos**2-3.0*(y-y0)**2/b**2
     $	-3.0*(z-z0)**2/cccos**2)

C     JLTYP＝1，表示为体热源
      JLTYP=1
	  IF(x .GE.(x0+d)) THEN
	  FLUX(1)=heat1*shape1
	  ELSE
	  FLUX(1)=heat2*shape2
	  ENDIF
      RETURN
      END