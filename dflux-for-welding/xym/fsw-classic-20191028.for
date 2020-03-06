C example for fsw_classic heating source
C fixed format for fortran77 (programming area 7-72)

C DFLUX for friction stir welding
C Systeme International
C welding speed 300 mm/min
C rotational velocity 600 rpm
C cylinder pin without extra topology

C ---------------------------------------
      SUBROUTINE DFLUX(FLUX,SOL,KSTEP,KINC,TIME,NOEL,NPT,COORDS,
     & JLTYP,TEMP,PRESS,SNAME)	
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION FLUX(2), TIME(2), COORDS(3)
      CHARACTER*80 SNAME
      real x,y,z ! coordinates
      real x0,y0,z0 ! initial position
      real x1,y1,z1 ! final position
      real t0,t1 ! dwelling time
      real v0_1,omega0_1 ! welding parameters 0_1
      real r_shoulder,r_pinr,l_pin ! shape parameters
      real P,miu,eff_heat,k ! internal parameters
      real x_c,y_c,z_c ! current position of welding tool
      real r_c ! radius between intergal point and tool
      real force ! frictional force
      real weld_t0_1 ! final time for weld 0_1
      REAL PI; PARAMETER ( PI= 3.14159 )

C ----------------------------------------	  
C definition of coordinates
      x=COORDS(1) ! coordinate of X
      y=COORDS(2) ! coordinate of Y
      z=COORDS(3) ! coordinate of Z
C initial position of welding tool
      x0=0.
      y0=0.
      z0=0.
C final position of welding tool
      x1=0.025
      y1=0.
      z1=0.
C dwelling time
      dwell_t0=3. ! initial
      dwell_t1=3. ! final
C welding speed and rotational velocity
      v0_1=0.005 ! welding speed
      omega0_1=600.*2.*PI/60. ! rotational velocity
C shape of welding tool
      r_shoulder=0.007 ! radius of shoulder
      r_pinr=0.002 ! radius of pin root
      l_pin=0.003 ! length of pin
C internal parameters
      P=15286242.0 ! axial pressure of welding tool
      miu=0.37 ! coefficient of friction
      eff_heat=0.98 ! coefficient of heat generation
      k=0.33 ! heating ratio of pin to shoulder
C frictional force
      if (SOL.lt.390) then
        force=miu*P
      else if (SOL.lt.450) then
        force=0.577*(98.1-(98.1-47)*(SOL-390)/(450-390))*1e6
      else if (SOL.lt.550) then
        force=0.577*(47-(47-5)*(SOL-450)/(550-450))*1e6
      else
        force=0.577*1e6
      end if

C ----------------------------------------
C initialize the flux value
      FLUX(1)=0.
C real-time position of welding tool
      weld_t0_1=dwell_t0+(sqrt((x1-x0)**2+(y1-y0)**2+
     & (z1-z0)**2)/v0_1)
      if (TIME(2).le.dwell_t0) then
        x_c=x0
        y_c=y0
        z_c=z0
      else if (TIME(2).gt.weld_t0_1) then
        x_c=x1
        y_c=y0
        z_c=z0
      else
        x_c=x0+(v0_1*(TIME(2)-dwell_t0))
        y_c=y0
        z_c=z0
      end if
C distance between tool and integral point
      r_c=sqrt((x-x_c)**2+(y-y_c)**2+(z-z_c)**2)
C fsw_classic heating source
      if ((JLTYP.eq.0).and.(r_c.le.r_shoulder)) then
      	FLUX(1)=force*r_c*omega0_1*eff_heat
      else if ((JLTYP.EQ.1)
     &            .and.(r_c.le.r_pinr)
     &            .and.(z.ge.(z0-l_pin))) then
      	FLUX(1)=(2.*k*force*omega0_1*eff_heat*
     &            (r_shoulder**3-r_pinr**3))/
     &            (3.*l_pin*r_pinr**2)
      else
      	FLUX(1)=0.
      end if

      return
      END

