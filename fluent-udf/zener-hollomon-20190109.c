#include "udf.h"
#include<math.h>
 
 DEFINE_PROPERTY(cell_viscosity,c,t)
 {
    real eta;
	real Z;
	real ZZ;
	real tempa;
    real temp = C_T(c,t);
	real sr=C_STRAIN_RATE_MAG(c,t);
	real alpha=1.6e-8;
	real Q=1.4888e5;
	real n=4.27;
	real R=8.314;
	real A=3.25e8;
	Z=sr*exp(Q/(R*(temp+273.15)));
	ZZ=log(pow(Z/A,1.0/n)+pow((1+pow(Z/A,2.0/n)),0.5));

	if(sr<1.e-5)
		sr=1.e-5;
	else
		sr=sr;
	
	eta=ZZ/(3.0*sr*alpha);


	if(eta>5.e12)
		eta=5.e12;
	else
		eta=eta;
    return eta;
 } 
 
  DEFINE_PROPERTY(cell_conduct,c,t)
 {
    real lambda;
	real temp = C_T(c,t);
	lambda=103.264+0.241*temp;
    return lambda;
 } 
 
DEFINE_SPECIFIC_HEAT(my_user_cp, T, Tref, h, yi)
 {
   real cp=754.08+0.3729*T+0.0012*T*T;
   *h = cp*(T-Tref);
   return cp;
 }
 
DEFINE_PROFILE(heat_flux,thread,index)
{
	real xc[ND_ND];
	face_t f;
	real t;
	real heat;
	real x,y,z,r,v,temp;
	real tempa,ina,inb,inc,sigma;
	real eta=1.6e-8;
	real Q=1.4888e5;
	real n=4.27;
	real R=8.314;
	real A=3.25e8;
	real delta;
	
	t=RP_Get_Real("flow-time");
	 begin_f_loop(f,thread)
	 {
		 F_CENTROID(xc,f,thread);
			x=xc[0];
			y=xc[1];
			z=xc[2];
			r=sqrt(x*x+z*z);
			v=83.77*r;
			temp=F_T(f,thread);
			
			ina=n*log(2.0)-log(A)+Q/((temp+273.15)*R);
			inb=eta*n;
			inc=ina/inb;
	
			if(temp>630.0)
				tempa=630.0;
			else
				tempa=temp;
		
			sigma=inc*(1-sqrt(tempa/630.85))+2.e7;
			
			delta=0.31*exp(83.77*r/1.87)-0.026;
			heat=sigma*(1.0-delta)*v*0.77/1.732;
	
		 F_PROFILE(f,thread,index)=heat;
	 }
	 end_f_loop(f,thread)
}

DEFINE_SOURCE(t_source,c,t,dS,eqn)
{        
	real source;
	real strain;
	real stress;
	real Z;
	real ZZ;
	real tempa;
    	real temp = C_T(c,t);
	real sr=C_STRAIN_RATE_MAG(c,t);
	real alpha=1.6e-8;
	real Q=1.4888e5;
	real n=4.27;
	real R=8.314;
	real A=3.25e8;
	Z=sr*exp(Q/(R*(temp+273.15)));
	ZZ=log(pow(Z/A,1.0/n)+pow((1+pow(Z/A,2.0/n)),0.5));
	stress=ZZ/alpha;
	strain=sr*CURRENT_TIMESTEP;
	source=stress*strain*0.5;
        return source;
       
}

DEFINE_PROFILE(velocity,thread,index)
{
	real xc[ND_ND];
	face_t f;
	real t;
	real velo;
	real x,y,z,r,v,delta,omega;
	
	t=RP_Get_Real("flow-time");
	 begin_f_loop(f,thread)
	 {

		F_CENTROID(xc,f,thread);
			x=xc[0];
			y=xc[1];
			z=xc[2];
		r=sqrt(x*x+z*z);
		delta=0.31*exp(83.77*r/1.87)-0.026;
		omega=delta*83.77;
		 F_PROFILE(f,thread,index)=omega;
	 }
	 end_f_loop(f,thread)
}