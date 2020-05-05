#include "udf.h"
#include <math.h>

DEFINE_PROPERTY(cell_viscosity, c, t)
{
	real visc;
	real cond1;
	real sellars;
	real temp = C_T(c, t);
	real sr = C_STRAIN_RATE_MAG(c, t);
	real time = RP_Get_Real("flow-time");
	real alpha = 1.6e-8;
	real Q = 1.4888e5;
	real n = 4.27;
	real R = 8.314;
	real A = 3.25e8;
	cond1 = sr * exp(Q / (R * (temp + 273.15)));
	sellars = log(pow(cond1 / A, 1.0 / n) + pow((1 + pow(cond1 / A, 2.0 / n)), 0.5));

	if (sr < 1.e-5)
		sr = 1.e-5;
	else
		sr = sr;

	if (time < 0.02)
		visc = 10.;
	else
		visc = sellars / (3.0 * sr * alpha);

	return visc;
}

DEFINE_PROPERTY(cell_conduct, c, t)
{
	real lambda;
	real temp = C_T(c, t);
	lambda = 103.264 + 0.241 * temp;
	return lambda;
}

DEFINE_SPECIFIC_HEAT(my_user_cp, T, Tref, h, yi)
{
	real cp = 754.08 + 0.3729 * T + 0.0012 * T * T;
	*h = cp * (T - Tref);
	return cp;
}

DEFINE_PROFILE(heat_flux, thread, index)
{
	real xc[ND_ND];
	face_t f;
	real t;
	real heat;
	real x, y, z, r, v, temp;
	real tempa, cond1, cond2, cond3, sigma, shear;
	real heat_pressure, heat_shear;
	real alpha = 1.6e-8;
	real Q = 1.4888e5;
	real n = 4.27;
	real R = 8.314;
	real A = 3.25e8;
	real delta;
	real fric_eff = 0.3;
	real beta2 = 0.3;
	real p = 5000. * 1e6 / (3.14 * 8. * 8.);

	t = RP_Get_Real("flow-time");
	begin_f_loop(f, thread)
	{
		F_CENTROID(xc, f, thread);
		x = xc[0];
		y = xc[1];
		z = xc[2];
		r = sqrt(x * x + z * z);
		v = 83.77 * r;
		temp = F_T(f, thread);

		cond1 = n * log(2.0) - log(A) + Q / ((temp + 273.15) * R);
		cond2 = alpha * n;
		cond3 = cond1 / cond2;

		if (temp > 630.0)
			tempa = 630.0;
		else
			tempa = temp;

		sigma = cond3 * (1 - sqrt(tempa / 630.85)) + 2.e7;
		shear = sigma / 1.732;

		delta = 0.31 * exp(83.77 * r / 1.87) - 0.026;
		heat_shear = shear * delta;
		heat_pressure = (1 - delta) * fric_eff * p;
		heat = beta2 * (heat_shear + heat_pressure) * v;

		F_PROFILE(f, thread, index) = heat;
	}
	end_f_loop(f, thread)
}

DEFINE_SOURCE(t_source, c, t, dS, eqn)
{
	real source;
	real miu = C_MU_EFF(c, t);
	real dudx = C_DUDX(c, t);
	real dudy = C_DUDY(c, t);
	real dudz = C_DUDZ(c, t);
	real dvdx = C_DVDX(c, t);
	real dvdy = C_DVDY(c, t);
	real dvdz = C_DVDZ(c, t);
	real dwdx = C_DWDX(c, t);
	real dwdy = C_DWDY(c, t);
	real dwdz = C_DWDZ(c, t);
	real cond1, cond2, cond3, cond4;
	real beta1 = 0.25;

	cond1 = (dudx * dudx + dvdy * dvdy + dwdz * dwdz) * 2.0;
	cond2 = (dudy + dvdx) * (dudy + dvdx);
	cond3 = (dudz + dwdx) * (dudz + dwdx);
	cond4 = (dvdz + dwdy) * (dvdz + dwdy);
	source = beta1 * miu * (cond1 + cond2 + cond3 + cond4);

	return source;
}

DEFINE_PROFILE(velocity, thread, index)
{
	real xc[ND_ND];
	face_t f;
	real t;
	real velo;
	real x, y, z, r, v, delta, omega;

	t = RP_Get_Real("flow-time");
	begin_f_loop(f, thread)
	{

		F_CENTROID(xc, f, thread);
		x = xc[0];
		y = xc[1];
		z = xc[2];
		r = sqrt(x * x + z * z);
		delta = 0.31 * exp(83.77 * r / 1.87) - 0.026;
		omega = delta * 83.77;
		F_PROFILE(f, thread, index) = omega;
	}
	end_f_loop(f, thread)
}