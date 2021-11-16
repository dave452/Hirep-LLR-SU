/***************************************************************************\
* Copyright (c) 2008, Agostino Patella                                      *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "suN.h"
#include "utils.h"
#include <math.h>

#if (NG!=4) 
#error : Mismatch between NG and ExpX!
#endif

/*
*  U = (1+i*y.sigma/4)(1-i*y.sigma/4)^{-1}
*  U = u[0] + i*u.sigma/4
*/
static void YtoU(double* u, double *y)
{
   double y2 = y[0]*y[0] + y[1]*y[1] +y[2]*y[2];
   double detY = 1.0 + y2/16.;
   u[0] = (1.0 - y2/16.)/detY;
   u[1] = y[0]/(2.*detY);
   u[2] = y[1]/(2.*detY);
   u[3] = y[2]/(2.*detY);
}


/*
*  Applica la rotazione di SU(2) definita da
*  U = s[0] + i*s.sigma/4
*  al vettore (v1,v2)
*/
static void su2_rotate(double *s,complex *v1,complex *v2)
{
   complex z1, z2;
   z1.re=
      s[0]*(*v1).re-s[1]*(*v2).im+s[2]*(*v2).re-s[3]*(*v1).im;
   z1.im=
      s[0]*(*v1).im+s[1]*(*v2).re+s[2]*(*v2).im+s[3]*(*v1).re;
   z2.re=
      s[0]*(*v2).re-s[1]*(*v1).im-s[2]*(*v1).re+s[3]*(*v2).im;
   z2.im=
      s[0]*(*v2).im+s[1]*(*v1).re-s[2]*(*v1).im-s[3]*(*v2).re;
   (*v1)=z1;
   (*v2)=z2;
}


/*
*  Approssima
*  V = exp(sum(A) iT(A)*h(A)*dt)
*  con una matrice unitaria secondo l'algoritmo di Luscher
*  e sostituisce
*  V.U -> U
*/
#ifdef GAUGE_SON
void ExpX(double dt, suNg_algebra_vector *h, suNg *r)
#else
void ExpX(double dt, suNg_algebra_vector *h, suNg *u)
#endif
{
#ifdef WITH_QUATERNIONS
suNg v_tmp, u_tmp;

u_tmp=*u;
_suNg_exp(dt,*h,v_tmp);
_suNg_times_suNg(*u,v_tmp,u_tmp);
#else //WITH_QUATERNIONS 

	double y[3];
	double s[6][4];
	y[0] = +5.000000000000000e-01*(*h).c[3]*(dt);
	y[1] = -5.000000000000000e-01*(*h).c[0]*(dt);
	y[2] = +2.500000000000000e-01*(*h).c[4]*(dt);
	YtoU(s[0],y);
	su2_rotate(s[0],&((*u).c[0]),&((*u).c[4]));
	su2_rotate(s[0],&((*u).c[1]),&((*u).c[5]));
	su2_rotate(s[0],&((*u).c[2]),&((*u).c[6]));
	su2_rotate(s[0],&((*u).c[3]),&((*u).c[7]));
	y[0] = +5.000000000000000e-01*(*h).c[7]*(dt);
	y[1] = -5.000000000000000e-01*(*h).c[1]*(dt);
	y[2] = +1.250000000000000e-01*(*h).c[4]*(dt)+2.165063509461096e-01*(*h).c[9]*(dt);
	YtoU(s[1],y);
	su2_rotate(s[1],&((*u).c[0]),&((*u).c[8]));
	su2_rotate(s[1],&((*u).c[1]),&((*u).c[9]));
	su2_rotate(s[1],&((*u).c[2]),&((*u).c[10]));
	su2_rotate(s[1],&((*u).c[3]),&((*u).c[11]));
	y[0] = +5.000000000000000e-01*(*h).c[8]*(dt);
	y[1] = -5.000000000000000e-01*(*h).c[5]*(dt);
	y[2] = -1.250000000000000e-01*(*h).c[4]*(dt)+2.165063509461096e-01*(*h).c[9]*(dt);
	YtoU(s[2],y);
	su2_rotate(s[2],&((*u).c[4]),&((*u).c[8]));
	su2_rotate(s[2],&((*u).c[5]),&((*u).c[9]));
	su2_rotate(s[2],&((*u).c[6]),&((*u).c[10]));
	su2_rotate(s[2],&((*u).c[7]),&((*u).c[11]));
	y[0] = +5.000000000000000e-01*(*h).c[11]*(dt);
	y[1] = -5.000000000000000e-01*(*h).c[2]*(dt);
	y[2] = +2.041241452319315e-01*(*h).c[14]*(dt)+1.250000000000000e-01*(*h).c[4]*(dt)+7.216878364870322e-02*(*h).c[9]*(dt);
	YtoU(s[3],y);
	su2_rotate(s[3],&((*u).c[0]),&((*u).c[12]));
	su2_rotate(s[3],&((*u).c[1]),&((*u).c[13]));
	su2_rotate(s[3],&((*u).c[2]),&((*u).c[14]));
	su2_rotate(s[3],&((*u).c[3]),&((*u).c[15]));
	y[0] = +5.000000000000000e-01*(*h).c[12]*(dt);
	y[1] = -5.000000000000000e-01*(*h).c[6]*(dt);
	y[2] = +2.041241452319315e-01*(*h).c[14]*(dt)-1.250000000000000e-01*(*h).c[4]*(dt)+7.216878364870322e-02*(*h).c[9]*(dt);
	YtoU(s[4],y);
	su2_rotate(s[4],&((*u).c[4]),&((*u).c[12]));
	su2_rotate(s[4],&((*u).c[5]),&((*u).c[13]));
	su2_rotate(s[4],&((*u).c[6]),&((*u).c[14]));
	su2_rotate(s[4],&((*u).c[7]),&((*u).c[15]));
	y[0] = +5.000000000000000e-01*(*h).c[13]*(dt);
	y[1] = -5.000000000000000e-01*(*h).c[10]*(dt);
	y[2] = +2.041241452319315e-01*(*h).c[14]*(dt)-1.443375672974064e-01*(*h).c[9]*(dt);
	YtoU(s[5],y);
	su2_rotate(s[5],&((*u).c[8]),&((*u).c[12]));
	su2_rotate(s[5],&((*u).c[9]),&((*u).c[13]));
	su2_rotate(s[5],&((*u).c[10]),&((*u).c[14]));
	su2_rotate(s[5],&((*u).c[11]),&((*u).c[15]));
	su2_rotate(s[5],&((*u).c[8]),&((*u).c[12]));
	su2_rotate(s[5],&((*u).c[9]),&((*u).c[13]));
	su2_rotate(s[5],&((*u).c[10]),&((*u).c[14]));
	su2_rotate(s[5],&((*u).c[11]),&((*u).c[15]));
	su2_rotate(s[4],&((*u).c[4]),&((*u).c[12]));
	su2_rotate(s[4],&((*u).c[5]),&((*u).c[13]));
	su2_rotate(s[4],&((*u).c[6]),&((*u).c[14]));
	su2_rotate(s[4],&((*u).c[7]),&((*u).c[15]));
	su2_rotate(s[3],&((*u).c[0]),&((*u).c[12]));
	su2_rotate(s[3],&((*u).c[1]),&((*u).c[13]));
	su2_rotate(s[3],&((*u).c[2]),&((*u).c[14]));
	su2_rotate(s[3],&((*u).c[3]),&((*u).c[15]));
	su2_rotate(s[2],&((*u).c[4]),&((*u).c[8]));
	su2_rotate(s[2],&((*u).c[5]),&((*u).c[9]));
	su2_rotate(s[2],&((*u).c[6]),&((*u).c[10]));
	su2_rotate(s[2],&((*u).c[7]),&((*u).c[11]));
	su2_rotate(s[1],&((*u).c[0]),&((*u).c[8]));
	su2_rotate(s[1],&((*u).c[1]),&((*u).c[9]));
	su2_rotate(s[1],&((*u).c[2]),&((*u).c[10]));
	su2_rotate(s[1],&((*u).c[3]),&((*u).c[11]));
	su2_rotate(s[0],&((*u).c[0]),&((*u).c[4]));
	su2_rotate(s[0],&((*u).c[1]),&((*u).c[5]));
	su2_rotate(s[0],&((*u).c[2]),&((*u).c[6]));
	su2_rotate(s[0],&((*u).c[3]),&((*u).c[7]));


#endif //WITH_QUATERNIONS
}



#define XG(m,a,b) ((m)+(a)*NG+(b))
#define XA(h,a) ((h)+(a))
#define TI(a,b) ((a)*NG+(b)-1)


void ExpX2(double dt, suNg_algebra_vector *h, suNg *u)
{
#ifdef WITH_QUATERNIONS

suNg v_tmp, u_tmp;

u_tmp=*u;
_suNg_exp(dt,*h,v_tmp);
_suNg_times_suNg(*u,v_tmp,u_tmp);

#else //WITH_QUATERNIONS 

	int i, j, k, n;
	double y[3];
	double d[NG];
	double s[NG*(NG-1)/2][4];
	double tmp;
	double *hf = (double*)h;
	complex *uf = (complex*)u;
	
	d[0] = 0.0;
	for(i = 1; i < NG; i++) {
		tmp = sqrt( 2./(i*(i+1)) ) * (*XA(hf,i));
		d[i] = -i *tmp;
		for(j = 0; j < i; j++)
			d[j] += tmp;
	}
	
	k = 0;
	for(j = 1; j<NG; j++) /* Aggiunto il j< testare. Claudio */
	for(i = 0; i < j; i++) {
		y[0] = dt * (*XA(hf,TI(j,i)));
		y[1] = -dt * (*XA(hf,TI(i,j)));
		y[2] = dt * (d[i]-d[j]) / NG;
		YtoU(s[k],y);
		for(n = 0; n < NG; n++)
			su2_rotate(s[k],XG(uf,i,n),XG(uf,j,n));
		k++;
	}

	k = NG*(NG-1)/2 - 1;
	for(j = NG-1; j >= 1; j--)
	for(i = j-1; i >= 0; i--) {
		for(n = 0; n < NG; n++)
			su2_rotate(s[k],XG(uf,i,n),XG(uf,j,n));
		k--;
	}

#endif //WITH_QUATERNIONS

}


