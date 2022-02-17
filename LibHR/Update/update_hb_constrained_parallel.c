/*************************************************************************** \
 * Copyright (c) 2008, Claudio Pica                                          *
 * All rights reserved.                                                      *
\***************************************************************************/

/*******************************************************************************
 *
 * File update.c
 *
 * Update programs
 *
 *******************************************************************************/

#define PROJECT_INTERVAL 10

#include "suN.h"
#include "utils.h"
#include "global.h"
#include "update.h"
#include "communications.h"
#include "logger.h"
#include "random.h"
#include "observables.h"
#define PI 3.141592653589793238462643383279502884197
#include <math.h>
#include <stdlib.h>

static int *dyn_gauge=NULL;

void project_gauge_field(void)
{
  _MASTER_FOR(&glattice,ix) {
    project_to_suNg(pu_gauge(ix,0));
    project_to_suNg(pu_gauge(ix,1));
    project_to_suNg(pu_gauge(ix,2));
    project_to_suNg(pu_gauge(ix,3));
  }

  start_gf_sendrecv(u_gauge);
}

#if defined(BASIC_SF) || defined(ROTATED_SF)
static void g_up_Dirichlet_BCs() {
  int ix,iy,iz,index;

  if(COORD[0] == NP_T-1) {
    for (ix=0;ix<X;++ix) for (iy=0;iy<Y;++iy) for (iz=0;iz<Z;++iz){
	  index=ipt(T-1,ix,iy,iz);
	  dyn_gauge[index*4]=dyn_gauge[index*4+1]=dyn_gauge[index*4+2]=dyn_gauge[index*4+3]=0;
	}
  }
}
#endif

#if defined(BASIC_SF) || defined(ROTATED_SF) || defined(BC_T_MIXED)
static void g_dn_Dirichlet_BCs() {
  int ix,iy,iz,index;

  if(COORD[0] == 0) {
    for (ix=0;ix<X;++ix) for (iy=0;iy<Y;++iy) for (iz=0;iz<Z;++iz){
	  index=ipt(0,ix,iy,iz);
	  dyn_gauge[index*4]=dyn_gauge[index*4+1]=dyn_gauge[index*4+2]=dyn_gauge[index*4+3]=0;
	}
    for (ix=0;ix<X;++ix) for (iy=0;iy<Y;++iy) for (iz=0;iz<Z;++iz){
	  index=ipt(1,ix,iy,iz);
	  dyn_gauge[index*4+1]=dyn_gauge[index*4+2]=dyn_gauge[index*4+3]=0;
	}
  }
}
#endif

#if defined(BC_T_OPEN) || defined(BC_T_MIXED)
static void g_up_open_BCs() {
  int ix,iy,iz,index;

  if(COORD[0] == NP_T-1) {
    for (ix=0;ix<X;++ix) for (iy=0;iy<Y;++iy) for (iz=0;iz<Z;++iz){
	  index=ipt(T-1,ix,iy,iz);
	  dyn_gauge[index*4]=0;
	}
  }
}
#endif

#if defined(BC_T_OPEN)
static void g_dn_open_BCs() {
  int ix,iy,iz,index;

  if(COORD[0] == 0) {
    for (ix=0;ix<X;++ix) for (iy=0;iy<Y;++iy) for (iz=0;iz<Z;++iz){
	  index=ipt(0,ix,iy,iz);
	  dyn_gauge[index*4]=dyn_gauge[index*4+1]=dyn_gauge[index*4+2]=dyn_gauge[index*4+3]=0;
	}
  }
}
#endif

static void free_hb_boundary() {
  if (dyn_gauge!=NULL) {
    free(dyn_gauge);
    dyn_gauge = NULL;
  }
}

static void init_hb_boundary() {
  dyn_gauge = malloc(sizeof(*dyn_gauge)*glattice.gsize_gauge*4);
  atexit(&free_hb_boundary); //register cleanup function at exit

  for(int i=0;i<glattice.gsize_gauge*4;i++) dyn_gauge[i]=1;
#if defined(BASIC_SF) || defined(ROTATED_SF)
  g_up_Dirichlet_BCs();
  g_dn_Dirichlet_BCs();
#endif
#ifdef BC_T_MIXED
  g_up_open_BCs();
  g_dn_Dirichlet_BCs();
#endif
#ifdef BC_T_OPEN
  g_up_open_BCs();
  g_dn_open_BCs();
#endif
}






static void update_all_constrained(double beta,int type, double * S, double Smin, double Smax)
{
  static int count=PROJECT_INTERVAL;

  if (count>=PROJECT_INTERVAL) {
    project_gauge_field();
    count=0;
  }
  ++count;
 //printf("boh= %d\n", glattice.local_master_pieces);
 //printf("boh= %d\n", glattice.master_start[0]);
 //printf("boh= %d\n", glattice.master_end[1]);
    suNg v;
//lprintf("MAIN",10,"Emin, Emax, E = %f, %f, %f\n", Smin, Smax, *S);
    for(int mu=0;mu<4;mu++){
      for(int i=0;i<glat_even.local_master_pieces;i++) {
        _OMP_PRAGMA ( _omp_for )
        for(int k=0;k<glat_even.local_master_pieces;k++) {
          for(int j=glat_even.master_start[k];j<=glat_even.master_end[k];j++){
            if(dyn_gauge[j*4+mu]!=0){
              staples(j,mu,&v);
	    //printf("E = %f, Emin= %f, Emax= %f, i=%d, j=%d\n", S,Smin,Smax,i,j);
              cabmar_constrained(beta,pu_gauge(j,mu),&v,(i == k)? type: 1, S,Smin , Smax);
            }
    //printf("E = %f, Emin= %f, Emax= %f\n", *S, Smin, Smax);
          }
        }
      }
    }
    _OMP_PRAGMA ( master )
          { complete_gf_sendrecv(u_gauge); }
    _OMP_PRAGMA ( barrier )
    for(int mu=0;mu<4;mu++){
      for(int i=0;i<glat_odd.local_master_pieces;i++) {
        _OMP_PRAGMA ( _omp_for )
        for(int k=0;k<glat_odd.local_master_pieces;k++) {
          for(int j=glat_odd.master_start[k];j<=glat_odd.master_end[k];j++){
            if(dyn_gauge[j*4+mu]!=0){
              staples(j,mu,&v);
	    //printf("E = %f, Emin= %f, Emax= %f, i=%d, j=%d\n", S,Smin,Smax,i,j);
              cabmar_constrained(beta,pu_gauge(j,mu),&v,(i == k)? type: 1, S,Smin , Smax);
            }
    //printf("E = %f, Emin= %f, Emax= %f\n", *S, Smin, Smax);
          }
        }
      }
    }
}


/*
for(int mu=0;mu<4;mu++){
  for(int i=0;i<glattice.local_master_pieces;i++) {
      if(CID == i){
        type = 0;
      }else{
        type = 1;
      }
      for(int j=glattice.master_start[i];j<=glattice.master_end[i];j++){
        if(dyn_gauge[j*4+mu]!=0){
          staples(j,mu,&v);
  //printf("E = %f, Emin= %f, Emax= %f, i=%d, j=%d\n", S,Smin,Smax,i,j);
          cabmar_constrained(beta,pu_gauge(j,mu),&v,type, S,Smin , Smax);
        }
//printf("E = %f, Emin= %f, Emax= %f\n", *S, Smin, Smax);
      }
  }
}
*/

void update_constrained_parallel(double beta,int nhb,int nor, double * S, double Smin, double Smax)
{
//  int mpiret; (void)mpiret;
  if(dyn_gauge==NULL ) init_hb_boundary();

  for (int n=0;n<nhb;n++){
    for (int i = 0; i<CART_SIZE;i++){
      if(CID == i)
      {
        update_all_constrained(beta,0,  S, Smin, Smax);
      }
      else
      {
        update_all_constrained(beta,1,  S, Smin, Smax);
      }
      MPI_Barrier(GLB_COMM);
      start_gf_sendrecv(u_gauge);
      bcast_from_rank(S, 1, i);
    }
  }

  for (int n=0;n<nor;n++){
    update_all_constrained(beta,1,  S, Smin, Smax);
  }
  MPI_Barrier(GLB_COMM);
  start_gf_sendrecv(u_gauge);

}


static void update_all_anneal(double beta, double *S)
{
  double Smin =  -1.*GLB_VOLUME*6.;
  double Smax = GLB_VOLUME*6.;
  static int count=PROJECT_INTERVAL;

  if (count>=PROJECT_INTERVAL) {
    project_gauge_field();
    count=0;
  }
  ++count;

_OMP_PRAGMA ( _omp_parallel )
  {
    suNg v;

    for(int mu=0;mu<4;mu++){
#ifdef WITH_MPI
_OMP_PRAGMA ( master )
      { start_gf_sendrecv(u_gauge); }
_OMP_PRAGMA ( barrier )
#endif
_OMP_PRAGMA ( _omp_for )
      for(int j=glat_even.master_start[0];j<=glat_even.master_end[0];j++){
        if(dyn_gauge[j*4+mu]!=0){
          staples(j,mu,&v);
          cabmar_constrained(beta,pu_gauge(j,mu),&v,0, S,Smin , Smax);
        }
      }
#ifdef WITH_MPI
_OMP_PRAGMA ( master )
      { complete_gf_sendrecv(u_gauge); }
_OMP_PRAGMA ( barrier )
#endif
      for(int i=1;i<glat_even.local_master_pieces;i++) {
_OMP_PRAGMA ( _omp_for )
        for(int j=glat_even.master_start[i];j<=glat_even.master_end[i];j++){
          if(dyn_gauge[j*4+mu]!=0){
            staples(j,mu,&v);
            cabmar_constrained(beta,pu_gauge(j,mu),&v,0, S,Smin , Smax);
          }
        }
      }
    }

    for(int mu=0;mu<4;mu++){
#ifdef WITH_MPI
_OMP_PRAGMA ( master )
      { start_gf_sendrecv(u_gauge); }
_OMP_PRAGMA ( barrier )
#endif
_OMP_PRAGMA ( _omp_for )
      for(int j=glat_odd.master_start[0];j<=glat_odd.master_end[0];j++) {
        if(dyn_gauge[j*4+mu]!=0){
          staples(j,mu,&v);
          cabmar_constrained(beta,pu_gauge(j,mu),&v,0, S,Smin , Smax);
        }
      }
#ifdef WITH_MPI
_OMP_PRAGMA ( master )
      { complete_gf_sendrecv(u_gauge); }
_OMP_PRAGMA ( barrier )
#endif
      for(int i=1;i<glat_odd.local_master_pieces;i++) {
_OMP_PRAGMA ( _omp_for )
        for(int j=glat_odd.master_start[i];j<=glat_odd.master_end[i];j++) {
          if(dyn_gauge[j*4+mu]!=0){
            staples(j,mu,&v);
            cabmar_constrained(beta,pu_gauge(j,mu),&v,0, S,Smin , Smax);
          }
        }
      }
    }

  }

}



int anneal_parallel(double beta, double dbeta, double *S, double S0, double dS){
  if(dyn_gauge==NULL ) init_hb_boundary();
  int k = 1;
  double S_old = *S;
  double Emin = S0 - dS*0.5;
  double Emax = S0 + dS*0.5;
  double Smin =  -1.*GLB_VOLUME*6.;
  double Smax = GLB_VOLUME*6.;
  while(((*S < Emin) ||  (*S > Emax)) && k < 10000){
    if(k%10==0)
    {
      if(((S_old - S0) > 0) && ((*S - S0) < 0) || ((S_old - S0) < 0) && ((*S - S0) > 0))
      {
          dbeta /= 2.;
      }
      if (*S < S0)
      {
        beta += dbeta;
      }else{
        beta -= dbeta;
      }

      S_old = *S;
      lprintf("ANNEAL",0,"beta = %f, S =  %f, %d ...\n", beta, S_old, k);
    }
//    update_all_anneal(beta, S);
    for (int i = 0; i<CART_SIZE;i++){
      //lprintf("ANNEAL",0,"%d ", i);
      if(CID == i)
      {
        update_all_constrained(beta,0,  S, Smin, Smax);
      }
      else
      {
        update_all_constrained(beta,1,  S, Smin, Smax);
      }
      MPI_Barrier(GLB_COMM);
      start_gf_sendrecv(u_gauge);
      lprintf("ANNEAL",0,"\n %d, Bef, CID=%d, S=%f \n", k, i, *S);
      bcast_from_rank(S, 1, i);
      lprintf("ANNEAL",0,"\n %d, Aft, CID=%d, S=%f \n", k, i, *S);
      lprintf("ANNEAL",0,"%d, Avr_plaquette, CID=%d, S=%f \n", k, i,  avr_plaquette()*GLB_VOLUME*6.);
    }
    k++;
  }
  S_old = *S;
  if(k > 10000){
    lprintf("ANNEAL",0,"Couldn't bring system to the interval beta = %f, S =  %f ...\n", beta, S_old);
    return 1;
  }
  lprintf("ANNEAL",0,"Finished S = %f, Smin = %f, Smax = %f, in %d steps ...\n",  S_old, Emin, Emax, k);
  lprintf("ANNEAL",0,"Check action S_an = %f, S_ap = %f",  S_old, avr_plaquette()*GLB_VOLUME*6.);
  MPI_Barrier(MPI_COMM_WORLD);
  lprintf("ANNEAL",0,"After barrier");
  return 0;
}




