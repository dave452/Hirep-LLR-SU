/*************************************************************************** \
 * Copyright (c) 2008, Claudio Pica                                          *   
 * All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
 *
 * File update_llr.c
 *
 * Update programs
 *
 *******************************************************************************/
#include "suN.h"
#include "utils.h"
#include "global.h"
#include "logger.h"
#include "random.h"
#include "communications.h"
#include <math.h>
#include "update.h"
#include "observables.h"
#include <stdlib.h>

int initHB = 0;

typedef struct {
  int nrm,nth;
  int it;
  double starta;
  double a;
  double S0;
  double dS;
#ifdef LLRHB
  double E;
#endif
} llrparams;

static llrparams llrp;

//reset it to 0 and rhoa to initial value
void restart_robbinsmonro(){
  llrp.it=8;
  llrp.a=llrp.starta;
}

void init_robbinsmonro(int nrm,int nth,double starta,int it,double dS,double S0){
  llrp.nrm=nrm;
  llrp.nth=nth;
  llrp.it=it;
  llrp.starta=starta;
  llrp.dS=dS;
  llrp.S0=S0;
#ifdef LLRHB
  if(!initHB){  
  llrp.E = avr_plaquette()*GLB_VOLUME*6.;
  lprintf("MAIN",0,"Bringing the system to the interval (S0,dS) = (%f, %f) ...\n", llrp.S0, llrp.dS);
  anneal(&(llrp.E), llrp.S0, llrp.dS);
  lprintf("MAIN",0,"System brought to the interval (S0,dS) = (%f, %f)\n", llrp.S0, llrp.dS);
	}
#endif
}
  
double get_llr_a(void){
  return llrp.a;	
}

double getS0(void){
  return llrp.S0;	
}

#ifdef LLRHB
double getE(void){
  return llrp.E;	
}
#endif
double getdS(void){
  return llrp.dS;	
}



void thermrobbinsmonro(void){
#ifdef LLRHB
  double Emin, Emax; 
  Emin = llrp.S0 - .5*llrp.dS;
  Emax = llrp.S0 + .5*llrp.dS;
  update_constrained(llrp.a, 1,0, &(llrp.E),Emin,Emax);
#else
  double S_llr,S_non_llr;
  update_llr_ghmc(&S_llr,&S_non_llr,1);
#endif
}

void llr_fixed_a_update(void){
#ifdef LLRHB
  double Emin, Emax; 
  Emin = llrp.S0 - .5*llrp.dS;
  Emax = llrp.S0 + .5*llrp.dS;
  for( int i=0; i<200 ; i++) {
	update_constrained(llrp.a, 1,0, &(llrp.E),Emin,Emax);
    	lprintf("ROBBINSMONRO",10,"Fixed a MC Step: %d E=%lf \n",i,llrp.E);
	polyakov();
	}
#ifdef WITH_UMBRELLA
  lprintf("ROBBINSMONRO",10,"Fixed a Swap : %d E=%lf \n",llrp.E);
  umbrella_swap(&(llrp.E),&llrp.S0,&llrp.a,&llrp.dS);
#endif
#else
  double S_llr,S_non_llr;
  update_llr_ghmc(&S_llr,&S_non_llr,0);
#ifdef WITH_UMBRELLA
  umbrella_swap(&S_llr,&llrp.S0,&llrp.a,&llrp.dS);
#endif
#endif

}


	
void robbinsmonro(void){
  
  int rmstep;
#ifdef LLRHB
  double Emin, Emax;
  Emin = llrp.S0 - .5*llrp.dS;
  Emax = llrp.S0 + .5*llrp.dS;
#else
  double S_llr;
  double S_non_llr;
#endif
    
  
  for(rmstep=0;rmstep<llrp.nth;rmstep++){
    //lprintf("llr",30,"Therm: %d\n",rmstep);
#ifdef LLRHB
  lprintf("llr",30,"Starting therm...\n");
  update_constrained(llrp.a, 1,0, &(llrp.E),Emin, Emax);
#else
  update_llr_ghmc(&S_llr,&S_non_llr,1);
#endif

  }
  
  
  
  
  double avr=0.;
  for(rmstep=0;rmstep<llrp.nrm;rmstep++){
#ifdef LLRHB
    update_constrained(llrp.a, 1,0,&(llrp.E),Emin,Emax);
    avr += llrp.E;
    //lprintf("ROBBINSMONRO",10,"RM Step: %d GMC Iter: %d E=%lf \n",llrp.it,rmstep,llrp.E);
    //if( rmstep%100 ==0 ) umbrella_swap(&(llrp.E),&llrp.S0,&llrp.a,&llrp.dS);
#else
    update_llr_ghmc(&S_llr,&S_non_llr,0);
    avr+=S_llr;
    //if( rmstep%100 ==0 ) umbrella_swap(&S_llr,&llrp.S0,&llrp.a,&llrp.dS);
    //lprintf("ROBBINSMONRO",10,"RM Step: %d GMC Iter: %d S_llr=%lf \n",llrp.it,rmstep,S_llr);
#endif
    
  }
  
  avr/=(double)llrp.nrm;
  lprintf("ROBBINSMONRO",10,"RM avr-S0: %lf delta_a: %lf \n",avr-llrp.S0,(avr-llrp.S0)*12./(llrp.dS*llrp.dS*llrp.it) );
  //llrp.a-=(avr-llrp.S0)*12./(llrp.dS*llrp.dS);
  llrp.a-=(avr-llrp.S0)*12./(llrp.dS*llrp.dS*llrp.it);
#ifdef WITH_UMBRELLA
#ifdef LLRHB
  //lprintf("ROBBINSMONRO",10,"Emin=%lf E=%lf Emax = %lf\n",Emin,llrp.E,Emax);
  if( llrp.it%2 == 0 ) umbrella_swap(&(llrp.E),&llrp.S0,&llrp.a,&llrp.dS);
 // if( rmstep%100 ==0 ) umbrella_swap(&(llrp.E),&llrp.S0,&llrp.a,&llrp.dS);
#else
  umbrella_swap(&S_llr,&llrp.S0,&llrp.a,&llrp.dS);
#endif
#endif
  llrp.it++;
}


typedef struct{
  double S_llr;
  double dS;
  double S0;
  double a;
  double deltaS;
  int rep;
  int repnext;

} reppar;


#ifdef WITH_UMBRELLA
static int compare_S0(const void *p, const void *q) {
  //lprintf("llr:compare",0,"inside compare \n");
  reppar x=*(reppar *)p;
  reppar y=*(reppar *)q;


  if (x.S0 < y.S0 )
    return -1;  // Return -1 if you want ascending, 1 if you want descending order. 
  else if (x.S0 > y.S0 )
    return 1;   // Return 1 if you want ascending, -1 if you want descending order. 

  return 0;
}


static int compare_deltaS(const void *p, const void *q) {
  //lprintf("llr:compare",0,"inside compare \n");
  reppar x=*(reppar *)p;
  reppar y=*(reppar *)q;


  if (x.deltaS < y.deltaS )
    return -1;  // Return -1 if you want ascending, 1 if you want descending order. 
  else if (x.deltaS > y.deltaS )
    return 1;   // Return 1 if you want ascending, -1 if you want descending order. 

  return 0;
}


void swap(double *data){
  reppar drep[N_REP];
  int toswap[N_REP];
  double rand,temp;
  int i,j;
  //for(i=0; i< 4*N_REP;  i++) printf("prova data[%d] = %f\n", i, data[i]);
#ifdef LLRHB
 toswap[N_REP-1]=0;
#else 
 toswap[N_REP-1]=1;
#endif

  //sorting the replicas with respect of their energy  S0 so we can swap nearest neighbourgh

  for(i=0;i<N_REP;i++){
    drep[i].S_llr=data[4*i];
    drep[i].S0=data[4*i+1];
    drep[i].a=data[4*i+2];
    drep[i].dS=data[4*i+3];
    drep[i].rep=i;
  }

  qsort(drep,N_REP,sizeof(drep[0]),compare_S0);
 
  //sorting with respect of the difference in hamiltonian between r and r+1 where r is the replica

  for(i=0;i<N_REP-1;i++){
    drep[i].repnext=drep[i+1].rep;
#ifdef LLRHB	
    //if( drep[i].S_llr > drep[i].S0 && drep[i+1].S_llr < drep[i+1].S0 ) {
    drep[i].deltaS = (drep[i+1].a - drep[i].a) * ( drep[i+1].S_llr - drep[i].S_llr) ;
#else
    toswap[i]=1;
    double S1=-drep[i+1].a*drep[i].S_llr-drep[i].a*drep[i+1].S_llr;
    double tmp1=(drep[i].S_llr-drep[i+1].S0)/drep[i+1].dS;
    double tmp2=(drep[i+1].S_llr-drep[i].S0)/drep[i].dS;
    double S2=(tmp1*tmp1+tmp2*tmp2)/2.;
    double S3=-drep[i+1].a*drep[i+1].S_llr-drep[i].a*drep[i].S_llr;
    tmp1=(drep[i].S_llr-drep[i].S0)/drep[i].dS;
    tmp2=(drep[i+1].S_llr-drep[i+1].S0)/drep[i+1].dS;
    double S4=(tmp1*tmp1+tmp2*tmp2)/2.;
    drep[i].deltaS = S1+S2-S3-S4;
#endif

  }
  
#ifndef LLRHB 
  qsort(drep,N_REP-1,sizeof(drep[0]),compare_deltaS);
#endif

#ifdef LLRHB
  for(i=0;i<N_REP-1;i++){
    if( drep[i].S_llr > drep[i].S0 && drep[i+1].S_llr < drep[i+1].S0 ) toswap[drep[i].rep] = 1;
    else toswap[ drep[i].rep ] = 0;
  //  lprintf("SWAP",10,"Replica %d: S0= %lf, E=%lf, next: %d, swap: %d\n", drep[i].rep, drep[i].S0, drep[i].S_llr, drep[i].repnext, toswap[drep[i].rep]);
  }
 // lprintf("SWAP",10,"Replica %d: S0= %lf, E=%lf\n", drep[N_REP-1].rep, drep[N_REP-1].S0, drep[N_REP-1].S_llr);
#endif

 
  //swap of replicas
  for(i=0;i<N_REP-1;i++){
#ifdef LLRHB
    if( toswap[drep[i].rep] ){
      toswap[drep[i].rep]=0;
#else
    if( toswap[drep[i].rep] && toswap[drep[i].repnext] ){
      toswap[drep[i].rep]=0;
      toswap[drep[i].repnext]=0;
#endif
      if(drep[i].deltaS<0){
	for(j=0;j<4;j++){ 
	  temp=data[4*drep[i].rep+j];
          data[4*drep[i].rep+j]=data[4*drep[i].repnext+j];
          data[4*drep[i].repnext+j]=temp;
	//lprintf("SWAP",10,"Replica %d with Emin E Emax a : %f, %f, %f, %f\n",drep[i].rep,drep[i].S0-0.5*drep[i].dS,drep[i].S_llr,drep[i].S0-0.5*drep[i].dS,drep[i].a);
	}
//	lprintf("SWAP", 10, "Replicas %d and %d swapped!\n",drep[i].rep,drep[i].repnext);
	
      }else{
	ranlxd(&rand,1);
	if(rand < exp(-drep[i].deltaS)){
	  for(j=0;j<4;j++){
            temp=data[4*drep[i].rep+j];
	    data[4*drep[i].rep+j]=data[4*drep[i].repnext+j];
	    data[4*drep[i].repnext+j]=temp; 
	  }	
	lprintf("SWAP", 10, "Replicas %d and %d swapped!\n",drep[i].rep,drep[i].repnext); 
	}
      }
    }
  }
}


void setreplica(double *data){
  lprintf("llr:setreplica",0,"Updating OLD LLR Param: S0 %lf,  a  %lf , dS %lf  \n",llrp.S0,llrp.a,llrp.dS);
  llrp.S0=data[1];
  llrp.dS=data[3];
  llrp.a=data[2];
  lprintf("llr:setreplica",0,"New LLR Param: S0 %lf,  a  %lf , dS %lf  \n",llrp.S0,llrp.a,llrp.dS);
  
}

#endif //WITH_UMBRELLA
//void swap_hb(double *data){
//  reppar drep[N_REP];
//  int toswap[N_REP];
//  double rand,temp;
//  int i,j;
//
//  toswap[N_REP-1]=1;
//
//  //sorting the replicas with respect of their energy  S0 so we can swap nearest neighbourgh
//
//  for(i=0;i<N_REP;i++){
//    
//    drep[i].S_llr=data[4*i];
//    drep[i].S0=data[4*i+1];
//    drep[i].a=data[4*i+2];
//    drep[i].dS=data[4*i+3];
//    drep[i].rep=i;
//
//  }
//
//  qsort(drep,N_REP,sizeof(drep[0]),compare_S0);
//
//  //sorting with respect of the difference in hamiltonian between r and r+1 where r is the replica
//
//  for(i=0;i<N_REP-1;i++){
//    toswap[i]=1;
//    drep[i].repnext=drep[i+1].rep;
//    drep[i].deltaS= (drep[i+1].a - drep[i].a) * ( drep[i+1].S_llr - drep[i].S_llr) ;
//  }
//  
//  qsort(drep,N_REP-1,sizeof(drep[0]),compare_deltaS);
//  
//  //swap of replicas
//
//  for(i=0;i<N_REP-1;i++){
//    if(toswap[drep[i].rep] && toswap[drep[i].repnext] ){
//      toswap[drep[i].rep]=0;
//      toswap[drep[i].repnext]=0;
//      if(drep[i].deltaS<0){
//	for(j=0;j<4;j++){                                                                                                                                                                                                                      
//	  temp=data[4*drep[i].rep+j];                                                                                                                                                                                                                      
//	  data[4*drep[i].rep+j]=data[4*drep[i].repnext+j];                                                                                                                                                                                                           
//	  data[4*drep[i].repnext+j]=temp;
//	}
//	lprintf("SWAP", 10, "Replicas %d and %d swapped!\n",drep[i].rep,drep[i].repnext); 
//	
//      }else{
//	ranlxd(&rand,1);
//	if(rand < exp(-drep[i].deltaS)){
//	  for(j=0;j<4;j++){                                                                                                                                                                                                                      
//	    temp=data[4*drep[i].rep+j];                                                                                                                                                                                                                      
//	    data[4*drep[i].rep+j]=data[4*drep[i].repnext+j];                                                                                                                                                                                                           
//	    data[4*drep[i].repnext+j]=temp; 
//	  }	
//	lprintf("SWAP", 10, "Replicas %d and %d swapped!\n",drep[i].rep,drep[i].repnext); 
//	}
//      }
//    }
//  }
//}
