/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica                                          *
 * All rights reserved.                                                      *
 \***************************************************************************/
/**************************************
Modified 2021 by David Mason and Davide Vadacchino
***************************************/
#include "global.h"
#include "communications.h"
#include "logger.h"
#include "error.h"
#include "geometry.h"
#include "spinor_field.h"
#include "suN_types.h"
#include "global.h"
#include "update.h"
#include "utils.h"
#include <string.h>
#include <stdio.h>
#ifdef WITH_MPI
#include <mpi.h>
#endif
#include "logger.h"

#ifdef MPI_TIMING
struct timeval gfstart, gfend, gfetime,sfstart, sfend, sfetime;
int gf_control=0,sf_control=0;
#endif

void global_sum(double *d, int n) {
#ifdef WITH_MPI
  int mpiret;(void)mpiret; // Remove warning of variable set but not used
  double pres[n];

#ifdef MPI_TIMING
  struct timeval start, end, etime;
  gettimeofday(&start,0);
#endif

  mpiret=MPI_Allreduce(d,pres,n,MPI_DOUBLE,MPI_SUM,GLB_COMM);


#ifdef MPI_TIMING
  gettimeofday(&end,0);
  timeval_subtract(&etime,&end,&start);
  lprintf("MPI TIMING",0,"global_sum " __FILE__ " %ld sec %ld usec\n",etime.tv_sec,etime.tv_usec);
#endif

#ifndef NDEBUG
  if (mpiret != MPI_SUCCESS) {
    char mesg[MPI_MAX_ERROR_STRING];
    int mesglen;
    MPI_Error_string(mpiret,mesg,&mesglen);
    lprintf("MPI",0,"ERROR: %s\n",mesg);
    error(1,1,"global_sum " __FILE__,"Cannot perform global_sum");
  }
#endif
  while(n>0) {
    --n;
    d[n]=pres[n];
  }
#else
  /* for non mpi do nothing */
  return;
#endif
}



void global_sum_int(int *d, int n) {
#ifdef WITH_MPI
  int mpiret; (void)mpiret; // Remove warning of variable set but not used
  int pres[n];

#ifdef MPI_TIMING
  struct timeval start, end, etime;
  gettimeofday(&start,0);
#endif

  mpiret=MPI_Allreduce(d,pres,n,MPI_INT,MPI_SUM,GLB_COMM);


#ifdef MPI_TIMING
  gettimeofday(&end,0);
  timeval_subtract(&etime,&end,&start);
  lprintf("MPI TIMING",0,"global_sum_int " __FILE__ " %ld sec %ld usec\n",etime.tv_sec,etime.tv_usec);
#endif

#ifndef NDEBUG
  if (mpiret != MPI_SUCCESS) {
    char mesg[MPI_MAX_ERROR_STRING];
    int mesglen;
    MPI_Error_string(mpiret,mesg,&mesglen);
    lprintf("MPI",0,"ERROR: %s\n",mesg);
    error(1,1,"global_sum_int " __FILE__,"Cannot perform global_sum");
  }
#endif
  while(n>0) {
    --n;
    d[n]=pres[n];
  }
#else
  /* for non mpi do nothing */
  return;
#endif
}



void global_max(double *d, int n) {
#ifdef WITH_MPI
  int mpiret;(void)mpiret; // Remove warning of variable set but not used
  double pres[n];

#ifdef MPI_TIMING
  struct timeval start, end, etime;
  gettimeofday(&start,0);
#endif

  mpiret=MPI_Allreduce(d,pres,n,MPI_DOUBLE,MPI_MAX,GLB_COMM);


#ifdef MPI_TIMING
  gettimeofday(&end,0);
  timeval_subtract(&etime,&end,&start);
  lprintf("MPI TIMING",0,"global_max " __FILE__ " %ld sec %ld usec\n",etime.tv_sec,etime.tv_usec);
#endif

#ifndef NDEBUG
  if (mpiret != MPI_SUCCESS) {
    char mesg[MPI_MAX_ERROR_STRING];
    int mesglen;
    MPI_Error_string(mpiret,mesg,&mesglen);
    lprintf("MPI",0,"ERROR: %s\n",mesg);
    error(1,1,"global_max " __FILE__,"Cannot perform global_sum");
  }
#endif
  while(n>0) {
    --n;
    d[n]=pres[n];
  }
#else
  /* for non mpi do nothing */
  return;
#endif
}


void bcast(double *d, int n) {
#ifdef WITH_MPI
  int mpiret; (void)mpiret; // Remove warning of variable set but not used

#ifdef MPI_TIMING
  struct timeval start, end, etime;
  gettimeofday(&start,0);
#endif

  mpiret=MPI_Bcast(d, n, MPI_DOUBLE, 0,GLB_COMM);

#ifdef MPI_TIMING
  gettimeofday(&end,0);
  timeval_subtract(&etime,&end,&start);
  lprintf("MPI TIMING",0,"bcast " __FILE__ " %ld sec %ld usec\n",etime.tv_sec,etime.tv_usec);
#endif

#ifndef NDEBUG
  if (mpiret != MPI_SUCCESS) {
    char mesg[MPI_MAX_ERROR_STRING];
    int mesglen;
    MPI_Error_string(mpiret,mesg,&mesglen);
    lprintf("MPI",0,"ERROR: %s\n",mesg);
    error(1,1,"bcast " __FILE__,"Cannot perform global_sum");
  }
#endif

#else
  /* for non mpi do nothing */
  return;
#endif
}

#ifdef LLRHB
void bcast_from_rank(double *d, int n, int r) {
#ifdef WITH_MPI
  int mpiret; (void)mpiret; // Remove warning of variable set but not used

#ifdef MPI_TIMING
  struct timeval start, end, etime;
  gettimeofday(&start,0);
#endif //MPI_TIMING

  mpiret=MPI_Bcast(d, n, MPI_DOUBLE, r,GLB_COMM);

#ifdef MPI_TIMING
  gettimeofday(&end,0);
  timeval_subtract(&etime,&end,&start);
  lprintf("MPI TIMING",0,"bcast " __FILE__ " %ld sec %ld usec\n",etime.tv_sec,etime.tv_usec);
#endif //MPI_TIMING
#ifndef NDEBUG
  if (mpiret != MPI_SUCCESS) {
    char mesg[MPI_MAX_ERROR_STRING];
    int mesglen;
    MPI_Error_string(mpiret,mesg,&mesglen);
    lprintf("MPI",0,"ERROR: %s\n",mesg);
    error(1,1,"bcast " __FILE__,"Cannot perform global_sum");
  }
#endif //NDEBUG

#else //WITH_MPI
  /* for non mpi do nothing */
  return;
#endif //WITH_MPI
}
#endif //LLRHB

void bcast_int(int *i, int n) {
#ifdef WITH_MPI
  int mpiret; (void)mpiret; // Remove warning of variable set but not used

#ifdef MPI_TIMING
  struct timeval start, end, etime;
  gettimeofday(&start,0);
#endif

  mpiret=MPI_Bcast(i, n, MPI_INT, 0,GLB_COMM);

#ifdef MPI_TIMING
  gettimeofday(&end,0);
  timeval_subtract(&etime,&end,&start);
  lprintf("MPI TIMING",0,"bcast " __FILE__ " %ld sec %ld usec\n",etime.tv_sec,etime.tv_usec);
#endif

#ifndef NDEBUG
  if (mpiret != MPI_SUCCESS) {
    char mesg[MPI_MAX_ERROR_STRING];
    int mesglen;
    MPI_Error_string(mpiret,mesg,&mesglen);
    lprintf("MPI",0,"ERROR: %s\n",mesg);
    error(1,1,"bcast " __FILE__,"Cannot perform global_sum");
  }
#endif

#else
  /* for non mpi do nothing */
  return;
#endif
}

/* functions for filling sending buffer */
#ifdef WITH_MPI
static void sync_gauge_field(suNg_field *gf) {
  int i;
  geometry_descriptor *gd=gf->type;
  /* int j, mu, x, y; */

  for(i=0; i<gd->ncopies_gauge; ++i) {
    /* this assumes that the 4 directions are contiguous in memory !!! */
    memcpy(((gf->ptr)+4*gd->copy_to[i]),((gf->ptr)+4*gd->copy_from[i]),4*(gd->copy_len[i])*sizeof(*(gf->ptr)));
    /*
         for(j=0; j<gd->copy_len[i]; j++) {
         x=gd->copy_from[i]+j;
         y=gd->copy_to[i]+j;
         for(mu=0; mu<4; mu++)
     *pu_gauge(y,mu) = *pu_gauge(x,mu);
     }
     */
  }
}

static void sync_spinor_field(spinor_field *p) {
  int i;
  /* int j, x, y; */
  geometry_descriptor *gd = p->type;

  for(i=0; i<gd->ncopies_spinor; ++i) {
    memcpy((p->ptr+gd->copy_to[i]-gd->master_shift),(p->ptr+gd->copy_from[i]-gd->master_shift),(gd->copy_len[i])*sizeof(*(p->ptr)));
    /*
       for(j=0; j<gd->copy_len[i]; j++) {
       x=gd->copy_from[i]+j;
       y=gd->copy_to[i]+j;
       p->ptr[y] = p->ptr[x];
       }
       */
  }
}

static void sync_gauge_transf(suNg_field *gf) {
  int i;
  geometry_descriptor *gd=gf->type;

  for(i=0; i<gd->ncopies_gauge; ++i) {
    memcpy(((gf->ptr)+gd->copy_to[i]),((gf->ptr)+gd->copy_from[i]),(gd->copy_len[i])*sizeof(*(gf->ptr)));
  }
}

#endif /* WITH_MPI */

/* This variable contains the information of the current status of communications
 * Values:
 * 0 => No communications pending
 *
 */
/* static unsigned int comm_status=0; */

void complete_gf_sendrecv(suNg_field *gf) {
#ifdef WITH_MPI
  int mpiret; (void)mpiret; // Remove warning of variable set but not used
  int nreq=2*gf->type->nbuffers_gauge;

  if(nreq>0) {
    MPI_Status status[nreq];

    mpiret=MPI_Waitall(nreq, gf->comm_req, status);

#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
      char mesg[MPI_MAX_ERROR_STRING];
      int mesglen, k;
      MPI_Error_string(mpiret,mesg,&mesglen);
      lprintf("MPI",0,"ERROR: %s\n",mesg);
      for (k=0; k<nreq; ++k) {
        if (status[k].MPI_ERROR != MPI_SUCCESS) {
          MPI_Error_string(status[k].MPI_ERROR,mesg,&mesglen);
          lprintf("MPI",0,"Req [%d] Source [%d] Tag [%] ERROR: %s\n",
              k,
              status[k].MPI_SOURCE,
              status[k].MPI_TAG,
              mesg);
        }
      }
      error(1,1,"complete_gf_sendrecv " __FILE__,"Cannot complete communications");
    }
#endif
  }

#ifdef MPI_TIMING
    if(gf_control>0)
      {
	gettimeofday(&gfend,0);
	timeval_subtract(&gfetime,&gfend,&gfstart);
	lprintf("MPI TIMING",0,"complete_gf_sendrecv" __FILE__ " %ld sec %ld usec\n",gfetime.tv_sec,gfetime.tv_usec);
	gf_control=0;
      }
#endif

#endif /* WITH_MPI */
}

void start_gf_sendrecv(suNg_field *gf) {
#ifdef WITH_MPI
  int i, mpiret;(void)mpiret; // Remove warning of variable set but not used
  geometry_descriptor *gd=gf->type;

  /* check communication status */
  /* questo credo che non sia il modo piu' efficiente!!! */
  /* bisognerebbe forse avere una variabile di stato nei campi?? */
  complete_gf_sendrecv(gf);


  /* fill send buffers */
  sync_gauge_field(gf);

#ifdef MPI_TIMING
  error(gf_control>0,1,"start_gf_sendrecv " __FILE__,"Multiple send without receive");
  gettimeofday(&gfstart,0);
  gf_control=1;
#endif

  for (i=0; i<(gd->nbuffers_gauge); ++i) {
    /* send ith buffer */
    mpiret=MPI_Isend((gf->ptr)+4*gd->sbuf_start[i], /* buffer */
        (gd->sbuf_len[i])*sizeof(suNg)/sizeof(double)*4, /* lenght in units of doubles */
        MPI_DOUBLE, /* basic datatype */
        gd->sbuf_to_proc[i], /* cid of destination */
        i, /* tag of communication */
        cart_comm, /* use the cartesian communicator */
        &(gf->comm_req[2*i]) /* handle to communication request */
        );
#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
      char mesg[MPI_MAX_ERROR_STRING];
      int mesglen;
      MPI_Error_string(mpiret,mesg,&mesglen);
      lprintf("MPI",0,"ERROR: %s\n",mesg);
      error(1,1,"start_gf_sendrecv " __FILE__,"Cannot start send buffer");
    }
#endif

    /* receive ith buffer */
    mpiret=MPI_Irecv((gf->ptr)+4*gd->rbuf_start[i], /* buffer */
        (gd->rbuf_len[i])*sizeof(suNg)/sizeof(double)*4, /* lenght in units of doubles */
        MPI_DOUBLE, /* basic datatype */
        gd->rbuf_from_proc[i], /* cid of origin */
        i, /* tag of communication */
        cart_comm, /* use the cartesian communicator */
        &(gf->comm_req[2*i+1]) /* handle to communication request */
        );
#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
      char mesg[MPI_MAX_ERROR_STRING];
      int mesglen;
      MPI_Error_string(mpiret,mesg,&mesglen);
      lprintf("MPI",0,"ERROR: %s\n",mesg);
      error(1,1,"start_gf_sendrecv " __FILE__,"Cannot start receive buffer");
    }
#endif


    /*
      if shifted
      force complete field transfer.

      shift field

      (for test write 0 in the remaining missed buffer)


(0,0,0,0) <= (t,x,y,z) < (T+2*T_BORDER,X+2*X_BORDER,Y+2*Y_BORDER,Z+2*Z_BORDER)


n[i] e' sul reticolo esteso

if(X_BORDER != 0 || Y_BORDER != 0 || Z_BORDER != 0)

for(n[1]=0;n[1]<X+2*X_BORDER;n[1]++)
for(n[2]=0;n[2]<Y+2*Y_BORDER;n[2]++)
for(n[3]=0;n[3]<Z+2*Z_BORDER;n[3]++)
n[0]=0

sono fuori se
n[1]+SHIFT_X==X+2*X_BORDER && X_BORDER!=0
or
n[2]+SHIFT_Y==Y+2*Y_BORDER && Y_BORDER!=0
or
n[3]+SHIFT_Z==Z+2*Z_BORDER && Z_BORDER!=0

shift_x*shift_y*shift_z*X_BORDER*Y_BORDER*Z_BORDER

shift_X*(X-1) shift_Y*shift_Z*Y_BORDER*Z_BORDER*(X-1)
shift_Y*(Y-1) shift_X*shift_Z*X_BORDER*Z_BORDER*(Y-1)
shift_Z*(Z-1) shift_X*shift_Y*X_BORDER*Y_BORDER(Z-1)

((Y-1)*(Z-1))*shift_X*BORDER_X



init shifted

n[0]=0
for(n[1]=0;n[1]<X+2*X_BORDER;n[1]++)
for(n[2]=0;n[2]<Y+2*Y_BORDER;n[2]++)
for(n[3]=0;n[3]<Z+2*Z_BORDER;n[3]++)
if(n[1]+shift_X==X+2*X_BORDER || n[2]+shift_Y==Y+2*Y_BORDER || n[3]+shift_Z==Z+2*Z_BORDER)
next=id
if(n[1]+shift_X==X+2*X_BORDER ) next=proc_up[next,1];
if(n[2]+shift_Y==Y+2*Y_BORDER ) next=proc_up[next,2];
if(n[3]+shift_Z==Z+2*Z_BORDER ) next=proc_up[next,3];

if(next!=id)
add_shift_mem[next]++;
map(n[1],n[2],n[3])=add_shift_mem[next]





      then send and complete again

*/

  }

#endif /* WITH_MPI */
}

void complete_sf_sendrecv(spinor_field *sf) {
#ifdef WITH_MPI
  int mpiret; (void)mpiret; // Remove warning of variable set but not used
  int nreq=2*sf->type->nbuffers_spinor;

  if(nreq>0) {
    MPI_Status status[nreq];

    mpiret=MPI_Waitall(nreq, sf->comm_req, status);

#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
      char mesg[MPI_MAX_ERROR_STRING];
      int mesglen, k;
      MPI_Error_string(mpiret,mesg,&mesglen);
      lprintf("MPI",0,"ERROR: %s\n",mesg);
      for (k=0; k<nreq; ++k) {
        if (status[k].MPI_ERROR != MPI_SUCCESS) {
          MPI_Error_string(status[k].MPI_ERROR,mesg,&mesglen);
          lprintf("MPI",0,"Req [%d] Source [%d] Tag [%] ERROR: %s\n",
              k,
              status[k].MPI_SOURCE,
              status[k].MPI_TAG,
              mesg);
        }
      }
      error(1,1,"complete_gf_sendrecv " __FILE__,"Cannot complete communications");
    }
#endif
  }

#ifdef MPI_TIMING
    if(sf_control>0)
      {
	gettimeofday(&sfend,0);
	timeval_subtract(&sfetime,&sfend,&sfstart);
	lprintf("MPI TIMING",0,"complete_sf_sendrecv" __FILE__ " %ld sec %ld usec\n",sfetime.tv_sec,sfetime.tv_usec);
	sf_control=0;
      }
#endif

#endif /* WITH_MPI */
}

void start_sf_sendrecv(spinor_field *sf) {
#ifdef WITH_MPI
  int i, mpiret; (void)mpiret; // Remove warning of variable set but not used
  geometry_descriptor *gd=sf->type;


  /* check communication status */
  /* questo credo che non sia il modo piu' efficiente!!! */
  /* bisognerebbe forse avere una variabile di stato nei campi?? */
  complete_sf_sendrecv(sf);

  /* fill send buffers */
  sync_spinor_field(sf);
#ifdef MPI_TIMING
  error(sf_control>0,1,"start_sf_sendrecv " __FILE__,"Multiple send without receive");
  gettimeofday(&sfstart,0);
  sf_control=1;
#endif

  for (i=0; i<(gd->nbuffers_spinor); ++i) {
    /* send ith buffer */
    mpiret=MPI_Isend((sf->ptr)+(gd->sbuf_start[i])-(gd->master_shift), /* buffer */
        (gd->sbuf_len[i])*(sizeof(suNf_spinor)/sizeof(double)), /* lenght in units of doubles */
        MPI_DOUBLE, /* basic datatype */
        gd->sbuf_to_proc[i], /* cid of destination */
        i, /* tag of communication */
        cart_comm, /* use the cartesian communicator */
        &(sf->comm_req[2*i]) /* handle to communication request */
        );
#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
      char mesg[MPI_MAX_ERROR_STRING];
      int mesglen;
      MPI_Error_string(mpiret,mesg,&mesglen);
      lprintf("MPI",0,"ERROR: %s\n",mesg);
      error(1,1,"start_gf_sendrecv " __FILE__,"Cannot start send buffer");
    }
#endif

    /* receive ith buffer */
    mpiret=MPI_Irecv((sf->ptr)+(gd->rbuf_start[i])-(gd->master_shift), /* buffer */
        (gd->rbuf_len[i])*(sizeof(suNf_spinor)/sizeof(double)), /* lenght in units of doubles */
        MPI_DOUBLE, /* basic datatype */
        gd->rbuf_from_proc[i], /* cid of origin */
        i, /* tag of communication */
        cart_comm, /* use the cartesian communicator */
        &(sf->comm_req[2*i+1]) /* handle to communication request */
        );
#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
      char mesg[MPI_MAX_ERROR_STRING];
      int mesglen;
      MPI_Error_string(mpiret,mesg,&mesglen);
      lprintf("MPI",0,"ERROR: %s\n",mesg);
      error(1,1,"start_gf_sendrecv " __FILE__,"Cannot start receive buffer");
    }
#endif

  }

#endif /* WITH_MPI */
}


void complete_gt_sendrecv(suNg_field *gf) {
#ifdef WITH_MPI
  int mpiret; (void)mpiret; // Remove warning of variable set but not used
  int nreq=2*gf->type->nbuffers_gauge;

  if(nreq>0) {
    MPI_Status status[nreq];

    mpiret=MPI_Waitall(nreq, gf->comm_req, status);

#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
      char mesg[MPI_MAX_ERROR_STRING];
      int mesglen, k;
      MPI_Error_string(mpiret,mesg,&mesglen);
      lprintf("MPI",0,"ERROR: %s\n",mesg);
      for (k=0; k<nreq; ++k) {
        if (status[k].MPI_ERROR != MPI_SUCCESS) {
          MPI_Error_string(status[k].MPI_ERROR,mesg,&mesglen);
          lprintf("MPI",0,"Req [%d] Source [%d] Tag [%] ERROR: %s\n",
              k,
              status[k].MPI_SOURCE,
              status[k].MPI_TAG,
              mesg);
        }
      }
      error(1,1,"complete_gt_sendrecv " __FILE__,"Cannot complete communications");
    }
#endif
  }

#endif /* WITH_MPI */
}

void start_gt_sendrecv(suNg_field *gf) {
#ifdef WITH_MPI
  int i, mpiret;
  (void)mpiret; // Remove warning of variable set but not used
  geometry_descriptor *gd=gf->type;

  /* check communication status */
  /* questo credo che non sia il modo piu' efficiente!!! */
  /* bisognerebbe forse avere una variabile di stato nei campi?? */
  complete_gt_sendrecv(gf);

  /* fill send buffers */
  sync_gauge_transf(gf);

  for (i=0; i<(gd->nbuffers_gauge); ++i) {
    /* send ith buffer */
    mpiret=MPI_Isend((gf->ptr)+gd->sbuf_start[i], /* buffer */
        (gd->sbuf_len[i])*sizeof(suNg)/sizeof(double), /* lenght in units of doubles */
        MPI_DOUBLE, /* basic datatype */
        gd->sbuf_to_proc[i], /* cid of destination */
        i, /* tag of communication */
        cart_comm, /* use the cartesian communicator */
        &(gf->comm_req[2*i]) /* handle to communication request */
        );
#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
      char mesg[MPI_MAX_ERROR_STRING];
      int mesglen;
      MPI_Error_string(mpiret,mesg,&mesglen);
      lprintf("MPI",0,"ERROR: %s\n",mesg);
      error(1,1,"start_gt_sendrecv " __FILE__,"Cannot start send buffer");
    }
#endif

    /* receive ith buffer */
    mpiret=MPI_Irecv((gf->ptr)+gd->rbuf_start[i], /* buffer */
        (gd->rbuf_len[i])*sizeof(suNg)/sizeof(double), /* lenght in units of doubles */
        MPI_DOUBLE, /* basic datatype */
        gd->rbuf_from_proc[i], /* cid of origin */
        i, /* tag of communication */
        cart_comm, /* use the cartesian communicator */
        &(gf->comm_req[2*i+1]) /* handle to communication request */
        );
#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
      char mesg[MPI_MAX_ERROR_STRING];
      int mesglen;
      MPI_Error_string(mpiret,mesg,&mesglen);
      lprintf("MPI",0,"ERROR: %s\n",mesg);
      error(1,1,"start_gt_sendrecv " __FILE__,"Cannot start receive buffer");
    }
#endif

  }

#endif /* WITH_MPI */
}


#ifdef WITH_UMBRELLA
void umbrella_swap(double* S_llr,double* S0, double* a, double* dS)
{


  int mpiret; (void)mpiret;

  lprintf("SWAP",10,"Starting Rep Par S0 = %f dS = %f a = %f \n",*S0,*dS,*a);

  /*wait on every processor*/
  mpiret=MPI_Barrier(MPI_COMM_WORLD);
  double data[4*N_REP];
  double locdata[4];

  locdata[0]=*S_llr;
  locdata[1]=*S0;
  locdata[2]=*a;
  locdata[3]=*dS;
  if(PID==0) {
    mpiret=MPI_Gather(locdata,4,MPI_DOUBLE,data,4,MPI_DOUBLE,0,UMB_WORLD);
#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
      char mesg[MPI_MAX_ERROR_STRING];
      int mesglen;
      MPI_Error_string(mpiret,mesg,&mesglen);
      lprintf("MPI",0,"ERROR: %s\n",mesg);
      error(1,1,"umbrella_swap " __FILE__,"Cannot complete gather");
    }
#endif
    if(UID==0) swap(data);
    mpiret=MPI_Scatter(data,4,MPI_DOUBLE,locdata,4,MPI_DOUBLE,0,UMB_WORLD);
#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
      char mesg[MPI_MAX_ERROR_STRING];
      int mesglen;
      MPI_Error_string(mpiret,mesg,&mesglen);
      lprintf("MPI",0,"ERROR: %s\n",mesg);
      error(1,1,"umbrella_swap " __FILE__,"Cannot complete scatter");
    }
#endif

  }
  bcast(locdata,4);
  setreplica(locdata);

  lprintf("SWAP",10,"New Rep Par S0 = %f dS = %f a = %f \n",*S0,*dS,*a);

}
#endif //WITH_UMBRELLA
