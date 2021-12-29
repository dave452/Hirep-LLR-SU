/***************************************************************************\
* Copyright (c) 2008, Agostino Patella, Claudio Pica                        *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef UPDATE_H
#define UPDATE_H

#include "suN.h"
#include "inverters.h"
#include "rational_functions.h"

/*Functions added for constrained heathbath 
 * and the related implementation of the LLR*/

void random_su2_creutz(double rho,double s[]);
void random_su2_constrained(double rho,double s[], double xmin, double xmax);

void cabmar_constrained(double beta,suNg *u, suNg *v,int type, double * E, double Emin, double Emax);
void update_constrained(double beta,int nhb,int nor, double *S, double Smin, double Smax);
//void update_constrained(double beta,int nhb,int nor);
void anneal(double *E, double S0, double dS);

void thermrobbinsmonro_hb(void);
void measrobbinsmonro_hb(void);
void robbinsmonro_hb(void);
void restart_robbinsmonro_hb();
void init_robbinsmonro_hb(int nrm,int nth,double starta,int it,double dS,double S0, int sfreq_RM, int sfreq_fxa);
/*================*/


void staples(int ix,int mu,suNg *v);
void test_staples();
void spatial_staples(int ix,int mu,suNg *v);
void cabmar(double beta,suNg *u, suNg *v,int type);

void project_gauge_field(void);

void update(double beta,int nhb,int nor);
void random_su2(double rho,double s[]);

/* functions and structures for the MRE algorithm */
typedef struct {
	spinor_field *s[2];
	int num[2];
	int max;
	int init;
} mre_par;

void mre_guess(mre_par*, int, spinor_field*, spinor_operator, spinor_field*);
void mre_store(mre_par*, int, spinor_field*);
void mre_init(mre_par*, int, double);

/* forces for the update */
void force0(double dt, suNg_av_field *force, void *par);
void force_obs_0pp(double dt, suNg_av_field *force, void *vpar);
typedef struct {
  int n_pf;
  spinor_field *pf;
  double mass;
  rational_app *ratio;
  double inv_err2;
} force_rhmc_par;

void init_force_rhmc();
void free_force_rhmc();
void force_rhmc(double dt, suNg_av_field *force, void *par);

typedef struct {
  int id;
  int n_pf;
  spinor_field *pf;
  int hasenbusch;
  double mass;
  double b;
  double mu;
  double inv_err2, inv_err2_flt;
  mre_par mpar;
} force_hmc_par;

typedef struct {
  int t0,t1;
  double shift;
} force_obs_0pp_par;



void force_fermion_core(spinor_field* Xs, spinor_field* Ys, suNg_av_field* force, double dt, double* forcestat, int type);
void force_hmc(double dt, suNg_av_field *force, void *par);
void force_hmc_tm(double dt, suNg_av_field *force, void *par);


void gaussian_momenta(suNg_av_field *momenta);
void gaussian_spinor_field(spinor_field *s);
void gaussian_spinor_field_flt(spinor_field_flt *s);
void z2_spinor_field(spinor_field *s);


/* For the fermion force ? */
void corret_pf_dist_hmc();
void calc_one_force(int n_force);


/* Action structures */
typedef enum {
  PureGauge,
  HMC,
  RHMC,
  TM,
  TM_alt,
  Hasenbusch,
  Hasenbusch_tm,
  Hasenbusch_tm_alt,
  LLRPureGauge,
  LLR_obs_0pp,
  LLR_HMC,
  NUM_MON_TYPE
} mon_type;

#ifdef MAIN_PROGRAM
int is_llr[NUM_MON_TYPE]={
0,//PureGauge
0,//HMC
0,//RHMC
0,//TM
0,//TM_alt
0,//Hasenbusch
0,//Hasenbusch_tm
0,//Hasenbusch_tm_alt
0,//LLRPureGauge
0,//LLR_obs_0pp
0};//LLR_HMC
#else
extern int is_llr[NUM_MON_TYPE];
#endif


typedef struct _mon_obs_0pp_par{
  force_obs_0pp_par fpar;
  int t0,t1;
  double shift;
} mon_obs_0pp_par;


typedef struct _mon_pg_par {
  double beta;
} mon_pg_par;

typedef struct _mon_hmc_par {
  double mass;
  int mre_past;
  force_hmc_par fpar;
  spinor_field *pf; /* pseudofermion field */
} mon_hmc_par;

typedef struct _mon_rhmc_par {
  double mass;
  rational_app ratio;
  force_rhmc_par fpar;
  spinor_field *pf; /* pseudofermion field */
} mon_rhmc_par;

typedef struct _mon_tm_par {
  double mass;
  double mu;
  int mre_past;
  force_hmc_par fpar;
  spinor_field *pf; /* pseudofermion field */
} mon_tm_par;


typedef struct _mon_hasenbusch_par {
  double mass;
  double dm;
  int mre_past;
  force_hmc_par fpar;
  spinor_field *pf; /* pseudofermion field */
} mon_hasenbusch_par;

typedef struct _mon_hasenbusch_tm_par {
  double mass;
  double mu;
  double dmu;
  int mre_past;
  force_hmc_par fpar;
  spinor_field *pf; /* pseudofermion field */
} mon_hasenbusch_tm_par;


typedef struct _monomial_data {
  int id; /* monomial id */
  mon_type type; /* type of monomial */
  void *par; /* parameters */
  double MT_prec; /* metropolis precision */
  double MD_prec; /* molecular dynamics precision */
  double force_prec; /* force precision */
} monomial_data;

typedef struct _monomial {
  monomial_data data;
  
  /* Functions */
  void (*free)(struct _monomial *m); /* free memory */
  
  void (*force_f)(double dt, suNg_av_field *force, void *par); /* force function */
  void *force_par; /* parameters for the force function */

  void (*init_traj)(const struct _monomial *m);
  void (*gaussian_pf)(const struct _monomial *m);
  void (*correct_pf)(const struct _monomial *m);
  void (*correct_la_pf)(const struct _monomial *m);
  const spinor_field *(*pseudofermion)(const struct _monomial *m); /* returns ps field pointer */
  void (*add_local_action)(const struct _monomial *m, scalar_field *loc_action);
  void (*add_llr_local_action)(const struct _monomial *m, scalar_field *loc_action);
  
} monomial;

struct _monomial* llr_obs_0pp_create(const monomial_data *data);
struct _monomial* llr_gauge_create(const monomial_data *data);
struct _monomial* llr_hmc_create(const monomial_data *data);
struct _monomial* pg_create(const monomial_data *data);
struct _monomial* hmc_create(const monomial_data *data);
struct _monomial* rhmc_create(const monomial_data *data);
struct _monomial* tm_create(const monomial_data *data);
struct _monomial* tm_alt_create(const monomial_data *data);
struct _monomial* hasen_create(const monomial_data *data);
struct _monomial* hasen_tm_create(const monomial_data *data);
struct _monomial* hasen_tm_alt_create(const monomial_data *data);

const monomial *add_mon(monomial_data *mon);
int num_mon();
const monomial *mon_n(int i);


typedef struct _integrator_par {
  int nsteps;
  int nmon;
  const monomial **mon_list;
  void (*integrator)(suNg_av_field*, double, struct _integrator_par*);
  struct _integrator_par *next;
  int level;
} integrator_par;

void gauge_integrator(suNg_av_field *momenta, double tlen, integrator_par *int_par);
void leapfrog_multistep(suNg_av_field *momenta, double tlen, integrator_par *int_par);
void O2MN_multistep(suNg_av_field *momenta, double tlen, integrator_par *int_par);
void O4MN_multistep(suNg_av_field *momenta, double tlen, integrator_par *int_par);


typedef struct _ghmc_par {
  
  /* integrator */
  integrator_par *integrator;
  double tlen;

  /* Fermion Theta angles */
  double theta[4];
  
  /* Probably not needed anymore */
  /* SF stuff */
  double SF_zf;
  double SF_ds;
  int SF_sign;
  double SF_ct;
  
} ghmc_par;

void init_ghmc(ghmc_par *par);
void free_ghmc();
void setstep();
int update_ghmc();
int update_ghmc_adapt();
int update_llr_ghmc(double *S_llr,double *S_non_llr,int therm);



/* local action */
typedef enum {
   NEW=1,
   DELTA=2
} local_action_type;

/*
 * compute the local action at every site for the HMC
 * H = | momenta |^2 + S_g + < phi1, phi2>
 */
void local_hmc_action(local_action_type type,
                      scalar_field *loc_action,
                      suNg_av_field *momenta);
void pf_local_action(scalar_field *loc_action,
                     spinor_field *pf);

void local_llr_hmc_action(double * S, double * S_llr, double * mom,scalar_field *loc_action,scalar_field *loc_llr_action,suNg_av_field *momenta);

void suNg_field_copy(suNg_field *g1, suNg_field *g2);
void suNf_field_copy(suNf_field *g1, suNf_field *g2);

/* find spectral interval using eva */
void find_spec_H2(double *max, double *min);

/*LLR Monomials*/

//int *llr_mons_idx;
//int *nonllr_mons_idx;
//void init_llr_mon();
//int num_llr_mon();

/* ROBBINS MONRO*/

void total_llr_action(double * S_llr);
void thermrobbinsmonro(void);
void measrobbinsmonro(void);
void robbinsmonro(void);
void restart_robbinsmonro();
void init_robbinsmonro(int nrm,int nth,double starta,int it,double dS,double S0, int sfreq_fxa, double Smin, double Smax, int nhb, int nor);
double getdS(void);
double get_llr_a(void);
double get_llr_a_hb(void);
double getS0(void);
double getS0_hb(void);
void llr_fixed_a_update(void);

#ifdef WITH_UMBRELLA
void swap(double *data);
//void swap_hb(double *data);
void setreplica(double *data);
//void setreplica_hb(double *data);
#endif


#endif
