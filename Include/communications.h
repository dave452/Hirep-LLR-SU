/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/
/**************************************
Modified 2021 by David Mason and Davide Vadacchino 
***************************************/
#ifndef COMMUNICATIONS_H
#define COMMUNICATIONS_H

void global_sum(double *d, int n);
void global_sum_int(int *d, int n);
void global_max(double *d, int n);
void bcast(double *d, int n);
void bcast_int(int *i, int n);

#include "spinor_field.h"
void complete_gf_sendrecv(suNg_field *gf);
void start_gf_sendrecv(suNg_field *gf);
void complete_sf_sendrecv(spinor_field *gf);
void start_sf_sendrecv(spinor_field *gf);

void complete_gt_sendrecv(suNg_field *gf);
void start_gt_sendrecv(suNg_field *gf);

void test_spinor_field(spinor_field *p);

/* Floating point sendrecv */
void complete_gf_sendrecv_flt(suNg_field_flt *gf);
void start_gf_sendrecv_flt(suNg_field_flt *gf);
void complete_sf_sendrecv_flt(spinor_field_flt *gf);
void start_sf_sendrecv_flt(spinor_field_flt *gf);
#ifdef WITH_UMBRELLA
#ifdef LLRHB
void bcast_from_rank(double *d, int n, int r);
#endif
#ifdef LLRHB_UM_BC
void umbrella_swap(double* S_llr,double* S0, double* a, double* dS);
void umbrella_swap_hb(double* S_llr,double* S0, double* a, double* dS);
#else
void umbrella_swap(double* S_llr,double* S0, double* a, double* dS, double* var_llr);
#endif
#endif

#endif /* COMMUNICATIONS_H */
