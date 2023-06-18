#ifndef GLOBAL_H
#define GLOBAL_H

#include "genotype.h"
#include "Goptions.hpp"
#include "pheno.hpp"

genotype g;
Goptions goptions;
pheno phe;

double **partialsums;
double *sum_op;

// Intermediate computations in E-step.
// Size = 3^(log_3(n)) * K
double **yint_e;
// n X K
double ***y_e;

// Intermediate computations in M-step. 
// Size = nthreads X 3^(log_3(n)) * K
double **yint_m;
// nthreads X log_3(n) X k
double ***y_m;

#endif
