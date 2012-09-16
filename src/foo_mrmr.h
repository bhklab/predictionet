#include <iostream>
#include <string>
#include <math.h> 
using namespace std;
#include <algorithm>
#include <ctime>
#include <cmath>
#include <vector>
#include <map>
#include "tree.h"
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

/* Entry points called from the R functions */
extern "C" 
{
	SEXP mrnet_adapted(SEXP data, SEXP namat, SEXP maxparents, SEXP nvar, SEXP nsample, SEXP predn, SEXP npredn, SEXP threshold);
	SEXP mrnet_adapted2(SEXP data, SEXP namat, SEXP prior, SEXP prior_weight, SEXP maxparents, SEXP nvar, SEXP nsample, SEXP predn, SEXP npredn, SEXP threshold);
	SEXP mrmr_ensemble(SEXP data, SEXP namat, SEXP maxparents, SEXP nvar, SEXP nsample, SEXP predn, SEXP npredn, SEXP rep_boot, SEXP maxnsol, SEXP threshold);
	SEXP mrmr_ensemble_remove(SEXP data, SEXP namat, SEXP maxparents, SEXP nvar, SEXP nsample, SEXP predn, SEXP npredn, SEXP rep_boot, SEXP maxnsol, SEXP threshold);

}

double get_correlation(double data [],int namat[],int ind_x, int ind_y, int size);
void build_mim_subset(double mim[],double data[], int namat [],int nvar,int nsample, int subset [],int size_subset);
int verify_equivalentset (tree<int>& tr, tree<int>::pre_order_iterator it, tree<int>::pre_order_iterator end, int maxnsol, int order_addition[], int res_vec[]);
double mrnet_onegene(double mim [], int size, int nbvar,int *var_ind,int var_target, int var_interest);
