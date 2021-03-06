#ifndef GUARD_bd_h
#define GUARD_bd_h

#include "info.h"
#include "tree.h"

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>


bool bd(std::vector<std::vector<double> >& X, tree& x, xinfo& xi, dinfo& di, pinfo& pi, size_t minobsnode, int *binaryX);

bool bd1(std::vector<std::vector<double> >& X, tree& x, xinfo& xi, dinfo& di, pinfo& pi, size_t minobsnode, int *binaryX, double* nLtD, double* pA, double* nN, double* nL, double* tD, std::vector<double>& iP);

#endif
