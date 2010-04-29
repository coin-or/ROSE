/*
** Name:       utils.cxx
** Author:     Leo Liberti
** Source:     GNU C++
** Purpose:    utilities for rose
** License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
** History:    050429 work started
*/

#include "utils.h"

#include <iostream>
#include <cassert>
#include <cstdlib>
#include <sys/time.h> 

/** Utility function for getting the CPU time
    @return the time
    @param UserSysTime: struct (see sys/resource.h)
*/
#include <sys/resource.h>
#include <cstdio>
double getcputime(double* UserSysTime) {
  struct rusage ru;
  if (getrusage(RUSAGE_SELF, &ru) == -1) {
    perror("ERROR: getrusage");
    return -1;
  } else {
    UserSysTime[0] = (double) ru.ru_utime.tv_sec
        + ((double) ru.ru_utime.tv_usec) / 1.0e6;
    UserSysTime[1] = (double) ru.ru_stime.tv_sec
        + ((double) ru.ru_stime.tv_usec) / 1.0e6;
    return (UserSysTime[0] + UserSysTime[1]);
  }
}

/** Utility function for generating a bounded random double
    @return the random number
    @param LB: lower bound
    @param UB: upper bound
*/
double boundedrandom(double LB, double UB) {
  assert (UB >= LB);
  double rc = (double) random() / (double) RAND_MAX;
  rc *= (UB - LB);
  rc += LB;
  return rc;
}

/** Utility function for generating a bounded random vector
    @param x the random vector (initialized with correct size)
    @param maxval the |bound| for vector values
*/
void randomvector(std::vector<double>& x, double maxval) {
  double rc;
  int n = x.size();
  assert(n > 0);
  for(int i = 0; i < n; i++) {
    rc = (double) random() / (double) RAND_MAX;
    rc = 2*rc - 1;
    rc *= maxval;
    x[i] = rc;
  }
}

/** Utility function for generating a bounded random vector in a box
    @param x the random vector (initialized with correct size)
    @param lb the box lower bounds
    @param ub the box upper bounds
*/
void randomvectorinbox(std::vector<double>& x, 
		       std::vector<double>& lb,
		       std::vector<double>& ub) {
  int n = x.size();
  assert(n > 0);
  assert(lb.size() == n);
  assert(ub.size() == n);
  // sample a random point in lb < x < ub
  double rc;
  double rv;
  double top;
  for(int i = 0; i < n; i++) {
    rc = (double) random() / (double) RAND_MAX;
    top = ub[i] - lb[i];
    rv = rc * top + lb[i];
    x[i] = rv;
  }
}  

/** Utility function for calculating an L2 distance
    @return the distance
    @param x: the first vector
    @param y: the second vector
*/
#include <vector>
#include <cmath>
double L2distance(std::vector<double>& x, std::vector<double>& y) {
  int n = x.size();
  assert(n == y.size());
  double val = 0;
  for(int i = 0; i < n; i++) {
    val += (x[i] - y[i])*(x[i] - y[i]);
  }
  val = sqrt(val);
  return val;
}

/** Utility function for calculating an L2 norm
    @return the distance
    @param x: the first vector
    @param y: the second vector
*/
double L2norm(std::vector<double>& x) {
  int n = x.size();
  double val = 0;
  for(int i = 0; i < n; i++) {
    val += x[i]*x[i];
  }
  val = sqrt(val);
  return val;
}

/** Utility function for calculating an Linf norm
    @return the distance
    @param x: the first vector
    @param y: the second vector
*/
double Linfnorm(std::vector<double>& x) {
  int n = x.size();
  double val = 0;
  for(int i = 0; i < n; i++) {
    if (fabs(x[i]) > val) {
      val = fabs(x[i]);
    }
  }
  return val;
}

/** Utility function for finding the type of a parameter value held in a string
    @return the type ({ IntType, BoolType, DoubleType, StringType })
    @param parmval the string parameter
    @param boolval the output boolean value
    @param intval the output integer value
    @param doubleval the output double value
*/
#include <string>
int parmtype(char* parmval, bool& boolval, int& intval, double& doubleval,
	     std::string& stringval) {

  using namespace std;

  char* ep;
  int paramtype;
  
  // is it an integer?
  intval = strtol(parmval, &ep, 0);
  paramtype = IntType;
  if (*ep != 0) {
    // no: is it a double?
    doubleval = strtod(parmval, &ep);
    paramtype = DoubleType;
    if (*ep != 0) {
      // no: it's a string		  
      stringval = parmval;
      paramtype = StringType;
      // test if bool
      if (stringval == "true" || stringval == "yes" ||
	  stringval == "TRUE" || stringval == "YES") {
	boolval = true;
	paramtype = BoolType;
      } else  if (stringval == "false" || stringval == "no" ||
		  stringval == "FALSE" || stringval == "NO") {
	boolval = false;
	paramtype = BoolType;
      } 
    } else {
      intval = (int) doubleval;
    }
  } else {
    // just in case we got the wrong type, initialize the compatible
    // types too
    doubleval = (double) intval;
    if (intval == 0 || intval == 1) {
      boolval = (bool) intval;
    }
  }

  return paramtype;

}

//// permutations of 3 elements

void permutation3old(int **ind)
{
  ind[0][0] = 0; ind[0][1] = 1; ind[0][2] = 2;
  ind[1][0] = 0; ind[1][1] = 2; ind[1][2] = 1;
  ind[2][0] = 1; ind[2][1] = 0; ind[2][2] = 2;
  ind[3][0] = 1; ind[3][1] = 2; ind[3][2] = 0;
  ind[4][0] = 2; ind[4][1] = 0; ind[4][2] = 1;
  ind[5][0] = 2; ind[5][1] = 1; ind[5][2] = 0;
}
void permutation3(int **ind,int *ibnd)
{
  ind[0][0] = ibnd[0]; ind[0][1] = ibnd[1]; ind[0][2] = ibnd[2];
  ind[1][0] = ibnd[0]; ind[1][1] = ibnd[2]; ind[1][2] = ibnd[1];
  ind[2][0] = ibnd[1]; ind[2][1] = ibnd[0]; ind[2][2] = ibnd[2];
  ind[3][0] = ibnd[1]; ind[3][1] = ibnd[2]; ind[3][2] = ibnd[0];
  ind[4][0] = ibnd[2]; ind[4][1] = ibnd[0]; ind[4][2] = ibnd[1];
  ind[5][0] = ibnd[2]; ind[5][1] = ibnd[1]; ind[5][2] = ibnd[0];
}



