/*
** Name:       utils.h
** Author:     Leo Liberti
** Source:     GNU C++
** Purpose:    utilities for rose - headers
** License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
** History:    050429 work started
*/

#ifndef ROSEUTILSH
#define ROSEUTILSH

#include <vector>
#include "parameters.h"

double getcputime(double* UserSysTime);
double boundedrandom(double LB, double UB);
void randomvector(std::vector<double>& x, double maxval);
void randomvectorinbox(std::vector<double>& x, 
		       std::vector<double>& lb,
		       std::vector<double>& ub);
double L2distance(std::vector<double>&, std::vector<double>&);
double L2norm(std::vector<double>& x);
double Linfnorm(std::vector<double>& x);
int parmtype(char* parmval, bool& boolval, int& intval, double& doubleval,
	     std::string& stringval);

void permutation3old(int **ind);
void permutation3(int **ind,int *ibnd);

int gcd(int a, int b);
std::string float2fraction(double a);
#endif
