/**********************************************************************
* Author:      Leo Liberti                                            *
* Name:        exprauxiliary.h                                        *
* Source:      GNU C++                                                *
* Purpose:     symbolic expression (auxiliary function header)        *
* History:     080928 0.0 work started                                *
* License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
***********************************************************************/

#ifndef __EV3EXPRAUXILIARYH__
#define __EV3EXPRAUXILIARYH__

bool is_integer(double a);
bool is_even(double a);
bool is_odd(double a);
// returns a very small positive value
double Ev3NearZero(void);
// returns a very large positive value
double Ev3Infinity(void);

double argmin(double a1, double a2, double a3, double a4);
double argmax(double a1, double a2, double a3, double a4);

// used by Rsmith and similar
void bilinearprodmkrange(double a, double b, double c, double d, 
			 double* wl, double *wu);
void fractionmkrange(double a, double b, double c, double d,
		     double* wl, double* wu);
void constpowermkrange(double a, double b, double c,
		       double* wl, double* wu);
void powermkrange(double a, double b, double c, double d,
		  double* wl, double* wu);

// used by Rsymmgroup
int computecolor(int symm, int level, int op, int cnst, int var, int ord, 
		 int all_lev_ops, int all_ops, int all_cnsts, 
		 int all_vars, int all_ords);

#endif
