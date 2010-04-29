/**********************************************************************
* Author:      Leo Liberti                                            *
* Name:        auxiliary.h                                            *
* Source:      GNU C++                                                *
* Purpose:     prototypes for auxiliary functions in expression.cxx   *
* History:     030121 0.0 work started                                *
* License:     (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
***********************************************************************/

bool is_integer(double);
bool is_even(double);
bool is_odd(double);
double Ev3NearZero(void);
double Ev3Infinity(void);
double argmin(double, double, double, double);
double argmax(double, double, double, double);
void bilinearprodmkrange(double, double, double, double, double*, double*);
void fractionmkrange(double, double, double, double, double*, double*);
void constpowermkrange(double, double, double, double*, double*);
void powermkrange(double, double, double, double, double*, double*);

