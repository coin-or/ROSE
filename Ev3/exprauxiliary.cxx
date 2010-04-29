/**********************************************************************
* Author:      Leo Liberti                                            *
* Name:        exprauxiliary.cxx                                      *
* Source:      GNU C++                                                *
* Purpose:     symbolic expression (auxiliary functions)              *
* History:     080928 0.0 work started                                *
* License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
***********************************************************************/

#include <iostream>
#include <cmath>
#include <cstdlib>
#include "common.h"
#include "exprauxiliary.h"

////////////// auxiliary functions ///////////////
	   
#define ISINTTOLERANCE 1e-8
bool is_integer(double a) {
  double b = fabs(a);
  int bi = (int) rint(b);
  if (fabs(b - bi) < ISINTTOLERANCE) {
    return true;
  } else {
    return false;
  }
}

bool is_even(double a) {
  if (is_integer(a)) {
    int ai = (int) rint(a);
    if (ai % 2 == 0) {
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}

bool is_odd(double a) {
  if (is_integer(a)) {
    int ai = (int) rint(a);
    if (ai % 2 == 1) {
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}

double Ev3NearZero(void) {
  // returns a very small positive value
  return 1/LARGE;
}

double Ev3Infinity(void) {
  // returns a very large positive value
  return LARGE;
}

double argmin(double a1, double a2, double a3, double a4) {
  double ret;
  ret = a1;
  if (a2 < ret)
    ret = a2;
  if (a3 < ret)
    ret = a3;
  if (a4 < ret)
    ret = a4;
  return ret;
}

double argmax(double a1, double a2, double a3, double a4) {
  double ret;
  ret = a1;
  if (a2 > ret)
    ret = a2;
  if (a3 > ret)
    ret = a3;
  if (a4 > ret)
    ret = a4;
  return ret;
}

void bilinearprodmkrange(double a, double b, double c, double d, 
			 double* wl, double *wu) {
  if ((a >= 0 && b >= 0 && c >= 0 && d >= 0) || \
      (a <  0 && b <  0 && c <  0 && d <  0)) {
    if(a*c <= b*d) {
      *wl = a * c;
      *wu = b * d;
    } else {
      *wl = b * d;
      *wu = a * c;
    }
  } else if (a <  0 && b <  0 && c >= 0 && d >= 0) {
    *wl = a * d;
    *wu = b * c;
  } else if (a >= 0 && b >= 0 && c <  0 && d <  0) {
    *wl = b * c;
    *wu = a * d;
  } else if (a <  0 && b >= 0 && c >= 0 && d >= 0) {
    *wl = a * d;
    *wu = b * d;
  } else if (a <  0 && b <  0 && c <  0 && d >= 0) {
    *wl = a * d;
    *wu = a * c;
  } else if (a >= 0 && b >= 0 && c <  0 && d >= 0) {
    *wl = b * c;
    *wu = b * d;
  } else if (a <  0 && b >= 0 && c <  0 && d <  0) {
    *wl = b * c;
    *wu = a * c;
  } else if (a <  0 && b >= 0 && c <  0 && d >= 0) {
    *wl = std::min(a * d, b * c);
    *wu = std::max(a * c, b * d);
  }
}
	
void fractionmkrange(double a, double b, double c, double d,
		     double* wl, double* wu) {
  if ((a >= 0 && b >= 0 && c >  0 && d >  0) || 
      (a <  0 && b <  0 && c <  0 && d <  0)) { 
    *wl = a / d;
    *wu = b / c;
  } else if ((a <  0 && b <  0 && c >  0 && d >  0) || 
	     (a >= 0 && b >= 0 && c <  0 && d <  0)) {
    *wl = b / c;
    *wu = a / d;
  } else if (a <  0 && b >= 0 && c <  0 && d <  0) {
    *wl = b / d;
    *wu = a / d;
  } else if (a <  0 && b >= 0 && c >  0 && d >  0) {
    *wl = a / c;
    *wu = b / c;
  } else if (c <  0 && d >  0) {
    *wl = -Ev3Infinity();
    *wu = Ev3Infinity();
  }
  if (c == 0 && d > 0) {
    if (a >= 0 && b >= 0) {
      *wl = a / d;
      *wu = Ev3Infinity();
    } else if (a < 0 && b >= 0) {
      *wl = -Ev3Infinity();
      *wu = Ev3Infinity();
    } else if (a < 0 && b < 0) {
      *wl = -Ev3Infinity();
      *wu = b / d;
    }
  } else if (c == 0 && d == 0) {
    *wl = -Ev3Infinity();
    *wu = Ev3Infinity();
  } else if (c < 0 && d == 0) {
    if (a >= 0 && b >= 0) {
      *wl = b / c;
      *wu = Ev3Infinity();
    } else if (a < 0 && b >= 0) {
      *wl = -Ev3Infinity();
      *wu = Ev3Infinity();
    } else if (a < 0 && b < 0) {
      *wl = -Ev3Infinity();
      *wu = a / c;
    }
  }
}   
	
void constpowermkrange(double a, double b, double c,
		       double* wl, double* wu) {
  double tmp = 0;
  double tmp1 = 0;
  if (is_integer(c)) {
    // exponent is integer
    tmp = std::pow(a, c);
    tmp1 = std::pow(b, c);
    if (is_even(c)) {
      if (a >= 0 && b >= 0) {
	*wl = tmp;
	*wu = tmp1;
      } else if (a < 0 && b < 0) {
	*wl = tmp1;
	*wu = tmp;
      } else {
	*wl = 0;
	*wu = std::max(tmp, tmp1);
      }
    } else {
      *wl = tmp;
      *wu = tmp1;
    }
  } else {
    // exponent not an integer: range can't go into negative numbers
    if (a < 0) {
      a = 0;
    }
    if (b < 0) {
      b = 0;
    }
    // calculate range
    tmp = std::pow(a, c);
    tmp1 = std::pow(b, c);
    if (tmp > tmp1) {
      double t = tmp;
      tmp = tmp1;
      tmp1 = t;
    }
    *wl = tmp;
    *wu = tmp1;
  }
  if (*wl <= -Ev3Infinity()) 
    *wl = -Ev3Infinity();
  if (*wu >= Ev3Infinity())
    *wu = Ev3Infinity();
}

void powermkrange(double a, double b, double c, double d,
		  double* wl, double* wu) {
  if (fabs(d-c) < 1/LARGE) {
    constpowermkrange(a, b, c, wl, wu);
    return;
  }
  // WARNING: don't have time to do it properly now
  double t1 = 0;
  double t2 = 0;
  double t3 = 0;
  double t4 = 0;
  t1 = std::pow(a, c);
  t2 = std::pow(a, d);
  t3 = std::pow(b, c);
  t4 = std::pow(b, d);
  *wl = argmin(t1, t2, t3, t4);
  *wu = argmax(t1, t2, t3, t4);

  std::cerr << "expression.cxx::powermkrange(): method not functional\n";
  exit(148);
}

int computecolor(int symm, int level, int op, int cnst, int var, int ord, 
		 int all_lev_ops, int all_ops, int all_cnsts, 
		 int all_vars, int all_ords) {
  int color = symm * (all_lev_ops * all_cnsts * all_vars * all_ords) +
    level * (all_ops * all_cnsts * all_vars * all_ords) +
    op * (all_cnsts * all_vars * all_ords) +
    cnst * (all_vars * all_ords) +
    var * all_ords + 
    ord;
  return color;
}
