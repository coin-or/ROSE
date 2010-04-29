/**********************************************************************
* Author:      Leo Liberti                                            *
* Name:        common.h                                               *
* Source:      GNU C++                                                *
* Purpose:     common stuff                                           *
* History:     050909 0.0 work started                                *
* License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
***********************************************************************/

#ifndef __EV3COMMONH__
#define __EV3COMMONH__

#define RCS12 "$Id: common.h,v 1.2 2006/07/30 05:36:44 liberti Exp liberti $"

#define NOVARIABLE -1
#define OPERANDSTRBUF 64
#define EXPRESSIONSTRBUF 2048
#define LARGE 1E10

#include <cassert>

typedef int Int;

// various operator types
enum OperatorType {  
  SUM, DIFFERENCE, PRODUCT, FRACTION, POWER,
  PLUS, MINUS, LOG, EXP, SIN, 
  COS, TAN, COT, SINH, COSH, 
  TANH, COTH, SQRT, VAR, CONST,
  ERROR
};

// utility functions
extern double Ev3NearZero(void);
extern double Ev3Infinity(void);

#endif
