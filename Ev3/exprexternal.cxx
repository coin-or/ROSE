/**********************************************************************
* Author:      Leo Liberti                                            *
* Name:        exprexternal.cxx                                       *
* Source:      GNU C++                                                *
* Purpose:     symbolic expression (external functions)               *
* History:     080928 0.0 work started                                *
* License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
***********************************************************************/

#include <string>
#if __GNUC__ >= 3
#include <sstream>
#else
#include <strstream>
#endif
#include <cmath>
#include <algorithm>
#include <cassert>
#include "expression.h"

using namespace std;

#define PEV3PI 3.1415926535897932385
#define NOTVARNAME "_var_not_found_"
#define VNAMEIDXCHAR "_"
//#define VNAMEIDXCHAR ""

/************** expression creation (no change to args) ***************/

// BIG FAT WARNING: when you change these operators, also please
// change their "-Link" counterparts!

// sums:
Expression operator + (Expression a, Expression b) {
  Expression ret;
  // make a preliminary check
  if (a->GetCoeff() == 0 || a->HasValue(0)) {
    ret.SetToCopyOf(b);
    return ret;
  }
  if (b->GetCoeff() == 0 || b->HasValue(0)) {
    ret.SetToCopyOf(a);
    return ret;  
  }
  if (!(a->IsConstant() && b->IsConstant()) && a->IsEqualToNoCoeff(b)) {
    a->SetCoeff(a->GetCoeff() + b->GetCoeff());
    if (fabs(a->GetCoeff()) < Ev3NearZero()) {
      // simplify to zero - for differences
      Expression zero(0.0);
      return zero;
    } else {
      ret.SetToCopyOf(a);
      return ret;
    }
  }
  // go for it
  if (a->IsLeaf() && a->GetOpType() == CONST &&
      b->IsLeaf() && b->GetOpType() == CONST) {
    // a, b are numbers - add them
    ret.SetToCopyOf(a);
    ret->SetValue(a->GetValue() + b->GetValue());
    ret->SetCoeff(1);
    ret->SetExponent(1);
    return ret;
  } else if (a->IsLeaf() && a->GetOpType() == VAR &&
	     b->IsLeaf() && b->GetOpType() == VAR &&
	     a->GetVarIndex() == b->GetVarIndex() &&
	     a->GetExponent() == b->GetExponent()) {
    // a, b are the same variable - add coefficients
    ret.SetToCopyOf(a);    
    ret->SetCoeff(a->GetCoeff() + b->GetCoeff());
    return ret;
  } else if (a->GetOpType() == SUM && b->GetOpType() != SUM) {    
    // a is a sum and b isn't - just add the term b
    ret.SetToCopyOf(a);    
    ret->DistributeCoeffOverSum();
    Int i = 0;
    bool couldsimplify = false;
    Expression tmp;
    // if t is a leaf and there is a like leaf in this,
    // just add it to the value/coefficient
    if (b->IsLeaf() && b->GetOpType() == CONST) {
      // b is a constant
      for (i = 0; i < ret->GetSize(); i++) {
	tmp = ret->GetNode(i);
	if (tmp->IsLeaf() && tmp->GetOpType() == CONST) {
	  tmp->SetValue(tmp->GetValue() + b->GetValue() / ret->GetCoeff());
	  tmp->SetCoeff(1);
	  tmp->SetExponent(1); // NB: changing tmp should also change a
	  couldsimplify = true;
	  break;
	}
      }
    } else if (b->IsLeaf() && b->GetOpType() == VAR) {
      // b is a variable
      for (i = 0; i < ret->GetSize(); i++) {
	if (ret->GetNode(i)->IsLeaf() && ret->GetNode(i)->GetOpType() == VAR &&
	    b->GetVarIndex() == ret->GetNode(i)->GetVarIndex() &&
	    b->GetExponent() == ret->GetNode(i)->GetExponent()) {
	  double tc = ret->GetNode(i)->GetCoeff() + 
	    b->GetCoeff() / ret->GetCoeff();
	  // warning: tc could be zero, but it would be cumbersome
	  // to simplify it here - do it in SimplifyConstant
	  ret->GetNode(i)->SetCoeff(tc);
	  couldsimplify = true;
	  break;
	}
      }
    } else if (!b->IsLeaf()) {
      // a is a sum, b is a nonleaf, look for a subnode of a similar to b
      for (i = 0; i < ret->GetSize(); i++) {
	if (ret->GetNode(i)->IsEqualTo(b)) {
	  // found one, add coefficients - notice, as above, coeff could
	  // be zero, but deal with that case in SimplifyConstant
	  ret->GetNode(i)->SetCoeff(ret->GetNode(i)->GetCoeff() 
				    + b->GetCoeff());
	  couldsimplify = true;
	  break;
	}
      }
    }
    if (!couldsimplify) {
      // either could not simplify in steps above, or b is an operator
      ret->AddCopyOfNode(b);
    }
    return ret;
  } else if (a->GetOpType() == SUM && b->GetOpType() == SUM) {    
    // a, b are sums - add terms of b to a
    b->DistributeCoeffOverSum();
    ret.SetToCopyOf(a);
    Int i = 0;
    Int s = b->GetSize();
    for (i = 0; i < s; i++) {
      ret = ret + b->GetNode(i);
    }
    return ret;
  } else if (a->GetOpType() != SUM && b->GetOpType() == SUM) {
    // a is not a sum but b is - transform this into a sum
    ret.SetToCopyOf(b);
    ret = ret + a;
    return ret;
  } else {
    // all other cases - make new node on top of the addends
    ret->SetOpType(SUM);
    ret->SetCoeff(1);
    ret->SetExponent(1);
    ret->AddCopyOfNode(a);
    ret->AddCopyOfNode(b);
    return ret;
  }
}


// product
Expression operator * (Expression a, Expression t) {
  Expression ret;
  // make a preliminary check  
  if (a->GetCoeff() == 0 || t->GetCoeff() == 0 || 
      a->HasValue(0) || t->HasValue(0)) {
    Expression zero(0.0);
    return zero;
  }
  if (a->HasValue(1)) {
    ret.SetToCopyOf(t);
    return ret;
  }
  if (t->HasValue(1)) {
    ret.SetToCopyOf(a);
    return ret; 
  }
  if (!(a->IsConstant() && t->IsConstant()) && a->IsEqualToNoCoeff(t)) {
    Expression two(2.0);
    ret.SetToCopyOf(a);
    ret->SetCoeff(1.0);
    ret = ret^two;
    ret->SetCoeff(a->GetCoeff() * t->GetCoeff());
    return ret;
  }
  // go for it
  if (a->IsLeaf() && a->GetOpType() == CONST &&
      t->IsLeaf() && t->GetOpType() == CONST) {
    // a, t are numbers - multiply them
    ret.SetToCopyOf(a);    
    ret->SetValue(a->GetValue() * t->GetValue());
    ret->SetCoeff(1);
    ret->SetExponent(1);
    return ret;
  } else if (a->IsLeaf() && a->GetOpType() == VAR &&
	     t->IsLeaf() && t->GetOpType() == VAR &&
	     a->GetVarIndex() == t->GetVarIndex()) {
    // a, t are the same variable - multiply coefficients
    // and add exponents
    ret.SetToCopyOf(a);
    ret->SetCoeff(a->GetCoeff() * t->GetCoeff());
    ret->SetExponent(a->GetExponent() + t->GetExponent());
    return ret;
  } else if (t->IsConstant()) {
    // t is constant, set coeff of a
    ret.SetToCopyOf(a);
    ret->SetCoeff(a->GetCoeff() * t->GetValue());
    ret->DistributeCoeffOverSum();
    return ret;
  } else if (a->IsConstant()) {
    // a is constant, set coeff of t
    ret.SetToCopyOf(t);
    ret->SetCoeff(t->GetCoeff() * a->GetValue());
    ret->DistributeCoeffOverSum();
    return ret;
  } else if (a->GetOpType() == PRODUCT && t->GetOpType() != PRODUCT) {
    // a is a product and t isn't - just multiply the term t
    ret.SetToCopyOf(a);
    Int i = 0;
    bool couldsimplify = false;
    if (t->IsLeaf() && t->GetOpType() == VAR) {
      // t is a variable
      Expression tmp;
      for (i = 0; i < ret->GetSize(); i++) {
	tmp = ret->GetNode(i);
	if (tmp->IsLeaf() && tmp->GetOpType() == VAR &&
	    t->GetVarIndex() == tmp->GetVarIndex()) {
	  // found same variable in a, multiply coeffs and add exponents
	  tmp->SetCoeff(tmp->GetCoeff() * t->GetCoeff());	
	  tmp->SetExponent(tmp->GetExponent() + t->GetExponent());
	  couldsimplify = true;
	  break;
	}
      }
    } 
    // here we shan't try to simplify f*f <-- f^2 (f nonleaf) 
    // because a product of nonleaves is easier to manipulate than 
    // a power (as it adds a depth level)
    if (!couldsimplify) {
      // either could not simplify in steps above, or t is an operator
      ret->AddCopyOfNode(t);
    }
    return ret;
  } else if (a->GetOpType() == PRODUCT && t->GetOpType() == PRODUCT) {
    // a, t are products - multiply terms of t to a
    t->DistributeCoeffOverProduct();
    ret.SetToCopyOf(a);
    Int i = 0;
    Int s = t->GetSize();
    for (i = 0; i < s; i++) {
      ret = ret * t->GetNode(i);
    }
    return ret;
  } else if (a->GetOpType() != PRODUCT && t->GetOpType() == PRODUCT) {
    // a is not a products but t is - transform this into a product
    ret.SetToCopyOf(t);
    ret = ret * a;
    return ret;
  } else {
    // all other cases - make new node on top of the addends
    ret->SetCoeff(1);
    ret->SetExponent(1);
    ret->SetOpType(PRODUCT);
    ret->AddCopyOfNode(a);
    ret->AddCopyOfNode(t);
    return ret;
  }
}

// fractions:
Expression operator / (Expression a, Expression t) throw(ErrDivideByZero) {
  Expression ret;
  // make a preliminary check
#ifdef DEBUG2
  cout << "----diff----" << endl;
  cout << a->PrintTree(2, 2) << endl;
  cout << "  --over--" << endl;
  cout << t->PrintTree(2, 2) << endl;
  cout << "------------" << endl;
#endif
  if (t->GetCoeff() == 0) {
    // divide by zero
    unsigned long mycode(0);
    string myif("Expression Building");
    string myscope("operator/");
    string myop("t.GetCoeff()==0");
    string mydesc("Divisor cannot be zero");
    string myinfo(HELPURL);
    string mydiv(NONE);
    throw ErrDivideByZero(mycode,myif,myscope,myop,mydesc,myinfo,mydiv);
  }
  if (a->GetCoeff() == 0 || a->HasValue(0)) {
    // dividend has zero coeff, return zero
    Expression zero(0.0);
    return zero;
  }
  if (t->HasValue(1)) {
    ret.SetToCopyOf(a);    
    return ret;
  }
  if (!(a->IsConstant() && t->IsConstant()) && a->IsEqualToNoCoeff(t)) {
    // dividend = divisor, return ratio of coefficients
    Expression one(1.0);
    one->SetCoeff(a->GetCoeff() / t->GetCoeff());
    return one;
  }
  // go for it
  if (a->IsLeaf() && a->GetOpType() == CONST &&
      t->IsLeaf() && t->GetOpType() == CONST) {
    // a, t are numbers - divide them
    if (t->GetValue() == 0) {
      unsigned long mycode(0);
      string myif("Expression Building");
      string myscope("operator/");
      string myop("t.GetValue()==0");
      string mydesc("Divisor cannot be zero");
      string myinfo(HELPURL);
      string mydiv(NONE);
      throw ErrDivideByZero(mycode,myif,myscope,myop,mydesc,myinfo,mydiv);
    } else {
      ret.SetToCopyOf(a);
      ret->SetValue(a->GetValue() / t->GetValue());
      ret->SetCoeff(1);
      ret->SetExponent(1);
      return ret;
    }
  } else if (t->HasValue(1)) {
    // divide by constant 1, don't do anything
    ret.SetToCopyOf(a);
    return ret;
  } else if (t->IsConstant()) {
    // t is constant, set coeff of a
    ret.SetToCopyOf(a);
    ret->SetCoeff(a->GetCoeff() / t->GetValue());
    ret->DistributeCoeffOverSum();
    return ret;
  } else if (a->IsVariable() && t->IsVariable() &&
	     a->GetVarIndex() == t->GetVarIndex()) {
    // cx^e / dx^f = (c/d)x^(e-f)
    ret.SetToCopyOf(a);
    double te = a->GetExponent() - t->GetExponent();
    double tc = a->GetCoeff() / t->GetCoeff();
    if (fabs(te) < Ev3NearZero()) {
      Expression c(tc);
      return tc;
    }
    ret->SetCoeff(tc);
    ret->SetExponent(te);
    return ret;
  } else if (a->IsVariable() && t->GetOpType() == PRODUCT) {
    // a is a variable, t is a product - see if a appears in t
    // and cancel common term
    // first simplify coeffs of divisor
    Expression at;
    at.SetToCopyOf(a);
    ret.SetToCopyOf(t);
    ret->ConsolidateProductCoeffs();
    // denominator
    if (fabs(ret->GetCoeff()) < Ev3NearZero()) {
      // divide by zero
      unsigned long mycode(22);
      string myif("Expression Building");
      string myscope("operator/");
      string myop("t->GetCoeff()");
      string mydesc("Divisor cannot be zero");
      string myinfo(HELPURL);
      string mydiv(NONE);
      throw ErrDivideByZero(mycode,myif,myscope,myop,mydesc,myinfo,mydiv);
    }
    if (fabs(at->GetCoeff()) < Ev3NearZero()) {
      Expression zero(0.0);
      return zero;
    }
    double accumulatedcoeff = at->GetCoeff() / ret->GetCoeff();
    at->SetCoeff(1.0);
    ret->SetCoeff(1.0);
    // now try simplification
    Int i;
    for (i = 0; i < ret->GetSize(); i++) {
      if (ret->GetNode(i)->GetOpType() == VAR && 
	  at->GetVarIndex() == ret->GetNode(i)->GetVarIndex()) {
	double te = at->GetExponent() - ret->GetNode(i)->GetExponent();
	if (fabs(te) < Ev3NearZero()) {
	  // exponents are the same, just cancel
	  at->One();
	  ret->DeleteNode(i);
	} else if (te > 0) {
	  // numerator remains, cancel denominator
	  at->SetExponent(te);
	  ret->DeleteNode(i);
	} else if (te < 0) {
	  // numerator goes to one, denominator remains
	  at->One();
	  ret->GetNode(i)->SetExponent(-te);
	}
	// exit loop
	break;
      }
    }
    // check that denominator (t) has more than one operand;
    // if not, bring up a rank level
    if (ret->GetSize() == 1) {
      ret = ret->GetNode(0);
    }
    // build ratio
    Expression ret2;
    ret2->SetOpType(FRACTION);
    ret2->SetCoeff(accumulatedcoeff);
    ret2->SetExponent(1);
    ret2->AddCopyOfNode(at);
    ret2->AddCopyOfNode(ret);  
    return ret2;
  } else if (t->IsVariable() && a->GetOpType() == PRODUCT) {
    // t is a variable, a is a product - see if t appears in a
    // and cancel common term
    // first simplify coeffs of divisor
    Expression bt;
    bt.SetToCopyOf(t);
    ret.SetToCopyOf(a);    
    ret->ConsolidateProductCoeffs();
    // denominator - already checked
    if (fabs(ret->GetCoeff()) < Ev3NearZero()) {
      Expression zero(0.0);
      return zero;
    }
    double accumulatedcoeff = ret->GetCoeff() / bt->GetCoeff();
    ret->SetCoeff(1.0);
    bt->SetCoeff(1.0);
    // now try simplification
    Int i;
    for (i = 0; i < ret->GetSize(); i++) {
      if (ret->GetNode(i)->GetOpType() == VAR && 
	  bt->GetVarIndex() == ret->GetNode(i)->GetVarIndex()) {
	double te = ret->GetNode(i)->GetExponent() - bt->GetExponent();
	if (fabs(te) < Ev3NearZero()) {
	  // exponents are the same, just cancel
	  bt->One();
	  ret->DeleteNode(i);
	} else if (te > 0) {
	  // numerator remains, cancel denominator
	  bt->One();
	  ret->GetNode(i)->SetExponent(te);
	} else if (te < 0) {
	  // numerator goes to one, denominator remains
	  bt->SetExponent(-te);
	  ret->DeleteNode(i);
	}
	// exit loop
	break;
      }
    }
    // check that numerator (a) has more than one operands;
    // if not, bring up a rank level
    if (ret->GetSize() == 1) {
      ret = ret->GetNode(0);
    }
    // build ratio
    Expression ret2;
    ret2->SetOpType(FRACTION);
    ret2->SetCoeff(accumulatedcoeff);
    ret2->SetExponent(1);
    ret2->AddCopyOfNode(ret);
    ret2->AddCopyOfNode(bt);  
    return ret2;
  } else if (a->GetOpType() == PRODUCT && t->GetOpType() == PRODUCT) {
    // a, t are products, try to cancel common terms
    Expression at;
    Expression bt;
    at.SetToCopyOf(a);
    bt.SetToCopyOf(t);
    Int i = 0, j = 0;
    double accumulatedcoeff;
    // first simplify coefficients of operands
    at->ConsolidateProductCoeffs();
    bt->ConsolidateProductCoeffs();
    // denominator
    if (fabs(bt->GetCoeff()) < Ev3NearZero()) {
      // divide by zero
      unsigned long mycode(21);
      string myif("Expression Building");
      string myscope("operator/");
      string myop("t->GetCoeff()");
      string mydesc("Divisor cannot be zero");
      string myinfo(HELPURL);
      string mydiv(NONE);
      throw ErrDivideByZero(mycode,myif,myscope,myop,mydesc,myinfo,mydiv);
    }
    if (fabs(at->GetCoeff()) < Ev3NearZero()) {
      Expression zero(0.0);
      return zero;
    }
    // save ratio of coeffs of products
    accumulatedcoeff = at->GetCoeff() / bt->GetCoeff();
    at->SetCoeff(1.0);
    bt->SetCoeff(1.0);
    // now try simplification
    i = 0;
    bool isnumeratorempty = false;
    bool isdenominatorempty = false;
    int szi = at->GetSize();
    int szj = bt->GetSize();
    while(!isnumeratorempty && !isdenominatorempty && i < szi) {
      j = 0;
      while(!isnumeratorempty && !isdenominatorempty && j < szj) {
	if (at->GetNode(i)->IsEqualTo(bt->GetNode(j))) {
	  // found like terms i and j
	  at->DeleteNode(i);
	  szi--;
	  if(szi == 0) {
	    isnumeratorempty = true;
	    at->One();
	  }
	  bt->DeleteNode(j);
	  szj--;
	  if (szj == 0) {
	    isdenominatorempty = true;
	    bt->One();
	  }
	  i--;   // cancel the effect of the later i++
	  break; // go to outer cycle
	} else {
	  j++;
	}
      }	
      i++;
    }
    if (bt->HasValue(1)) {
      // denominator is 1, return a
      at->SetCoeff(accumulatedcoeff);
      return at;
    }
    // now construct fraction
    // check that numerator, denominator have more than one operands;
    // if not, bring up a rank level
    if (at->GetSize() == 1) {
      at = at->GetNode(0);
    }
    if (bt->GetSize() == 1) {
      bt = bt->GetNode(0);
    }
    ret->SetCoeff(accumulatedcoeff); // already contains coeffs of a, t
    ret->SetExponent(1);
    ret->SetOpType(FRACTION);
    ret->AddCopyOfNode(at);
    ret->AddCopyOfNode(bt);
    return ret;
  } else {
    Expression at;
    Expression bt;
    at.SetToCopyOf(a);
    bt.SetToCopyOf(t);
    ret->SetCoeff(at->GetCoeff() / bt->GetCoeff());
    at->SetCoeff(1);
    bt->SetCoeff(1);
    ret->SetExponent(1);
    ret->SetOpType(FRACTION);
    ret->AddCopyOfNode(at);
    ret->AddCopyOfNode(bt);
    return ret;
  }
}    

// unary minus:
Expression operator - (Expression a) {
  Expression ret;
  ret.SetToCopyOf(a);
  if (ret->IsLeaf() && ret->GetOpType() == CONST) {
    ret->SetValue(- ret->GetValue());
    ret->SetCoeff(1);
    ret->SetExponent(1);    
  } else {
    ret->SetCoeff(- ret->GetCoeff());
  }
  return ret;
}

// binary minus:
Expression operator - (Expression a, Expression b) {
  Expression ret;
  if (a->HasValue(0)) 
    return -b;
  if (b->HasValue(0)) {
    ret.SetToCopyOf(a);
    return a;
  }
  ret = a + (-b);
  return ret;
}

// power:
Expression operator ^ (Expression a, Expression t) {
  // make a preliminary check
  Expression ret;
  if (a->GetCoeff() == 0) {
    // *this is zero, just return zero
    Expression zero(0.0);
    return zero;
  }
  if (t->HasValue(0.0)) {
    // exponent is 0, just return 1
    Expression one(1.0);
    return one;
  } else if (t->HasValue(1.0)) {
    // exponent is 1, just return a
    ret.SetToCopyOf(a);
    return ret;
  } 
  if (a->HasValue(0.0)) {
    // base is zero, return 0
    Expression zero(0.0);
    return zero;
  } else if (a->HasValue(1.0)) {
    // base is one, return 1
    Expression one(1.0);
    return one;
  }
  // go for it
  if (a->IsLeaf() && a->GetOpType() == CONST &&
      t->IsLeaf() && t->GetOpType() == CONST) {
    // constant to constant
    ret.SetToCopyOf(a);
    ret->SetValue(pow(ret->GetValue(), t->GetValue()));
    ret->SetCoeff(1);
    ret->SetExponent(1);
    return ret;
  } else if (a->IsLeaf() && a->GetOpType() == VAR &&
	     t->IsLeaf() && t->GetOpType() == CONST) {
    // variable to constant
    ret.SetToCopyOf(a);
    ret->SetCoeff(pow(ret->GetCoeff(), t->GetValue()));
    ret->SetExponent(ret->GetExponent() * t->GetValue());
    return ret;
  } else {
    ret->SetCoeff(1);
    ret->SetExponent(1);
    ret->SetOpType(POWER);
    ret->AddCopyOfNode(a);
    ret->AddCopyOfNode(t);
    return ret;
  }
}

Expression Log(Expression a) throw(ErrNotPermitted) {
  // make a preliminary check
  if (a->IsZero()) {
    // *this is zero, can't do
    unsigned long mycode(0);
    string myif("Expression Building");
    string myscope("Log");
    string myop("IsZero()");
    string mydesc("log(0) is undefined");
    string myinfo(HELPURL);
    throw ErrNotPermitted(mycode,myif,myscope,myop,mydesc,myinfo);
  }
  if (a->IsLessThan(0)) {
    // argument is < zero, can't do
    unsigned long mycode(0);
    string myif("Expression Building");
    string myscope("Log");
    string myop("value <= 0");
    string mydesc("log(<=0) is undefined");
    string myinfo(HELPURL);
    throw ErrNotPermitted(mycode,myif,myscope,myop,mydesc,myinfo);
  } 
  // go for it
  if (a->IsLeaf() && a->GetOpType() == CONST) {
    Expression ret;
    ret.SetToCopyOf(a);
    double t = ret->GetValue();
    assert(t >= 0);
    ret->SetCoeff(1);    
    ret->SetValue(log(t));
    ret->SetExponent(1);
    ret->SetOpType(CONST);
    return ret;
  } else {
    Expression ret;
    ret->SetCoeff(1);
    ret->SetExponent(1);
    ret->SetOpType(LOG);
    ret->AddCopyOfNode(a);
    return ret;
  }
}

Expression Exp(Expression a) {
  // go for it
  if (a->IsLeaf() && a->GetOpType() == CONST) {
    Expression ret;
    ret.SetToCopyOf(a);
    ret->SetCoeff(1);
    ret->SetValue(exp(a->GetValue()));
    ret->SetExponent(1);
    ret->SetOpType(CONST);
    return ret;
  } else {
    Expression ret;
    ret->SetCoeff(1);
    ret->SetExponent(1);
    ret->SetOpType(EXP);
    ret->AddCopyOfNode(a);
    return ret;
  }
}

Expression Sin(Expression a) {
  // go for it
  if (a->IsLeaf() && a->GetOpType() == CONST) {
    Expression ret;
    ret.SetToCopyOf(a);
    ret->SetCoeff(1);
    ret->SetValue(sin(a->GetValue()));
    ret->SetExponent(1);
    ret->SetOpType(CONST);
    return ret;
  } else {
    Expression ret;
    ret->SetCoeff(1);
    ret->SetExponent(1);
    ret->SetOpType(SIN);
    ret->AddCopyOfNode(a);
    return ret;
  }
}

Expression Cos(Expression a) {
  // go for it
  if (a->IsLeaf() && a->GetOpType() == CONST) {
    Expression ret;
    ret.SetToCopyOf(a);
    ret->SetCoeff(1);
    ret->SetValue(cos(a->GetValue()));
    ret->SetExponent(1);
    ret->SetOpType(CONST);
    return ret;
  } else {
    Expression ret;
    ret->SetCoeff(1);
    ret->SetExponent(1);
    ret->SetOpType(COS);
    ret->AddCopyOfNode(a);
    return ret;
  }
}

Expression Tan(Expression a) {
  // go for it
  if (a->IsLeaf() && a->GetOpType() == CONST) {
    Expression ret;
    ret.SetToCopyOf(a);
    ret->SetCoeff(1);
    ret->SetValue(tan(a->GetValue()));
    ret->SetExponent(1);
    ret->SetOpType(CONST);
    return ret;
  } else {
    Expression ret;
    ret->SetCoeff(1);
    ret->SetExponent(1);
    ret->SetOpType(TAN);
    ret->AddCopyOfNode(a);
    return ret;
  }
}

Expression Cot(Expression a)  throw(ErrNotPermitted) {
  // make a preliminary check
  if (a->IsZero()) {
    // *this is zero, can't do
    unsigned long mycode(0);
    string myif("Expression Building");
    string myscope("Cot");
    string myop("IsZero()");
    string mydesc("cot(0) is undefined");
    string myinfo(HELPURL);
    throw ErrNotPermitted(mycode,myif,myscope,myop,mydesc,myinfo);
  }
  // go for it
  if (a->IsLeaf() && a->GetOpType() == CONST) {
    Expression ret;
    ret.SetToCopyOf(a);
    double t = tan(a->GetValue());
    assert(t != 0);
    ret->SetCoeff(1);
    ret->SetValue(1 / t);
    ret->SetExponent(1);
    ret->SetOpType(CONST);
    return ret;
  } else {
    Expression ret;
    ret->SetCoeff(1);
    ret->SetExponent(1);
    ret->SetOpType(COT);
    ret->AddCopyOfNode(a);
    return ret;
  }
}

Expression Sinh(Expression a) {
  // go for it
  if (a->IsLeaf() && a->GetOpType() == CONST) {
    Expression ret;
    ret.SetToCopyOf(a);
    ret->SetCoeff(1);
    ret->SetValue(sinh(a->GetValue()));
    ret->SetExponent(1);
    ret->SetOpType(CONST);
    return ret;
  } else {
    Expression ret;
    ret->SetCoeff(1);
    ret->SetExponent(1);
    ret->SetOpType(SINH);
    ret->AddCopyOfNode(a);
    return ret;
  }
}

Expression Cosh(Expression a) {
  // go for it
  if (a->IsLeaf() && a->GetOpType() == CONST) {
    Expression ret;
    ret.SetToCopyOf(a);
    ret->SetCoeff(1);
    ret->SetValue(cosh(a->GetValue()));
    ret->SetExponent(1);
    ret->SetOpType(CONST);
    return ret;
  } else {
    Expression ret;
    ret->SetCoeff(1);
    ret->SetExponent(1);
    ret->SetOpType(COSH);
    ret->AddCopyOfNode(a);
    return ret;
  }
}

Expression Tanh(Expression a) {
  // go for it
  if (a->IsLeaf() && a->GetOpType() == CONST) {
    Expression ret;
    ret.SetToCopyOf(a);
    ret->SetCoeff(1);
    ret->SetValue(tanh(a->GetValue()));
    ret->SetExponent(1);
    ret->SetOpType(CONST);
    return ret;
  } else {
    Expression ret;
    ret->SetCoeff(1);
    ret->SetExponent(1);
    ret->SetOpType(TANH);
    ret->AddCopyOfNode(a);
    return ret;
  }
}

Expression Coth(Expression a)  throw(ErrNotPermitted) {
  // make a preliminary check
  if (a->IsZero()) {
    // *this is zero, can't do
    unsigned long mycode(0);
    string myif("Expression Building");
    string myscope("Coth");
    string myop("IsZero()");
    string mydesc("coth(0) is undefined");
    string myinfo(HELPURL);
    throw ErrNotPermitted(mycode,myif,myscope,myop,mydesc,myinfo);
  }
  // go for it
  if (a->IsLeaf() && a->GetOpType() == CONST) {
    Expression ret;
    ret.SetToCopyOf(a);
    double t = tanh(a->GetValue());
    assert(t != 0);
    ret->SetCoeff(1);
    ret->SetValue(1 / t);
    ret->SetExponent(1);
    ret->SetOpType(CONST);
    return ret;
  } else {
    Expression ret;
    ret->SetCoeff(1);
    ret->SetExponent(1);
    ret->SetOpType(COTH);
    ret->AddCopyOfNode(a);
    return ret;
  }
}

Expression Sqrt(Expression a) throw(ErrNotPermitted) {
  // make a preliminary check
  if (a->IsLessThan(0) && !a->HasValue(0)) {
    // argument is < zero, can't do
    unsigned long mycode(0);
    string myif("Expression Building");
    string myscope("Sqrt");
    string myop("value < 0");
    string mydesc("sqrt(<0) is complex, can't do");
    string myinfo(HELPURL);
    throw ErrNotPermitted(mycode,myif,myscope,myop,mydesc,myinfo);
  } 
  // go for it
  if (a->IsLeaf() && a->GetOpType() == CONST) {
    Expression ret;
    ret.SetToCopyOf(a);
    double t = a->GetValue();
    assert(t >= 0);
    ret->SetCoeff(1);
    ret->SetValue(sqrt(t));
    ret->SetExponent(1);
    ret->SetOpType(CONST);
    return ret;
  } else {
    Expression ret;
    ret->SetCoeff(1);
    ret->SetExponent(1);
    ret->SetOpType(SQRT);
    ret->AddCopyOfNode(a);
    return ret;
  }
}

/***************** expression creation (affects arguments) ***********/

// sums:
Expression SumLink(Expression a, Expression b) {
  // make a preliminary check
  if (a->GetCoeff() == 0 || a->HasValue(0))
    return b;
  if (b->GetCoeff() == 0 || b->HasValue(0))
    return a;  
  if (!(a->IsConstant() && b->IsConstant()) && a->IsEqualToNoCoeff(b)) {
    a->SetCoeff(a->GetCoeff() + b->GetCoeff());
    if (fabs(a->GetCoeff()) < Ev3NearZero()) {
      // simplify to zero - for differences
      Expression zero(0.0);
      return zero;
    } else
      return a;
  }
  // go for it
  if (a->IsLeaf() && a->GetOpType() == CONST &&
      b->IsLeaf() && b->GetOpType() == CONST) {
    // a, b are numbers - add them
    a->SetValue(a->GetValue() + b->GetValue());
    a->SetCoeff(1);
    a->SetExponent(1);
    return a;
  } else if (a->IsLeaf() && a->GetOpType() == VAR &&
	     b->IsLeaf() && b->GetOpType() == VAR &&
	     a->GetVarIndex() == b->GetVarIndex() &&
	     a->GetExponent() == b->GetExponent()) {
    // a, b are the same variable - add coefficients
    a->SetCoeff(a->GetCoeff() + b->GetCoeff());
    return a;
  } else if (a->GetOpType() == SUM && b->GetOpType() != SUM) {    
    // a is a sum and b isn't - just add the term b
    a->DistributeCoeffOverSum();
    Int i = 0;
    bool couldsimplify = false;
    Expression tmp;
    // if t is a leaf and there is a like leaf in this,
    // just add it to the value/coefficient
    if (b->IsLeaf() && b->GetOpType() == CONST) {
      // b is a constant
      for (i = 0; i < a->GetSize(); i++) {
	tmp = a->GetNode(i);
	if (tmp->IsLeaf() && tmp->GetOpType() == CONST) {
	  tmp->SetValue(tmp->GetValue() + b->GetValue() / a->GetCoeff());
	  tmp->SetCoeff(1);
	  tmp->SetExponent(1); // NB: changing tmp should also change a
	  couldsimplify = true;
	  break;
	}
      }
    } else if (b->IsLeaf() && b->GetOpType() == VAR) {
      // b is a variable
      for (i = 0; i < a->GetSize(); i++) {
	if (a->GetNode(i)->IsLeaf() && a->GetNode(i)->GetOpType() == VAR &&
	    b->GetVarIndex() == a->GetNode(i)->GetVarIndex() &&
	    b->GetExponent() == a->GetNode(i)->GetExponent()) {
	  double tc = a->GetNode(i)->GetCoeff() + b->GetCoeff()/a->GetCoeff();
	  // warning: tc could be zero, but it would be cumbersome
	  // to simplify it here - do it in SimplifyConstant
	  a->GetNode(i)->SetCoeff(tc);
	  couldsimplify = true;
	  break;
	}
      }
    } else if (!b->IsLeaf()) {
      // a is a sum, b is a nonleaf, look for a subnode of a similar to b
      for (i = 0; i < a->GetSize(); i++) {
	if (a->GetNode(i)->IsEqualTo(b)) {
	  // found one, add coefficients - notice, as above, coeff could
	  // be zero, but deal with that case in SimplifyConstant
	  a->GetNode(i)->SetCoeff(a->GetNode(i)->GetCoeff() + b->GetCoeff());
	  couldsimplify = true;
	  break;
	}
      }
    }
    if (!couldsimplify) {
      // either could not simplify in steps above, or b is an operator
      a->AddNode(b);
    }
    return a;
  } else if (a->GetOpType() == SUM && b->GetOpType() == SUM) {    
    // a, b are sums - add terms of b to a
    b->DistributeCoeffOverSum();
    Int i = 0;
    Int s = b->GetSize();
    for (i = 0; i < s; i++) {
      a = SumLink(a, b->GetNode(i));
    }
    return a;
  } else if (a->GetOpType() != SUM && b->GetOpType() == SUM) {
    // a is not a sum but b is - transform this into a sum
    b = SumLink(b, a);
    return b;
  } else {
    // all other cases - make new node on top of the addends
    Expression ret;
    

    ret->SetOpType(SUM);
    ret->SetCoeff(1);
    ret->SetExponent(1);
    ret->AddNode(a);
    ret->AddNode(b);

    return ret;
  }
}


// product
Expression ProductLink(Expression a, Expression t) {
  // make a preliminary check
  if (a->GetCoeff() == 0 || t->GetCoeff() == 0 || 
      a->HasValue(0) || t->HasValue(0)) {
    Expression zero(0.0);
    return zero;
  }
  if (a->HasValue(1))
    return t;
  if (t->HasValue(1))
    return a; 
  if (!(a->IsConstant() && t->IsConstant()) && a->IsEqualToNoCoeff(t)) {
    Expression two(2.0);
    a->SetCoeff(a->GetCoeff() * t->GetCoeff());
    return a^two;
  }
  // go for it
  if (a->IsLeaf() && a->GetOpType() == CONST &&
      t->IsLeaf() && t->GetOpType() == CONST) {
    // a, t are numbers - multiply them
    a->SetValue(a->GetValue() * t->GetValue());
    a->SetCoeff(1);
    a->SetExponent(1);
    return a;
  } else if (a->IsLeaf() && a->GetOpType() == VAR &&
	     t->IsLeaf() && t->GetOpType() == VAR &&
	     a->GetVarIndex() == t->GetVarIndex()) {
    // a, t are the same variable - multiply coefficients
    // and add exponents
    a->SetCoeff(a->GetCoeff() * t->GetCoeff());
    a->SetExponent(a->GetExponent() + t->GetExponent());
    return a;
  } else if (t->IsConstant()) {
    // t is constant, set coeff of a
    a->SetCoeff(a->GetCoeff() * t->GetValue());
    a->DistributeCoeffOverSum();
    return a;
  } else if (a->IsConstant()) {
    // a is constant, set coeff of t
    t->SetCoeff(t->GetCoeff() * a->GetValue());
    t->DistributeCoeffOverSum();
    return t;
  } else if (a->GetOpType() == PRODUCT && t->GetOpType() != PRODUCT) {
    // a is a product and t isn't - just multiply the term t
    Int i = 0;
    bool couldsimplify = false;
    if (t->IsLeaf() && t->GetOpType() == VAR) {
      // t is a variable
      Expression tmp;
      for (i = 0; i < a->GetSize(); i++) {
	tmp = a->GetNode(i);
	if (tmp->IsLeaf() && tmp->GetOpType() == VAR &&
	    t->GetVarIndex() == tmp->GetVarIndex()) {
	  // found same variable in a, multiply coeffs and add exponents
	  tmp->SetCoeff(tmp->GetCoeff() * t->GetCoeff());	
	  tmp->SetExponent(tmp->GetExponent() + t->GetExponent());
	  couldsimplify = true;
	  break;
	}
      }
    } 
    // here we shan't try to simplify f*f <-- f^2 (f nonleaf) 
    // because a product of nonleaves is easier to manipulate than 
    // a power (as it adds a depth level)
    if (!couldsimplify) {
      // either could not simplify in steps above, or t is an operator
      a->AddNode(t);
    }
    return a;
  } else if (a->GetOpType() == PRODUCT && t->GetOpType() == PRODUCT) {
    // a, t are products - multiply terms of t to a
    t->DistributeCoeffOverProduct();
    Int i = 0;
    Int s = t->GetSize();
    for (i = 0; i < s; i++) {
      a = ProductLink(a, t->GetNode(i));
    }
    return a;
  } else if (a->GetOpType() != PRODUCT && t->GetOpType() == PRODUCT) {
    // a is not a products but t is - transform this into a product
    t = ProductLink(t, a);
    return t;
  } else {
    // all other cases - make new node on top of the addends
    Expression ret;
    ret->SetCoeff(1);
    ret->SetExponent(1);
    ret->SetOpType(PRODUCT);
    ret->AddNode(a);
    ret->AddNode(t);
    return ret;
  }
}

// fractions:
Expression FractionLink(Expression a, Expression t) 
  throw(ErrDivideByZero) {
  // make a preliminary check
  if (t->GetCoeff() == 0) {
    // divide by zero
    unsigned long mycode(0);
    string myif("Expression Building");
    string myscope("FractionLink");
    string myop("t.GetCoeff()==0");
    string mydesc("Divisor cannot be zero");
    string myinfo(HELPURL);
    string mydiv(NONE);
    throw ErrDivideByZero(mycode,myif,myscope,myop,mydesc,myinfo,mydiv);
  }
  if (a->GetCoeff() == 0 || a->HasValue(0)) {
    // dividend has zero coeff, return zero
    Expression zero(0.0);
    return zero;
  }
  if (t->HasValue(1))
    return a;
  if (!(a->IsConstant() && t->IsConstant()) && a->IsEqualToNoCoeff(t)) {
    // dividend = divisor, return ratio of coefficients
    Expression one(1.0);
    one->SetCoeff(a->GetCoeff() / t->GetCoeff());
    return one;
  }
  // go for it
  if (a->IsLeaf() && a->GetOpType() == CONST &&
      t->IsLeaf() && t->GetOpType() == CONST) {
    // a, t are numbers - divide them
    if (t->GetValue() == 0) {
      unsigned long mycode(0);
      string myif("Expression Building");
      string myscope("FractionLink");
      string myop("t.GetValue()==0");
      string mydesc("Divisor cannot be zero");
      string myinfo(HELPURL);
      string mydiv(NONE);
      throw ErrDivideByZero(mycode,myif,myscope,myop,mydesc,myinfo,mydiv);
    } else {
      a->SetValue(a->GetValue() / t->GetValue());
      a->SetCoeff(1);
      a->SetExponent(1);
      return a;
    }
  } else if (t->HasValue(1)) {
    // divide by constant 1, don't do anything
    return a;
  } else if (t->IsConstant()) {
    // t is constant, set coeff of a
    a->SetCoeff(a->GetCoeff() / t->GetValue());
    a->DistributeCoeffOverSum();
    return a;
  } else if (a->IsVariable() && t->IsVariable() &&
	     a->GetVarIndex() == t->GetVarIndex()) {
    // cx^e / dx^f = (c/d)x^(e-f)
    double te = a->GetExponent() - t->GetExponent();
    double tc = a->GetCoeff() / t->GetCoeff();
    if (fabs(te) < Ev3NearZero()) {
      Expression c(tc);
      return tc;
    }
    a->SetCoeff(tc);
    a->SetExponent(te);
    return a;
  } else if (a->IsVariable() && t->GetOpType() == PRODUCT) {
    // a is a variable, t is a product - see if a appears in t
    // and cancel common term
    // first simplify coeffs of divisor
    t->ConsolidateProductCoeffs();
    // denominator
    if (fabs(t->GetCoeff()) < Ev3NearZero()) {
      // divide by zero
      unsigned long mycode(22);
      string myif("Expression Building");
      string myscope("FractionLink");
      string myop("t->GetCoeff()");
      string mydesc("Divisor cannot be zero");
      string myinfo(HELPURL);
      string mydiv(NONE);
      throw ErrDivideByZero(mycode,myif,myscope,myop,mydesc,myinfo,mydiv);
    }
    if (fabs(a->GetCoeff()) < Ev3NearZero()) {
      Expression zero(0.0);
      return zero;
    }
    double accumulatedcoeff = a->GetCoeff() / t->GetCoeff();
    a->SetCoeff(1.0);
    t->SetCoeff(1.0);
    // now try simplification
    Int i;
    for (i = 0; i < t->GetSize(); i++) {
      if (t->GetNode(i)->GetOpType() == VAR && 
	  a->GetVarIndex() == t->GetNode(i)->GetVarIndex()) {
	double te = a->GetExponent() - t->GetNode(i)->GetExponent();
	if (fabs(te) < Ev3NearZero()) {
	  // exponents are the same, just cancel
	  a->One();
	  t->DeleteNode(i);
	} else if (te > 0) {
	  // numerator remains, cancel denominator
	  a->SetExponent(te);
	  t->DeleteNode(i);
	} else if (te < 0) {
	  // numerator goes to one, denominator remains
	  a->One();
	  t->GetNode(i)->SetExponent(-te);
	}
	// exit loop
	break;
      }
    }
    // check that denominator (t) has more than one operands;
    // if not, bring up a rank level
    if (t->GetSize() == 1) {
      t = t->GetNode(0);
    }
    // build ratio
    Expression ret;
    ret->SetOpType(FRACTION);
    ret->SetCoeff(accumulatedcoeff);
    ret->SetExponent(1);
    ret->AddNode(a);
    ret->AddNode(t);  
    return ret;
  } else if (t->IsVariable() && a->GetOpType() == PRODUCT) {
    // t is a variable, a is a product - see if t appears in a
    // and cancel common term
    // first simplify coeffs of divisor
    a->ConsolidateProductCoeffs();
    // denominator - already checked
    if (fabs(a->GetCoeff()) < Ev3NearZero()) {
      Expression zero(0.0);
      return zero;
    }
    double accumulatedcoeff = a->GetCoeff() / t->GetCoeff();
    a->SetCoeff(1.0);
    t->SetCoeff(1.0);
    // now try simplification
    Int i;
    for (i = 0; i < a->GetSize(); i++) {
      if (a->GetNode(i)->GetOpType() == VAR && 
	  t->GetVarIndex() == a->GetNode(i)->GetVarIndex()) {
	double te = a->GetNode(i)->GetExponent() - t->GetExponent();
	if (fabs(te) < Ev3NearZero()) {
	  // exponents are the same, just cancel
	  t->One();
	  a->DeleteNode(i);
	} else if (te > 0) {
	  // numerator remains, cancel denominator
	  t->One();
	  a->GetNode(i)->SetExponent(te);
	} else if (te < 0) {
	  // numerator goes to one, denominator remains
	  t->SetExponent(-te);
	  a->DeleteNode(i);
	}
	// exit loop
	break;
      }
    }
    // check that numerator (a) has more than one operands;
    // if not, bring up a rank level
    if (a->GetSize() == 1) {
      a = a->GetNode(0);
    }
    // build ratio
    Expression ret;
    ret->SetOpType(FRACTION);
    ret->SetCoeff(accumulatedcoeff);
    ret->SetExponent(1);
    ret->AddNode(a);
    ret->AddNode(t);  
    return ret;
  } else if (a->GetOpType() == PRODUCT && t->GetOpType() == PRODUCT) {
    // a, t are products, try to cancel common terms
    Int i = 0, j = 0;
    double accumulatedcoeff;
    // first simplify coefficients of operands
    a->ConsolidateProductCoeffs();
    t->ConsolidateProductCoeffs();
    // denominator
    if (fabs(t->GetCoeff()) < Ev3NearZero()) {
      // divide by zero
      unsigned long mycode(21);
      string myif("Expression Building");
      string myscope("FractionLink");
      string myop("t->GetCoeff()");
      string mydesc("Divisor cannot be zero");
      string myinfo(HELPURL);
      string mydiv(NONE);
      throw ErrDivideByZero(mycode,myif,myscope,myop,mydesc,myinfo,mydiv);
    }
    if (fabs(a->GetCoeff()) < Ev3NearZero()) {
      Expression zero(0.0);
      return zero;
    }
    // save ratio of coeffs of products
    accumulatedcoeff = a->GetCoeff() / t->GetCoeff();
    a->SetCoeff(1.0);
    t->SetCoeff(1.0);
    // now try simplification
    i = 0;
    bool isnumeratorempty = false;
    bool isdenominatorempty = false;
    int szi = a->GetSize();
    int szj = t->GetSize();
    while(!isnumeratorempty && !isdenominatorempty && i < szi) {
      j = 0;
      while(!isnumeratorempty && !isdenominatorempty && j < szj) {
	if (a->GetNode(i)->IsEqualTo(t->GetNode(j))) {
	  // found like terms i and j
	  a->DeleteNode(i);
	  szi--;
	  if(szi == 0) {
	    isnumeratorempty = true;
	    a->One();
	  }
	  t->DeleteNode(j);
	  szj--;
	  if (szj == 0) {
	    isdenominatorempty = true;
	    t->One();
	  }
	  i--;   // cancel the effect of the later i++
	  break; // go to outer cycle
	} else {
	  j++;
	}
      }	
      i++;
    }
    if (t->HasValue(1)) {
      // denominator is 1, return a
      a->SetCoeff(accumulatedcoeff);
      return a;
    }
    // now construct fraction
    // check that numerator, denominator have more than one operands;
    // if not, bring up a rank level
    if (a->GetSize() == 1) {
      a = a->GetNode(0);
    }
    if (t->GetSize() == 1) {
      t = t->GetNode(0);
    }
    Expression ret;
    ret->SetCoeff(accumulatedcoeff); // already contains coeffs of a, t
    ret->SetExponent(1);
    ret->SetOpType(FRACTION);
    ret->AddNode(a);
    ret->AddNode(t);
    return ret;
  } else {
    Expression ret;
    ret->SetCoeff(a->GetCoeff() / t->GetCoeff());
    a->SetCoeff(1);
    t->SetCoeff(1);
    ret->SetExponent(1);
    ret->SetOpType(FRACTION);
    ret->AddNode(a);
    ret->AddNode(t);
    return ret;
  }
}    

// unary minus:
Expression MinusLink(Expression a) {
  if (a->IsLeaf() && a->GetOpType() == CONST) {
    a->SetValue(- a->GetValue());
    a->SetCoeff(1);
    a->SetExponent(1);    
  } else {
    a->SetCoeff(- a->GetCoeff());
  }
  return a;
}

// binary minus:
Expression DifferenceLink(Expression a, Expression b) {
  if (a->HasValue(0))
    return MinusLink(b);
  if (b->HasValue(0))
    return a;
  return SumLink(a, MinusLink(b));
}

// power:
Expression PowerLink(Expression a, Expression t) {
  // make a preliminary check
  if (a->GetCoeff() == 0) {
    // *this is zero, just return zero
    Expression zero(0.0);
    return zero;
  }
  if (t->HasValue(0.0)) {
    // exponent is 0, just return 1
    Expression one(1.0);
    return one;
  } else if (t->HasValue(1.0)) {
    // exponent is 1, just return a
    return a;
  } 
  if (a->HasValue(0.0)) {
    // base is zero, return 0
    Expression zero(0.0);
    return zero;
  } else if (a->HasValue(1.0)) {
    // base is one, return 1
    Expression one(1.0);
    return one;
  }
  // go for it
  if (a->IsLeaf() && a->GetOpType() == CONST &&
      t->IsLeaf() && t->GetOpType() == CONST) {
    // constant to constant
    a->SetValue(pow(a->GetValue(), t->GetValue()));
    a->SetCoeff(1);
    a->SetExponent(1);
    return a;
  } else if (a->IsLeaf() && a->GetOpType() == VAR &&
	     t->IsLeaf() && t->GetOpType() == CONST) {
    // variable to constant
    a->SetExponent(a->GetExponent() * t->GetValue());
    return a;
  } else {
    Expression ret;
    ret->SetCoeff(1);
    ret->SetExponent(1);
    ret->SetOpType(POWER);
    ret->AddNode(a);
    ret->AddNode(t);
    return ret;
  }
}

Expression LogLink(Expression a) throw(ErrNotPermitted) {
  // make a preliminary check
  if (a->IsZero()) {
    // *this is zero, can't do
    unsigned long mycode(0);
    string myif("Expression Building");
    string myscope("LogLink");
    string myop("IsZero()");
    string mydesc("log(0) is undefined");
    string myinfo(HELPURL);
    throw ErrNotPermitted(mycode,myif,myscope,myop,mydesc,myinfo);
  }
  if (a->IsLessThan(0)) {
    // argument is < zero, can't do
    unsigned long mycode(0);
    string myif("Expression Building");
    string myscope("LogLink");
    string myop("value <= 0");
    string mydesc("log(<=0) is undefined");
    string myinfo(HELPURL);
    throw ErrNotPermitted(mycode,myif,myscope,myop,mydesc,myinfo);
  } 
  // go for it
  if (a->IsLeaf() && a->GetOpType() == CONST) {
    double t = a->GetValue();
    assert(t >= 0);
    a->SetCoeff(1);    
    a->SetValue(log(t));
    a->SetExponent(1);
    a->SetOpType(CONST);
    return a;
  } else {
    Expression ret;
    ret->SetCoeff(1);
    ret->SetExponent(1);
    ret->SetOpType(LOG);
    ret->AddNode(a);
    return ret;
  }
}

Expression ExpLink(Expression a) {
  // go for it
  if (a->IsLeaf() && a->GetOpType() == CONST) {
    a->SetCoeff(1);
    a->SetValue(exp(a->GetValue()));
    a->SetExponent(1);
    a->SetOpType(CONST);
    return a;
  } else {
    Expression ret;
    ret->SetCoeff(1);
    ret->SetExponent(1);
    ret->SetOpType(EXP);
    ret->AddNode(a);
    return ret;
  }
}

Expression SinLink(Expression a) {
  // go for it
  if (a->IsLeaf() && a->GetOpType() == CONST) {
    a->SetCoeff(1);
    a->SetValue(sin(a->GetValue()));
    a->SetExponent(1);
    a->SetOpType(CONST);
    return a;
  } else {
    Expression ret;
    ret->SetCoeff(1);
    ret->SetExponent(1);
    ret->SetOpType(SIN);
    ret->AddNode(a);
    return ret;
  }
}

Expression CosLink(Expression a) {
  // go for it
  if (a->IsLeaf() && a->GetOpType() == CONST) {
    a->SetCoeff(1);
    a->SetValue(cos(a->GetValue()));
    a->SetExponent(1);
    a->SetOpType(CONST);
    return a;
  } else {
    Expression ret;
    ret->SetCoeff(1);
    ret->SetExponent(1);
    ret->SetOpType(COS);
    ret->AddNode(a);
    return ret;
  }
}

Expression TanLink(Expression a) {
  // go for it
  if (a->IsLeaf() && a->GetOpType() == CONST) {
    a->SetCoeff(1);
    a->SetValue(tan(a->GetValue()));
    a->SetExponent(1);
    a->SetOpType(CONST);
    return a;
  } else {
    Expression ret;
    ret->SetCoeff(1);
    ret->SetExponent(1);
    ret->SetOpType(TAN);
    ret->AddNode(a);
    return ret;
  }
}

Expression CotLink(Expression a)  throw(ErrNotPermitted) {
  // make a preliminary check
  if (a->IsZero()) {
    // *this is zero, can't do
    unsigned long mycode(0);
    string myif("Expression Building");
    string myscope("CotLink");
    string myop("IsZero()");
    string mydesc("cot(0) is undefined");
    string myinfo(HELPURL);
    throw ErrNotPermitted(mycode,myif,myscope,myop,mydesc,myinfo);
  }
  // go for it
  if (a->IsLeaf() && a->GetOpType() == CONST) {
    double t = tan(a->GetValue());
    assert(t != 0);
    a->SetCoeff(1);
    a->SetValue(1 / t);
    a->SetExponent(1);
    a->SetOpType(CONST);
    return a;
  } else {
    Expression ret;
    ret->SetCoeff(1);
    ret->SetExponent(1);
    ret->SetOpType(COT);
    ret->AddNode(a);
    return ret;
  }
}

Expression SinhLink(Expression a) {
  // go for it
  if (a->IsLeaf() && a->GetOpType() == CONST) {
    a->SetCoeff(1);
    a->SetValue(sinh(a->GetValue()));
    a->SetExponent(1);
    a->SetOpType(CONST);
    return a;
  } else {
    Expression ret;
    ret->SetCoeff(1);
    ret->SetExponent(1);
    ret->SetOpType(SINH);
    ret->AddNode(a);
    return ret;
  }
}

Expression CoshLink(Expression a) {
  // go for it
  if (a->IsLeaf() && a->GetOpType() == CONST) {
    a->SetCoeff(1);
    a->SetValue(cosh(a->GetValue()));
    a->SetExponent(1);
    a->SetOpType(CONST);
    return a;
  } else {
    Expression ret;
    ret->SetCoeff(1);
    ret->SetExponent(1);
    ret->SetOpType(COSH);
    ret->AddNode(a);
    return ret;
  }
}

Expression TanhLink(Expression a) {
  // go for it
  if (a->IsLeaf() && a->GetOpType() == CONST) {
    a->SetCoeff(1);
    a->SetValue(tanh(a->GetValue()));
    a->SetExponent(1);
    a->SetOpType(CONST);
    return a;
  } else {
    Expression ret;
    ret->SetCoeff(1);
    ret->SetExponent(1);
    ret->SetOpType(TANH);
    ret->AddNode(a);
    return ret;
  }
}

Expression CothLink(Expression a)  throw(ErrNotPermitted) {
  // make a preliminary check
  if (a->IsZero()) {
    // *this is zero, can't do
    unsigned long mycode(0);
    string myif("Expression Building");
    string myscope("CothLink");
    string myop("IsZero()");
    string mydesc("coth(0) is undefined");
    string myinfo(HELPURL);
    throw ErrNotPermitted(mycode,myif,myscope,myop,mydesc,myinfo);
  }
  // go for it
  if (a->IsLeaf() && a->GetOpType() == CONST) {
    double t = tanh(a->GetValue());
    assert(t != 0);
    a->SetCoeff(1);
    a->SetValue(1 / t);
    a->SetExponent(1);
    a->SetOpType(CONST);
    return a;
  } else {
    Expression ret;
    ret->SetCoeff(1);
    ret->SetExponent(1);
    ret->SetOpType(COTH);
    ret->AddNode(a);
    return ret;
  }
}

Expression SqrtLink(Expression a) throw(ErrNotPermitted) {
  // make a preliminary check
  if (a->IsLessThan(0) && !a->HasValue(0)) {
    // argument is < zero, can't do
    unsigned long mycode(0);
    string myif("Expression Building");
    string myscope("SqrtLink");
    string myop("value < 0");
    string mydesc("sqrt(<0) is complex, can't do");
    string myinfo(HELPURL);
    throw ErrNotPermitted(mycode,myif,myscope,myop,mydesc,myinfo);
  } 
  // go for it
  if (a->IsLeaf() && a->GetOpType() == CONST) {
    double t = a->GetValue();
    assert(t >= 0);
    a->SetCoeff(1);
    a->SetValue(sqrt(t));
    a->SetExponent(1);
    a->SetOpType(CONST);
    return a;
  } else {
    Expression ret;
    ret->SetCoeff(1);
    ret->SetExponent(1);
    ret->SetOpType(SQRT);
    ret->AddNode(a);
    return ret;
  }
}


/***************** differentiation ******************/

// differentiate w.r.t. variable varindex
Expression Diff(const Expression& ac, Int vi) {
  Expression ret(DiffNoSimplify(ac, vi));
  Simplify(&ret);
  return ret;
}

Expression DiffNoSimplify(const Expression& ac, Int vi) {
  Expression a;
  a.SetToCopyOf(ac);
  Expression zero(0.0);
  Expression c(1.0);

  if (a->DependsOnVariable(vi)) {
    if (a->IsLeaf()) {
      
      if (a->GetOpType() == CONST || a->GetVarIndex() != vi) {
	// safety check
	cerr << "Expression::Diff: warning: this node should "
	     << "not diff to zero\n";
	return zero;
      } else {
	// variable vi, check exponent
	if (a->GetExponent() == 0) {
	  // exponent is zero, node is actually constant
	  return zero;
	} else if (a->GetExponent() == 1) {
	  // exponent is one, node is variable
	  c->SetValue(a->GetCoeff());

	  return c;
	} else {
	  // for all other cases, apply rule x^c = c*x^(c-1)
	  double expon = a->GetExponent();
	  Expression ret(a.Copy());
	  ret->SetExponent(expon - 1);
	  ret->SetCoeff(ret->GetCoeff() * expon);
	  
	  return ret;
	}
      }
    } else {

      // non-leaf node. build derivative.
      int op = a->GetOpType();
      Int sz = a->GetSize();
      double opcoeff = a->GetCoeff();
      if (sz == 0) {
	throw ErrNotPermitted(10, "Expression", "Diff", "GetSize() == 0",
			      "non-leaf node can't have size 0", HELPURL);
      }
      Int i, j;
      Expression ret(0.0);
      Expression tmp(1.0);
      Expression two(2.0);
      switch(op) {
      case SUM:
	ret = Diff(a->GetNode(0), vi); // f_0'
	for(i = 1; i < sz; i++) {
	  ret = ret + Diff(a->GetNode(i), vi); // ... + g_i'
	}
	break;
      case DIFFERENCE:
	ret = Diff(a->GetNode(0), vi);  // f_0'
	for(i = 1; i < sz; i++) {
	  ret = ret - Diff(a->GetNode(i), vi);  // ... - g_i'
	}
	break;
      case PRODUCT:
	if (sz == 1) {
	  // warn about product with one operand 
	  cerr << "Expression::Diff: warning: product with 1 operand "
	       << "should not occur\n";
	} 
	ret = Diff(a->GetNode(0), vi);  // f_0'
	for(j = 1; j < sz; j++) {
	  // get copies, not references
	  ret = ret * a->GetCopyOfNode(j); // ... * f_i[i!=0]
	}
	tmp->One(); // reset temporary to 1.0
	for(i = 1; i < sz; i++) {
	  tmp = Diff(a->GetNode(i), vi); // tmp = f_i'
	  for(j = 0; j < sz; j++) {
	    if (j != i) 
	      // get references, and copy later (in sum)
	      tmp = tmp * a->GetNode(j); // ... * f_j[i!=i]
	  }
	  ret = ret + tmp.Copy();  // ... + tmp
	  tmp->One(); // reset temporary to 1.0 
	}
	break;
      case FRACTION:
	if (sz != 2) {
	  // check that there are exactly two operands
	  throw ErrNotPermitted(11, "Expression", "Diff", "GetSize() != 2",
				"fraction must have exactly 2 operands", 
				HELPURL);
	}
	if (a->GetNode(1)->IsZero()) {
	  // check denominator is not zero
	  throw ErrDivideByZero(20, "Expression", "Diff", 
				"GetNode(1)->IsZero()", 
				"cannot divide by zero", HELPURL, 
				a->GetNode(1)->ToString());
	}
	tmp->One();
	ret = Diff(a->GetNode(0), vi); // f'
	ret = ret * a->GetCopyOfNode(1);  // f'g
	// can dispense from using GetCopyOf because tmp gets copied anyway
	tmp = a->GetNode(0);           // tmp = f
	tmp = tmp * Diff(a->GetNode(1), vi);  // tmp = fg'
	ret = ret - tmp.Copy();    // f'g - fg'
	tmp->One();
	tmp = a->GetNode(1);  // tmp = g
	tmp = tmp ^ two; // g^2
	// can dispense from using copy here - tmp is not used thereafter
	// and when tmp is deleted, its subnodes are not automatically
	// deleted unless reference counter is zero - which won't be.
	ret = ret / tmp;   // (f'g - fg')/g^2
	break;
      case POWER:
	if (sz != 2) {
	  // check that there are exactly two operands
	  throw ErrNotPermitted(12, "Expression", "Diff", "GetSize() != 2",
				"power must have exactly 2 operands", 
				HELPURL);
	}
	// check exponent
	if (a->GetNode(1)->IsZero()) {
	  // exponent is zero, whole node is one, diff is zero
	  ret->Zero();
	} else if (a->GetNode(1)->HasValue(1.0)) {
	  // exponent is one, whole node is first operand, diff
	  // is diff of first operand
	  ret = Diff(a->GetNode(0), vi);
	} else if (a->GetNode(1)->HasValue(2.0)) { 
	  // exponent is two, diff is 2 * first op * diff of first operand
	  ret = Diff(a->GetNode(0), vi);  // f'
	  ret = ret * a->GetCopyOfNode(0);   // f'f
	  ret->SetCoeff(ret->GetCoeff() * 2.0);  // 2f'f
	} else {
	  // all other cases
	  if (a->GetNode(1)->IsConstant()) {
	    Expression tmp2;

	    // numeric exponent != 0,1,2
	    
	    ret = Diff(a->GetNode(0), vi); // f'
	    tmp = a->GetCopyOfNode(0);     // f
	    tmp2 = a->GetCopyOfNode(1);
	    tmp2->ConsolidateValue();
	    tmp->SetCoeff(tmp->GetCoeff() * tmp2->GetValue());
	    tmp->SetExponent(tmp2->GetValue()-1);

	    // can dispense from using copy here - Diff returns copies anyway.
	    // when temporary is deleted, its subnodes are not automatically
	    // deleted unless their reference counter is zero - which won't be.
	    ret = ret * tmp; // f'(cf^(c-1))
	  } else {
	    // symbolic exponent f^g    
	    ret = a->GetCopyOfNode(0); // f
	    ret = Log(ret); // log(f)
	    // can dispense from using copy here - Diff returns copies anyway.
	    // when temporary is deleted, its subnodes are not automatically
	    // deleted unless their reference counter is zero - which won't be.
	    ret = ret * Diff(a->GetNode(1), vi); // g' log(f)
	    tmp = Diff(a->GetNode(0), vi);  // f'
	    tmp = tmp * a->GetCopyOfNode(1);   // gf'
	    tmp = tmp / a->GetCopyOfNode(0);    // gf'/f
	    // can dispense from using copy here - tmp is not used thereafter
	    // and when tmp is deleted, its subnodes are not automatically
	    // deleted unless their reference counter is zero - which won't be.
	    ret = ret / tmp;           // g'log(f) + gf'/f
	    ret = ret * a.Copy();        // (g'log(f) + gf'/f)(f^g)
	  }
	}
	break;
      case MINUS:
	if (sz != 1) {
	  // check that there is exactly one operand
	  throw ErrNotPermitted(13, "Expression", "Diff", "GetSize() != 1",
				"unary minus must have exactly 1 operand", 
				HELPURL);
	}
	ret = Diff(a->GetNode(0), vi);
	ret->SetCoeff(- ret->GetCoeff());
	break;
      case LOG:
	if (sz != 1) {
	  // check that there is exactly one operand
	  throw ErrNotPermitted(14, "Expression", "Diff", "GetSize() != 1",
				"log must have exactly 1 operand", 
				HELPURL);
	}
	if (a->GetNode(0)->IsLessThan(0)) {
	  throw ErrNotPermitted(15, "Expression", "Diff", "arg <= 0",
				"log argument must be symbolic or positive",
				HELPURL);
	}
	ret = Diff(a->GetNode(0), vi);  // f'
	ret = ret / a->GetCopyOfNode(0);  // f'/f       
	break;
      case EXP:
	if (sz != 1) {
	  // check that there is exactly one operand
	  throw ErrNotPermitted(16, "Expression", "Diff", "GetSize() != 1",
				"exp must have exactly 1 operand", 
				HELPURL);
	}
	ret = Diff(a->GetNode(0), vi);  // f'
	ret = ret * a.Copy();  // f' e^f
	break;
      case SIN:
	if (sz != 1) {
	  // check that there is exactly one operand
	  throw ErrNotPermitted(17, "Expression", "Diff", "GetSize() != 1",
				"sin must have exactly 1 operand", 
				HELPURL);
	}
	ret = Diff(a->GetNode(0), vi) * Cos(a->GetCopyOfNode(0));  // f' cos(f)
	break;
      case COS:
	if (sz != 1) {
	  // check that there is exactly one operand
	  throw ErrNotPermitted(18, "Expression", "Diff", "GetSize() != 1",
				"cos must have exactly 1 operand", 
				HELPURL);
	}
	ret = -Diff(a->GetNode(0), vi) * Sin(a->GetCopyOfNode(0)); // -f'sin(f)
	break;
      case TAN:
	if (sz != 1) {
	  // check that there is exactly one operand
	  throw ErrNotPermitted(19, "Expression", "Diff", "GetSize() != 1",
				"tan must have exactly 1 operand", 
				HELPURL);
	}
	ret = a.Copy();  // tan(f)
	ret = ret ^ two;      // tan(f)^2
	c->One();
	ret = ret + c;         // tan(f)^2 + 1
	ret = ret * Diff(a->GetNode(0), vi);    // f' * (tan(f)^2 + 1)
	break;
      default:
	// missing cases: COT, SINH, COSH, TANH, COTH, SQRT
	break;
      }
      ret->SetCoeff(ret->GetCoeff() * opcoeff);
      
      return ret;
    }
  } else {
    return zero;
  }
  return zero;
}

/************************ simplifications **********************/

bool TrigSimp(Expression a) {
  Int i = 0;
  Int ret = 0;
  bool bret = false;
  bool ischanged = false;
  // first, recurse over all subnodes of a
  for(i = 0; i < a->GetSize(); i++) {
    ischanged = TrigSimp(a->GetNode(i));
    if (!bret && ischanged) {
      bret = true;
    }
  }
  // now try simplification on a
  if (a->GetOpType() == SUM && a->GetSize() > 1) {    
    // try to look for a sin^2 and a cos^2
    Int sinpos = -1;
    Int cospos = -1;
    Int sinpossimple = -1;
    Int cospossimple = -1;
    for (i = 0; i < a->GetSize(); i++) {
      // cycle over subnodes
      if (a->GetNode(i)->GetOpType() == POWER) {
	if (a->GetNode(i)->GetNode(0)->GetOpType() == SIN &&
	    a->GetNode(i)->GetNode(1)->HasValue(2))
	  sinpos = i;
      }
      if (a->GetNode(i)->GetOpType() == POWER) {
	if (a->GetNode(i)->GetNode(0)->GetOpType() == COS &&
	    a->GetNode(i)->GetNode(1)->HasValue(2))
	  cospos = i;
      }
      if (a->GetNode(i)->GetOpType() == SIN &&
	  a->GetNode(i)->GetExponent() == 2) {
	sinpossimple = i;
      } 
      if (a->GetNode(i)->GetOpType() == COS &&
	  a->GetNode(i)->GetExponent() == 2) {
	cospossimple = i;
      }
    }
    if (sinpos != -1 && cospos != -1) {
      // found both, check their arguments
      if (a->GetNode(sinpos)->GetNode(0)->GetNode(0)->IsEqualTo
	  (a->GetNode(cospos)->GetNode(0)->GetNode(0))) {
	ret++; // augment simplification counter
	bret = true;
	// arguments are equal, can simplify
	Int f = sinpos < cospos ? sinpos : cospos; // last to delete
	Int l = sinpos > cospos ? sinpos : cospos; // first to delete
	a->DeleteNode(l); 
	a->DeleteNode(f); 
	// verify that there are still some nodes left
	if (a->GetSize() == 0) {
	  // there aren't any, set a to one
	  a->One();
	} else {
	  // there are some, check whether there is a constant part
	  // we can add the 1 to
	  bool addflag = false;
	  for (i = 0; i < a->GetSize(); i++) {
	    if (a->GetNode(i)->IsConstant()) {
	      // yes there is
	      a->GetNode(i)->SetValue(a->GetNode(i)->GetSimpleValue() + 1);
	      addflag = true;
	      break;
	    }
	  }
	  if (!addflag) {
	    // no there wasn't, add it as a symbolic node
	    Expression one(1.0);
	    a->AddNode(one);
	  }
	  // check that there is more than just one node
	  if (a->GetSize() == 1) {
	    // only one node, shift everything one rank level up
	    a = a->GetNode(0);
	  }
	}
      }
    }
    if (sinpossimple != -1 && cospossimple != -1) {
      if (a->GetNode(sinpossimple)->GetNode(0)->IsEqualTo
	  (a->GetNode(cospossimple)->GetNode(0))) {
	ret++;
	bret = true;
	// arguments are equal, can simplify
	Int f = sinpossimple < cospossimple ? sinpossimple : cospossimple; 
	Int l = sinpossimple > cospossimple ? sinpossimple : cospossimple; 
	a->DeleteNode(l); 
	a->DeleteNode(f); 
	// verify that there are still some nodes left
	if (a->GetSize() == 0) {
	  // there aren't any, set a to one
	  a->One();
	} else {
	  // there are some, check whether there is a constant part
	  // we can add the 1 to
	  bool addflag = false;
	  for (i = 0; i < a->GetSize(); i++) {
	    if (a->GetNode(i)->IsConstant()) {
	      // yes there is
	      a->GetNode(i)->SetValue(a->GetNode(i)->GetSimpleValue() + 1);
	      addflag = true;
	      break;
	    }
	  }
	  if (!addflag) {
	    // no there wasn't, add it as a symbolic node
	    Expression one(1.0);
	    a->AddNode(one);
	  }
	  // check that there is more than just one node
	  if (a->GetSize() == 1) {
	    // only one node, shift everything one rank level up
	    a = a->GetNode(0);
	  }
	}
      }
    }
  } 
  if (ret > 0)
    bret = true;
  return bret;
}

// generic simplification with modification of the expression
bool Simplify(Expression* a) {
  bool changedflag = false;
  bool ret = false;
  bool goonflag = true;
  while(goonflag) {
    goonflag = false;
    (*a)->ConsolidateProductCoeffs();
    (*a)->DistributeCoeffOverSum();
    ret = DifferenceToSum(a);
    if (ret) {
      changedflag = true;
      goonflag = true;
    }
    ret = SimplifyConstant(a);
    if (ret) {
      changedflag = true;
      goonflag = true;
    }
    ret = CompactProducts(a);
    if (ret) {
      changedflag = true;  
      goonflag = true;
    }
    ret = CompactLinearPart(a);
    if (ret) {
      changedflag = true;  
      goonflag = true;
    }    
    ret = SimplifyRecursive(a);
    if (ret) {
      changedflag = true;
      goonflag = true;
    }
    ret = TrigSimp(*a);
    if (ret) {
      changedflag = true;
      goonflag = true;
    }
  }
  return changedflag;
}

// call after DifferenceToSum
bool SimplifyConstant(Expression* a) {
  bool ret = false;
  Expression one(1.0);
  Expression zero(0.0);
  if ((*a)->GetExponent() == 0) {
    RecursiveDestroy(a);
    a->SetTo(one);
    ret = true;
  } else if ((*a)->GetCoeff() == 0) {
    RecursiveDestroy(a);
    a->SetTo(zero);
    ret = true;
  } else {
    int i = 0;
    int op, ops;
    op = (*a)->GetOpType();
    bool ischanged = false;
    int sz = (*a)->GetSize();
    while (i < sz) {
      // simplify each of the terms
      ischanged = SimplifyConstant((*a)->GetNodePtr(i));
      if (!ret && ischanged) {
	ret = true;
      }
      i++;
    }
    i = 0;
    while (i < sz) {
      // simplify this operand as a whole
      ops = (*a)->GetNode(i)->GetOpType();
      switch(op) {
      case SUM:
	if (ops == CONST) {
	  if ((*a)->GetNode(i)->GetValue() == 0) {
	    (*a)->DeleteNode(i);
	    ret = true;
	    sz--;
	    if (sz == 1) {
	      a->SetTo((*a)->GetNode(0));
	      i = 0;
	      sz = (*a)->GetSize();
	    }
	  } else {
	    i++;
	  }
	} else {
	  i++;
	}
	break;
      case PRODUCT:
	if (ops == CONST) {
	  if ((*a)->GetNode(i)->GetValue() == 1) {
	    (*a)->DeleteNode(i);
	    ret = true;
	    sz--;
	    if (sz == 1) {
	      a->SetTo((*a)->GetNode(0));
	      i = 0;
	      sz = (*a)->GetSize();
	    }
	  } else if ((*a)->GetNode(i)->GetValue() == 0) {
	    RecursiveDestroy(a);
	    ret = true;
	    a->SetTo(zero);
	    sz = 0;
	  } else {
	    i++;
	  }
	} else {
	  i++;
	}
	break;
      case FRACTION:
	if (ops == CONST && i == 1) {
	  a->SetTo((*a)->GetNode(0));
	  ret = true;
	  sz--;
	} else {
	  i++;
	}
	if (sz >= 2 && (*a)->GetNode(0)->IsConstant() && 
	    (*a)->GetNode(1)->IsConstant()) {
	  double d = (*a)->GetNode(1)->GetValue();
	  if (d == 0) {
	    throw ErrDivideByZero(23, "Expression", "SimplifyConstant", 
				  "d==0", 
				  "cannot divide by zero", HELPURL, 
				  (*a)->ToString());
	  } 
	  ret = true;
	  a->SetTo((*a)->GetNode(0));
	  (*a)->SetValue((*a)->GetValue() / d);
	  (*a)->SetCoeff(1);
	  (*a)->SetExponent(1);
	  sz = 0;
	}
	break;
      case POWER:
	if (sz >= 2 && (*a)->GetNode(0)->IsConstant() && 
	    (*a)->GetNode(1)->IsConstant()) {
	  double d = (*a)->GetNode(1)->GetValue();
	  ret = true;
	  a->SetTo((*a)->GetNode(0));
	  (*a)->SetValue(pow((*a)->GetValue(), d));
	  (*a)->SetCoeff(1);
	  (*a)->SetExponent(1);
	  sz = 0;
	} else {
	  i++;
	}
	break;
      case LOG:
	if ((*a)->GetNode(0)->IsConstant()) {
	  double d = (*a)->GetNode(0)->GetValue();
	  if (d <= 0) {
	    throw ErrNotPermitted(24, "Expression", "SimplifyConstant", 
				  "d<=0", 
				  "log of nonpositive not allowed", HELPURL);
	  }
	  ret = true;
	  a->SetTo((*a)->GetNode(0));
	  (*a)->SetValue(log(d));
	  (*a)->SetCoeff(1);
	  (*a)->SetExponent(1);
	  sz = 0;
	} else {
	  i++;
	}
	break;
      case EXP:
	if ((*a)->GetNode(0)->IsConstant()) {
	  double d = (*a)->GetNode(0)->GetValue();
	  ret = true;
	  a->SetTo((*a)->GetNode(0));
	  (*a)->SetValue(exp(d));
	  (*a)->SetCoeff(1);
	  (*a)->SetExponent(1);
	  sz = 0;
	} else {
	  i++;
	}
	break;
      case SIN:
	if ((*a)->GetNode(0)->IsConstant()) {
	  double d = (*a)->GetNode(0)->GetValue();
	  ret = true;
	  a->SetTo((*a)->GetNode(0));
	  (*a)->SetValue(sin(d));
	  (*a)->SetCoeff(1);
	  (*a)->SetExponent(1);
	  sz = 0;
	} else {
	  i++;
	}
	break;
      case COS:
	if ((*a)->GetNode(0)->IsConstant()) {
	  double d = (*a)->GetNode(0)->GetValue();
	  ret = true;
	  a->SetTo((*a)->GetNode(0));
	  (*a)->SetValue(cos(d));
	  (*a)->SetCoeff(1);
	  (*a)->SetExponent(1);
	  sz = 0;
	} else {
	  i++;
	}
	break;
      case TAN:
	if ((*a)->GetNode(0)->IsConstant()) {
	  double d = (*a)->GetNode(0)->GetValue();
	  ret = true;
	  a->SetTo((*a)->GetNode(0));
	  (*a)->SetValue(tan(d));
	  (*a)->SetCoeff(1);
	  (*a)->SetExponent(1);
	  sz = 0;
	} else {
	  i++;
	}
	break;
      case COT:
	if ((*a)->GetNode(0)->IsConstant()) {
	  double d = (*a)->GetNode(0)->GetValue();
	  ret = true;
	  a->SetTo((*a)->GetNode(0));
	  (*a)->SetValue(1/tan(d));
	  (*a)->SetCoeff(1);
	  (*a)->SetExponent(1);
	  sz = 0;
	} else {
	  i++;
	}
	break;
      case SINH:
	if ((*a)->GetNode(0)->IsConstant()) {
	  double d = (*a)->GetNode(0)->GetValue();
	  ret = true;
	  a->SetTo((*a)->GetNode(0));
	  (*a)->SetValue(sinh(d));
	  (*a)->SetCoeff(1);
	  (*a)->SetExponent(1);
	  sz = 0;
	} else {
	  i++;
	}
	break;
      case COSH:
	if ((*a)->GetNode(0)->IsConstant()) {
	  double d = (*a)->GetNode(0)->GetValue();
	  ret = true;
	  a->SetTo((*a)->GetNode(0));
	  (*a)->SetValue(cosh(d));
	  (*a)->SetCoeff(1);
	  (*a)->SetExponent(1);
	  sz = 0;
	} else {
	  i++;
	}
	break;
      case TANH:
	if ((*a)->GetNode(0)->IsConstant()) {
	  double d = (*a)->GetNode(0)->GetValue();
	  ret = true;
	  a->SetTo((*a)->GetNode(0));
	  (*a)->SetValue(tanh(d));
	  (*a)->SetCoeff(1);
	  (*a)->SetExponent(1);
	  sz = 0;
	} else {
	  i++;
	}
	break;
      case COTH:
	if ((*a)->GetNode(0)->IsConstant()) {
	  double d = (*a)->GetNode(0)->GetValue();
	  ret = true;
	  a->SetTo((*a)->GetNode(0));
	  (*a)->SetValue(1/tanh(d));
	  (*a)->SetCoeff(1);
	  (*a)->SetExponent(1);
	  sz = 0;
	} else {
	  i++;
	}
	break;
      case SQRT:
	if ((*a)->GetNode(0)->IsConstant()) {
	  double d = (*a)->GetNode(0)->GetValue();
	  if (d <= 0) {
	    throw ErrNotPermitted(25, "Expression", "SimplifyConstant", 
				  "d<=0", 
				  "sqrt of nonpositive not allowed", HELPURL);
	  }
	  ret = true;
	  a->SetTo((*a)->GetNode(0));
	  (*a)->SetValue(sqrt(d));
	  (*a)->SetCoeff(1);
	  (*a)->SetExponent(1);
	  sz = 0;
	} else {
	  i++;
	}
	break;
      default:
	i++;
	break;
      }
    }
  }
  return ret;
}

bool SimplifyRecursive(Expression* a) {
  bool ret = false;
  bool ischanged = false;
  // signals whether we're dealing with 
  // -1 : nothing
  // 0 : constants
  // 1 : linear variables
  // 2 : a variable raised to an exponent different from 1
  // 3 : any other more complex node
  if (!(*a)->IsLeaf()) {
    int op = (*a)->GetOpType();
    Expression t1, t2;
    int i, j;
    for(i = 0; i < (*a)->GetSize(); i++) {
      ischanged = SimplifyRecursive((*a)->GetNodePtr(i));
      if (!ret && ischanged)
	ret = true;
    }
    int status = -1;
    int prestatus = -1;
    double consolidated[4] = {0, 0, 0, 0};
    double expon = 0;
    double preexpon = 0;
    double c = 0;
    int prevarindex = -1;
    int varindex = -1;
    int firstvarindex = -1;
    int firstconstindex = -1;
    int sz = (*a)->GetSize();
    Expression one(1.0);

    switch(op) {
    case SUM:
      i = 0;
      while(i < sz) {
	// work out which status we're in
	if ((*a)->GetNode(i)->IsConstant()) {
	  if (status == -1 || firstconstindex == -1) {
	    firstconstindex = i;
	  }
	  // constant
	  status = 0;
	} else if ((*a)->GetNode(i)->IsVariable() && 
		   (*a)->GetNode(i)->GetExponent() == 1) {
	  // variable
	  status = 1;
	} else if ((*a)->GetNode(i)->IsVariable() && 
		   (*a)->GetNode(i)->GetExponent() != 1) {
	  // variable raised to power
	  status = 2;
	} else {
	  // other term
	  status = 3;
	}
	if (status == 0) {
	  // constant
	  consolidated[status] += (*a)->GetNode(i)->GetValue();
	  (*a)->GetNode(firstconstindex)->SetValue(consolidated[status]);
	  (*a)->GetNode(firstconstindex)->SetCoeff(1);
	  (*a)->GetNode(firstconstindex)->SetExponent(1);
	  if (prestatus == 0) {
	    (*a)->DeleteNode(i);
	    ret = true;
	    sz--;
	    if (sz == 1) {
	      a->SetTo((*a)->GetNode(0));
	      i = 0;
	      sz = (*a)->GetSize();
	    }
	  } else {
	    i++;
	  }
	} else if (status == 1) {
	  // variable
	  varindex = (*a)->GetNode(i)->GetVarIndex();
	  c = (*a)->GetNode(i)->GetCoeff();
	  if (varindex != prevarindex) {
	    firstvarindex = i;
	    consolidated[status] = c;
	    i++;
	  } else {
	    consolidated[status] += c;
	    (*a)->GetNode(firstvarindex)->SetCoeff(consolidated[status]);
	    ret = true;
	    (*a)->DeleteNode(i);
	    sz--;
	    if (sz == 1) {
	      a->SetTo((*a)->GetNode(0));
	      i = 0;
	      sz = (*a)->GetSize();
	    }
	  }
	  prevarindex = varindex;
	} else if (status == 2) {
	  // variable raised to power
	  varindex = (*a)->GetNode(i)->GetVarIndex();
	  expon = (*a)->GetNode(i)->GetExponent();
	  c = (*a)->GetNode(i)->GetCoeff();
	  if (expon != preexpon || varindex != prevarindex) {
	    firstvarindex = i;
	    consolidated[status] = c;
	    i++;
	  } else {
	    consolidated[status] += c;
	    (*a)->GetNode(firstvarindex)->SetCoeff(consolidated[status]);
	    ret = true;
	    (*a)->DeleteNode(i);
	    sz--;
	    if (sz == 1) {
	      a->SetTo((*a)->GetNode(0));
	      i = 0;
	      sz = (*a)->GetSize();
	    }
	  }
	  preexpon = expon;
	  prevarindex = varindex;
	} else if (status == 3) {
	  // other term
	  c = (*a)->GetNode(i)->GetCoeff();
	  firstvarindex = i;
	  consolidated[status] = c;
	  j = i + 1;
	  while(j < sz) {
	    if ((*a)->GetNode(i)->IsEqualToNoCoeff((*a)->GetNode(j))) {
	      consolidated[status] += c;
	      ret = true;
	      (*a)->GetNode(firstvarindex)->SetCoeff(consolidated[status]);
	      RecursiveDestroy((*a)->GetNodePtr(j));
	      (*a)->DeleteNode(j);
	      sz--;
	      if (sz == 1) {
		a->SetTo((*a)->GetNode(0));
		j = 0;
		sz = (*a)->GetSize();
	      }
	    } else {
	      j++;
	    }
	  }
	  i++;
	} else {
	  // should never happen, but anyway...
	  i++;
	}
	// update status of last iteration
	prestatus = status;
      }
      break;
    case PRODUCT:
      i = 0;
      prevarindex = -1;
      consolidated[0] = 1;
      expon = 0;
      while(i < sz) {
	if ((*a)->GetNode(i)->IsVariable()) {
	  varindex = (*a)->GetNode(i)->GetVarIndex();
	  if (varindex != prevarindex) {
	    firstvarindex = i;
	    consolidated[0] = (*a)->GetNode(i)->GetCoeff();
	    expon = (*a)->GetNode(i)->GetExponent();
	    i++;
	  } else {
	    consolidated[0] *= (*a)->GetNode(i)->GetCoeff();
	    expon += (*a)->GetNode(i)->GetExponent();
	    (*a)->GetNode(firstvarindex)->SetCoeff(consolidated[0]);
	    (*a)->GetNode(firstvarindex)->SetExponent(expon);
	    (*a)->DeleteNode(i);
	    ret = true;
	    sz--;
	    if (sz == 1) {
	      a->SetTo((*a)->GetNode(0));
	      i = 0;
	      sz = (*a)->GetSize();
	    }
	  }
	} else if (!(*a)->GetNode(i)->IsLeaf()) {
	  // WARNING: work to be done
	  // not going to do the same as in sum just yet, maybe future
	  // work - transform expr * expr in expr^2 when expr not a variable
	  i++;
	}
      }
      break;
    case FRACTION:
      if ((*a)->GetNode(0)->IsEqualTo((*a)->GetNode(1))) {
	// f(x)/f(x)
	RecursiveDestroy(a);
	a->SetTo(one);
	ret = true;
	sz = 0;
      } else {
	// try to simplify denominator by one of numerator factors
	if ((*a)->GetNode(0)->GetOpType() == PRODUCT) {
	  for(j = 0; j < (*a)->GetNode(0)->GetSize(); j++) {
	    if ((*a)->GetNode(1)->IsEqualTo((*a)->GetNode(0)->GetNode(j))) {
	      a->SetTo((*a)->GetNode(0));
	      (*a)->DeleteNode(j);	      
	      ret = true;
	      sz = 0;
	      break;
	    } 
	  }
	}
	// do the opposite
	if (sz > 0 && (*a)->GetNode(1)->GetOpType() == PRODUCT) {
	  for(j = 0; j < (*a)->GetNode(1)->GetSize(); j++) {
	    if ((*a)->GetNode(0)->IsEqualTo((*a)->GetNode(1)->GetNode(j))) {
	      (*a)->GetNode(0).SetTo(one);
	      (*a)->GetNode(1)->DeleteNode(j);
	      ret = true;
	      sz = 0;
	      break;
	    } 
	  }
	} 
	if (sz > 0 && (*a)->GetNode(0)->GetOpType() == PRODUCT &&
	    (*a)->GetNode(1)->GetOpType() == PRODUCT) {
	  // both num and denom. are products, try and find common factors
	  int k = 0;
	  int sz1, sz2;
	  j = 0;
	  sz1 = (*a)->GetNode(0)->GetSize();
	  sz2 = (*a)->GetNode(1)->GetSize();	
	  while (j < sz1) {
	    k = 0;
	    while (k < sz2) {
	      if ((*a)->GetNode(0)->GetNode(j)->IsEqualTo
		  ((*a)->GetNode(1)->GetNode(k))) {
		(*a)->GetNode(0)->DeleteNode(j);
		(*a)->GetNode(1)->DeleteNode(k);
		ret = true;
		sz1--;
		if (sz1 == 0) {
		  // numerator empty, replace with 1
		  (*a)->GetNode(0)->One();
		}
		sz2--;
		if (sz2 == 0) {
		  // denominator empty, node becomes num.
		  a->SetTo((*a)->GetNode(0));
		}
		if (sz1 == 0 && sz2 == 0) {
		  // 1/1, simplify
		  (*a)->One();
		}
		if (sz1 == 0 || sz2 == 0) {
		  // either num. or den. 1, exit loop
		  sz1 = 0;
		  sz2 = 0;
		  break;
		}
		j--;
	      } else {
		k++;
	      }
	    }
	    j++;
	  }
	}
      }
      sz = 0;
      break;
    case POWER:
      if (sz == 2 && 
	  (*a)->GetNode(0)->IsVariable() &&
	  (*a)->GetNode(1)->IsConstant()) {
	// case var^const, transform in variable with an exponent
	(*a)->GetNode(0)->SetExponent((*a)->GetNode(1)->GetValue());
	(*a)->DeleteNode(1);
	a->SetTo((*a)->GetNode(0));
      }
      break;
    default:
      break;
    }
  }
  return ret;
}

bool DifferenceToSum(Expression* a) {
  bool ret = false;
  double d, e;
  if (!(*a)->IsLeaf()) {
    if (((*a)->GetOpType() == SUM || (*a)->GetOpType() == DIFFERENCE) && 
	(*a)->GetSize() == 1) {
      DifferenceToSum((*a)->GetNodePtr(0));
      // replace a with its child
      a->SetTo((*a)->GetNode(0));
      ret = true;
    }
    if ((*a)->GetOpType() == DIFFERENCE) {
      int i;
      (*a)->SetOpType(SUM);
      for(i = 1; i < (*a)->GetSize(); i++) { 
	// start from node 1 not 0 because a difference is +op0 -op1 -op2 ...
	(*a)->GetNode(i)->SetCoeff(- (*a)->GetNode(i)->GetCoeff());
      }
    } else if ((*a)->GetOpType() == MINUS) {
      d = (*a)->GetCoeff();
      e = (*a)->GetExponent();
      if (is_even(e)) {
	// replace a with its child and adjust coeff
	a->SetTo((*a)->GetNode(0));
	(*a)->SetCoeff((*a)->GetCoeff() * d); // since exponent is even, +
	(*a)->SetExponent((*a)->GetExponent() * e);
	ret = true;
      } else if (is_odd(e)) {
	// replace a with its child and adjust coeff
	a->SetTo((*a)->GetNode(0));
	(*a)->SetCoeff(- (*a)->GetCoeff() * d); // since exponent is odd, -
	(*a)->SetExponent((*a)->GetExponent() * e);
	ret = true;
      }
    } else if ((*a)->GetOpType() == PLUS) {
      // replace a with its child
      a->SetTo((*a)->GetNode(0));      
      (*a)->SetCoeff((*a)->GetCoeff() * d); // since exponent is even, +
      (*a)->SetExponent((*a)->GetExponent() * e);
      ret = true;
    }
  } 
  return ret;
}    

// standard order for a set of subnodes of a sum is:
// constants + linear variables + vars^{rising powers} + complex operands
class NodeOrderSum {
public:
  int operator() (const Expression& a, const Expression& b) {
    if (a->IsConstant() && !b->IsConstant()) {
      return true;
    } else if (a->IsVariable() && b->IsVariable()) {
      if (a->GetExponent() == 1 && b->GetExponent() != 1) {
	return true;
      } else if (a->GetExponent() < b->GetExponent()) {
	return true;
      } else if (a->GetExponent() > b->GetExponent()) {
	return false;
      } else {
	if (a->GetVarIndex() < b->GetVarIndex()) {
	  return true;
	} else {
	  return false;
	}
      }
    } else if (a->IsLeaf() && !b->IsLeaf()) {
      return true;
    } else {
      return false;
    }
  }
};

// standard order for a set of subnodes is:
// constants + vars^{rising powers} + complex operands
class NodeOrder {
public:
  int operator() (const Expression& a, const Expression& b) {
    if (a->IsConstant() && !b->IsConstant()) {
      return true;
    } else if (a->IsVariable() && b->IsVariable()) {
      if (a->GetExponent() < b->GetExponent()) {
	return true;
      } else if (a->GetExponent() > b->GetExponent()) {
	return false;
      } else {
	if (a->GetVarIndex() < b->GetVarIndex()) {
	  return true;
	} else {
	  return false;
	}
      }
    } else if (a->IsLeaf() && !b->IsLeaf()) {
      return true;
    } else {
      return false;
    }
  }
};

bool ReorderNodes(Expression* a) {
  bool ret = false;
  if (!(*a)->IsLeaf() && (*a)->GetSize() > 1 &&
      ((*a)->GetOpType() == SUM || (*a)->GetOpType() == PRODUCT)) {
    int i;
    for(i = 0; i < (*a)->GetSize(); i++) {
      ReorderNodes((*a)->GetNodePtr(i));
    }
    // how do I get this sort to tell me whether it did anything or not?
    // at the moment this function returns false by definition, incorrect
    if ((*a)->GetOpType() == SUM) {
      sort(((*a)->GetNodeVectorPtr())->begin(), 
	   ((*a)->GetNodeVectorPtr())->end(), NodeOrderSum());
    } else {
      sort(((*a)->GetNodeVectorPtr())->begin(), 
	   ((*a)->GetNodeVectorPtr())->end(), NodeOrder());
    }
  }
  return ret;
}

bool CompactLinearPart(Expression* a) {
  bool ret = false;
  bool ischanged = false;
  ischanged = SimplifyConstant(a);
  if (!ret && ischanged)
    ret = true;
  ischanged = DifferenceToSum(a);
  if (!ret && ischanged)
    ret = true;
  ischanged = CompactLinearPartRecursive(a);
  if (!ret && ischanged) 
    ret = true;
  ischanged = ReorderNodes(a);
  return ret;
}

bool CompactLinearPartRecursive(Expression* a) {
  bool ret = false;
  bool ischanged = false;
  int i, j;
  i = 0;
  int sz = (*a)->GetSize();
  if ((*a)->GetOpType() == SUM) {
    while(i < sz) {
      ischanged = CompactLinearPartRecursive((*a)->GetNodePtr(i));
      if (!ret && ischanged)
	ret = true;
      if ((*a)->GetNode(i)->GetOpType() == SUM) {
	ret = true;
	for(j = 0; j < (*a)->GetNode(i)->GetSize(); j++) {
	  (*a)->AddNode((*a)->GetNode(i)->GetNode(j));
	}	
	(*a)->DeleteNode(i);
	sz--; // we have deleted a node
	// we don't need to increase i, since we have deleted the i-th node
	// the next node is still the i-th node
	if (sz == 1) {
	  a->SetTo((*a)->GetNode(0));
	  i = 0;
	  sz = (*a)->GetSize();
	}
      } else {
	i++; 
      }
    }
  }
  return ret;
}


// deals with cases like *(*(x,y), z) --> *(x,y,z)
bool CompactProducts(Expression* a) {
  bool ret = false;
  bool ischanged = false;
  int i, j;
  if ((*a)->GetOpType() == PRODUCT) {
    // if a node is a product, recurse until no more
    for(i = 0; i < (*a)->GetSize(); i++) {
      ischanged = CompactProducts((*a)->GetNodePtr(i));
      if (!ret && ischanged) {
	ret = true;
      }
      if ((*a)->GetNode(i)->GetOpType() == PRODUCT) {
	ret = true;
	for(j = 0; j < (*a)->GetNode(i)->GetSize(); j++) {
	  (*a)->AddNode((*a)->GetNode(i)->GetNode(j));
	}
	(*a)->DeleteNode(i);
      }
    } 
    if ((*a)->GetSize() == 1) {
      // product with just one operand, up a level
      a->SetTo((*a)->GetNode(0));
      ret = true;
    }
  } else {
    // not a product, recurse
    for(i = 0; i < (*a)->GetSize(); i++) {
      ischanged = CompactProducts((*a)->GetNodePtr(i));
      if (!ret && ischanged) {
	ret = true;
      }
    }
  }
  (*a)->ConsolidateProductCoeffs();
  return ret;
}

// generic simplification on a copy of the expression
Expression SimplifyCopy(Expression* a, bool& ischanged) {
  Expression b;
  b = (*a).Copy();
  ischanged = Simplify(&b);
  return b;
}

// recursive destroy - delete all the tree and nodes in a tree - use
// with caution
void RecursiveDestroy(Expression* a) {
  int i;
  for(i = 0; i < (*a)->GetSize(); i++) {
    RecursiveDestroy((*a)->GetNodePtr(i));
  }
  a->Destroy();
}
