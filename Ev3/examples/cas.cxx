/**********************************************************************
* Author:      Leo Liberti                                            *
* Name:        test.cxx                                               *
* Source:      GNU C++                                                *
* Purpose:     main driver program for Expression v3                  *
* History:     021218 0.0 work started                                *
* License:    (C) Leo Liberti, all rights reserved. Code published under the
               Common Public License.
***********************************************************************/

#include "expression.h"
#include "parser.h"

#include <stdlib.h>
#include <iostream>
#include <map>
#include <string>
#include <cstring>

using namespace std;

#define BUFSIZE 8192
#define MAXVARS 1000

enum { EVALUATE, DIFFERENTIATE, SIMPLIFY, ASSIGN, UNKNOWN };

int main(int argc, char** argv) {

  char buffer[BUFSIZE];
  char bufsave[BUFSIZE];
  char *tok;
  char *nexttok;
  int cmd;

  double varvalues[MAXVARS];

  ExpressionParser p;
  Expression e, f, d;
  map<string, Expression> store;
  int nerr;

  while(true) {
    cout << "Ev3> ";
    cin >> buffer;
    strncpy(bufsave, buffer, BUFSIZE);
    tok = buffer;
    // parse the command
    if (strncmp(buffer, "quit", 4) == 0) {
      break;
    } else if (strncmp(tok, "eval", 4) == 0) {
      cmd = EVALUATE;
      nexttok = strchr(tok, '(');
      if (!nexttok) {
	cerr << "syntax error: was expecting a (\n";
	continue;
      } else {
	tok = nexttok + 1;
	nexttok = strchr(tok, ',');
	if (!nexttok) {
	  cerr << "syntax error: was expecting a ,\n";
	  continue;
	} else {
	  // delimit tok
	  *nexttok = '\0';
	  string t(tok);
	  if (store.find(t) != store.end()) {
	    // string t (i.e. tok) was in the store
	    e = store[t];
	  } else {
	    // not in the store, parse it in
	    e = p.Parse(tok, nerr);
	  }
	  int vv = e->NumberOfVariables();
	  // "undelimit"
	  *nexttok = ',';
	  tok = nexttok + 1;
	  // read in variable values
	  int i = 0;
	  bool errflag = false;
	  while(nexttok) {
	    if (i >= vv) {
	      break;
	    }
	    tok = nexttok + 1;
	    if (i == vv - 1) {
	      nexttok = strchr(tok, ')');
	    } else {
	      nexttok = strchr(tok, ',');
	    }
	    if (!nexttok) {
	      cerr << "can't evaluate, not enough values: " << i + 1
		   << " out of " << vv << " needed\n";
	      errflag = true;
	      break;
	    }
	    *nexttok = '\0';
	    varvalues[i] = atof(tok);
	    i++;
	    *nexttok = ',';
	  }
	  if (errflag) {
	    continue;
	  }
	  cout << "  " << e->Eval(varvalues, vv) << endl;
	}
      }
    } else if (strncmp(tok, "diff", 4) == 0) {
      cmd = DIFFERENTIATE;
    } else if (strncmp(tok, "simp", 4) == 0) {
      cmd = SIMPLIFY;
    } else if (strncmp(tok, "let", 3) == 0) {
      cmd = ASSIGN;
    } else {
      cmd = UNKNOWN;
    }

    /*
    bool ischanged;
    cout << "D(" << e->ToString();
    f = SimplifyCopy(&e, ischanged);
    cout << " --> " << f->ToString();
    d = Diff(f, 1);
    cout << ", var1) = " << d->ToString() << endl;
    */
  }

  return 0;
}
