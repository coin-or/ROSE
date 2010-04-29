/**********************************************************************
 * Name:       nl2ros.cxx                                           
 * Authors:    Stefano Galli / Leo Liberti                            
 * Purpose:    part of ROSE - AMPL .nl to .ROS translator (external) 
 * Source:     GNU C++                                                
 * License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
 * History:    040328 v.0.0 derived from Stefano Galli's project      
 **********************************************************************/

#include <iostream>
#include <cstdio>
#include <vector>
#include <map>

#include "asl.h"
#include "nlp.h"
#include "getstub.h"
#include "r_opn.hd" /* for N_OPS */
#include "opcode.hd"

#define CHR (char*)  
      /* for suppressing "String literal to char*" warnings */

#define R_OPS ((ASL_fg*)asl)->I.r_ops_
#define OBJ_DE ((ASL_fg*)asl)->I.obj_de_
#define CON_DE ((ASL_fg*)asl)->I.con_de_

#define EPSILONTOLERANCE 1e-30
#define ROSEINFINITY 1e30

fint timing = 0;

static keyword keywds[] = { /* must be alphabetical */
   KW(CHR"timing", L_val, &timing, CHR"display timings for the run"),
   };

static Option_Info Oinfo = { CHR"testampl", CHR"ANALYSIS TEST",
   CHR"concert_options", keywds, nkeywds, 0, CHR"ANALYSIS TEST" };

bool is_expr_zero(expr* e) {
  efunc* op = e->op;
  int opnum = Intcast op;
  if (opnum == OPNUM && fabs(e->dL) < EPSILONTOLERANCE) {
    return true;
  } else {
    return false;
  }
}

void display_expr(expr *e) {

  using namespace std;

  efunc *op;
  expr **ep;
  int opnum;
  float number;
  float  y, z, k, ind;
  char buffer[50];
  char * pEnd;
  double t;
  
   op = e->op;
   opnum = Intcast op;

  /* Printf ("op %d  optype %d  ", opnum, optype[opnum]);*/

   switch(opnum) {

      case OPPLUS:
	 printf("(");
         display_expr (e->L.e);
	 printf(" + ");
         display_expr (e->R.e);
	 printf(")");
         return;

      case OPMINUS:
	 printf("(");
         display_expr (e->L.e);
	 printf(" - ");
         display_expr (e->R.e);
	 printf(")");
         return;

      case OPMULT:
	 printf("(");
         display_expr (e->L.e);
	 printf(" * ");
         display_expr (e->R.e);
	 printf(")");
         return;

      case OPDIV:
         display_expr (e->L.e);
	 printf(" / ");
         display_expr (e->R.e);
         return;

      case OPREM:
         Printf ("remainder\n");
         display_expr (e->L.e);
         display_expr (e->R.e);
         return;

      case OPPOW:
	 printf("(");
         display_expr (e->L.e);
	 printf(")^(");
         display_expr (e->R.e);
	 printf(")");
         return;

      case OPLESS:
         Printf ("less\n");
         display_expr (e->L.e);
         display_expr (e->R.e);
         return;

      case MINLIST:
         Printf ("min -- not implemented\n");
         exit(1);

      case MAXLIST:
         Printf ("max -- not implemented\n");
         exit(1);

      case FLOOR:
         Printf ("floor -- not implemented\n");
         exit(1);

      case CEIL:
         Printf ("ceil -- not implemented\n");
         exit(1);

      case ABS:
         Printf ("abs");
         display_expr (e->L.e);
         return;

      case OPUMINUS:
         Printf ("(-(");
         display_expr (e->L.e);
	 printf("))");
         return;

      case OPIFnl:
         Printf ("if nl -- not implemented\n");
         exit(1);

      case OP_tanh:
         Printf ("tanh");
         display_expr (e->L.e);
         return;

      case OP_tan:
         Printf ("tan");
         display_expr (e->L.e);
         return;

      case OP_sqrt:
         Printf ("sqrt(");
         display_expr (e->L.e);
	 printf(")");
         return;

      case OP_sinh:
         Printf ("sinh");
         display_expr (e->L.e);
         return;

      case OP_sin:
         Printf ("sin");
         display_expr (e->L.e);
         return;

      case OP_log10:
         Printf ("log10");
         display_expr (e->L.e);
         return;

      case OP_log:
         Printf ("log(");
         display_expr (e->L.e);
	 printf(")");
         return;

      case OP_exp:
         Printf ("exp");
         display_expr (e->L.e);
         return;

      case OP_cosh:
         Printf ("cosh");
         display_expr (e->L.e);
         return;

      case OP_cos:
         Printf ("cos");
         display_expr (e->L.e);
         return;

      case OP_atanh:
         Printf ("atanh");
         display_expr (e->L.e);
         return;

      case OP_atan2:
         Printf ("atan2");
         display_expr (e->L.e);
         return;

      case OP_atan:
         Printf ("atan");
         display_expr (e->L.e);
         return;

      case OP_asinh:
         Printf ("asinh");
         display_expr (e->L.e);
         return;

      case OP_asin:
         Printf ("asin");
         display_expr (e->L.e);
         return;

      case OP_acosh:
         Printf ("acosh");
         display_expr (e->L.e);
         return;

      case OP_acos:
         Printf ("acos");
         display_expr (e->L.e);
         return;

      case OPSUMLIST:
	 printf("(");
         for (ep = e->L.ep; ep < e->R.ep; *ep++) {
	   display_expr (*ep);
	   if (ep < e->R.ep - 1) {
	     printf(" + ");
	   }
	 }
	 printf(")");
         return;

      case OPintDIV:
         Printf ("int division");
         display_expr (e->L.e);
         display_expr (e->R.e);

      case OPprecision:
         Printf ("precision");
         display_expr (e->L.e);
         return;

      case OPround:
         Printf ("round");
         display_expr (e->L.e);
         return;

      case OPtrunc:
         Printf ("trunc");
         display_expr (e->L.e);
         return;

      case OP1POW:
	 cout << "(";
         display_expr (e->L.e);
	 t = e->R.en->v;
	 if (t >= 0) {
	   cout << ")^" << e->R.en->v;
	 } else {
	   cout << ")^(" << e->R.en->v << ")";
	 }
         return;
	 
      case OP2POW:
         cout << "(";
	 display_expr (e->L.e);
         cout << ")^2";
         return;

      case OPCPOW:
	 t = e->L.en->v;
	 if (t >= 0) {
	   cout << e->L.en->v << "^(";
	 } else {
	   cout << "(" << e->L.en->v << ")^(";
	 }
         display_expr (e->R.e);
	 cout << ")";
         return;

      case OPFUNCALL:
         Printf ("function call -- not implemented\n");
         exit(1);

      case OPNUM:
	t = ((expr_n*)e)->v;
	if (t >= 0) {
	  cout << t;
	} else {
	  cout << "(" << t << ")";
	}
	return;

      case OPPLTERM:
         Printf ("pl term -- not implemented\n");
         exit(1);

      case OPIFSYM:
         Printf ("if sym -- not implemented\n");
         exit(1);

      case OPHOL:
         Printf ("string argument -- not implemented\n");
         exit(1);

      case OPVARVAL:
	 // Printf ("(x%d)", e->a + 1);
         Printf ("x%d", e->a + 1);
         return;

      default:
         Printf ("other -- not implemented\n");
         exit(1);
   }
}

/*----------------------------------------------------------------------
  Main Program
----------------------------------------------------------------------*/

int main(int argc, char **argv) {

  using namespace std;

   FILE *nl;
   ASL *asl;

   int i, j, k;
   char *stub, *getenvErr;

   cgrad *cg;
   ograd *og;

   efunc *r_ops_int[N_OPS];

   /*** Get name of .nl file; read problem sizes ***/

   if (argc < 2) {
     cout << "nl2ros --- Stefano Galli / Leo Liberti\n";
     cout << "  This utility translates AMPL .nl stub files into .ros files\n";
     cout << "  that can be read by Rose. In order to obtain an .nl file from\n";
     cout << "  an AMPL .mod/.dat problem pair, call AMPL with the option '-ogproblem'.\n";
     cout << "  This will produce a file 'problem.nl'. Then run nl2ros as follows:\n";
     cout << "    " << argv[0] << " problem.nl > problem.ros\n";
     cout << "  to produce the .ros file\n";
     return 1;
   }

   asl = ASL_alloc(ASL_read_fg);
   stub = getstub(&argv, &Oinfo);
   nl = jac0dim(stub, (fint)strlen(stub));

   /*** Read coefficients, bounds, and op tree from .nl file ***/

   Uvx = (real *)Malloc(n_var*sizeof(real));
   Urhsx = (real *)Malloc(n_con*sizeof(real));
   A_vals = (real *)Malloc(nzc*sizeof(real));

   
   for(i = 0; i < N_OPS; i++)
      r_ops_int[i] = (efunc*)(unsigned long)i;
   R_OPS = r_ops_int;
   want_derivs = 0;
   fg_read(nl,0);
   R_OPS = 0;

   /*** Get and record options ***/

   if (getopts(argv, &Oinfo)) exit(1);

   int n_badvals = 0;

   switch(timing) {
      case 0:
         break;
      case 1:
         break;
      default:
         cout << "Invalid value " << timing 
            << " for directive timing" << endl;
         n_badvals++;
         }

   if (n_badvals)
      exit(1);

   /*-------------------------------------------------------------------
     For a start, write out some statistics (using printf or <<)
   -------------------------------------------------------------------*/

   cout << "# .ros file automatically translated by nl2ros on .nl input\n";

   // integrality of variables
   int thenlv = max(nlvc, nlvo);
   int theother = n_var - thenlv - nwv - nbv - niv;
   vector<bool> integrality(n_var, false);
   i = 0;
   for(; i < nlvb - nlvbi; i++) {
     integrality[i] = false;
   }
   for(; i < nlvb; i++) {
     integrality[i] = true;
   }
   for(; i < nlvc - nlvci; i++) {
     integrality[i] = false;
   }
   for(; i < nlvc; i++) {
     integrality[i] = true;
   }
   for(; i < nlvo - nlvoi; i++) {
     integrality[i] = false;
   } 
   for(; i < nlvo; i++) {
     integrality[i] = true;
   }
   for(; i < thenlv + nwv; i++) {
     integrality[i] = false;
   }
   for(; i < thenlv + nwv + theother; i++) {
     integrality[i] = false;
   }
   for(; i < n_var - niv; i++) {
     integrality[i] = true;
     LUv[i] = 0;
     Uvx[i] = 1;
   }
   for(; i < n_var; i++) {
     integrality[i] = true;
   }

#ifdef NL2ROSDEBUG
   cout << "# n_var = " << n_var << endl;
   cout << "# thenlv = " << thenlv << endl;
   cout << "# nlvb = " << nlvb << endl;
   cout << "# nlvbi = " << nlvbi << endl;
   cout << "# nlvc = " << nlvc << endl;
   cout << "# nlvci = " << nlvci << endl;
   cout << "# nlvo = " << nlvo << endl;
   cout << "# nlvoi = " << nlvoi << endl;
   cout << "# nwv = " << nwv << endl;
   cout << "# theother = " << theother << endl;
   cout << "# nbv = " << nbv << endl;
   cout << "# niv = " << niv << endl;
#endif

   // variables section
   cout << "variables =\n";
   for(i = 0; i < n_var; i++) {
     if (LUv[i] < -ROSEINFINITY) {
       cout << "MinusInfinity < ";
     } else {
       cout << LUv[i] << " < ";
     }
     cout << "x" << i + 1;
     if (Uvx[i] > ROSEINFINITY) {
       cout << " < PlusInfinity";
     } else {
       cout << " < " << Uvx[i];
     }
     if (integrality[i]) {
       cout << " / Integer";
     } else {
       cout << " / Continuous";
     }
     if (i < n_var - 1) {
       cout << ",\n";
     } else {
       cout << ";\n";
     }
   }
     
   // objective function
   if (n_obj > 1) {
     cout << "nl2ros: cannot deal with more than 1 objective function\n";
     exit(2);
   }
   bool firstflag = false;
   for (i = 0; i < n_obj; i++) {
     firstflag = true;
     printf("objfun = [");
     for(og = Ograd[i]; og; og = og->next) {
       if (fabs(og->coef) < EPSILONTOLERANCE) {
	 og->coef = 0;
       }
       if (og->coef >= 0) {
	 if (og->coef == 1) {
	   if (!firstflag) {
	     cout << " + ";
	   }
	   cout << "x" << og->varno + 1;
	   firstflag = false;
	 } else if (og->coef != 0) {
	   if (!firstflag) {
	     cout << " + ";
	   }
	   cout << og->coef << "*x" << og->varno + 1;
	   firstflag = false;
	 }
       } else {
	 if (!firstflag) {
	   cout << " + ";
	 }
	 cout << "(" << og->coef << ")*x" << og->varno + 1;
	 firstflag = false;
       }
     }

     /*
     ASL_fg* theaslfg = (ASL_fg*) asl;
     cout << (int) *(theaslfg->i.objtype_) << endl; // if 0, min, if 1, max
     */
     if (!is_expr_zero((OBJ_DE + i)->e)) {
       if (!firstflag) {
	 cout << " + ";
       }
       display_expr((OBJ_DE + i)->e);
       firstflag = false;
     }
     printf("];\n");
   }
	
   // constraints
   if (n_con > 0) {
     // read linear jacobian in sparse form
     map<pair<int,int>,double> A;
     map<pair<int,int>,double>::iterator mit;
     pair<int,int> p;
     for(j = 0; j < n_var; j++) {
       p.second = j;
       for(i = A_colstarts[j]; i < A_colstarts[j + 1]; i++) {
	 p.first = A_rownos[i];
	 if (fabs(A_vals[i]) < EPSILONTOLERANCE) {
	   A_vals[i] = 0;
	 }
	 A[p] = A_vals[i];
       }
     }
     // print constraints
     cout << "constraints = \n";
     for (i = 0; i < n_con; i++) {
       firstflag = true;
       cout << "[ ";
       if (LUrhs[i] < -ROSEINFINITY) {
	 cout << "MinusInfinity";
       } else {
	 cout << LUrhs[i];
       } 
       cout << " < ";
       p.first = i;
       for(j = 0; j < n_var; j++) {
	 p.second = j;
	 mit = A.find(p);
	 if (mit != A.end()) {
	   if (mit->second >= 0) {
	     if (mit->second == 1) {
	       if (!firstflag) {
		 cout << " + ";
	       }
	       cout << "x" << j + 1;
	       firstflag = false;
	     } else if (mit->second != 0) {	     
	       if (!firstflag) {
		 cout << " + ";
	       }
	       cout << mit->second << "*x" << j + 1;
	       firstflag = false;
	     }
	   } else {
	     if (!firstflag) {
	       cout << " + ";
	     }
	     cout << "(" << mit->second << ")*x" << j + 1;
	     firstflag = false;
	   }
	 }
       }
       if (!is_expr_zero((CON_DE + i)->e)) {
	 if (!firstflag) {
	   cout << " + ";
	 }
	 display_expr((CON_DE + i)->e);
	 firstflag = false;
       }
       cout << " < ";
       if (Urhsx[i] > ROSEINFINITY) {
	 cout << "PlusInfinity";
       } else {
	 cout << Urhsx[i];
       }
       cout << " ]";
       if (i == n_con - 1) {
	 cout << ";\n";
       } else {
	 cout << ",\n";
       }
     }
   } else {
     cout << "constraints = [ ];" << endl;
   }

   // starting point
   cout << "startingpoint = ";
   for(i = 0; i < n_var; i++) {
     cout << "0";
     if (i < n_var - 1) {
       cout << ",";
     } else {
       cout << ";\n";
     }
   }

   // options
   cout << "options = MainSolver sBB;\n";
   return 0;
}
