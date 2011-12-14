/*
** Name:       solver_Rconvexifiermod.cxx
** Author:     Leo Liberti & Sonia Cafieri
** Source:     GNU C++
** Purpose:    Solver Rconvexifiermod implementation
** License:    (C) Leo Liberti, all rights reserved. Code published under the
               Common Public License.
** History:    081122 work started
*/

#include <cstdio>
#include <cstring>
#include <cassert>
#include <sys/time.h>
#include <fstream>
#include <sstream>
#include "newsolver.h"
#include "solver.h"
#include "solver_Rconvexifiermod.h"
#include "utils.h"
#include "Ev3/auxiliary.h"
#include <iomanip>

#define FIRSTOBJ 1
#define EPSILON 1e-4
#define WDEFNAMEBUFSIZE 256
#define SMALL 1e-12

// RCONVEXIFIERMODSOLVER CLASS METHODS

RconvexifiermodSolver::RconvexifiermodSolver() {
  TheName = "Rconvexifiermod";
  NumberOfVariables = 0;
  NumberOfConstraints = 0;
  NumberOfObjectives = 0;
  IsSolved = false;
  OptDir = Minimization;
  OptDirCoeff = 1;
  ManualOptDir = false;
  OptimalObjVal = ROSEINFINITY;
  Epsilon = EPSILON;
  LocalSolverName = "Rsmith";
  MaxRunningTime = 0;
  OutFile = "Rconvexifiermod_out.mod";
  OutFile2 = "Rconvexifiermod_out.dat";
}

RconvexifiermodSolver::~RconvexifiermodSolver() {
  if (LocalSolver) {
    DeleteSolver(LocalSolver);
  }
}

void RconvexifiermodSolver::SetProblem(Problem* p) {
  InProb = p;
  InitProb = false;
}

bool RconvexifiermodSolver::CanSolve(int problemtype) {
  // fill in here with information about solver capabilities
}

void RconvexifiermodSolver::Initialize(bool force = false) {

  if (InitProb) {
    ;
  } else {
    InitProb = true;
  }

  // check the input problem is set to something
  if (InProb == NULL) {
    cerr << "RconvexifiermodSolver::Initialize(): error: InProb is NULL\n";
    assert(InProb != NULL);
  }

  NumberOfVariables = InProb->GetNumberOfVariables();
  NumberOfConstraints = InProb->GetNumberOfConstraints();
  NumberOfObjectives = InProb->GetNumberOfObjectives();

  // parameters
  ReplaceParams(InProb->GetParamsRef());
  int np = ParameterBlob.GetNumberOfParameters();
  for(int i = 1; i <= np; i++) {
    if (ParameterBlob.GetParameterName(i) == "RconvexifiermodLocalSolver") {
      LocalSolverName = ParameterBlob.GetParameterStringValue(i);
    } else if (ParameterBlob.GetParameterName(i) == "ConvexifiermodMaxTime") {
      MaxRunningTime = (double) ParameterBlob.GetParameterIntValue(i);
      assert(MaxRunningTime >= 0);
    } else if (ParameterBlob.GetParameterName(i) == "ConvexifiermodEpsilon") {
      Epsilon = ParameterBlob.GetParameterDoubleValue(i);
      assert(Epsilon > 0);
    } else if (ParameterBlob.GetParameterName(i) == "ConvexifiermodQuiet" ||
               ParameterBlob.GetParameterName(i) == "Quiet") {
      Quiet = ParameterBlob.GetParameterBoolValue(i);
    } else if (ParameterBlob.GetParameterName(i) == "MainSolver") {
      MainSolver = ParameterBlob.GetParameterStringValue(i);
    } else if (ParameterBlob.GetParameterName(i) == "AmplFlag") {
      AmplFlag = ParameterBlob.GetParameterBoolValue(i);
    }
  }

  if (MainSolver != TheName && AmplFlag) {
    Quiet = true;
  }

  // obj fun dir
  OptDir = InProb->GetObjectiveLI(1)->OptDir;

  // variable integrality
  for(int i = 1; i <= NumberOfVariables; i++) {
    integrality.push_back(InProb->GetVariableLI(i)->IsIntegral);
  }

  f = ROSEINFINITY;
  fstar = ROSEINFINITY;

  // create and configure local solver
  LocalSolver = NewSolver(LocalSolverName);
  LocalSolver->SetProblem(InProb);
  LocalSolver->ReplaceParams(ParameterBlob);

  f = ROSEINFINITY;
  fstar = ROSEINFINITY;
}

string RconvexifiermodSolver::VarIndexToString(int vi) const {
  stringstream outbuf;
  if (vi <= nvar_orig) {
    if (integrality[vi-1]) {
      outbuf << "Rcvx_xInt[" << vi << "]";
    } else {
      outbuf << "Rcvx_x[" << vi << "]";
    }
  } else {
    outbuf << "Rcvx_w[" << vi << "]";
  }
  return outbuf.str();
}

string RconvexifiermodSolver::VarIndexToString(double cf, int vi) const {
  stringstream outbuf;
  if (fabs(cf - 1) > SMALL) {
    outbuf << "(" << cf << ")*";
  }
  if (vi <= nvar_orig) {
    if (integrality[vi-1]) {
      outbuf << "Rcvx_xInt[" << vi << "]";
    } else {
      outbuf << "Rcvx_x[" << vi << "]";
    }
  } else {
    outbuf << "Rcvx_w[" << vi << "]";
  }

  return outbuf.str();
}

int RconvexifiermodSolver::Solve(bool reinitialize) {

  int ret = 0;
  if (reinitialize || !InitProb) {
    Initialize();
  }

  // do it
  nvar_orig = NumberOfVariables;
  nconstr_orig = NumberOfConstraints;

  // Standardize the problem using Rsmith
  //LocalSolver->Initialize(false);
  LocalSolver->Solve();

  NumberOfVariables = InProb->GetNumberOfVariables();
  NumberOfConstraints = InProb->GetNumberOfConstraints();

  map<int,vector<int> > addedvar2constrLI;
  int addedvarcounter = nvar_orig;
  int addedconcounter = nconstr_orig + 1;
  int stdconstr = 0;

  // re-set variable data
  for(int i = 0; i < NumberOfVariables; i++) {
    vlb.push_back(InProb->GetVariableLI(i + 1)->LB) ;
    vub.push_back(InProb->GetVariableLI(i + 1)->UB) ;
    x.push_back(0);
    xstar.push_back(0);
    StartingPoint.push_back(x[i]);
  }

  vector<Expression> defcons;
  string cn = "c";

  int frac = 0;
  int v1f, v2f;

  ofstream out(OutFile.c_str());

  out << "# convexified .mod file produced by ROSE (Rconvexifiermod solver)"
      << endl;
  out << "# cardinalities" << endl;
  out << "param Rcvx_NumberOfVars := " << nvar_orig << ";\n";
  out << "param Rcvx_NumberOfAddedVars := " << NumberOfVariables - nvar_orig
      << ";\n" << endl;
  out << "# sets of variable indices" << endl;
  out << "set Rcvx_Vars := 1 .. Rcvx_NumberOfVars;" << endl;
  out << "set Rcvx_AddedVars := \n"
      << "  Rcvx_NumberOfVars+1 .. Rcvx_NumberOfVars+Rcvx_NumberOfAddedVars;\n";
  out << "set Rcvx_AllVars := Rcvx_Vars union Rcvx_AddedVars;\n" << endl;
  out << "# lower/upper variable bounds" << endl;
  out << "param Rcvx_L{Rcvx_AllVars};" << endl;
  out << "param Rcvx_U{Rcvx_AllVars};\n" << endl;
  out << "# original variables (continuous)" << endl;
  if (InProb->IsProblemContinuous()) {
    out << "var Rcvx_x{j in Rcvx_Vars} <= Rcvx_U[j], >= Rcvx_L[j];" << endl;
  } else {
    out << "set  Rcvx_contVars;" << endl;
    out << "var Rcvx_x{j in Rcvx_contVars} <= Rcvx_U[j], >= Rcvx_L[j];" << endl;
  }
  if (!InProb->IsProblemContinuous()) {
    out << "# original variables (integer)" << endl;
    out << "set  Rcvx_intVars;" << endl;
    out << 
      "var Rcvx_xInt{j in Rcvx_intVars} <= Rcvx_U[j], >= Rcvx_L[j], integer;" 
	<< endl;

    int n_varcont = 0;
    for(int i = 1; i <= nvar_orig; i++) {
       if (!InProb->GetVariableLI(i)->IsIntegral) {
         n_varcont ++;
       }
    }
    int n_varint = nvar_orig - n_varcont;
    out << "let  Rcvx_contVars := {";
    int n_varcont1 = 0;
    for(int i = 1; i <= nvar_orig; i++) {
       if (!InProb->GetVariableLI(i)->IsIntegral) {
         n_varcont1 ++;
         out << i;
         if (n_varcont1 < n_varcont) out << ",";
       }
    }
    out << "};" << endl;
    out << "let  Rcvx_intVars := {";
    int n_varint1 = 0;
    for(int i = 1; i <= nvar_orig; i++) {
       if (InProb->GetVariableLI(i)->IsIntegral) {
         n_varint1 ++;
         out << i;
         if (n_varint1 < n_varint) out << ",";
       }
    }
    out << "};" << endl;

  }
  out << "# added variables" << endl;
  out << "var Rcvx_w{j in Rcvx_AddedVars} <= Rcvx_U[j], >= Rcvx_L[j];\n"
      << endl;
  int objvarindex = InProb->GetObjective(1)->Function->GetVarIndex();
  if (OptDir == Minimization) {
    out << "minimize ";
  } else {
    out << "maximize ";
  }
//new 090720
  int sizeobj = InProb->GetObjective(1)->Function->GetSize();
  int vobj;
  double coeffobj;
  if (sizeobj > 1) {
    // linear objective in original problem, Smith reformulator does not standardize:
    // more than 1 node in the obj expression tree is possible 
    out << "Rcvx_obj : " ;
    for(int k = 0; k < sizeobj; k++) {
      coeffobj = InProb->GetObjective(1)->Function->GetNode(k)->GetCoeff();
      vobj = InProb->GetObjective(1)->Function->GetNode(k)->GetVarIndex();
      out << VarIndexToString(coeffobj, vobj);
      if (k < sizeobj - 1) {
	  out << " + ";
      }
    }
    out << ";\n\n";
  } else {
    coeffobj = InProb->GetObjective(1)->Function->GetCoeff();
    out << "Rcvx_obj : " << coeffobj << "*" <<VarIndexToString(objvarindex) << ";\n\n";
  }

  int addedvarindex;
  int size;
  double coeff;
  int v1, v2;
  int lcc = 1;
  double l, u;
  stringstream cname;
  stringstream cnamelist;

  for(int j = 1; j <= NumberOfConstraints; j++) {

    if (InProb->GetConstraintLI(j)->Function->IsLinear()) {
      // linear constraint, just print
      out << "subject to ";
      cname << "Rcvx_lincon" << lcc;
      out << cname.str() << " : ";
      cnamelist << ", " << cname.str();
      cname.str("");
      size = InProb->GetConstraintLI(j)->Function->GetSize();
      for(int k = 0; k < size; k++) {
	coeff = InProb->GetConstraintLI(j)->Function->GetNode(k)->GetCoeff();
	v1 = InProb->GetConstraintLI(j)->Function->GetNode(k)->GetVarIndex();
	out << VarIndexToString(coeff, v1);
	if (k < size - 1) {
	  out << " + ";
	}
      }
      l = InProb->GetConstraintLI(j)->LB;
      u = InProb->GetConstraintLI(j)->UB;
      if (fabs(u-l) < SMALL) {
	// equality constraint
	out << " = " << u;
      } else if (l > -ROSEINFINITY && u >= ROSEINFINITY) {
	// >= constr
	out << " >= " << l;
      } else if (l <= -ROSEINFINITY && u < ROSEINFINITY) {
	// <= constr
	out << " <= " << u;
      } else if (l > -ROSEINFINITY && u < ROSEINFINITY) {
	// two-sided constraint
	out << " <= " << u << ";\n";
	out << "subject to ";
	cname << "Rcvx_lincon" << lcc << "_2";
	cnamelist << ", " << cname.str();
	out << cname.str() << " : ";
	cname.str("");
	for(int k = 0; k < size; k++) {
	  coeff = InProb->GetConstraintLI(j)->Function->GetNode(k)->GetCoeff();
	  v1 = InProb->GetConstraintLI(j)->Function->GetNode(k)->GetVarIndex();
	  out << VarIndexToString(coeff, v1);
	  if (k < size - 1) {
	    out << " + ";
	  }
	}
	out << " >= " << l;
      }
      out << ";" << endl;
      lcc++;
    } else {
      // nonlinear defining constraint, relax
      size = InProb->GetConstraintLI(j)->Function->GetNode(1)->GetSize();
      int op = InProb->GetConstraintLI(j)->Function->GetNode(1)->GetOpType();
      double expo =
	InProb->GetConstraintLI(j)->Function->GetNode(1)->GetExponent();
      coeff=InProb->GetConstraintLI(j)->Function->GetNode(1)->GetCoeff();
      Expression e = - InProb->GetConstraintLI(j)->Function->GetLinearPart();
      addedvarindex = InProb->GetConstraintLI(j)->Function->
	GetNode(0)->GetVarIndex();

      if(op == FRACTION) {
	// find variables involved in fraction
	vector<int> vidx;
	Expression xs(1, 1, "x"), ys(1, 2, "y");
	Expression schema = xs/ys;
	Constraint* theCon = InProb->GetConstraintLI(j);
	if (!theCon->Function->IsLinear()) {
	  theCon->Function->GetVarIndicesInSchema(vidx, schema);
	}
	// TO DO
	cerr << "Rconvexifiermod: error: fractional relaxation not implemented"
	     << endl;
	exit(132);
      }

      if(op == PRODUCT && size == 2) {

	// find variables involved in bilinear terms
	vector<int> vidx;
	Expression xs(1, 1, "x"), ys(1, 2, "y");
	Expression schema = xs*ys;
	Constraint* theCon = InProb->GetConstraintLI(j);
	if (!theCon->Function->IsLinear()) {
	  theCon->Function->GetVarIndicesInSchema(vidx, schema);
	}
	v1 = vidx[0];
	v2 = vidx[1];

	double cf = 1;
	double coeff1 = InProb->
	  GetConstraintLI(j)->Function->GetNode(1)->GetNode(0)->GetCoeff();
	if (coeff != 1) {
	  cf = cf*coeff;    // case  coeff*(xx*yy)
	}
	if (coeff1 != 1) {
	  cf = cf*coeff1; // case  ((coeff1*xx)*yy)
	}

	out << "subject to ";
	cname << "Rcvx_prod_" << v1 << "_" << v2 << "_a";
	cnamelist << ", " << cname.str();
	out << cname.str() << " : "
	    << VarIndexToString(1/cf, addedvarindex) << " >= "
	    << "Rcvx_L[" << v1 << "]*" << VarIndexToString(v2) << " + "
	    << "Rcvx_L[" << v2 << "]*" << VarIndexToString(v1) << " - "
	    << "Rcvx_L[" << v1 << "]*Rcvx_L[" << v2 << "];" << endl;
	cname.str("");
	out << "subject to ";
	cname << "Rcvx_prod_" << v1 << "_" << v2 << "_b";
	cnamelist << ", " << cname.str();
	out << cname.str() << " : "
	    << VarIndexToString(1/cf, addedvarindex) << " >= "
	    << "Rcvx_U[" << v1 << "]*" << VarIndexToString(v2) << " + "
	    << "Rcvx_U[" << v2 << "]*" << VarIndexToString(v1) << " - "
	    << "Rcvx_U[" << v1 << "]*Rcvx_U[" << v2 << "];" << endl;
	cname.str("");
	out << "subject to ";
	cname << "Rcvx_prod_" << v1 << "_" << v2 << "_c";
	cnamelist << ", " << cname.str();
	out << cname.str() << " : "
	    << VarIndexToString(1/cf, addedvarindex) << " <= "
	    << "Rcvx_L[" << v1 << "]*" << VarIndexToString(v2) << " + "
	    << "Rcvx_U[" << v2 << "]*" << VarIndexToString(v1) << " - "
	    << "Rcvx_L[" << v1 << "]*Rcvx_U[" << v2 << "];" << endl;
	cname.str("");
	out << "subject to ";
	cname << "Rcvx_prod_" << v1 << "_" << v2 << "_d";
	cnamelist << ", " << cname.str();
	out << cname.str() << " : "
	    << VarIndexToString(1/cf, addedvarindex) << " <= "
	    << "Rcvx_U[" << v1 << "]*" << VarIndexToString(v2) << " + "
	    << "Rcvx_L[" << v2 << "]*" << VarIndexToString(v1) << " - "
	    << "Rcvx_U[" << v1 << "]*Rcvx_L[" << v2 << "];" << endl;
	cname.str("");

      } else if(op == POWER || (op == VAR && (expo!=0 && expo!=1))) {

	// find variable involved in power
	vector<int> vidx;
	Expression xs(1, 1, "x");
	Expression schema = xs^expo;
	Constraint* theCon = InProb->GetConstraintLI(j);
	if (!theCon->Function->IsLinear()) {
	  theCon->Function->GetVarIndicesInSchema(vidx, schema);
	}

	int v1 = vidx[0];

	// convex powers
	if(int(expo)%2==0 ||
	   (int(expo)%2==1 && vlb[v1-1]>=0 && vub[v1-1]>=0) ) {

	  out << "subject to ";
	  cname << "Rcvx_cvxpow_" << v1 << "_" << expo << "_secant";
	  cnamelist << ", " << cname.str();
	  out << cname.str() << " : "
	      << VarIndexToString(1/coeff, addedvarindex) << " <= "
	      << "Rcvx_L[" << v1 << "]^" << expo << " + "
	      << "((Rcvx_U[" << v1 << "]^" << expo << " - "
	      << "Rcvx_L[" << v1 << "]^" << expo << ") / ("
	      << "Rcvx_U[" << v1 << "] - Rcvx_L[" << v1 << "]))*("
	      << VarIndexToString(v1) << " - Rcvx_L[" << v1 << "]);" << endl;
	  cname.str("");

	  out << "subject to ";
	  cname << "Rcvx_cvxpow_" << v1 << "_" << expo << "_tanlow";
	  cnamelist << ", " << cname.str();
	  out << cname.str() << " : "
	      << VarIndexToString(1/coeff, addedvarindex) << " >= "
	      << "(" << expo << "*Rcvx_L[" << v1 << "]";
	  if (expo > 2) {
	    out << "^" << expo - 1;
	  }
	  out << ")*(" << VarIndexToString(v1) << " - Rcvx_L[" << v1 << "])"
	      << " + Rcvx_L[" << v1 << "]^" << expo << ";" << endl;
	  cname.str("");

	  out << "subject to ";
	  cname << "Rcvx_cvxpow_" << v1 << "_" << expo << "_tanupp";
	  cnamelist << ", " << cname.str();
	  out << cname.str() << " : "
	      << VarIndexToString(1/coeff, addedvarindex) << " >= "
	      << "(" << expo << "*Rcvx_U[" << v1 << "]";
	  if (expo > 2) {
	    out << "^" << expo - 1;
	  }
	  out << ")*(" << VarIndexToString(v1) << " - Rcvx_U[" << v1 << "])"
	      << " + Rcvx_U[" << v1 << "]^" << expo << ";" << endl;
	  cname.str("");

	} else if (int(expo)%2==1 && vlb[v1-1]<=0 && vub[v1-1]<=0) {
	  // concave powers
	  out << "subject to ";
	  cname << "Rcvx_cncpow_" << v1 << "_" << expo << "_secant";
	  cnamelist << ", " << cname.str();
	  out << cname.str() << " : "
	      << VarIndexToString(1/coeff, addedvarindex) << " >= "
	      << "Rcvx_L[" << v1 << "]^" << expo << " + "
	      << "((Rcvx_U[" << v1 << "]^" << expo << " - "
	      << "Rcvx_L[" << v1 << "]^" << expo << ") / ("
	      << "Rcvx_U[" << v1 << "] - Rcvx_L[" << v1 << "]))*("
	      << VarIndexToString(v1) << " - Rcvx_L[" << v1 << "]);" << endl;
	  cname.str("");

	  out << "subject to ";
	  cname << "Rcvx_cncpow_" << v1 << "_" << expo << "_tanlow";
	  cnamelist << ", " << cname.str();
	  out << cname.str() << " : "
	      << VarIndexToString(1/coeff, addedvarindex) << " <= "
	      << "(" << expo << "*Rcvx_L[" << v1 << "]";
	  if (expo > 2) {
	    out << "^" << expo - 1;
	  }
	  out << ")*(" << VarIndexToString(v1) << " - Rcvx_L[" << v1 << "])"
	      << " + Rcvx_L[" << v1 << "]^" << expo << ";" << endl;
	  cname.str("");

	  out << "subject to ";
	  cname << "Rcvx_cncpow_" << v1 << "_" << expo << "_tanupp";
	  cnamelist << ", " << cname.str();
	  out << cname.str() << " : "
	      << VarIndexToString(1/coeff, addedvarindex) << " <= "
	      << "(" << expo << "*Rcvx_U[" << v1 << "]";
	  if (expo > 2) {
	    out << "^" << expo - 1;
	  }
	  out << ")*(" << VarIndexToString(v1) << " - Rcvx_U[" << v1 << "])"
	      << " + Rcvx_U[" << v1 << "]^" << expo << ";" << endl;
	  cname.str("");

	} else if(int(expo)%2==1 && vlb[v1-1]<=0 && vub[v1-1] >= 0) {
	  // x^(2k+1) with range including 0

	  int k = int((expo-1)/2);
	  if (k > 10) {
	    cerr << "RconvexifiermodSolver::Solve(): error: "
		 << "exponent greater than 21 - case not implemented \n";
	    exit(0);
	  }
	  double rk[10];
	  rk[0] = -0.5; rk[1] = -0.6058295862; rk[2] = -0.6703320476;
	  rk[3] = -0.7145377272; rk[4] = -0.7470540749; rk[5] = -0.7721416355;
	  rk[6] = -0.7921778546; rk[7] = -0.8086048979; rk[8] = -0.8223534102;
	  rk[9] = -0.8340533676;

	  double one = 1.;
	  double rrk = (pow(rk[k-1],expo)-1)/(rk[k-1]-1);
	  double cpnt = rk[k-1]*vlb[v1-1];
	  double dpnt = rk[k-1]*vub[v1-1];

	  out << "subject to ";
	  cname << "Rcvx_oddpow_" << v1 << "_" << expo << "_1";
	  cnamelist << ", " << cname.str();
	  out << cname.str() << " : "
	      << "Rcvx_L[" << v1 << "]^" << expo << "*(1 + " << rk[k-1]
	      << "*(" << VarIndexToString(v1) << "/Rcvx_L[" << v1
	      << "] - 1)) >= " << VarIndexToString(1/coeff, addedvarindex)
	      << ";" << endl;
	  cname.str("");

	  out << "subject to ";
	  cname << "Rcvx_oddpow_" << v1 << "_" << expo << "_2";
	  cnamelist << ", " << cname.str();
	  out << cname.str() << " : "
	      << "Rcvx_U[" << v1 << "]^" << expo << "*(1 + " << rk[k-1]
	      << "*(" << VarIndexToString(v1) << "/Rcvx_U[" << v1
	      << "] - 1)) <= " << VarIndexToString(1/coeff, addedvarindex)
	      << ";" << endl;
	  cname.str("");

	  out << "subject to ";
	  cname << "Rcvx_oddpow_" << v1 << "_" << expo << "_3";
	  cnamelist << ", " << cname.str();
	  out << cname.str() << " : "
	      << expo << "*Rcvx_U[" << v1 << "]^" << expo-1
	      << "*" << VarIndexToString(v1) << " - 2*" << expo
	      << "*Rcvx_U[" << v1 << "]^" << expo << " <= "
	      << VarIndexToString(1/coeff, addedvarindex) << ";" << endl;
	  cname.str("");

	  out << "subject to ";
	  cname << "Rcvx_oddpow_" << v1 << "_" << expo << "_4";
	  cnamelist << ", " << cname.str();
	  out << cname.str() << " : "
	      << expo << "*Rcvx_L[" << v1 << "]^" << expo-1
	      << "*" << VarIndexToString(v1) << " - 2*" << expo
	      << "*Rcvx_L[" << v1 << "]^" << expo << " >= "
	      << VarIndexToString(1/coeff, addedvarindex) << ";" << endl;
	  cname.str();
	}
      } // end if(op)
    }
  }

  out << "problem Rcvx : Rcvx_x, ";
  if (!InProb->IsProblemContinuous()) {
    out << "Rcvx_xInt, ";
  }
  out << "Rcvx_w, Rcvx_obj" << cnamelist.str();
  out << ";" << endl;

  out.close();

  // output data file
  ofstream datout(OutFile2.c_str());
  datout << "# convexified .dat file produced by ROSE (Rconvexifiermod solver)"
	 << endl;
  datout << "param : Rcvx_L Rcvx_U :=" << endl;
  for(int i = 0; i < NumberOfVariables; i++) {
    datout << " " << i+1 << " " << vlb[i] << " " << vub[i] << endl;
  }
  datout << ";" << endl;
  datout.close();

  ////////////////////////////////
  // after solution: set optimal values in problem
  IsSolved = true;
  for(int i = 0; i < NumberOfVariables; i++) {
    InProb->SetCurrentVariableValueLI(i + 1, xstar[i]);
  }
  return ret;
}

void RconvexifiermodSolver::SetOptimizationDirection(int theoptdir) {
  assert(theoptdir == Minimization || theoptdir == Maximization);
  ManualOptDir = true;
  bool changedir = false;
  if (OptDir != theoptdir) {
    changedir = true;
  }
  OptDir = theoptdir;
  if (OptDir == Minimization) {
    OptDirCoeff = 1;
  } else {
    OptDirCoeff = -1;
  }
  if (changedir && InitProb) {
    // set opt dir in the local solver
  }
}

void RconvexifiermodSolver::GetSolution(map<int,double>& objfunval,
			     map<int,double>& solution) {
  if (IsSolved) {
    objfunval[1] = OptimalObjVal;
    for(int i = 0; i < NumberOfVariables; i++) {
      solution[InProb->GetVariableID(i + 1)] = xstar[i];
    }
  } else {
    objfunval[1] = ROSEINFINITY;
    for(int i = 0; i < NumberOfVariables; i++) {
      solution[InProb->GetVariableID(i + 1)] = ROSEINFINITY;
    }
  }
}

void RconvexifiermodSolver::GetSolutionLI(vector<double>& objfunval,
			       vector<double>& solution) {
  if (IsSolved) {
    objfunval[1] = OptimalObjVal;
    for(int i = 0; i < NumberOfVariables; i++) {
      solution[i] = xstar[i];
    }
  } else {
    objfunval[1] = ROSEINFINITY;
    for(int i = 0; i < NumberOfVariables; i++) {
      solution[i] = ROSEINFINITY;
    }
  }
}

void RconvexifiermodSolver::GetSolutionLI(double& objfunval,
				   vector<double>& solution) {
  if (IsSolved) {
    objfunval = OptimalObjVal;
    for(int i = 0; i < NumberOfVariables; i++) {
      solution[i] = xstar[i];
    }
  } else {
    objfunval = ROSEINFINITY;
    for(int i = 0; i < NumberOfVariables; i++) {
      solution[i] = ROSEINFINITY;
    }
  }
}


void RconvexifiermodSolver::SetStartingPoint(int varID, double sp) {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

void RconvexifiermodSolver::SetStartingPointLI(int localindex, double sp) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

double RconvexifiermodSolver::GetStartingPoint(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

double RconvexifiermodSolver::GetStartingPointLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

void RconvexifiermodSolver::SetVariableLB(int varID, double LB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableLBLI(localindex, LB);
}

double RconvexifiermodSolver::GetVariableLB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void RconvexifiermodSolver::SetVariableUB(int varID, double UB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableUBLI(localindex, UB);
}

double RconvexifiermodSolver::GetVariableUB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}

void RconvexifiermodSolver::SetVariableLBLI(int localindex, double LB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vlb[localindex - 1] = LB;
}

double RconvexifiermodSolver::GetVariableLBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void RconvexifiermodSolver::SetVariableUBLI(int localindex, double UB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vub[localindex - 1] = UB;
}

double RconvexifiermodSolver::GetVariableUBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}
