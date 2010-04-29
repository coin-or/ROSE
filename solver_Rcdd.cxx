/*
** Name:       solver_Rcdd.cxx
** Author:     Leo Liberti & Sonia Cafieri
** Source:     GNU C++
** Purpose:    print out CDD .ine file
** License:    (C) Leo Liberti, all rights reserved. Code published under the
               Common Public License.
** History:    080924 work started
*/

#include <fstream>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <cassert>
#include <sys/time.h>
#include "newsolver.h"
#include "solver.h"
#include "solver_Rcdd.h"
#include "utils.h"

#define FIRSTOBJ 1
#define EPSILON 1e-4

// RCDDSOLVER CLASS METHODS

RcddSolver::RcddSolver() {
  TheName = "Rcdd";
  NumberOfVariables = 0;
  NumberOfConstraints = 0;
  NumberOfObjectives = 0;
  IsSolved = false;
  OptDir = Minimization;
  OptDirCoeff = 1;
  ManualOptDir = false;
  OptimalObjVal = ROSEINFINITY;
  Epsilon = EPSILON;
  MaxRunningTime = 0;
  //OutFile = "Rcdd_out.ine";
}

RcddSolver::~RcddSolver() {
}

void RcddSolver::SetProblem(Problem* p) {
  InProb = p;
  InitProb = false;
}

bool RcddSolver::CanSolve(int problemtype) {
  // fill in here with information about solver capabilities
}

void RcddSolver::Initialize(bool force = false) {


  if (InitProb) {
    ;
  } else {
    InitProb = true;
  }

  // check the input problem is set to something
  if (InProb == NULL) {
    cerr << "RcddSolver::Initialize(): error: InProb is NULL\n";
    assert(InProb != NULL);
  }

  NumberOfVariables = InProb->GetNumberOfVariables();
  NumberOfConstraints = InProb->GetNumberOfConstraints();
  NumberOfObjectives = InProb->GetNumberOfObjectives();

  // parameters
  ReplaceParams(InProb->GetParamsRef());
  int np = ParameterBlob.GetNumberOfParameters();
  for(int i = 1; i <= np; i++) {
    if (ParameterBlob.GetParameterName(i) == "RcddLocalSolver") {
      LocalSolverName = ParameterBlob.GetParameterStringValue(i);
    } else if (ParameterBlob.GetParameterName(i) == "RcddMaxTime") {
      MaxRunningTime = (double) ParameterBlob.GetParameterIntValue(i);
      assert(MaxRunningTime >= 0);
    } else if (ParameterBlob.GetParameterName(i) == "RcddEpsilon") {
      Epsilon = ParameterBlob.GetParameterDoubleValue(i);
      assert(Epsilon > 0);
    } else if (ParameterBlob.GetParameterName(i) == "RcddQuiet" ||
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

  // variable integrality
  for(int i = 1; i <= NumberOfVariables; i++) {
    integrality.push_back(InProb->GetVariableLI(i)->IsIntegral);
  }


  // current, best solution, bounds
  for(int i = 0; i < NumberOfVariables; i++) {
    x.push_back(InProb->GetStartingPointLI(i + 1));
    xstar.push_back(x[i]);
    vlb.push_back(InProb->GetVariableLI(i + 1)->LB);
    vub.push_back(InProb->GetVariableLI(i + 1)->UB);
    StartingPoint.push_back(x[i]);
  }
  f = ROSEINFINITY;
  fstar = ROSEINFINITY;



  // constraints  bounds
  cl = new double [NumberOfConstraints];
  cu = new double [NumberOfConstraints];
  for(int i = 0; i < NumberOfConstraints; i++) {
    cl[i] = InProb->GetConstraintLI(i + 1)->LB;
    cu[i] = InProb->GetConstraintLI(i + 1)->UB;
  }


  // get the rhs
  rhs = new double [NumberOfConstraints];
  int sgn [NumberOfConstraints];
  for(int i = 0; i < NumberOfConstraints; i++) {
    if (cl[i] < cu[i]) {
      if (cl[i] <= -ROSEINFINITY) {
	// case g(x) <= cu
	rhs[i] = cu[i];
        sgn[i] = 1;
      } else if (cu[i] >= ROSEINFINITY) {
	// case g(x) >= cl  => -g(x) <= - cl
	rhs[i] = -cl[i];
        sgn[i] = -1;
      }
    }
  }



  vector<Expression> ConstrLin;
  for(int i = 1; i <= NumberOfConstraints; i++) {
    ConstrLin.push_back(InProb->GetConstraintLI(i)->Function->GetLinearPart());
  }

  // calculate number of nonzeroes in constraint matrix
  int nzmat = 0;
  for(int i = 1; i <= NumberOfConstraints; i++) {
    if (ConstrLin[i - 1]->IsVariable()) {
      nzmat++;
    } else {
      nzmat += ConstrLin[i - 1]->GetSize();
    }
  }

  rowidx = new int [nzmat + 1];
  colidx = new int [nzmat + 1];
  a = new double [nzmat + 1];

  //int acnt = 0;
  // get nonzeroes in constraint matrix
  for(int i = 0; i < NumberOfConstraints; i++) {
    if (ConstrLin[i]->IsVariable()) {
      rowidx[acnt] = i;
      colidx[acnt] = ConstrLin[i]->GetVarIndex();
      colidx[acnt] = InProb->GetVarLocalIndex(colidx[acnt]);
      colidx[acnt] = colidx[acnt]-1;
      a[acnt] = sgn[i] * ConstrLin[i]->GetCoeff();
      acnt++;
    } else {
      int thecsize = ConstrLin[i]->GetSize();
      for(int j = 0; j < thecsize; j++) {
	rowidx[acnt] = i;
	colidx[acnt] = ConstrLin[i]->GetNode(j)->GetVarIndex();
	colidx[acnt] = InProb->GetVarLocalIndex(colidx[acnt]);
        colidx[acnt] = colidx[acnt]-1;
	a[acnt] = sgn[i]* ConstrLin[i]->GetNode(j)->GetCoeff();
	acnt++;
      }
    }
  }
/*
   cout << " acnt = " << acnt<< endl;
   cout << " a= " << endl;
   for(int i=0; i < acnt; i++) {
         cout << a[i] << " ";
   }
    cout << " \n";
   cout << " rowidx= " << endl;
   for(int i=0; i < acnt; i++) {
         cout << rowidx[i] << " ";
   }
    cout << " \n";
   cout << " colidx= " << endl;
   for(int i=0; i < acnt; i++) {
         cout << colidx[i] << " ";
   }
    cout << " \n\n";
*/

  // count the number of finite bounds
  nb = 0;
  for(int i = 0; i < NumberOfVariables; i++) {
      if(vlb[i] != -ROSEINFINITY) {
        nb ++;
      }
      if(vub[i] != ROSEINFINITY) {
        nb ++;
      }
  }


}

int RcddSolver::Solve(bool reinitialize) {

  int ret = 0;
  if (reinitialize || !InitProb) {
    Initialize();
  }

  // do it

  ///ofstream out(OutFile.c_str());
  ofstream out("incdd.ine");
  out << "CDD .ine file, automatic translation provided by rose" << endl;
  out << "H-representation" << endl;
  out << "begin" << endl;
  out << NumberOfConstraints+nb << " " << NumberOfVariables+1 << " rational" << endl ;

  for(int i = 0; i < NumberOfConstraints; i++) {

      if((int)rhs[i]==rhs[i]) {
         out << rhs[i] << " " ;
      } else {
         int dd = decimal(rhs[i]);
         out << rhs[i]*pow(10,(double)dd) <<"/" << pow(10, (double)dd) << " ";
      }
      for(int j = 0; j < NumberOfVariables; j++) {
          int flag = 0;
          for(int k = 0; k < acnt; k++) {
              if(rowidx[k]==i && colidx[k]==j) {
                 if((int)a[k]==a[k]) {
                    out << -a[k] << " ";
                 } else {
                    int dd = decimal(a[k]);
                    out << -a[k]*pow(10,(double)dd) <<"/" << pow(10, (double)dd) << " ";
                 }
                 flag = 1;
              }
          }
          if(flag == 0) {
             out << 0 << " ";
          }
      }

      out <<" \n" ;
  }

  for(int i = 0; i < NumberOfVariables; i++) {
      if(vlb[i] != -ROSEINFINITY) {
         if((int)vlb[i]==vlb[i]) {
             out << -vlb[i] << " ";
         } else {
             int dd = decimal(vlb[i]);
             out << -vlb[i]*pow(10,(double)dd) <<"/" << pow(10, (double)dd) << " ";
         }
         for(int j = 0; j < NumberOfVariables; j++) {
             if (j!=i) {
                 out << 0 << " ";
             } else {
                 out << 1 << " ";
             }
         }
         out <<" \n";
      }
      if(vub[i] != ROSEINFINITY) {
         if((int)vub[i]==vub[i]) {
             out << vub[i] << " ";
         } else {
             int dd = decimal(vub[i]);
             out << vub[i]*pow(10,(double)dd) <<"/" << pow(10, (double)dd) << " ";
         }
         for(int j = 0; j < NumberOfVariables; j++) {
             if (j!=i) {
                 out << 0 << " ";
             } else {
                 out << -1 << " ";
             }
         }
         out <<" \n";
      }
  }


  out << "end" << endl;

  out.close();



  // after solution: set optimal values in problem
  IsSolved = true;
  for(int i = 0; i < NumberOfVariables; i++) {
    InProb->SetCurrentVariableValueLI(i + 1, xstar[i]);
  }
  double discrepancy = 0;
  bool constrfeas = InProb->TestConstraintsFeasibility(Epsilon, discrepancy);
  bool varfeas = InProb->TestVariablesFeasibility(Epsilon, discrepancy);
  if (constrfeas && varfeas) {
    IsFeasible = 1;
  }
  else {
    IsFeasible = 0;
  }
  for(int i = 0; i < NumberOfVariables; i++) {
    InProb->SetOptimalVariableValueLI(i + 1, xstar[i]);
  }
  OptimalObjVal = InProb->EvalObj(FIRSTOBJ);
  if (!Quiet) {
    cout << "RcddSolver: best solution f* = " << OptimalObjVal << endl;
  }
  InProb->SetOptimalObjectiveValue(FIRSTOBJ, OptimalObjVal);
  InProb->SetSolved(true);
  if (IsFeasible) {
    InProb->SetFeasible(1);
  } else {
    if (!Quiet) {
      cout << "RcddSolver: this solution is infeasible by at least "
           << discrepancy << endl;
    }
    InProb->SetFeasible(-1);
  }

  return ret;
}

void RcddSolver::SetOptimizationDirection(int theoptdir) {
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

void RcddSolver::GetSolution(map<int,double>& objfunval,
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

void RcddSolver::GetSolutionLI(vector<double>& objfunval,
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

void RcddSolver::GetSolutionLI(double& objfunval,
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


void RcddSolver::SetStartingPoint(int varID, double sp) {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

void RcddSolver::SetStartingPointLI(int localindex, double sp) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

double RcddSolver::GetStartingPoint(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

double RcddSolver::GetStartingPointLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

void RcddSolver::SetVariableLB(int varID, double LB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableLBLI(localindex, LB);
}

double RcddSolver::GetVariableLB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void RcddSolver::SetVariableUB(int varID, double UB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableUBLI(localindex, UB);
}

double RcddSolver::GetVariableUB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}

void RcddSolver::SetVariableLBLI(int localindex, double LB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vlb[localindex - 1] = LB;
}

double RcddSolver::GetVariableLBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void RcddSolver::SetVariableUBLI(int localindex, double UB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vub[localindex - 1] = UB;
}

double RcddSolver::GetVariableUBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}


int RcddSolver::decimal(double n) {

  std::string s;
  // convert double n to string s
  { std::ostringstream ss;
    ss << n;
    s = ss.str();
  }

  const char *s1 = s.c_str(); // get const char * representation
  int len = strlen(s1) -1;

  int ii=0;
  while(s[ii]!='.') {
     ii++;
  }

  return len-ii;
}
