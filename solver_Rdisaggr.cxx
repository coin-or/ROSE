/*
** Name:       solver_Rdisaggr.cxx
** Author:     Leo Liberti 
** Source:     GNU C++
** Purpose:    disaggregating reformulation
** License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
** History:    080228 work started
*/

#include <cstdio>
#include <cstring>
#include <cassert>
#include <sys/time.h>
#include "newsolver.h"
#include "solver.h"
#include "solver_Rdisaggr.h"
#include "utils.h"

#define FIRSTOBJ 1
#define EPSILON 1e-4

DisaggrSolver::DisaggrSolver() {
  TheName = "disaggr";
  NumberOfVariables = 0;
  NumberOfConstraints = 0;
  NumberOfObjectives = 0;
  IsSolved = false;
  OptDir = Minimization;
  OptDirCoeff = 1;
  ManualOptDir = false;
  OptimalObjVal = ROSEINFINITY;
  Epsilon = EPSILON;
  LocalSolverName = "snopt";
  MaxRunningTime = 0;
  K = 5;
}

DisaggrSolver::~DisaggrSolver() {
}

void DisaggrSolver::SetProblem(Problem* p) {
  InProb = p; 
  InitProb = false;
}

bool DisaggrSolver::CanSolve(int problemtype) {
  // fill in here with information about solver capabilities
}

void DisaggrSolver::Initialize(bool force = false) {

  if (InitProb) {
    ;
  } else {
    InitProb = true;
  }

  // check the input problem is set to something
  if (InProb == NULL) {
    cerr << "DisaggrSolver::Initialize(): error: InProb is NULL\n";
    assert(InProb != NULL);
  }

  NumberOfVariables = InProb->GetNumberOfVariables();
  NumberOfConstraints = InProb->GetNumberOfConstraints();
  NumberOfObjectives = InProb->GetNumberOfObjectives();

  // parameters
  ReplaceParams(InProb->GetParamsRef());
  int np = ParameterBlob.GetNumberOfParameters();
  for(int i = 1; i <= np; i++) {
    if (ParameterBlob.GetParameterName(i) == "DisaggrLocalSolver") {
      LocalSolverName = ParameterBlob.GetParameterStringValue(i);
    } else if (ParameterBlob.GetParameterName(i) == "DisaggrMaxTime") {
      MaxRunningTime = (double) ParameterBlob.GetParameterIntValue(i);
      assert(MaxRunningTime >= 0);
    } else if (ParameterBlob.GetParameterName(i) == "DisaggrEpsilon") {
      Epsilon = ParameterBlob.GetParameterDoubleValue(i);
      assert(Epsilon > 0);
    } else if (ParameterBlob.GetParameterName(i) == "DisaggrQuiet" ||
               ParameterBlob.GetParameterName(i) == "Quiet") {
      Quiet = ParameterBlob.GetParameterBoolValue(i);
    } else if (ParameterBlob.GetParameterName(i) == "MainSolver") {
      MainSolver = ParameterBlob.GetParameterStringValue(i);
    } else if (ParameterBlob.GetParameterName(i) == "AmplFlag") {
      AmplFlag = ParameterBlob.GetParameterBoolValue(i);
    } else if (ParameterBlob.GetParameterName(i) == "DisaggrK") {
      K = ParameterBlob.GetParameterIntValue(i);
    }
  }

  if (MainSolver != TheName && AmplFlag) {
    Quiet = true;
  }

  // variable integrality
  for(int i = 1; i <= NumberOfVariables; i++) {
    integrality.push_back(InProb->GetVariableLI(i)->IsIntegral);
  }

  // create and configure local solver
  LocalSolver = NewSolver(LocalSolverName);
  LocalSolver->SetProblem(InProb);
  LocalSolver->ReplaceParams(ParameterBlob);

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

}

int DisaggrSolver::Solve(bool reinitialize) {

  int ret = 0;
  if (reinitialize || !InitProb) {
    Initialize();
  }

  // do it

  // find variables involved in bilinear terms
  string vn;
  vector<int> vidx;
  Expression x(1.0, 1, "x");
  Expression y(1.0, 2, "y");
  Expression schema = x * y;
  for(int i = 1; i <= NumberOfObjectives; i++) {
    Objective* theObj = InProb->GetObjectiveLI(i);
    if (!theObj->Function->IsLinear()) {
      theObj->Function->GetVarIndicesInSchema(vidx, schema);
    }
  }
  for(int i = 1; i <= NumberOfConstraints; i++) {
    Constraint* theCon = InProb->GetConstraintLI(i);
    if (!theCon->Function->IsLinear()) {
      theCon->Function->GetVarIndicesInSchema(vidx, schema);
    }
  }

  // replace appropriate vars with disaggregation
  for(int i = 0; i < vidx.size(); i ++) {
    Variable* theVar = InProb->GetVariable(vidx[i]);
    if (!theVar->IsIntegral) {
      // only if var is continuous
      Expression needle(1.0, vidx[i], "x");
      // prepare sum of disaggregation
      Expression replace(theVar->LB);
      double chi = (theVar->UB - theVar->LB) / (double) K;
      for(int k = 0; k < K; k++) {
	vn = "d";
	InProb->AddVariable(vn, false, false, 0.0, chi, chi/2);
	Expression d(1.0, InProb->GetNumberOfVariables(), "d");
	replace = replace + d;
      }
      // replace obj objs and constrs
      for(int i = 1; i <= NumberOfObjectives; i++) {
	Objective* theObj = InProb->GetObjectiveLI(i);
	if (!theObj->Function->IsLinear()) {
	  theObj->Function->ReplaceSubexpression(needle, replace);
	  while(theObj->Function->DistributeProductsOverSums());
	}
      }
      for(int i = 1; i <= NumberOfConstraints; i++) {
	Constraint* theCon = InProb->GetConstraintLI(i);
	if (!theCon->Function->IsLinear()) {
	  theCon->Function->ReplaceSubexpression(needle, replace);
	  while(theCon->Function->DistributeProductsOverSums());
	}
      }
      // add constraint replace - needle = 0
      Expression needledef = replace - needle;
      vn = "c";
      InProb->AddConstraint(vn, needledef, 0.0, 0.0);
    }
  }

  if (!Quiet) {
    cout << *InProb;
  }
  
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
    cout << "DisaggrSolver: best solution f* = " << OptimalObjVal << endl;
  }
  InProb->SetOptimalObjectiveValue(FIRSTOBJ, OptimalObjVal);
  InProb->SetSolved(true);    
  if (IsFeasible) {
    InProb->SetFeasible(1);
  } else {
    if (!Quiet) {
      cout << "DisaggrSolver: this solution is infeasible by at least " 
           << discrepancy << endl;
    }
    InProb->SetFeasible(-1);
  }

  return ret;
}

void DisaggrSolver::SetOptimizationDirection(int theoptdir) {
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

void DisaggrSolver::GetSolution(map<int,double>& objfunval, 
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

void DisaggrSolver::GetSolutionLI(vector<double>& objfunval, 
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

void DisaggrSolver::GetSolutionLI(double& objfunval, 
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


void DisaggrSolver::SetStartingPoint(int varID, double sp) {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

void DisaggrSolver::SetStartingPointLI(int localindex, double sp) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

double DisaggrSolver::GetStartingPoint(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

double DisaggrSolver::GetStartingPointLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

void DisaggrSolver::SetVariableLB(int varID, double LB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableLBLI(localindex, LB);
}

double DisaggrSolver::GetVariableLB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void DisaggrSolver::SetVariableUB(int varID, double UB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableUBLI(localindex, UB);
}

double DisaggrSolver::GetVariableUB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}

void DisaggrSolver::SetVariableLBLI(int localindex, double LB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vlb[localindex - 1] = LB;
}

double DisaggrSolver::GetVariableLBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void DisaggrSolver::SetVariableUBLI(int localindex, double UB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vub[localindex - 1] = UB;
}

double DisaggrSolver::GetVariableUBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}
