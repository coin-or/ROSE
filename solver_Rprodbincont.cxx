/*
** Name:       solver_Rprodbincont.cxx
** Author:     Leo Liberti & Fabien Tarissan
** Source:     GNU C++
** Purpose:    ProdBinCont reformulation solver
** License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
** History:    071121 work started
*/

#include <cstdio>
#include <cstring>
#include <cassert>
#include <sys/time.h>
#include "newsolver.h"
#include "solver.h"
#include "solver_Rprodbincont.h"
#include "utils.h"

#define FIRSTOBJ 1
#define EPSILON 1e-4
#define WDEFNAMEBUFSIZE 256

// SNOPT6SOLVER CLASS METHODS

ProdBinContSolver::ProdBinContSolver() {
  TheName = "prodbincont";
  NumberOfVariables = 0;
  NumberOfConstraints = 0;
  NumberOfObjectives = 0;
  IsSolved = false;
  OptDir = Minimization;
  OptDirCoeff = 1;
  ManualOptDir = false;
  OptimalObjVal = ROSEINFINITY;
  Epsilon = EPSILON;
  LocalSolverName = "null";
  MaxRunningTime = 0;
}

ProdBinContSolver::~ProdBinContSolver() {
}

void ProdBinContSolver::SetProblem(Problem* p) {
  InProb = p; 
  InitProb = false;
}

bool ProdBinContSolver::CanSolve(int problemtype) {
  // fill in here with information about solver capabilities
}

void ProdBinContSolver::Initialize(bool force = false) {

  if (InitProb) {
    ;
  } else {
    InitProb = true;
  }

  // check the input problem is set to something
  if (InProb == NULL) {
    cerr << "ProdBinContSolver::Initialize(): error: InProb is NULL\n";
    assert(InProb != NULL);
  }

  NumberOfVariables = InProb->GetNumberOfVariables();
  NumberOfConstraints = InProb->GetNumberOfConstraints();
  NumberOfObjectives = InProb->GetNumberOfObjectives();

  // parameters
  ReplaceParams(InProb->GetParamsRef());
  int np = ParameterBlob.GetNumberOfParameters();
  for(int i = 1; i <= np; i++) {
    if (ParameterBlob.GetParameterName(i) == "ProdBinContLocalSolver") {
      LocalSolverName = ParameterBlob.GetParameterStringValue(i);
    } else if (ParameterBlob.GetParameterName(i) == "ProdBinContMaxTime") {
      MaxRunningTime = (double) ParameterBlob.GetParameterIntValue(i);
      assert(MaxRunningTime >= 0);
    } else if (ParameterBlob.GetParameterName(i) == "ProdBinContEpsilon") {
      Epsilon = ParameterBlob.GetParameterDoubleValue(i);
      assert(Epsilon > 0);
    } else if (ParameterBlob.GetParameterName(i) == "ProdBinContQuiet" ||
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

  /*
  // create and configure local solver
  LocalSolver = NewSolver(LocalSolverName);
  LocalSolver->SetProblem(InProb);
  LocalSolver->ReplaceParams(ParameterBlob);
  */

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

int ProdBinContSolver::Solve(bool reinitialize) {

  int ret = 0;
  if (reinitialize || !InitProb) {
    Initialize();
  }

  // do it



  // maps containing integrality, variable bounds and var names 
  map<int,double> mlb, mub;
  map<int,bool> mint;
  map<int,string> mstr;
  for(int i = 1; i <= NumberOfVariables; i++) {
    Variable* theVar = InProb->GetVariableLI(i);
    mlb[theVar->ID] = theVar->LB;
    mub[theVar->ID] = theVar->UB;
    mint[theVar->ID] = theVar->IsIntegral;
    mstr[theVar->ID] = theVar->Name;
  }

#undef sprintf
  // find variables involved in bilinear terms
  char wdefnbuf[WDEFNAMEBUFSIZE];
  string vn = "w";
  string cn = "c";
  int vi = NumberOfVariables + 1;
  int found;
  double wlb, wub;
  vector<Expression> defcons;
  map<int,pair<double,double> > addvarbounds;
  for(int i = 1; i <= NumberOfObjectives; i++) {
    Objective* theObj = InProb->GetObjectiveLI(i);
    if (!theObj->Function->IsLinear()) {
      vi += theObj->Function->ProdBinCont(vi, vn, mint, mlb, mub, mstr, 
					  addvarbounds, defcons);
      // add variables
      for(map<int,pair<double,double> >::iterator mi = addvarbounds.begin();
	  mi != addvarbounds.end(); mi++) {
	int i = mi->first;
	wlb = mi->second.first;
	wub = mi->second.second;
	//	sprintf(wdefnbuf, "w%d", vi + i);
	sprintf(wdefnbuf, "w%d", i);
	vn = wdefnbuf;
	InProb->AddVariable(vn, false, false, wlb, wub, (wlb+wub)/2);
      }
      // add constraints
      for(int i = 0; i < defcons.size(); i++) {
	InProb->AddConstraint(cn, defcons[i], -ROSEINFINITY, 0);
      }
    }
  vn = "w";    
  }
  for(int i = 1; i <= NumberOfConstraints; i++) {
    Constraint* theCon = InProb->GetConstraintLI(i);
    if (!theCon->Function->IsLinear()) {
      vi += theCon->Function->ProdBinCont(vi, vn, mint, mlb, mub, mstr,
					  addvarbounds, defcons);
      for(map<int,pair<double,double> >::iterator mi = addvarbounds.begin();
	  mi != addvarbounds.end(); mi++) {
	int i = mi->first;
	wlb = mi->second.first;
	wub = mi->second.second;
	//	sprintf(wdefnbuf, "w%d", vi + i);
	sprintf(wdefnbuf, "w%d", i);
	vn = wdefnbuf;
	InProb->AddVariable(vn, false, false, wlb, wub, (wlb+wub)/2);
      }
      for(int i = 0; i < defcons.size(); i++) {
	InProb->AddConstraint(cn, defcons[i], -ROSEINFINITY, 0);
      }
    }
    vn = "w";
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

  return ret;
}

void ProdBinContSolver::SetOptimizationDirection(int theoptdir) {
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

void ProdBinContSolver::GetSolution(map<int,double>& objfunval, 
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

void ProdBinContSolver::GetSolutionLI(vector<double>& objfunval, 
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

void ProdBinContSolver::GetSolutionLI(double& objfunval, 
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


void ProdBinContSolver::SetStartingPoint(int varID, double sp) {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

void ProdBinContSolver::SetStartingPointLI(int localindex, double sp) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

double ProdBinContSolver::GetStartingPoint(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

double ProdBinContSolver::GetStartingPointLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

void ProdBinContSolver::SetVariableLB(int varID, double LB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableLBLI(localindex, LB);
}

double ProdBinContSolver::GetVariableLB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void ProdBinContSolver::SetVariableUB(int varID, double UB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableUBLI(localindex, UB);
}

double ProdBinContSolver::GetVariableUB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}

void ProdBinContSolver::SetVariableLBLI(int localindex, double LB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vlb[localindex - 1] = LB;
}

double ProdBinContSolver::GetVariableLBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void ProdBinContSolver::SetVariableUBLI(int localindex, double UB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vub[localindex - 1] = UB;
}

double ProdBinContSolver::GetVariableUBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}
