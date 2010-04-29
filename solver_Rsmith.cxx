/*
** Name:       solver_Rsmith.cxx
** Author:     Leo Liberti
** Source:     GNU C++
** Purpose:    Smith standard form reformulation solver
** License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
** History:    080225 work started
*/

#include <cstdio>
#include <cstring>
#include <cassert>
#include <sys/time.h>
#include "newsolver.h"
#include "solver.h"
#include "solver_Rsmith.h"
#include "utils.h"

#define FIRSTOBJ 1
#define EPSILON 1e-4
#define WDEFNAMEBUFSIZE 256

// SNOPT6SOLVER CLASS METHODS

SmithSolver::SmithSolver() {
  TheName = "smith";
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
  InProb = NULL;
}

SmithSolver::~SmithSolver() {
  // localsolver not used by this solver
  /*if (LocalSolver) {
    DeleteSolver(LocalSolver);
  }*/
}

void SmithSolver::SetProblem(Problem* p) {
  InProb = p; 
  InitProb = false;
}

bool SmithSolver::CanSolve(int problemtype) {
  // fill in here with information about solver capabilities
}

void SmithSolver::Initialize(bool force = false) {
  using namespace std;
  if (InitProb) {
    ;
  } else {
    InitProb = true;
  }

  // check the input problem is set to something
  if (InProb == NULL) {
    cerr << "SmithSolver::Initialize(): error: InProb is NULL\n";
    assert(InProb != NULL);
  }

  NumberOfVariables = InProb->GetNumberOfVariables();
  NumberOfConstraints = InProb->GetNumberOfConstraints();
  NumberOfObjectives = InProb->GetNumberOfObjectives();

  // parameters
  ReplaceParams(InProb->GetParamsRef());
  int np = ParameterBlob.GetNumberOfParameters();
  for(int i = 1; i <= np; i++) {
    if (ParameterBlob.GetParameterName(i) == "SmithMaxTime") {
      MaxRunningTime = (double) ParameterBlob.GetParameterIntValue(i);
      assert(MaxRunningTime >= 0);
    } else if (ParameterBlob.GetParameterName(i) == "SmithEpsilon") {
      Epsilon = ParameterBlob.GetParameterDoubleValue(i);
      assert(Epsilon > 0);
    } else if (ParameterBlob.GetParameterName(i) == "SmithQuiet" ||
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

}

int SmithSolver::Solve(bool reinitialize) {
  using namespace std;
  int ret = 0;
  if (reinitialize || !InitProb) {
    Initialize();
  }
  // maps containing the variable bounds (used for interval arithmetics)
  map<int,double> mlb, mub;

  // do it
  // get problem ranges
  for(int i = 1; i <= NumberOfVariables; i++) {
    Variable* theVar = InProb->GetVariableLI(i);
    mlb[theVar->ID] = theVar->LB;
    mub[theVar->ID] = theVar->UB;
  }
  
  // schemata dealt with by Smith's reformulation
  vector<Expression> schemata;
  Expression x(1,1,"x");
  x->SetExponent(2);
  schemata.push_back(x);
  // operators dealt with by Smith's reformulation
  vector<int> oplabels;
  oplabels.push_back(SUM);
  oplabels.push_back(PRODUCT);
  oplabels.push_back(FRACTION);
  oplabels.push_back(POWER);
  oplabels.push_back(LOG);
  oplabels.push_back(EXP);
  oplabels.push_back(SQRT);
  // stores the intermediate expressions
  vector<Expression> defcons;
  // added variable starting index
  int wi = NumberOfVariables + 1;
  // default added variable name
  string wn;
  // default definining constraint name
  char wdefnbuf[WDEFNAMEBUFSIZE];
  string wdefname;
  // store added variable bounds
  double wlb, wub;
  // added variable counter
  int waddcounter = NumberOfVariables + 1;
  int wfirstindex = waddcounter;
  int defconcurrent = 0;
  int addconstr = 0;
#undef sprintf
  // standardize objectives
  for(int i = 1; i <= NumberOfObjectives; i++) {
    Objective* theObj = InProb->GetObjectiveLI(i);
    if (!theObj->Function->IsLinear()) {
      // obj fun is nonlinear, standardize
      wn = "w";
      theObj->Function->SmithStandardForm(wi, wfirstindex, 
					  wn, oplabels, schemata, defcons);
      // for each operator subexpression
      for(int d = defconcurrent; d < defcons.size(); d++) {
	// compute added var bounds by interval arithmetic
	defcons[d]->Interval(mlb, mub, wlb, wub);
	if(wub > ROSEINFINITY) {
	  wub = ROSEINFINITY;
	}
	if(wlb < -ROSEINFINITY) {
	  wlb = -ROSEINFINITY;
	}
	// add new bounds to v bound maps
	mlb[wi] = wlb;
	mub[wi] = wub;
	sprintf(wdefnbuf, "w%d", waddcounter);
	wn = wdefnbuf;
	// add the new variable
	InProb->AddVariable(wn, false, false, wlb, wub, 0);
	sprintf(wdefnbuf, "wdef%d", wi);
	wdefname = wdefnbuf;
	// create the expression for the new defining constraint
	Expression wvar(1.0, wi, wn);
	Expression wdefexpr = defcons[d] - wvar;
	// add the defining constraint
	InProb->AddConstraint(wdefname, wdefexpr, 0, 0); 
	addconstr++;
	// increase counters
	wi++;
	waddcounter++;
      }
      defconcurrent = defcons.size();
    }
  }

  // standardize constraints
  for(int i = 1; i <= NumberOfConstraints; i++) {
    Constraint* theCon = InProb->GetConstraintLI(i);
    if (!theCon->Function->IsLinear()) {
      wn = "w";
      theCon->Function->SmithStandardForm(wi, wfirstindex, 
					  wn, oplabels, schemata, defcons);
      for(int d = defconcurrent; d < defcons.size(); d++) {
	defcons[d]->Interval(mlb, mub, wlb, wub);
	if(wub > ROSEINFINITY) {
	  wub = ROSEINFINITY;
	}
	if(wlb < -ROSEINFINITY) {
	  wlb = -ROSEINFINITY;
	}
	mlb[wi] = wlb;
	mub[wi] = wub;
	sprintf(wdefnbuf, "w%d", waddcounter);
	wn = wdefnbuf;
	InProb->AddVariable(wn, false, false, wlb, wub, 0);
	sprintf(wdefnbuf, "wdef%d", wi);
	wdefname = wdefnbuf;
	Expression wvar(1.0, wi, wn);
	Expression wdefexpr = defcons[d] - wvar;
	InProb->AddConstraint(wdefname, wdefexpr, 0, 0); 
	addconstr++;
	wi++;
	waddcounter++;
      }
      defconcurrent = defcons.size();
    }
  }
  InProb->Simplifier(true);

  //if (!Quiet) {
    cout << *InProb;
  //}

  return ret;
}

void SmithSolver::SetOptimizationDirection(int theoptdir) {
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

void SmithSolver::GetSolution(std::map<int,double>& objfunval, 
			      std::map<int,double>& solution) {
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

void SmithSolver::GetSolutionLI(std::vector<double>& objfunval, 
				std::vector<double>& solution) {
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

void SmithSolver::GetSolutionLI(double& objfunval, 
				std::vector<double>& solution) {
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

void SmithSolver::SetStartingPoint(int varID, double sp) {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

void SmithSolver::SetStartingPointLI(int localindex, double sp) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

double SmithSolver::GetStartingPoint(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

double SmithSolver::GetStartingPointLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

void SmithSolver::SetVariableLB(int varID, double LB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableLBLI(localindex, LB);
}

double SmithSolver::GetVariableLB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void SmithSolver::SetVariableUB(int varID, double UB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableUBLI(localindex, UB);
}

double SmithSolver::GetVariableUB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}

void SmithSolver::SetVariableLBLI(int localindex, double LB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vlb[localindex - 1] = LB;
}

double SmithSolver::GetVariableLBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void SmithSolver::SetVariableUBLI(int localindex, double UB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vub[localindex - 1] = UB;
}

double SmithSolver::GetVariableUBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}
