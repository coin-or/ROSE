/*
** Name:       solver_gomory.cxx
** Author:     Leo Liberti
** Source:     GNU C++
** Purpose:    Solver Gomory implementation
** License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
** History:    051003 work started (Gomory cutting plane algorithm for MILPs)
*/

#include <cstdio>
#include <cstring>
#include <cassert>
#include <sys/time.h>
#include "newsolver.h"
#include "solver.h"
#include "solver_gomory.h"
#include "utils.h"

#define FIRSTOBJ 1
#define EPSILON 1e-4

GomorySolver::GomorySolver() {
  TheName = "gomory";
  NumberOfVariables = 0;
  NumberOfConstraints = 0;
  IsSolved = false;
  OptDir = Minimization;
  OptDirCoeff = 1;
  ManualOptDir = false;
  OptimalObjVal = ROSEINFINITY;
  Epsilon = EPSILON;
  LocalSolverName = "glpk";
  LocalSolver = NULL;
  MaxRunningTime = 0;
  Relaxed = false;
}

GomorySolver::~GomorySolver() {
  if (LocalSolver) {
    DeleteSolver(LocalSolver);
  }
}

void GomorySolver::SetProblem(Problem* p) {
  InProb = p; 
  InitProb = false;
}

bool GomorySolver::CanSolve(int problemtype) {
  if (problemtype == ROSE_LP || problemtype == ROSE_MILP) {
    return true;
  } else {
    return false;
  }
}

void GomorySolver::Initialize(bool force = false) {

  if (InitProb) {
    ;
  } else {
    InitProb = true;
  }

  // check the input problem is set to something
  if (InProb == NULL) {
    cerr << "GomorySolver::Initialize(): error: InProb is NULL\n";
    assert(InProb != NULL);
  }

  NumberOfVariables = InProb->GetNumberOfVariables();
  NumberOfConstraints = InProb->GetNumberOfConstraints();

  // parameters
  ReplaceParams(InProb->GetParamsRef());
  int np = ParameterBlob.GetNumberOfParameters();
  for(int i = 1; i <= np; i++) {
    if (ParameterBlob.GetParameterName(i) == "GomoryLocalSolver") {
      LocalSolverName = ParameterBlob.GetParameterStringValue(i);
    } else if (ParameterBlob.GetParameterName(i) == "GomoryMaxTime") {
      MaxRunningTime = (double) ParameterBlob.GetParameterIntValue(i);
      assert(MaxRunningTime >= 0);
    } else if (ParameterBlob.GetParameterName(i) == "GomoryEpsilon") {
      Epsilon = ParameterBlob.GetParameterDoubleValue(i);
      assert(Epsilon > 0);
    } else if (ParameterBlob.GetParameterName(i) == "RelaxInteger") {
      Relaxed = true;
    } else if (ParameterBlob.GetParameterName(i) == "GomoryQuiet" ||
               ParameterBlob.GetParameterName(i) == "Quiet") {
      Quiet = ParameterBlob.GetParameterBoolValue(i);
    }
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

  // variable integrality
  bool ismilp = false;
  bool tmpb;
  for(int i = 1; i <= NumberOfVariables; i++) {
    tmpb = InProb->GetVariableLI(i)->IsIntegral;
    integrality.push_back(tmpb);
    if (Relaxed) {
      integrality[i - 1] = false;
    }
    if (!ismilp && integrality[i - 1]) {
      ismilp = true;
    }
  }
  if (!ismilp) {
    cout << "GomorySolver: WARNING: problem is an LP, Gomory is overkill\n";
  }
}

int GomorySolver::Solve(bool reinitialize) {

  int ret = 0;
  if (reinitialize || !InitProb) {
    Initialize();
  }

  // timings
  double tstart[2], tend[2];
  double starttime, stoptime, usertimelap;
  starttime = getcputime(tstart);

  // do it
  bool termination = false;
  bool feasibility = false;
  vector<pair<int, double> > cutinfo;
  pair<int, double> p;
  double b;
  int cisize;
  int NumberOfCuts = 0;
  bool tmpbas, tmpint;
  while(!termination) {
    ret = LocalSolver->Solve();
    LocalSolver->GetSolutionLI(f, x);
    if (!Quiet) {
      cout << "GomorySolver: calling local solver, ret = " << ret 
	   << ", f = " << f << endl;
    }
    termination = true;
    for(int i = 0; i < NumberOfVariables; i++) {
      tmpbas = LocalSolver->IsBasicLI(i + 1);
      tmpint = integrality[i]; 
      if (tmpbas && tmpint && fabs(x[i] - rint(x[i])) > Epsilon) {
	// found non-integer variable: generate cut
	termination = false;
	b = LocalSolver->GetSimplexTableauRow(i + 1, cutinfo);
	// make cut
	cisize = cutinfo.size();
	for(int j = 0; j < cisize; j++) {
	  cutinfo[j].second = floor(cutinfo[j].second) - cutinfo[j].second;
	}
	// add the new slack variable for the cut
	b = floor(b) - b;
	if (!Quiet) {
	  cout << "GomorySolver: Gomory from basic " 
	       << i + 1 << ": ";
	  for(int j = 0; j < cisize; j++) {
	    cout << "(" << cutinfo[j].second << ")*";
	    if (cutinfo[j].first <= NumberOfVariables) {
	      cout << "x" << cutinfo[j].first;
	    } else {
	      cout << "c" 
		   << cutinfo[j].first - NumberOfVariables;
	    }
	    if (j < cisize - 1) {
	      cout << " + ";
	    }
	  }
	  cout << " <= " << b << endl;
	}
	LocalSolver->AddCut(cutinfo, -ROSEINFINITY, b);
	NumberOfCuts++;
	break;
      }      
    }
    if (MaxRunningTime > 0) {
      stoptime = getcputime(tend);
      usertimelap = tend[0] - tstart[0];
      if (usertimelap > (double) MaxRunningTime) {
	termination = true;
	if (!Quiet) {
	  cout << "GomorySolver: running time exceeds limit, exiting\n";
	}
      }
    }
    if (termination) {
      // means all variables are integer, update optimal solution
      fstar = f;
      for(int i = 0; i < NumberOfVariables; i++) {
	xstar[i] = x[i];
      }
    }
  }
  
  // final stopwatch
  stoptime = getcputime(tend);
  usertimelap = tend[0] - tstart[0];
  if (!Quiet) {
    cout << "GomorySolver: finished search after " << usertimelap
	 << " seconds of user CPU, f* = " << fstar << endl;
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
    cout << "GomorySolver: best solution f* = " << OptimalObjVal << endl;
  }
  InProb->SetOptimalObjectiveValue(FIRSTOBJ, OptimalObjVal);
  InProb->SetSolved(true);    
  if (IsFeasible) {
    InProb->SetFeasible(1);
  } else {
    if (!Quiet) {
      cout << "GomorySolver: this solution is infeasible by at least " 
           << discrepancy << endl;
    }
    InProb->SetFeasible(-1);
  }

  return ret;
}

void GomorySolver::SetOptimizationDirection(int theoptdir) {
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

void GomorySolver::GetSolution(map<int,double>& objfunval, 
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

void GomorySolver::GetSolutionLI(vector<double>& objfunval, 
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

void GomorySolver::GetSolutionLI(double& objfunval, 
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


void GomorySolver::SetStartingPoint(int varID, double sp) {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

void GomorySolver::SetStartingPointLI(int localindex, double sp) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

double GomorySolver::GetStartingPoint(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

double GomorySolver::GetStartingPointLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

void GomorySolver::SetVariableLB(int varID, double LB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableLBLI(localindex, LB);
}

double GomorySolver::GetVariableLB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void GomorySolver::SetVariableUB(int varID, double UB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableUBLI(localindex, UB);
}

double GomorySolver::GetVariableUB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}

void GomorySolver::SetVariableLBLI(int localindex, double LB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vlb[localindex - 1] = LB;
}

double GomorySolver::GetVariableLBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void GomorySolver::SetVariableUBLI(int localindex, double UB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vub[localindex - 1] = UB;
}

double GomorySolver::GetVariableUBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}
