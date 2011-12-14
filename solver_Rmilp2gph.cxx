/*
** Name:       solver_Rmilp2gph.cxx
** Author:     Leo Liberti
** Source:     GNU C++
** Purpose:    Solver Rmilp2gph implementation 
** License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
** History:    111214 work started
*/

#include <cstdio>
#include <cstring>
#include <cassert>
#include <sys/time.h>
#include "newsolver.h"
#include "solver.h"
#include "solver_Rmilp2gph.h"
#include "utils.h"

#define FIRSTOBJ 1
#define EPSILON 1e-4

// SNOPT6SOLVER CLASS METHODS

Rmilp2gphSolver::Rmilp2gphSolver() {
  TheName = "Rmilp2gph";
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

Rmilp2gphSolver::~Rmilp2gphSolver() {
}

void Rmilp2gphSolver::SetProblem(Problem* p) {
  InProb = p; 
  InitProb = false;
}

bool Rmilp2gphSolver::CanSolve(int problemtype) {
  // fill in here with information about solver capabilities
}

void Rmilp2gphSolver::Initialize(bool force = false) {

  if (InitProb) {
    ;
  } else {
    InitProb = true;
  }

  // check the input problem is set to something
  if (InProb == NULL) {
    cerr << "Rmilp2gphSolver::Initialize(): error: InProb is NULL\n";
    assert(InProb != NULL);
  }

  NumberOfVariables = InProb->GetNumberOfVariables();
  NumberOfConstraints = InProb->GetNumberOfConstraints();
  NumberOfObjectives = InProb->GetNumberOfObjectives();

  // parameters
  ReplaceParams(InProb->GetParamsRef());
  int np = ParameterBlob.GetNumberOfParameters();
  for(int i = 1; i <= np; i++) {
    if (ParameterBlob.GetParameterName(i) == "Rmilp2gphLocalSolver") {
      LocalSolverName = ParameterBlob.GetParameterStringValue(i);
    } else if (ParameterBlob.GetParameterName(i) == "Rmilp2gphMaxTime") {
      MaxRunningTime = (double) ParameterBlob.GetParameterIntValue(i);
      assert(MaxRunningTime >= 0);
    } else if (ParameterBlob.GetParameterName(i) == "Rmilp2gphEpsilon") {
      Epsilon = ParameterBlob.GetParameterDoubleValue(i);
      assert(Epsilon > 0);
    } else if (ParameterBlob.GetParameterName(i) == "Rmilp2gphQuiet" ||
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

  /* ACTIVATE IN CASE OF A COMPLEX SOLVER WITH SUBSOLVER(s)
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

int Rmilp2gphSolver::Solve(bool reinitialize) {

  int ret = 0;
  if (reinitialize || !InitProb) {
    Initialize();
  }

  // do it
  
  // sparse constraint matrix
  map<int, vector<pair<int, double> > > wstar;
  map<int,int> vi2idx;
  int ncols = 1;
  int nrows = 0;
  int h;
  vector<double> lincoeff;
  vector<int> linvi;
  vector<string> linvn;
  double c;
  int linconstrsize = 0;
  // check that there are constraints
  if (NumberOfConstraints > 0) {
    for(int i = 1; i <= NumberOfConstraints; i++) {
      Constraint* theCon = InProb->GetConstraintLI(i);
      // count linear constraints
      if (theCon->Function->IsLinear()) {
	linconstrsize++;
      }
    }
    // check that there are some linear constraints
    if (linconstrsize > 0) {
      // read sparse linear constraint matrix
      for(int i = 1; i <= NumberOfConstraints; i++) {
	Constraint* theCon = InProb->GetConstraintLI(i);
	if (theCon->Function->IsLinear()) {
	  nrows++;
	  theCon->Function->GetLinearInfo(lincoeff, linvi, linvn, c);
	  for(int j = 0; j < linvi.size(); j++) {
	    // map varindices on another index from 1 to nvars
	    h = vi2idx[linvi[j]];
	    if (h == 0) {
	      vi2idx[linvi[j]] = ncols;
	      h = ncols;
	      ncols++;
	    }
	    pair<int,double> p(i - 1, lincoeff[j]);
	    wstar[h - 1].push_back(p);
	  }
	}
      }
      // write .gph file to Rmilp2gph_out.gph
      OutFile = "Rmilp2gph_out.gph";
      ofstream out(OutFile.c_str());
      out << "Undirected;" << endl;
      out << "Vertices:" << endl;
      for(int i = 0; i < ncols + nrows; i++) {
	out << i << endl;
      }
      out << "Edges:" << endl;
      for(map<int, vector<pair<int, double> > >::iterator mi = wstar.begin();
	  mi != wstar.end(); mi++) {
	for(vector<pair<int,double> >::iterator vi = mi->second.begin();
	    vi != mi->second.end(); vi++) {
	  out << mi->first << " ";
	  out << vi->first + ncols << " ";
	  out << vi->second << " 1" << endl;
	}
      }
      out.close();
    } else {
      if (!Quiet) {
	cout << "Rmilp2gphSolver: no linear constraints, nothing to do\n";
      }
    }
  } else {
    if (!Quiet) {
      cout << "Rmilp2gphSolver: unconstrained problem, nothing to do\n";
    }
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
    cout << "Rmilp2gphSolver: best solution f* = " << OptimalObjVal << endl;
  }
  InProb->SetOptimalObjectiveValue(FIRSTOBJ, OptimalObjVal);
  InProb->SetSolved(true);    
  if (IsFeasible) {
    InProb->SetFeasible(1);
  } else {
    if (!Quiet) {
      cout << "Rmilp2gphSolver: this solution is infeasible by at least " 
           << discrepancy << endl;
    }
    InProb->SetFeasible(-1);
  }

  return ret;
}

void Rmilp2gphSolver::SetOptimizationDirection(int theoptdir) {
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

void Rmilp2gphSolver::GetSolution(map<int,double>& objfunval, 
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

void Rmilp2gphSolver::GetSolutionLI(vector<double>& objfunval, 
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

void Rmilp2gphSolver::GetSolutionLI(double& objfunval, 
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


void Rmilp2gphSolver::SetStartingPoint(int varID, double sp) {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

void Rmilp2gphSolver::SetStartingPointLI(int localindex, double sp) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

double Rmilp2gphSolver::GetStartingPoint(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

double Rmilp2gphSolver::GetStartingPointLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

void Rmilp2gphSolver::SetVariableLB(int varID, double LB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableLBLI(localindex, LB);
}

double Rmilp2gphSolver::GetVariableLB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void Rmilp2gphSolver::SetVariableUB(int varID, double UB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableUBLI(localindex, UB);
}

double Rmilp2gphSolver::GetVariableUB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}

void Rmilp2gphSolver::SetVariableLBLI(int localindex, double LB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vlb[localindex - 1] = LB;
}

double Rmilp2gphSolver::GetVariableLBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void Rmilp2gphSolver::SetVariableUBLI(int localindex, double UB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vub[localindex - 1] = UB;
}

double Rmilp2gphSolver::GetVariableUBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}
