/*
** Name:       solver_Rprintmod.cxx
** Author:     Leo Liberti 
** Source:     GNU C++
** Purpose:    print flat form in AMPL .mod format
** License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
** History:    080225 work started
*/

#include <fstream>
#include <cstdio>
#include <cstring>
#include <cassert>
#include <sys/time.h>
#include "newsolver.h"
#include "solver.h"
#include "solver_Rprintmod.h"
#include "utils.h"

#define FIRSTOBJ 1
#define EPSILON 1e-6

PrintmodSolver::PrintmodSolver() {
  TheName = "printmod";
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
  OutFile = "Rprintmod_out.mod";
}

PrintmodSolver::~PrintmodSolver() {
  if (LocalSolver) {
    DeleteSolver(LocalSolver);
  }
}

void PrintmodSolver::SetProblem(Problem* p) {
  InProb = p; 
  InitProb = false;
}

bool PrintmodSolver::CanSolve(int problemtype) {
  // fill in here with information about solver capabilities
}

void PrintmodSolver::Initialize(bool force = false) {

  if (InitProb) {
    ;
  } else {
    InitProb = true;
  }

  // check the input problem is set to something
  if (InProb == NULL) {
    cerr << "PrintmodSolver::Initialize(): error: InProb is NULL\n";
    assert(InProb != NULL);
  }

  NumberOfVariables = InProb->GetNumberOfVariables();
  NumberOfConstraints = InProb->GetNumberOfConstraints();
  NumberOfObjectives = InProb->GetNumberOfObjectives();

  // parameters
  ReplaceParams(InProb->GetParamsRef());
  int np = ParameterBlob.GetNumberOfParameters();
  for(int i = 1; i <= np; i++) {
    if (ParameterBlob.GetParameterName(i) == "PrintLocalSolver") {
      LocalSolverName = ParameterBlob.GetParameterStringValue(i);
    } else if (ParameterBlob.GetParameterName(i) == "PrintMaxTime") {
      MaxRunningTime = (double) ParameterBlob.GetParameterIntValue(i);
      assert(MaxRunningTime >= 0);
    } else if (ParameterBlob.GetParameterName(i) == "PrintEpsilon") {
      Epsilon = ParameterBlob.GetParameterDoubleValue(i);
      assert(Epsilon > 0);
    } else if (ParameterBlob.GetParameterName(i) == "PrintQuiet" ||
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

int PrintmodSolver::Solve(bool reinitialize) {

  int ret = 0;
  if (reinitialize || !InitProb) {
    Initialize();
  }

  // do it
  if (!Quiet) {
  //  cout << "SymmgroupSolver: writing AMPL .dat file " << OutFile 
  //	 << " for findsymm.mod" << endl;
  }
  ofstream out(OutFile.c_str());
  out << "# flat form AMPL .mod file, automatic translation provided by rose" 
      << endl;
  out << "# P has " << NumberOfVariables << " variables, "
       << NumberOfObjectives << " objectives, "
       << NumberOfConstraints << " constraints" << endl << endl;
  for(int i = 1; i <= NumberOfVariables; i++) {
    out << "var x" << InProb->GetVariableID(i) << " >= " 
	<< vlb[i-1] << ", <= " << vub[i-1]; 
    if (integrality[i-1]) {
      if (fabs(vlb[i-1]) < EPSILON && fabs(vub[i-1] - 1) < EPSILON) {
	out << ", binary";
      } else {
	out << ", integer";
      }
    }
    out << ";" << endl;
  }
  out << endl;
  for(int i = 1; i <= NumberOfObjectives; i++) {
    Objective* theObj = InProb->GetObjectiveLI(i);
    int maxvi = theObj->Function->NumberOfVariables();
    theObj->Function->ResetVarNames("x", 1, maxvi);
    int optdir = InProb->GetOptimizationDirectionLI(i);
    if (optdir == 0) {
      // minimization
      out << "minimize ";
    } else if (optdir == 1) {
      // minimization
      out << "maximize ";
    }
    out << "obj" << i << ": " << theObj->Function->ToString() << ";" << endl;
  }
  out << endl;
  for(int i = 1; i <= NumberOfConstraints; i++) {
    Constraint* theCon = InProb->GetConstraintLI(i);
    int maxvi = theCon->Function->NumberOfVariables();
    theCon->Function->ResetVarNames("x", 1, maxvi);
    if (theCon->LB <= -ROSEINFINITY && theCon->UB >= ROSEINFINITY) {
      continue;
    } else {
      out << "subject to c" << i << ": ";
      out << theCon->Function->ToString();
      if (theCon->LB > -ROSEINFINITY && theCon->UB >= ROSEINFINITY) {
	// \ge
	out << " >= " << theCon->LB << ";";
      } else if (theCon->LB <= -ROSEINFINITY && theCon->UB < ROSEINFINITY) {
	// \le
	out << " <= " << theCon->UB << ";";
      } else {
	// ranged
	if (fabs(theCon->UB - theCon->LB) < EPSILON) {
	  out << " = " << theCon->LB << ";";
	} else {
	  out << " >= " << theCon->LB << ";" << endl;
	  out << "subject to c" << i << "b: ";
	  out << theCon->Function->ToString();
	  out << " <= " << theCon->UB << ";";
	}
      }
      out << endl;
    }
  }
  out << endl;
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

  return ret;
}

void PrintmodSolver::SetOptimizationDirection(int theoptdir) {
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

void PrintmodSolver::GetSolution(map<int,double>& objfunval, 
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

void PrintmodSolver::GetSolutionLI(vector<double>& objfunval, 
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

void PrintmodSolver::GetSolutionLI(double& objfunval, 
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


void PrintmodSolver::SetStartingPoint(int varID, double sp) {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

void PrintmodSolver::SetStartingPointLI(int localindex, double sp) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

double PrintmodSolver::GetStartingPoint(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

double PrintmodSolver::GetStartingPointLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

void PrintmodSolver::SetVariableLB(int varID, double LB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableLBLI(localindex, LB);
}

double PrintmodSolver::GetVariableLB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void PrintmodSolver::SetVariableUB(int varID, double UB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableUBLI(localindex, UB);
}

double PrintmodSolver::GetVariableUB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}

void PrintmodSolver::SetVariableLBLI(int localindex, double LB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vlb[localindex - 1] = LB;
}

double PrintmodSolver::GetVariableLBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void PrintmodSolver::SetVariableUBLI(int localindex, double UB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vub[localindex - 1] = UB;
}

double PrintmodSolver::GetVariableUBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}
