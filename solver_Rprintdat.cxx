/*
** Name:       solver_Rprintdat.cxx
** Author:     Leo Liberti & Sonia Cafieri
** Source:     GNU C++
** Purpose:    Solver Rprintdat implementation / change Rprintdat with solvername
** License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
** History:    080917 work started
*/

#include <fstream>
#include <cstdio>
#include <cstring>
#include <cassert>
#include <sys/time.h>
#include "newsolver.h"
#include "solver.h"
#include "solver_Rprintdat.h"
#include "utils.h"

#define FIRSTOBJ 1
#define EPSILON 1e-4

// RPRINTDATSOLVER CLASS METHODS

RprintdatSolver::RprintdatSolver() {
  LocalSolver = NULL;
  TheName = "Rprintdat";
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
  OutFile = "Rprintdat_out.dat";
}

RprintdatSolver::~RprintdatSolver() {
  if (LocalSolver) {
    DeleteSolver(LocalSolver);
  }
}

void RprintdatSolver::SetProblem(Problem* p) {
  InProb = p; 
  InitProb = false;
}

bool RprintdatSolver::CanSolve(int problemtype) {
  // fill in here with information about solver capabilities
}

void RprintdatSolver::Initialize(bool force = false) {

  if (InitProb) {
    ;
  } else {
    InitProb = true;
  }

  // check the input problem is set to something
  if (InProb == NULL) {
    cerr << "RprintdatSolver::Initialize(): error: InProb is NULL\n";
    assert(InProb != NULL);
  }

  NumberOfVariables = InProb->GetNumberOfVariables();
  NumberOfConstraints = InProb->GetNumberOfConstraints();
  NumberOfObjectives = InProb->GetNumberOfObjectives();

  // parameters
  ReplaceParams(InProb->GetParamsRef());
  int np = ParameterBlob.GetNumberOfParameters();
  for(int i = 1; i <= np; i++) {
    if (ParameterBlob.GetParameterName(i) == "RprintdatMaxTime") {
      MaxRunningTime = (double) ParameterBlob.GetParameterIntValue(i);
      assert(MaxRunningTime >= 0);
    } else if (ParameterBlob.GetParameterName(i) == "RprintdatEpsilon") {
      Epsilon = ParameterBlob.GetParameterDoubleValue(i);
      assert(Epsilon > 0);
    } else if (ParameterBlob.GetParameterName(i) == "RprintdatQuiet" ||
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

int RprintdatSolver::Solve(bool reinitialize) {

  int ret = 0;
  if (reinitialize || !InitProb) {
    Initialize();
  }

  // do it
  
  ofstream out1("Rprintdat.mod");
  out1 << "param n >=1, integer;" << endl;
  out1 << "param m >=1, integer;\n" << endl;
  out1 << "set I := 1..n;" << endl;
  out1 << "set J := 1..m; \n" << endl;
  out1 << "param xL{i in I};" << endl;
  out1 << "param xU{i in I}; \n" << endl;
  out1 << "param c{i in I} default 0;" << endl;
  out1 << "param bl{j in J};" << endl;
  out1 << "param bu{j in J};" << endl;
  out1 << "param a{j in J, i in I} default 0;\n\n" << endl;
  out1 << "var x{i in I} >= xL[i], <= xU[i]; \n" << endl;
  out1 << "minimize Rpd_obj: sum{i in I: c[i]<>0} c[i]*x[i]; \n" << endl;
  out1 << "subject to con {j in J}: bl[j] <= sum{i in I: a[j,i]<>0} a[j,i]*x[i] <= bu[j];\n" << endl;
  out1.close(); 
  
  ofstream out("Rprintdat.dat");
  out << "# AMPL .dat file, automatic translation provided by rose\n" << endl;
  out << "param n := " << NumberOfVariables << ";" << endl;
  out << "param m := " << NumberOfConstraints << ";" << endl;
  out << "\n";
  out << "param xL := " << endl;
  for(int i = 1; i <= NumberOfVariables; i++) {
    out << i << " " << vlb[i-1] << endl;
  }
  out << ";" << endl;
  out << "\n";
  out << "param xU := " << endl;
  for(int i = 1; i <= NumberOfVariables; i++) {
    out << i << " " << vub[i-1] << endl;
  }
  out << ";" << endl;
  out << "\n";
  
  vector<double> lincoeff;
  vector<int> linvi;
  vector<string> linvn;
  double c = 0;
  out << "param c :=" << endl;
  for(int i = 1; i <= NumberOfObjectives; i++) {
    Objective* theObj = InProb->GetObjectiveLI(i);
    theObj->Function->GetLinearInfo(lincoeff, linvi, linvn, c);
    for(int j = 0; j < lincoeff.size(); j++) {
      out << linvi[j] << " " << lincoeff[j] << endl; 
    }
  }
  out << ";" << endl;

  // constraints  bounds
  double cl;
  double cu;
  out << "param : bl bu := " << endl;
  for(int i = 1; i <= NumberOfConstraints; i++) {
    cl = InProb->GetConstraintLI(i)->LB;
    cu = InProb->GetConstraintLI(i)->UB;
    if(cl < -1.e+18) {
      out << i << " " << "-Infinity";
    } else {
      out << i << " " << cl;
    }
    if(cu > 1.e+18) {
      out << " " << "Infinity" << endl;
    } else {
      out << " " << cu << endl; 
    }
  }
  out << ";" << endl;

  // constraints
  out << "param a := " << endl;
  for(int i = 1; i <= NumberOfConstraints; i++) {
    Constraint* theCon = InProb->GetConstraintLI(i);
    theCon->Function->GetLinearInfo(lincoeff, linvi, linvn, c);
    for(int j = 0; j < lincoeff.size(); j++) {
      out << i << " " << linvi[j] << " " << lincoeff[j] << endl; 
    }
  }
  out << ";" << endl;
  out.close();

  return ret;
}

void RprintdatSolver::SetOptimizationDirection(int theoptdir) {
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

void RprintdatSolver::GetSolution(map<int,double>& objfunval, 
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

void RprintdatSolver::GetSolutionLI(vector<double>& objfunval, 
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

void RprintdatSolver::GetSolutionLI(double& objfunval, 
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

void RprintdatSolver::SetStartingPoint(int varID, double sp) {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

void RprintdatSolver::SetStartingPointLI(int localindex, double sp) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

double RprintdatSolver::GetStartingPoint(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

double RprintdatSolver::GetStartingPointLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

void RprintdatSolver::SetVariableLB(int varID, double LB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableLBLI(localindex, LB);
}

double RprintdatSolver::GetVariableLB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void RprintdatSolver::SetVariableUB(int varID, double UB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableUBLI(localindex, UB);
}

double RprintdatSolver::GetVariableUB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}

void RprintdatSolver::SetVariableLBLI(int localindex, double LB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vlb[localindex - 1] = LB;
}

double RprintdatSolver::GetVariableLBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void RprintdatSolver::SetVariableUBLI(int localindex, double UB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vub[localindex - 1] = UB;
}

double RprintdatSolver::GetVariableUBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}
