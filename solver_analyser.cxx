/*
** Name:       solver_analyser.cxx
** Author:     Leo Liberti
** Source:     GNU C++
** Purpose:    Analyser solver - kind of "dr ampl"
** License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
** History:    071120 work started
*/

#include <cstdio>
#include <cstring>
#include <cassert>
#include <sys/time.h>
#include "newsolver.h"
#include "solver.h"
#include "solver_analyser.h"
#include "utils.h"

#define FIRSTOBJ 1
#define EPSILON 1e-3



AnalyserSolver::AnalyserSolver() {
  TheName = "analyser";
  NumberOfVariables = 0;
  NumberOfConstraints = 0;
  IsSolved = false;
  OptDir = Minimization;
  OptDirCoeff = 1;
  ManualOptDir = false;
  OptimalObjVal = ROSEINFINITY;
  Epsilon = EPSILON;
  LocalSolverName = "none";
  MaxRunningTime = 0;
}

AnalyserSolver::~AnalyserSolver() {
}

void AnalyserSolver::SetProblem(Problem* p) {
  InProb = p; 
  InitProb = false;
}

bool AnalyserSolver::CanSolve(int problemtype) {
  // fill in here with information about solver capabilities
}

void AnalyserSolver::Initialize(bool force = false) {

  if (InitProb) {
    ;
  } else {
    InitProb = true;
  }

  // check the input problem is set to something
  if (InProb == NULL) {
    cerr << "AnalyserSolver::Initialize(): error: InProb is NULL\n";
    assert(InProb != NULL);
  }

  NumberOfVariables = InProb->GetNumberOfVariables();
  NumberOfConstraints = InProb->GetNumberOfConstraints();
  NumberOfObjectives = InProb->GetNumberOfObjectives();

  // parameters
  ReplaceParams(InProb->GetParamsRef());
  int np = ParameterBlob.GetNumberOfParameters();
  for(int i = 1; i <= np; i++) {
    if (ParameterBlob.GetParameterName(i) == "AnalyserLocalSolver") {
      LocalSolverName = ParameterBlob.GetParameterStringValue(i);
    } else if (ParameterBlob.GetParameterName(i) == "AnalyserMaxTime") {
      MaxRunningTime = (double) ParameterBlob.GetParameterIntValue(i);
      assert(MaxRunningTime >= 0);
    } else if (ParameterBlob.GetParameterName(i) == "AnalyserEpsilon") {
      Epsilon = ParameterBlob.GetParameterDoubleValue(i);
      assert(Epsilon > 0);
    } else if (ParameterBlob.GetParameterName(i) == "AnalyserQuiet" ||
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

int AnalyserSolver::Solve(bool reinitialize) {

  int ret = 0;
  if (reinitialize || !InitProb) {
    Initialize();
  }

  // do it

  // add an AMPL suffix var.integer
  Ampl::SufDesc* vtype = suf_get_ASL(Ampl::asl, "integer", 
				     Ampl::ASL_Sufkind_var);
  // add an AMPL suffix constr.linear
  Ampl::SufDesc* ctype = suf_get_ASL(Ampl::asl, "linear", 
				     Ampl::ASL_Sufkind_con); 
  // add an AMPL suffix constr.feasible
  Ampl::SufDesc* cfeas = suf_get_ASL(Ampl::asl, "feasible", 
				     Ampl::ASL_Sufkind_con); 
  // add an AMPL suffix constr.evidently_convex
  Ampl::SufDesc* cevconv = suf_get_ASL(Ampl::asl, "evconvex", 
				     Ampl::ASL_Sufkind_con); 
  // add an AMPL suffix constr.evidently_convex
  Ampl::SufDesc* cevconc = suf_get_ASL(Ampl::asl, "evconcave", 
				     Ampl::ASL_Sufkind_con); 
  // add an AMPL suffix obj.linear
  Ampl::SufDesc* otype = suf_get_ASL(Ampl::asl, "linear", 
				     Ampl::ASL_Sufkind_obj); 
  // add an AMPL suffix constr.infeasquant
  Ampl::SufDesc* cinfeasquant = suf_get_ASL(Ampl::asl, "infeasquant", 
				     Ampl::ASL_Sufkind_con); 

  // store space for suffixes
  vtype->u.i = (int*) M1zapalloc_ASL(&Ampl::asl->i, 
				     NumberOfVariables * sizeof(int));
  ctype->u.i = (int*) M1zapalloc_ASL(&Ampl::asl->i, 
				     NumberOfConstraints * sizeof(int));
  cfeas->u.i = (int*) M1zapalloc_ASL(&Ampl::asl->i, 
				     NumberOfConstraints * sizeof(int));
  cevconv->u.i = (int*) M1zapalloc_ASL(&Ampl::asl->i, 
				     NumberOfConstraints * sizeof(int));
  cevconc->u.i = (int*) M1zapalloc_ASL(&Ampl::asl->i, 
				     NumberOfConstraints * sizeof(int));
  otype->u.i = (int*) M1zapalloc_ASL(&Ampl::asl->i, 
				     NumberOfObjectives * sizeof(int));
  cinfeasquant->u.r = (real*) M1zapalloc_ASL(&Ampl::asl->i, 
				     NumberOfConstraints * sizeof(real));


  if (!Quiet) {
    cout << "AnalyserSolver: performing problem analysis" << endl;
  }
  // return variable integrality
  for(int i = 0; i < NumberOfVariables; i++) {
    if (integrality[i]) {
      vtype->u.i[i] = 1;
    } else {
      vtype->u.i[i] = 0;
    }
  }

  // return constraint linearity/nonlinearity
  for(int i = 1; i <= NumberOfConstraints; i++) {
    if (InProb->GetConstraintLI(i)->Function->IsLinear()) {
      ctype->u.i[i-1] = 1;
    } else {
      ctype->u.i[i-1] = 0;
    }
  }


  // return constraint evident convexity
  for(int i = 1; i <= NumberOfConstraints; i++) {
    if (InProb->GetConstraintLI(i)->Function->IsEvidentlyConvex()) {
      cevconv->u.i[i-1] = 1;
    } else {
      cevconv->u.i[i-1] = 0;
    }
  }

  // return constraint evident concavity
  for(int i = 1; i <= NumberOfConstraints; i++) {
    if (InProb->GetConstraintLI(i)->Function->IsEvidentlyConcave()) {
      cevconc->u.i[i-1] = 1;
    } else {
      cevconc->u.i[i-1] = 0;
    }
  }

  // return objective linearity/nonlinearity
  for(int i = 1; i <= NumberOfObjectives; i++) {
    if (InProb->GetObjectiveLI(i)->Function->IsLinear()) {
      otype->u.i[i-1] = 1;
    } else {
      otype->u.i[i-1] = 0;
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

  // updating the suffix .feasible for each constraint
  discrepancy = 0;
  // return constraint feasibility
  for(int i = 1; i <= NumberOfConstraints; i++) {
    if (InProb->TestConstraintsFeasibility(i,Epsilon,discrepancy)) {
      cfeas->u.i[i-1] = 1;
    } else {
      cfeas->u.i[i-1] = 0;
    }
  }

  for(int i = 0; i < NumberOfVariables; i++) {
    InProb->SetOptimalVariableValueLI(i + 1, xstar[i]);
  }
  OptimalObjVal = InProb->EvalObj(FIRSTOBJ);
  // updating the suffix .infeasquant for each constraint
  discrepancy = 0;
  // return constraint infeasibility
  for(int i = 1; i <= NumberOfConstraints; i++) {
    cinfeasquant->u.r[i-1] = (InProb->TestConstraintsInfeasibilityQuantity(i,Epsilon,discrepancy)); 
  }

  /*
  if (!Quiet) {
    cout << "AnalyserSolver: best solution f* = " << OptimalObjVal << endl;
  }
  */
  InProb->SetOptimalObjectiveValue(FIRSTOBJ, OptimalObjVal);
  InProb->SetSolved(true);    
  if (IsFeasible) {
    InProb->SetFeasible(1);
  } else {
    if (!Quiet) {
      /*
      cout << "AnalyserSolver: this solution is infeasible by at least " 
           << discrepancy << endl;
      */
    }
    InProb->SetFeasible(-1);
  }

  return ret;
}

void AnalyserSolver::SetOptimizationDirection(int theoptdir) {
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

void AnalyserSolver::GetSolution(map<int,double>& objfunval, 
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

void AnalyserSolver::GetSolutionLI(vector<double>& objfunval, 
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

void AnalyserSolver::GetSolutionLI(double& objfunval, 
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


void AnalyserSolver::SetStartingPoint(int varID, double sp) {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

void AnalyserSolver::SetStartingPointLI(int localindex, double sp) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

double AnalyserSolver::GetStartingPoint(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

double AnalyserSolver::GetStartingPointLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

void AnalyserSolver::SetVariableLB(int varID, double LB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableLBLI(localindex, LB);
}

double AnalyserSolver::GetVariableLB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void AnalyserSolver::SetVariableUB(int varID, double UB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableUBLI(localindex, UB);
}

double AnalyserSolver::GetVariableUB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}

void AnalyserSolver::SetVariableLBLI(int localindex, double LB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vlb[localindex - 1] = LB;
}

double AnalyserSolver::GetVariableLBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void AnalyserSolver::SetVariableUBLI(int localindex, double UB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vub[localindex - 1] = UB;
}

double AnalyserSolver::GetVariableUBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}
