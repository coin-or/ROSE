/*
** Name:       solver_oa.cxx
** Author:     Leo Liberti & Claudia D'Ambrosio
** Source:     GNU C++
** Purpose:    Solver Oa implementation
** License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
** History:    080707 work started
*/

#include <cstdio>
#include <cstring>
#include <cassert>
#include <fstream>
#include <sys/time.h>
#include "newsolver.h"
#include "solver.h"
#include "solver_Roa.h"
#include "utils.h"

#define FIRSTOBJ 1
#define EPSILON 1e-4

#define CNAMEBUFSIZE 64

// ROASOLVER CLASS METHODS

OaSolver::OaSolver() {
  TheName = "oa";
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
  NumberOfReformulations = 0;
  OutFile = "Roa.out";
}

OaSolver::~OaSolver() {
}

void OaSolver::SetProblem(Problem* p) {
  InProb = p; 
  InitProb = false;
}

bool OaSolver::CanSolve(int problemtype) {
  // fill in here with information about solver capabilities
}

void OaSolver::Initialize(bool force = false) {

  if (InitProb) {
    ;
  } else {
    InitProb = true;
  }

  // check the input problem is set to something
  if (InProb == NULL) {
    cerr << "OaSolver::Initialize(): error: InProb is NULL\n";
    assert(InProb != NULL);
  }

  NumberOfVariables = InProb->GetNumberOfVariables();
  NumberOfConstraints = InProb->GetNumberOfConstraints();
  NumberOfObjectives = InProb->GetNumberOfObjectives();

  // parameters
  ReplaceParams(InProb->GetParamsRef());
  int np = ParameterBlob.GetNumberOfParameters();
  for(int i = 1; i <= np; i++) {
    if (ParameterBlob.GetParameterName(i) == "OaLocalSolver") {
      LocalSolverName = ParameterBlob.GetParameterStringValue(i);
    } else if (ParameterBlob.GetParameterName(i) == "OaMaxTime") {
      MaxRunningTime = (double) ParameterBlob.GetParameterIntValue(i);
      assert(MaxRunningTime >= 0);
    } else if (ParameterBlob.GetParameterName(i) == "OaEpsilon") {
      Epsilon = ParameterBlob.GetParameterDoubleValue(i);
      assert(Epsilon > 0);
    } else if (ParameterBlob.GetParameterName(i) == "OaQuiet" ||
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

int OaSolver::Solve(bool reinitialize) {

  int ret = 0;
  if (reinitialize || !InitProb) {
    Initialize();
  }

  // do it

  // add an AMPL suffix var.amplsolverindex
  Ampl::SufDesc* vamplindex = suf_get_ASL(Ampl::asl, "amplsolverindex", 
				     Ampl::ASL_Sufkind_var);
  vamplindex->u.i = (int*) M1zapalloc_ASL(&Ampl::asl->i, 
				     NumberOfVariables * sizeof(int));

  // set initial variable values to current variable values
  for(int i = 1; i <= NumberOfVariables; i++) {
    InProb->SetCurrentVariableValueLI(i, x[i-1]);
  }
  ofstream out(OutFile.c_str());

  // check for evidently convex constraints
  double discrepancy;
  int cID;
  int vID;
  bool isactive;
  int uplow;
  char cname[CNAMEBUFSIZE];
#undef sprintf
  int constr_counter=0;
  int temp = NumberOfVariables;
  for(int i = 1; i <= NumberOfConstraints; i++) {
    Constraint* theCon = InProb->GetConstraintLI(i);
    if (!theCon->Function->IsLinear()) {
      if ((theCon->Function->IsEvidentlyConvex() && 
	   theCon->UB < ROSEINFINITY) ||
	  (theCon->Function->IsEvidentlyConcave() &&
	   theCon->LB > -ROSEINFINITY)) {
	cID = InProb->GetConstraintID(i);
	isactive = InProb->IsConstrActive(cID, Epsilon, uplow);
//	if (isactive && uplow == 1) {
        if(InProb->TestConstraintsFeasibility(InProb->GetConstrLocalIndex(cID),Epsilon,discrepancy)){
           constr_counter++;
	}	
      }
    }
  }

  double **coeff;
  coeff = new double*[constr_counter];
  for(int i=0; i<constr_counter; i++) {
    coeff[i] = new double[temp];
  }
  double * bnd = new double[constr_counter];
  int temp_counter=0;

  out << "param card_OA := " << constr_counter << ";\n\n";
  out << "param coeff_OA : ";
  for(int i = 1; i <= NumberOfVariables; i++) {
    out << " " << i;
  }
  out << " :=\n";
  for(int i = 1; i <= NumberOfConstraints; i++) {
    Constraint* theCon = InProb->GetConstraintLI(i);
    if (!theCon->Function->IsLinear()) {
      if ((theCon->Function->IsEvidentlyConvex() && 
	   theCon->UB < ROSEINFINITY) ||
	  (theCon->Function->IsEvidentlyConcave() &&
	   theCon->LB > -ROSEINFINITY)) {
	cID = InProb->GetConstraintID(i);
	isactive = InProb->IsConstrActive(cID, Epsilon, uplow);
//	if (isactive && uplow == 1) {
        if(InProb->TestConstraintsFeasibility(InProb->GetConstrLocalIndex(cID),Epsilon,discrepancy)){
          temp_counter++;
	  bnd[temp_counter-1] = 0; // theCon->UB - InProb->EvalConstr(cID); // always =
	  // generate OA constraint at x
	  for(int j = 1; j <= NumberOfVariables; j++) {
	    vID = InProb->GetVarLocalIndex(j);
	    out << "# " << j << "\t" << InProb->GetVariableID(j) << "\t" << InProb->GetVarLocalIndex(j) <<"\n";
	    coeff[temp_counter-1][vID-1] = InProb->EvalConstrDiff(cID, vID);
	    out << "# " << j << "\t" << x[j-1] << "\t" << coeff[temp_counter-1][vID-1]<<"\n";
	    bnd[temp_counter-1] += coeff[temp_counter-1][vID-1]*x[j-1]; // a '+' because it's -coeff*x[j] on LHS
//            if(fabs(coeff)>1/ROSEINFINITY) 
//	      out << "+(" << coeff << ")*_FP_MILP_x" << vID << " ";
	  }
	}	
      }
    }
  }
  for(int i = 0; i < constr_counter; i++) {
    out << i+1 << " ";
    for(int j = 0; j < NumberOfVariables; j++) {
      out << coeff[i][j] << " ";
    }
    out << "\n";
  }
  out << ";\n";
  out << "\nparam: rhs_OA := ";
  for(int i = 0; i < constr_counter; i++) {
    out << i+1 << " ";
    out << bnd[i] << "\n";
  }
  out << ";\n";

  out.close();
  free(coeff);
  free(bnd);

  // return variable ampl solver index
  for(int i = 0; i < NumberOfVariables; i++) {
      vamplindex->u.i[i] = i+1;
  }

  // after solution: set optimal values in problem
  IsSolved = true;
  IsFeasible = false;
  discrepancy = ROSEINFINITY;
  for(int i = 0; i < NumberOfVariables; i++) {
    InProb->SetCurrentVariableValueLI(i+1, i+1); // this is NOT a mistake!
    InProb->SetOptimalVariableValueLI(i+1, i+1); // this is NOT a mistake!
  }
  OptimalObjVal = InProb->EvalObj(FIRSTOBJ);
  if (!Quiet) {
    cout << "OaSolver: best solution f* = " << OptimalObjVal << endl;
  }
  InProb->SetOptimalObjectiveValue(FIRSTOBJ, OptimalObjVal);
  InProb->SetSolved(true);    
  InProb->SetFeasible(1);

  return ret;
}

int OaSolver::GetNumberOfReformulations(void) const {
  return NumberOfReformulations;
}

Problem* OaSolver::GetReformulation(int ridx) {
  assert(ridx >= 0 && ridx <= 1);
  return InProb;
}

void OaSolver::SetOptimizationDirection(int theoptdir) {
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

void OaSolver::GetSolution(map<int,double>& objfunval, 
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

void OaSolver::GetSolutionLI(vector<double>& objfunval, 
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

void OaSolver::GetSolutionLI(double& objfunval, 
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


void OaSolver::SetStartingPoint(int varID, double sp) {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

void OaSolver::SetStartingPointLI(int localindex, double sp) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

double OaSolver::GetStartingPoint(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

double OaSolver::GetStartingPointLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

void OaSolver::SetVariableLB(int varID, double LB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableLBLI(localindex, LB);
}

double OaSolver::GetVariableLB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void OaSolver::SetVariableUB(int varID, double UB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableUBLI(localindex, UB);
}

double OaSolver::GetVariableUB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}

void OaSolver::SetVariableLBLI(int localindex, double LB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vlb[localindex - 1] = LB;
}

double OaSolver::GetVariableLBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void OaSolver::SetVariableUBLI(int localindex, double UB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vub[localindex - 1] = UB;
}

double OaSolver::GetVariableUBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}
