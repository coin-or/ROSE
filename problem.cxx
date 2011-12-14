/*
** Name:       problem.cxx
** Author:     Leo Liberti
** Source:     GNU C++
** Purpose:    implementation of the Problem class
** License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
** History:    050215 work started
*/

#include "problem.h"
#include <cassert>
#include <iomanip>

#define PARSEBUFFER 256000
#ifndef NULL
#define NULL 0
#endif
#define SYMSIMPTOLERANCE 1e-20

Problem::Problem() {
  Name = "My problem";
  IsProbContinuous = true;
  IsProbLinear = false;
  Parent = NULL;
  FormulationName = "original";
  NumberOfVariables = 0;
  NumberOfObjectives = 0;
  NumberOfConstraints = 0;
  x = NULL;
  Obj1AdditiveConstant = 0;
  DeletionFlag = false;
  ComputedSymbolicDiff = false;
  ComputedSymbolicDiff2 = false;
  Solved = false;
  Feasible = 0;
  NoSimplifier = false;
  NextFreeObjFunID = 0;
  NextFreeVarID = 0;
  NextFreeConstrID = 0;
}

Problem::Problem(bool thenosimpl) {
  Name = "My problem";
  IsProbContinuous = true;
  IsProbLinear = false;
  Parent = NULL;
  FormulationName = "original";
  NumberOfVariables = 0;
  NumberOfObjectives = 0;
  NumberOfConstraints = 0;
  x = NULL;
  Obj1AdditiveConstant = 0;
  DeletionFlag = false;
  ComputedSymbolicDiff = false;
  ComputedSymbolicDiff2 = false;
  Solved = false;
  Feasible = 0;
  NoSimplifier = thenosimpl;
  NextFreeObjFunID = 0;
  NextFreeVarID = 0;
  NextFreeConstrID = 0;
}

Problem::~Problem() {
  int s;
  s = Var.size();
  if (s > 0) {
    for(int i = 0; i < s; i++) {
      delete Var[i];
    }
  }
  s = Obj.size();
  if (s > 0) {
    for(int i = 0; i < s; i++) {
      delete Obj[i];
    }
  }
  s = Constr.size();
  if (s > 0) {
    for(int i = 0; i < s; i++) {
      delete Constr[i];
    }
  }
  if (x) {
    delete [] x;
  }
}

int Problem::GetNumberOfIntegerVariables(void) {
  int ret = 0;
  for(int i = 1; i <= NumberOfVariables; i++) {
    if (GetVariableLI(i)->IsIntegral) {
      ret++;
    }
  }      
  return ret;
}

int Problem::NewVariableID(void) {
  NextFreeVarID++;
  return NextFreeVarID;
}

int Problem::NewObjectiveID(void) {
  NextFreeObjFunID++;
  return NextFreeObjFunID;
}

int Problem::NewConstraintID(void) {
  NextFreeConstrID++;
  return NextFreeConstrID;
}

void Problem::AddVariable(std::string& myName, bool myIntegral, 
			  bool myPersistent, 
			  double myLB, double myUB, double myOptimum) {
  Variable* theVar = new Variable;
  theVar->IsIntegral = myIntegral;
  theVar->Persistent = myPersistent;
  theVar->Optimum = myOptimum;
  theVar->ID = NewVariableID();
  theVar->Name = myName;
  theVar->LB = myLB;
  theVar->UB = myUB;
  NumberOfVariables++;
  VarIDToLocalIndex[theVar->ID] = NumberOfVariables;
  Var.push_back(theVar);
  StartingPoint.push_back(myOptimum);
 
  double *new_x = new double[NumberOfVariables];

  for (int i = 0; i < NumberOfVariables - 1; i++) {
    new_x[i] = x[i];
  }

  delete[] x;

  x = new_x;
}

void Problem::AddVariable(std::string& myName, bool myIntegral, 
			  bool myPersistent,
			  bool myConstant, double myLB, 
			  double myUB, double myOptimum) {
  Variable* theVar = new Variable;
  theVar->IsIntegral = myIntegral;
  theVar->Persistent = myPersistent;
  theVar->IsConstant = myConstant;
  theVar->Optimum = myOptimum;
  theVar->ID = NewVariableID();
  theVar->Name = myName;
  theVar->LB = myLB;
  theVar->UB = myUB;
  NumberOfVariables++;
  VarIDToLocalIndex[theVar->ID] = NumberOfVariables;
  Var.push_back(theVar);
  StartingPoint.push_back(myOptimum);

  double *new_x = new double[NumberOfVariables];

  for (int i = 0; i < NumberOfVariables - 1; i++) {
    new_x[i] = x[i];
  }

  delete[] x;

  x = new_x;
}

void Problem::AddObjective(std::string& myName, Expression myExpr, int myOptDir, 
			   double myOptimum) {
  Objective* theObj = new Objective;
  theObj->ID = NewObjectiveID();
  theObj->OptDir = myOptDir;
  theObj->Name = myName;
  theObj->Optimum = myOptimum;
  theObj->Function.SetTo(myExpr);
  theObj->NonlinearPart.SetTo(theObj->Function->GetNonlinearPart());
  theObj->FunctionFET = NULL;
  theObj->NonlinearPartFET = NULL;
  NumberOfObjectives++;
  ObjIDToLocalIndex[theObj->ID] = NumberOfObjectives;
  Obj.push_back(theObj);
}

void Problem::AddConstraint(std::string& myName, Expression myExpr, 
			    double myLB, double myUB) {
  Constraint* theConstr = new Constraint;
  theConstr->ID = NewConstraintID();
  theConstr->Name = myName;
  theConstr->LB = myLB;
  theConstr->UB = myUB;
  theConstr->Function.SetTo(myExpr);
  theConstr->NonlinearPart.SetTo(theConstr->Function->GetNonlinearPart());
  theConstr->FunctionFET = NULL;
  theConstr->NonlinearPartFET = NULL;
  NumberOfConstraints++;
  ConstrIDToLocalIndex[theConstr->ID] = NumberOfConstraints;
  Constr.push_back(theConstr);
}

int Problem::GetOptimizationDirection(int objID) {
  int localindex = ObjIDToLocalIndex[objID];
  assert(localindex >= 1 && localindex <= NumberOfObjectives);
  return Obj[localindex - 1]->OptDir;
}

int Problem::GetOptimizationDirectionLI(int localindex) {
  assert(localindex >= 1 && localindex <= NumberOfObjectives);
  return Obj[localindex - 1]->OptDir;
}

void Problem::SetOptimizationDirection(int objID, int theOptDir) {
  int localindex = ObjIDToLocalIndex[objID];
  assert(localindex >= 1 && localindex <= NumberOfObjectives);
  assert(theOptDir == Minimization || theOptDir == Maximization);
  Obj[localindex - 1]->OptDir = theOptDir;
}

void Problem::SetOptimizationDirectionLI(int localindex, int theOptDir) {
  assert(localindex >= 1 && localindex <= NumberOfObjectives);
  assert(theOptDir == Minimization || theOptDir == Maximization);
  Obj[localindex - 1]->OptDir = theOptDir;
}

Problem* Problem::GetChild(int childindex) { 
  assert(childindex <= Children.size() && childindex >= 1);
  return Children[childindex - 1]; 
}

Variable* Problem::GetVariable(int varID) {
  int localindex = VarIDToLocalIndex[varID];  
  assert(localindex <= NumberOfVariables && localindex >= 1);
  return Var[localindex - 1];
}

Variable* Problem::GetVariableLI(int localindex) {
  assert(localindex <= NumberOfVariables && localindex >= 1);
  return Var[localindex - 1];
}

Objective* Problem::GetObjective(int objID) {
  int localindex = ObjIDToLocalIndex[objID];  
  assert(localindex <= NumberOfObjectives && localindex >= 1);
  return Obj[localindex - 1];
}

Objective* Problem::GetObjectiveLI(int localindex) {
  assert(localindex <= NumberOfObjectives && localindex >= 1);
  return Obj[localindex - 1];
}

Constraint* Problem::GetConstraint(int constrID) {
  int localindex = ConstrIDToLocalIndex[constrID];  
  assert(localindex <= NumberOfConstraints && localindex >= 1);
  return Constr[localindex - 1];
}

Constraint* Problem::GetConstraintLI(int localindex) {
  assert(localindex <= NumberOfConstraints && localindex >= 1);
  return Constr[localindex - 1];
}

int Problem::GetVariableID(int localindex) {
  assert(localindex <= NumberOfVariables && localindex >= 1);
  return Var[localindex - 1]->ID;
}

int Problem::GetObjectiveID(int localindex) {
  assert(localindex <= NumberOfObjectives && localindex >= 1);
  return Obj[localindex - 1]->ID;
}

int Problem::GetConstraintID(int localindex) {
  assert(localindex <= NumberOfConstraints && localindex >= 1);
  return Constr[localindex - 1]->ID;
}

int Problem::GetVarLocalIndex(int varID) {
  int localindex = VarIDToLocalIndex[varID];
  return localindex;
}

int Problem::GetObjLocalIndex(int objID) {
  int localindex = ObjIDToLocalIndex[objID];
  return localindex;
}

int Problem::GetConstrLocalIndex(int constrID) {
  int localindex = ConstrIDToLocalIndex[constrID];
  return localindex;
}

double Problem::GetOptimalVariableValue(int varID) {
  int localindex = VarIDToLocalIndex[varID];
  assert(localindex <= NumberOfVariables && localindex >= 1);
  return Var[localindex - 1]->Optimum;
}
double Problem::GetOptimalVariableValueLI(int localindex) {
  assert(localindex <= NumberOfVariables && localindex >= 1);
  return Var[localindex - 1]->Optimum;
}

void Problem::SetOptimalVariableValue(int varID, double value) {
  int localindex = VarIDToLocalIndex[varID];
  assert(localindex <= NumberOfVariables && localindex >= 1);
  Var[localindex - 1]->Optimum = value;
}
void Problem::SetOptimalVariableValueLI(int localindex, double value) {
  assert(localindex <= NumberOfVariables && localindex >= 1);
  Var[localindex - 1]->Optimum = value;
}

void Problem::SetVariableLBValue(int varID, double value) {
  int localindex = VarIDToLocalIndex[varID];
  assert(localindex <= NumberOfVariables && localindex >= 1);
  Var[localindex-1]->LB = value;
}

void Problem::SetVariableLBValueLI(int localindex, double value) {
  assert(localindex <= NumberOfVariables && localindex >= 1);
  Var[localindex-1]->LB = value;
}

void Problem::SetVariableUBValue(int varID, double value) {
  int localindex = VarIDToLocalIndex[varID];
  assert(localindex <= NumberOfVariables && localindex >= 1);
  Var[localindex-1]->UB = value;
}

void Problem::SetVariableUBValueLI(int localindex, double value) {
  assert(localindex <= NumberOfVariables && localindex >= 1);
  Var[localindex-1]->UB = value;
}

void Problem::SetConstantVariable(int varID, bool constant) {
  int localindex = VarIDToLocalIndex[varID];
  assert(localindex <= NumberOfVariables && localindex >= 1);
  Var[localindex-1]->IsConstant = constant;
}

 
void Problem::SetConstantVariableLI(int localindex, bool constant) {
  assert(localindex <= NumberOfVariables && localindex >= 1);
  Var[localindex-1]->IsConstant = constant;
}

double Problem::GetCurrentVariableValue(int varID) {
  int localindex = VarIDToLocalIndex[varID];
  assert(localindex <= NumberOfVariables && localindex >= 1);
  return x[localindex - 1];
}
double Problem::GetCurrentVariableValueLI(int localindex) {
  assert(localindex <= NumberOfVariables && localindex >= 1);
  return x[localindex - 1];
}

void Problem::SetCurrentVariableValue(int varID, double value) {
  int localindex = VarIDToLocalIndex[varID];
  assert(localindex <= NumberOfVariables && localindex >= 1);
  x[localindex - 1] = value;
}
void Problem::SetCurrentVariableValueLI(int localindex, double value) {
  assert(localindex <= NumberOfVariables && localindex >= 1);
  x[localindex - 1] = value;
}

bool Problem::TestConstraintsFeasibility(double tolerance, double& disc) {
  bool ret = true;
  double t, l, u;
  for(int i = 1; i <= NumberOfConstraints; i++) {
    // if problematic, try taking away "Constr[i-1]->Function->"
    
    Constr[i-1]->FunctionFET = Constr[i-1]->Function->GetFastEvalTree();
    t = Constr[i-1]->Function->FastEval(Constr[i-1]->FunctionFET, x, 
					VarIDToLocalIndex, NumberOfVariables);
    l = Constr[i - 1]->LB;
    u = Constr[i - 1]->UB;
    if (!(t >= l - tolerance && t <= u + tolerance)) {
#ifdef VERBOSEDEBUG
      std::cerr << "c[" << i << "]: l < " 
		<< Constr[i - 1]->Function->ToString() 
		<< " < u (" << l << " < " << t << " < " 
		<< u << ")" << std::endl;
#endif
      ret = false;
      disc = std::max(t - l + tolerance, u + tolerance - t);
      break;
    }
  }
  return ret;
}

bool Problem::TestConstraintsFeasibility(int cID, 
					 double tolerance, double& disc) {
  bool ret = true;
  double t, l, u;
//  int i = ConstrIDToLocalIndex[cID];
  int i = cID;
  // if problematic, try taking away "Constr[i-1]->Function->"
  //t = FastEval(Constr[i-1]->FunctionFET, x, 
  t = Constr[i-1]->Function->FastEval(Constr[i-1]->FunctionFET, x, 
				      VarIDToLocalIndex, NumberOfVariables);
  l = Constr[i - 1]->LB;
  u = Constr[i - 1]->UB;
  if (!(t >= l - tolerance && t <= u + tolerance)) {
#ifdef VERBOSEDEBUG
    std::cerr << setprecision(24) << "c[" << i << "]: l < " 
	      << Constr[i - 1]->Function->ToString() 
	      << " < u (" << l-tolerance << " < " << t << " < " 
	      << u+tolerance << ")" << std::endl;
#endif
    ret = false;
    disc = std::max(t - l + tolerance, u + tolerance - t);
  }
  return ret;
}

double Problem::TestConstraintsInfeasibilityQuantity(int cID, 
					 double tolerance, double& disc) {
  double t, l, u;
  //int i = ConstrIDToLocalIndex[cID];
  int i = cID;
  // if problematic, try taking away "Constr[i-1]->Function->"
  //t = FastEval(Constr[i-1]->FunctionFET, x, 
  t = Constr[i-1]->Function->FastEval(Constr[i-1]->FunctionFET, x, 
				      VarIDToLocalIndex, NumberOfVariables);
  l = Constr[i - 1]->LB;
  u = Constr[i - 1]->UB;
  int prova;
//  std::cout << l << " <= " << t << " <= " << u << std::endl;
  if (!(t >= l - tolerance)) {
#ifdef VERBOSEDEBUG
    std::cerr << setprecision(24) << "c[" << i << "]: l < " 
	      << Constr[i - 1]->Function->ToString() 
	      << " < u (" << l-tolerance << " < " << t << " < " 
	      << u+tolerance << ")" << std::endl;
#endif
    disc = std::max(t - l + tolerance, u + tolerance - t);
    return - abs(l - t);
  }
  if (!(t <= u + tolerance)) {
#ifdef VERBOSEDEBUG
    std::cerr << setprecision(24) << "c[" << i << "]: l < " 
	      << Constr[i - 1]->Function->ToString() 
	      << " < u (" << l-tolerance << " < " << t << " < " 
	      << u+tolerance << ")" << std::endl;
#endif
    disc = std::max(t - l + tolerance, u + tolerance - t);
    return abs(t - u);
  }

  return 0;
}

bool Problem::TestVariablesFeasibility(double tolerance, double& disc) {
  bool ret = true;
  for(int i = 0; i < NumberOfVariables; i++) {
    if (Var[i]->IsIntegral) {
      disc = fabs(x[i] - rint(x[i]));
      if (disc > tolerance) {
	ret = false;
	break;
      }
    }
    if (x[i] < Var[i]->LB - tolerance || x[i] > Var[i]->UB + tolerance) {
      disc = std::max(Var[i]->LB - x[i], x[i] - Var[i]->UB);
      ret = false;
      break;
    }
  }
  return ret;
}

bool Problem::TestVariablesFeasibility(int vID, 
				       double tolerance, double& disc) {
  bool ret = true;
  int i = VarIDToLocalIndex[vID];
  if (Var[i]->IsIntegral) {
    disc = fabs(x[i] - rint(x[i]));
    if (disc > tolerance) {
      ret = false;
    }
  } else if (x[i] < Var[i]->LB - tolerance || x[i] > Var[i]->UB + tolerance) {
    disc = std::max(Var[i]->LB - x[i], x[i] - Var[i]->UB);
    ret = false;
  }
  return ret;
}

double Problem::GetStartingPoint(int varID) {
  int localindex = VarIDToLocalIndex[varID];
  assert(localindex <= NumberOfVariables && localindex >= 1);
  return StartingPoint[localindex - 1];
}
double Problem::GetStartingPointLI(int localindex) {
  assert(localindex <= NumberOfVariables && localindex >= 1);
  return StartingPoint[localindex - 1];
}

double Problem::GetOptimalObjectiveValue(int objID) {
  int localindex = ObjIDToLocalIndex[objID];
  assert(localindex <= NumberOfObjectives && localindex >= 1);
  return Obj[localindex - 1]->Optimum;
}
void Problem::SetOptimalObjectiveValue(int objID, double value) {
  int localindex = ObjIDToLocalIndex[objID];
  assert(localindex <= NumberOfObjectives && localindex >= 1);
  Obj[localindex - 1]->Optimum = value;
}

void Problem::GetSolution(std::map<int,double>& objfunval, 
			  std::map<int,double>& solution) {
  int nof = GetNumberOfObjectives();
  int nov = GetNumberOfVariables();
  if (IsSolved()) {
    objfunval.erase(objfunval.begin(), objfunval.end());
    for(int i = 1; i <= nof; i++) {
      objfunval[GetObjectiveID(i)] = GetObjectiveLI(i)->Optimum;
    }
    solution.erase(solution.begin(), solution.end());
    for(int i = 1; i <= nov; i++) {
      solution[GetVariableID(i)] = GetVariableLI(i)->Optimum;
    }
  } else if (IsFeasible() == -1) {
    // if infeasible, set everything to INFINITY
    objfunval.erase(objfunval.begin(), objfunval.end());
    for(int i = 1; i <= nof; i++) {
      objfunval[GetObjectiveID(i)] = ROSEINFINITY;
    }
    solution.erase(solution.begin(), solution.end());
    for(int i = 1; i <= nov; i++) {
      solution[GetVariableID(i)] = ROSEINFINITY;
    }
  } else if (IsSolved() == false) {
    // if still unsolved, set everything to NaN
    objfunval.erase(objfunval.begin(), objfunval.end());
    for(int i = 1; i <= nof; i++) {
      objfunval[GetObjectiveID(i)] = 0.0/0.0;
    }
    solution.erase(solution.begin(), solution.end());
    for(int i = 1; i <= nov; i++) {
      solution[GetVariableID(i)] = 0.0/0.0;
    }
  }    
}

void Problem::GetSolutionLI(std::vector<double>& objfunval, 
			    std::vector<double>& solution) {
  int nof = GetNumberOfObjectives();
  int nov = GetNumberOfVariables();
  if (IsSolved()) {
    if (objfunval.size() != nof) {
      objfunval.erase(objfunval.begin(), objfunval.end());
      for(int i = 1; i <= nof; i++) {
	objfunval.push_back(0.0);
      }
    }
    for(int i = 1; i <= nof; i++) {
      objfunval[i - 1] = GetObjectiveLI(i)->Optimum;
    }
    if (solution.size() != nov) {
      solution.erase(solution.begin(), solution.end());
      for(int i = 1; i <= nov; i++) {
	solution.push_back(0.0);
      }
    }
    for(int i = 1; i <= nov; i++) {
      solution[i - 1] = GetVariableLI(i)->Optimum;
    }
  } else if (IsFeasible() == -1) {
    // if infeasible, set everything to INFINITY
    if (objfunval.size() != nof) {
      objfunval.erase(objfunval.begin(), objfunval.end());
      for(int i = 1; i <= nof; i++) {
	objfunval.push_back(ROSEINFINITY);
      }
    }
    if (solution.size() != nov) {
      solution.erase(solution.begin(), solution.end());
      for(int i = 1; i <= nov; i++) {
	solution.push_back(ROSEINFINITY);
      }
    }
  } else if (IsSolved() == false) {
    // if still unsolved, set everything to NaN
    if (objfunval.size() != nof) {
      objfunval.erase(objfunval.begin(), objfunval.end());
      for(int i = 1; i <= nof; i++) {
	objfunval.push_back(0.0/0.0);
      }
    }
    if (solution.size() != nov) {
      solution.erase(solution.begin(), solution.end());
      for(int i = 1; i <= nov; i++) {
	solution.push_back(0.0/0.0);
      }
    }
  }    
}

double Problem::EvalObj(int objID) {
  int localindex = ObjIDToLocalIndex[objID];
  assert(localindex <= NumberOfObjectives && localindex >= 1);

  Obj[localindex - 1]->FunctionFET = Obj[localindex - 1]->Function->GetFastEvalTree();

  double t = FastEval(Obj[localindex - 1]->FunctionFET, x, VarIDToLocalIndex, 
		      NumberOfVariables);
  if (objID == FIRSTOBJ) {
    t += Obj1AdditiveConstant;
  }
  return t;
}

double Problem::EvalNLObj(int objID) {
  int localindex = ObjIDToLocalIndex[objID];
  assert(localindex <= NumberOfObjectives && localindex >= 1);

  Obj[localindex-1]->NonlinearPartFET = Obj[localindex-1]->NonlinearPart->GetFastEvalTree();
  double t = FastEval(Obj[localindex-1]->NonlinearPartFET, x, 
		      VarIDToLocalIndex, NumberOfVariables);
  return t;
}

double Problem::EvalObjDiff(int objID, int varID) {
  int loindex = ObjIDToLocalIndex[objID];
  int lvindex = VarIDToLocalIndex[varID];
  assert(loindex <= NumberOfObjectives && loindex >= 1);
  assert(lvindex <= NumberOfVariables && lvindex >= 1);
  if (!ComputedSymbolicDiff) {
    ComputeSymbolicDiffs();
  }
  double t = FastEval(Obj[loindex - 1]->DiffFET[lvindex - 1], x, 
		      VarIDToLocalIndex, NumberOfVariables);
  return t;
}

double Problem::EvalObjDiffNoConstant(int objID, int varID) {
  int loindex = ObjIDToLocalIndex[objID];
  int lvindex = VarIDToLocalIndex[varID];
  assert(loindex <= NumberOfObjectives && loindex >= 1);
  assert(lvindex <= NumberOfVariables && lvindex >= 1);
  if (!ComputedSymbolicDiff) {
    ComputeSymbolicDiffs();
  }
  double kt = Obj[loindex-1]->Diff[lvindex-1]->GetConstantPart();
  double t = FastEval(Obj[loindex-1]->DiffFET[lvindex-1], x, 
		      VarIDToLocalIndex, NumberOfVariables) - kt; 
  return t;
}

double Problem::EvalObjDiff2(int objID, int varID1, int varID2) {
  int loindex = ObjIDToLocalIndex[objID];
  int lvidx1 = VarIDToLocalIndex[varID1];
  int lvidx2 = VarIDToLocalIndex[varID2];
  assert(loindex <= NumberOfObjectives && loindex >= 1);
  assert(lvidx1 <= NumberOfVariables && lvidx1 >= 1);
  assert(lvidx2 <= NumberOfVariables && lvidx2 >= 1);
  if (!ComputedSymbolicDiff2) {
    ComputeSymbolicDiffs2();
  }
  double t = FastEval(Obj[loindex-1]->Diff2FET[lvidx1-1][lvidx2-1], x, 
		      VarIDToLocalIndex, NumberOfVariables);
  return t;
}

double Problem::EvalConstr(int constrID) {
  int localindex = ConstrIDToLocalIndex[constrID];
  assert(localindex <= NumberOfConstraints && localindex >= 1);

  Constr[localindex - 1]->FunctionFET =   Constr[localindex - 1]->Function->GetFastEvalTree();
  double t = FastEval(Constr[localindex - 1]->FunctionFET, x, 
		      VarIDToLocalIndex, NumberOfVariables);
  return t;
}

double Problem::EvalNLConstr(int constrID) {
  int localindex = ConstrIDToLocalIndex[constrID];
  assert(localindex <= NumberOfConstraints && localindex >= 1);

  Constr[localindex-1]->NonlinearPartFET = Constr[localindex-1]->NonlinearPart->GetFastEvalTree();
  double t = FastEval(Constr[localindex - 1]->NonlinearPartFET, x, 
		      VarIDToLocalIndex, NumberOfVariables);
  return t;
}

double Problem::EvalConstrDiff(int constrID, int varID) {
  int lcindex = ConstrIDToLocalIndex[constrID];
  int lvindex = VarIDToLocalIndex[varID];
  assert(lcindex <= NumberOfConstraints && lcindex >= 1);
  assert(lvindex <= NumberOfVariables && lvindex >= 1);
  if (!ComputedSymbolicDiff) {
    ComputeSymbolicDiffs();
  }
  double t = FastEval(Constr[lcindex - 1]->DiffFET[lvindex - 1], x, 
		      VarIDToLocalIndex, NumberOfVariables);
  return t;
}

double Problem::EvalConstrDiffNoConstant(int constrID, int varID) {
  int lcindex = ConstrIDToLocalIndex[constrID];
  int lvindex = VarIDToLocalIndex[varID];
  assert(lcindex <= NumberOfConstraints && lcindex >= 1);
  assert(lvindex <= NumberOfVariables && lvindex >= 1);
  if (!ComputedSymbolicDiff) {
    ComputeSymbolicDiffs();
  }
  double kt = Constr[lcindex-1]->Diff[lvindex-1]->GetConstantPart();
  double t = FastEval(Constr[lcindex-1]->DiffFET[lvindex-1], x, 
		      VarIDToLocalIndex, NumberOfVariables) - kt;
  return t;
}

double Problem::EvalConstrDiff2(int constrID, int varID1, int varID2){
  int lcindex = ConstrIDToLocalIndex[constrID];
  int lv1 = VarIDToLocalIndex[varID1];
  int lv2 = VarIDToLocalIndex[varID2];
  assert(lcindex <= NumberOfConstraints && lcindex >= 1);
  assert(lv1 <= NumberOfVariables && lv1 >= 1);
  assert(lv2 <= NumberOfVariables && lv2 >= 1);
  if (!ComputedSymbolicDiff2) {
    ComputeSymbolicDiffs2();
  }
  double t = FastEval(Constr[lcindex-1]->Diff2FET[lv1-1][lv2-1], x, 
		      VarIDToLocalIndex, NumberOfVariables);
  return t;
}

bool Problem::IsObjConstant(int objID) {
  int loindex = ObjIDToLocalIndex[objID];
  assert(loindex <= NumberOfObjectives && loindex >= 1);
  if (!ComputedSymbolicDiff) {
    ComputeSymbolicDiffs();
  }
  return Obj[loindex - 1]->Function->IsConstant();
}

bool Problem::IsObjDiffConstant(int objID, int varID) {
  int loindex = ObjIDToLocalIndex[objID];
  int lvindex = VarIDToLocalIndex[varID];
  assert(loindex <= NumberOfObjectives && loindex >= 1);
  assert(lvindex <= NumberOfVariables && lvindex >= 1);
  if (!ComputedSymbolicDiff) {
    ComputeSymbolicDiffs();
  }
  return Obj[loindex - 1]->Diff[lvindex - 1]->IsConstant();
}

bool Problem::IsObjDiff2Constant(int objID, int varID1, int varID2) {
  int loindex = ObjIDToLocalIndex[objID];
  int lv1 = VarIDToLocalIndex[varID1];
  int lv2 = VarIDToLocalIndex[varID2];
  assert(loindex <= NumberOfObjectives && loindex >= 1);
  assert(lv1 <= NumberOfVariables && lv1 >= 1);
  assert(lv2 <= NumberOfVariables && lv2 >= 1);
  if (!ComputedSymbolicDiff) {
    ComputeSymbolicDiffs();
  }
  return Obj[loindex - 1]->Diff2[lv1-1][lv2-1]->IsConstant();
}

bool Problem::IsConstrConstant(int constrID) {
  int lcindex = ConstrIDToLocalIndex[constrID];
  assert(lcindex <= NumberOfConstraints && lcindex >= 1);
  if (!ComputedSymbolicDiff) {
    ComputeSymbolicDiffs();
  }
  return Constr[lcindex - 1]->Function->IsConstant();
}

bool Problem::IsConstrDiffConstant(int constrID, int varID) {
  int lcindex = ConstrIDToLocalIndex[constrID];
  int lvindex = VarIDToLocalIndex[varID];
  assert(lcindex <= NumberOfConstraints && lcindex >= 1);
  assert(lvindex <= NumberOfVariables && lvindex >= 1);
  if (!ComputedSymbolicDiff) {
    ComputeSymbolicDiffs();
  }
  return Constr[lcindex - 1]->Diff[lvindex - 1]->IsConstant();
}

bool Problem::IsConstrDiff2(int constrID, int varID1, int varID2) {
  int lcindex = ConstrIDToLocalIndex[constrID];
  int lv1 = VarIDToLocalIndex[varID1];
  int lv2 = VarIDToLocalIndex[varID2];
  assert(lcindex <= NumberOfConstraints && lcindex >= 1);
  assert(lv1 <= NumberOfVariables && lv1 >= 1);
  assert(lv2 <= NumberOfVariables && lv2 >= 1);
  if (!ComputedSymbolicDiff) {
    ComputeSymbolicDiffs();
  }
  return Constr[lcindex - 1]->Diff2[lv1-1][lv2-1]->IsConstant();
}

bool Problem::IsConstrActive(int cID, double tolerance, int& UpperLower){
  int i = ConstrIDToLocalIndex[cID];
  assert(i <= NumberOfConstraints && i >= 1);
  bool ret = false;
  double discrepancy;
  if(TestConstraintsFeasibility(i,tolerance,discrepancy))
  {
    // if problematic, try taking away "Constr[i-1]->Function->"
    double t = Constr[i-1]->Function->FastEval(Constr[i-1]->FunctionFET, x, 
					     VarIDToLocalIndex, 
					     NumberOfVariables);
    double l = Constr[i - 1]->LB;
    double u = Constr[i - 1]->UB;
    if (u - t < tolerance) {
      UpperLower = 1;
      ret = true;
    } 
    if (t - l < tolerance) {
      UpperLower = -1;
      ret = true;
    }
    if (l <= t && t <= u && u-l < tolerance) {
      UpperLower = 0;
      ret = true;
    }
    if (!ret) {
      UpperLower = 0;
    }
  }
  return ret;
}

void Problem::ComputeSymbolicDiffs(void) {
  if (!ComputedSymbolicDiff) {
    ComputedSymbolicDiff = true;
    for(int i = 0; i < NumberOfObjectives; i++) {
      for(int j = 1; j <= NumberOfVariables; j++) {
	Obj[i]->Diff.push_back(Diff(Obj[i]->Function, GetVariableID(j)));
	Obj[i]->DiffFET.push_back(Obj[i]->Diff[j - 1]->GetFastEvalTree());
      }
    }
    for(int i = 0; i < NumberOfConstraints; i++) {
//std::cout << "qui -i= " << i << std::endl;
      for(int j = 1; j <= NumberOfVariables; j++) {
	Constr[i]->Diff.push_back(Diff(Constr[i]->Function, GetVariableID(j)));
	Constr[i]->DiffFET.push_back(Constr[i]->Diff[j-1]->GetFastEvalTree());
      }
    }
  }
}

void Problem::ComputeSymbolicDiffs2(void) {
  if (!ComputedSymbolicDiff2) {
    ComputedSymbolicDiff2 = true;
    for(int i = 0; i < NumberOfObjectives; i++) {
      for(int j = 0; j < NumberOfVariables; j++) {
	std::vector<Expression> tmpve;
	std::vector<FastEvalTree*> tmpvf;
	Obj[i]->Diff.push_back(Diff(Obj[i]->Function, GetVariableID(j+1)));
	for(int k = 1; k <= NumberOfVariables; k++) {
	  tmpve.push_back(Diff(Obj[i]->Diff[j], GetVariableID(k)));
	  tmpvf.push_back(tmpve[k - 1]->GetFastEvalTree());
	}
	Obj[i]->Diff2.push_back(tmpve);
	Obj[i]->Diff2FET.push_back(tmpvf);
      }
    }
    for(int i = 0; i < NumberOfConstraints; i++) {
      for(int j = 0; j < NumberOfVariables; j++) {
	std::vector<Expression> tmpve;
	std::vector<FastEvalTree*> tmpvf;
	for(int k = 1; k <= NumberOfVariables; k++) {
	  tmpve.push_back(Diff(Constr[i]->Diff[j], GetVariableID(k)));
	  tmpvf.push_back(tmpve[k - 1]->GetFastEvalTree());
	}
	Constr[i]->Diff2.push_back(tmpve);
	Constr[i]->Diff2FET.push_back(tmpvf);
      }
    }
  }
}

void Problem::DebugSymbolicDiffs(void) {
  if (!ComputedSymbolicDiff) {
    ComputeSymbolicDiffs();
  }
  if (!ComputedSymbolicDiff2) {
    ComputeSymbolicDiffs2();
  }
  std::cerr << "* Derivatives debug point (only nonzeroes displayed)\n";
  std::cerr << "** First order derivatives:\n";
  for(int i = 0; i < NumberOfObjectives; i++) {
    std::cerr << "**** Objective " << i + 1 << ":\n";
    for(int j = 0; j < NumberOfVariables; j++) {
      if (!Obj[i]->Diff[j]->IsZero()) {
	std::cerr << "      (d/dx_" << GetVariableID(j+1) << ")("
	     << Obj[i]->Function->ToString() << ") = "
	     << Obj[i]->Diff[j]->ToString() << std::endl;
      }
    }
  }
  for(int i = 0; i < NumberOfConstraints; i++) {
    std::cerr << "**** Constraint " << i + 1 << ":\n";
    for(int j = 0; j < NumberOfVariables; j++) {
      if (!Constr[i]->Diff[j]->IsZero()) {
	std::cerr << "      (d/dx_" << GetVariableID(j+1) << ")("
	     << Constr[i]->Function->ToString() << ") = "
	     << Constr[i]->Diff[j]->ToString() << std::endl;
      }
    }
  }
  std::cerr << "** Second order derivatives:\n";
  for(int i = 0; i < NumberOfObjectives; i++) {
    std::cerr << "**** Objective " << i + 1 << ":\n";
    for(int j = 0; j < NumberOfVariables; j++) {
      for(int k = 0; k < NumberOfVariables; k++) {
	if (!Obj[i]->Diff2[j][k]->IsZero()) {
	  std::cerr << "      (d^2/dx_" << GetVariableID(j+1) 
	       << "x_" << GetVariableID(k+1) << ")("
	       << Obj[i]->Function->ToString() << ") = "
	       << Obj[i]->Diff2[j][k]->ToString() << std::endl;
	}
      }
    }
  }
  for(int i = 0; i < NumberOfConstraints; i++) {
    std::cerr << "**** Constraint " << i + 1 << ":\n";
    for(int j = 0; j < NumberOfVariables; j++) {
      for(int k = 0; k < NumberOfVariables; k++) {
	if (!Constr[i]->Diff2[j][k]->IsZero()) {
	  std::cerr << "      (d^2/dx_" << GetVariableID(j+1) 
	       << "x_" << GetVariableID(k+1) << ")("
	       << Constr[i]->Function->ToString() << ") = "
	       << Constr[i]->Diff2[j][k]->ToString() << std::endl;
	}
      }
    }
  }
  std::cerr << "** End derivatives\n";
  std::cerr << "* End derivatives debug point\n";
}

void Problem::DeleteVariable(int varID) {
  std::vector<Variable*>::iterator vi = 
    find_if(Var.begin(), Var.end(), IsIDEqual<Variable*>(varID));
  if (vi != Var.end()) {
    Var.erase(vi);
  }
  NumberOfVariables--;

  for (int i = 1; i <= NumberOfVariables; i++) {
    VarIDToLocalIndex[Var[i-1]->ID] = i;
  }
}

void Problem::DeleteObjective(int objID) {
  std::vector<Objective*>::iterator oi = 
    find_if(Obj.begin(), Obj.end(), IsIDEqual<Objective*>(objID));
  if (oi != Obj.end()) {
    Obj.erase(oi);
  }
  NumberOfObjectives--;
  for (int i = objID; i <= NumberOfObjectives; i++) {
      Obj[i-1]->NonlinearPart.SetTo(Obj[i-1]->Function->GetNonlinearPart());
      Obj[i-1]->FunctionFET = Obj[i-1]->Function->GetFastEvalTree();
      Obj[i-1]->NonlinearPartFET = Obj[i-1]->NonlinearPart->GetFastEvalTree();
  }

  for (int i = 1; i <= NumberOfObjectives; i++) {
    ObjIDToLocalIndex[Obj[i-1]->ID] = i;
  }
}

void Problem::DeleteConstraint(int constrID) {
  std::vector<Constraint*>::iterator ci = 
    find_if(Constr.begin(), Constr.end(), IsIDEqual<Constraint*>(constrID));
  if (ci != Constr.end()) {
    Constr.erase(ci);
  }
  NumberOfConstraints--;
  for (int i = constrID; i <= NumberOfConstraints; i++) {
      Constr[i-1]->NonlinearPart.SetTo(Constr[i-1]->Function->GetNonlinearPart());
      DeleteFastEvalTree(Constr[i-1]->FunctionFET);
      DeleteFastEvalTree(Constr[i-1]->NonlinearPartFET);
      Constr[i-1]->FunctionFET = Constr[i-1]->Function->GetFastEvalTree();
      Constr[i-1]->NonlinearPartFET = Constr[i-1]->NonlinearPart->GetFastEvalTree();
  }

  int i = 1;
  for(ci = Constr.begin(); ci != Constr.end(); ci++) {
    ConstrIDToLocalIndex[(*ci)->ID] = i;
    i++;
  }
}

void Problem::Parse(char* thefilename) {
  FILE* Ev3f;
  char buffer[PARSEBUFFER * 2];
  char element[PARSEBUFFER * 2];

  char* t1;
  char* t2;
  char* buf;

  bool getbufferfromfileflag = true;
  bool changesectionflag = false;
  bool bufferishalffullflag = false;

  if ((Ev3f = fopen(thefilename, "r")) == NULL) {
    std::cerr << "rose::Parse: can't open " << thefilename << "\n";
    exit(128);
  }

  // initialization
  Name = thefilename;
  NumberOfVariables = 0;
  NumberOfConstraints = 0;
  int status = NoSectionType;

  int lineno = 0;

  while(1) {
  top:
    if (getbufferfromfileflag) {
      // get line from file
      if (bufferishalffullflag) {
	if (strlen(buffer) > 0) {
	  buf = buffer + strlen(buffer) - 1;
	} else {
	  buf = buffer;
	}
      } else {
	buf = buffer;
      }
      lineno++;
      fgets(buf, PARSEBUFFER, Ev3f);
      if (feof(Ev3f)) {
	break;
      }
      if (bufferishalffullflag) {
	bufferishalffullflag = false;
	while (*buf == ' ' || *buf == '\t') {
	  *buf = ' ';
	  buf++;
	}
	buf = buffer;
      }
    }
    // if not an empty line or a comment
    if (strlen(buf) > 1 && buf[0] != '#') {
      // skip spaces
      while (*buf == ' ' || *buf == '\t') buf++;
      // if status is undefined, find out
      if (status == NoSectionType) {
	if (strstr(buf, "variables") != NULL) {
	  // variable section
	  status = VariableSectionType;
	} else if (strstr(buf, "objfun") != NULL) {
	  // objective function section
	  status = ObjfunSectionType;
	} else if (strstr(buf, "constraints") != NULL) {
	  // constraints section
	  status = ConstraintSectionType;
	} else if (strstr(buf, "startingpoint") != NULL) {
	  // starting point section
	  status = StartingpointSectionType;
	} else if (strstr(buf, "parameters") != NULL || 
		   strstr(buf, "options") != NULL) {
	  // options section
	  status = ParameterSectionType;
	} else {
	  // file is invalid
	  std::cerr << "rose::Parse: file does not seem to be in Rose "
	       << "format\n";
	  exit(129);
	}
      }
      if (status != NoSectionType) {
	// look for the "=" sign
	t1 = strchr(buf, '=');
	// if not found, set at beginning of buf
	if (t1 == NULL) {
	  t1 = buf;
	} else {
	  // found, augment by one
	  t1++;
	}
	// skip spaces
	while(*t1 == ' ' || *t1 == '\t') t1++;
	// look for a comma
	t2 = strchr(t1, ',');
	if (t2 == NULL) {
	  // no comma, look for a semicolon
	  t2 = strchr(t1, ';');
	  if (t2 == NULL) {
	    // no semicolon, look for a newline
	    t2 = strchr(t1, '\n');
	    if (t2 == NULL) {
	      // file not valid
	      std::cerr << "rose::Parse: newline not found\n";
	      assert(false);
	    } else {
	      // found a newline, go back to top of while to read new line
	      getbufferfromfileflag = true;
	      changesectionflag = false;
	      bufferishalffullflag = true;
	      goto top;
	    }
	  } else {
	    // found a semicolon, go on to next section
	    changesectionflag = true;
	    getbufferfromfileflag = true; // after a semicolon, change line
	  }
	} else {
	  // comma found, see if next info is on this line or next line
	  buf = t2 + 1;
	  while (*buf == ' ' || *buf == '\t') buf++;
	  if (*buf == '\n') {
	    // next info is on next line
	    getbufferfromfileflag = true;
	  } else {
	    // next info is on the same line after the comma
	    getbufferfromfileflag = false;
	  }
	}
	// terminate buffer at comma or semicolon
	*t2 = '\0';
	// copy the interesting bit into the element array
	strncpy(element, t1, t2 - t1 + 1);
	element[strlen(t1) + 1] = '\0';
	ParseEv3String(element, status, lineno);
	if (changesectionflag) {
	  changesectionflag = false;
	  status = NoSectionType;
	}
      }
    }
  }
  fclose(Ev3f);

  // size up the problem
  int i, j;

  // variables
  if (NumberOfVariables <= 0) {
    std::cerr << "rose::Parse: number of variables must be positive\n";
    assert(false);
  }

  if (NumberOfVariables != (int) StartingPoint.size()) {
    std::cerr << "rose::Parse: number of vars != components "
	 << "of starting point\n";
    assert(false);
  }

  x = new double [NumberOfVariables];

  // simplify the problem if possible
  Simplifier(false);
  
  // set formulation information 
  SetProblemLinear();
}

void Problem::SetProblemLinear(void) {
  // is problem linear (MILP)? 
  IsProbLinear = true;
  for(int i = 1; i <= NumberOfObjectives; i++) {
    if (!GetObjectiveLI(i)->Function->IsLinear()) {
      IsProbLinear = false;
      cerr << "WARNING: problem is nonlinear" << endl;
      break;
    }
  }
  if (IsProbLinear) {
    for(int i = 1; i <= NumberOfConstraints; i++) {
      if (!GetConstraintLI(i)->Function->IsLinear()) {
	IsProbLinear = false;
	cerr << "WARNING: problem is nonlinear" << endl;
	break;
      }
    }
  }
}

void Problem::ParseEv3String(char* element, int status, int lineno) {

  // temporary character space
  char* t1;
  char* t2;
  char* t3;

  // parameters
  std::string pn;
  char *ep;
  char backupchar;
  bool boolval;
  int intval;
  double doubleval;
  std::string stringval;
  int paramtype;

  // etc
  int i;
  int nerr;
  double tmpvlb, tmpvub;

  // parse the blasted string
  switch(status) {
  case VariableSectionType:
    {
      // create a new variable
      Variable* theVar = new Variable;
      theVar->Optimum = 0;
      theVar->Persistent = false;
      theVar->IsIntegral = false;
      theVar->ID = NewVariableID();
      // read in lower bound
      t1 = strchr(element, '<');
      if (t1 == NULL) {
	std::cerr << "rose::ParseEv3String(Var): '<' not found on line " 
	     << lineno << "\n";
	assert(false);
      }
      *t1 = '\0';
      t1++;
      if (strstr(element, "MinusInfinity") != NULL) {
	tmpvlb = -ROSEINFINITY;
      } else {
	sscanf(element, "%lf", &tmpvlb);
      }
      theVar->LB = tmpvlb;
      // read in variable name
      t2 = strchr(t1, '<');
      if (t2 == NULL) {
	std::cerr << "rose::ParseEv3String(Var): '<' not found on line " << lineno
	     << "\n";
	assert(false);
      }
      // skip white space from end of t1 token
      t3 = t2 - 1;
      while(*t3 == ' ' || *t3 == '\t') t3--;
      t3++;
      *t3 = '\0';
      t2++;
      // skip white space at the beginning of t1 token
      while(*t1 == ' ' || *t1 == '\t') t1++;
      theVar->Name = t1;
      // read in upper bound
      t3 = strchr(t2, '/');
      if (t3 == NULL) {
      // rest of string is upper bound
	if (strstr(t2, "PlusInfinity") != NULL) {
	  tmpvub = ROSEINFINITY;
	} else {
	  sscanf(t2, "%lf", &tmpvub);
	}
	theVar->UB = tmpvub;
	theVar->IsIntegral = false;	
      } else {
	// there is also a "variable type" field after the '/' char
	*t3 = '\0';
	t3++;
	if (strstr(t2, "PlusInfinity") != NULL) {
	  tmpvub = ROSEINFINITY;
	} else {
	  sscanf(t2, "%lf", &tmpvub);
	}
	theVar->UB = tmpvub;
	// read in variable type
	if (strstr(t3, "Int") != NULL || strstr(t3, "int") != NULL) {
	  // integer variable
	  theVar->IsIntegral = true;
	  IsProbContinuous = false;
	} else {
	  theVar->IsIntegral = false;
	}
      }
      // increase the variable (local) counter
      NumberOfVariables++;
      VarIDToLocalIndex[theVar->ID] = NumberOfVariables;
      Var.push_back(theVar);
    }
    break;
  case ObjfunSectionType:
    {
      Objective* theObj = new Objective;
      theObj->ID = NewObjectiveID();
      // configure expression parser
      for(i = 1; i <= NumberOfVariables; i++) {
	exprparser.SetVariableID(Var[i-1]->Name, i);
      }
      // strip leading and closing square brackets
      t1 = strchr(element, '[');
      if (t1 == NULL) {
	std::cerr << "rose::ParseEv3String(Obj): '[' not found on line " << lineno
	     << "\n";
      }
      // look for min or max
      *t1 = '\0';
      if (strstr(element, "max") != NULL) {
	theObj->OptDir = Maximization;
      } else {
	theObj->OptDir = Minimization;
      }
      *t1 = '[';
      // start parsing
      t1++;
      while(*t1 == ' ' || *t1 == '\t' || *t1 == '|') t1++;
      t2 = strchr(t1, ']');
      if (t2 == NULL) {
	std::cerr << "rose::ParseEv3String(Obj): ']' not found on line " << lineno
	     << "\n";
      }
      t2--;
      while(*t2 == ' ' || *t2 == '\t') t2--;
      t2++;
      *t2 = '\0';    
      theObj->Name = t1;
      theObj->Function.SetTo(exprparser.Parse(t1, nerr));
      theObj->NonlinearPart.SetTo(theObj->Function->GetNonlinearPart());
      theObj->FunctionFET = NULL;
      theObj->NonlinearPartFET = NULL;
      theObj->OptDir = Minimization;
      NumberOfObjectives++;
      ObjIDToLocalIndex[theObj->ID] = NumberOfObjectives;
      Obj.push_back(theObj);
    }
    break;
  case ConstraintSectionType:
    {
      bool arethereconstrs = true;
      if ((t1 = strchr(element, '0')) != NULL) {
	t1++;
	while(*t1 == ' ' || *t1 == '\t') t1++;
	if (*t1 == '\0') {
	  arethereconstrs = false;
	}
      }
      if ((t1 = strchr(element, '[')) != NULL) {
	t1++;
	while(*t1 == ' ' || *t1 == '\t') t1++;
	if (*t1 == ']') {
	  arethereconstrs = false;
	} 
      }
      if (arethereconstrs) {
	Constraint* theConstr = new Constraint;
	theConstr->ID = NewConstraintID();
	// strip leading and closing square brackets
	t1 = strchr(element, '[');
	if (t1 == NULL) {
	  std::cerr << "rose::ParseEv3String(Constr): '[' not found on line " 
	       << lineno << "\n";
	  assert(false);
	}
	t1++;
	while(*t1 == ' ' || *t1 == '\t') t1++;
	t2 = strchr(t1, ']');
	if (t2 == NULL) {
	  std::cerr << "rose::ParseEv3String(Constr): ']' not found on line " 
	       << lineno << "\n";
	  assert(false);
	}
	t2--;
	while(*t2 == ' ' || *t2 == '\t') t2--;
	t2++;
	*t2 = '\0';    
	element = t1;
	// read in lower bound
	t1 = strchr(element, '<');
	if (t1 == NULL) {
	  std::cerr << "rose::ParseEv3String(Constr): '<' not found on line " 
	       << lineno << "\n";
	  assert(false);
	}
	*t1 = '\0';
	t1++;
	if (strstr(element, "MinusInfinity") != NULL) {
	  tmpvlb = -ROSEINFINITY;
	} else {
	  sscanf(element, "%lf", &tmpvlb);
	}
	theConstr->LB = tmpvlb;
	// read in upper bound
	t2 = strchr(t1, '<');
	if (t2 == NULL) {
	  std::cerr << "rose::ParseEv3String(Constr): '<' not found on line " 
	       << lineno << "\n";
	  assert(false);
	}
	*t2 = '\0';
	t2++;
	if (strstr(t2, "PlusInfinity") != NULL) {
	  tmpvub = ROSEINFINITY;
	} else {
	  sscanf(t2, "%lf", &tmpvub);
	}
	theConstr->UB = tmpvub;
	// read in constraint expression
	while (*t1 == '|' || *t1 == ' ' || *t1 == '\t') {
	  t1++;
	}
	theConstr->Name = t1;
	theConstr->Function.SetTo(exprparser.Parse(t1, nerr));
	theConstr->FunctionFET = NULL;
	theConstr->NonlinearPartFET = NULL;
	theConstr->NonlinearPart.SetTo(theConstr->Function->
				       GetNonlinearPart());
	NumberOfConstraints++;
	ConstrIDToLocalIndex[theConstr->ID] = NumberOfConstraints;
	Constr.push_back(theConstr);
      }
    }
    break;
  case StartingpointSectionType:
    sscanf(element, "%lf", &tmpvlb);
    StartingPoint.push_back(tmpvlb);
    break;
  case ParameterSectionType:
    t1 = element;
    // skip blanks
    while(*t1 == ' ' || *t1 == '\t') {
      t1++;
    }
    // look for param name
    t2 = t1;
    while(*t2 != ' ' && *t2 != '\t' && *t2 != '\n') {
      t2++;
    }
    backupchar = *t2;
    *t2 = '\0';
    pn = t1;
    *t2 = backupchar;
    t1 = t2;
    // skip blanks
    while(*t1 == ' ' || *t1 == '\t') {
      t1++;
    }
    // now value from t1 to end of element is param value
    paramtype = parmtype(t1, boolval, intval, doubleval, stringval);
    if (pn == "NoSimplifier" && boolval == true) {
      NoSimplifier = true;
    }
    if (paramtype == BoolType) {
      ParameterBlob.SetBoolParameter(pn, boolval);
    } else if (paramtype == IntType) {
      ParameterBlob.SetIntParameter(pn, intval);
    } else if (paramtype == DoubleType) {
      ParameterBlob.SetDoubleParameter(pn, doubleval);
    } else if (paramtype == StringType) {
      ParameterBlob.SetStringParameter(pn, stringval);
    }
    break;
  default:
    break;
  }
}

Ampl::ASL* Problem::ParseAMPL(char** argv, int argc) {

  using namespace Ampl;
  using namespace std;

  // AMPL stuff
  FILE *nl;
  char *stub, *getenvErr;
  ASL *asl;
  cgrad *cg;
  ograd *og;
  efunc *r_ops_int[N_OPS];
  asl = ASL_alloc(ASL_read_fg);
  stub = getstub(&argv, TheOInfo);
  nl = jac0dim(stub, (fint)strlen(stub));
  Uvx = (real *)Malloc(n_var*sizeof(real));
  Urhsx = (real *)Malloc(n_con*sizeof(real));
  A_vals = (real *)Malloc(nzc*sizeof(real));   
  for(int i = 0; i < N_OPS; i++) {
    r_ops_int[i] = (efunc*)(unsigned long)i;
  }
  R_OPS = r_ops_int;
  want_derivs = 0;
  want_xpi0 = 3;
  X0 = (real*) Malloc(n_var*sizeof(real));
  for(int i = 0; i < n_var; i++) {
    X0[i] = 0;
  }
  fg_read(nl,0);
  R_OPS = 0;
  int number_nonlinear_vars = std::max(nlvc, nlvo);
  int number_linear_arcs = nwv;
  int number_linear_vars = n_var - number_nonlinear_vars - 
    niv - nbv - nwv;
  int number_binary_vars = nbv;
  int number_integer_vars = niv;
  int count_nonlinear_vars = number_nonlinear_vars;
  int count_linear_arcs = count_nonlinear_vars + number_linear_arcs;
  int count_linear_vars = count_linear_arcs + number_linear_vars;
  int count_binary_vars = count_linear_vars + number_binary_vars;
  int count_integer_vars = count_binary_vars + number_integer_vars;
  int count_continuous_vars = count_linear_vars;

#ifdef VERBOSEDEBUG
  std::cout << "number_nonlinear_vars: " << number_nonlinear_vars << std::endl;
  std::cout << "number_linear_arcs: " << number_linear_arcs << std::endl;
  std::cout << "number_linear_vars = n_var-number_nonlinear_vars-niv-nbv-nwv = " 
       << n_var << "-" << number_nonlinear_vars << "-" << niv << "-" << nbv 
       << "-" << nwv << " = " << number_linear_vars << std::endl;
  std::cout << "number_binary_vars: " << nbv << std::endl;
  std::cout << "number_integer_vars: " << niv << std::endl;
  std::cout << "continuous variables: x1, ..., x" << count_linear_vars << std::endl;
  std::cout << "binary variables: x" << count_linear_vars + 1 << ", ..., x"
       << count_binary_vars << std::endl;
  std::cout << "integer variables: x" << count_binary_vars + 1 << ", ..., x"
       << count_integer_vars << std::endl;
  assert (count_integer_vars == n_var);
  std::cout << "\n";
#endif

  int ret = 0;
  bool Quiet = false;

  // parse the problem
  ASL_fg* theaslfg = (ASL_fg*) asl;
  // variables
  char vnamebuf[16];
  strcpy(vnamebuf, "x");
  char* vthename = vnamebuf + 1;
  for(int i = 0; i < n_var; i++) {
    Variable *theVar = new Variable;
    theVar->Optimum = 0;
    theVar->Persistent = false;
    theVar->ID = NewVariableID();
    if (LUv[i] < -ROSEINFINITY) {
      theVar->LB = -ROSEINFINITY;
    } else {
      theVar->LB = LUv[i];
    }
    sprintf(vthename, "%d", theVar->ID);
    theVar->Name = vnamebuf;
    if (Uvx[i] > ROSEINFINITY) {
      theVar->UB = ROSEINFINITY;
    } else {
      theVar->UB = Uvx[i];
    }
    // setting the integrality flag of variables
    if (i < nlvb - nlvbi)  {
      theVar->IsIntegral = false;
    } else if (i < nlvb) {
      theVar->IsIntegral = true;
    } else if (i < nlvc - nlvci) {
      theVar->IsIntegral = false;
    } else if (i < nlvc) {
      theVar->IsIntegral = true;
    } else if (i < nlvo - nlvoi) {
      theVar->IsIntegral = false;
    } else if (i < nlvo) {
      theVar->IsIntegral = true;
    } else if (i < count_linear_arcs) {
      theVar->IsIntegral = false;
    } else if (i < count_linear_vars) {
      theVar->IsIntegral = false;
    } else if (i < n_var - number_integer_vars) {
      theVar->IsIntegral = true;
      theVar->LB = 0;
      theVar->UB = 1;
    } else {
      theVar->IsIntegral = true;
    }
    if (theVar->IsIntegral) {
      IsProbContinuous = false;
    }
    NumberOfVariables++;
    VarIDToLocalIndex[theVar->ID] = NumberOfVariables;
    Var.push_back(theVar);
#ifdef VERBOSEDEBUG
    std::cout << "variable " << i + 1 << ": " << theVar->LB << " <= " 
	 << theVar->Name << " <= " << theVar->UB;
    if (theVar->IsIntegral) {
      std::cout << ", integer";
    }
    std::cout << std::endl;
#endif
  }

  // objective functions
  char onamebuf[16];
  char* othename = onamebuf + 1;
  strcpy(onamebuf, "o");
  for(int i = 0; i < n_obj; i++) {
    Objective* theObj = new Objective;
    theObj->ID = NewObjectiveID();
    int minmax = (int) (*(theaslfg->i.objtype_));
    if (minmax == 0) {
      theObj->OptDir = Minimization;
    } else {
      theObj->OptDir = Maximization;
    }
    sprintf(othename, "%d", theObj->ID);
    theObj->Name = onamebuf;
    // linear part
    for(og = Ograd[i]; og; og = og->next) {
      if (fabs(og->coef) > EPSILONTOLERANCE) {
	sprintf(vthename, "%d", og->varno + 1);
	std::string thevn = vnamebuf;
	Expression thevar(og->coef, og->varno + 1, thevn);
	theObj->Function = SumLink(theObj->Function, thevar);
      }
    }
    // nonlinear part
    Ampl::expr* theamplexpr = ((OBJ_DE)+i)->e;
    if (!IsAMPLExpressionZero(theamplexpr)) {
      theObj->Function = SumLink(theObj->Function, 
				 ParseAMPLExpressionTree(theamplexpr));
    }
    theObj->NonlinearPart.SetTo(theObj->Function->GetNonlinearPart());
    theObj->FunctionFET = NULL;
    theObj->NonlinearPartFET = NULL; 
    NumberOfObjectives++;
    ObjIDToLocalIndex[theObj->ID] = NumberOfObjectives;
    Obj.push_back(theObj);
#ifdef VERBOSEDEBUG
    std::cout << "objective: " << theObj->Function->ToString() << std::endl;
#endif
  }

  // constraints
  char cnamebuf[16];
  char* cthename = cnamebuf + 1;
  strcpy(cnamebuf, "c");
  if (n_con > 0) {
    // read linear jacobian in sparse form
    std::map<std::pair<int,int>,double> A;
    std::map<std::pair<int,int>,double>::iterator mit;
    std::pair<int,int> p;
    for(int j = 0; j < n_var; j++) {
      p.second = j;
      for(int i = A_colstarts[j]; i < A_colstarts[j + 1]; i++) {
	p.first = A_rownos[i];
	if (fabs(A_vals[i]) < EPSILONTOLERANCE) {
	  A_vals[i] = 0;
	}
	A[p] = A_vals[i];
      }
    }
    for(int i = 0; i < n_con; i++) {
      Constraint* theConstr = new Constraint;
      theConstr->ID = NewConstraintID();
      theConstr->LB = LUrhs[i];
      theConstr->UB = Urhsx[i];
      sprintf(cthename, "%d", theConstr->ID);
      theConstr->Name = cnamebuf;
      p.first = i;
      for(int j = 0; j < n_var; j++) {
	p.second = j;
	mit = A.find(p);
	if (mit != A.end()) {
	  if (fabs(mit->second) > EPSILONTOLERANCE) {
	    sprintf(vthename, "%d", j + 1);
	    std::string thevn = vnamebuf;
	    Expression thevar(mit->second, j + 1, thevn);
	    theConstr->Function = SumLink(theConstr->Function, thevar);
	  }
	}
      }
      Ampl::expr* theamplexpr = ((CON_DE)+i)->e;
      if (!IsAMPLExpressionZero(theamplexpr)) {
	theConstr->Function = SumLink(theConstr->Function, 
				      ParseAMPLExpressionTree(theamplexpr));
      }
      theConstr->NonlinearPart.SetTo(theConstr->Function->GetNonlinearPart());
      theConstr->FunctionFET = NULL; 
      theConstr->NonlinearPartFET = NULL; 
      NumberOfConstraints++;
      ConstrIDToLocalIndex[theConstr->ID] = NumberOfConstraints;
      Constr.push_back(theConstr);
#ifdef VERBOSEDEBUG
      std::cout << "constraint " << i + 1 << ": " 
	   << theConstr->Function->ToString() << std::endl;
#endif 
    }
  }

  // starting point
  for(int i = 0; i < n_var; i++) {
    if (X0) {
      StartingPoint.push_back(X0[i]);
    } else {
      StartingPoint.push_back(0.0);
    }
  }

  // AMPL parameters
  // to add an AMPL option, please also change rose.cxx, problem.h
  getopts(argv, Ampl::TheOInfo);
  // main solver
  std::string mainsolvername;
  switch(*Ampl::TheMainSolver) {
  case 0: mainsolvername = "snopt"; break;
  case 1: mainsolvername = "tabu"; break;
  case 2: mainsolvername = "glpk"; break;
  case 3: mainsolvername = "vns"; break;
  case 4: mainsolvername = "gomory"; break;
  case 5: mainsolvername = "limitedbranch"; break;
  case 6: mainsolvername = "localbranch"; break;
  case 7: mainsolvername = "analyser"; break;
  case 8: mainsolvername = "prodbincont"; break;
  case 9: mainsolvername = "print"; break;
  case 10: mainsolvername = "smith"; break;
  case 11: mainsolvername = "copy"; break;
  case 12: mainsolvername = "disaggr"; break;
  case 13: mainsolvername = "oa"; break;
  case 14: mainsolvername = "symmgroup"; break;
  case 15: mainsolvername = "printmod"; break;
  case 16: mainsolvername = "ipopt"; break;
  case 17: mainsolvername = "rvinci"; break;
  case 18: mainsolvername = "rporta"; break;
  case 19: mainsolvername = "rquarticconvex"; break;
  case 20: mainsolvername = "rconvexifier"; break;
  case 21: mainsolvername = "rprintdat"; break;
  case 22: mainsolvername = "rcdd"; break;
  case 23: mainsolvername = "rrelaxation"; break;
  case 24: mainsolvername = "subgradient"; break;
  case 25: mainsolvername = "rconvexifiermod"; break;
  case 26: mainsolvername = "rprintncvxdiscrepancy"; break;
  case 27: mainsolvername = "rfbbtfp"; break;
  case 28: mainsolvername = "rmilp2gph"; break;
  default: mainsolvername = "print"; break;
  }
  ParameterBlob.SetStringParameter("MainSolver", mainsolvername);
  // quiet
  if (*Ampl::TheQuiet == 0) {
    Quiet = false;
  } else {
    Quiet = true;
  }  
  ParameterBlob.SetBoolParameter("Quiet", Quiet);
  // relaxinteger
  bool relaxed;
  if (*Ampl::TheRelaxInteger == 0) {
    relaxed = false;
  } else {
    relaxed = true;
  }
  ParameterBlob.SetBoolParameter("RelaxInteger", relaxed);
  // no automatic simplification
  if (*Ampl::TheNoSimplifier == 0) {
    NoSimplifier = false;
  } else {
    NoSimplifier = true;
  }
  ParameterBlob.SetBoolParameter("NoSimplifier", NoSimplifier);
  // glpk parameters
  std::string glpksolution;
  switch(*Ampl::TheGLPKSolutionMethod) {
  case 0: glpksolution = "simplex"; break;
  case 1: glpksolution = "interior"; break;
  default: glpksolution = "simplex"; break;
  }
  ParameterBlob.SetStringParameter("GLPKSolutionMethod", glpksolution);
  // snopt6 parameters
  ParameterBlob.SetIntParameter("Snopt6MajorIterations", 
				*Ampl::TheSnopt6MajorIterations);
  ParameterBlob.SetIntParameter("Snopt6MinorIterations", 
				*Ampl::TheSnopt6MinorIterations);
  // Rconvexifier parameters
  if (*Ampl::TheConvexifierOutAmpl == 1) {
    ParameterBlob.SetBoolParameter("ConvexifierOutAmpl", true);
  } else {
    ParameterBlob.SetBoolParameter("ConvexifierOutAmpl", false);
  }
  // Rsymmgroup parameters
  ParameterBlob.SetIntParameter("SymmgroupOutType", *Ampl::TheSymmgroupOutType);
  // Rfbbtfp parameters
  // [none]
  // Rmilp2gph parameters
  // [none]
  // vns parameters
  ParameterBlob.SetDoubleParameter("VNSEpsilon", *Ampl::TheVNSEpsilon);
  ParameterBlob.SetIntParameter("VNSKmax", *Ampl::TheVNSKmax);
  ParameterBlob.SetIntParameter("VNSKmin", *Ampl::TheVNSKmin);
  ParameterBlob.SetIntParameter("VNSNeighbourhood", 
				*Ampl::TheVNSNeighbourhood);
  ParameterBlob.SetBoolParameter("VNSWarmUp", *Ampl::TheVNSWarmUp);
  std::string vnslocalsolvername;
  switch(*Ampl::TheVNSLocalSolver) {
  case 0: vnslocalsolvername = "snopt"; break;
  case 1: vnslocalsolvername = "tabu"; break;
  case 2: vnslocalsolvername = "glpk"; break;
  case 4: vnslocalsolvername = "gomory"; break;
  case 5: vnslocalsolvername = "limitedbranch"; break;
  case 6: vnslocalsolvername = "localbranch"; break;
  default: vnslocalsolvername = "NONE"; break;
  }
  if (vnslocalsolvername != "NONE") {
    ParameterBlob.SetStringParameter("VNSLocalSolver", vnslocalsolvername);
  }
  ParameterBlob.SetIntParameter("VNSMaxTime", *Ampl::TheVNSMaxTime);
  ParameterBlob.SetIntParameter("VNSSamples", *Ampl::TheVNSSamples);  

  // limited branch parameters
  ParameterBlob.SetDoubleParameter("LimitedBranchArbitrary", 
				   *Ampl::TheLimitedBranchArbitrary);
  std::string limitedbranchlocalsolvername;
  switch(*Ampl::TheLimitedBranchLocalSolver) {
  case 0: limitedbranchlocalsolvername = "snopt"; break;
  case 1: limitedbranchlocalsolvername = "tabu"; break;
  case 2: limitedbranchlocalsolvername = "glpk"; break;
  case 4: limitedbranchlocalsolvername = "gomory"; break;
  case 5: limitedbranchlocalsolvername = "limitedbranch"; break;
  case 6: limitedbranchlocalsolvername = "localbranch"; break;
  default: limitedbranchlocalsolvername = "NONE"; break;
  }

  // local branch parameters
  ParameterBlob.SetIntParameter("LocalBranchKmax", 
				*Ampl::TheLocalBranchKmax); 
  ParameterBlob.SetBoolParameter("LocalBranchWarmUp", 
				 *Ampl::TheLocalBranchWarmUp); 

  // disaggr parameters
  ParameterBlob.SetIntParameter("DisaggrK", *Ampl::TheDisaggrK);

  // quartic convex reformulator parameters
  ParameterBlob.SetIntParameter("RQuarticConvexType", 
				*Ampl::TheRQuarticConvexType); 

  // initialise variable values
  x = new double [NumberOfVariables];

  // simplify the problem if possible
  Simplifier(true);

  // set formulation information 
  SetProblemLinear();

  return asl;
}

Expression Problem::ParseAMPLExpressionTree(Ampl::expr* e) {
  using namespace Ampl;
  using namespace std;

  Expression ret(0.0);
  Expression kval(0.0);

  char vnamebuf[16];
  strcpy(vnamebuf, "x");
  char* vthename = vnamebuf + 1;

  efunc *op;
  expr **ep;
  int opnum;
  float number;
  float  y, z, k, ind;
  char buffer[50];
  char * pEnd;
  double t;
  
  op = e->op;
  opnum = Intcast op;

  switch(opnum) {
    
  case OPPLUS:
    ret = ParseAMPLExpressionTree(e->L.e);
    ret = SumLink(ret, ParseAMPLExpressionTree(e->R.e));
    break;
    
  case OPMINUS:
    ret = ParseAMPLExpressionTree(e->L.e);
    ret = DifferenceLink(ret, ParseAMPLExpressionTree(e->R.e));
    break;
    
  case OPMULT:
    ret = ParseAMPLExpressionTree(e->L.e);
    ret = ProductLink(ret, ParseAMPLExpressionTree(e->R.e));
    break;
    
  case OPDIV:
    ret = ParseAMPLExpressionTree(e->L.e);
    ret = FractionLink(ret, ParseAMPLExpressionTree(e->R.e));
    break;
    
  case OPREM:
    printf("Problem::ParseAMPLExpressionTree: remainder -- not implemented\n");
    exit(1);
    break;
    
  case OPPOW:
    ret = ParseAMPLExpressionTree(e->L.e);
    ret = PowerLink(ret, ParseAMPLExpressionTree(e->R.e));
    break;
    
  case OPLESS:
    printf("Problem::ParseAMPLExpressionTree: less -- not implemented\n");
    exit(1);
    break;
    
  case MINLIST:
    printf("Problem::ParseAMPLExpressionTree: min -- not implemented\n");
    exit(1);
    break;
    
  case MAXLIST:
    printf("Problem::ParseAMPLExpressionTree: max -- not implemented\n");
    exit(1);
    break;

  case FLOOR:
    printf("Problem::ParseAMPLExpressionTree: floor -- not implemented\n");
    exit(1);
    break;

  case CEIL:
    printf("Problem::ParseAMPLExpressionTree: ceil -- not implemented\n");
    exit(1);
    break;

  case ABS:
    printf("Problem::ParseAMPLExpressionTree: abs -- not implemented\n");
    exit(1);
    break;
    
  case OPUMINUS:
    ret = MinusLink(ParseAMPLExpressionTree(e->L.e));
    break;
    
  case OPIFnl:
    printf("Problem::ParseAMPLExpressionTree: if_nl -- not implemented\n");
    exit(1);
    break;
        
  case OP_tanh:
    ret = TanhLink(ParseAMPLExpressionTree(e->L.e));
    break;
    
  case OP_tan:
    ret = TanLink(ParseAMPLExpressionTree(e->L.e));
    break;
    
  case OP_sqrt:
    ret = SqrtLink(ParseAMPLExpressionTree(e->L.e));
    break;
    
  case OP_sinh:
    ret = SinhLink(ParseAMPLExpressionTree(e->L.e));
    break;
    
  case OP_sin:
    ret = SinLink(ParseAMPLExpressionTree(e->L.e));
    break;
    
  case OP_log10:
    printf("Problem::ParseAMPLExpressionTree: log10 -- not implemented\n");
    exit(1);
    break;
    
  case OP_log:
    ret = LogLink(ParseAMPLExpressionTree(e->L.e));
    break;
    
  case OP_exp:
    ret = ExpLink(ParseAMPLExpressionTree(e->L.e));
    break;
    
  case OP_cosh:
    ret = CoshLink(ParseAMPLExpressionTree(e->L.e));
    break;
    
  case OP_cos:
    ret = CosLink(ParseAMPLExpressionTree(e->L.e));
    break;
    
  case OP_atanh:
    printf("Problem::ParseAMPLExpressionTree: atanh -- not implemented\n");
    exit(1);
    break;
    
  case OP_atan2:
    printf("Problem::ParseAMPLExpressionTree: atan2 -- not implemented\n");
    exit(1);
    break;
    
  case OP_atan:
    printf("Problem::ParseAMPLExpressionTree: atan -- not implemented\n");
    exit(1);
    break;
    
  case OP_asinh:
    printf("Problem::ParseAMPLExpressionTree: asinh -- not implemented\n");
    exit(1);
    break;
    
  case OP_asin:
    printf("Problem::ParseAMPLExpressionTree: asin -- not implemented\n");
    exit(1);
    break;
    
  case OP_acosh:
    printf("Problem::ParseAMPLExpressionTree: acosh -- not implemented\n");
    exit(1);
    break;
    
  case OP_acos:
    printf("Problem::ParseAMPLExpressionTree: acos -- not implemented\n");
    exit(1);
    break;
    
  case OPSUMLIST:
    for (ep = e->L.ep; ep < e->R.ep; *ep++) {
      ret = SumLink(ret, ParseAMPLExpressionTree(*ep));
    }
    break;
    
  case OPintDIV:
    printf("Problem::ParseAMPLExpressionTree: int_div -- not implemented\n");
    exit(1);
    break;
    
  case OPprecision:
    printf("Problem::ParseAMPLExpressionTree: precision -- not implemented\n");
    exit(1);
    break;

  case OPround:
    printf("Problem::ParseAMPLExpressionTree: round -- not implemented\n");
    exit(1);
    break;
    
  case OPtrunc:
    printf("Problem::ParseAMPLExpressionTree: trunc -- not implemented\n");
    exit(1);
    break;
    
  case OP1POW:
    kval->SetValue(e->R.en->v);
    ret = PowerLink(ParseAMPLExpressionTree(e->L.e), kval);
    break;
    
  case OP2POW:
    kval->SetValue(2.0);
    ret = PowerLink(ParseAMPLExpressionTree(e->L.e), kval);
    break;

  case OPCPOW:
    kval->SetValue(e->L.en->v);
    ret = PowerLink(kval, ParseAMPLExpressionTree(e->R.e));
    break;
    
  case OPFUNCALL:
    printf("Problem::ParseAMPLExpressionTree: func_call -- not implemented\n");
    exit(1);
    break;
    
  case OPNUM:
    t = ((expr_n*)e)->v;
    ret->SetValue(t);
    break;
    
  case OPPLTERM:
    printf("Problem::ParseAMPLExpressionTree: pl_term -- not implemented\n");
    exit(1);
    break;
    
  case OPIFSYM:
    printf("Problem::ParseAMPLExpressionTree: if_sym -- not implemented\n");
    exit(1);
    break;
    
  case OPHOL:
    printf("Problem::ParseAMPLExpressionTree: string -- not implemented\n");
    exit(1);
    break;

  case OPVARVAL:
    ret->SetOpType(VAR);
    ret->SetExponent(1.0);
    ret->SetCoeff(1.0);
    sprintf(vthename, "%d", e->a + 1);
    ret->SetVarName(vnamebuf);
    ret->SetVarIndex(e->a + 1);
    break;
    
  default:
    printf("Problem::ParseAMPLExpressionTree: unknownop -- not implemented\n");
    exit(1);
  }
  
  return ret;
}
 
bool Problem::IsAMPLExpressionZero(Ampl::expr* e) {
  using namespace Ampl;
  efunc* op = e->op;
  int opnum = Intcast op;
  if (opnum == OPNUM && fabs(e->dL) < EPSILONTOLERANCE) {
    return true;
  } else {
    return false;
  }
}

void Problem::Simplifier(bool theamplflag) {
  int i, j, v1, v2;
  double c, c1, c2, c3, e;
  bool goonflag = true;
  bool ischanged = false;
  std::vector<bool> IsConstrDeleted(NumberOfConstraints, false);
  std::vector<bool> IsVarDeleted(NumberOfVariables, false);

  if (!NoSimplifier) {
    while(goonflag) {
      goonflag = false;
      // simplify expressions
      for(i = 0; i < NumberOfObjectives; i++) {
	ischanged = Simplify(&Obj[i]->Function);
      }
      if (!goonflag && ischanged) {
	goonflag = true;
      }
      for (i = 0; i < NumberOfConstraints; i++) {
	if (!IsConstrDeleted[i]) {
	  ischanged = Simplify(&Constr[i]->Function);
	}
	if (!goonflag && ischanged) {
	  goonflag = true;
	}
      }
      // check in constraints for -inf < constr < inf situations and 
      // just delete them
      for (i = 0; i < NumberOfConstraints; i++) {
	if (!IsConstrDeleted[i]) {
	  c1 = Constr[i]->LB;
	  c2 = Constr[i]->UB;
	  if (c1 <= -ROSEINFINITY && c2 >= ROSEINFINITY) {
	    IsConstrDeleted[i] = true;
	  }
	}
      }
      // check in constraints for l < constant < u situations and 
      // just delete them
      for (i = 0; i < NumberOfConstraints; i++) {
	if (!IsConstrDeleted[i]) {
	  c1 = Constr[i]->LB;
	  c2 = Constr[i]->UB;
	  if (Constr[i]->Function->IsConstant()) {
	    c3 = Constr[i]->Function->GetValue();
	    IsConstrDeleted[i] = true;
	    // check feasibility!
	    if (c1 > c3 || c2 < c3) {
	      std::cerr << "rose::Simplify: constraint " << c1 << " <= "
		   << c3 << " <= " << c2 << " infeasible\n";
	      assert(c1 <= c3 && c3 <= c2);
	    }
	  }
	}
      }
      // check in constraints for l < x < u situations and move them in ranges
      for (i = 0; i < NumberOfConstraints; i++) {
	if (!IsConstrDeleted[i]) {
	  if (Constr[i]->Function->IsLeaf()) {
	    c1 = Constr[i]->LB;
	    c2 = Constr[i]->UB;
	    if (Constr[i]->Function->IsVariable()) {
	      v1 = Constr[i]->Function->GetVarIndex();
	      c = Constr[i]->Function->GetCoeff();
	      e = Constr[i]->Function->GetExponent();
	      if (c < 0) {
		// coefficient < 0, invert lower and upper bound
                double tmp = c2;
		c2 = c1 / c;
		c1 = tmp / c;
	      } else if (c > 0) {
		// coefficient > 0, just divide bounds by it
		c1 = c1 / c;
		c2 = c2 / c;
	      } else {
		// c = 0, just delete the constr
		IsConstrDeleted[i] = true;
		ischanged = true;
	      }
	      if (c != 0 && e != 0) {
		// WARNING: check that maths below also works with e negative,
		// e not integer, e odd.
		c1 = pow(c1, 1 / e);
		c2 = pow(c2, 1 / e);
		// update var range if needed
		if (c1 > GetVariable(v1)->LB) {
		  GetVariable(v1)->LB = c1;
		} 
		if (c2 < GetVariable(v1)->UB) {
		  GetVariable(v1)->UB = c2;
		}
                if (GetVariable(v1)->UB < GetVariable(v1)->LB) {
		  std::cerr << "rose::Simplify: constraint " << GetVariable(v1)->LB << " <= "
		    << c << " <= " << GetVariable(v1)->UB << " infeasible\n";
                }
		// delete the constr
		IsConstrDeleted[i] = true;
		ischanged = true;
	      }
	    } else {
	      c = Constr[i]->Function->GetValue();
	      if (! (c1 <= c && c <= c2)) {
		std::cerr << "rose::Simplify: constraint " << c1 << " <= "
		     << c << " <= " << c2 << " infeasible\n";
		assert(c1 <= c && c <= c2);
	      }
	    }
	  }
	}
      }
      if (!goonflag && ischanged)
	goonflag = true;
      // check in variable ranges for variables that can be replaced by a const
      if (NumberOfVariables > 1) {
	for(i = 0; i < NumberOfVariables; i++) {
	  if (Var[i]->Persistent = false && 
	      fabs(Var[i]->UB - Var[i]->LB) < SYMSIMPTOLERANCE) {
	    // range is c < var < c, replace variable with constant
	    c = (Var[i]->UB + Var[i]->LB) / 2;
	    Var[i]->UB = c;
	    Var[i]->LB = c;
	    ischanged = true;
	    for(j = 0; j < NumberOfObjectives; j++) {
	      Obj[j]->Function->VariableToConstant(GetVariableID(i + 1), c);
	    }
	    for(j = 0; j < NumberOfConstraints; j++) {
	      Constr[j]->Function->VariableToConstant(GetVariableID(i + 1), c);
	    }
	    IsVarDeleted[i] = true;
	  }
	}
      }
      if (!goonflag && ischanged) {
	goonflag = true;
      }
      if (!theamplflag) {
	// check in constraints for 0 < x - y < 0 situations
	for(i = 0; i < NumberOfConstraints; i++) {
	  if (!IsConstrDeleted[i]) {
	    if (Constr[i]->Function->GetSize() == 2 &&
		Constr[i]->Function->GetNode(0)->IsVariable() &&
		Constr[i]->Function->GetNode(1)->IsVariable() &&
		fabs(Constr[i]->LB) < SYMSIMPTOLERANCE &&
		fabs(Constr[i]->UB) < SYMSIMPTOLERANCE) {
	      // two subnodes, both of which variables, and constr bounds at 0
	      if ((Constr[i]->Function->GetOpType() == DIFFERENCE ||
		   Constr[i]->Function->GetOpType() == SUM) &&
		  Constr[i]->Function->GetNode(0)->GetExponent() == 1 &&
		  Constr[i]->Function->GetNode(1)->GetExponent() == 1 &&
		  Constr[i]->Function->GetNode(0)->GetCoeff() != 0 &&
		  Constr[i]->Function->GetNode(1)->GetCoeff() != 0) {
		// difference or sum with exponents set to 1 and
		// nonzero coefficients
		v1 = Constr[i]->Function->GetNode(0)->GetVarIndex();
		v2 = Constr[i]->Function->GetNode(1)->GetVarIndex();
		c1 = Constr[i]->Function->GetNode(0)->GetCoeff();
		c2 = Constr[i]->Function->GetNode(1)->GetCoeff();
		if (Constr[i]->Function->GetOpType() == SUM) {
		  c2 = -c2;
		}
		// situation is c1*v1 - c2*v2 = 0, i.e. v2=(c1/c2)v1
		if (v1 == v2) {
		  if (c1 != -c2) {
		    // only solution has x_v1 = 0
		    // TODO100524: search and replace x_v1 by 0 throughout
		  } else {
		    // constraint is 0 < x - x < 0, just delete it
		    IsConstrDeleted[i] = true;
		    ischanged = true;
		  }
		} else {
		  // constraint is a potential candidate, check that
		  // at least one var is not immutable
		  std::string vn;
		  double theb1, theb2;
		  bool theint;
		  if (Var[VarIDToLocalIndex[v2] - 1]->Persistent == false) {
		    // substitute v2 with v1 throughout
		    IsConstrDeleted[i] = true;
		    IsVarDeleted[VarIDToLocalIndex[v2] - 1] = true;
		    ischanged = true;		
		    theb1 = Var[VarIDToLocalIndex[v1] - 1]->UB;
		    theb2 = Var[VarIDToLocalIndex[v2] - 1]->UB;
		    if (theb2 < theb1) {
		      Var[VarIDToLocalIndex[v1] - 1]->UB = theb2;
		    }
		    theb1 = Var[VarIDToLocalIndex[v1] - 1]->LB;
		    theb2 = Var[VarIDToLocalIndex[v2] - 1]->LB;
		    if (theb2 > theb1) {
		      Var[VarIDToLocalIndex[v1] - 1]->LB = theb2;
		    }
		    theint = Var[VarIDToLocalIndex[v1] - 1]->IsIntegral;
		    if (!theint) {
		      // if v1 continuous, copy integrality property from v2
		      Var[VarIDToLocalIndex[v1] - 1]->IsIntegral = 
			Var[VarIDToLocalIndex[v2] - 1]->IsIntegral;
		    }
		    vn = Var[VarIDToLocalIndex[v1] - 1]->Name;
		    for(j = 0; j < NumberOfObjectives; j++) {
		      Obj[j]->Function->ReplaceVariable(v2, v1, vn, c1/c2);
		    }
		    for(j = 0; j < NumberOfConstraints; j++) {
		      if (!IsConstrDeleted[j]) {
			Constr[j]->Function->ReplaceVariable(v2, v1, vn, 
							     c1/c2);
		      }
		    }
		  } else if (Var[VarIDToLocalIndex[v1]-1]->Persistent == 
			     false) {
		    // substitute v1 with v2 throughout
		    IsConstrDeleted[i] = true;
		    IsVarDeleted[VarIDToLocalIndex[v1] - 1] = true;
		    ischanged = true;
		    theb1 = Var[VarIDToLocalIndex[v2] - 1]->UB;
		    theb2 = Var[VarIDToLocalIndex[v1] - 2]->UB;
		    if (theb2 < theb1) {
		      Var[VarIDToLocalIndex[v2] - 1]->UB = theb2;
		    }
		    theb1 = Var[VarIDToLocalIndex[v2] - 1]->LB;
		    theb2 = Var[VarIDToLocalIndex[v1] - 1]->LB;
		    if (theb2 > theb1) {
		      Var[VarIDToLocalIndex[v2] - 1]->LB = theb2;
		    }
		    theint = Var[VarIDToLocalIndex[v2] - 1]->IsIntegral;
		    if (!theint) {
		      // if v2 continuous, copy integrality property from v1
		      Var[VarIDToLocalIndex[v2] - 1]->IsIntegral = 
			Var[VarIDToLocalIndex[v1] - 1]->IsIntegral;
		    }
		    vn = Var[VarIDToLocalIndex[v2] - 1]->Name;
		    for(j = 0; j < NumberOfObjectives; j++) {
		      Obj[j]->Function->ReplaceVariable(v1, v2, vn, c2/c1);
		    }
		    for(j = 0; j < NumberOfConstraints; j++) {
		      if (!IsConstrDeleted[j]) {
			Constr[j]->Function->ReplaceVariable(v1, v2, vn, 
							     c2/c1);
		      }
		    }
		  }
		} 
	      }
	    }
	  }
	}
      }
      if (!goonflag && ischanged) {
	goonflag = true;
      }
    }
  }

  // extract additive constant from 1st objective (single-objective
  // NLP solvers usually require additive constants separated from func)
  double d = Obj[0]->Function->RemoveAdditiveConstant();
  Obj1AdditiveConstant += d;
  // extract the additive constants from constrs and put them in ranges  
  for(i = 0; i < NumberOfConstraints; i++) {
    if (!IsConstrDeleted[i]) {
      d = Constr[i]->Function->RemoveAdditiveConstant();
      if (d != 0) {
	Constr[i]->LB -= d;
	Constr[i]->UB -= d;
      }
    }
  }

  // now perform the actual deletion.
  DeletionFlag = false;

  // delete variables:
  std::vector<int> VarDelIndex;
  for(i = 0; i < NumberOfVariables; i++) {
    if (IsVarDeleted[i]) {
      VarDelIndex.push_back(Var[i]->ID);
    }
  }
  int delvars = VarDelIndex.size();
  std::vector<Variable*>::iterator vi;
  for(i = 0; i < delvars; i++) {
    for(vi = Var.begin(); vi != Var.end(); vi++) {
      if ((*vi)->ID == VarDelIndex[i]) {
	Var.erase(vi);
	VarIDToLocalIndex[VarDelIndex[i]] = 0;
	NumberOfVariables--;
	DeletionFlag = true;
	break;
      }
    }
  }

  // delete constraints:
  std::vector<int> ConstrDelIndex;
  for(i = 0; i < NumberOfConstraints; i++) {
    if (IsConstrDeleted[i]) {
      ConstrDelIndex.push_back(Constr[i]->ID);
    }
  }
  int delconstrs = ConstrDelIndex.size();
  std::vector<Constraint*>::iterator ci;
  for(i = 0; i < delconstrs; i++) {
    for(ci = Constr.begin(); ci != Constr.end(); ci++) {
      if ((*ci)->ID == ConstrDelIndex[i]) {
	Constr.erase(ci);
	ConstrIDToLocalIndex[ConstrDelIndex[i]] = 0;
	NumberOfConstraints--;
	DeletionFlag = true;
	break;
      }
    }
  }
  
  
  // remap VarIDToLocalIndex, ConstrIDToLocalIndex 
  // (...and ObjIDToLocalIndex if used)
  i = 1;
  for(vi = Var.begin(); vi != Var.end(); vi++) {
    VarIDToLocalIndex[(*vi)->ID] = i;
    i++;
  }

#ifdef OBJECTIVEDELETION
  i = 1;
  std::vector<Objective*>::iterator oi;
  for(oi = Obj.begin(); oi != Obj.end(); ci++) {
    ConstrIDToLocalIndex[(*oi)->ID] = i;
    i++;
  }
#endif

  i = 1;
  for(ci = Constr.begin(); ci != Constr.end(); ci++) {
    ConstrIDToLocalIndex[(*ci)->ID] = i;
    i++;
  }

  // re-build nonlinear parts
  for(i = 0; i < NumberOfObjectives; i++) {  
    Objective* theObj = Obj[i];
    theObj->NonlinearPart.SetTo(theObj->Function->GetNonlinearPart());
  }
  for(i = 0; i < NumberOfConstraints; i++) {
    Constraint* theConstr = Constr[i];
    theConstr->NonlinearPart.SetTo(theConstr->Function->GetNonlinearPart());
  }

  // re-build FET functions
  for(i = 0; i < NumberOfObjectives; i++) {  
    Objective* theObj = Obj[i];
    DeleteFastEvalTree(theObj->FunctionFET);
    DeleteFastEvalTree(theObj->NonlinearPartFET);
    theObj->FunctionFET = theObj->Function->GetFastEvalTree();
    theObj->NonlinearPartFET = theObj->NonlinearPart->GetFastEvalTree();
  }
  for(i = 0; i < NumberOfConstraints; i++) {
    Constraint* theConstr = Constr[i];
    DeleteFastEvalTree(theConstr->FunctionFET);
    DeleteFastEvalTree(theConstr->NonlinearPartFET);
    theConstr->FunctionFET = theConstr->Function->GetFastEvalTree();
    theConstr->NonlinearPartFET = theConstr->NonlinearPart->GetFastEvalTree();
  }  
}

// output to std::cout
std::ostream& operator<<(std::ostream& s, Problem& gp) {
  int i;
  int v = gp.GetNumberOfVariables();
  int c = gp.GetNumberOfConstraints();
  s << "# ROSE problem: " << gp.GetName() << std::endl;
  s << "# Problem has " << v << " variables and " << c << " constraints\n";
  s << "# Variables: " << std::endl;
  if (gp.HasDeleted()) {    
    s << "# Reformulation process deleted some variables or constraints\n";
  }
  s << std::endl;
  s << "variables = ";
  for(i = 1; i <= v - 1; i++) {
    s << gp.GetVariableLI(i)->LB << " < ";
    s << gp.GetVariableLI(i)->Name << " < ";
    s << gp.GetVariableLI(i)->UB << " / ";
    if (gp.GetVariableLI(i)->IsIntegral) {
      s << "Integer" << "," << std::endl;
    } else {
      s << "Continuous" << "," << std::endl;
    }
  }
  s << gp.GetVariableLI(v)->LB << " < ";
  s << gp.GetVariableLI(v)->Name << " < ";
  s << gp.GetVariableLI(v)->UB << " / ";
  if (gp.GetVariableLI(v)->IsIntegral) {
    s << "Integer" << ";" << std::endl;
  } else {
    s << "Continuous" << ";" << std::endl;
  }
  s << std::endl;
  s << "# Objective Function:" << std::endl;
  s << "objfun = ";
  if (gp.GetOptimizationDirection(1) == Minimization) {
    s << "min";
  } else {
    s << "max";
  }
  s << " [ " << gp.GetObjectiveLI(1)->Function->ToString() 
    << " ];" << std::endl;
  s << std::endl;
  s << "# Constraints: " << std::endl;
  s << "constraints = ";
  for(i = 1; i <= c - 1; i++) {
    s << "[ " << gp.GetConstraintLI(i)->LB << " < ";
    s << gp.GetConstraintLI(i)->Function->ToString() << " < ";
    s << gp.GetConstraintLI(i)->UB << " ]," << std::endl;
  }
  if (c != 0) {
    s << "[ " << gp.GetConstraintLI(c)->LB << " < ";
    s << gp.GetConstraintLI(c)->Function->ToString() << " < ";
    s << gp.GetConstraintLI(c)->UB << " ];" << std::endl;
    s << std::endl;
  } else {
    s << "0;" << std::endl << std::endl;
  }
  s << "# Starting Point: " << std::endl;
  s << "startingpoint = ";
  for(i = 1; i <= v - 1; i++) {
    s << gp.GetStartingPointLI(i) << ", ";
  }
  s << gp.GetStartingPointLI(v) << ";" << std::endl;
  s << std::endl;
  s << "# end of problem " << gp.GetName() << std::endl;
  return s;
}
