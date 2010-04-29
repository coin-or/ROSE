/*
** Name:       solver_limitedbranch.cxx
** Author:     Leo Liberti
** Source:     GNU C++
** Purpose:    Solver limitedbranch implementation 
** License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
** History:    0511-- work started
*/

#include <cstdio>
#include <cstring>
#include <cassert>
#include <sys/time.h>
#include "newsolver.h"
#include "solver.h"
#include "solver_limitedbranch.h"
#include "utils.h"

#define FIRSTOBJ 1
#define EPSILON 1e-4
#define LIMITEDBRANCHMAXNUMBEROFCUTS 500
#define LIMITEDBRANCHARBITRARYCHOICE 0.2
//#define LIMITEDBRANCHFIX

// SNOPT6SOLVER CLASS METHODS

LimitedBranchSolver::LimitedBranchSolver() {
  TheName = "limitedbranch";
  NumberOfVariables = 0;
  NumberOfConstraints = 0;
  IsSolved = false;
  OptDir = Minimization;
  OptDirCoeff = 1;
  ManualOptDir = false;
  OptimalObjVal = ROSEINFINITY;
  Epsilon = EPSILON;
  LocalSolverName = "snopt";
  LocalSolver = NULL;
  MaxRunningTime = 0;
  Relaxed = false;
  Quiet = false;
  ArbitraryChoice = LIMITEDBRANCHARBITRARYCHOICE;
  MaxNumberOfCuts = LIMITEDBRANCHMAXNUMBEROFCUTS;
  NumberOfCuts = 0;
  NeedsCuts = false;
  MINLPSols = 0;
  MainSolver = "";
  AmplFlag = false;
}

LimitedBranchSolver::~LimitedBranchSolver() {
  if (LocalSolver) {
    DeleteSolver(LocalSolver);
  }
}

void LimitedBranchSolver::SetProblem(Problem* p) {
  InProb = p; 
  InitProb = false;
}

bool LimitedBranchSolver::CanSolve(int problemtype) {
  return true;
}

void LimitedBranchSolver::SetMaxNumberOfCuts(int mncuts) {
  // this is only possible if the problem is uninitialized
  assert(!InitProb);
  MaxNumberOfCuts = mncuts;
}

void LimitedBranchSolver::Initialize(bool force = false) {

  if (!InitProb || force) {
    InitProb = true;

    // check the input problem is set to something
    if (InProb == NULL) {
      cerr << "LimitedBranchSolver::Initialize(): error: InProb is NULL\n";
      assert(InProb != NULL);
    }

    NumberOfVariables = InProb->GetNumberOfVariables();
    NumberOfConstraints = InProb->GetNumberOfConstraints();

    // parameters
    ReplaceParams(InProb->GetParamsRef());
    int np = ParameterBlob.GetNumberOfParameters();
    for(int i = 1; i <= np; i++) {
      if (ParameterBlob.GetParameterName(i) == "LimitedBranchLocalSolver") {
	LocalSolverName = ParameterBlob.GetParameterStringValue(i);
      } else if (ParameterBlob.GetParameterName(i) == "LimitedBranchMaxTime") {
	MaxRunningTime = (double) ParameterBlob.GetParameterIntValue(i);
	assert(MaxRunningTime >= 0);
      } else if (ParameterBlob.GetParameterName(i) == "LimitedBranchEpsilon") {
	Epsilon = ParameterBlob.GetParameterDoubleValue(i);
	assert(Epsilon > 0);
      } else if (ParameterBlob.GetParameterName(i)=="LimitedBranchArbitrary") {
	ArbitraryChoice = ParameterBlob.GetParameterDoubleValue(i);
	assert(ArbitraryChoice >= 0 && ArbitraryChoice <= 1);
      } else if (ParameterBlob.GetParameterName(i)=="LimitedBranchNeedsCuts") {
	int needscuts = ParameterBlob.GetParameterIntValue(i);
	assert(needscuts == 0 || needscuts == 1);
	NeedsCuts = (bool) needscuts;
      } else if (ParameterBlob.GetParameterName(i) == "RelaxInteger") {
	Relaxed = true;
      } else if (ParameterBlob.GetParameterName(i) == "LimitedBranchQuiet" ||
		 ParameterBlob.GetParameterName(i) == "Quiet") {
	Quiet = ParameterBlob.GetParameterBoolValue(i);
      } else if (ParameterBlob.GetParameterName(i) == "AmplSolver") {
	AmplFlag = ParameterBlob.GetParameterBoolValue(i);
      } else if (ParameterBlob.GetParameterName(i) == "MainSolver") {
	MainSolver = ParameterBlob.GetParameterStringValue(i);
      }
    }

    if (MainSolver != TheName && AmplFlag) {
      Quiet = true;
    }

    // variable integrality
    for(int i = 1; i <= NumberOfVariables; i++) {
      integrality.push_back(InProb->GetVariableLI(i)->IsIntegral);
      if (Relaxed) {
	integrality[i - 1] = false;
      }
    }

    // create and configure local solver
    LocalSolver = NewSolver(LocalSolverName);
    LocalSolver->SetProblem(InProb);
    LocalSolver->ReplaceParams(ParameterBlob);
    if (NeedsCuts) {
      LocalSolver->SetMaxNumberOfCuts(MaxNumberOfCuts);
    }
    LocalSolver->Initialize(false);

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
}

int LimitedBranchSolver::Solve(bool reinitialize) {

  int ret = 0;
  if (reinitialize || !InitProb) {
    Initialize(reinitialize);
  }
  bool isfeas1, isfeas2;
  double disc = 0;
  double themaxdisc = 0;
  double tmp;
  int themaxdiscindex = -1;

  // initialize randomizer
  struct timeval theTV;
  struct timezone theTZ;
  gettimeofday(&theTV, &theTZ);
  srandom(theTV.tv_usec);

  // re-initialize some stuff for warm starts
  fstar = ROSEINFINITY;
  for(int i = 0; i < NumberOfVariables; i++) {
    vlb[i] = InProb->GetVariableLI(i + 1)->LB;
    vub[i] = InProb->GetVariableLI(i + 1)->UB;
  }

  // initial point
  for(int i = 0; i < NumberOfVariables; i++) {
    LocalSolver->SetStartingPointLI(i + 1, x[i]);
  }

  // timings
  starttime = getcputime(tstart);

  // do it
  ret = LocalSolver->Solve();
  LocalSolver->GetSolutionLI(f, x);
  for(int i = 0; i < NumberOfVariables; i++) {
    InProb->SetCurrentVariableValueLI(i + 1, x[i]);
  }
  isfeas1 = InProb->TestConstraintsFeasibility(Epsilon, disc);
  if (!isfeas1) {
    if (!Quiet) {
      cout << "LimitedBranchSolver: WARNING: "
	   << "infeasible initial solution (LSret = " << ret << ")" << endl;
    }
  } else {
    // do NOT set f* here
    for(int i = 0; i < NumberOfVariables; i++) {
      xstar[i] = x[i];
    }
    if (!Quiet) {
      cout << "LimitedBranchSolver: starting with objective f* = " << fstar 
	   << " (LSret = " << ret << ")" << endl;
    }
  }  
  isfeas2 = InProb->TestVariablesFeasibility(Epsilon, disc);
  if (!isfeas2) {
    for(int i = 0; i < NumberOfVariables; i++) {
      if (integrality[i]) {
	tmp = fabs(x[i] - rint(x[i]));
	if (tmp > Epsilon && tmp > themaxdisc) {
	  themaxdisc = tmp;
	  themaxdiscindex = i + 1;
	}
      }
    }
    // choose biggest discrepancy for branching
    if (themaxdiscindex != -1) {
      Branch(themaxdiscindex, vlb, vub);
    }
  } else {
    // all is integral, OK
    ;
  }

  // final stopwatch
  stoptime = getcputime(tend);
  usertimelap = tend[0] - tstart[0];
  if (!Quiet) {
    cout << "LimitedBranchSolver: finished search after " << usertimelap
	 << " seconds of user CPU, f* = " << fstar << endl;
  }

  // after solution: set optimal values in problem
  IsSolved = true;
  for(int i = 0; i < NumberOfVariables; i++) {
    InProb->SetCurrentVariableValueLI(i + 1, xstar[i]);
    InProb->SetOptimalVariableValueLI(i + 1, xstar[i]);
  }
  OptimalObjVal = InProb->EvalObj(FIRSTOBJ);

  double constrdiscr = 0;
  double vardiscr = 0;
  IsFeasible = InProb->TestConstraintsFeasibility(Epsilon, constrdiscr);
  bool constrfeas = IsFeasible;
  IsFeasible = InProb->TestVariablesFeasibility(Epsilon, vardiscr);
  bool varfeas = IsFeasible;
  for(int i = 0; i < NumberOfVariables; i++) {
    InProb->SetOptimalVariableValueLI(i + 1, xstar[i]);
  }
  OptimalObjVal = InProb->EvalObj(FIRSTOBJ);
  if (!Quiet) {
    cout << "LimitedBranchSolver: best solution f* = " 
	 << OptimalObjVal << endl;
  }
  InProb->SetOptimalObjectiveValue(FIRSTOBJ, OptimalObjVal);
  InProb->SetSolved(true);    
  if (constrfeas && varfeas) {
    InProb->SetFeasible(1);
    ret = 0;
  } else {
    if (!Quiet) {
      if (!constrfeas) {
	cout << "LimitedBranchSolver: infeasible by at least " 
	     << constrdiscr << " in constraints\n";
      }
      if (!varfeas) {
	cout << "LimitedBranchSolver: infeasible by at least " 
	     << vardiscr << " in integrality/ranges\n";
      }
    }
    InProb->SetFeasible(-1);
    ret = 1;
  }

  return ret;
}

int LimitedBranchSolver::Branch(int localvarindex,
				std::vector<double> vlb, 
				std::vector<double> vub) {
  int ret = 0;
  bool isfeas;
  double disc = 0;
  double tmp;
  double themaxdisc = 0;
  int themaxdiscindex = -1;
  double saveval;
  double leftobj;
  double rightobj;
  double r;
  int optdir;
  double optdircoeff;

  // verify cpu time constrained termination
  if (MaxRunningTime > 0) {
    stoptime = getcputime(tend);
    usertimelap = tend[0] - tstart[0];
    if (usertimelap > (double) MaxRunningTime) {
      if (!Quiet) {
	cout << "LimitedBranchSolver: running time exceeds limit, exiting\n";
      }
      return ret;
    }
  }

  // branch on localvarindex: first, try to fix variables to integer values
  double leftub = floor(x[localvarindex - 1]);
#ifdef LIMITEDBRANCHFIX
  double leftlb = leftub;
#else
  double leftlb = vlb[localvarindex - 1];
#endif
  if (fabs(leftub - leftlb) < Epsilon) {
    leftub = leftlb;
  }
  double rightlb = ceil(x[localvarindex - 1]);
#ifdef LIMITEDBRANCHFIX
  double rightub = rightlb;
#else
  double rightub = vub[localvarindex - 1];
#endif
  if (fabs(rightub - rightlb) < Epsilon) {
    rightlb = rightub;
  }
  if (!Quiet) {
    cout << "LimitedBranchSolver: branching on x" << localvarindex 
	 << "=" << x[localvarindex - 1] << "; left: [" 
	 << leftlb << "," << leftub << "], right: ["
	 << rightlb << "," << rightub << "]" << endl;
  }

  // choose L or R branch according to the choice which decreases 
  // (min) the objective  
  saveval = InProb->GetCurrentVariableValueLI(localvarindex);
  // try with leftub
  InProb->SetCurrentVariableValueLI(localvarindex, leftub);
  leftobj = InProb->EvalObj(FIRSTOBJ);
  // now with rightlb
  InProb->SetCurrentVariableValueLI(localvarindex, rightlb);
  rightobj = InProb->EvalObj(FIRSTOBJ);
  // reset val
  InProb->SetCurrentVariableValueLI(localvarindex, saveval);
  r = boundedrandom(0,1);
  if (fabs(leftobj - rightobj) < Epsilon || r < ArbitraryChoice) {
    // if objfun does not depend significantly on x_localvarindex,
    // and in any case sometimes, 
    // decide randomly whether to look at left or right first
    r = boundedrandom(0,1);
    if (r <= 0.5) {
      // do right first
      double tmp = leftlb;
      leftlb = rightlb;
      rightlb = tmp;
      tmp = leftub;
      leftub = rightub;
      rightub = tmp;
    }
  } else {
    optdir = InProb->GetOptimizationDirection(FIRSTOBJ);
    if (optdir == Minimization) {
      optdircoeff = 1;
    } else {
      optdircoeff = -1;
    }
    leftobj *= optdircoeff;
    rightobj *= optdircoeff;
    if (leftobj > rightobj) {
      // do right first
      double tmp = leftlb;
      leftlb = rightlb;
      rightlb = tmp;
      tmp = leftub;
      leftub = rightub;
      rightub = tmp;      
    }
  }

  for(int LR = 0; LR < 2; LR++) {
    themaxdisc = 0;
    themaxdiscindex = -1;
    if (LR == 0) {
      // first branch (can be left or right according to previous choice)
      vlb[localvarindex - 1] = leftlb;
      vub[localvarindex - 1] = leftub;
    } else if (LR == 1) {
      // second branch
      vlb[localvarindex - 1] = rightlb;
      vub[localvarindex - 1] = rightub;
    }
    // solve locally
    for(int i = 1; i <= NumberOfVariables; i++) {
      LocalSolver->SetVariableLBLI(i, vlb[i - 1]);
      LocalSolver->SetVariableUBLI(i, vub[i - 1]);
    }
    ret = LocalSolver->Solve();
    LocalSolver->GetSolutionLI(f, x);
    for(int i = 0; i < NumberOfVariables; i++) {
      InProb->SetCurrentVariableValueLI(i + 1, x[i]);
    }

    // is local solution feasible in constraints?
    isfeas = InProb->TestConstraintsFeasibility(Epsilon, disc);
    if (isfeas) {
      // recurse only if feasible, otherwise prune search at this branch

      // is local solution integer feasible?
      isfeas = InProb->TestVariablesFeasibility(Epsilon, disc);
      if (isfeas) {

	// feasible integer, we have a candidate solution
	if (f < fstar) {
	  if (!Quiet) {
	    cout << "LimitedBranchSolver: improved f* = " << f 
		 << " integer feasible (LSret = " << ret << ")" << endl;
	  }
	  fstar = f;
	  for(int i = 0; i < NumberOfVariables; i++) {
	    xstar[i] = x[i];
	  }
	  ret++;
	  MINLPSols++;
	}
	
      } else {
	// not integer feasible, but may be near - try rounding

	for(int i = 0; i < NumberOfVariables; i++) {
	  if (integrality[i]) {
	    InProb->SetCurrentVariableValueLI(i + 1, rint(x[i]));
	  }
	}

	// is rounding feasible in the constraints?
	isfeas = InProb->TestConstraintsFeasibility(Epsilon, disc);	
	if (!isfeas) {
	  // no it isn't, re-set the old variable values
	  for(int i = 0; i < NumberOfVariables; i++) {
	    if (integrality[i]) {
	      InProb->SetOptimalVariableValueLI(i + 1, x[i]);
	    }
	  }
	  if (!Quiet) {
	    cout << "LimitedBranchSolver: local solver returned integer "
		 << "infeasible " << "(LSret = " << ret << ")" << endl;
	  }
	  
	} else {	  
	  // rounding worked, feasible integer candidate solution

	  f = InProb->EvalObj(FIRSTOBJ);
	  if (f < fstar) {
	    if (!Quiet) {
	      cout << "LimitedBranchSolver: improved f* = " << f 
		   << " integer feasible (rounding)" << endl;
	    }
	    fstar = f;
	    for(int i = 0; i < NumberOfVariables; i++) {
	      x[i] = rint(x[i]);
	      xstar[i] = x[i];
	    }
	    ret++;
	    MINLPSols++;
	  }
	}
      }
      
      // look for branching variable (largest discrepancy)
      for(int i = 0; i < NumberOfVariables; i++) {
	if (integrality[i] && i != localvarindex - 1) {
	  tmp = fabs(x[i] - rint(x[i]));
	  if (tmp > Epsilon && tmp > themaxdisc) {
	    themaxdisc = tmp;
	    themaxdiscindex = i + 1;
	  }
	}
      }
      
      // if themaxdiscindex == -1, all is integral, OK, otherwise, 
      // choose biggest discrepancy for branching -
      // branch only if no solutions found yet
      if (themaxdiscindex != -1 && MINLPSols == 0) {
	Branch(themaxdiscindex, vlb, vub);
      }
    }
  }
}

int LimitedBranchSolver::AddCut(std::vector<std::pair<int, double> >& cutinfo, 
				double LB, double UB) {
  assert(LocalSolver);
  NumberOfCuts++;
  return LocalSolver->AddCut(cutinfo, LB, UB);
}
int LimitedBranchSolver::AddCut(Expression e, double LB, double UB) {
  assert(LocalSolver);
  NumberOfCuts++;
  return LocalSolver->AddCut(e, LB, UB);
}
void LimitedBranchSolver::EnableCut(int cutID) {
  assert(LocalSolver);
  LocalSolver->EnableCut(cutID);
}
void LimitedBranchSolver::DisableCut(int cutID) {
  assert(LocalSolver);
  LocalSolver->DisableCut(cutID);
}
void LimitedBranchSolver::SetCutLB(int cutID, double LB) {
  assert(LocalSolver);
  LocalSolver->SetCutLB(cutID, LB);
}
double LimitedBranchSolver::GetCutLB(int cutID) {
  assert(LocalSolver);
  return LocalSolver->GetCutLB(cutID);
}
void LimitedBranchSolver::SetCutUB(int cutID, double UB) {
  assert(LocalSolver);
  LocalSolver->SetCutUB(cutID, UB);
}
double LimitedBranchSolver::GetCutUB(int cutID) {
  assert(LocalSolver);
  return LocalSolver->GetCutUB(cutID);
}

void LimitedBranchSolver::SetOptimizationDirection(int theoptdir) {
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

void LimitedBranchSolver::GetSolution(map<int,double>& objfunval, 
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

void LimitedBranchSolver::GetSolutionLI(vector<double>& objfunval, 
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

void LimitedBranchSolver::GetSolutionLI(double& objfunval, 
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


void LimitedBranchSolver::SetStartingPoint(int varID, double sp) {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

void LimitedBranchSolver::SetStartingPointLI(int localindex, double sp) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

double LimitedBranchSolver::GetStartingPoint(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

double LimitedBranchSolver::GetStartingPointLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

void LimitedBranchSolver::SetVariableLB(int varID, double LB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableLBLI(localindex, LB);
}

double LimitedBranchSolver::GetVariableLB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void LimitedBranchSolver::SetVariableUB(int varID, double UB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableUBLI(localindex, UB);
}

double LimitedBranchSolver::GetVariableUB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}

void LimitedBranchSolver::SetVariableLBLI(int localindex, double LB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vlb[localindex - 1] = LB;
}

double LimitedBranchSolver::GetVariableLBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void LimitedBranchSolver::SetVariableUBLI(int localindex, double UB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vub[localindex - 1] = UB;
}

double LimitedBranchSolver::GetVariableUBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}
