/*
** Name:       solver_VNS.cxx
** Author:     Leo Liberti
** Source:     GNU C++
** Purpose:    Solver VNS implementation
** License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
** History:    050921 work started
*/

#include <cstdio>
#include <cstring>
#include <cassert>
#include <cassert>
#include <sys/time.h>
#include "newsolver.h"
#include "solver_vns.h"
#include "utils.h"

#define FIRSTOBJ 1
#define EPSILON 1e-4
#define VNSWARMUPTOLERANCE 0.5
#define VNSWARMUPMAXITN 6
#define VNSLSMAXTIME 3
#define VNSMAXNUMBEROFCUTS 500
#define VNSOUTBOUNDSMIDPOINTPROB 0.3

// SNOPT6SOLVER CLASS METHODS

VNSSolver::VNSSolver() {
  TheName = "vns";
  NumberOfVariables = 0;
  NumberOfConstraints = 0;
  IsSolved = false;
  OptDir = Minimization;
  OptDirCoeff = 1;
  ManualOptDir = false;
  OptimalObjVal = ROSEINFINITY;
  LocalSolverName = "snopt";
  LocalSolver = NULL;
  MINLPLocalSolverName = "limitedbranch";
  MaxRunningTime = 0;
  Kmin = 0;
  Kmax = 10;
  Samples = 1;
  Quiet = false;
  f = ROSEINFINITY;
  fstar = ROSEINFINITY;
  IsFeasible = false;
  Epsilon = EPSILON;
  NeighbourhoodType = 0;
  IsMINLP = false;
  Relaxed = false;
  MaxNumberOfCuts = VNSMAXNUMBEROFCUTS;
  NumberOfCuts = 0;
  NeedsCuts = false;
  MainSolver = "";
  AmplFlag = false;
  WarmUp = true;
}

VNSSolver::~VNSSolver() {
  if (LocalSolver) {
    DeleteSolver(LocalSolver);
  }
}

void VNSSolver::SetMaxNumberOfCuts(int mncuts) {
  // this is only possible if the problem is uninitialized
  assert(!InitProb);
  MaxNumberOfCuts = mncuts;
}

void VNSSolver::SetProblem(Problem* p) {
  InProb = p; 
  InitProb = false;
}

bool VNSSolver::CanSolve(int problemtype) {
  // can solve everything
  return true;
}

void VNSSolver::Initialize(bool force = false) {

  if (!InitProb || force) {
    InitProb = true;

    // check the input problem is set to something
    if (InProb == NULL) {
      cerr << "VNSSolver::Initialize(): error: InProb is NULL\n";
      assert(InProb != NULL);
    }

    NumberOfVariables = InProb->GetNumberOfVariables();
    NumberOfConstraints = InProb->GetNumberOfConstraints();

    // parameters
    ReplaceParams(InProb->GetParamsRef());
    int np = ParameterBlob.GetNumberOfParameters();
    bool forcemethod = false;
    bool forcequiet = false;
    for(int i = 1; i <= np; i++) {
      if (ParameterBlob.GetParameterName(i) == "VNSLocalSolver") {
	LocalSolverName = ParameterBlob.GetParameterStringValue(i);
	forcemethod = true;
      } else if (ParameterBlob.GetParameterName(i) == "VNSLocalMINLPSolver") {
	MINLPLocalSolverName = ParameterBlob.GetParameterStringValue(i);
      } else if (ParameterBlob.GetParameterName(i) == "VNSMaxTime") {
	MaxRunningTime = (double) ParameterBlob.GetParameterIntValue(i);
	assert(MaxRunningTime >= 0);
      } else if (ParameterBlob.GetParameterName(i) == "VNSKmin") {
	Kmin = ParameterBlob.GetParameterIntValue(i);
	assert(Kmin >= 0);
      } else if (ParameterBlob.GetParameterName(i) == "VNSKmax") {
	Kmax = ParameterBlob.GetParameterIntValue(i);
	assert(Kmin >= 0);
      } else if (ParameterBlob.GetParameterName(i) == "VNSNeighbourhood") {
	NeighbourhoodType = ParameterBlob.GetParameterIntValue(i);
      } else if (ParameterBlob.GetParameterName(i) == "VNSSamples" ||
		 ParameterBlob.GetParameterName(i) == "VNSTrials") {
	Samples = ParameterBlob.GetParameterIntValue(i);
	assert(Samples > 0);
      } else if (ParameterBlob.GetParameterName(i) == "VNSEpsilon") {
	Epsilon = ParameterBlob.GetParameterDoubleValue(i);
	assert(Epsilon > 0);
      } else if (ParameterBlob.GetParameterName(i) == "RelaxInteger") {
	if (ParameterBlob.GetParameterIntValue(i) == 1) {
	  Relaxed = true;
	}
      } else if (ParameterBlob.GetParameterName(i) == "VNSQuiet" ||
		 ParameterBlob.GetParameterName(i) == "Quiet") {
	Quiet = ParameterBlob.GetParameterBoolValue(i);
	forcequiet = true;
      } else if (ParameterBlob.GetParameterName(i) == "VNSWarmUp") {
	WarmUp = ParameterBlob.GetParameterBoolValue(i);
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
      if (!IsMINLP && integrality[i - 1]) {
	IsMINLP = true;
      }
    }

    // create and configure local solver
    LocalSolver = NewSolver(LocalSolverName);
    if (IsMINLP && !LocalSolver->CanSolve(ROSE_MINLP)) {
      if (forcemethod) {
	cout << "VNSSolver: WARNING: local solver can't deal with integer "
	     << "variables, relaxing\n";
      } else {
	DeleteSolver(LocalSolver);
	LocalSolverName = "limitedbranch";
	LocalSolver = NewSolver(LocalSolverName);
      }
    }
    if (!Quiet) {
      cout << "VNSSolver: using " << LocalSolverName << " as local solver\n";
    }
    LocalSolver->SetProblem(InProb);
    savequiet = ParameterBlob.GetBoolParameter("Quiet");
    if (!forcequiet) {
      ParameterBlob.SetBoolParameter("Quiet", true);
    } else {
      ParameterBlob.SetBoolParameter("Quiet", Quiet);    
    }
    double lsmaxtime = VNSLSMAXTIME;
    ParameterBlob.SetIntParameter("LimitedBranchMaxTime", (int) lsmaxtime);
    InProb->ReplaceParams(ParameterBlob);
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

int VNSSolver::Solve(bool reinitialize) {

  int ret = 0;
  if (reinitialize || !InitProb) {
    Initialize(reinitialize);
  }

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

  // do it
  bool termination = false;

  // timings
  double tstart[2], tend[2];
  double starttime, stoptime, usertimelap;
  starttime = getcputime(tstart);
  int samplecount = 0;
  int majoriteration = 0;
  int k = Kmin;

  // warm up local solver first
  if (WarmUp) {
    double feastolerance = VNSWARMUPTOLERANCE;
    int maxfeasitn = VNSWARMUPMAXITN;
    int feasitn = 0;
    bool isfeas;
    while(!termination) {
      GetRandomPointInNeighbourhood(x, xstar, vlb, vub, Kmax, 1);
      for(int i = 0; i < NumberOfVariables; i++) {
	LocalSolver->SetStartingPointLI(i + 1, x[i]);
      }
      ret = LocalSolver->Solve();
      LocalSolver->GetSolutionLI(f, x);
      for(int i = 0; i < NumberOfVariables; i++) {
	InProb->SetCurrentVariableValueLI(i + 1, x[i]);
      }
      double disc = 0;
      isfeas = InProb->TestConstraintsFeasibility(feastolerance, disc);
      if (isfeas) {
	termination = true;
      }
      if (!Quiet) {
	cout << "VNSSolver: preprocessing: f* = " << f << ", LSret = " 
	     << ret << ", feasible = ";
	if (isfeas) {
	  cout << "true";
	} else {
	  cout << "false";
	}
	cout << " (tol=" << feastolerance << ")" << endl;
      }
      feasitn++;
      feastolerance = (feasitn + 1)*VNSWARMUPTOLERANCE;
      if (feasitn > maxfeasitn) {
	termination = true;
      }
    }
    if (!isfeas) {
      // try one last time with x=0
      for(int i = 0; i < NumberOfVariables; i++) {
	LocalSolver->SetStartingPointLI(i + 1, 0.0);
      }
      ret = LocalSolver->Solve();
      LocalSolver->GetSolutionLI(f, x);    
      if (!Quiet) {
	cout << "VNSSolver: preprocessing (last): f* = " << f << ", LSret = " 
	     << ret << ", feasible = ";
	if (isfeas) {
	  cout << "true";
	} else {
	  cout << "false";
	}
	cout << " (tol=" << feastolerance << ")" << endl;
      }
    }
    for(int i = 0; i < NumberOfVariables; i++) {
      xstar[i] = x[i];
      // fstar = f;
      fstar = ROSEINFINITY; // don't keep this as a valid solution
    }
  }

  // start VNS proper
  if (!Quiet) {
    cout << "VNSSolver: starting VNS proper\n";
  }
  termination = false;
  while(!termination) {
    // find local solution
    for(samplecount = 0; samplecount < Samples; samplecount++) {
      GetRandomPointInNeighbourhood(x, xstar, vlb, vub, k, NeighbourhoodType);
      for(int i = 0; i < NumberOfVariables; i++) {
	LocalSolver->SetStartingPointLI(i + 1, x[i]);
      }
      ret = LocalSolver->Solve();
      int saveret = ret;
      if (LocalSolverName == "snopt" && (ret == 3 || ret == 4 || ret == 9)) {
	// using SNOPT and returned fishy code, check feasibility manually
	LocalSolver->GetSolutionLI(f, x);
	for(int i = 0; i < NumberOfVariables; i++) {
	  InProb->SetCurrentVariableValueLI(i + 1, x[i]);
	}
	double disc = 0;
	bool isfeas = InProb->TestConstraintsFeasibility(Epsilon, disc);
	if (isfeas) {
	  ret = 0;
	}
      }
      if (ret == 0) {
	// retrieve current variable values
	LocalSolver->GetSolutionLI(f, x);
	// see if current solution better than current best
	if (f < fstar) {
	  if (!Quiet) {
	    cout << "VNSSolver: improved f* = " << f
		 << " at itn. " << majoriteration + 1 << ", neighb. " << k + 1
		 << ", sample " << samplecount + 1 
		 << ", LSret=" << saveret << endl;
	  }
	  fstar = f;
	  for(int j = 0; j < NumberOfVariables; j++) {
	    xstar[j] = x[j];
	  }
	  k = Kmin - 1;
	  break;
	}
      }
    }
    // termination conditions
    if (MaxRunningTime > 0) {
      stoptime = getcputime(tend);
      usertimelap = tend[0] - tstart[0];
      if (usertimelap > (double) MaxRunningTime) {
	termination = true;
	if (!Quiet) {
	  cout << "VNSSolver: running time exceeds limit, exiting\n";
	}
      }
    }
    if (k == Kmax) {
      termination = true;
      if (!Quiet) {
	cout << "VNSSolver: k exceeds Kmax, exiting\n";
      }
    }
    // increase neighbourhood
    k++;
    majoriteration++;
  }

  // final stopwatch
  stoptime = getcputime(tend);
  usertimelap = tend[0] - tstart[0];
  if (fstar >= ROSEINFINITY) {
    fstar = f;
  }
  if (!Quiet) {
    cout << "VNSSolver: finished search after " << usertimelap
	 << " seconds of user CPU, f* = " << fstar << endl;
  }

  // after solution: set optimal values in problem
  IsSolved = true;
  for(int i = 0; i < NumberOfVariables; i++) {
    InProb->SetCurrentVariableValueLI(i + 1, xstar[i]);
  }
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
    cout << "VNSSolver: best solution f* = " << OptimalObjVal << endl;
  }
  InProb->SetOptimalObjectiveValue(FIRSTOBJ, OptimalObjVal);
  InProb->SetSolved(true);    
  if (constrfeas && varfeas) {
    InProb->SetFeasible(1);
    ret = 0;
  } else {
    if (!constrfeas && varfeas) {
      ret = 2;
    } else if (constrfeas && !varfeas) {
      ret = 3;
    } else if (!constrfeas && !varfeas) {
      ret = 1;
    }
      
    if (!Quiet) {
      if (!constrfeas) {
	cout << "VNSSolver: infeasible by at least " 
	     << constrdiscr << " in constraints\n";
      }
      if (!varfeas) {
	cout << "VNSSolver: infeasible by at least " 
	     << vardiscr << " in integrality/ranges\n";
      }
    }
    InProb->SetFeasible(-1);
  }

  // re-set the original problem flags
  ParameterBlob.SetBoolParameter("Quiet", savequiet);
  InProb->ReplaceParams(ParameterBlob);

  return ret;
}

void VNSSolver::SetOptimizationDirection(int theoptdir) {
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

void VNSSolver::GetSolution(map<int,double>& objfunval, 
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

void VNSSolver::GetSolutionLI(vector<double>& objfunval, 
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

void VNSSolver::GetSolutionLI(double& objfunval, 
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


void VNSSolver::SetStartingPoint(int varID, double sp) {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

void VNSSolver::SetStartingPointLI(int localindex, double sp) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

double VNSSolver::GetStartingPoint(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

double VNSSolver::GetStartingPointLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

void VNSSolver::SetVariableLB(int varID, double LB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableLBLI(localindex, LB);
}

double VNSSolver::GetVariableLB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void VNSSolver::SetVariableUB(int varID, double UB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableUBLI(localindex, UB);
}

double VNSSolver::GetVariableUB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}

void VNSSolver::SetVariableLBLI(int localindex, double LB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vlb[localindex - 1] = LB;
}

double VNSSolver::GetVariableLBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void VNSSolver::SetVariableUBLI(int localindex, double UB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vub[localindex - 1] = UB;
}

double VNSSolver::GetVariableUBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}

int VNSSolver::AddCut(std::vector<std::pair<int, double> >& cutinfo, 
				double LB, double UB) {
  assert(LocalSolver);
  NumberOfCuts++;
  return LocalSolver->AddCut(cutinfo, LB, UB);
}
int VNSSolver::AddCut(Expression e, double LB, double UB) {
  assert(LocalSolver);
  NumberOfCuts++;
  return LocalSolver->AddCut(e, LB, UB);
}
void VNSSolver::EnableCut(int cutID) {
  assert(LocalSolver);
  LocalSolver->EnableCut(cutID);
}
void VNSSolver::DisableCut(int cutID) {
  assert(LocalSolver);
  LocalSolver->DisableCut(cutID);
}
void VNSSolver::SetCutLB(int cutID, double LB) {
  assert(LocalSolver);
  LocalSolver->SetCutLB(cutID, LB);
}
double VNSSolver::GetCutLB(int cutID) {
  assert(LocalSolver);
  return LocalSolver->GetCutLB(cutID);
}
void VNSSolver::SetCutUB(int cutID, double UB) {
  assert(LocalSolver);
  LocalSolver->SetCutUB(cutID, UB);
}
double VNSSolver::GetCutUB(int cutID) {
  assert(LocalSolver);
  return LocalSolver->GetCutUB(cutID);
}

void VNSSolver::GetRandomPointInNeighbourhood(std::vector<double>& x,
					      std::vector<double>& c,
					      std::vector<double>& l,
					      std::vector<double>& u,
					      int k, 
					      int neighbourhoodtype) {
  int n = x.size();
  assert(n > 0);
  assert(n == l.size());
  assert(n == u.size());
  double rc;
  double rv;
  double top;
  double thevlb, thevub;
  if (neighbourhoodtype == 0) {
    // hyperrectangular shell-like neighbourhoods
    // get random direction in (-1,1)
    vector<double> d(n, 0.0);
    randomvector(d, 1);
    // scale by inf norm
    double dinfnorm = Linfnorm(d);
    assert(dinfnorm != 0);
    for(int i = 0; i < n; i++) {
      d[i] /= dinfnorm;
    }
    // get random radius in [r_k,r_k+1]  
    double tmp = 1.0 / (double) Kmax;
    double r = boundedrandom(tmp * k, tmp * (k+1));
    for(int i = 0; i < n; i++) {
      x[i] = r * d[i];
    }
    // now transform x affinely to be centered at c and bounded by l,u
    for(int i = 0; i < n; i++) {
      if (integrality[i] == false) {
	if (x[i] > 0) {
	  x[i] *= u[i] - c[i];
	} else {
	  x[i] *= c[i] - l[i];
	}
	x[i] += c[i];
      }
    }
  } else if (neighbourhoodtype == 1) {
    // hyperrectangular (filled) neighbourhoods
    for(int i = 0; i < n; i++) {
      if (integrality[i] == false) {
	thevlb = c[i] - ((double) k) * (c[i] - l[i]) / ((double) Kmax);
	thevub = c[i] + ((double) k) * (u[i] - c[i]) / ((double) Kmax);
	rc = (double) random() / (double) RAND_MAX;
	top = thevub - thevlb;
	rv = rc * top + thevlb;
	x[i] = rv;
      }
    } 
  }

  // deal with integer variables
  double theproblim = (double) k / (double) Kmax;
  if (IsMINLP) {
    for(int i = 0; i < n; i++) {
      if (integrality[i]) {
	thevlb = c[i] - ((double) k) * (c[i] - l[i]) / ((double) Kmax);
	thevub = c[i] + ((double) k) * (u[i] - c[i]) / ((double) Kmax);
	top = thevub - thevlb;
	rv = rc * top + thevlb;
	rv = rint(rv);
	x[i] = rint(rv);
	if (fabs(x[i] - rv) < Epsilon) {
	  // x value hasn't changed, change it artificially with probability
	  // k / k_max	  
	  double theprob = boundedrandom(0, 1);
	  if (theprob < theproblim) {
	    if (fabs(x[i] - l[i]) < Epsilon) {
	      x[i] = x[i] + k;
	    } else if (fabs(u[i] - x[i]) < Epsilon) {
	      x[i] = x[i] - k;
	    }
	    if (x[i] < l[i]) {
	      x[i] = l[i];
	    } else if (x[i] > u[i]) {
	      x[i] = u[i];
	    }
	  }
	}
      }
    }
  }
  
  for(int i = 0; i < NumberOfVariables; i++) {
    if (x[i] < l[i]) {
      x[i] = l[i];
      rc = boundedrandom(0,1);
      if (rc < VNSOUTBOUNDSMIDPOINTPROB) {
	x[i] = (l[i] + u[i]) / 2;
      }
    } else if (x[i] > u[i]) {
      x[i] = u[i];
      if (rc < VNSOUTBOUNDSMIDPOINTPROB) {
	x[i] = (l[i] + u[i]) / 2;
      }
    }
  }   
  
}
