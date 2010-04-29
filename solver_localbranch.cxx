/*
** Name:       solver_localbranch.cxx
** Author:     Leo Liberti
** Source:     GNU C++
** Purpose:    Solver localbranch implementation
** License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
** History:    051122 (flying west of Morocco's coast) work started
*/

#include <cstdio>
#include <cstring>
#include <cassert>
#include <sys/time.h>
#include "newsolver.h"
#include "solver.h"
#include "solver_localbranch.h"
#include "utils.h"

#define FIRSTOBJ 1
#define EPSILON 1e-4
#define LOCALBRANCHRECURSIONLIMIT 20
// for solution method = "localbranching"
#define LOCALBRANCHDEFAULTK 10
// for solution method = "vns"
#define LOCALBRANCHKMAX 10
#define LOCALBRANCHKMIN 1
#define LOCALBRANCHSAMPLES 1
#define LOCALBRANCHWARMUPTOLERANCE 0.5
#define LOCALBRANCHWARMUPMAXITN 5
#define LOCALBRANCHWARMUPFACTOR 0.2
#define LOCALBRANCHLSMAXTIME 3

LocalBranchSolver::LocalBranchSolver() {
  TheName = "localbranch";
  NumberOfVariables = 0;
  NumberOfConstraints = 0;
  IsSolved = false;
  OptDir = Minimization;
  OptDirCoeff = 1;
  ManualOptDir = false;
  OptimalObjVal = ROSEINFINITY;
  Epsilon = EPSILON;
  LocalSolverName = "vns";
  LocalSolver = NULL;
  SolutionMethod = "vns";
  MaxRunningTime = 0;
  Kmax = LOCALBRANCHKMAX;
  Kmin = LOCALBRANCHKMIN;
  Samples = LOCALBRANCHSAMPLES;
  RecursionLimit = LOCALBRANCHRECURSIONLIMIT;
  MaxRank = 0;
  MainSolver = "";
  AmplFlag = false;
  WarmUp = true;
}

LocalBranchSolver::~LocalBranchSolver() {
  if (LocalSolver) {
    DeleteSolver(LocalSolver);
  }
}

void LocalBranchSolver::SetProblem(Problem* p) {
  InProb = p; 
  InitProb = false;
}

bool LocalBranchSolver::CanSolve(int problemtype) {
  // can solve everything
  return true;
}

void LocalBranchSolver::Initialize(bool force = false) {

  if (!InitProb || force) {
    InitProb = true;

    // check the input problem is set to something
    if (InProb == NULL) {
      cerr << "LocalBranchSolver::Initialize(): error: InProb is NULL\n";
      assert(InProb != NULL);
    }

    NumberOfVariables = InProb->GetNumberOfVariables();
    NumberOfConstraints = InProb->GetNumberOfConstraints();

    // parameters
    bool forcemethod = false;
    bool forcesolver = false;
    bool forcequiet = false;
    ReplaceParams(InProb->GetParamsRef());
    int np = ParameterBlob.GetNumberOfParameters();
    for(int i = 1; i <= np; i++) {
      if (ParameterBlob.GetParameterName(i) == "LocalBranchLocalSolver") {
	LocalSolverName = ParameterBlob.GetParameterStringValue(i);
	forcesolver = true;
      } else if (ParameterBlob.GetParameterName(i) == 
		 "LocalBranchSolutionMethod") {
	SolutionMethod = ParameterBlob.GetParameterStringValue(i);
	forcemethod = true;
      } else if (ParameterBlob.GetParameterName(i) == "LocalBranchMaxTime") {
	MaxRunningTime = (double) ParameterBlob.GetParameterIntValue(i);
	assert(MaxRunningTime >= 0);
      } else if (ParameterBlob.GetParameterName(i) == "LocalBranchEpsilon") {
	Epsilon = ParameterBlob.GetParameterDoubleValue(i);
	assert(Epsilon > 0);
      } else if (ParameterBlob.GetParameterName(i) == "LocalBranchKmax") {
	Kmax = ParameterBlob.GetParameterIntValue(i);
	assert(Kmax > 0);
      } else if (ParameterBlob.GetParameterName(i) == "LocalBranchKmin") {
	Kmin = ParameterBlob.GetParameterIntValue(i);
	assert(Kmin > 0);
      } else if (ParameterBlob.GetParameterName(i) == 
		 "LocalBranchRecursionLimit") {
	RecursionLimit = ParameterBlob.GetParameterIntValue(i);
	assert(RecursionLimit > 0);
      } else if (ParameterBlob.GetParameterName(i) == "LocalBranchDefaultK") {
	Kmin = ParameterBlob.GetParameterIntValue(i);
	assert(Kmin > 0);
      } else if (ParameterBlob.GetParameterName(i) == "LocalBranchSamples") {
	Samples = ParameterBlob.GetParameterIntValue(i);
	assert(Samples > 0);
      } else if (ParameterBlob.GetParameterName(i) == "LocalBranchWarmUp") {
	WarmUp = ParameterBlob.GetParameterBoolValue(i);
      } else if (ParameterBlob.GetParameterName(i) == "LocalBranchQuiet" ||
		 ParameterBlob.GetParameterName(i) == "Quiet") {
	Quiet = ParameterBlob.GetParameterBoolValue(i);
	forcequiet = true;
      } else if (ParameterBlob.GetParameterName(i) == "AmplSolver") {
	AmplFlag = ParameterBlob.GetParameterBoolValue(i);
      } else if (ParameterBlob.GetParameterName(i) == "MainSolver") {
	MainSolver = ParameterBlob.GetParameterStringValue(i);
      }
    }

    if (MainSolver != TheName && AmplFlag) {
      Quiet = true;
    }

    // set kmin at default k for solution method = "localbranching"
    if (SolutionMethod != "vns" && Kmin == LOCALBRANCHKMIN) {
      Kmin = LOCALBRANCHDEFAULTK;
    }

    // variable integrality
    for(int i = 1; i <= NumberOfVariables; i++) {
      integrality.push_back(InProb->GetVariableLI(i)->IsIntegral);
      if (!IsMINLP && integrality[i - 1]) {
	IsMINLP = true;
      }
    }

    // create and configure local solver
    LocalSolver = NewSolver(LocalSolverName);
    LocalSolver->SetProblem(InProb);
    LocalSolver->ReplaceParams(ParameterBlob);
    if (IsMINLP && !LocalSolver->CanSolve(ROSE_MINLP)) {
      if (forcesolver) {
	cerr << "LocalBranchSolver: local solver " << LocalSolverName 
	     << " can't deal with integer variables\n";
	exit(120);
      } else {
	DeleteSolver(LocalSolver);
	LocalSolverName = "limitedbranch";
	LocalSolver = NewSolver(LocalSolverName);
      }
    } 

    if (!Quiet) {
      cout << "LocalBranchSolver: using ";
      if (IsMINLP) {
	cout << "(mixed-integer) ";
      } else {
	cout << "(continuous) ";
      }
      cout << LocalSolverName << " as local solver\n";
    }  
    LocalSolver->SetProblem(InProb);
    savequiet = ParameterBlob.GetBoolParameter("Quiet");
    if (!forcequiet) {
      ParameterBlob.SetBoolParameter("Quiet", true);
    } else {
      ParameterBlob.SetBoolParameter("Quiet", Quiet);
    }
    double lsmaxtime = LOCALBRANCHLSMAXTIME;
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

int LocalBranchSolver::LocalBranching(int k,
				      double currval, 
				      std::vector<double> currsol,
				      double& newval,
				      std::vector<double>& newsol, 
				      int rank) {
  // perform local branching with constant k
  using namespace std;
  int ret = 0;
  if (rank > MaxRank) {
    MaxRank = rank;
  }
  // create local branching constraint extended to general integer variables
  // (see Fischetti and Lodi, Local Branching, last page before References)
  double lhsterm = 0;
  vector<pair<int,double> > cutinfo;
  pair<int,double> p;
  double mu;
  double LB, UB;
  int constrtype = 0; // 0=Local branching constr, 1=Eq.(1) in fisch/lodi
  for(int i = 0; i < NumberOfVariables; i++) {
    if ((IsMINLP && integrality[i]) || !IsMINLP) {
      mu = 1 / (vub[i] - vlb[i]);
      p.first = i + 1;
      if (fabs(vub[i] - currsol[i]) < Epsilon) {
	lhsterm = lhsterm + mu * vub[i];
	p.second = -mu;
	cutinfo.push_back(p);
      } else if (fabs(currsol[i] - vlb[i]) < Epsilon) {
	lhsterm = lhsterm - mu * vlb[i];
	p.second = mu;
	cutinfo.push_back(p);
      } else {
        // only add this term if variable is nonnegative
	if (vlb[i] >= 0) {
	  p.second = mu;
	  cutinfo.push_back(p);
	}
      }
    }
  }
  LB = -ROSEINFINITY;
  UB = (double) k - lhsterm;
  if (IsMINLP) {
    UB = floor(UB);
  }
  if (cutinfo.size() == 0) {
    // cut is empty, it may happen if current sol is far from var
    // bounds and variables may be negative; in this case, 
    // use a generalization of cut (1) in Fischetti and Lodi's paper
    constrtype = 1;
    lhsterm = 0;
    for(int i = 0; i < NumberOfVariables; i++) {
      if ((IsMINLP && integrality[i]) || !IsMINLP) {
	lhsterm += fabs(currsol[i]);
	p.first = i + 1;
	p.second = fabs(currsol[i]);
	cutinfo.push_back(p);
      }
    }
    lhsterm = (1.0 - ((double) k / (double) Kmax)) * lhsterm;
    LB = lhsterm;
    UB = ROSEINFINITY;
    if (IsMINLP) {
      LB = ceil(LB);
    }    
  }
  if (cutinfo.size() == 0) {
    cerr << "LocalBranchSolver: can't identify a valid LocBra cut at itn. "
	 << k << endl;
    exit(235);
  }
  // add cut to solver
  int thecutid = LocalSolver->AddCut(cutinfo, LB, UB);
  // try to find new incumbent
  LocalSolver->Solve();
  LocalSolver->GetSolutionLI(newval, newsol);
  for(int i = 0; i < NumberOfVariables; i++) {
    InProb->SetCurrentVariableValueLI(i + 1, newsol[i]);
  }
  double disc = 0;
  bool isfeas1 = InProb->TestConstraintsFeasibility(Epsilon, disc);
  bool isfeas2 = InProb->TestVariablesFeasibility(Epsilon, disc);
  if (isfeas1 && isfeas2 && newval < currval) {
    // new solution feasible and improved, update
    currval = newval;
    for(int i = 0; i < NumberOfVariables; i++) {
      currsol[i] = newsol[i];
    }
    // reverse k-constraint sense
    if (constrtype == 0) {
      if (IsMINLP) {
	LB = (double) (k + 1 - lhsterm);
	LB = ceil(LB);
      } else {
	LB = (double) (k - lhsterm + Epsilon);
      }
      UB = ROSEINFINITY;
    } else if (constrtype == 1) {
      if (IsMINLP) {
	UB = floor(lhsterm);
      } else {
	UB = lhsterm + Epsilon;
      }
      LB = ROSEINFINITY;
    }
    LocalSolver->SetCutLB(thecutid, LB);
    LocalSolver->SetCutUB(thecutid, UB);
    if (rank < RecursionLimit) {
      // recurse using the new solution
      LocalBranching(k, currval, currsol, newval, newsol, rank + 1);
    }
  } else {
    // new solution is feasible but not improved, or infeasible:
    // original original incumbent
    newval = currval;
    for(int i = 0; i < NumberOfVariables; i++) {
      newsol[i] = currsol[i];
    }
  }  
  return ret;
}

int LocalBranchSolver::Solve(bool reinitialize) {

  int ret = 0;
  if (reinitialize || !InitProb) {
    Initialize(reinitialize);
  }

  // initialize randomizer
  struct timeval theTV;
  struct timezone theTZ;
  gettimeofday(&theTV, &theTZ);
  srandom(theTV.tv_usec);

  // timings
  double tstart[2], tend[2];
  double starttime, stoptime, usertimelap;
  starttime = getcputime(tstart);
  int samplecount = 0;
  int majoriteration = 0;
  int k = Kmin;
  int saveret;

  // re-initialize some stuff for warm starts
  fstar = ROSEINFINITY;
  for(int i = 0; i < NumberOfVariables; i++) {
    vlb[i] = InProb->GetVariableLI(i + 1)->LB;
    vub[i] = InProb->GetVariableLI(i + 1)->UB;
  }

  double disc = 0;
  bool isfeas1, isfeas2;

  // warm up local solver first
  if (WarmUp) {
    double feastolerance = LOCALBRANCHWARMUPTOLERANCE;
    int maxfeasitn = LOCALBRANCHWARMUPMAXITN;
    int feasitn = 0;
    bool isfeas;
    bool wutermination = false;
    while(!wutermination) {
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
	wutermination = true;
      }
      if (!Quiet) {
	cout << "LocalBranchSolver: preprocessing: f* = " << f << ", LSret = " 
	     << ret << ", feasible = ";
	if (isfeas) {
	  cout << "true";
	} else {
	  cout << "false";
	}
	cout << " (tol=" << feastolerance << ")" << endl;
      }
      feastolerance *= LOCALBRANCHWARMUPFACTOR;
      feasitn++;
      if (feasitn > maxfeasitn) {
	wutermination = true;
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
	cout << "LocalBranchSolver: preprocessing (last): f* = " 
	     << f << ", LSret = " << ret << ", feasible = ";
	if (isfeas) {
	  cout << "true";
	} else {
	  cout << "false";
	}
	cout << " (tol=" << feastolerance << ")" << endl;
      }
    }
    if (!Quiet) {
      cout << "LocalBranchSolver: starting Local Branching proper\n";
    }
    for(int i = 0; i < NumberOfVariables; i++) {
      xstar[i] = x[i];
      // fstar = f;
      fstar = ROSEINFINITY; // don't keep this as a valid solution
    }
  }

  // do it
  if (SolutionMethod == "vns") {
    // VNS on k
    bool termination = false;
    while(!termination) {
      for(samplecount = 0; samplecount < Samples; samplecount++) {
	ret = LocalBranching(k, fstar, xstar, f, x, 0);
	if (ret == 0) {
	  // feasible solution
	  if (f < fstar) {
	    if (!Quiet) {
	      cout << "LocalBranchSolver: improved f* = " << f
		   << " with k = " << k << ", sample " << samplecount + 1
		   << ", LSret = " << ret << endl;
	    }
	    fstar = f;
	    for(int j = 0; j < NumberOfVariables; j++) {
	      xstar[j] = x[j];
	    }
	    k = Kmin;
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
	    cout << "LocalBranchSolver: running time exceeds limit, exiting\n";
	  }
	}
      }
      if (k == Kmax) {
	termination = true;
	if (!Quiet) {
	  cout << "LocalBranchSolver: k exceeds Kmax, exiting\n";
	}
      }
      // increase neighbourhood
      k++;
    }

  } else {
    // do one local branching iteration with a fixed k 
    ret = LocalBranching(k, fstar, xstar, f, x, 0);
  }

  // final stopwatch
  stoptime = getcputime(tend);
  usertimelap = tend[0] - tstart[0];
  if (!Quiet) {
    cout << "LocalBranchSolver: used " << LocalSolver->GetNumberOfCuts()
	 << " local branching cuts, max recursive rank = "
	 << MaxRank << endl;
    cout << "LocalBranchSolver: finished search after " << usertimelap
         << " seconds of user CPU, f* = " << fstar << endl;
  }
  
  // after solution: set optimal values in problem
  IsSolved = true;
  for(int i = 0; i < NumberOfVariables; i++) {
    InProb->SetCurrentVariableValueLI(i + 1, xstar[i]);
  }
  double constrdiscr = 0;
  double vardiscr = 0;
  bool constrfeas = false;
  bool varfeas = false;
  IsFeasible = InProb->TestConstraintsFeasibility(Epsilon, constrdiscr);
  constrfeas = IsFeasible;
  IsFeasible = InProb->TestVariablesFeasibility(Epsilon, vardiscr);
  varfeas = IsFeasible;
  for(int i = 0; i < NumberOfVariables; i++) {
    InProb->SetOptimalVariableValueLI(i + 1, xstar[i]);
  }
  OptimalObjVal = InProb->EvalObj(FIRSTOBJ);
  if (!Quiet) {
    cout << "LocalBranchSolver: best solution f* = " << OptimalObjVal << endl;
  }
  InProb->SetOptimalObjectiveValue(FIRSTOBJ, OptimalObjVal);
  InProb->SetSolved(true);    
  if (constrfeas && varfeas) {
    InProb->SetFeasible(1);
  } else {
    if (!Quiet) {
      if (!constrfeas) {
	cout << "LocalBranchSolver: infeasible by at least " 
	     << constrdiscr << " in constraints\n";
      }
      if (!varfeas) {
	cout << "LocalBranchSolver: infeasible by at least " 
	     << vardiscr << " in integrality/ranges\n";
      }
    }
    InProb->SetFeasible(-1);
  }

  ParameterBlob.SetBoolParameter("Quiet", savequiet);
  InProb->ReplaceParams(ParameterBlob);

  return ret;
}

void LocalBranchSolver::SetOptimizationDirection(int theoptdir) {
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

void LocalBranchSolver::GetSolution(map<int,double>& objfunval, 
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

void LocalBranchSolver::GetSolutionLI(vector<double>& objfunval, 
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

void LocalBranchSolver::GetSolutionLI(double& objfunval, 
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

void LocalBranchSolver::SetStartingPoint(int varID, double sp) {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

void LocalBranchSolver::SetStartingPointLI(int localindex, double sp) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

double LocalBranchSolver::GetStartingPoint(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

double LocalBranchSolver::GetStartingPointLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

void LocalBranchSolver::SetVariableLB(int varID, double LB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableLBLI(localindex, LB);
}

double LocalBranchSolver::GetVariableLB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void LocalBranchSolver::SetVariableUB(int varID, double UB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableUBLI(localindex, UB);
}

double LocalBranchSolver::GetVariableUB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}

void LocalBranchSolver::SetVariableLBLI(int localindex, double LB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vlb[localindex - 1] = LB;
}

double LocalBranchSolver::GetVariableLBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void LocalBranchSolver::SetVariableUBLI(int localindex, double UB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vub[localindex - 1] = UB;
}

double LocalBranchSolver::GetVariableUBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}

void LocalBranchSolver::GetRandomPointInNeighbourhood(std::vector<double>& x,
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
	  // x value hasn't changed, change it artificially
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
