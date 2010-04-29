/*
** Name:       solver_tabu.cxx
** Author:     Leo Liberti
** Source:     GNU C++
** Purpose:    GO tabu search solver - implementation (LETSGO)
** License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
** History:    050429 work started
*/

#include <deque>
#include <cstdio>
#include <cstring>
#include <cassert>
#include <sys/time.h>
#include "newsolver.h"
#include "solver.h"
#include "solver_tabu.h"
#include "utils.h"

#define FIRSTOBJ 1
#define EPSILON 1e-4
#define DISTURBEPSILON 1e-1
#define LAGMULTEPSILON 1e-3
#define STPTDISTFACTOR1 1.5
#define STPTDISTFACTOR2 2
#define MAXLOCALSOLUTIONTRIES 10

TabuSolver::TabuSolver() {
  LocalSolverName = "snopt";
  LocalSolver = NULL;
  NumberOfVariables = 0;
  NumberOfConstraints = 0;
  IsSolved = false;
  OptDir = Minimization;
  OptDirCoeff = 1;
  ManualOptDir = false;
  OptimalObjVal = ROSEINFINITY;
  MaxLocalSolutions = 100;
  IsFeasible = false;
  MaxRunningTime = 0; // disabled
  ImprovingTolerance = EPSILON;
  MaxNonImprovingIterations = 10;
  MaxListLength = 5;
  MinRadius = 1e-1;
  MaxRadius = 1;
  int NumberOfIntervals = 10;
  DeltaRadius = (MaxRadius - MinRadius) / NumberOfIntervals;
  Quiet = false;
  TheName = "tabu";
  Epsilon = EPSILON;
  MainSolver = "";
  AmplFlag = false;
}

TabuSolver::~TabuSolver() {
  if (LocalSolver) {
    DeleteSolver(LocalSolver);
  }
}

void TabuSolver::SetProblem(Problem* p) {
  InProb = p; 
  InitProb = false;
}

bool TabuSolver::CanSolve(int problemtype) {
  if (ROSE_LP || ROSE_NLP) {
    return true;
  } else {
    return false;
  }
}

void TabuSolver::Initialize(bool force = false) {

  if (!InitProb || force) {
    InitProb = true;

    // initialize randomizer
    struct timeval theTV;
    struct timezone theTZ;
    gettimeofday(&theTV, &theTZ);
    srandom(theTV.tv_usec);

    // check the input problem is set to something
    if (InProb == NULL) {
      cerr << "TabuSolver::Initialize(): error: InProb is NULL\n";
      assert(InProb != NULL);
    }

    NumberOfVariables = InProb->GetNumberOfVariables();
    NumberOfConstraints = InProb->GetNumberOfConstraints();

    // maximum number of stored local solutions, corresponds to maxaddedcuts
    MaxLocalSolutions = max(100, 10 * NumberOfVariables);
  
    // parameters
    ReplaceParams(InProb->GetParamsRef());
    int np = ParameterBlob.GetNumberOfParameters();
    int tp;
    bool bpv;
    int ipv;
    double dpv;
    string spv;
    for(int i = 1; i <= np; i++) {
      if (ParameterBlob.GetParameterName(i) == "TabuLocalSolver") {
	LocalSolverName = ParameterBlob.GetParameterStringValue(i);
      } else if (ParameterBlob.GetParameterName(i) == "TabuMaxTime") {
	MaxRunningTime = (double) ParameterBlob.GetParameterIntValue(i);
	assert(MaxRunningTime >= 0);
      } else if (ParameterBlob.GetParameterName(i)=="TabuMaxLocalSolutions") {
	MaxLocalSolutions = ParameterBlob.GetParameterIntValue(i);
	assert(MaxLocalSolutions > 0);
      } else if (ParameterBlob.GetParameterName(i)=="TabuImprovingTolerance") {
	ImprovingTolerance = ParameterBlob.GetParameterDoubleValue(i);
	assert(ImprovingTolerance >= 0);
      } else if (ParameterBlob.GetParameterName(i) == "TabuMaxNonImproving") {
	MaxNonImprovingIterations = ParameterBlob.GetParameterIntValue(i);
	assert(MaxNonImprovingIterations >= 0);
      } else if (ParameterBlob.GetParameterName(i) == "TabuMaxListLength") {
	MaxListLength = ParameterBlob.GetParameterIntValue(i);
	assert(MaxListLength > 0);
      } else if (ParameterBlob.GetParameterName(i) == "TabuMinRadius") {
	MinRadius = ParameterBlob.GetParameterDoubleValue(i);
	assert(MinRadius >= 0);
      } else if (ParameterBlob.GetParameterName(i) == "TabuMaxRadius") {
	MaxRadius = ParameterBlob.GetParameterDoubleValue(i);
	assert(MaxRadius >= 0);
      } else if (ParameterBlob.GetParameterName(i) == "TabuDeltaRadius") {
	DeltaRadius = ParameterBlob.GetParameterDoubleValue(i);
      } else if (ParameterBlob.GetParameterName(i) == "TabuEpsilon") {
	Epsilon = ParameterBlob.GetParameterDoubleValue(i);
      } else if (ParameterBlob.GetParameterName(i) == "TabuQuiet" ||
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

    assert(MinRadius <= MaxRadius);

    // create and configure local solver
    LocalSolver = NewSolver(LocalSolverName);
    LocalSolver->SetMaxNumberOfCuts(MaxLocalSolutions);
    LocalSolver->SetProblem(InProb);
    LocalSolver->ReplaceParams(ParameterBlob);
    LocalSolver->Initialize(false);    

    // current solution
    for(int i = 0; i < NumberOfVariables; i++) {
      x.push_back(InProb->GetStartingPointLI(i + 1));
      xstar.push_back(x[i]);
      lb.push_back(InProb->GetVariableLI(i + 1)->LB);
      ub.push_back(InProb->GetVariableLI(i + 1)->UB);
      StartingPoint.push_back(x[i]);
    }
    f = ROSEINFINITY;
    fstar = ROSEINFINITY;
  }
}

int TabuSolver::Solve(bool reinitialize) {

  int ret = 0;
  if (reinitialize || !InitProb) {
    Initialize(reinitialize);
  }

  // re-initialize some stuff for warm starts
  fstar = ROSEINFINITY;
  for(int i = 0; i < NumberOfVariables; i++) {
    lb[i] = InProb->GetVariableLI(i + 1)->LB;
    ub[i] = InProb->GetVariableLI(i + 1)->UB;
  }

  // do it!
  bool termination = false;
  
  // timings
  double tstart[2], tend[2];
  double starttime, stoptime, usertimelap;
  starttime = getcputime(tstart);

  // outer iteration counter
  int outcounter = 0;
  
  // non improving iterations counter
  int nonimproving = 0;

  // tabu list on the cut indices
  deque<int> tabulist;

  // remember all pools seen up to now: map cut index to center, radius
  map<int,pair<vector<double>,double> > pool;
  map<int,pair<vector<double>,double> >::iterator poolit;
  pair<vector<double>,double> poolitem;

  // the local maximum
  vector<double> xmax(NumberOfVariables, 0);  
  double fmax;

  // initial pool radius
  double radius = DeltaRadius < 0 ? MaxRadius : MinRadius;
  // current radius
  double currentradius = radius;

  // current cut to add
  Expression e;
  double constcut;

  // etc
  int tlsize;
  int ci;
  double clm;
  bool radiuschange;
  double thedistance;
  vector<double> d(NumberOfVariables, 0.0);
  double startingpointdistance = radius;

  // outer loop: repeat until termination condition
  while(!termination) {

    if (!Quiet) {
      stoptime = getcputime(tend);  
      usertimelap = tend[0] - tstart[0];
      cout << "TabuSolver: starting outer iteration " << outcounter 
	   << ", elapsed user time " << usertimelap << endl;
    }

    // 1. find random point in hypercube, place it in x[] and localsolver
    RandomLocalSolverStartingPoint();
      
    // inner loop: TS until cuts can be added, then have to re-initialize
    int incounter = 0;
    while(incounter < MaxLocalSolutions && !termination) {
      
      if (!Quiet) {
	stoptime = getcputime(tend);  
	usertimelap = tend[0] - tstart[0];
	cout << "TabuSolver:   starting inner iteration " << incounter 
	     << ", elapsed user time " << usertimelap << endl;
      }

      // 2. local solution (min). 
      bool localsolutionok = false;
      int localsolutioncounter = 0;
      thedistance = ROSEINFINITY;
      while(!localsolutionok) {
	if (!Quiet) {
	  cout << "TabuSolver:     local solver... ";
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
	  // solution is feasible, find values
	  LocalSolver->GetSolutionLI(f, x);
	  localsolutionok = true;
	  IsFeasible = true;
	  thedistance = L2distance(StartingPoint, x);
	  if (!Quiet) {
	    cout << "LSret=" << saveret 
		 << "; f=" << f << "; ||x_0-x||=" << thedistance << endl;
	  }
	  break;
	} else {	  
	  // can't find a solution, reduce radius of cuts in tabu list
	  radiuschange = false;
	  tlsize = tabulist.size();
	  for(int k = 0; k < tlsize; k++) {
	    ci = tabulist[k];
	    currentradius = pool[ci].second - DeltaRadius;
	    if (currentradius >= MinRadius && currentradius <= MaxRadius) {
	      // new radius ok, change it
	      pool[ci].second = currentradius;
	      LocalSolver->SetCutLB(ci, LocalSolver->GetCutLB(ci)+DeltaRadius);
	      radiuschange = true;
	    }
	  }
	  if (!Quiet) {
	    cout << "[no sol, LSret=" << ret << "]" << endl;
	  }

	  if (!radiuschange) {
	    // no change in radius, simply change starting point and repeat
	    RandomLocalSolverStartingPoint();
	    if (!Quiet) {
	      cout << "TabuSolver:     initial point: can't update radiuses, "
		   << "generating random point\n";
	    }
	  } else {
	    if (!Quiet) {
	      cout << "TabuSolver:     initial point: updating radiuses" 
		   << endl;
	    }
	  }

	  localsolutioncounter++;
	  if (localsolutioncounter > MAXLOCALSOLUTIONTRIES) {
	    if (!Quiet) {
	      cout << "TabuSolver:     no feasible local solution found, ";
	      cout << "terminating" << endl;
	    }
	    termination = true;
	    break;
	  }
	}
      }

      // count the non-improving iterations
      if ((fstar - f)*OptDirCoeff < ImprovingTolerance) {
	nonimproving++;
      } else {
	nonimproving = 0;
      }
      if (nonimproving > MaxNonImprovingIterations) {
	if (!Quiet) {
	  cout << "TabuSolver:     too many non-improving iterations, "
	       << "terminating";
	  cout << endl;
	}
	termination = true;
	break;
      }

      if ((f - fstar)*OptDirCoeff < 0) {

	// new solution is better than previous best
	fstar = f;
	for(int i = 0; i < NumberOfVariables; i++) {
	  xstar[i] = x[i];
	}
	if (!Quiet) {
	  cout << "TabuSolver:     new best solution, f* = " << fstar << endl;
	}

	// 6. add cut ||x - min|| > radius and pool to list

	// 6a. see if there is a suitable cut already in pool
	int cutindex = -1;
	for(poolit = pool.begin(); poolit != pool.end(); poolit++) {
	  vector<double>& c = poolit->second.first;
	  double& r = poolit->second.second;
	  double l2d = L2distance(x, c);
	  if (l2d <= r) {
	    cutindex = poolit->first;
	    break;
	  }
	}

	if (cutindex == -1) {
	  // not already in solver

	  // 6b. cut is (sum x_i^2 - 2 sum a_ix_i + sum a_i^2) > radius^2
	  //     hence sum (x_i^2 - 2 a_ix_i) > radius^2 - sum a_i^2, 
	  //     where a is the center of the ball (called 'min' in comments 
	  //     above, and in fact 'x', i.e. current solution, in the code)

	  //     Compute constant
	  constcut = radius;
	  for(int i = 0; i < NumberOfVariables; i++) {
	    constcut -= x[i]*x[i];
	  }
	  
	  //     Make expression
	  e->Zero();
	  for(int i = 0; i < NumberOfVariables; i++) {
	    int varid = InProb->GetVariableID(i + 1);
	    // makes a variable with vid varid
	    Expression ex(1, varid, InProb->GetVariableLI(i + 1)->Name);  
	    Expression square = ex * ex; // makes a square
	    Expression linear;           
	    linear.SetToCopyOf(ex);      // makes a copy of ex
	    linear->SetCoeff(-2*x[i]);   // makes -2*a*ex
	    e = e + square + linear;     // add these terms to expression
	  }
	  Simplify(&e);

	  //     Add cut to solver
	  cutindex = LocalSolver->AddCut(e, constcut, ROSEINFINITY);
	  
	  //     Add cut to pool
	  poolitem.first = x;
	  poolitem.second = radius;
	  pool[cutindex] = poolitem;

	  // starting point outside this tabu zone
	  startingpointdistance *= STPTDISTFACTOR1;

	  if (!Quiet) {
	    cout << "TabuSolver:     adding cut " << cutindex;
#ifdef TABUPRINT
	    cout << ": " << e->ToString() << " >= " << constcut << endl;
#else
	    cout << endl;
#endif
	  }

	} else {

	  // found a relevant cut in the pool, make sure it's not
	  // active already (means local solver has failed, for a 
	  // number of numerical instability-related reasons)
	  bool poolcutok = true;
	  for(int i = 0; i < tabulist.size(); i++) {
	    if (tabulist[i] == cutindex) {
	      poolcutok = false;
	      break;
	    }
	  }

	  if (poolcutok) {
	    // a similar cut is already in pool, activate it
	    LocalSolver->EnableCut(cutindex);
	    if (!Quiet) {
	      cout << "TabuSolver:     enabling cut " << cutindex << endl;
	    }
	    startingpointdistance = STPTDISTFACTOR1 * pool[cutindex].second;
	  } else {
	    // current point is tabu: move initial point well outside
	    // this tabu ball
	    startingpointdistance = STPTDISTFACTOR2 * pool[cutindex].second;
	    // don't add this to the tabu list, it's already there
	    cutindex = -1;
	    if (!Quiet) {
	      cout << "TabuSolver:     current solution is tabu, "
		   << "moving next starting point" << endl;
	    }
	  }

	  // this inner iteration doesn't really count, we have space 
	  // in the solver for more added cuts, enabling cuts does
	  // not use memory
	  incounter--; 
	}
	
	// 6c. Add to tabu list
	if (cutindex != -1) {
	  tabulist.push_back(cutindex);
	}

	// 7. if more than maxtabulist cuts, disable oldest
	if (tabulist.size() > MaxListLength) {
	  cutindex = tabulist[0];
	  tabulist.pop_front();
	  LocalSolver->DisableCut(cutindex);
	  if (!Quiet) {
	    cout << "TabuSolver:     disabling cut " << cutindex << endl;
	  }
	}

	// set next starting point at distance radius from current sol
	// find a random (even infeasible) direction
	FindRandomUnitDirection(d);
	for(int i = 0; i < NumberOfVariables; i++) {
	  d[i] = startingpointdistance * d[i];
	  StartingPoint[i] = x[i] + d[i];
	}

      } else {

	// current solution worse than previous one

	tlsize = tabulist.size();

	// if solution is in tabu zone, local solver has failed
	bool xtabuzone = false;
	int ksave;
	vector<double>* c;
	double* r;
	for(int k = 0; k < tlsize; k++) {
	  c = &(pool[tabulist[k]].first);
	  r = &(pool[tabulist[k]].second);
	  if (L2distance(x, *c) <= *r) {
	    xtabuzone = true;
	    ksave = k;
	    break;
	  }
	}

	if (xtabuzone) {

	  // local solver has failed, likeliest reason is that 
	  // it never managed to move out of the current minimum,
	  // set next starting point well away
	  startingpointdistance *= STPTDISTFACTOR2;
	  FindRandomUnitDirection(d);
	  for(int i = 0; i < NumberOfVariables; i++) {
	    d[i] = startingpointdistance * d[i];
	    StartingPoint[i] = (((*c)[i] + x[i]) / 2) + d[i];
	  }

	  if (!Quiet) {
	    cout << "TabuSolver:     current solution is tabu, "
		 << "moving next starting point" << endl;
	  }	  
	  
	} else {

	  // check whether any added cut is active at the current solution
	  // by testing lagrange multipliers
	  clm = 0;
	  ci = 0;
	  radiuschange = false;
	  for(int k = 0; k < tlsize; k++) {
	    ci = tabulist[k];
	    clm = LocalSolver->GetCutLagrangeMultiplier(ci);
	    if (fabs(clm) > LAGMULTEPSILON) {
	      // lagrange mult for cut is nonzero, constr is active
	      currentradius = pool[ci].second + DeltaRadius;
	      if (currentradius >= MinRadius && currentradius <= MaxRadius) {
		// increase its radius in the pool
		pool[ci].second = currentradius;
		// increase its radius in the solver
		LocalSolver->SetCutLB(ci, 
				      LocalSolver->GetCutLB(ci)+DeltaRadius);
		radiuschange = true;
		if (!Quiet) {
		  cout << "TabuSolver:     cut " << ci 
		       << " active at current solution, new radius "
		       << currentradius << endl;
		}
	      }
	    }
	  }
	  
	  if (!radiuschange) {
	    // no cut is active - increase all radiuses in tabu list then
	    for(int k = 0; k < tlsize; k++) {
	      ci = tabulist[k];
	      currentradius = pool[ci].second + DeltaRadius;
	      if (currentradius >= MinRadius && currentradius <= MaxRadius) {
		// increase its radius in the pool
		pool[ci].second = currentradius;
		// increase its radius in the solver
		LocalSolver->SetCutLB(ci, 
				      LocalSolver->GetCutLB(ci)+DeltaRadius);
		radiuschange = true;
		if (!Quiet) {
		  cout << "TabuSolver:     cut " << ci 
		       << " inactive at current solution, new radius " 
		       << currentradius << endl;
		}
	      }
	    }
	  }
	  
	  if (!radiuschange) {
	    // none of the cuts have changed radius - discard solution
	    // and generate another random one
	    RandomLocalSolverStartingPoint();
	    if (!Quiet) {
	      cout << "TabuSolver:     can't update radiuses, current solution"
		   << " discarded" << endl;
	    }
	  } else {
	    if (!Quiet) {
	      cout << "TabuSolver:     moving to current solution" << endl;
	    }
	  }
	}
      }

      // check the time
      if (MaxRunningTime > 0) {
	stoptime = getcputime(tend);  
	usertimelap = tend[0] - tstart[0];
	if (usertimelap > MaxRunningTime) {
	  if (!Quiet) {
	    cout << "TabuSolver:     maximum time exceeded, terminating\n";
	  }
	  termination = true;
	  break;
	}
      }

      incounter++;
    }

    // finished space within the local solver, re-initialize
    if (!Quiet) {
      cout << "TabuSolver:   deleting structures\n";
    }
    DeleteSolver(LocalSolver);
    LocalSolver = NewSolver(LocalSolverName);
    LocalSolver->SetMaxNumberOfCuts(MaxLocalSolutions);
    LocalSolver->SetProblem(InProb);
    LocalSolver->ReplaceParams(ParameterBlob);
    tabulist.erase(tabulist.begin(), tabulist.end());
    pool.erase(pool.begin(), pool.end());
    startingpointdistance = radius;

    outcounter++;
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
    cout << "TabuSolver: best solution f* = " << OptimalObjVal << endl;
  }
  InProb->SetOptimalObjectiveValue(FIRSTOBJ, OptimalObjVal);
  InProb->SetSolved(true);    
  if (constrfeas && varfeas) {
    InProb->SetFeasible(1);
  } else {
    if (!Quiet) {
      if (!constrfeas) {
	cout << "TabuSolver: infeasible by at least " 
	     << constrdiscr << " in constraints\n";
      }
      if (!varfeas) {
	cout << "TabuSolver: infeasible by at least " 
	     << vardiscr << " in integrality/ranges\n";
      }
    }
    InProb->SetFeasible(-1);
  }

  return ret;
} 

void TabuSolver::SetOptimizationDirection(int theoptdir) {
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

void TabuSolver::GetSolution(map<int,double>& objfunval, 
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

void TabuSolver::GetSolutionLI(vector<double>& objfunval, 
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

void TabuSolver::GetSolutionLI(double& objfunval, 
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


void TabuSolver::SetStartingPoint(int varID, double sp) {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

void TabuSolver::SetStartingPointLI(int localindex, double sp) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

double TabuSolver::GetStartingPoint(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

double TabuSolver::GetStartingPointLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

void TabuSolver::RandomLocalSolverStartingPoint(void) {
  for(int i = 0; i < NumberOfVariables; i++) {
    x[i] = boundedrandom(lb[i], ub[i]);
    StartingPoint[i] = x[i];
    LocalSolver->SetStartingPointLI(i + 1, x[i]);
  }
}

void TabuSolver::DisturbLocalSolverStartingPoint(void) {
  for(int i = 0; i < NumberOfVariables; i++) {
    x[i] = LocalSolver->GetStartingPointLI(i + 1);
    x[i] = boundedrandom(x[i] - DISTURBEPSILON, x[i] + DISTURBEPSILON);
    StartingPoint[i] = x[i];
    LocalSolver->SetStartingPointLI(i + 1, x[i]);
  }
}

void TabuSolver::FindRandomUnitDirection(std::vector<double>& d) {
  int n = d.size();
  double thenorm = 0;
  while(thenorm == 0) {
    for(int i = 0; i < n; i++) {
      d[i] = boundedrandom(-1, 1);
    }
    thenorm = L2norm(d);
  }
  for(int i = 0; i < n; i++) {
    d[i] /= thenorm;
  }
}

////////////////// OBLIVION ////////////////////////////

#ifdef OBLIVION

	// 3. find pool radius:
	if (FindRadiusByInverseOptimization) {
	  // find radius by finding distance to nearest maximum 
	  // (stupid and expensive idea)
	  // 3a. invert opt direction
	  if (OptDir == Minimization) {
	    LocalSolver->SetOptimizationDirection(Maximization);
	  } else {
	    LocalSolver->SetOptimizationDirection(Minimization);
	  }
	  // 4. local solution (max) from (min)
	  DisturbLocalSolverStartingPoint();
	  LocalSolver->Solve();
	  LocalSolver->GetSolutionLI(fmax, xmax);
	  // 5. |max - min| = radius of pool of attraction
	  radius = L2distance(x, xmax);
	  // 5a. re-invert direction
	  LocalSolver->SetOptimizationDirection(OptDir);
	  if (radius < MinRadius) {
	    radius = MinRadius;
	  }
	  if (radius > MaxRadius) {
	    radius = MaxRadius;
	  }
	} else {
	  // 3b. use heuristic
	  radius += DeltaRadius;
	  
#ifdef RANDOMRADIUS
	  // random radius (stupid idea)
	  radius = boundedrandom(MinRadius, MaxRadius);
#endif
	}


#endif
