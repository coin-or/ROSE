/*
** Name:       solver_Rfbbtfp.cxx
** Author:     Leo Liberti
** Source:     GNU C++
** Purpose:    Solver Rfbbtfp implementation 
** License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
** History:    100524 work started
*/

#include <cstdio>
#include <cstring>
#include <cassert>
#include <sys/time.h>
#include <fstream>
#include "newsolver.h"
#include "solver.h"
#include "solver_Rfbbtfp.h"
#include "utils.h"

#define FIRSTOBJ 1
#define EPSILON 1e-4

// SNOPT6SOLVER CLASS METHODS

RfbbtfpSolver::RfbbtfpSolver() {
  TheName = "Rfbbtfp";
  NumberOfVariables = 0;
  NumberOfConstraints = 0;
  NumberOfObjectives = 0;
  IsSolved = false;
  OptDir = Minimization;
  OptDirCoeff = 1;
  ManualOptDir = false;
  OptimalObjVal = ROSEINFINITY;
  Epsilon = EPSILON;
  LocalSolverName = "glpk";
  MaxRunningTime = 0;
  OutFile = "Rfbbtfp_out.dat";
}

RfbbtfpSolver::~RfbbtfpSolver() {
}

void RfbbtfpSolver::SetProblem(Problem* p) {
  InProb = p; 
  InitProb = false;
}

bool RfbbtfpSolver::CanSolve(int problemtype) {
  // fill in here with information about solver capabilities
}

void RfbbtfpSolver::Initialize(bool force = false) {

  if (InitProb) {
    ;
  } else {
    InitProb = true;
  }

  // check the input problem is set to something
  if (InProb == NULL) {
    cerr << "RfbbtfpSolver::Initialize(): error: InProb is NULL\n";
    assert(InProb != NULL);
  }

  NumberOfVariables = InProb->GetNumberOfVariables();
  NumberOfConstraints = InProb->GetNumberOfConstraints();
  NumberOfObjectives = InProb->GetNumberOfObjectives();

  // parameters
  ReplaceParams(InProb->GetParamsRef());
  int np = ParameterBlob.GetNumberOfParameters();
  for(int i = 1; i <= np; i++) {
    if (ParameterBlob.GetParameterName(i) == "RfbbtfpLocalSolver") {
      LocalSolverName = ParameterBlob.GetParameterStringValue(i);
    } else if (ParameterBlob.GetParameterName(i) == "RfbbtfpMaxTime") {
      MaxRunningTime = (double) ParameterBlob.GetParameterIntValue(i);
      assert(MaxRunningTime >= 0);
    } else if (ParameterBlob.GetParameterName(i) == "RfbbtfpEpsilon") {
      Epsilon = ParameterBlob.GetParameterDoubleValue(i);
      assert(Epsilon > 0);
    } else if (ParameterBlob.GetParameterName(i) == "RfbbtfpQuiet" ||
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

int RfbbtfpSolver::Solve(bool reinitialize) {

  int ret = 0;
  if (reinitialize || !InitProb) {
    Initialize();
  }

  // do it
  int n = NumberOfVariables;
  int m = NumberOfConstraints;
  int counter = 0;

  double bigU = 1 / Epsilon;

  // read linear constraint expressions with bounds
  vector<Expression*> eptr;
  vector<double> elb;
  vector<double> eub;
  for(int i = 0; i < m; i++) {
    Expression& e = InProb->GetConstraintLI(i + 1)->Function;
    if (e->IsLinear()) {
      eptr.push_back(&e);
      elb.push_back(InProb->GetConstraintLI(i + 1)->LB);
      eub.push_back(InProb->GetConstraintLI(i + 1)->UB);
      if (elb[counter] < -ROSEINFINITY) {
	elb[counter] = -ROSEINFINITY;
      } 
      if (eub[counter] > ROSEINFINITY) {
	eub[counter] = ROSEINFINITY;
      }
      counter++;
    }
  }
  m = counter;

  // read variable bounds and cap at bigU
  map<int,double> xL;
  map<int,double> xU;
  for(int j = 1; j <= n; j++) {
    xL[j] = vlb[j - 1];
    if (xL[j] < -bigU) {
      xL[j] = -bigU;
    } else if (xL[j] > bigU) {
      xL[j] = bigU;
    }
    xU[j] = vub[j - 1];
    if (xU[j] > bigU) {
      xU[j] = bigU;
    } else if (xU[j] < -bigU) {
      xU[j] = -bigU;
    }
  }  

  map<int,double> XL(xL);
  map<int,double> XU(xU);
  map<int,double> ZL(xL);
  map<int,double> ZU(xU);

  ofstream out(OutFile.c_str());

  out << "## .dat instance for fbbt.mod, produced by ROSE ###" << endl;
  out << "## WARNING: read (x) as (1*x)" << endl;
  for (int i = 0; i < m; i++) {
    out << "# " 
         << elb[i] << " <= " << (*eptr[i])->ToString() << " <= " << eub[i] 
	 << endl; 
  }
  for(int j = 1; j <= n; j++) {
    out << "#  x" << j << " in [" << ZL[j] << "," << ZU[j] << "]" << endl;
  }
  out << "###################################################" << endl;

  out.precision(6);
  out << fixed;

  // print out AMPL .dat instance
  vector<std::pair<int,int> > edgeList;
  map<int,int> nodeType; 
  map<int,int> varindex;
  map<int,double> coeff;
  vector<int> nNodes;
  map<int,pair<double,double> > bnds;
  out << "param bigU := " << bigU << ";" << endl;
  out << "param n := " << n << ";" << endl;
  out << "param m := " << m << ";" << endl;
  for(int i = 0; i < m; i++) {
    out << "set V[" << i+1 << "] := ";
    nNodes.push_back
      ((*eptr[i])->TreeData(1, edgeList, nodeType, varindex, coeff, bnds));
    for(int v = 1; v <= nNodes[i] - 1; v++) {
      out << v << ", ";
    }
    out << nNodes[i] << ";" << endl;
    out << "set A[" << i+1 << "] := ";
    int a;
    for(a = 0; a < (int) edgeList.size(); a++) {
      out << "(" << edgeList[a].first << "," << edgeList[a].second << ")";
      if (a < (int) edgeList.size() - 1) {
        out << ", ";
      }
    }
    out << ";" << endl;
  }

  out << "param lambda := " << endl;
  for(int i = 0; i < m; i++) {
    (*eptr[i])->TreeData(1, edgeList, nodeType, varindex, coeff, bnds);
    for(int v = 1; v <= nNodes[i]; v++) {
      out << "  " << i+1 << " " << v << "  " << nodeType[v]+1 << endl;
    }
  }
  out << ";" << endl;

  out << "param index := " << endl;
  for(int i = 0; i < m; i++) {
    (*eptr[i])->TreeData(1, edgeList, nodeType, varindex, coeff, bnds);
    for(int v = 0; v < nNodes[i]; v++) {
      out << "  " << i+1 << " " << v+1 << "  " << varindex[v+1] << endl;
    }
  }
  out << ";" << endl;

  out << "param a := " << endl;
  for(int i = 0; i < m; i++) {
    (*eptr[i])->TreeData(1, edgeList, nodeType, varindex, coeff, bnds);
    for(int v = 0; v < nNodes[i]; v++) {
      out << "  " << i+1 << " " << v+1 << "  " << coeff[v+1] << endl;
    }
  }
  out << ";" << endl;

  out << "param r := " << endl;
  for(int i = 0; i < m; i++) {
    out << "  " << i+1 << "  1" << endl;
  }
  out << ";" << endl;
  
  out << "param gL := " << endl;
  for(int i = 0; i < m; i++) {
    out << "  " << i+1 << "  " << elb[i] << endl;
  }
  out << ";" << endl;

  out << "param gU := " << endl;
  for(int i = 0; i < m; i++) {
    out << "  " << i+1 << "  " << eub[i] << endl;
  }
  out << ";" << endl;

  out << "param x0L := " << endl;
  for(int j = 1; j <= n; j++) {
    out << "  " << j << "  " << xL[j] << endl;
  }
  out << ";" << endl;

  out << "param x0U := " << endl;
  for(int j = 1; j <= n; j++) {
    out << "  " << j << "  " << xU[j] << endl;
  }
  out << ";" << endl;

  out << "param y00L := " << endl;
  for(int i = 0; i < m; i++) {
    (*eptr[i])->TreeData(1, edgeList, nodeType, varindex, coeff, bnds);
    for(int v = 0; v < nNodes[i]; v++) {
       if (nodeType[v+1]+1 == 19) {
         bnds[v+1].first = xL[varindex[v+1]];
         bnds[v+1].second = xU[varindex[v+1]];
       }
     }
     for(int v = 0; v < nNodes[i]; v++) {
       int noderel[2]; int k = 0;
       if (nodeType[v+1]+1 == 3) {
         for(int a = 0; a < (int) edgeList.size(); a++) {
            if (edgeList[a].first == v+1) {
              k++;
              noderel[k] = edgeList[a].second;
            }
         }
         bilinearprodmkrange(bnds[noderel[1]].first, bnds[noderel[1]].second, 
                             bnds[noderel[2]].first, bnds[noderel[2]].second,
                             &bnds[v+1].first, &bnds[v+1].second);
       }
     }
     for(int v = 0; v < nNodes[i]; v++) {
       int noderel[2]; int k = 0;
       if (nodeType[v+1]+1 == 1) {
         for(int a = 0; a < (int) edgeList.size(); a++) {
            if (edgeList[a].first == v+1) {
              k++;
              noderel[k] = edgeList[a].second;
            }
         }
         bnds[v+1].first = bnds[noderel[1]].first+bnds[noderel[2]].first;
         bnds[v+1].second = bnds[noderel[1]].second+bnds[noderel[2]].second;
       }
     }
    for(int v = 0; v < nNodes[i]; v++) {
      out << "  " << i+1 << " " << v+1 << "  " << bnds[v+1].first << endl;
    }
  }
  out << ";" << endl;
  
  out << "param y00U := " << endl;
  for(int i = 0; i < m; i++) {
    (*eptr[i])->TreeData(1, edgeList, nodeType, varindex, coeff, bnds);
    for(int v = 0; v < nNodes[i]; v++) {
       if (nodeType[v+1]+1 == 19) {
         bnds[v+1].first = xL[varindex[v+1]];
         bnds[v+1].second = xU[varindex[v+1]];
       }
     }
     for(int v = 0; v < nNodes[i]; v++) {
       int noderel[2]; int k = 0;
       if (nodeType[v+1]+1 == 3) {
         for(int a = 0; a < (int) edgeList.size(); a++) {
            if (edgeList[a].first == v+1) {
              k++;
              noderel[k] = edgeList[a].second;
            }
         }
         bilinearprodmkrange(bnds[noderel[1]].first, bnds[noderel[1]].second,
                             bnds[noderel[2]].first, bnds[noderel[2]].second,
                             &bnds[v+1].first, &bnds[v+1].second);
       }
     }
     for(int v = 0; v < nNodes[i]; v++) {
       int noderel[2]; int k = 0;
       if (nodeType[v+1]+1 == 1) {
         for(int a = 0; a < (int) edgeList.size(); a++) {
            if (edgeList[a].first == v+1) {
              k++;
              noderel[k] = edgeList[a].second;
            }
         }
         bnds[v+1].first = bnds[noderel[1]].first+bnds[noderel[2]].first;
         bnds[v+1].second = bnds[noderel[1]].second+bnds[noderel[2]].second;
       }
     }

    for(int v = 0; v < nNodes[i]; v++) {
      out << "  " << i+1 << " " << v+1 << "  " << bnds[v+1].second << endl;
    }
  }
  out << ";" << endl;

  out << "## FBBT with epsilon = " << Epsilon << " gives: ########" << endl;
  out << "# FBBT start" << endl;
  for(int j = 1; j <= n; j++) {
    out << "#  x" << j << " in [" << ZL[j] << "," << ZU[j] << "]" << endl;
  }

  // timings
  double tstart[2], tend[2];
  double ttot1, ttot2, tpartial;
  ttot1 = getcputime(tstart);

  // FBBT
  int index;
  int iteration = 1;
  map<int,double>::iterator mi;
  double fbbt_discrepancy = 2*Epsilon;
  while(fbbt_discrepancy > Epsilon) {
    // FBBT Up/Down iteration
    for (int h = 0; h < m; h++) {
      (*eptr[h])->FBBTUpDown(XL, XU, elb[h], eub[h]);
    }
    // compute fbbt_discrepancy
    fbbt_discrepancy = 0;
    for(mi = XL.begin(); mi != XL.end(); mi++) {
      index = mi->first;
      fbbt_discrepancy += fabs(mi->second - ZL[index]);
    } 
    for(mi = XU.begin(); mi != XU.end(); mi++) {
      index = mi->first;
      fbbt_discrepancy += fabs(ZU[index] - mi->second);
    } 
    // timings
    ttot2 = getcputime(tend);
    tpartial = tend[0] - tstart[0];
    // print out current intervals
    if (fbbt_discrepancy > Epsilon) {
      out << "# FBBT iteration " << iteration << ": userCPU = " << tpartial
           << ", discrepancy = " << fbbt_discrepancy << endl;
      for(index = 1; index <= n; index++) {
        out << "#  x" << index << " in [" << XL[index] << "," << XU[index] 
             << "]" << endl;
      }
      // update Z
      for(mi = XL.begin(); mi != XL.end(); mi++) {
        index = mi->first;
        ZL[index] = mi->second;
      }
      for(mi = XU.begin(); mi != XU.end(); mi++) {
        index = mi->first;
        ZU[index] = mi->second;
      }
    }
    iteration++;
  }

  // compute final width
  double totalsum = 0;
  for(index = 1; index <= n; index++) {
    totalsum += XU[index] - XL[index];
  }

  // compute final CPU time
  ttot2 = getcputime(tend);
  tpartial = tend[0] - tstart[0];
  out << "# sum of all variable intervals = " << totalsum << endl;
  out << "# final discrepancy (termination when <= " << Epsilon << ") = "
      << fbbt_discrepancy << endl;
  out << "# Total CPU time = " << ttot2 << ", user CPU time = " << tpartial
       << endl;
  out << "###################################################" << endl;

  out.close();

  if (!Quiet) {
    cout << "# sum of all variable intervals = " << totalsum << endl;
    cout << "# Total CPU time = " << ttot2 << ", user CPU time = " << tpartial
	 << endl;
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
    cout << "RfbbtfpSolver: best solution f* = " << OptimalObjVal << endl;
  }
  InProb->SetOptimalObjectiveValue(FIRSTOBJ, OptimalObjVal);
  InProb->SetSolved(true);    
  if (IsFeasible) {
    InProb->SetFeasible(1);
  } else {
    if (!Quiet) {
      cout << "RfbbtfpSolver: this solution is infeasible by at least " 
           << discrepancy << endl;
    }
    InProb->SetFeasible(-1);
  }

  return ret;
}

void RfbbtfpSolver::SetOptimizationDirection(int theoptdir) {
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

void RfbbtfpSolver::GetSolution(map<int,double>& objfunval, 
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

void RfbbtfpSolver::GetSolutionLI(vector<double>& objfunval, 
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

void RfbbtfpSolver::GetSolutionLI(double& objfunval, 
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


void RfbbtfpSolver::SetStartingPoint(int varID, double sp) {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

void RfbbtfpSolver::SetStartingPointLI(int localindex, double sp) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

double RfbbtfpSolver::GetStartingPoint(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

double RfbbtfpSolver::GetStartingPointLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

void RfbbtfpSolver::SetVariableLB(int varID, double LB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableLBLI(localindex, LB);
}

double RfbbtfpSolver::GetVariableLB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void RfbbtfpSolver::SetVariableUB(int varID, double UB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableUBLI(localindex, UB);
}

double RfbbtfpSolver::GetVariableUB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}

void RfbbtfpSolver::SetVariableLBLI(int localindex, double LB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vlb[localindex - 1] = LB;
}

double RfbbtfpSolver::GetVariableLBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void RfbbtfpSolver::SetVariableUBLI(int localindex, double UB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vub[localindex - 1] = UB;
}

double RfbbtfpSolver::GetVariableUBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}
