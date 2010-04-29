/*
** Name:       solver_Rprintncvxdiscr.cxx
** Author:     Sonia Cafieri
** Source:     GNU C++
** Purpose:    Solver Rprintncvxdiscr implementation
** License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
** History:    090917 work started
*/

#include <cstdio>
#include <cstring>
#include <cassert>
#include <sys/time.h>
#include <fstream>
#include <sstream>
#include "newsolver.h"
#include "solver.h"
#include "solver_Rprintncvxdiscr.h"
#include "utils.h"

#define FIRSTOBJ 1
#define EPSILON 1e-4

// RPRINTNCVXDISCRSOLVER CLASS METHODS

RprintncvxdiscrSolver::RprintncvxdiscrSolver() {
  TheName = "Rprintncvxdiscr";
  NumberOfVariables = 0;
  NumberOfConstraints = 0;
  NumberOfObjectives = 0;
  IsSolved = false;
  OptDir = Minimization;
  OptDirCoeff = 1;
  ManualOptDir = false;
  OptimalObjVal = ROSEINFINITY;
  Epsilon = EPSILON;
  LocalSolverName = "smith";
  MaxRunningTime = 0;
  OutFile = "Rprintncvxdiscr.out";
}

RprintncvxdiscrSolver::~RprintncvxdiscrSolver() {
}

void RprintncvxdiscrSolver::SetProblem(Problem* p) {
  InProb = p; 
  InitProb = false;
}

bool RprintncvxdiscrSolver::CanSolve(int problemtype) {
  // fill in here with information about solver capabilities
}

void RprintncvxdiscrSolver::Initialize(bool force = false) {

  if (InitProb) {
    ;
  } else {
    InitProb = true;
  }

  // check the input problem is set to something
  if (InProb == NULL) {
    cerr << "RprintncvxdiscrSolver::Initialize(): error: InProb is NULL\n";
    assert(InProb != NULL);
  }

  NumberOfVariables = InProb->GetNumberOfVariables();
  NumberOfConstraints = InProb->GetNumberOfConstraints();
  NumberOfObjectives = InProb->GetNumberOfObjectives();

  // parameters
  ReplaceParams(InProb->GetParamsRef());
  int np = ParameterBlob.GetNumberOfParameters();
  for(int i = 1; i <= np; i++) {
    if (ParameterBlob.GetParameterName(i) == "RprintncvxdiscrLocalSolver") {
      LocalSolverName = ParameterBlob.GetParameterStringValue(i);
    } else if (ParameterBlob.GetParameterName(i) == "RprintncvxdiscrMaxTime") {
      MaxRunningTime = (double) ParameterBlob.GetParameterIntValue(i);
      assert(MaxRunningTime >= 0);
    } else if (ParameterBlob.GetParameterName(i) == "RprintncvxdiscrEpsilon") {
      Epsilon = ParameterBlob.GetParameterDoubleValue(i);
      assert(Epsilon > 0);
    } else if (ParameterBlob.GetParameterName(i) == "RprintncvxdiscrQuiet" ||
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
cout << "integr " << integrality[i-1] << endl;
  }

  // create and configure local solver
  LocalSolver = NewSolver(LocalSolverName);
  LocalSolver->SetProblem(InProb);
  LocalSolver->ReplaceParams(ParameterBlob);

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

int RprintncvxdiscrSolver::Solve(bool reinitialize) {

  int ret = 0;
  if (reinitialize || !InitProb) {
    Initialize();
  }

  // do it
  nvar_orig = NumberOfVariables;
  nconstr_orig = NumberOfConstraints;

  // Standardize the problem using Rsmith
  //LocalSolver->Initialize(false);
  LocalSolver->Solve();

  NumberOfVariables = InProb->GetNumberOfVariables();
  NumberOfConstraints = InProb->GetNumberOfConstraints();

  // re-set variable data
  vlb.erase(vlb.begin(),vlb.end());
  vub.erase(vub.begin(),vub.end());
  for(int i = 0; i < NumberOfVariables; i++) {
    vlb.push_back(InProb->GetVariableLI(i + 1)->LB) ;
    vub.push_back(InProb->GetVariableLI(i + 1)->UB) ;
    x.push_back(0);
    xstar.push_back(0);
    StartingPoint.push_back(x[i]);
  }

  int addedvarindex;
  int size;
  double coeff;
  vector<int> vidx;
  int ncvxdefcon = 0;
  int maxsize = 1;

  ofstream out(OutFile.c_str());

  for(int j = 1; j <= NumberOfConstraints; j++) {

    if (!InProb->GetConstraintLI(j)->Function->IsLinear()) {
      // nonlinear defining constraint
      ncvxdefcon ++;
      size = InProb->GetConstraintLI(j)->Function->GetNode(1)->GetSize();
      if (size > maxsize) maxsize = size;
    }
  }
  out << "let sBB_ncvxdefcon := " << ncvxdefcon <<";" << endl; 
  out << "let sBB_ncvxsize := " << maxsize <<";" << endl; 

  
  ncvxdefcon = 0;
  for(int j = 1; j <= NumberOfConstraints; j++) {

    if (!InProb->GetConstraintLI(j)->Function->IsLinear()) {

      // nonlinear defining constraint
      ncvxdefcon ++;

      size = InProb->GetConstraintLI(j)->Function->GetNode(1)->GetSize();
      int op = InProb->GetConstraintLI(j)->Function->GetNode(1)->GetOpType();
      double expo =
        InProb->GetConstraintLI(j)->Function->GetNode(1)->GetExponent();
      coeff=InProb->GetConstraintLI(j)->Function->GetNode(1)->GetCoeff();

      //Expression e = - InProb->GetConstraintLI(j)->Function->GetLinearPart();
      addedvarindex = InProb->GetConstraintLI(j)->Function->
        GetNode(0)->GetVarIndex();
      
      InProb->GetConstraintLI(j)->Function->
        GetNode(1)->GetVarIndices(vidx);

      // discrepancy between added linearizing variable and original nonconvex expression
      out << " " << endl;
      out << "let sBB_ncvx_discrepancy["<<ncvxdefcon<<"]:= abs(";
      out << "Rcvx_w["<<addedvarindex<<"] - ";

      if(op == PRODUCT) {
        out <<"(";
        for (int k=0; k<size-1; k++) {
           if(vidx[k]<= nvar_orig) {
             if (integrality[vidx[k]-1]) {
               out <<"Rcvx_xInt["<<vidx[k];
             } else {
               out <<"Rcvx_x["<<vidx[k];
             }
           } else {
            out <<"Rcvx_w["<<vidx[k];
           }
           out<<"]*";
        }
        if(vidx[size-1]<= nvar_orig) {
          if (integrality[vidx[size-1]-1]) {
            out <<"Rcvx_xInt["<<vidx[size-1];
          } else {
            out <<"Rcvx_x["<<vidx[size-1];
          }
        } else {
           out <<"Rcvx_w["<<vidx[size-1];
        }
        out <<"])";

      } else if(op == POWER || (op == VAR && (expo!=0 && expo!=1))) {
        if(vidx[0]<= nvar_orig) {
           if (integrality[vidx[0]-1]) {
             out <<"Rcvx_xInt["<<vidx[0]<<"]^"<<expo;
           } else {
             out <<"Rcvx_x["<<vidx[0]<<"]^"<<expo;
           }
        } else {
           out <<"Rcvx_w["<<vidx[0]<<"]^"<<expo;
        }
        size = 1;
      }
      out <<");" << endl;
      vidx.erase(vidx.begin(),vidx.end());
    
 
      // for each nonlinear defining constraint,
      // indices of (original) variables involved in the nonconvex expression 
      // and inverse map (for each variable index, gives its position in the 
      // corresponding sBB_ncvx)
      for (int k=0; k<maxsize; k++) {
        
         out << "let sBB_ncvx["<<ncvxdefcon<<","<<k+1<<"]:= ";
         if (k <= size) {
            out << vidx[k] << ";" << endl;
         } else {
            out << 0 << ";" << endl;
         }
         out << "let sBB_ncvx_inv["<<ncvxdefcon<<",";
         if (k <= size) {
            out << vidx[k] << "]:= ";
         } else {
            out << 0 <<"]:= ";
         }
         out << k+1  << ";" << endl;

      }


    }

  }

  // after solution: set optimal values in problem
  IsSolved = true;
  for(int i = 0; i < NumberOfVariables; i++) {
    InProb->SetCurrentVariableValueLI(i + 1, xstar[i]);
  }

  return ret;
}

void RprintncvxdiscrSolver::SetOptimizationDirection(int theoptdir) {
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

void RprintncvxdiscrSolver::GetSolution(map<int,double>& objfunval, 
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

void RprintncvxdiscrSolver::GetSolutionLI(vector<double>& objfunval, 
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

void RprintncvxdiscrSolver::GetSolutionLI(double& objfunval, 
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


void RprintncvxdiscrSolver::SetStartingPoint(int varID, double sp) {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

void RprintncvxdiscrSolver::SetStartingPointLI(int localindex, double sp) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

double RprintncvxdiscrSolver::GetStartingPoint(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

double RprintncvxdiscrSolver::GetStartingPointLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

void RprintncvxdiscrSolver::SetVariableLB(int varID, double LB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableLBLI(localindex, LB);
}

double RprintncvxdiscrSolver::GetVariableLB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void RprintncvxdiscrSolver::SetVariableUB(int varID, double UB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableUBLI(localindex, UB);
}

double RprintncvxdiscrSolver::GetVariableUB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}

void RprintncvxdiscrSolver::SetVariableLBLI(int localindex, double LB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vlb[localindex - 1] = LB;
}

double RprintncvxdiscrSolver::GetVariableLBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void RprintncvxdiscrSolver::SetVariableUBLI(int localindex, double UB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vub[localindex - 1] = UB;
}

double RprintncvxdiscrSolver::GetVariableUBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}
