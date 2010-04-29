/*
** Name:       solver_Rvinci.cxx
** Author:     Leo Liberti & Sonia Cafieri
** Source:     GNU C++
** Purpose:    print out VINCI .ine file
** License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
** History:    080503 work started
*/

#include <fstream>
#include <cstdio>
#include <cstring>
#include <cassert>
#include <sys/time.h>
#include "newsolver.h"
#include "solver.h"
#include "solver_Rvinci.h"
#include "utils.h"

#define FIRSTOBJ 1
#define EPSILON 1e-4

// RVINCISOLVER CLASS METHODS

RvinciSolver::RvinciSolver() {
  TheName = "Rvinci";
  NumberOfVariables = 0;
  NumberOfConstraints = 0;
  NumberOfObjectives = 0;
  IsSolved = false;
  OptDir = Minimization;
  OptDirCoeff = 1;
  ManualOptDir = false;
  OptimalObjVal = ROSEINFINITY;
  Epsilon = EPSILON;
  LocalSolverName = "snopt6";
  MaxRunningTime = 0;
  //OutFile = "Rvinci_out.ine";
}

RvinciSolver::~RvinciSolver() {
}

void RvinciSolver::SetProblem(Problem* p) {
  InProb = p; 
  InitProb = false;
}

bool RvinciSolver::CanSolve(int problemtype) {
  // fill in here with information about solver capabilities
}

void RvinciSolver::Initialize(bool force = false) {


  if (InitProb) {
    ;
  } else {
    InitProb = true;
  }

  // check the input problem is set to something
  if (InProb == NULL) {
    cerr << "RvinciSolver::Initialize(): error: InProb is NULL\n";
    assert(InProb != NULL);
  }

  NumberOfVariables = InProb->GetNumberOfVariables();
  NumberOfConstraints = InProb->GetNumberOfConstraints();
  NumberOfObjectives = InProb->GetNumberOfObjectives();

  // parameters
  ReplaceParams(InProb->GetParamsRef());
  int np = ParameterBlob.GetNumberOfParameters();
  for(int i = 1; i <= np; i++) {
    if (ParameterBlob.GetParameterName(i) == "RvinciLocalSolver") {
      LocalSolverName = ParameterBlob.GetParameterStringValue(i);
    } else if (ParameterBlob.GetParameterName(i) == "RvinciMaxTime") {
      MaxRunningTime = (double) ParameterBlob.GetParameterIntValue(i);
      assert(MaxRunningTime >= 0);
    } else if (ParameterBlob.GetParameterName(i) == "RvinciEpsilon") {
      Epsilon = ParameterBlob.GetParameterDoubleValue(i);
      assert(Epsilon > 0);
    } else if (ParameterBlob.GetParameterName(i) == "RvinciQuiet" ||
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



  // constraints  bounds
  cl = new double [NumberOfConstraints];
  cu = new double [NumberOfConstraints];
  for(int i = 0; i < NumberOfConstraints; i++) {
    cl[i] = InProb->GetConstraintLI(i + 1)->LB;
    cu[i] = InProb->GetConstraintLI(i + 1)->UB;
  }


  // get the rhs
  rhs = new double [NumberOfConstraints];
  int sgn [NumberOfConstraints];
  for(int i = 0; i < NumberOfConstraints; i++) {
    if (cl[i] < cu[i]) {
      if (cl[i] <= -ROSEINFINITY) {
	// case g(x) <= cu
	rhs[i] = cu[i];
        sgn[i] = 1;
      } else if (cu[i] >= ROSEINFINITY) {
	// case g(x) >= cl  => -g(x) <= - cl 
	rhs[i] = -cl[i];
        sgn[i] = -1;
      }
    }
  }



  vector<Expression> ConstrLin;
  for(int i = 1; i <= NumberOfConstraints; i++) {
    ConstrLin.push_back(InProb->GetConstraintLI(i)->Function->GetLinearPart());
  }

  // calculate number of nonzeroes in constraint matrix
  int nzmat = 0;
  for(int i = 1; i <= NumberOfConstraints; i++) {
    if (ConstrLin[i - 1]->IsVariable()) {
      nzmat++;
    } else {
      nzmat += ConstrLin[i - 1]->GetSize();
    }
  }

  rowidx = new int [nzmat + 1];
  colidx = new int [nzmat + 1];
  a = new double [nzmat + 1];  

  // get nonzeroes in constraint matrix
  for(int i = 0; i < NumberOfConstraints; i++) {
    //int acnt = 0;
    if (ConstrLin[i]->IsVariable()) {
      rowidx[acnt] = i;
      colidx[acnt] = ConstrLin[i]->GetVarIndex();
      colidx[acnt] = InProb->GetVarLocalIndex(colidx[acnt]);
      colidx[acnt] = colidx[acnt]-1;
      a[acnt] = sgn[i] * ConstrLin[i]->GetCoeff();
      acnt++;
    } else {
      int thecsize = ConstrLin[i]->GetSize();
      for(int j = 0; j < thecsize; j++) {
	rowidx[acnt] = i;
	colidx[acnt] = ConstrLin[i]->GetNode(j)->GetVarIndex();
	colidx[acnt] = InProb->GetVarLocalIndex(colidx[acnt]);
        colidx[acnt] = colidx[acnt]-1;
	a[acnt] = sgn[i]* ConstrLin[i]->GetNode(j)->GetCoeff();
	acnt++;
      }	
    }
  }
/*
   cout << " acnt = " << acnt<< endl;
   cout << " a= " << endl;
   for(int i=0; i < acnt; i++) {
         cout << a[i] << " ";
   }
    cout << " \n";
   cout << " rowidx= " << endl;
   for(int i=0; i < acnt; i++) {
         cout << rowidx[i] << " ";
   }
    cout << " \n";
   cout << " colidx= " << endl;
   for(int i=0; i < acnt; i++) {
         cout << colidx[i] << " ";
   }
    cout << " \n\n";
*/

  // count the number of finite bounds
  nb = 0;
  for(int i = 0; i < NumberOfVariables; i++) {
      if(vlb[i] != -ROSEINFINITY) {
        nb ++;
      }
      if(vub[i] != ROSEINFINITY) {
        nb ++;
      }
  }


}

int RvinciSolver::Solve(bool reinitialize) {

  int ret = 0;
  if (reinitialize || !InitProb) {
    Initialize();
  }

  // do it

  ///ofstream out(OutFile.c_str());
  ofstream out("invinci.ine");
  out << "VINCI .ine file, automatic translation provided by rose" << endl;
  out << "H-representation" << endl;
  out << "begin" << endl;
  out << NumberOfConstraints+nb << " " << NumberOfVariables+1 << " real" << endl ;

  for(int i = 0; i < NumberOfConstraints; i++) {

      out << rhs[i] << " " ;
      for(int j = 0; j < NumberOfVariables; j++) {
          int flag = 0;
          for(int k = 0; k < acnt; k++) {
              if(rowidx[k]==i && colidx[k]==j) {
                 out << -a[k] << " ";
                 flag = 1;
              } 
          }
          if(flag == 0) {  
             out << 0 << " ";
          }
      }
      
      out <<" \n" ;
  }

  for(int i = 0; i < NumberOfVariables; i++) {
      if(vlb[i] != -ROSEINFINITY) {
         out << -vlb[i] << " ";
         for(int j = 0; j < NumberOfVariables; j++) {
             if (j!=i) { 
                 out << 0 << " ";
             } else {
                 out << 1 << " ";
             }
         }
         out <<" \n";
      }
      if(vub[i] != ROSEINFINITY) {
         out << vub[i] << " ";
         for(int j = 0; j < NumberOfVariables; j++) {
             if (j!=i) { 
                 out << 0 << " ";
             } else {
                 out << -1 << " ";
             }
         }
         out <<" \n";
      }
  }


  out << "end" << endl;

  out.close();


/*
    vector<int> constrnodes;
    vector<Expression> Constr;
    for(int i = 1; i <= NumberOfConstraints; i++) {
      Constr.push_back(InProb->GetConstraintLI(i);
      if (Constr[i - 1]->IsVariable()) {
	constrnodes.push_back(1);
      } else {
	constrnodes.push_back(Constr[i - 1]->GetSize());
      }
    }

  for(int i = 0; i < NumberOfConstraints; i++) { 

     out << cl[i] << -1* InProb->Constr(i)->GetCoeff();
//     out << cl[i] << -1* InProb->GetConstraintLI(i+1)->GetCoeff();

  }

*/


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
    cout << "RvinciSolver: best solution f* = " << OptimalObjVal << endl;
  }
  InProb->SetOptimalObjectiveValue(FIRSTOBJ, OptimalObjVal);
  InProb->SetSolved(true);    
  if (IsFeasible) {
    InProb->SetFeasible(1);
  } else {
    if (!Quiet) {
      cout << "RvinciSolver: this solution is infeasible by at least " 
           << discrepancy << endl;
    }
    InProb->SetFeasible(-1);
  }

  return ret;
}

void RvinciSolver::SetOptimizationDirection(int theoptdir) {
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

void RvinciSolver::GetSolution(map<int,double>& objfunval, 
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

void RvinciSolver::GetSolutionLI(vector<double>& objfunval, 
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

void RvinciSolver::GetSolutionLI(double& objfunval, 
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


void RvinciSolver::SetStartingPoint(int varID, double sp) {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

void RvinciSolver::SetStartingPointLI(int localindex, double sp) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

double RvinciSolver::GetStartingPoint(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

double RvinciSolver::GetStartingPointLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

void RvinciSolver::SetVariableLB(int varID, double LB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableLBLI(localindex, LB);
}

double RvinciSolver::GetVariableLB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void RvinciSolver::SetVariableUB(int varID, double UB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableUBLI(localindex, UB);
}

double RvinciSolver::GetVariableUB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}

void RvinciSolver::SetVariableLBLI(int localindex, double LB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vlb[localindex - 1] = LB;
}

double RvinciSolver::GetVariableLBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void RvinciSolver::SetVariableUBLI(int localindex, double UB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vub[localindex - 1] = UB;
}

double RvinciSolver::GetVariableUBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}
