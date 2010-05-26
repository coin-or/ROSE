/*
** Name:       solver_Rporta.cxx
** Author:     Leo Liberti & Sonia Cafieri
** Source:     GNU C++
** Purpose:    print out PORTA .ieq file
** License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
** History:    080517 work started
*/

#include <fstream>
#include <cstdio>
#include <cstring>
#include <cassert>
#include <sys/time.h>
#include "newsolver.h"
#include "solver.h"
#include "solver_Rporta.h"
#include <string>
#include <sstream>
#include "utils.h"

#define FIRSTOBJ 1
#define EPSILON 1e-4

// RPORTASOLVER CLASS METHODS

RportaSolver::RportaSolver() {
  TheName = "Rporta";
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
  //OutFile = "Rporta.ieq";
}

RportaSolver::~RportaSolver() {
}

void RportaSolver::SetProblem(Problem* p) {
  InProb = p; 
  InitProb = false;
}

bool RportaSolver::CanSolve(int problemtype) {
  // fill in here with information about solver capabilities
}

void RportaSolver::Initialize(bool force = false) {


  if (InitProb) {
    ;
  } else {
    InitProb = true;
  }

  // check the input problem is set to something
  if (InProb == NULL) {
    cerr << "RportaSolver::Initialize(): error: InProb is NULL\n";
    assert(InProb != NULL);
  }

  NumberOfVariables = InProb->GetNumberOfVariables();
  NumberOfConstraints = InProb->GetNumberOfConstraints();
  NumberOfObjectives = InProb->GetNumberOfObjectives();

  // parameters
  ReplaceParams(InProb->GetParamsRef());
  int np = ParameterBlob.GetNumberOfParameters();
  for(int i = 1; i <= np; i++) {
    if (ParameterBlob.GetParameterName(i) == "Rporta") {
      LocalSolverName = ParameterBlob.GetParameterStringValue(i);
    } else if (ParameterBlob.GetParameterName(i) == "Rporta") {
      MaxRunningTime = (double) ParameterBlob.GetParameterIntValue(i);
      assert(MaxRunningTime >= 0);
    } else if (ParameterBlob.GetParameterName(i) == "Rporta") {
      Epsilon = ParameterBlob.GetParameterDoubleValue(i);
      assert(Epsilon > 0);
    } else if (ParameterBlob.GetParameterName(i) == "Rporta" ||
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

}

int RportaSolver::Solve(bool reinitialize) {

  int ret = 0;
  if (reinitialize || !InitProb) {
    Initialize();
  }

  // do it

  // number of var. to be eliminated by Porta - change if necessary
  int nel = 1;

  ///ofstream out(OutFile.c_str());
  ofstream out("inporta.ieq");
  out << "DIM = ";
  out << NumberOfVariables << endl ;
  out << "\nCOMMENT\n";
  out << "PORTA .ieq file, automatic translation provided by rose" << endl;
  out << "\nELIMINATION_ORDER" << endl;
  //out << "0 0 0 0 0 2 1 " << endl ;
  for (int i = 0; i < NumberOfVariables-nel; i++) {
     out << "0 ";
  }
  for (int i = nel; i > 0; i--) {
     out << i << " " ;
  }
  out << "\n " << endl;
  out << "\nINEQUALITIES_SECTION" << endl;

  for(int i = 0; i < NumberOfConstraints; i++) {
      
      int cnt = 0;
      for(int j = 0; j < NumberOfVariables; j++) {
          for(int k = 0; k < acnt; k++) {
              if(rowidx[k]==i && colidx[k]==j) {
                 cnt ++;
                 if(a[k]>0 && cnt>1) out <<"+";
                 if((int)a[k]==a[k]) {
                    out << a[k] << "x"<< j+1 << " ";
                 } else {
                    int dd = decimal(a[k]);
                    out << a[k]*pow(10,(double)dd) <<"/" << pow(10, (double)dd) <<"x"<< j+1 << " ";
                 }    
              } 
          }
      }
      out <<" <=";
      if((int)rhs[i]==rhs[i]) {
         out << rhs[i];
      } else {
         int dd = decimal(rhs[i]);
         out << rhs[i]*pow(10,(double)dd) <<"/" << pow(10, (double)dd);
      }
      out <<" \n" ;
  }

  for(int i = 0; i < NumberOfVariables; i++) {
      if((int)vlb[i]==vlb[i]) {
        out << "-x" << i+1 << " <= " << -vlb[i] << endl;
      } else {
        int dd = decimal(vlb[i]);
        out << "-x" << i+1 << " <= " << -vlb[i]*pow(10,(double)dd) <<"/" << pow(10, (double)dd) << endl;
      }
      if((int)vub[i]==vub[i]) {
        out << "x" << i+1 << " <= " << vub[i] << endl;
      } else {
        int dd = decimal(vub[i]);
        out << "x" << i+1 << " <= " << vub[i]*pow(10,(double)dd) <<"/" << pow(10, (double)dd) << endl;
      }
  }

  out << "\nEND" << endl;

  out.close();


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
    cout << "Rporta: best solution f* = " << OptimalObjVal << endl;
  }
  InProb->SetOptimalObjectiveValue(FIRSTOBJ, OptimalObjVal);
  InProb->SetSolved(true);    
  if (IsFeasible) {
    InProb->SetFeasible(1);
  } else {
    if (!Quiet) {
      cout << "Rporta: this solution is infeasible by at least " 
           << discrepancy << endl;
    }
    InProb->SetFeasible(-1);
  }

  return ret;
}

void RportaSolver::SetOptimizationDirection(int theoptdir) {
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

void RportaSolver::GetSolution(map<int,double>& objfunval, 
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

void RportaSolver::GetSolutionLI(vector<double>& objfunval, 
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

void RportaSolver::GetSolutionLI(double& objfunval, 
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


void RportaSolver::SetStartingPoint(int varID, double sp) {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

void RportaSolver::SetStartingPointLI(int localindex, double sp) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

double RportaSolver::GetStartingPoint(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

double RportaSolver::GetStartingPointLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

void RportaSolver::SetVariableLB(int varID, double LB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableLBLI(localindex, LB);
}

double RportaSolver::GetVariableLB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void RportaSolver::SetVariableUB(int varID, double UB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableUBLI(localindex, UB);
}

double RportaSolver::GetVariableUB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}

void RportaSolver::SetVariableLBLI(int localindex, double LB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vlb[localindex - 1] = LB;
}

double RportaSolver::GetVariableLBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void RportaSolver::SetVariableUBLI(int localindex, double UB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vub[localindex - 1] = UB;
}

double RportaSolver::GetVariableUBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}



int RportaSolver::decimal(double n) {

  std::string s;
  // convert double n to string s
  { std::ostringstream ss;
    ss << n;
    s = ss.str();
  }

  const char *s1 = s.c_str(); // get const char * representation
  int len = strlen(s1) -1;

  int ii=0;
  while(s[ii]!='.') {
     ii++;
  }

  return len-ii;
}
