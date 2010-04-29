/*
** Name:       solver_RQuarticConvex.cxx
** Author:     Leo Liberti & Sonia Cafieri
** Source:     GNU C++
** Purpose:    Solver RQuarticConvex implementation
** License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
** History:    080709 work started
*/

#include <cstdio>
#include <cstring>
#include <cassert>
#include <sys/time.h>
#include <fstream>
#include "newsolver.h"
#include "solver.h"
#include "solver_RQuarticConvex.h"
#include "utils.h"
#include "Ev3/auxiliary.h"

#define FIRSTOBJ 1
#define EPSILON 1e-4
#define WDEFNAMEBUFSIZE 256

// RQUARTICCONVEX-SOLVER CLASS METHODS

RQuarticConvexSolver::RQuarticConvexSolver() {
  TheName = "RQuarticConvex";
  LocalSolver = NULL;
  NumberOfVariables = 0;
  NumberOfConstraints = 0;
  NumberOfObjectives = 0;
  IsSolved = false;
  OptDir = Minimization;
  OptDirCoeff = 1;
  ManualOptDir = false;
  OptimalObjVal = ROSEINFINITY;
  Epsilon = EPSILON;
  LocalSolverName = "Rsmith";
  MaxRunningTime = 0;
  //OutFile = "RQuarticConvex_out.ros";
  QuarticConvexType = 0;
}

RQuarticConvexSolver::~RQuarticConvexSolver() {
  if (LocalSolver) {
    DeleteSolver(LocalSolver);
  }
}

void RQuarticConvexSolver::SetProblem(Problem* p) {
  InProb = p; 
  InitProb = false;
}

bool RQuarticConvexSolver::CanSolve(int problemtype) {
  // fill in here with information about solver capabilities
}

void RQuarticConvexSolver::Initialize(bool force = false) {

  if (InitProb) {
    ;
  } else {
    InitProb = true;
  }

  // check the input problem is set to something
  if (InProb == NULL) {
    cerr << "RQuarticConvexSolver::Initialize(): error: InProb is NULL\n";
    assert(InProb != NULL);
  }

  NumberOfVariables = InProb->GetNumberOfVariables();
  NumberOfConstraints = InProb->GetNumberOfConstraints();
  NumberOfObjectives = InProb->GetNumberOfObjectives();

  // parameters
  ReplaceParams(InProb->GetParamsRef());
  int np = ParameterBlob.GetNumberOfParameters();
  for(int i = 1; i <= np; i++) {
    if (ParameterBlob.GetParameterName(i) == "RQuarticConvexLocalSolver") {
      LocalSolverName = ParameterBlob.GetParameterStringValue(i);
    } else if (ParameterBlob.GetParameterName(i) == "RQuarticConvexMaxTime") {
      MaxRunningTime = (double) ParameterBlob.GetParameterIntValue(i);
      assert(MaxRunningTime >= 0);
    } else if (ParameterBlob.GetParameterName(i) == "RQuarticConvexEpsilon") {
      Epsilon = ParameterBlob.GetParameterDoubleValue(i);
      assert(Epsilon > 0);
    } else if (ParameterBlob.GetParameterName(i) == "RQuarticConvexQuiet" ||
               ParameterBlob.GetParameterName(i) == "Quiet") {
      Quiet = ParameterBlob.GetParameterBoolValue(i);
    } else if (ParameterBlob.GetParameterName(i) == "RQuarticConvexType") {
      QuarticConvexType = ParameterBlob.GetParameterIntValue(i);
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

}

int RQuarticConvexSolver::Solve(bool reinitialize) {

  int ret = 0;
  if (reinitialize || !InitProb) {
    Initialize();
  }

  // do it
////////////////////
  double addconst = InProb->GetObj1AdditiveConstant();

  int nvar_orig = NumberOfVariables;
  int nconstr_orig = NumberOfConstraints;

  // Standardize the problem using RSMITH

  //cout << "\nSMITH: " << endl;
  LocalSolver->Initialize();
  LocalSolver->Solve();

  NumberOfVariables = InProb->GetNumberOfVariables();
  NumberOfConstraints = InProb->GetNumberOfConstraints();
  int nconstr_smith = NumberOfConstraints;

//  cout << *InProb;

  for(int i = nvar_orig; i < NumberOfVariables; i++) {
  //for(int i = 0; i < NumberOfVariables; i++) {
     vlb.push_back(InProb->GetVariableLI(i + 1)->LB) ;
     vub.push_back(InProb->GetVariableLI(i + 1)->UB) ;
     x.push_back(0);
     xstar.push_back(0);
     StartingPoint.push_back(x[i]);
  }


  // added variable name
  string wn = "w";
  char wdefnbuf[WDEFNAMEBUFSIZE];
  // added variable bounds
  double wlb, wub;
  double wlb1, wub1;
  double wlb2, wub2;
  // added variable counter
  int waddcounter = NumberOfVariables + 1;


  vector<Expression> cc;
  vector<Expression> defcons;
  string cn = "c";
  int v1, v2, v3, v4;

 
  // note: here NumberOfConstraints is that computed after Rsmith

  for(int j = 1; j <= NumberOfConstraints; j++) {

    int op = InProb->GetConstraintLI(j)->Function->GetNode(1)->GetOpType();
    double expo = InProb->GetConstraintLI(j)->Function->GetNode(1)->GetExponent();
    double coeff = InProb->GetConstraintLI(j)->Function->GetNode(1)->GetCoeff();
    Expression e = - InProb->GetConstraintLI(j)->Function->GetLinearPart();
    //if (coeff < 0) e = -e;

    int v1, v2, v3, v4;

    int size = InProb->GetConstraintLI(j)->Function->GetNode(1)->GetSize();

    if(op == PRODUCT && size == 4) {
//      cout << "Product" << endl;

      // find variables involved in quadrilinear terms
      vector<int> vidx;
      Expression xs1(1, 1, "x1"), xs2(1, 2, "x2"), xs3(1, 3, "x3"), xs4(1, 4, "x4");
      Expression schema = xs1*xs2*xs3*xs4;
      InProb->GetConstraintLI(j)->Function->GetVarIndicesInSchema(vidx, schema);

      v1 = vidx[0];
      v2 = vidx[1];
      v3 = vidx[2];
      v4 = vidx[3];
 
      string x1name = InProb->GetVariableLI(v1)->Name; 
      string x2name = InProb->GetVariableLI(v2)->Name; 
      string x3name = InProb->GetVariableLI(v3)->Name; 
      string x4name = InProb->GetVariableLI(v4)->Name; 

      double cf = 1;
      double coeff1 = InProb->GetConstraintLI(j)->Function->GetNode(1)->GetNode(0)->GetCoeff();
      if (coeff != 1 && coeff != -1) cf = cf*coeff;    // case  coeff*(x1*x2*x3*x4)
      if (coeff1 != 1 && coeff1 != -1) cf = cf*coeff1; // case  ((coeff1*x1)*x2*x3*x4)      

      Expression x1(cf, v1, x1name), x2(1, v2, x2name), x3(1, v3, x3name), x4(1, v4, x4name);


      // choice of the relaxation type
      int bilin1 = 0, bilin2 = 0, trilin = 0, trilin2 = 0; 
      
      if (QuarticConvexType == 0) {
	trilin2 = 1;
      } else if (QuarticConvexType == 1) {
	trilin = 1;
      } else if (QuarticConvexType == 2) {
	bilin1 = 1;
      } else if (QuarticConvexType == 3) {
	bilin2 = 1;
      }

#undef sprintf
      sprintf(wdefnbuf, "w%d", waddcounter);
      wn = wdefnbuf;  // name
      Expression w1(1,waddcounter,wn);
      int storeidx;
      if(trilin2 ==1 ) storeidx = waddcounter;

      if(bilin1 == 1 || bilin2 == 1 || trilin2 == 1) {

#undef sprintf
        // add constraint w - x1*x2 = 0
        //sprintf(wdefnbuf, "w%d", waddcounter);
        //wn = wdefnbuf;  // name
        //Expression w1(1,waddcounter,wn);
        InProb->AddConstraint(cn, -w1+x1*x2, 0, 0);
        // add variable w = x1*x2
        w1 = x1*x2;
        bilinearprodmkrange(vlb[vidx[0]-1], vub[vidx[0]-1], vlb[vidx[1]-1], vub[vidx[1]-1], &wlb, &wub); // bounds
        if(cf != 1.) bilinearprodmkrange(wlb, wub, cf, cf, &wlb, &wub); // bounds - to take into account cf
        InProb->AddVariable(wn, false, false, wlb, wub, 0);
        vlb.push_back(wlb) ;
        vub.push_back(wub) ;
        x.push_back(0);
        xstar.push_back(0);
        StartingPoint.push_back(0);
        
        // replace by schema 
        Expression schema1 = x1*x2;
        Expression ww1 = w1->ReplaceBySchema(waddcounter, wn, schema1);

        if(bilin1 == 1) {

          // add constraint w - x3*x4 = 0
          sprintf(wdefnbuf, "w%d", waddcounter+1);
          wn = wdefnbuf;  // name
          Expression w2(1,waddcounter+1,wn);
          InProb->AddConstraint(cn, -w2+x3*x4, 0, 0);
          // add variable w = x3*x4
          w2 = x3*x4;
          bilinearprodmkrange(vlb[vidx[2]-1], vub[vidx[2]-1], vlb[vidx[3]-1], vub[vidx[3]-1], &wlb1, &wub1); // bounds
          InProb->AddVariable(wn, false, false, wlb1, wub1, 0);
          vlb.push_back(wlb1) ;
          vub.push_back(wub1) ;
          x.push_back(0);
          xstar.push_back(0);
          StartingPoint.push_back(0);
 
          // replace by schema 
          Expression schema2 = x3*x4;
          Expression ww2 = w2->ReplaceBySchema(waddcounter+1, wn, schema2);

          // add constraint e - w1*w2 = 0, where e is the variable w added by Rsmith
          Expression wx(1, waddcounter+1, "wx");
          wx = w1*w2;

          //Expression cw = -e + wx;
          Expression cw;
          if(coeff == -1) cw = e + wx;
          else cw = -e + wx;
          InProb->AddConstraint(cn, cw, 0, 0);

          waddcounter = waddcounter + 2;

        } // endif(bilin1=1)

        if(bilin2 == 1) {

          // add constraint w - (x1*x2)*x3 = 0
          sprintf(wdefnbuf, "w%d", waddcounter+1);
          wn = wdefnbuf;  // name
          Expression w3(1,waddcounter+1,wn);
          InProb->AddConstraint(cn, -w3+w1*x3, 0, 0);
          // add variable w = w1*x3
          w3 = w1*x3;
          bilinearprodmkrange(wlb, wub, vlb[vidx[2]-1], vub[vidx[2]-1], &wlb1, &wub1); // bounds
          InProb->AddVariable(wn, false, false, wlb1, wub1, 0);
          vlb.push_back(wlb1) ;
          vub.push_back(wub1) ;
          x.push_back(0);
          xstar.push_back(0);
          StartingPoint.push_back(0);
 
          // replace by schema 
          Expression schema2 = w1*x3;
          Expression ww3 = w3->ReplaceBySchema(waddcounter+1, wn, schema2);

          // add constraint e - w3*x4 = 0, where e is the variable w added by Rsmith
          Expression wx(1, waddcounter+1, "wx");
          wx = w3*x4;

          //Expression cw = -e + wx;
          Expression cw;
          if(coeff == -1) cw = e + wx;
          else cw = -e + wx;
          InProb->AddConstraint(cn, cw, 0, 0);

          waddcounter = waddcounter + 2;

        } // endif(bilin2=1)
  //cout << *InProb;

      } // endif (bilin1=1 || bilin2=1 || trilin2=1)

      if(trilin == 1 || trilin2 == 1) {

        sprintf(wdefnbuf, "w%d", waddcounter);
        wn = wdefnbuf;  // name
        Expression w(1,waddcounter,wn);

      if(trilin == 1) {
#undef sprintf
      // add constraint w - x1*x2*x3 = 0
      //sprintf(wdefnbuf, "w%d", waddcounter);
      //wn = wdefnbuf;  // name
      //Expression w(1,waddcounter,wn);
      InProb->AddConstraint(cn, -w+x1*x2*x3, 0, 0);
      // add variable w = x1*x2*x3
      w = x1*x2*x3;
      bilinearprodmkrange(vlb[vidx[0]-1], vub[vidx[0]-1], vlb[vidx[1]-1], vub[vidx[1]-1], &wlb1, &wub1);
      bilinearprodmkrange(wlb1, wub1, vlb[vidx[2]-1], vub[vidx[2]-1], &wlb, &wub); // bounds
      if(cf != 1.) bilinearprodmkrange(wlb, wub, cf, cf, &wlb, &wub); // bounds - to take into account cf
      InProb->AddVariable(wn, false, false, wlb, wub, 0);
      vlb.push_back(wlb) ;
      vub.push_back(wub) ;
      x.push_back(0);
      xstar.push_back(0);
      StartingPoint.push_back(0);

      // replace by schema 
      Expression schema1 = x1*x2*x3;
      Expression ww = w->ReplaceBySchema(waddcounter, wn, schema1);

      // add constraint e - w*x4 = 0, where e is the variable w added by Rsmith
      Expression wx(1, waddcounter+1, "wx");
      wx = w*x4;

      //Expression cw = -(e - wx);
      Expression cw;
      if(coeff == -1) cw = e + wx;
      else cw = -e + wx;
      InProb->AddConstraint(cn, cw, 0, 0);

      waddcounter++;

  //cout << *InProb;

      }

      if(trilin2 == 1) {
#undef sprintf
      // add constraint e - w1*x3*x4 = 0, where e is the variable w added by Rsmith
      Expression cw;
      if(coeff == -1) cw = e + w1*x3*x4;
      else cw = -e + w1*x3*x4;
      InProb->AddConstraint(cn, cw, 0, 0);

      waddcounter++;

  //cout << *InProb;

      }


      int **ind;
      ind = new int*[6];
      for(int i=0; i<6; i++) {
        ind[i] = new int[6];
      }
     
      int *ibnd; 
      ibnd = new int[3];
      if(trilin == 1) {
         ibnd[0] = v1-1; ibnd[1] = v2-1; ibnd[2] = v3-1;
      } 
      if(trilin2 == 1) {
         ibnd[0] = storeidx-1; ibnd[1] = v3-1; ibnd[2] = v4-1;
      }
 
      // compute the 6 permutations of the 3 variables 
      permutation3(ind,ibnd);

      int i, flag=0, idx=0;
      i = 0;
      while(i < 6 && flag == 0) {
         if(vlb[ind[i][0]] >=0 && vlb[ind[i][1]] >=0 && vlb[ind[i][2]] <=0 && vub[ind[i][2]] >=0) {
            idx = i;   // store the index of the permutation satisfying the condition
            flag = 1;  // this is case 1
         }
         i++; 
      }
      i = 0;
      while(i < 6 && flag == 0) {
         if(vlb[ind[i][0]] >=0 && vlb[ind[i][1]] <=0 && vlb[ind[i][2]] <=0 
            && vub[ind[i][1]] >=0 && vub[ind[i][2]] >=0) {
            idx = i;   // store the index of the permutation satisfying the condition
            flag = 2;  // this is case 2
         }
         i++; 
      }
      i = 0;
      while(i < 6 && flag == 0) {
         if(vlb[ind[i][0]] <=0 && vlb[ind[i][1]] <=0 && vlb[ind[i][2]] <=0 
            && vub[ind[i][0]] >=0 && vub[ind[i][1]] >=0 && vub[ind[i][2]] >=0) {
            idx = i;   // store the index of the permutation satisfying the condition
            flag = 3;  // this is case 3
         }
         i++; 
      }
      i = 0;
      while(i < 6 && flag == 0) {
         if(vlb[ind[i][0]] >=0 && vlb[ind[i][1]] <=0 && vlb[ind[i][2]] <=0 
            && vub[ind[i][1]] >=0 && vub[ind[i][2]] <=0) {
            idx = i;   // store the index of the permutation satisfying the condition
            flag = 4;  // this is case 4
         }
         i++; 
      }
      i = 0;
      while(i < 6 && flag == 0) {
         if(vlb[ind[i][0]] <=0 && vlb[ind[i][1]] <=0 && vlb[ind[i][2]] <=0 
            && vub[ind[i][0]] >=0 && vub[ind[i][1]] >=0 && vub[ind[i][2]] <=0) {
            idx = i;   // store the index of the permutation satisfying the condition
            flag = 5;  // this is case 5
         }
         i++; 
      }
      i = 0;
      while(i < 6 && flag == 0) {
         if(vlb[ind[i][0]] <=0 && vub[ind[i][0]] >=0 && vub[ind[i][1]] <=0 && vub[ind[i][2]] <=0) {
            idx = i;   // store the index of the permutation satisfying the condition
            flag = 6;  // this is case 6
         }
         i++; 
      }
      i = 0;
      while(i < 6 && flag == 0) {
         if(vlb[ind[i][0]] >=0 && vlb[ind[i][1]] >=0 && vlb[ind[i][2]] >=0) {
            idx = i;   // store the index of the permutation satisfying the condition
            flag = 7;  // this is case 7
         }
         i++; 
      }
      i = 0;
      while(i < 6 && flag == 0) {
         if(vlb[ind[i][0]] >=0 && vlb[ind[i][1]] >=0 && vlb[ind[i][2]] <=0 && vub[ind[i][2]] <=0) {
            idx = i;   // store the index of the permutation satisfying the condition
            flag = 8;  // this is case 8
         }
         i++; 
      }
      i = 0;
      while(i < 6 && flag == 0) {
         if(vlb[ind[i][0]] >=0 && vlb[ind[i][1]] <=0 && vub[ind[i][1]] <=0 
            && vlb[ind[i][2]] <=0 && vub[ind[i][2]] <=0) {
            idx = i;   // store the index of the permutation satisfying the condition
            flag = 9;  // this is case 9
         }
         i++; 
      }
      i = 0;
      while(i < 6 && flag == 0) {
         if(vlb[ind[i][0]] <=0 && vub[ind[i][0]] <=0 && vlb[ind[i][1]] <=0 && vub[ind[i][1]] <=0 
            && vlb[ind[i][2]] <=0 && vub[ind[i][2]] <=0) {
            idx = i;   // store the index of the permutation satisfying the condition
            flag = 10;  // this is case 10
         }
         i++; 
      }

      if(flag==0) {
         cout << "RQuarticConvex: error: case not implemented" << endl;
         exit(0);
      }

      v1 = ind[idx][0]; v2 = ind[idx][1]; v3 = ind[idx][2];

      Expression xL1(cf*vlb[v1]); Expression xU1(cf*vub[v1]);
      Expression xL2(vlb[v2]); Expression xU2(vub[v2]);
      Expression xL3(vlb[v3]); Expression xU3(vub[v3]);

      Expression xx1(1, v1+1, InProb->GetVariableLI(v1+1)->Name);
      Expression xx2(1, v2+1, InProb->GetVariableLI(v2+1)->Name);
      Expression xx3(1, v3+1, InProb->GetVariableLI(v3+1)->Name);
      //cout << "v1 =" << v1 << "  v2 =" << v2 << "  v3=" << v3 << endl;
      //cout << "xx1 :" << xx1->ToString() << "  xx2 :" << xx2->ToString() << "  xx3: " << xx3->ToString() << endl;

      if(trilin2 == 1) {
          w=e;
      } 


      // case 1
      if(flag == 1) {
      cout << " -- case 1 --" << endl;

         for(int ii = 0; ii < 12; ii++) {
             cc.push_back(0.);
         }

         Expression theta  = xL1*xU2*xU3 - xU1*xU2*xL3 - xL1*xL2*xU3 + xU1*xL2*xU3; 
         Expression theta1 = xU1*xL2*xL3 - xU1*xU2*xU3 - xL1*xL2*xL3 + xL1*xU2*xL3; 

         cc[0]  = xU2*xU3*xx1 + xU1*xU3*xx2 + xU1*xU2*xx3 - 2.*xU1*xU2*xU3;  cc[0] = w- cc[0];
         cc[1]  = xU2*xL3*xx1 + xL1*xU3*xx2 + xL1*xU2*xx3 - xL1*xU2*xL3 - xL1*xU2*xU3;  cc[1] = w- cc[1];
         cc[2]  = xU2*xL3*xx1 + xL1*xL3*xx2 + xL1*xL2*xx3 - xL1*xU2*xL3 - xL1*xL2*xL3;  cc[2] = w- cc[2];
         cc[3]  = xL2*xU3*xx1 + xU1*xL3*xx2 + xU1*xL2*xx3 - xU1*xL2*xU3 - xU1*xL2*xL3;  cc[3] = w- cc[3];
         cc[4]  = xL2*xL3*xx1 + xU1*xL3*xx2 + xL1*xL2*xx3 - xU1*xL2*xL3 - xL1*xL2*xL3;  cc[4] = w- cc[4];
         cc[5]  = xL2*xU3*xx1 + xL1*xU3*xx2 + (theta/(xU3-xL3))*xx3 + (-(theta*xL3)/(xU3-xL3))
                  - xL1*xU2*xU3 - xU1*xL2*xU3 + xU1*xU2*xL3;  cc[5] = w- cc[5];

         cc[6]  = xU2*xL3*xx1 + xU1*xL3*xx2 + xU1*xU2*xx3 - 2.*xU1*xU2*xL3;  cc[6] = w- cc[6];
         cc[7]  = xL2*xL3*xx1 + xU1*xU3*xx2 + xU1*xL2*xx3 - xU1*xL2*xU3 - xU1*xL2*xL3;  cc[7] = w- cc[7];
         cc[8]  = xU2*xU3*xx1 + xL1*xU3*xx2 + xL1*xL2*xx3 - xL1*xU2*xU3 - xL1*xL2*xU3;  cc[8] = w- cc[8];
         cc[9]  = xU2*xU3*xx1 + xL1*xL3*xx2 + xL1*xU2*xx3 - xL1*xU2*xU3 - xL1*xU2*xL3;  cc[9] = w- cc[9];
         cc[10] = xL2*xU3*xx1 + xU1*xU3*xx2 + xL1*xL2*xx3 - xU1*xL2*xU3 - xL1*xL2*xU3;  cc[10] = w- cc[10];
         cc[11] = xL2*xL3*xx1 + xL1*xL3*xx2 + (theta1/(xL3-xU3))*xx3 + (-(theta1*xU3)/(xL3-xU3))
                  - xU1*xL2*xL3 - xL1*xU2*xL3 + xU1*xU2*xU3;  cc[11] = w- cc[11];
      } // end if case 1

      // case 2
      if(flag == 2) {
      cout << " -- case 2 --" << endl;

         for(int ii = 0; ii < 12; ii++) {
             cc.push_back(0.);
         }

         // compute the 6 permutations of the 3 variables 
         ibnd[0] = v1; ibnd[1] = v2; ibnd[2] = v3; 
         permutation3(ind,ibnd);
         int i, flagg=0, idx=0;
         i = 0;
         while(i < 6 && flagg == 0) {
            if(vub[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]]
                <= vlb[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vub[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]])
            {
               idx = i;   // store the index of the permutation satisfying the condition
               flagg = 1;  // condition is satisfied
            }
            i++; 
         }
         if (flagg==0) {
            cout << "ERROR!!!" << endl; exit(0);
         }
         v1 = ind[idx][0]; v2 = ind[idx][1]; v3 = ind[idx][2];
         //cout << "vlb[v1] =" << vlb[v1] << "  vlb[v2] =" << vlb[v2] << "  vlb[v3]=" << vlb[v3] << endl;
         //cout << "vub[v1] =" << vub[v1] << "  vub[v2] =" << vub[v2] << "  vub[v3]=" << vub[v3] << endl;
         Expression xL1(cf*vlb[v1]); Expression xU1(cf*vub[v1]);
         Expression xL2(vlb[v2]); Expression xU2(vub[v2]);
         Expression xL3(vlb[v3]); Expression xU3(vub[v3]);

         Expression xx1(1, v1+1, InProb->GetVariableLI(v1+1)->Name);
         Expression xx2(1, v2+1, InProb->GetVariableLI(v2+1)->Name);
         Expression xx3(1, v3+1, InProb->GetVariableLI(v3+1)->Name);
         cout << "v1 =" << v1 << "  v2 =" << v2 << "  v3=" << v3 << endl;
         cout << "xx1 :" << xx1->ToString() << "  xx2 :" << xx2->ToString() << "  xx3: " << xx3->ToString() << endl;

       //if(vub[v1]*vub[v2]*vlb[v3] + vlb[v1]*vlb[v2]*vub[v3]
         // <= vlb[v1]*vub[v2]*vlb[v3] + vub[v1]*vlb[v2]*vub[v3]) {

         Expression theta1 = xL1*xL2*xL3 - xU1*xU2*xL3 - xL1*xL2*xU3 + xU1*xL2*xU3;
         Expression theta2 = xU1*xL2*xU3 - xU1*xU2*xL3 - xL1*xL2*xU3 + xL1*xU2*xU3;

         cc[0]  = xU2*xU3*xx1 + xU1*xU3*xx2 + xU1*xU2*xx3 - 2.*xU1*xU2*xU3;  cc[0] = w- cc[0];
         cc[1]  = xL2*xL3*xx1 + xU1*xL3*xx2 + xU1*xL2*xx3 - 2.*xU1*xL2*xL3;  cc[1] = w- cc[1];
         cc[2]  = xU2*xL3*xx1 + xL1*xU3*xx2 + xL1*xU2*xx3 - xL1*xU2*xL3 - xL1*xU2*xU3;  cc[2] = w- cc[2];
         cc[3]  = xU2*xL3*xx1 + xL1*xL3*xx2 + xL1*xL2*xx3 - xL1*xU2*xL3 - xL1*xL2*xL3;  cc[3] = w- cc[3];
         cc[4]  = xL2*xU3*xx1 + (theta1/(xL2-xU2))*xx2 + xL1*xL2*xx3 + (-(theta1*xU2)/(xL2-xU2))
                  - xL1*xL2*xL3 - xU1*xL2*xU3 + xU1*xU2*xL3;  cc[4] = w- cc[4];
         cc[5]  = xL2*xU3*xx1 + xL1*xU3*xx2 + (theta2/(xU3-xL3))*xx3 + (-(theta2*xL3)/(xU3-xL3))
                  - xU1*xL2*xU3 - xL1*xU2*xU3 + xU1*xU2*xL3;  cc[5] = w- cc[5];

      

       if(vlb[v1]*vlb[v2]*vlb[v3] + vub[v1]*vub[v2]*vub[v3]
          >= vub[v1]*vlb[v2]*vlb[v3] + vlb[v1]*vub[v2]*vub[v3]) {

         Expression theta1c = xU1*xL2*xL3 - xU1*xU2*xU3 - xL1*xL2*xL3 + xL1*xU2*xL3;
         Expression theta2c = xU1*xL2*xL3 - xU1*xU2*xU3 - xL1*xL2*xL3 + xL1*xL2*xU3;

         cc[6]  = xL2*xU3*xx1 + xU1*xU3*xx2 + xU1*xL2*xx3 - 2.*xU1*xL2*xU3;  cc[6] = w- cc[6];
         cc[7]  = xU2*xL3*xx1 + xU1*xL3*xx2 + xU1*xU2*xx3 - 2.*xU1*xU2*xL3;  cc[7] = w- cc[7];
         cc[8]  = xU2*xU3*xx1 + xL1*xU3*xx2 + xL1*xL2*xx3 - xL1*xU2*xU3 - xL1*xL2*xU3;  cc[8] = w- cc[8];
         cc[9]  = xU2*xU3*xx1 + xL1*xL3*xx2 + xL1*xU2*xx3 - xL1*xU2*xU3 - xL1*xU2*xL3;  cc[9] = w- cc[9];
         cc[10] = xL2*xL3*xx1 + xL1*xL3*xx2 + (theta1c/(xL3-xU3))*xx3 + (-(theta1c*xU3)/(xL3-xU3))
                  - xU1*xL2*xL3 - xL1*xU2*xL3 + xU1*xU2*xU3;  cc[10] = w- cc[10];
         cc[11] = xL2*xL3*xx1 + (theta2c/(xL2-xU2))*xx2 + xL1*xL2*xx3 + (-(theta2c*xU2)/(xL2-xU2))
                  - xU1*xL2*xL3 - xL1*xL2*xU3 + xU1*xU2*xU3;  cc[11] = w- cc[11];

       } else {
         Expression theta1c = xU1*xU2*xU3 - xU1*xL2*xL3 - xL1*xU2*xU3 + xL1*xL2*xU3;
         Expression theta2c = xU1*xU2*xU3 - xU1*xL2*xL3 - xL1*xU2*xU3 + xL1*xU2*xL3;

         cc[6]  = xL2*xU3*xx1 + xU1*xU3*xx2 + xU1*xL2*xx3 - 2.*xU1*xL2*xU3;  cc[6] = w- cc[6];
         cc[7]  = xU2*xL3*xx1 + xU1*xL3*xx2 + xU1*xU2*xx3 - 2.*xU1*xU2*xL3;  cc[7] = w- cc[7];
         cc[8]  = xL2*xL3*xx1 + xL1*xL3*xx2 + xL1*xU2*xx3 - xL1*xL2*xL3 - xL1*xU2*xL3;  cc[8] = w- cc[8];
         cc[9]  = xL2*xL3*xx1 + xL1*xU3*xx2 + xL1*xL2*xx3 - xL1*xL2*xL3 - xL1*xL2*xU3;  cc[9] = w- cc[9];
         cc[10] = xU2*xU3*xx1 + xL1*xU3*xx2 + (theta1c/(xU3-xL3))*xx3 + (-(theta1c*xL3)/(xU3-xL3))
                  - xU1*xU2*xU3 - xL1*xL2*xU3 + xU1*xL2*xL3;  cc[10] = w- cc[10];
         cc[11] = xU2*xU3*xx1 + (theta2c/(xU2-xL2))*xx2 + xL1*xU2*xx3 + (-(theta2c*xL2)/(xU2-xL2))
                  - xU1*xU2*xU3 - xL1*xU2*xL3 + xU1*xL2*xL3;  cc[11] = w- cc[11];
       }
      
      } // end if case 2

      // case 3
      if(flag == 3) {
      cout << " -- case 3 --" << endl;

         for(int ii = 0; ii < 10; ii++) {
             cc.push_back(0.);
         }
         int last;

       if(vub[v1]*vub[v2]*vlb[v3] + vub[v1]*vlb[v2]*vub[v3] + vlb[v1]*vub[v2]*vub[v3]
          <= vlb[v1]*vlb[v2]*vlb[v3] + 2.*vub[v1]*vub[v2]*vub[v3] &&
          vlb[v1]*vlb[v2]*vlb[v3] + vub[v1]*vub[v2]*vlb[v3] + vub[v1]*vlb[v2]*vub[v3]
          <= vlb[v1]*vub[v2]*vub[v3] + 2.*vub[v1]*vlb[v2]*vlb[v3] &&
          vlb[v1]*vlb[v2]*vlb[v3] + vub[v1]*vub[v2]*vlb[v3] + vlb[v1]*vub[v2]*vub[v3]
          <= vub[v1]*vlb[v2]*vub[v3] + 2.*vlb[v1]*vub[v2]*vlb[v3] &&
          vlb[v1]*vlb[v2]*vlb[v3] + vub[v1]*vlb[v2]*vub[v3] + vlb[v1]*vub[v2]*vub[v3]
          <= vub[v1]*vub[v2]*vlb[v3] + 2.*vlb[v1]*vlb[v2]*vub[v3] ) {
    
         Expression theta3x = 0.5*(xL1*xU2*xU3 + xL1*xL2*xL3 - xU1*xU2*xL3 - xU1*xL2*xU3)/(xL1-xU1);
         Expression theta3y = 0.5*(xU1*xL2*xU3 + xL1*xL2*xL3 - xU1*xU2*xL3 - xL1*xU2*xU3)/(xL2-xU2);
         Expression theta3z = 0.5*(xU1*xU2*xL3 + xL1*xL2*xL3 - xU1*xL2*xU3 - xL1*xU2*xU3)/(xL3-xU3);
         Expression theta3c = xL1*xL2*xL3 - theta3x*xL1 - theta3y*xL2 - theta3z*xL3;

         cc[0]  = xU2*xL3*xx1 + xL1*xL3*xx2 + xL1*xU2*xx3 - 2.*xL1*xU2*xL3;  cc[0] = w- cc[0];
         cc[1]  = xL2*xU3*xx1 + xL1*xU3*xx2 + xL1*xL2*xx3 - 2.*xL1*xL2*xU3;  cc[1] = w- cc[1];
         cc[2]  = xU2*xU3*xx1 + xU1*xU3*xx2 + xU1*xU2*xx3 - 2.*xU1*xU2*xU3;  cc[2] = w- cc[2];
         cc[3]  = xL2*xL3*xx1 + xU1*xL3*xx2 + xU1*xL2*xx3 - 2.*xU1*xL2*xL3;  cc[3] = w- cc[3];
         cc[4]  = theta3x*xx1 + theta3y*xx2 + theta3z*xx3 + theta3c;  cc[4] = w- cc[4];

         last=4;

       } else if (vub[v1]*vub[v2]*vlb[v3] + vub[v1]*vlb[v2]*vub[v3] + vlb[v1]*vub[v2]*vub[v3] 
                  >= vlb[v1]*vlb[v2]*vlb[v3] + 2.*vub[v1]*vub[v2]*vub[v3]) {
         
         Expression theta1 = xU1*xU2*xL3 - xL1*xL2*xL3 - xU1*xU2*xU3 + xU1*xL2*xU3;
         Expression theta2 = xU1*xL2*xU3 - xL1*xL2*xL3 - xU1*xU2*xU3 + xL1*xU2*xU3;
         Expression theta3 = xU1*xU2*xL3 - xL1*xL2*xL3 - xU1*xU2*xU3 + xL1*xU2*xU3;

         cc[0]  = xU2*xL3*xx1 + xL1*xL3*xx2 + xL1*xU2*xx3 - 2.*xL1*xU2*xL3;  cc[0] = w- cc[0];
         cc[1]  = xL2*xU3*xx1 + xL1*xU3*xx2 + xL1*xL2*xx3 - 2.*xL1*xL2*xU3;  cc[1] = w- cc[1];
         cc[2]  = xL2*xL3*xx1 + xU1*xL3*xx2 + xU1*xL2*xx3 - 2.*xU1*xL2*xL3;  cc[2] = w- cc[2];
         cc[3]  = (theta1/(xU1-xL1))*xx1 + xU1*xU3*xx2 + xU1*xU2*xx3 + (-(theta1*xL1)/(xU1-xL1))
                  - xU1*xU2*xL3 - xU1*xL2*xU3 + xL1*xL2*xL3;  cc[3] = w- cc[3];
         cc[4]  = xU2*xU3*xx1 + xU1*xU3*xx2 + (theta2/(xU3-xL3))*xx3 + (-(theta2*xL3)/(xU3-xL3))
                  - xU1*xL2*xU3 - xL1*xU2*xU3 + xL1*xL2*xL3;  cc[4] = w- cc[4];
         cc[5]  = xU2*xU3*xx1 + (theta3/(xU2-xL2))*xx2 + xU1*xU2*xx3 + (-(theta3*xL2)/(xU2-xL2))
                  - xU1*xU2*xL3 - xL1*xU2*xU3 + xL1*xL2*xL3;  cc[5] = w- cc[5];

         last=5;

       } else {
         // compute the 6 permutations of the 3 variables 
         ibnd[0] = v1; ibnd[1] = v2; ibnd[2] = v3; 
         permutation3(ind,ibnd);
         int i, flagg=0, idx=0;
         i = 0;
         while(i < 6 && flagg == 0) {
            if ((vlb[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] + vub[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vub[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]] 
                 >= vlb[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]] + 2.*vub[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]])     && 
              (vub[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]]
                >= vub[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]] + 2.*vlb[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] &&
               vlb[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]] + vub[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]]
                >= vub[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] + 2.*vlb[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]] &&
               vub[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]] + vub[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]]
                >= vlb[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + 2.*vub[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]] &&
               vub[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vub[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]]
                >= vlb[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]] + 2.*vub[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]])  || 
              (vub[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vlb[ind[i][2]]*vub[ind[i][2]]
                <= vub[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]] + 2.*vlb[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]])  ||
              (vlb[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]] + vub[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]]
                <= vub[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] + 2.*vlb[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]]) )
            {
               idx = i;   // store the index of the permutation satisfying the condition
               flagg = 1;  // condition is satisfied
            }
            i++; 
         }
         if (flagg==0) {
            cout << "ERROR!!!" << endl; exit(0);
         }
         v1 = ind[idx][0]; v2 = ind[idx][1]; v3 = ind[idx][2];
         Expression xL1(cf*vlb[v1]); Expression xU1(cf*vub[v1]);
         Expression xL2(vlb[v2]); Expression xU2(vub[v2]);
         Expression xL3(vlb[v3]); Expression xU3(vub[v3]);

         Expression xx1(1, v1+1, InProb->GetVariableLI(v1+1)->Name);
         Expression xx2(1, v2+1, InProb->GetVariableLI(v2+1)->Name);
         Expression xx3(1, v3+1, InProb->GetVariableLI(v3+1)->Name);
         cout << "v1 =" << v1 << "  v2 =" << v2 << "  v3=" << v3 << endl;
         cout << "xx1 :" << xx1->ToString() << "  xx2 :" << xx2->ToString() << "  xx3: " << xx3->ToString() << endl;

       //} else if (vlb[v1]*vlb[v2]*vlb[v3] + vub[v1]*vub[v2]*vlb[v3] + vub[v1]*vlb[v2]*vub[v3] 
         //         >= vlb[v1]*vub[v2]*vub[v3] + 2.*vub[v1]*vlb[v2]*vlb[v3]) {

         Expression theta1 = xL1*xL2*xL3 - xL1*xU2*xU3 - xU1*xL2*xL3 + xU1*xL2*xU3;
         Expression theta2 = xU1*xU2*xL3 - xL1*xU2*xU3 - xU1*xL2*xL3 + xU1*xL2*xU3;
         Expression theta3 = xL1*xL2*xL3 - xL1*xU2*xU3 - xU1*xL2*xL3 + xU1*xU2*xL3;

         cc[0]  = xU2*xL3*xx1 + xL1*xL3*xx2 + xL1*xU2*xx3 - 2.*xL1*xU2*xL3;  cc[0] = w- cc[0];
         cc[1]  = xL2*xU3*xx1 + xL1*xU3*xx2 + xL1*xL2*xx3 - 2.*xL1*xL2*xU3;  cc[1] = w- cc[1];
         cc[2]  = xU2*xU3*xx1 + xU1*xU3*xx2 + xU1*xU2*xx3 - 2.*xU1*xU2*xU3;  cc[2] = w- cc[2];
         cc[3]  = xL2*xL3*xx1 + (theta1/(xL2-xU2))*xx2 + xU1*xL2*xx3 + (-(theta1*xU2)/(xL2-xU2))
                  - xL1*xL2*xL3 - xU1*xL2*xU3 + xL1*xU2*xU3;  cc[3] = w- cc[3];
         cc[4]  = (theta2/(xU1-xL1))*xx1 + xU1*xL3*xx2 + xU1*xL2*xx3 + (-(theta2*xL1)/(xU1-xL1))
                  - xU1*xU2*xL3 - xU1*xL2*xU3 + xL1*xU2*xU3;  cc[4] = w- cc[4];
         cc[5]  = xL2*xL3*xx1 + xU1*xL3*xx2 + (theta3/(xL3-xU3))*xx3 + (-(theta3*xU3)/(xL3-xU3))
                  - xL1*xL2*xL3 - xU1*xU2*xL3 + xL1*xU2*xU3;  cc[5] = w- cc[5];

         last=5;
       }
       if(last > 4) cc.push_back(0.);  

       if(vub[v1]*vlb[v2]*vlb[v3] + vlb[v1]*vub[v2]*vlb[v3] + vlb[v1]*vlb[v2]*vub[v3]
          >= vub[v1]*vub[v2]*vub[v3] + 2.*vlb[v1]*vlb[v2]*vlb[v3] &&
          vlb[v1]*vub[v2]*vlb[v3] + vlb[v1]*vlb[v2]*vub[v3] + vub[v1]*vub[v2]*vub[v3]
          >= vub[v1]*vlb[v2]*vlb[v3] + 2.*vlb[v1]*vub[v2]*vub[v3] &&
          vub[v1]*vlb[v2]*vlb[v3] + vlb[v1]*vlb[v2]*vub[v3] + vub[v1]*vub[v2]*vub[v3]
          >= vlb[v1]*vub[v2]*vlb[v3] + 2.*vub[v1]*vlb[v2]*vub[v3] &&
          vub[v1]*vlb[v2]*vlb[v3] + vlb[v1]*vub[v2]*vlb[v3] + vub[v1]*vub[v2]*vub[v3]
          >= vlb[v1]*vlb[v2]*vub[v3] + 2.*vub[v1]*vub[v2]*vlb[v3] ) {

      cout << "if 1c" << endl;
         Expression theta3x = 0.5*(xU1*xL2*xL3 + xU1*xU2*xU3 - xL1*xL2*xU3 - xL1*xU2*xL3)/(xU1-xL1);
         Expression theta3y = 0.5*(xL1*xU2*xL3 + xU1*xU2*xU3 - xL1*xL2*xU3 - xU1*xL2*xL3)/(xU2-xL2);
         Expression theta3z = 0.5*(xL1*xL2*xU3 + xU1*xU2*xU3 - xL1*xU2*xL3 - xU1*xL2*xL3)/(xU3-xL3);
         Expression theta3c = xU1*xU2*xU3 - theta3x*xU1 - theta3y*xU2 - theta3z*xU3;

         cc[last+1]  = xL2*xL3*xx1 + xL1*xL3*xx2 + xL1*xL2*xx3 - 2.*xL1*xL2*xL3;  cc[last+1] = w- cc[last+1];
         cc[last+2]  = xU2*xL3*xx1 + xU1*xL3*xx2 + xU1*xU2*xx3 - 2.*xU1*xU2*xL3;  cc[last+2] = w- cc[last+2];
         cc[last+3]  = xL2*xU3*xx1 + xU1*xU3*xx2 + xU1*xL2*xx3 - 2.*xU1*xL2*xU3;  cc[last+3] = w- cc[last+3];
         cc[last+4]  = xU2*xU3*xx1 + xL1*xU3*xx2 + xL1*xU2*xx3 - 2.*xL1*xU2*xU3;  cc[last+4] = w- cc[last+4];
         cc[last+5]  = theta3x*xx1 + theta3y*xx2 + theta3z*xx3 + theta3c;  cc[last+5] = w- cc[last+5];

       } else if (vub[v1]*vlb[v2]*vlb[v3] + vlb[v1]*vub[v2]*vlb[v3] + vlb[v1]*vlb[v2]*vub[v3] 
                  <= vub[v1]*vub[v2]*vub[v3] + 2.*vlb[v1]*vlb[v2]*vlb[v3]) {

      cout << "if 2c" << endl;
         Expression theta1 = xU1*xL2*xL3 - xU1*xU2*xU3 - xL1*xL2*xL3 + xL1*xU2*xL3;
         Expression theta2 = xL1*xL2*xU3 - xU1*xU2*xU3 - xL1*xL2*xL3 + xL1*xU2*xL3;
         Expression theta3 = xL1*xL2*xU3 - xU1*xU2*xU3 - xL1*xL2*xL3 + xU1*xL2*xL3;

         cc.push_back(0.);  
         cc[last+1]  = xU2*xL3*xx1 + xU1*xL3*xx2 + xU1*xU2*xx3 - 2.*xU1*xU2*xL3;  cc[last+1] = w- cc[last+1];
         cc[last+2]  = xL2*xU3*xx1 + xU1*xU3*xx2 + xU1*xL2*xx3 - 2.*xU1*xL2*xU3;  cc[last+2] = w- cc[last+2];
         cc[last+3]  = xU2*xU3*xx1 + xL1*xU3*xx2 + xL1*xU2*xx3 - 2.*xL1*xU2*xU3;  cc[last+3] = w- cc[last+3];
         cc[last+4]  = xL2*xL3*xx1 + xL1*xL3*xx2 + (theta1/(xL3-xU3))*xx3 + (-(theta1*xU3)/(xL3-xU3))
                       - xU1*xL2*xL3 - xL1*xU2*xL3 + xU1*xU2*xU3;  cc[last+4] = w- cc[last+4];
         cc[last+5]  = (theta2/(xL1-xU1))*xx1 + xL1*xL3*xx2 + xL1*xL2*xx3 + (-(theta2*xU1)/(xL1-xU1))
                       - xL1*xL2*xU3 - xL1*xU2*xL3 + xU1*xU2*xU3;  cc[last+5] = w- cc[last+5];
         cc[last+6]  = xL2*xL3*xx1 + (theta3/(xL2-xU2))*xx2 + xL1*xL2*xx3 + (-(theta3*xU2)/(xL2-xU2))
                       - xL1*xL2*xU3 - xU1*xL2*xL3 + xU1*xU2*xU3;  cc[last+6] = w- cc[last+6];

       } else if (vlb[v1]*vub[v2]*vlb[v3] + vlb[v1]*vlb[v2]*vub[v3] + vub[v1]*vub[v2]*vub[v3]
                  <= vub[v1]*vlb[v2]*vlb[v3] + 2.*vlb[v1]*vub[v2]*vub[v3]) {

      cout << "if 3c" << endl;
         Expression theta1 = xL1*xL2*xU3 - xU1*xL2*xL3 - xL1*xU2*xU3 + xL1*xU2*xL3;
         Expression theta2 = xU1*xU2*xU3 - xU1*xL2*xL3 - xL1*xU2*xU3 + xL1*xU2*xL3;
         Expression theta3 = xL1*xL2*xU3 - xU1*xL2*xL3 - xL1*xU2*xU3 + xU1*xU2*xU3;

         cc.push_back(0.);  
         cc[last+1]  = xL2*xL3*xx1 + xL1*xL3*xx2 + xL1*xL2*xx3 - 2.*xL1*xL2*xL3;  cc[last+1] = w- cc[last+1];
         cc[last+2]  = xU2*xL3*xx1 + xU1*xL3*xx2 + xU1*xU2*xx3 - 2.*xU1*xU2*xL3;  cc[last+2] = w- cc[last+2];
         cc[last+3]  = xL2*xU3*xx1 + xU1*xU3*xx2 + xU1*xL2*xx3 - 2.*xU1*xL2*xU3;  cc[last+3] = w- cc[last+3];
         cc[last+4]  = (theta1/(xL1-xU1))*xx1 + xL1*xU3*xx2 + xL1*xU2*xx3 + (-(theta1*xU3)/(xL3-xU3))
                       - xL1*xL2*xU3 - xL1*xU2*xL3 + xU1*xL2*xL3;  cc[last+4] = w- cc[last+4];
         cc[last+5]  = xU2*xU3*xx1 + (theta2/(xU2-xL2))*xx2 + xL1*xU2*xx3 + (-(theta2*xL2)/(xU2-xL2))
                       - xU1*xU2*xU3 - xL1*xU2*xL3 + xU1*xL2*xL3;  cc[last+5] = w- cc[last+5];
         cc[last+6]  = xU2*xU3*xx1 + xL1*xU3*xx2 + (theta3/(xU3-xL3))*xx3 + (-(theta3*xL3)/(xU3-xL3))
                       - xL1*xL2*xU3 - xU1*xU2*xU3 + xU1*xL2*xL3;  cc[last+6] = w- cc[last+6];
       }
      } // end if case 3

      // case 4
      if(flag == 4) {
      cout << " -- case 4 --" << endl;

         for(int ii = 0; ii < 12; ii++) {
             cc.push_back(0.);
         }

         Expression theta  = xU1*xL2*xU3 - xU1*xU2*xL3 - xL1*xL2*xU3 + xL1*xL2*xL3;
         Expression theta1 = xL1*xU2*xL3 - xU1*xL2*xL3 - xL1*xU2*xU3 + xU1*xU2*xU3;

         cc[0]  = xL2*xL3*xx1 + xU1*xL3*xx2 + xU1*xL2*xx3 - 2.*xU1*xL2*xL3;  cc[0] = w- cc[0];
         cc[1]  = xU2*xL3*xx1 + xL1*xU3*xx2 + xL1*xU2*xx3 - xL1*xU2*xL3 - xL1*xU2*xU3;  cc[1] = w- cc[1];
         cc[2]  = xU2*xL3*xx1 + xL1*xL3*xx2 + xL1*xL2*xx3 - xL1*xU2*xL3 - xL1*xL2*xL3;  cc[2] = w- cc[2];
         cc[3]  = xL2*xU3*xx1 + xU1*xU3*xx2 + xU1*xU2*xx3 - xU1*xL2*xU3 - xU1*xU2*xU3;  cc[3] = w- cc[3];
         cc[4]  = xU2*xU3*xx1 + xL1*xU3*xx2 + xU1*xU2*xx3 - xL1*xU2*xU3 - xU1*xU2*xU3;  cc[4] = w- cc[4];
         cc[5]  = xL2*xU3*xx1 + (theta/(xL2-xU2))*xx2 + xL1*xL2*xx3 + (-(theta*xU2)/(xL2-xU2))
                  - xU1*xL2*xU3 - xL1*xL2*xL3 + xU1*xU2*xL3;  cc[5] = w- cc[5];

         cc[6]  = xU2*xL3*xx1 + xU1*xL3*xx2 + xU1*xU2*xx3 - 2.*xU1*xU2*xL3;  cc[6] = w- cc[6];
         cc[7]  = xU2*xU3*xx1 + xU1*xU3*xx2 + xU1*xL2*xx3 - xU1*xU2*xU3 - xU1*xL2*xU3;  cc[7] = w- cc[7];
         cc[8]  = xL2*xL3*xx1 + xL1*xL3*xx2 + xL1*xU2*xx3 - xL1*xU2*xL3 - xL1*xL2*xL3;  cc[8] = w- cc[8];
         cc[9]  = xL2*xL3*xx1 + xL1*xU3*xx2 + xL1*xL2*xx3 - xL1*xL2*xU3 - xL1*xL2*xL3;  cc[9] = w- cc[9];
         cc[10] = xL2*xU3*xx1 + xL1*xU3*xx2 + xU1*xL2*xx3 - xL1*xL2*xU3 - xU1*xL2*xU3;  cc[10] = w- cc[10];
         cc[11] = xU2*xU3*xx1 + (theta1/(xU2-xL2))*xx2 + xL1*xU2*xx3 + (-(theta1*xL2)/(xU2-xL2))
                  - xL1*xU2*xL3 - xU1*xU2*xU3 + xU1*xL2*xL3;  cc[11] = w- cc[11];

      } // end if case 4

      // case 5
      if(flag == 5) {
      cout << " -- case 5 --" << endl;

         for(int ii = 0; ii < 12; ii++) {
             cc.push_back(0.);
         }

         // compute the 6 permutations of the 3 variables 
         ibnd[0] = v1; ibnd[1] = v2; ibnd[2] = v3; 
         permutation3(ind,ibnd);
         int i, flagg=0, idx=0;
         i = 0;
         while(i < 6 && flagg == 0) {
            if(vub[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]]
                >= vlb[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vub[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]]) 
            {
               idx = i;   // store the index of the permutation satisfying the condition
               flagg = 1;  // condition is satisfied
            }
            i++; 
         }
         if (flagg==0) {
            cout << "ERROR!!!" << endl; exit(0);
         }
         v1 = ind[idx][0]; v2 = ind[idx][1]; v3 = ind[idx][2];
         Expression xL1(cf*vlb[v1]); Expression xU1(cf*vub[v1]);
         Expression xL2(vlb[v2]); Expression xU2(vub[v2]);
         Expression xL3(vlb[v3]); Expression xU3(vub[v3]);

         Expression xx1(1, v1+1, InProb->GetVariableLI(v1+1)->Name);
         Expression xx2(1, v2+1, InProb->GetVariableLI(v2+1)->Name);
         Expression xx3(1, v3+1, InProb->GetVariableLI(v3+1)->Name);
         cout << "v1 =" << v1 << "  v2 =" << v2 << "  v3=" << v3 << endl;
         cout << "xx1 :" << xx1->ToString() << "  xx2 :" << xx2->ToString() << "  xx3: " << xx3->ToString() << endl;


       if(vlb[v1]*vlb[v2]*vlb[v3] + vub[v1]*vub[v2]*vub[v3]
          <= vub[v1]*vub[v2]*vlb[v3] + vlb[v1]*vlb[v2]*vub[v3]) {

         Expression theta1 = xU1*xU2*xL3 - xL1*xL2*xL3 - xU1*xU2*xU3 + xU1*xL2*xU3;
         Expression theta2 = xU1*xU2*xL3 - xL1*xL2*xL3 - xU1*xU2*xU3 + xL1*xU2*xU3;

         cc[0]  = xU2*xL3*xx1 + xL1*xL3*xx2 + xL1*xU2*xx3 - 2.*xL1*xU2*xL3;  cc[0] = w- cc[0];
         cc[1]  = xL2*xL3*xx1 + xU1*xL3*xx2 + xU1*xL2*xx3 - 2.*xU1*xL2*xL3;  cc[1] = w- cc[1];
         cc[2]  = xU2*xU3*xx1 + xL1*xU3*xx2 + xL1*xL2*xx3 - xL1*xU2*xU3 - xL1*xL2*xU3;  cc[2] = w- cc[2];
         cc[3]  = xL2*xU3*xx1 + xU1*xU3*xx2 + xL1*xL2*xx3 - xU1*xL2*xU3 - xL1*xL2*xU3;  cc[3] = w- cc[3];
         cc[4]  = (theta1/(xU1-xL1))*xx1 + xU1*xU3*xx2 + xU1*xU2*xx3 + (-(theta1*xL1)/(xU1-xL1))
                  - xU1*xU2*xL3 - xU1*xL2*xU3 + xL1*xL2*xL3;  cc[4] = w- cc[4];
         cc[5]  = xU2*xU3*xx1 + (theta2/(xU2-xL2))*xx2 + xU1*xU2*xx3 + (-(theta2*xL2)/(xU2-xL2))
                  - xU1*xU2*xL3 - xL1*xU2*xU3 + xL1*xL2*xL3;  cc[5] = w- cc[5];

       } else {
         Expression theta1 = xL1*xL2*xL3 - xU1*xU2*xL3 - xL1*xL2*xU3 + xL1*xU2*xU3;
         Expression theta2 = xL1*xL2*xL3 - xU1*xU2*xL3 - xL1*xL2*xU3 + xU1*xL2*xU3;

         cc[0]  = xU2*xL3*xx1 + xL1*xL3*xx2 + xL1*xU2*xx3 - 2.*xL1*xU2*xL3;  cc[0] = w- cc[0];
         cc[1]  = xL2*xL3*xx1 + xU1*xL3*xx2 + xU1*xL2*xx3 - 2.*xU1*xL2*xL3;  cc[1] = w- cc[1];
         cc[2]  = xL2*xU3*xx1 + xU1*xU3*xx2 + xU1*xU2*xx3 - xU1*xL2*xU3 - xU1*xU2*xU3;  cc[2] = w- cc[2];
         cc[3]  = xU2*xU3*xx1 + xL1*xU3*xx2 + xU1*xU2*xx3 - xL1*xU2*xU3 - xU1*xU2*xU3;  cc[3] = w- cc[3];
         cc[4]  = (theta1/(xL1-xU1))*xx1 + xL1*xU3*xx2 + xL1*xL2*xx3 + (-(theta1*xU1)/(xL1-xU1))
                  - xL1*xL2*xL3 - xL1*xU2*xU3 + xU1*xU2*xL3;  cc[4] = w- cc[4];
         cc[5]  = xL2*xU3*xx1 + (theta2/(xL2-xU2))*xx2 + xL1*xL2*xx3 + (-(theta2*xU2)/(xL2-xU2))
                  - xL1*xL2*xL3 - xU1*xL2*xU3 + xU1*xU2*xL3;  cc[5] = w- cc[5];
       }

       //if(vub[v1]*vlb[v2]*vlb[v3] + vlb[v1]*vub[v2]*vub[v3]
         // >= vlb[v1]*vub[v2]*vlb[v3] + vub[v1]*vlb[v2]*vub[v3]) {

         Expression theta1c = xL1*xU2*xL3 - xU1*xL2*xL3 - xL1*xU2*xU3 + xL1*xL2*xU3;
         Expression theta2c = xU1*xU2*xU3 - xU1*xL2*xL3 - xL1*xU2*xU3 + xL1*xU2*xL3;

         cc[6]  = xL2*xL3*xx1 + xL1*xL3*xx2 + xL1*xL2*xx3 - 2.*xL1*xL2*xL3;  cc[6] = w- cc[6];
         cc[7]  = xU2*xL3*xx1 + xU1*xL3*xx2 + xU1*xU2*xx3 - 2.*xU1*xU2*xL3;  cc[7] = w- cc[7];
         cc[8]  = xL2*xU3*xx1 + xL1*xU3*xx2 + xU1*xL2*xx3 - xL1*xL2*xU3 - xU1*xL2*xU3;  cc[8] = w- cc[8];
         cc[9]  = xU2*xU3*xx1 + xU1*xU3*xx2 + xU1*xL2*xx3 - xU1*xU2*xU3 - xU1*xL2*xU3;  cc[9] = w- cc[9];
         cc[10] = (theta1c/(xL1-xU1))*xx1 + xL1*xU3*xx2 + xL1*xU2*xx3 + (-(theta1c*xU1)/(xL1-xU1))
                  - xL1*xU2*xL3 - xL1*xL2*xU3 + xU1*xL2*xL3;  cc[10] = w- cc[10];
         cc[11] = xU2*xU3*xx1 + (theta2c/(xU2-xL2))*xx2 + xL1*xU2*xx3 + (-(theta2c*xL2)/(xU2-xL2))
                  - xU1*xU2*xU3 - xL1*xU2*xL3 + xU1*xL2*xL3;  cc[11] = w- cc[11];
       //}

      } // end if case 5

      // case 6
      if(flag == 6) {
      cout << " -- case 6 --" << endl;

         for(int ii = 0; ii < 12; ii++) {
             cc.push_back(0.);
         }

         Expression theta = xU1*xU2*xL3 - xL1*xL2*xL3 - xU1*xU2*xU3 + xU1*xL2*xU3;
         Expression theta1 = xL1*xU2*xL3 - xU1*xL2*xL3 - xL1*xU2*xU3 + xL1*xL2*xU3;

         cc[0]  = xL2*xL3*xx1 + xU1*xL3*xx2 + xU1*xL2*xx3 - 2.*xU1*xL2*xL3;  cc[0] = w- cc[0];
         cc[1]  = xU2*xU3*xx1 + xL1*xL3*xx2 + xL1*xU2*xx3 - xL1*xU2*xU3 - xL1*xU2*xL3;  cc[1] = w- cc[1];
         cc[2]  = xU2*xU3*xx1 + xL1*xU3*xx2 + xL1*xL2*xx3 - xL1*xU2*xU3 - xL1*xL2*xU3;  cc[2] = w- cc[2];
         cc[3]  = xL2*xU3*xx1 + xU1*xU3*xx2 + xL1*xL2*xx3 - xU1*xL2*xU3 - xL1*xL2*xU3;  cc[3] = w- cc[3];
         cc[4]  = xU2*xL3*xx1 + xL1*xL3*xx2 + xU1*xU2*xx3 - xL1*xU2*xL3 - xU1*xU2*xL3;  cc[4] = w- cc[4];
         cc[5]  = (theta/(xU1-xL1))*xx1 + xU1*xU3*xx2 + xU1*xU2*xx3 + (-(theta*xL1)/(xU1-xL1))
                  - xU1*xU2*xL3 - xU1*xL2*xU3 + xL1*xL2*xL3;  cc[5] = w- cc[5];

         cc[6]  = xL2*xL3*xx1 + xL1*xL3*xx2 + xL1*xL2*xx3 - 2.*xL1*xL2*xL3;  cc[6] = w- cc[6];
         cc[7]  = xL2*xU3*xx1 + xL1*xU3*xx2 + xU1*xL2*xx3 - xL1*xL2*xU3 - xU1*xL2*xU3;  cc[7] = w- cc[7];
         cc[8]  = xU2*xU3*xx1 + xU1*xU3*xx2 + xU1*xL2*xx3 - xU1*xU2*xU3 - xU1*xL2*xU3;  cc[8] = w- cc[8];
         cc[9]  = xU2*xU3*xx1 + xU1*xL3*xx2 + xU1*xU2*xx3 - xU1*xU2*xU3 - xU1*xU2*xL3;  cc[9] = w- cc[9];
         cc[10] = xU2*xL3*xx1 + xU1*xL3*xx2 + xL1*xU2*xx3 - xL1*xU2*xL3 - xU1*xU2*xL3;  cc[10] = w- cc[10];
         cc[11] = (theta1/(xL1-xU1))*xx1 + xL1*xU3*xx2 + xL1*xU2*xx3 + (-(theta1*xU1)/(xL1-xU1))
                  - xL1*xU2*xL3 - xL1*xL2*xU3 + xU1*xL2*xL3;  cc[11] = w- cc[11];
      } // end if case 6


      // case 7
      if(flag == 7) {
      cout << " -- case 7 --" << endl;

         for(int ii = 0; ii < 12; ii++) {
             cc.push_back(0.);
         }

         // compute the 6 permutations of the 3 variables 
         ibnd[0] = v1; ibnd[1] = v2; ibnd[2] = v3; 
         permutation3(ind,ibnd);
         int i, flagg=0, idx=0;
         i = 0;
         while(i < 6 && flagg == 0) {
            if(vub[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]]
                <= vlb[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vub[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]] &&
               vub[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]]
                <= vub[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]]) 
            {
               idx = i;   // store the index of the permutation satisfying the condition
               flagg = 1;  // condition is satisfied
            }
            i++; 
         }
         if (flagg==0) {
            cout << "ERROR!!!" << endl; exit(0);
         }
         v1 = ind[idx][0]; v2 = ind[idx][1]; v3 = ind[idx][2];
         Expression xL1(cf*vlb[v1]); Expression xU1(cf*vub[v1]);
         Expression xL2(vlb[v2]); Expression xU2(vub[v2]);
         Expression xL3(vlb[v3]); Expression xU3(vub[v3]);

         Expression xx1(1, v1+1, InProb->GetVariableLI(v1+1)->Name);
         Expression xx2(1, v2+1, InProb->GetVariableLI(v2+1)->Name);
         Expression xx3(1, v3+1, InProb->GetVariableLI(v3+1)->Name);
         cout << "v1 =" << v1 << "  v2 =" << v2 << "  v3=" << v3 << endl;
         cout << "xx1 :" << xx1->ToString() << "  xx2 :" << xx2->ToString() << "  xx3: " << xx3->ToString() << endl;

       //if(vub[v1]*vlb[v2]*vlb[v3] + vlb[v1]*vub[v2]*vub[v3]
         // <= vlb[v1]*vub[v2]*vlb[v3] + vub[v1]*vlb[v2]*vub[v3] &&
         // vub[v1]*vlb[v2]*vlb[v3] + vlb[v1]*vub[v2]*vub[v3]
         // <= vub[v1]*vub[v2]*vlb[v3] + vlb[v1]*vlb[v2]*vub[v3]) {

         Expression theta1 = xU1*xU2*xL3 - xL1*xU2*xU3 - xU1*xL2*xL3 + xU1*xL2*xU3;
         Expression theta2 = xL1*xL2*xU3 - xU1*xL2*xL3 - xL1*xU2*xU3 + xL1*xU2*xL3;

         cc[0]  = xL2*xL3*xx1 + xL1*xL3*xx2 + xL1*xL2*xx3 - 2.*xL1*xL2*xL3;  cc[0] = w- cc[0];
         cc[1]  = xU2*xU3*xx1 + xU1*xU3*xx2 + xU1*xU2*xx3 - 2.*xU1*xU2*xU3;  cc[1] = w- cc[1];
         cc[2]  = xL2*xU3*xx1 + xL1*xU3*xx2 + xU1*xL2*xx3 - xL1*xL2*xU3 - xU1*xL2*xU3;  cc[2] = w- cc[2];
         cc[3]  = xU2*xL3*xx1 + xU1*xL3*xx2 + xL1*xU2*xx3 - xU1*xU2*xL3 - xL1*xU2*xL3;  cc[3] = w- cc[3];
         cc[4]  = (theta1/(xU1-xL1))*xx1 + xU1*xL3*xx2 + xU1*xL2*xx3 + (-(theta1*xL1)/(xU1-xL1))
                  - xU1*xU2*xL3 - xU1*xL2*xU3 + xL1*xU2*xU3;  cc[4] = w- cc[4];
         cc[5]  = (theta2/(xL1-xU1))*xx1 + xL1*xU3*xx2 + xL1*xU2*xx3 + (-(theta2*xU1)/(xL1-xU1))
                  - xL1*xL2*xU3 - xL1*xU2*xL3 + xU1*xL2*xL3;  cc[5] = w- cc[5];
       //}
         cc[6]  = xL2*xL3*xx1 + xU1*xL3*xx2 + xU1*xU2*xx3 - xU1*xU2*xL3 - xU1*xL2*xL3;  cc[6] = w- cc[6];
         cc[7]  = xU2*xL3*xx1 + xL1*xL3*xx2 + xU1*xU2*xx3 - xU1*xU2*xL3 - xL1*xU2*xL3;  cc[7] = w- cc[7];
         cc[8]  = xL2*xL3*xx1 + xU1*xU3*xx2 + xU1*xL2*xx3 - xU1*xL2*xU3 - xU1*xL2*xL3;  cc[8] = w- cc[8];
         cc[9]  = xU2*xU3*xx1 + xL1*xL3*xx2 + xL1*xU2*xx3 - xL1*xU2*xU3 - xL1*xU2*xL3;  cc[9] = w- cc[9];
         cc[10] = xL2*xU3*xx1 + xU1*xU3*xx2 + xL1*xL2*xx3 - xU1*xL2*xU3 - xL1*xL2*xU3;  cc[10] = w- cc[10];
         cc[11] = xU2*xU3*xx1 + xL1*xU3*xx2 + xL1*xL2*xx3 - xL1*xU2*xU3 - xL1*xL2*xU3;  cc[11] = w- cc[11];
      } // end if case 7
      
      // case 8
      if(flag == 8) {
      cout << " -- case 8 --" << endl;

         for(int ii = 0; ii < 12; ii++) {
             cc.push_back(0.);
         }

         // compute the 6 permutations of the 3 variables 
         ibnd[0] = v1; ibnd[1] = v2; ibnd[2] = v3; 
         permutation3(ind,ibnd);
         int i, flagg=0, idx=0;
         i = 0;
         while(i < 6 && flagg == 0) {
            if(vub[ind[i][2]] <=0  &&
               ((vlb[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] + vub[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]]
                 >= vub[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]] &&
                vlb[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] + vub[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]]
                 >= vlb[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vub[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]]) ||
                (vub[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]]
                  >= vlb[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vub[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]])) ) 
            {
               idx = i;   // store the index of the permutation satisfying the condition
               flagg = 1;  // condition is satisfied
            }
            i++; 
         }
         if (flagg==0) {
            cout << "ERROR!!!" << endl; exit(0);
         }
         v1 = ind[idx][0]; v2 = ind[idx][1]; v3 = ind[idx][2];
         Expression xL1(cf*vlb[v1]); Expression xU1(cf*vub[v1]);
         Expression xL2(vlb[v2]); Expression xU2(vub[v2]);
         Expression xL3(vlb[v3]); Expression xU3(vub[v3]);

         Expression xx1(1, v1+1, InProb->GetVariableLI(v1+1)->Name);
         Expression xx2(1, v2+1, InProb->GetVariableLI(v2+1)->Name);
         Expression xx3(1, v3+1, InProb->GetVariableLI(v3+1)->Name);
         cout << "v1 =" << v1 << "  v2 =" << v2 << "  v3=" << v3 << endl;
         cout << "xx1 :" << xx1->ToString() << "  xx2 :" << xx2->ToString() << "  xx3: " << xx3->ToString() << endl;

       //if(vub[v3]<=0) {

         cc[0]  = xU2*xL3*xx1 + xL1*xL3*xx2 + xL1*xL2*xx3 - xL1*xU2*xL3 - xL1*xL2*xL3;  cc[0] = w- cc[0];
         cc[1]  = xU2*xL3*xx1 + xL1*xU3*xx2 + xL1*xU2*xx3 - xL1*xU2*xL3 - xL1*xU2*xU3;  cc[1] = w- cc[1];
         cc[2]  = xL2*xU3*xx1 + xU1*xL3*xx2 + xU1*xL2*xx3 - xU1*xL2*xU3 - xU1*xL2*xL3;  cc[2] = w- cc[2];
         cc[3]  = xL2*xU3*xx1 + xU1*xU3*xx2 + xU1*xU2*xx3 - xU1*xL2*xU3 - xU1*xU2*xU3;  cc[3] = w- cc[3];
         cc[4]  = xL2*xL3*xx1 + xU1*xL3*xx2 + xL1*xL2*xx3 - xU1*xL2*xL3 - xL1*xL2*xL3;  cc[4] = w- cc[4];
         cc[5]  = xU2*xU3*xx1 + xL1*xU3*xx2 + xU1*xU2*xx3 - xU1*xU2*xU3 - xL1*xU2*xU3;  cc[5] = w- cc[5];

         if(vlb[v1]*vlb[v2]*vlb[v3] + vub[v1]*vub[v2]*vub[v3]
            >= vub[v1]*vlb[v2]*vlb[v3] + vlb[v1]*vub[v2]*vub[v3] &&
            vlb[v1]*vlb[v2]*vlb[v3] + vub[v1]*vub[v2]*vub[v3]
            >= vlb[v1]*vub[v2]*vlb[v3] + vub[v1]*vlb[v2]*vub[v3]) {

           Expression theta1c = xU1*xL2*xL3 - xU1*xU2*xU3 - xL1*xL2*xL3 + xL1*xU2*xL3;
           Expression theta2c = xU1*xL2*xU3 - xL1*xL2*xL3 - xU1*xU2*xU3 + xL1*xU2*xU3;

         cc[6]  = xL2*xU3*xx1 + xL1*xU3*xx2 + xL1*xL2*xx3 - 2.*xL1*xL2*xU3;  cc[6] = w- cc[6];
         cc[7]  = xU2*xL3*xx1 + xU1*xL3*xx2 + xU1*xU2*xx3 - 2.*xU1*xU2*xL3;  cc[7] = w- cc[7];
         cc[8]  = xL2*xL3*xx1 + xU1*xU3*xx2 + xU1*xL2*xx3 - xU1*xL2*xU3 - xU1*xL2*xL3;  cc[8] = w- cc[8];
         cc[9]  = xU2*xU3*xx1 + xL1*xL3*xx2 + xL1*xU2*xx3 - xL1*xU2*xU3 - xL1*xU2*xL3;  cc[9] = w- cc[9];
         cc[10] = xL2*xL3*xx1 + xL1*xL3*xx2 + (theta1c/(xL3-xU3))*xx3 + (-(theta1c*xU3)/(xL3-xU3))
                  - xU1*xL2*xL3 - xL1*xU2*xL3 + xU1*xU2*xU3;  cc[10] = w- cc[10];
         cc[11] = xU2*xU3*xx1 + xU1*xU3*xx2 + (theta2c/(xU3-xL3))*xx3 + (-(theta2c*xL3)/(xU3-xL3))
                  - xU1*xL2*xU3 - xL1*xU2*xU3 + xL1*xL2*xL3;  cc[11] = w- cc[11];

         } else if(vub[v1]*vlb[v2]*vlb[v3] + vlb[v1]*vub[v2]*vub[v3]
                    >= vlb[v1]*vub[v2]*vlb[v3] + vub[v1]*vlb[v2]*vub[v3]) {

           Expression theta1c = xL1*xL2*xL3 - xL1*xU2*xU3 - xU1*xL2*xL3 + xU1*xL2*xU3;
           Expression theta2c = xL1*xU2*xL3 - xU1*xL2*xL3 - xL1*xU2*xU3 + xU1*xU2*xU3;

         cc[6]  = xL2*xU3*xx1 + xL1*xU3*xx2 + xL1*xL2*xx3 - 2.*xL1*xL2*xU3;  cc[6] = w- cc[6];
         cc[7]  = xU2*xL3*xx1 + xU1*xL3*xx2 + xU1*xU2*xx3 - 2.*xU1*xU2*xL3;  cc[7] = w- cc[7];
         cc[8]  = xL2*xL3*xx1 + xL1*xL3*xx2 + xL1*xU2*xx3 - xL1*xL2*xL3 - xL1*xU2*xL3;  cc[8] = w- cc[8];
         cc[9]  = xU2*xU3*xx1 + xU1*xU3*xx2 + xU1*xL2*xx3 - xU1*xL2*xU3 - xU1*xU2*xU3;  cc[9] = w- cc[9];
         cc[10] = xL2*xL3*xx1 + (theta1c/(xL2-xU2))*xx2 + xU1*xL2*xx3 + (-(theta1c*xU2)/(xL2-xU2))
                  - xL1*xL2*xL3 - xU1*xL2*xU3 + xL1*xU2*xU3;  cc[10] = w- cc[10];
         cc[11] = xU2*xU3*xx1 + (theta2c/(xU2-xL2))*xx2 + xL1*xU2*xx3 + (-(theta2c*xL2)/(xU2-xL2))
                  - xL1*xU2*xL3 - xU1*xU2*xU3 + xU1*xL2*xL3;  cc[11] = w- cc[11];
         }
      //} 
      } // end if case 8

      // case 9
      if(flag == 9) {
      cout << " -- case 9 --" << endl;

         for(int ii = 0; ii < 12; ii++) {
             cc.push_back(0.);
         }

         // compute the 6 permutations of the 3 variables 
         ibnd[0] = v1; ibnd[1] = v2; ibnd[2] = v3; 
         permutation3(ind,ibnd);
         int i, flagg=0, idx=0;
         i = 0;
         while(i < 6 && flagg == 0) {
            if(vlb[ind[i][0]] >=0  &&
               ((vlb[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] + vub[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]]
                 <= vub[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]] &&
                 vlb[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] + vub[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]]
                 <= vlb[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vub[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]]) || 
                (vub[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]]
                 <= vlb[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] + vub[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]] &&
                 vub[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]]
                 <= vlb[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vub[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]]) ))
            {
               idx = i;   // store the index of the permutation satisfying the condition
               flagg = 1;  // condition is satisfied
            }
            i++; 
         }
         if (flagg==0) {
            cout << "ERROR!!!" << endl; exit(0);
         }
         v1 = ind[idx][0]; v2 = ind[idx][1]; v3 = ind[idx][2];
         Expression xL1(cf*vlb[v1]); Expression xU1(cf*vub[v1]);
         Expression xL2(vlb[v2]); Expression xU2(vub[v2]);
         Expression xL3(vlb[v3]); Expression xU3(vub[v3]);

         Expression xx1(1, v1+1, InProb->GetVariableLI(v1+1)->Name);
         Expression xx2(1, v2+1, InProb->GetVariableLI(v2+1)->Name);
         Expression xx3(1, v3+1, InProb->GetVariableLI(v3+1)->Name);
         cout << "v1 =" << v1 << "  v2 =" << v2 << "  v3=" << v3 << endl;
         cout << "xx1 :" << xx1->ToString() << "  xx2 :" << xx2->ToString() << "  xx3: " << xx3->ToString() << endl;


       //if(vlb[v1]>=0) {
       if(vlb[v1]*vlb[v2]*vlb[v3] + vub[v1]*vub[v2]*vub[v3]
          <= vub[v1]*vub[v2]*vlb[v3] + vlb[v1]*vlb[v2]*vub[v3] &&
          vlb[v1]*vlb[v2]*vlb[v3] + vub[v1]*vub[v2]*vub[v3]
          <= vlb[v1]*vub[v2]*vlb[v3] + vub[v1]*vlb[v2]*vub[v3]) {

         Expression theta1 = xL1*xL2*xU3 - xU1*xU2*xU3 - xL1*xL2*xL3 + xL1*xU2*xL3;
         Expression theta2 = xU1*xL2*xU3 - xL1*xL2*xL3 - xU1*xU2*xU3 + xU1*xU2*xL3;

         cc[0]  = xU2*xU3*xx1 + xL1*xU3*xx2 + xL1*xU2*xx3 - 2.*xL1*xU2*xU3;  cc[0] = w- cc[0];
         cc[1]  = xL2*xL3*xx1 + xU1*xL3*xx2 + xU1*xL2*xx3 - 2.*xU1*xL2*xL3;  cc[1] = w- cc[1];
         cc[2]  = xL2*xU3*xx1 + xU1*xU3*xx2 + xL1*xL2*xx3 - xU1*xL2*xU3 - xL1*xL2*xU3;  cc[2] = w- cc[2];
         cc[3]  = xU2*xL3*xx1 + xL1*xL3*xx2 + xU1*xU2*xx3 - xU1*xU2*xL3 - xL1*xU2*xL3;  cc[3] = w- cc[3];
         cc[4]  = (theta1/(xL1-xU1))*xx1 + xL1*xL3*xx2 + xL1*xL2*xx3 + (-(theta1*xU1)/(xL1-xU1))
                  - xL1*xL2*xU3 - xL1*xU2*xL3 + xU1*xU2*xU3;  cc[4] = w- cc[4];
         cc[5]  = (theta2/(xU1-xL1))*xx1 + xU1*xU3*xx2 + xU1*xU2*xx3 + (-(theta2*xL1)/(xU1-xL1))
                  - xU1*xL2*xU3 - xU1*xU2*xL3 + xL1*xL2*xL3;  cc[5] = w- cc[5];

       } else if(vub[v1]*vub[v2]*vlb[v3] + vlb[v1]*vlb[v2]*vub[v3]
                 <= vlb[v1]*vlb[v2]*vlb[v3] + vub[v1]*vub[v2]*vub[v3] &&
                 vub[v1]*vub[v2]*vlb[v3] + vlb[v1]*vlb[v2]*vub[v3]
                 <= vlb[v1]*vub[v2]*vlb[v3] + vub[v1]*vlb[v2]*vub[v3]) {

         Expression theta1 = xU1*xU2*xU3 - xL1*xL2*xU3 - xU1*xU2*xL3 + xL1*xU2*xL3;
         Expression theta2 = xL1*xL2*xL3 - xU1*xU2*xL3 - xL1*xL2*xU3 + xU1*xL2*xU3;

         cc[0]  = xU2*xU3*xx1 + xL1*xU3*xx2 + xL1*xU2*xx3 - 2.*xL1*xU2*xU3;  cc[0] = w- cc[0];
         cc[1]  = xL2*xL3*xx1 + xU1*xL3*xx2 + xU1*xL2*xx3 - 2.*xU1*xL2*xL3;  cc[1] = w- cc[1];
         cc[2]  = xL2*xU3*xx1 + xU1*xU3*xx2 + xU1*xU2*xx3 - xU1*xL2*xU3 - xU1*xU2*xU3;  cc[2] = w- cc[2];
         cc[3]  = xU2*xL3*xx1 + xL1*xL3*xx2 + xL1*xL2*xx3 - xL1*xL2*xL3 - xL1*xU2*xL3;  cc[3] = w- cc[3];
         cc[4]  = xU2*xL3*xx1 + (theta1/(xU2-xL2))*xx2 + xU1*xU2*xx3 + (-(theta1*xL2)/(xU2-xL2))
                  - xU1*xU2*xU3 - xL1*xU2*xL3 + xL1*xL2*xU3;  cc[4] = w- cc[4];
         cc[5]  = xL2*xU3*xx1 + (theta2/(xL2-xU2))*xx2 + xL1*xL2*xx3 + (-(theta2*xU2)/(xL2-xU2))
                  - xL1*xL2*xL3 - xU1*xL2*xU3 + xU1*xU2*xL3;  cc[5] = w- cc[5];
       }

         cc[6]  = xL2*xL3*xx1 + xL1*xU3*xx2 + xL1*xL2*xx3 - xL1*xL2*xU3 - xL1*xL2*xL3;  cc[6] = w- cc[6];
         cc[7]  = xL2*xL3*xx1 + xL1*xL3*xx2 + xL1*xU2*xx3 - xL1*xU2*xL3 - xL1*xL2*xL3;  cc[7] = w- cc[7];
         cc[8]  = xL2*xU3*xx1 + xL1*xU3*xx2 + xU1*xL2*xx3 - xL1*xL2*xU3 - xU1*xL2*xU3;  cc[8] = w- cc[8];
         cc[9]  = xU2*xU3*xx1 + xU1*xU3*xx2 + xU1*xL2*xx3 - xU1*xU2*xU3 - xU1*xL2*xU3;  cc[9] = w- cc[9];
         cc[10] = xU2*xU3*xx1 + xU1*xL3*xx2 + xU1*xU2*xx3 - xU1*xU2*xU3 - xU1*xU2*xL3;  cc[10] = w- cc[10];
         cc[11] = xU2*xL3*xx1 + xU1*xL3*xx2 + xL1*xU2*xx3 - xL1*xU2*xL3 - xU1*xU2*xL3;  cc[11] = w- cc[11];
      //}
      } // end if case 9

      // case 10
      if(flag == 10) {
      cout << " -- case 10 --" << endl;

         for(int ii = 0; ii < 12; ii++) {
             cc.push_back(0.);
         }

         // compute the 6 permutations of the 3 variables 
         ibnd[0] = v1; ibnd[1] = v2; ibnd[2] = v3; 
         permutation3(ind,ibnd);
         int i, flagg=0, idx=0;
         i = 0;
         while(i < 6 && flagg == 0) {
            if(vub[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]]
                >= vlb[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vub[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]] &&
               vub[ind[i][0]]*vlb[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vub[ind[i][1]]*vub[ind[i][2]]
                >= vub[ind[i][0]]*vub[ind[i][1]]*vlb[ind[i][2]] + vlb[ind[i][0]]*vlb[ind[i][1]]*vub[ind[i][2]]) {
            {
               idx = i;   // store the index of the permutation satisfying the condition
               flagg = 1;  // condition is satisfied
            }
            i++; 
         }
         if (flagg==0) {
            cout << "ERROR!!!" << endl; exit(0);
         }
         v1 = ind[idx][0]; v2 = ind[idx][1]; v3 = ind[idx][2];
         Expression xL1(cf*vlb[v1]); Expression xU1(cf*vub[v1]);
         Expression xL2(vlb[v2]); Expression xU2(vub[v2]);
         Expression xL3(vlb[v3]); Expression xU3(vub[v3]);

         Expression xx1(1, v1+1, InProb->GetVariableLI(v1+1)->Name);
         Expression xx2(1, v2+1, InProb->GetVariableLI(v2+1)->Name);
         Expression xx3(1, v3+1, InProb->GetVariableLI(v3+1)->Name);
         cout << "v1 =" << v1 << "  v2 =" << v2 << "  v3=" << v3 << endl;
         cout << "xx1 :" << xx1->ToString() << "  xx2 :" << xx2->ToString() << "  xx3: " << xx3->ToString() << endl;

         cc[0]  = xL2*xL3*xx1 + xU1*xL3*xx2 + xU1*xU2*xx3 - xU1*xU2*xL3 - xU1*xL2*xL3;  cc[0] = w- cc[0];
         cc[1]  = xU2*xU3*xx1 + xL1*xL3*xx2 + xL1*xU2*xx3 - xL1*xU2*xU3 - xL1*xU2*xL3;  cc[1] = w- cc[1];
         cc[2]  = xU2*xL3*xx1 + xL1*xL3*xx2 + xU1*xU2*xx3 - xU1*xU2*xL3 - xL1*xU2*xL3;  cc[2] = w- cc[2];
         cc[3]  = xU2*xU3*xx1 + xL1*xU3*xx2 + xL1*xL2*xx3 - xL1*xU2*xU3 - xL1*xL2*xU3;  cc[3] = w- cc[3];
         cc[4]  = xL2*xU3*xx1 + xU1*xU3*xx2 + xL1*xL2*xx3 - xU1*xL2*xU3 - xL1*xL2*xU3;  cc[4] = w- cc[4];
         cc[5]  = xL2*xL3*xx1 + xU1*xU3*xx2 + xU1*xL2*xx3 - xU1*xL2*xU3 - xU1*xL2*xL3;  cc[5] = w- cc[5];

       //if(vub[v1]*vlb[v2]*vlb[v3] + vlb[v1]*vub[v2]*vub[v3]
         // >= vlb[v1]*vub[v2]*vlb[v3] + vub[v1]*vlb[v2]*vub[v3] &&
         // vub[v1]*vlb[v2]*vlb[v3] + vlb[v1]*vub[v2]*vub[v3]
         // >= vub[v1]*vub[v2]*vlb[v3] + vlb[v1]*vlb[v2]*vub[v3]) {

           Expression theta1c = xL1*xU2*xL3 - xU1*xL2*xL3 - xL1*xU2*xU3 + xL1*xL2*xU3;
           Expression theta2c = xU1*xU2*xL3 - xL1*xU2*xU3 - xU1*xL2*xL3 + xU1*xL2*xU3;

         cc[6]  = xL2*xL3*xx1 + xL1*xL3*xx2 + xL1*xL2*xx3 - 2.*xL1*xL2*xL3;  cc[6] = w- cc[6];
         cc[7]  = xU2*xU3*xx1 + xU1*xU3*xx2 + xU1*xU2*xx3 - 2.*xU1*xU2*xU3;  cc[7] = w- cc[7];
         cc[8]  = xL2*xU3*xx1 + xL1*xU3*xx2 + xU1*xL2*xx3 - xL1*xL2*xU3 - xU1*xL2*xU3;  cc[8] = w- cc[8];
         cc[9]  = xU2*xL3*xx1 + xU1*xL3*xx2 + xL1*xU2*xx3 - xL1*xU2*xL3 - xU1*xU2*xL3;  cc[9] = w- cc[9];
         cc[10] = (theta1c/(xL1-xU1))*xx1 + xL1*xU3*xx2 + xL1*xU2*xx3 + (-(theta1c*xU1)/(xL1-xU1))
                  - xL1*xU2*xL3 - xL1*xL2*xU3 + xU1*xL2*xL3;  cc[10] = w- cc[10];
         cc[11] = (theta2c/(xU1-xL1))*xx1 + xU1*xL3*xx2 + xU1*xL2*xx3 + (-(theta2c*xL1)/(xU1-xL1))
                  - xU1*xU2*xL3 - xU1*xL2*xU3 + xL1*xU2*xU3;  cc[11] = w- cc[11];
       }

      } // end if case 10

      // get the bound on the constraints to be added 
      // simplify constraints
      double bnd[cc.size()];    
      for(int ii = 0; ii < cc.size(); ii++) {
          bnd[ii] = cc[ii]->GetConstantPart();  cc[ii] = cc[ii] - bnd[ii];
          Simplify(&cc[ii]);
          defcons.push_back(cc[ii]);
      }

      // add constraints
      for(int ii = 0; ii < defcons.size(); ii++) {
          if(ii < defcons.size()/2) {
            if(cf > 0) InProb->AddConstraint(cn, defcons[ii], - bnd[ii], ROSEINFINITY);
            if(cf < 0) InProb->AddConstraint(cn, defcons[ii], -ROSEINFINITY, - bnd[ii]);
          } else {
            if(cf > 0) InProb->AddConstraint(cn, defcons[ii], -ROSEINFINITY, - bnd[ii]);
            if(cf < 0) InProb->AddConstraint(cn, defcons[ii], - bnd[ii], ROSEINFINITY);
          }
      }
      defcons.erase(defcons.begin(), defcons.end());
      cc.erase(cc.begin(), cc.end());

     //cout << *InProb;

      } // end if (trilin==1 || trilin2=1)
      vidx.erase(vidx.begin(),vidx.end());

    //} // end if(op==product)

    } else if(op == POWER || (op == VAR && (expo!=0 && expo!=1))) {
//      cout << "Power " << endl;

      // find variable involved in power
      vector<int> vidx;
      Expression xs(1, 1, "x");
      Expression schema = xs^expo;

      InProb->GetConstraintLI(j)->Function->GetVarIndicesInSchema(vidx, schema);

      int v1 = vidx[0];
      string xname = InProb->GetVariableLI(v1)->Name; 
      Expression xx(1, v1, xname);

      Expression xL1(vlb[v1-1]);
      Expression xU1(vub[v1-1]);

      // x^2k, x^(2k+1) with range not including 0
      if(int(expo)%2==0 || 
         (int(expo)%2==1 && vlb[v1-1]<=0 && vub[v1-1]<=0)|| (int(expo)%2==1 && vlb[v1-1]>=0 && vub[v1-1]>=0) ) {

        Expression xLUm = (xU1+xL1)/2.;
        Expression xLm = (xL1+xLUm)/2.;
        Expression xUm = (xU1+xLUm)/2.;
        Expression dr1 = coeff*expo*(xL1^(expo-1));
        Expression c1_1 = coeff*(xL1^expo) + dr1*(xx-xL1);
        Expression dr2 = coeff*expo*(xU1^(expo-1));
        Expression c1_2 = coeff*(xU1^expo) + dr2*(xx-xU1);
        Expression dr3 = coeff*expo*(xLm^(expo-1));
        Expression c1_3 = coeff*(xLm^expo) + dr3*(xx-xLm);
        Expression dr4 = coeff*expo*(xUm^(expo-1));
        Expression c1_4 = coeff*(xUm^expo) + dr4*(xx-xUm);
        Expression dr5 = coeff*expo*(xLUm^(expo-1));
        Expression c1_5 = coeff*(xLUm^expo) + dr5*(xx-xLUm);

        c1_1 = e - c1_1;
        c1_2 = e - c1_2;
        c1_3 = e - c1_3;
        c1_4 = e - c1_4;
        c1_5 = e - c1_5;

        Expression c2;
        Expression cc1 = coeff*(xU1^expo);
        Expression cc2 = coeff*(xL1^expo);
        Expression cc3 = cc1-cc2;
        Expression cc4 = cc3/(xU1-xL1);
        Expression cc5 = cc4*(xx-xL1);
        c2 = cc2 + cc5;

        c2 = e - c2;

        double bnd[6];    
        bnd[0] = c1_1->GetConstantPart();  c1_1 = c1_1 -bnd[0];
        bnd[1] = c1_2->GetConstantPart();  c1_2 = c1_2 -bnd[1];
        bnd[2] = c1_3->GetConstantPart();  c1_3 = c1_3 -bnd[2];
        bnd[3] = c1_4->GetConstantPart();  c1_4 = c1_4 -bnd[3];
        bnd[4] = c1_5->GetConstantPart();  c1_5 = c1_5 -bnd[4];
        bnd[5] = c2->GetConstantPart();  c2 = c2 -bnd[5];

        Simplify(&c1_1);
        Simplify(&c1_2);
        Simplify(&c1_3);
        Simplify(&c1_4);
        Simplify(&c1_5);
        Simplify(&c2);
        defcons.push_back(c1_1);
        defcons.push_back(c1_2);
        defcons.push_back(c1_3);
        defcons.push_back(c1_4);
        defcons.push_back(c1_5);
        defcons.push_back(c2);
    
        int sgn;
        if(vlb[v1-1]>=0 && vub[v1-1]>=0) sgn = 1; 
        if(vlb[v1-1]<=0 && vub[v1-1]<=0) sgn = -1; 

        // add constraints
        if(int(expo)%2==0) {
        for(int ii = 0; ii < defcons.size(); ii++) {
           if(coeff >0) {
              if(ii < defcons.size()-1) {
	         InProb->AddConstraint(cn, defcons[ii], - bnd[ii], ROSEINFINITY);
               } else {
	         InProb->AddConstraint(cn, defcons[ii], -ROSEINFINITY, - bnd[ii]);
               } 
           } else { 
              if(ii < defcons.size()-1) {
	         InProb->AddConstraint(cn, defcons[ii], -ROSEINFINITY, - bnd[ii]);
               } else {
	         InProb->AddConstraint(cn, defcons[ii], - bnd[ii], ROSEINFINITY);
               } 
           }
        }
        } else { // int(expo%2==1)---- x^3--------
        for(int ii = 0; ii < defcons.size(); ii++) {
           if((coeff>0 && sgn==1)||(coeff<0 && sgn==-1)) { // convex 
              if(ii < defcons.size()-1) {
	         InProb->AddConstraint(cn, defcons[ii], - bnd[ii], ROSEINFINITY);
               } else {
	         InProb->AddConstraint(cn, defcons[ii], -ROSEINFINITY, - bnd[ii]);
               } 
           } else { // concave
              if(ii < defcons.size()-1) {
	         InProb->AddConstraint(cn, defcons[ii], -ROSEINFINITY, - bnd[ii]);
               } else {
	         InProb->AddConstraint(cn, defcons[ii], - bnd[ii], ROSEINFINITY);
               } 
           }
        }
        }
        defcons.erase(defcons.begin(), defcons.end());


      // x^(2k+1) with range including 0
      //} else if(int(expo)%2==1 && vlb[v1-1]<=0 && vub[v1-1]>=0) {
      } else if(expo==3 && vlb[v1-1]<=0 && vub[v1-1]>=0) {

        int k = 1;
        double rk = -0.5;
        double one = 1.;
        double rrk = (pow(rk,expo)-1)/(rk-1);
        double cpnt = rk*vlb[v1-1];
        double dpnt = rk*vub[v1-1];
      
        Expression c1, c2, c3, c4;
        Expression cc1 = (xL1^expo);
        Expression cc2 = (xU1^expo);
        Expression cc3 = xx/xL1 - one;
        Expression cc4 = xx/xU1 - one;
        Expression cc5 = one + rrk*cc3;
        Expression cc6 = one + rrk*cc4;

        Expression ccc1 = xL1^(expo-1);
        Expression ccc2 = xU1^(expo-1);

        c1 = coeff*cc1*cc5;
        c3 = coeff*cc2*cc6;
     
        c2 = coeff*(expo*ccc2*xx - (expo-1)*(xU1^expo));
        c4 = coeff*(expo*ccc1*xx - (expo-1)*(xL1^expo));
   
        c1 = e - c1;
        c2 = e - c2;
        c3 = e - c3;
        c4 = e - c4;

        double bnd[5];    
        bnd[0] = c1->GetConstantPart();  c1 = c1 -bnd[0];
        bnd[1] = c2->GetConstantPart();  c2 = c2 -bnd[1];
        bnd[2] = c3->GetConstantPart();  c3 = c3 -bnd[2];
        bnd[3] = c4->GetConstantPart();  c4 = c4 -bnd[3];

        Simplify(&c1);
        Simplify(&c2);
        Simplify(&c3);
        Simplify(&c4);
        defcons.push_back(c1);
        defcons.push_back(c2);
        defcons.push_back(c3);
        defcons.push_back(c4);

        if(coeff*cpnt > coeff*vub[v1-1] && coeff*dpnt > coeff*vlb[v1-1] ||
           coeff*cpnt < coeff*vub[v1-1] && coeff*dpnt < vlb[v1-1]) {
           Expression c5 = cc1 + ((cc2 - cc1)/(xU1-xL1))*(xx-xL1);
           c5 = c5*coeff;
           c5 = e - c5;
           bnd[4] = c5->GetConstantPart();  c5 = c5 -bnd[4];
           Simplify(&c5);
           defcons.push_back(c5);
        }


        // add constraints
        if(cpnt <= vub[v1-1] && dpnt >= vlb[v1-1]) {
          for(int ii = 0; ii < defcons.size(); ii++) {
             if(ii < defcons.size()/2) {
	        if(coeff > 0) InProb->AddConstraint(cn, defcons[ii], - bnd[ii], ROSEINFINITY);
	        if(coeff < 0) InProb->AddConstraint(cn, defcons[ii], -ROSEINFINITY, - bnd[ii]);
             } else {
	        if(coeff > 0) InProb->AddConstraint(cn, defcons[ii], -ROSEINFINITY, - bnd[ii]);
	        if(coeff < 0) InProb->AddConstraint(cn, defcons[ii], - bnd[ii], ROSEINFINITY);
             }
          }
        } else if (coeff*cpnt > coeff*vub[v1-1] && coeff*dpnt > coeff*vlb[v1-1]) {
	      InProb->AddConstraint(cn, defcons[4], - bnd[4], ROSEINFINITY);
	      InProb->AddConstraint(cn, defcons[2], -ROSEINFINITY, - bnd[2]);
	      InProb->AddConstraint(cn, defcons[3], -ROSEINFINITY, - bnd[3]);
        } else if (coeff*cpnt < coeff*vub[v1-1] && coeff*dpnt < coeff*vlb[v1-1]) {
	      InProb->AddConstraint(cn, defcons[0], - bnd[0], ROSEINFINITY);
	      InProb->AddConstraint(cn, defcons[4], -ROSEINFINITY, - bnd[4]);
	      InProb->AddConstraint(cn, defcons[1], - bnd[1], ROSEINFINITY);
        }              
        defcons.erase(defcons.begin(), defcons.end());

      } // end if(expo%2)
      //vidx.erase(vidx.begin(),vidx.end());
    } // end if (op)

  } // end for j 


  // simplify the problem if possible
  // this is necessary to simplify constraints like w < 1
  // before bilinear relax. where Function->GetNode(1) is required
  InProb->Simplifier(false);

  //cout << "Simplifier" << endl;
 // cout << *InProb;


  // Bilinear relaxation
  //

  NumberOfConstraints = InProb->GetNumberOfConstraints(); 
  for(int j = 1; j <= NumberOfConstraints; j++) {

    int op = InProb->GetConstraintLI(j)->Function->GetNode(1)->GetOpType();
    double expo = InProb->GetConstraintLI(j)->Function->GetNode(1)->GetExponent();
    double coeff = InProb->GetConstraintLI(j)->Function->GetNode(1)->GetCoeff();
    Expression e = - InProb->GetConstraintLI(j)->Function->GetLinearPart();

    int size = InProb->GetConstraintLI(j)->Function->GetNode(1)->GetSize();
    if(op == PRODUCT && size == 2) {
//      cout << "Product - bilinear" << endl;

        // find variables involved in bilinear terms
        vector<int> vidx;
        Expression xs(1, 1, "x"), ys(1, 2, "y");
        Expression schema = xs*ys;
        
        InProb->GetConstraintLI(j)->Function->GetVarIndicesInSchema(vidx, schema);

        string xname = InProb->GetVariableLI(vidx[0])->Name; 
        string yname = InProb->GetVariableLI(vidx[1])->Name; 

        double cf = 1;
        double coeff1 = InProb->GetConstraintLI(j)->Function->GetNode(1)->GetNode(0)->GetCoeff();
        if (coeff != 1) cf = cf*coeff;    // case  coeff*(xx*yy)
        if (coeff1 != 1) cf = cf*coeff1; // case  ((coeff1*xx)*yy)      

        Expression xx(cf, vidx[0], xname), yy(1, vidx[1], yname);

        Expression xL1(cf*vlb[vidx[0]-1]);
        Expression xU1(cf*vub[vidx[0]-1]);
        Expression xL2(vlb[vidx[1]-1]);
        Expression xU2(vub[vidx[1]-1]);

        Expression c1 = xL1*yy + xL2*xx - xL1*xL2;  c1 = e- c1;
        Expression c2 = xU1*yy + xU2*xx - xU1*xU2;  c2 = e- c2;
        Expression c3 = xL1*yy + xU2*xx - xL1*xU2;  c3 = e- c3;
        Expression c4 = xU1*yy + xL2*xx - xU1*xL2;  c4 = e- c4;

        double bbnd[4];    
        bbnd[0] = c1->GetConstantPart();  c1 = c1 -bbnd[0];
        bbnd[1] = c2->GetConstantPart();  c2 = c2 -bbnd[1];
        bbnd[2] = c3->GetConstantPart();  c3 = c3 -bbnd[2];
        bbnd[3] = c4->GetConstantPart();  c4 = c4 -bbnd[3];

        Simplify(&c1);
        Simplify(&c2);
        Simplify(&c3);
        Simplify(&c4);
        defcons.push_back(c1);
        defcons.push_back(c2);
        defcons.push_back(c3);
        defcons.push_back(c4);


        // add constraints
        for(int ii = 0; ii < defcons.size(); ii++) {
            if(ii < defcons.size()/2) {
	       if(cf > 0) InProb->AddConstraint(cn, defcons[ii], - bbnd[ii], ROSEINFINITY);
	       if(cf < 0) InProb->AddConstraint(cn, defcons[ii], -ROSEINFINITY, - bbnd[ii]);
            } else {
	       if(cf > 0) InProb->AddConstraint(cn, defcons[ii], -ROSEINFINITY, - bbnd[ii]);
	       if(cf < 0) InProb->AddConstraint(cn, defcons[ii], - bbnd[ii], ROSEINFINITY);
            }
        }
        defcons.erase(defcons.begin(), defcons.end());
    
        vidx.erase(vidx.begin(),vidx.end());
    }// end if bilinear
  }
  //  cout << "end of Bilinear:" << endl;  
    //cout << *InProb;
  

  // delete constraints added by Rsmith
  for(int i = NumberOfConstraints; i>=1; i--) {
     int idc = InProb->GetConstraintID(i);
     if (!InProb->GetConstraintLI(i)->Function->IsLinear()) {
       //if(idc > nconstr_orig) InProb->DeleteConstraint(idc);
       InProb->DeleteConstraint(idc);
     }
  }
 // cout << "\nDeleteConstraint:" << endl;
//  cout << *InProb;


  // simplify the problem if possible
  InProb->Simplifier(false);
  
  //cout << "\nSimplifier:" << endl;
  //cout << *InProb;

  //cout << "addconst = " << addconst << endl;

  // output in .ros format
  //ofstream out(OutFile.c_str());
  //out << *InProb;
  //out.close();

  NumberOfVariables = InProb->GetNumberOfVariables();
  NumberOfConstraints = InProb->GetNumberOfConstraints();

  // output in .mod format
  LocalSolverName = "Rprintdat";
  //LocalSolverName = "Rprintmod";
  Solver* OtherLocalSolver = NewSolver(LocalSolverName);
  OtherLocalSolver->SetProblem(InProb);
  OtherLocalSolver->ReplaceParams(ParameterBlob);
  OtherLocalSolver->Initialize(true);
  OtherLocalSolver->Solve();
  DeleteSolver(OtherLocalSolver);

  return ret;
}

void RQuarticConvexSolver::SetOptimizationDirection(int theoptdir) {
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

void RQuarticConvexSolver::GetSolution(map<int,double>& objfunval, 
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

void RQuarticConvexSolver::GetSolutionLI(vector<double>& objfunval, 
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

void RQuarticConvexSolver::GetSolutionLI(double& objfunval, 
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


void RQuarticConvexSolver::SetStartingPoint(int varID, double sp) {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

void RQuarticConvexSolver::SetStartingPointLI(int localindex, double sp) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

double RQuarticConvexSolver::GetStartingPoint(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

double RQuarticConvexSolver::GetStartingPointLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

void RQuarticConvexSolver::SetVariableLB(int varID, double LB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableLBLI(localindex, LB);
}

double RQuarticConvexSolver::GetVariableLB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void RQuarticConvexSolver::SetVariableUB(int varID, double UB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableUBLI(localindex, UB);
}

double RQuarticConvexSolver::GetVariableUB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}

void RQuarticConvexSolver::SetVariableLBLI(int localindex, double LB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vlb[localindex - 1] = LB;
}

double RQuarticConvexSolver::GetVariableLBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void RQuarticConvexSolver::SetVariableUBLI(int localindex, double UB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vub[localindex - 1] = UB;
}

double RQuarticConvexSolver::GetVariableUBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}


