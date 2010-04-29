/*
** Name:       solver_Rconvexifier.cxx
** Author:     Leo Liberti & Sonia Cafieri
** Source:     GNU C++
** Purpose:    Solver Rconvexifier implementation
** License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
** History:    080606 work started
*/

#include <cstdio>
#include <cstring>
#include <cassert>
#include <sys/time.h>
#include <fstream>
#include "newsolver.h"
#include "solver.h"
#include "solver_Rconvexifier.h"
#include "utils.h"
#include "Ev3/auxiliary.h"

#define FIRSTOBJ 1
#define EPSILON 1e-4
#define WDEFNAMEBUFSIZE 256

// RCONVEXIFIERSOLVER CLASS METHODS

RconvexifierSolver::RconvexifierSolver() {
  TheName = "Rconvexifier";
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
  OutFile = "Rconvexifier_out.ros";
  OutFileMaps = "Rconvexifier_out.run";
  OutAmpl = false;
}

RconvexifierSolver::~RconvexifierSolver() {
  if (LocalSolver) {
    DeleteSolver(LocalSolver);
  }
}

void RconvexifierSolver::SetProblem(Problem* p) {
  InProb = p; 
  InitProb = false;
}

bool RconvexifierSolver::CanSolve(int problemtype) {
  // fill in here with information about solver capabilities
}

Problem *RconvexifierSolver::GetProblem(double *dlb, double *dub, map<int, int> var_map) {
  //Global indices used

  //Replace all constant variables with constant values
  //But do all variables go from 1 to vlb.size(),thats a local index!!not global
  for (int i = 1; i <= var_map.size(); i++) {
    map <int, int>::iterator lbiter = lbconstmap.find(var_map[i]);
    map <int, int>::iterator ubiter = ubconstmap.find(var_map[i]);

    cout << "Attempting to replace VarIndex " << i << endl;

    if (lbiter != lbconstmap.end()) {
      cout << "LB Found at " << lbiter->second << " replacing with " << dlb[i-1] << endl;

      SetVariableLB(lbiter->second, dlb[i-1]);
      SetVariableUB(lbiter->second, dlb[i-1]);
      SetStartingPoint(lbiter->second, dlb[i-1]);
      
      InProb->SetOptimalVariableValue(lbiter->second, dlb[i-1]);
      InProb->SetVariableLBValue(lbiter->second, dlb[i-1]);
      InProb->SetVariableUBValue(lbiter->second, dlb[i-1]);
    }

    if (ubiter != ubconstmap.end()) {
      cout << "UB Found at " << ubiter->second << " replacing with " << dub[i-1] << endl;
      SetVariableLB(lbiter->second, dub[i-1]);
      SetVariableUB(lbiter->second, dub[i-1]);
      SetStartingPoint(lbiter->second, dub[i-1]);
      
      InProb->SetOptimalVariableValue(ubiter->second, dub[i-1]);
      InProb->SetVariableLBValue(ubiter->second, dub[i-1]);
      InProb->SetVariableUBValue(ubiter->second, dub[i-1]);
    }
  }

  return InProb;
}

void RconvexifierSolver::Initialize(bool force = false) {

  if (InitProb) {
    ;
  } else {
    InitProb = true;
  }

  // check the input problem is set to something
  if (InProb == NULL) {
    cerr << "RconvexifierSolver::Initialize(): error: InProb is NULL\n";
    assert(InProb != NULL);
  }

  NumberOfVariables = InProb->GetNumberOfVariables();
  NumberOfConstraints = InProb->GetNumberOfConstraints();
  NumberOfObjectives = InProb->GetNumberOfObjectives();

  // parameters
  ReplaceParams(InProb->GetParamsRef());
  int np = ParameterBlob.GetNumberOfParameters();
  for(int i = 1; i <= np; i++) {
    if (ParameterBlob.GetParameterName(i) == "RconvexifierLocalSolver") {
      LocalSolverName = ParameterBlob.GetParameterStringValue(i);
    } else if (ParameterBlob.GetParameterName(i) == "ConvexifierOutAmpl") {
      OutAmpl = ParameterBlob.GetParameterBoolValue(i);
    } else if (ParameterBlob.GetParameterName(i) == "ConvexifierMaxTime") {
      MaxRunningTime = (double) ParameterBlob.GetParameterIntValue(i);
      assert(MaxRunningTime >= 0);
    } else if (ParameterBlob.GetParameterName(i) == "ConvexifierEpsilon") {
      Epsilon = ParameterBlob.GetParameterDoubleValue(i);
      assert(Epsilon > 0);
    } else if (ParameterBlob.GetParameterName(i) == "ConvexifierQuiet" ||
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

  // current, best solution, boundsa
/*
  for(int i = 0; i < NumberOfVariables; i++) {
    x.push_back(InProb->GetStartingPointLI(i + 1));
    xstar.push_back(x[i]);
    vlb.push_back(InProb->GetVariableLI(i + 1)->LB);
    vub.push_back(InProb->GetVariableLI(i + 1)->UB);
    StartingPoint.push_back(x[i]);
  }
*/
  f = ROSEINFINITY;
  fstar = ROSEINFINITY;
}

Expression RconvexifierSolver::addLbConstVar(double coeff, int index, 
					     vector <double> *vlb, vector <double> *vub) {
  //varindex is the local variable index
  double lb = (*vlb)[index-1];
  double ub = (*vub)[index-1];
  stringstream varnamess;
  string varname;
  int newvarindex;

  varnamess << InProb->GetVariableLI(index)->Name << "_lb";
  varname = varnamess.str();

  InProb->AddVariable(varname,false,false,true,lb,ub,lb);
  vlb->push_back(lb);
  vub->push_back(ub);
  x.push_back(lb);
  xstar.push_back(lb);
  StartingPoint.push_back(lb);

  newvarindex = InProb->GetNumberOfVariables();

  lbconstmap.insert(pair <int, int> (index, newvarindex));
  nconvar++;

  Expression ret(coeff,newvarindex,varname);

  return ret;
}

Expression RconvexifierSolver::addUbConstVar(double coeff, int index, 
					     vector <double> *vlb, vector <double> *vub) {
  double lb = (*vlb)[index-1];
  double ub = (*vub)[index-1];
  stringstream varnamess;
  string varname;
  int newvarindex;

  varnamess << InProb->GetVariableLI(index)->Name << "_ub";
  varname = varnamess.str();

  InProb->AddVariable(varname,false,false,true,lb,ub,ub);
  vlb->push_back(lb);
  vub->push_back(ub);
  x.push_back(ub);
  xstar.push_back(ub);
  StartingPoint.push_back(ub);

  newvarindex = InProb->GetNumberOfVariables();
  ubconstmap.insert(pair <int, int> (index, newvarindex));
  nconvar++;

  Expression ret(coeff,newvarindex,varname);
  
  return ret;
}

int RconvexifierSolver::Solve(bool reinitialize) {

  int ret = 0;
  if (reinitialize || !InitProb) {
    Initialize();
  }

  // do it

  int nvar_orig = NumberOfVariables;
  int nconstr_orig = NumberOfConstraints;

  // Standardize the problem using Rsmith
  //LocalSolver->Initialize(false);
  LocalSolver->Solve();

  NumberOfVariables = InProb->GetNumberOfVariables();
  NumberOfConstraints = InProb->GetNumberOfConstraints();

  map<int,vector<int> > addedvar2constrLI;
  int addedvarcounter = nvar_orig;
  int addedconcounter = nconstr_orig + 1;
  int stdconstr = 0;

  for(int i = 0; i < NumberOfVariables; i++) {
    vlb.push_back(InProb->GetVariableLI(i + 1)->LB) ;
    vub.push_back(InProb->GetVariableLI(i + 1)->UB) ;
    x.push_back(0);
    xstar.push_back(0);
    StartingPoint.push_back(x[i]);
  }

  vector<Expression> defcons;
  string cn = "c";
 
  int frac = 0;
  int v1f, v2f;


  for(int j = 1; j <= NumberOfConstraints; j++) {

    if (!InProb->GetConstraintLI(j)->Function->IsLinear()) {
      addedvarcounter++;
      // only convexify nonlinear defining constraints
      int size = InProb->GetConstraintLI(j)->Function->GetNode(1)->GetSize();
      int op = InProb->GetConstraintLI(j)->Function->GetNode(1)->GetOpType();
      double expo = 
	InProb->GetConstraintLI(j)->Function->GetNode(1)->GetExponent();
      double coeff=InProb->GetConstraintLI(j)->Function->GetNode(1)->GetCoeff();
      Expression e = - InProb->GetConstraintLI(j)->Function->GetLinearPart();
      
      int v1, v2;

      if(op == FRACTION) {
	//cout << "Fraction" << endl;

	// find variables involved in fraction
	vector<int> vidx;
	Expression xs(1, 1, "x"), ys(1, 2, "y");
	Expression schema = xs/ys;
	Constraint* theCon = InProb->GetConstraintLI(j);
	if (!theCon->Function->IsLinear()) {
	  theCon->Function->GetVarIndicesInSchema(vidx, schema);
	}

	// add variable z = 1/y
#undef sprintf
	string zn = "z";
	char wdefnbuf[WDEFNAMEBUFSIZE];
	int waddcounter = NumberOfVariables + 1;
	sprintf(wdefnbuf, "z%d", waddcounter);
	zn = wdefnbuf;  // name
	Expression z(1,waddcounter,zn);
	double zlb,zub;
	fractionmkrange(1.0, 1.0, vlb[vidx[1]-1], vub[vidx[1]-1], 
			&zlb, &zub); //bounds

	InProb->AddVariable(zn,false,false,zlb,zub,0);
	vlb.push_back(zlb);
	vub.push_back(zub);
	x.push_back(0);
	xstar.push_back(0);
	StartingPoint.push_back(0);

	NumberOfVariables = InProb->GetNumberOfVariables();
	int v4 = NumberOfVariables;

	double cf = 1;
	for (int jj = 1; jj<= j; jj++) {
	  for (int i = 0; i < InProb->GetConstraintLI(jj)->Function->GetSize(); 
	       i++) {
	    if(InProb->GetConstraintLI(jj)->Function->GetNode(i)->GetOpType() 
	       == VAR && 
	       InProb->GetConstraintLI(jj)->Function->GetNode(i)->GetVarIndex()
	       == vidx[1]) {
	    
	      Expression con =InProb->GetConstraintLI(jj)->Function->GetNode(i);
	      double val = 
		InProb->GetConstraintLI(jj)->Function->GetNode(i)->GetValue();
	      cf=InProb->GetConstraintLI(jj)->Function->GetNode(i)->GetCoeff(); 
	      Expression zz(1,v4,"z");
	      Expression cons = cf*(1./zz);
	    
	      InProb->GetConstraintLI(jj)->Function->
		GetNode(i)->SetValue(1/val);
	      InProb->GetConstraintLI(jj)->Function->GetNode(i)->
		ReplaceWithExpression(cons);
	      InProb->GetConstraintLI(jj)->Function->GetNode(i)->
		CreateFastEvalTree();
	      InProb->GetConstraintLI(jj)->FunctionFET = 
		InProb->GetConstraintLI(jj)->Function->GetFastEvalTree();
	    }
	  }
	}

	if(cf !=1 && !OutAmpl) {
	  InProb->Simplifier(false);
	}
     
	// simplify the problem if possible
	if (!OutAmpl) {
	  //InProb->Simplifier(false);
	}

	op = PRODUCT;
	frac = 1;
	v1f = vidx[0];
	v2f = v4;
      }

      if(op == PRODUCT && size == 2) {
	//cout << "Product" << endl;

	int v1, v2;
      
	if (frac == 0) {
	  // find variables involved in bilinear terms
	  vector<int> vidx;
	  Expression xs(1, 1, "x"), ys(1, 2, "y");
	  Expression schema = xs*ys;
	  Constraint* theCon = InProb->GetConstraintLI(j);
	  if (!theCon->Function->IsLinear()) {
	    theCon->Function->GetVarIndicesInSchema(vidx, schema);
	  }
	  v1 = vidx[0];
	  v2 = vidx[1];
	}
	else {
	  v1 = v1f;
	  v2 = v2f;
	}  
     
	string xname = InProb->GetVariableLI(v1)->Name; 
	string yname = InProb->GetVariableLI(v2)->Name; 
      
	double cf = 1;
	double coeff1 = InProb->
	  GetConstraintLI(j)->Function->GetNode(1)->GetNode(0)->GetCoeff();
	if (coeff != 1) {
	  cf = cf*coeff;    // case  coeff*(xx*yy)
	}
	if (coeff1 != 1) {
	  cf = cf*coeff1; // case  ((coeff1*xx)*yy)      
	}
      
	Expression xx(cf, v1, xname), yy(1, v2, yname);
      
	Expression xL1 = addLbConstVar(cf, v1, &vlb, &vub);
	Expression xU1 = addUbConstVar(cf, v1, &vlb, &vub);
	Expression xL2 = addLbConstVar(1.0, v2, &vlb, &vub);
	Expression xU2 = addUbConstVar(1.0, v2, &vlb, &vub);
      
	Expression c1 = xL1*yy + xL2*xx - xL1*xL2;  c1 = e- c1;
	Expression c2 = xU1*yy + xU2*xx - xU1*xU2;  c2 = e- c2;
	Expression c3 = xL1*yy + xU2*xx - xL1*xU2;  c3 = e- c3;
	Expression c4 = xU1*yy + xL2*xx - xU1*xL2;  c4 = e- c4;
      
	double bnd[4];    
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
      

	// add constraints
	for(int ii = 0; ii < defcons.size(); ii++) {
	  if(ii < defcons.size()/2) {
	    if(cf > 0) {
	      InProb->AddConstraint(cn, defcons[ii], -bnd[ii], ROSEINFINITY);
	      addedvar2constrLI[addedvarcounter].push_back(addedconcounter);
	      addedconcounter++;
	    }
	    if(cf < 0) {
	      InProb->AddConstraint(cn, defcons[ii], -ROSEINFINITY,-bnd[ii]);
	      addedvar2constrLI[addedvarcounter].push_back(addedconcounter);
	      addedconcounter++;
	    }
	  } else {
	    if(cf > 0) {
	      InProb->AddConstraint(cn, defcons[ii],-ROSEINFINITY,-bnd[ii]);
	      addedvar2constrLI[addedvarcounter].push_back(addedconcounter);
	      addedconcounter++;
	    }
	    if(cf < 0) {
	      InProb->AddConstraint(cn, defcons[ii], -bnd[ii],ROSEINFINITY);
	      addedvar2constrLI[addedvarcounter].push_back(addedconcounter);
	      addedconcounter++;
	    }
	  }
	}
	defcons.erase(defcons.begin(), defcons.end());
      
      } else if(op == POWER || (op == VAR && (expo!=0 && expo!=1))) {
      
	// find variable involved in power
	vector<int> vidx;
	Expression xs(1, 1, "x");
	Expression schema = xs^expo;
	Constraint* theCon = InProb->GetConstraintLI(j);
	if (!theCon->Function->IsLinear()) {
	  theCon->Function->GetVarIndicesInSchema(vidx, schema);
	}
      
	int v1 = vidx[0];
	string xname = InProb->GetVariableLI(v1)->Name; 
	Expression xx(1, v1, xname);
      
	Expression xL1 = addLbConstVar(1.0, v1, &vlb, &vub);
	Expression xU1 = addUbConstVar(1.0, v1, &vlb, &vub);
      

	// x^2k, x^(2k+1) with range not including 0
	if(int(expo)%2==0 || 
	   (int(expo)%2==1 && vlb[v1-1]<=0 && vub[v1-1]<=0) || 
	   (int(expo)%2==1 && vlb[v1-1]>=0 && vub[v1-1]>=0) ) {
	        
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
    
	  // add constraints
	  for(int ii = 0; ii < defcons.size(); ii++) {
	    if(coeff > 0) {
	      if(ii < defcons.size()-1) {
		InProb->AddConstraint(cn, defcons[ii], - bnd[ii], ROSEINFINITY);
		addedvar2constrLI[addedvarcounter].push_back(addedconcounter);
		addedconcounter++;
	      } else {
		InProb->AddConstraint(cn, defcons[ii],-ROSEINFINITY, -bnd[ii]);
		addedvar2constrLI[addedvarcounter].push_back(addedconcounter);
		addedconcounter++;
	      } 
	    } else {
	      if(ii < defcons.size()-1) {
		InProb->AddConstraint(cn, defcons[ii],-ROSEINFINITY, -bnd[ii]);
		addedvar2constrLI[addedvarcounter].push_back(addedconcounter);
		addedconcounter++;
	      } else {
		InProb->AddConstraint(cn, defcons[ii], -bnd[ii], ROSEINFINITY);
		addedvar2constrLI[addedvarcounter].push_back(addedconcounter);
		addedconcounter++;
	      } 
	    }
	  }
	  defcons.erase(defcons.begin(), defcons.end());
        
	} else if(int(expo)%2==1 && vlb[v1-1]<=0 && vub[v1-1]>=0) {
	  // x^(2k+1) with range including 0

	  int k = int((expo-1)/2);
	  if (k > 10) {
	    cerr << "RconvexifierSolver::Solve(): error: exponent greater than 21 - case not implemented \n";
	    exit(0);
	  }
	  double rk[10];
	  rk[0] = -0.5; rk[1] = -0.6058295862; rk[2] = -0.6703320476;
	  rk[3] = -0.7145377272; rk[4] = -0.7470540749; rk[5] = -0.7721416355;
	  rk[6] = -0.7921778546; rk[7] = -0.8086048979; rk[8] = -0.8223534102; 
	  rk[9] = -0.8340533676; 

	  double one = 1.;
	  double rrk = (pow(rk[k-1],expo)-1)/(rk[k-1]-1);
	  double cpnt = rk[k-1]*vlb[v1-1];
	  double dpnt = rk[k-1]*vub[v1-1];
      
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
		if(coeff > 0) {
		  InProb->AddConstraint(cn, defcons[ii],-bnd[ii],ROSEINFINITY);
		  addedvar2constrLI[addedvarcounter].push_back(addedconcounter);
		  addedconcounter++;
		}
		if(coeff < 0) {
		  InProb->AddConstraint(cn, defcons[ii],-ROSEINFINITY,-bnd[ii]);
		  addedvar2constrLI[addedvarcounter].push_back(addedconcounter);
		  addedconcounter++;
		}
	      } else {
		if(coeff > 0) {
		  InProb->AddConstraint(cn, defcons[ii],-ROSEINFINITY,-bnd[ii]);
		  addedvar2constrLI[addedvarcounter].push_back(addedconcounter);
		  addedconcounter++;
		}
		if(coeff < 0) {
		  InProb->AddConstraint(cn, defcons[ii],-bnd[ii], ROSEINFINITY);
		  addedvar2constrLI[addedvarcounter].push_back(addedconcounter);
		  addedconcounter++;
		}
	      }
	    }
	  } else if (coeff*cpnt > coeff*vub[v1-1] && 
		     coeff*dpnt > coeff*vlb[v1-1]) {
	    InProb->AddConstraint(cn, defcons[4], - bnd[4], ROSEINFINITY);
	    addedvar2constrLI[addedvarcounter].push_back(addedconcounter);
	    addedconcounter++;

	    InProb->AddConstraint(cn, defcons[2], -ROSEINFINITY, - bnd[2]);
	    addedvar2constrLI[addedvarcounter].push_back(addedconcounter);
	    addedconcounter++;

	    InProb->AddConstraint(cn, defcons[3], -ROSEINFINITY, - bnd[3]);
	    addedvar2constrLI[addedvarcounter].push_back(addedconcounter);
	    addedconcounter++;
	  } else if (coeff*cpnt < coeff*vub[v1-1] && 
		     coeff*dpnt < coeff*vlb[v1-1]) {
	    InProb->AddConstraint(cn, defcons[0], - bnd[0], ROSEINFINITY);
	    addedvar2constrLI[addedvarcounter].push_back(addedconcounter);
	    addedconcounter++;

	    InProb->AddConstraint(cn, defcons[4], -ROSEINFINITY, - bnd[4]);
	    addedvar2constrLI[addedvarcounter].push_back(addedconcounter);
	    addedconcounter++;

	    InProb->AddConstraint(cn, defcons[1], - bnd[1], ROSEINFINITY);
	    addedvar2constrLI[addedvarcounter].push_back(addedconcounter);
	    addedconcounter++;
	  }              
	  defcons.erase(defcons.begin(), defcons.end());
	
	} // end if(expo%2)
          
      } else {
	addedvarcounter--;
      } // end if(op)
    }
  }
  
  // delete constraints added by Rsmith
  for(int i = NumberOfConstraints; i>=1; i--) {
     int idc = InProb->GetConstraintID(i);
     if (!InProb->GetConstraintLI(i)->Function->IsLinear()) {
       if(idc > nconstr_orig) {
	 InProb->DeleteConstraint(idc);
	 stdconstr++;
       }
     }
  }

  // simplify the problem if possible
  if (!OutAmpl) {
    InProb->Simplifier(false);
  }
  
  // update NumberOfVariables and NumberOfConstraints
  int oldNC = NumberOfConstraints;
  NumberOfVariables = InProb->GetNumberOfVariables();
  NumberOfConstraints = InProb->GetNumberOfConstraints();

  // output in .ros format
  ofstream out(OutFile.c_str());
  out << *InProb;
  out.close();

  if (OutAmpl) {
    // output maps
    ofstream outmaps(OutFileMaps.c_str());
    outmaps << "# Rconvexifier_out.run -- "
	    << "maps added variables to constraint indices in AMPL" << endl;
    outmaps << "set addedvars := " << nvar_orig+1 << ".." << addedvarcounter
	    << ";" << endl;
    outmaps << "set addedvar2constr{addedvars};" << endl;
    for(int j = nvar_orig+1; j <= addedvarcounter; j++) {
      outmaps << "let addedvar2constr[" << j << "] :=";
      for(vector<int>::iterator vi = addedvar2constrLI[j].begin();
	  vi != addedvar2constrLI[j].end(); vi++) {
	outmaps << " " << *vi + (oldNC - nconstr_orig - stdconstr);
      }
      outmaps << ";" << endl;
    }
    outmaps.close();

    // output in .mod format
    LocalSolverName = "Rprintmod";
    Solver* OtherLocalSolver = NewSolver(LocalSolverName);
    OtherLocalSolver->SetProblem(InProb);
    OtherLocalSolver->ReplaceParams(ParameterBlob);
    OtherLocalSolver->Initialize();
    OtherLocalSolver->Solve();
    DeleteSolver(OtherLocalSolver);
  }
  ////////////////////////////////
  // after solution: set optimal values in problem
  IsSolved = true;
  for(int i = 0; i < NumberOfVariables; i++) {
    InProb->SetCurrentVariableValueLI(i + 1, xstar[i]);
  }
  return ret;
}

void RconvexifierSolver::SetOptimizationDirection(int theoptdir) {
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

void RconvexifierSolver::GetSolution(map<int,double>& objfunval, 
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

void RconvexifierSolver::GetSolutionLI(vector<double>& objfunval, 
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

void RconvexifierSolver::GetSolutionLI(double& objfunval, 
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


void RconvexifierSolver::SetStartingPoint(int varID, double sp) {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

void RconvexifierSolver::SetStartingPointLI(int localindex, double sp) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

double RconvexifierSolver::GetStartingPoint(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

double RconvexifierSolver::GetStartingPointLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

void RconvexifierSolver::SetVariableLB(int varID, double LB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableLBLI(localindex, LB);
}

double RconvexifierSolver::GetVariableLB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void RconvexifierSolver::SetVariableUB(int varID, double UB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableUBLI(localindex, UB);
}

double RconvexifierSolver::GetVariableUB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}

void RconvexifierSolver::SetVariableLBLI(int localindex, double LB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vlb[localindex - 1] = LB;
}

double RconvexifierSolver::GetVariableLBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void RconvexifierSolver::SetVariableUBLI(int localindex, double UB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vub[localindex - 1] = UB;
}

double RconvexifierSolver::GetVariableUBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}
