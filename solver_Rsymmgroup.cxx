/*
** Name:       solver_Rsymmgroup.cxx
** Author:     Leo Liberti
** Source:     GNU C++
** Purpose:    Solver Rsymmgroup implementation
**             This solver outputs a representation of the problem structure
**             so that the symmetry group of the problem can be computed
**             either through Nauty (expressions DAG, for general MINLPs),
**             or through the method given in the COCOA08 paper 
**             (AMPL, for linear problems only)
** License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
** History:    080708 work started
**             10---- spits out MILP description in GAP format
**             111229 outputs MINLP DAG in .gph format
*/

#include <cstdio>
#include <cstring>
#include <cassert>
#include <fstream>
#include <sys/time.h>
#include "newsolver.h"
#include "solver.h"
#include "solver_Rsymmgroup.h"
#include "utils.h"

#define FIRSTOBJ 1
#define EPSILON 1e-4
#define MAXWIDTH 50

// Rsymmgroup CLASS METHODS

SymmgroupSolver::SymmgroupSolver() {
  TheName = "symmgroup";
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
  // 0=nauty(MINLP in DAG format)
  // 1=ampl(MILP for COCOA08 paper)
  // 2=gap(MILP in GAP format)
  // 3=gph(MINLP in DAG format)
  OutType = 0; 
}

SymmgroupSolver::~SymmgroupSolver() { }

void SymmgroupSolver::SetProblem(Problem* p) {
  InProb = p; 
  InitProb = false;
}

bool SymmgroupSolver::CanSolve(int problemtype) {
  // fill in here with information about solver capabilities
}

void SymmgroupSolver::Initialize(bool force = false) {

  if (InitProb) {
    ;
  } else {
    InitProb = true;
  }

  // check the input problem is set to something
  if (InProb == NULL) {
    cerr << "SymmgroupSolver::Initialize(): error: InProb is NULL\n";
    assert(InProb != NULL);
  }

  NumberOfVariables = InProb->GetNumberOfVariables();
  NumberOfConstraints = InProb->GetNumberOfConstraints();
  NumberOfObjectives = InProb->GetNumberOfObjectives();

  // parameters
  ReplaceParams(InProb->GetParamsRef());
  int np = ParameterBlob.GetNumberOfParameters();
  for(int i = 1; i <= np; i++) {
    if (ParameterBlob.GetParameterName(i) == "SymmgroupLocalSolver") {
      LocalSolverName = ParameterBlob.GetParameterStringValue(i);
    } else if (ParameterBlob.GetParameterName(i) == "SymmgroupMaxTime") {
      MaxRunningTime = (double) ParameterBlob.GetParameterIntValue(i);
      assert(MaxRunningTime >= 0);
    } else if (ParameterBlob.GetParameterName(i) == "SymmgroupOutType") {
      OutType = (int) ParameterBlob.GetParameterIntValue(i);
    } else if (ParameterBlob.GetParameterName(i) == "SymmgroupEpsilon") {
      Epsilon = ParameterBlob.GetParameterDoubleValue(i);
      assert(Epsilon > 0);
    } else if (ParameterBlob.GetParameterName(i) == "SymmgroupQuiet" ||
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

  if (OutType < 0 || OutType > 3) {
    cerr << "SymmgroupSolver::Initialize(): error: OutType (" 
	 << OutType << ") out of bounds\n";
    exit(20);
  }

  if (OutType == 0) {
    OutFile = "Rsymmgroup_out.nauty";
  } else if (OutType == 1) {
    OutFile = "Rsymmgroup_out.dat";
  } else if (OutType == 2) {
    OutFile = "Rsymmgroup_out.gap";
  } else if (OutType == 3) {
    OutFile = "Rsymmgroup_out.gph";
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

int SymmgroupSolver::Solve(bool reinitialize) {

  int ret = 0;
  if (reinitialize || !InitProb) {
    Initialize();
  }

  // do it
  if (!Quiet) {
    cout << "SymmgroupSolver: writing ";
    if (OutType == 0) {
      cout << "MINLP DAG in nauty format (dreadnaut < " << OutFile << ")";
    } else if (OutType == 1) {
      cout << "MILP as AMPL .dat file (use with fpp.mod [COCOA08]):" 
	   << OutFile;
    } else if (OutType == 2) {
      cout << "MILP data in GAP format: " << OutFile;
    } else if (OutType == 3) {
      cout << "MINLP DAG in .gph format: " << OutFile;
    }
    cout << "\n";
  }
  ofstream out(OutFile.c_str());
  Constraint* theCon = NULL;

  if (OutType == 0 || OutType == 3) {

    // MINLP DAG output

    if (OutType == 0) {
      // nauty output to execute with nauty (dreadnaut shell)
      out << "! expression DAG automatically generated by ROSE's Rsymmgroup\n";
      out << "! P instance representation for symmgroup" << endl;
      out << "! P has " << NumberOfVariables << " vars, "
	  << NumberOfObjectives << " objs, "
	  << NumberOfConstraints << " constrs" << endl << endl;
    } else if (OutType == 3) {
      // .gph output
      out << "# .gph file automatically generated by ROSE (MINLP DAG)" << endl;
    }
    // the DAG in adjacency lists
    map<int,vector<int> > DAG;
    // vertex oplabels
    map<int,int> VOp;
    // vertex colors
    map<int,set<int> > VColor;
    // startnode = maxvarindex
    int startnode = 0;
    for(int i = 1; i <= NumberOfVariables; i++) {
      if (startnode < InProb->GetVariableID(i)) {
	startnode = InProb->GetVariableID(i);
      }
    }
    int topnode = 1;
    int totalops = ERROR + 1; // ERROR is the last OP ID, see Ev3/common.h
    // compute maximum *explicit* recursion level across all expressions
    // (explicit means: nodes like "coeff*node^expon" are reformulated to
    //  "simple" nodes in this reckoning)
    int level = 0;
    int maxlevel = 0; 
    for(int i = 1; i <= NumberOfObjectives; i++) {
      Objective* theObj = InProb->GetObjectiveLI(i);
      level = theObj->Function->GetNumberOfExplicitLevels();
      if (level > maxlevel) {
	maxlevel = level;
      }
    }
    for(int i = 1; i <= NumberOfConstraints; i++) {
      theCon = InProb->GetConstraintLI(i);
      level = theCon->Function->GetNumberOfExplicitLevels();
      if (level > maxlevel) {
	maxlevel = level;
      }
    }
    // find all constants across expressions
    map<double,int,DoubleLessThan> ConstantSymm;
    for(int i = 1; i <= NumberOfObjectives; i++) {
      Objective* theObj = InProb->GetObjectiveLI(i);
      theObj->Function->GetNumberOfDistinctConstantsRecursive(ConstantSymm);
    }
    for(int i = 1; i <= NumberOfConstraints; i++) {
      theCon = InProb->GetConstraintLI(i);
      theCon->Function->GetNumberOfDistinctConstantsRecursive(ConstantSymm);
    }
    // constant symmetry classes (colours) according to value equality
    int constcounter = 0;
    for(map<double,int>::iterator mi = ConstantSymm.begin();
	mi != ConstantSymm.end(); mi++) {
      constcounter++;
      mi->second = constcounter;
    }
    // variables symmetry classes (colours) according to limits and integrality
    map<pair<bool,pair<double,double> >, set<int> > VC;
    pair<bool,pair<double,double> > p2;
    for(int i = 1; i <= NumberOfVariables; i++) {
      p2.first = integrality[i-1];
      p2.second.first = vlb[i-1];
      p2.second.second = vub[i-1];
      VC[p2].insert(i);
    }
    int varsymmclasses = VC.size() + 1;
    map<int,string> VarNameInv;
    int vcounter;
    for(vcounter = 1; vcounter <= NumberOfVariables; vcounter++) {
      VarNameInv[vcounter] = InProb->GetVariable(vcounter)->Name;
    }
    map<int,int> VarSymm;
    vcounter = 1;
    for(map<pair<bool,pair<double,double> >, set<int> >::iterator mit =
	  VC.begin(); mit != VC.end(); mit++) {
      for(set<int>::iterator sit = mit->second.begin(); 
	  sit != mit->second.end(); sit++) {
	VarSymm[InProb->GetVariableID(*sit)] = vcounter;
      }
      vcounter++;
    }
    // maximum color index on each expression
    int maxcolor = maxlevel * (ERROR + 1);
    // objectives symmetry classes based on number of nodes
    map<int,int> ObjSymm;
    map<int,set<int> > OS2;
    // compute map symmclass -> objid
    for(int i = 1; i <= NumberOfObjectives; i++) {
      Objective* theObj = InProb->GetObjectiveLI(i);
      OS2[theObj->Function->GetNumberOfNodes()].insert(i);
    }
    // now invert the map
    int symmclasscounter = 1;
    for(map<int,set<int> >::iterator mit = OS2.begin(); 
	mit != OS2.end(); mit++) {
      for(set<int>::iterator si = mit->second.begin(); si != mit->second.end();
	  si++) {
	ObjSymm[*si] = symmclasscounter;
      }
      symmclasscounter++;
    }
    // do objectives
    level = 1;
    for(int i = 1; i <= NumberOfObjectives; i++) {
      Objective* theObj = InProb->GetObjectiveLI(i);
      int thetopnode = startnode + topnode;
      theObj->Function->AddExpressionTreeToDAG(DAG, VOp, VColor, 
					       ObjSymm, i,
					       ConstantSymm, 
					       varsymmclasses, 
					       VarSymm, 
					       maxcolor, 
					       startnode, topnode, 
					       totalops, level);
      if (DAG[thetopnode].size() == 0 && !(theObj->Function->IsConstant())) {
	// this is necessary in order to deal with single-variable 
	// objectives --- it's not necessary for constraints
	// because those are already colored by direction and RHS
	VColor[maxcolor + 1].insert(thetopnode);
      }
      topnode++;
    }

    // constraints symmetry classes based on number of nodes and LHS/RHS consts
    map<int,int> ConstrSymm;
    map<pair<int,pair<double,double> >, set<int> > CS2;
    pair<int,pair<double,double> > p;
    // compute map symmclass -> constrid
    for(int i = 1; i <= NumberOfConstraints; i++) {
      theCon = InProb->GetConstraintLI(i);
      p.first = theCon->Function->GetNumberOfNodes();
      p.second.first = theCon->LB;
      p.second.second = theCon->UB;
      CS2[p].insert(i);
    }
    // now invert the map
    for(map<pair<int,pair<double,double> >, set<int> >::iterator 
	  mit = CS2.begin(); mit != CS2.end(); mit++) {
      for(set<int>::iterator si = mit->second.begin(); si != mit->second.end();
	  si++) {
	ConstrSymm[*si] = symmclasscounter;
      }
      symmclasscounter++;
    }

    // do constraints
    set<int> ConstraintNodeID;
    level = 1;
    for(int i = 1; i <= NumberOfConstraints; i++) {
      theCon = InProb->GetConstraintLI(i);
      ConstraintNodeID.insert(startnode + topnode);
      theCon->Function->AddExpressionTreeToDAG(DAG, VOp, VColor, 
					       ConstrSymm, i, 
					       ConstantSymm,  
					       varsymmclasses, VarSymm, 
					       maxcolor, startnode, topnode, 
					       totalops, level);
      topnode++;
    }

    // output Rsymmgroup_out.[nauty|gph]
    int Vsize = DAG.size();
    int vid;
    int vertices = 0;
    map<int,vector<int> >::iterator mit;
    vector<int>::iterator vit;
    // output varnames
    for(map<int,string>::iterator it = VarNameInv.begin(); 
	it != VarNameInv.end(); it++) {
      if (OutType == 0) {
	out << "!";
      } else if (OutType == 3) {
	out << "#";
      }
      out << " varname " << it->first << " = " << it->second << endl;
    }
    // output constraint topmost operator node IDs
    if (OutType == 0) {
      out << "!";
    } else if (OutType == 3) {
      out << "#";
    }
    out << " ConstrNodeID := [";
    set<int>::iterator theit = ConstraintNodeID.begin();
    if (theit != ConstraintNodeID.end()) {
      out << *theit;
      theit++;
      for( ; theit != ConstraintNodeID.end(); theit++) { 
	out << "," << *theit;
      }
    }
    out << "];\n";
    // output ops
    for(map<int,int>::iterator it = VOp.begin(); it != VOp.end(); it++) {
      if (OutType == 0) {
	out << "!";
      } else if (OutType == 3) {
	out << "#";
      }
      out << " op(" << it->first << ") = " << it->second << endl;
    }
    // compute number of vertices in DAG
    for(mit = DAG.begin(); mit != DAG.end(); mit++) {
      if (mit->first > vertices) {
	vertices = mit->first;
      }
      for(vit = mit->second.begin(); vit != mit->second.end(); vit++) {
	if (*vit > vertices) {
	  vertices = *vit;
	}
      }
    }
    if (OutType == 0) {
      // output adjacency lists in nauty format
      out << "n=" << vertices << endl << "$=1\ng\n";
      for(int i = 1; i <= vertices; i++) {
	if (i > 1) {
	  out << " ; !" << i-1 << "\n";
	}
	for(vit = DAG[i].begin(); vit != DAG[i].end(); vit++) {
	  out << " " << *vit;
	}
      }
      // output colors in nauty format
      out << ".\n! colors\nf=[\n ";
      map<int,set<int> >::iterator mit2;
      set<int>::iterator vit2;
      for(mit2 = VColor.begin(); mit2 != VColor.end(); mit2++) {
	if (mit2 != VColor.begin()) {
	  // mit2; // why ever did I insert such a null instruction? LEO111229
	  out << " |\n ";
	}
	for(vit2 = mit2->second.begin(); vit2 != mit2->second.end(); vit2++) {
	  if (vit2 != mit2->second.begin()) {
	    out << ",";
	  }
	  out << " " << *vit2;
	}
      }
      out << "\n]\n! execute nauty\nx\n! type orbits\no\n";
    } else if (OutType == 3) {
      // output MINLP dag in .gph format
      out << "Directed;" << endl;
      out << "Vertices:" << endl;
      // invert the colour map
      map<int,int> v2c;
      map<int,set<int> >::iterator mit2;
      set<int>::iterator vit2;
      for(mit2 = VColor.begin(); mit2 != VColor.end(); mit2++) {
	for(vit2 = mit2->second.begin(); vit2 != mit2->second.end(); vit2++) {
	  v2c[*vit2] = mit2->first;
	}
      }
      // out vertices and colours
      for(int i = 0; i < vertices; i++) {
	out << i << " " << v2c[i] << endl;
      }
      // out edges
      out << "Edges:" << endl;
      for(int i = 1; i <= vertices; i++) {
	for(vit = DAG[i].begin(); vit != DAG[i].end(); vit++) {
	  out << i-1 << " " << (*vit) - 1 << " 1 1" << endl;
	}
      }
    }
    // end MINLP dag (nauty/.gph) output

  } else if (OutType == 1) {

    // AMPL .dat output
    if (!InProb->IsProblemLinear()) {
      cerr << "solver_Rsymmgroup.cxx::Solve(): ERROR: symmgroupouttype=1 on (MI)NLP, abort" << endl;
      exit(154);
    }

    // AMPL .dat output for use with fpp.mod (see COCOA08 paper)
    // generalities
    out << "# P instance representation for symmgroup" << endl;
    out << "# P has " << NumberOfVariables << " variables, "
	<< NumberOfObjectives << " objectives, "
	<< NumberOfConstraints << " constraints" << endl << endl;
    out << "param f := " << NumberOfObjectives << ";" << endl;
    out << "param m := " << NumberOfConstraints << ";" << endl;
    out << "param n := " << NumberOfVariables << ";" << endl << endl;
    // variable info
    out << "# variable symmetry" << endl;
    out << "param : xL xU I := " << endl;
    for(int i = 1; i <= NumberOfVariables; i++) {
      out << " " << i << "  " << vlb[i - 1] << " " << vub[i - 1] << " " 
	  << integrality[i - 1] << endl;
    }
    out << ";" << endl << endl;
    // objective info
    out << "# objective(s) symmetry" << endl;
    vector<double> lincoeff;
    vector<int> linvi;
    vector<string> linvn;
    double linc;
    map<pair<int,int>,double> linpart;
    vector<double> linrhs;
    pair<int,int> p;
    for(int i = 1; i <= NumberOfObjectives; i++) {
      Objective* theObj = InProb->GetObjectiveLI(i);
      theObj->Function->GetLinearInfo(lincoeff, linvi, linvn, linc);
      p.first = i;
      for(int j = 0; j < lincoeff.size(); j++) {
	p.second = InProb->GetVarLocalIndex(linvi[j]);
	linpart[p] = lincoeff[j];
      }
      linrhs.push_back(linc);
    }
    out << "param F :";
    for(int j = 1; j <= NumberOfVariables; j++) {
      out << " " << j;
    }
    out << " :=" << endl;
    for(int i = 1; i <= NumberOfObjectives; i++) {
      out << " " << i << " ";
      p.first = i;
      for(int j = 1; j <= NumberOfVariables; j++) {
	p.second = j;
	out << " " << linpart[p];
      }
      out << endl;
    }
    out << ";" << endl;
    out << "param : Fdir Fconst := " << endl;
    for(int i = 1; i <= NumberOfObjectives; i++) {
      out << " " << i << " ";
      out << " " << InProb->GetOptimizationDirectionLI(i);
      out << " " << linrhs[i-1] << endl;
    }
    out << ";" << endl;
    // constraint info
    out << "# constraint(s) symmetry" << endl;
    linpart.erase(linpart.begin(),linpart.end());
    linrhs.erase(linrhs.begin(),linrhs.end());
    for(int i = 1; i <= NumberOfConstraints; i++) {
      theCon = InProb->GetConstraintLI(i);
      theCon->Function->GetLinearInfo(lincoeff, linvi, linvn, linc);
      p.first = i;
      for(int j = 0; j < lincoeff.size(); j++) {
	p.second = InProb->GetVarLocalIndex(linvi[j]);
	linpart[p] = lincoeff[j];
      }
      linrhs.push_back(linc);
    }
    out << "param G :";
    for(int j = 1; j <= NumberOfVariables; j++) {
      out << " " << j;
    }
    out << " :=" << endl;
    for(int i = 1; i <= NumberOfConstraints; i++) {
      out << " " << i << " ";
      p.first = i;
      for(int j = 1; j <= NumberOfVariables; j++) {
	p.second = j;
	out << " " << linpart[p];
      }
      out << endl;
    }
    out << ";" << endl;
    out << "param : L U := " << endl;

    // first pass: represent 1e30 with something else
    double infrep = 0;
    for(int i = 1; i <= NumberOfConstraints; i++) {
      theCon = InProb->GetConstraintLI(i);
      if (fabs(theCon->LB + linrhs[i-1]) < 1e30 && 
	  fabs(theCon->LB + linrhs[i-1]) > infrep) {
	infrep = fabs(theCon->LB + linrhs[i-1]) + 1;
      }
      if (fabs(theCon->UB + linrhs[i-1]) < 1e30 && 
	  fabs(theCon->UB + linrhs[i-1]) > infrep) {
	infrep = fabs(theCon->UB + linrhs[i-1]) + 1;
      }
    }
    // second pass
    for(int i = 1; i <= NumberOfConstraints; i++) {
      theCon = InProb->GetConstraintLI(i);
      out << " " << i << "  ";
      if (theCon->LB + linrhs[i-1] <= -1e30) {
	out << -infrep;
      } else {
	out << theCon->LB + linrhs[i - 1];
      }
      out << " ";
      if (theCon->UB + linrhs[i-1] >= 1e30) { 
	out << infrep << endl;
      } else {
	out << theCon->UB + linrhs[i - 1] << endl;
      }
    }
    out << ";" << endl;

  } else if (OutType == 2) {

    // output a MILP as Gap code: an m x n  constraint matrix A, 
    // a cost n-vector c, lhs/rhs m-vectors bL, bU, lower/upper var bound
    // n-vectors xL, xU, integrality n-vector I
    // PRECONDITION: all data must be integer
    if (!InProb->IsProblemLinear()) {
      cerr << "solver_Rsymmgroup.cxx::Solve(): ERROR: symmgroupouttype=2 on (MI)NLP, abort" << endl;
      exit(154);
    }

    // variable info
    out << "# P instance representation in Gap format" << endl;
    out << "xL := [";
    out << float2fraction(vlb[0]);
    for(int i = 1; i <= NumberOfVariables - 1; i++) {
      out << ",";
      if (i % MAXWIDTH == 0) {
	out << "\n";
      }
      out << float2fraction(vlb[i]);
    }
    out << "];\n";
    out << "xU := [";
    out << float2fraction(vub[0]);
    for(int i = 1; i <= NumberOfVariables - 1; i++) {
      out << ",";
      if (i % MAXWIDTH == 0) {
	out << "\n";
      }
      out << float2fraction(vub[i]);
    }
    out << "];\n";
    out << "I := [";
    out << integrality[0];
    for(int i = 1; i <= NumberOfVariables - 1; i++) {
      out << ",";
      if (i % MAXWIDTH == 0) {
	out << "\n";
      }
      out << integrality[i];
    }
    out << "];\n";

    // objective info
    vector<double> lincoeff;
    vector<int> linvi;
    vector<string> linvn;
    double linc;
    map<pair<int,int>,double> linpart;
    vector<double> linrhs;
    pair<int,int> p;
    for(int i = 1; i <= NumberOfObjectives; i++) {
      Objective* theObj = InProb->GetObjectiveLI(i);
      theObj->Function->GetLinearInfo(lincoeff, linvi, linvn, linc);
      p.first = i;
      for(int j = 0; j < lincoeff.size(); j++) {
	p.second = InProb->GetVarLocalIndex(linvi[j]);
	linpart[p] = lincoeff[j];
      }
      linrhs.push_back(linc);
    }
    out << "c := [";
    p.first = 1;
    p.second = 1;
    out << float2fraction(linpart[p]);
    for(int j = 2; j <= NumberOfVariables; j++) {
      p.second = j;
      out << ",";
      if (j % MAXWIDTH == 0) {
	out << "\n";
      }
      out << float2fraction(linpart[p]);
    }
    out << "];\n" << endl;

    // constraints info
    linpart.erase(linpart.begin(),linpart.end());
    linrhs.erase(linrhs.begin(),linrhs.end());
    for(int i = 1; i <= NumberOfConstraints; i++) {
      theCon = InProb->GetConstraintLI(i);
      theCon->Function->GetLinearInfo(lincoeff, linvi, linvn, linc);
      p.first = i;
      for(int j = 0; j < lincoeff.size(); j++) {
	p.second = InProb->GetVarLocalIndex(linvi[j]);
	linpart[p] = lincoeff[j];
      }
      linrhs.push_back(linc);
    }
    out << "A := [";
    p.first = 1;
    p.second = 1;
    out << "[" << float2fraction(linpart[p]);
    for(int j = 2; j <= NumberOfVariables; j++) {
      p.second = j;
      out << ",";
      if (j % MAXWIDTH == 0) {
	out << "\n";
      }
      out << float2fraction(linpart[p]);
    }
    out << "]";
    for(int i = 2; i <= NumberOfConstraints; i++) {
      p.first = i;
      p.second = 1;
      out << ",\n [" << linpart[p];
      for(int j = 2; j <= NumberOfVariables; j++) {
	p.second = j;
	out << ","; 
	if (j % MAXWIDTH == 0) {
	  out << "\n";
	}
	out << float2fraction(linpart[p]);
      }
      out << "]";
    }
    out << "];\n" << endl;
    // first pass: represent 1e30 with something smaller
    double infrep = 0;
    for(int i = 1; i <= NumberOfConstraints; i++) {
      theCon = InProb->GetConstraintLI(i);
      if (fabs(theCon->LB + linrhs[i-1]) < 1e30 && 
	  fabs(theCon->LB + linrhs[i-1]) > infrep) {
	infrep = fabs(theCon->LB + linrhs[i-1]) + 1;
      }
      if (fabs(theCon->UB + linrhs[i-1]) < 1e30 && 
	  fabs(theCon->UB + linrhs[i-1]) > infrep) {
	infrep = fabs(theCon->UB + linrhs[i-1]) + 1;
      }
    }
    out << "bL := [";
    // second pass
    vector<int> coninfL;
    theCon = InProb->GetConstraintLI(1);
    if (theCon->LB + linrhs[0] <= -1e30) {
      coninfL.push_back(-1);
      out << float2fraction(-infrep);
    } else {
      coninfL.push_back(0);
      out << float2fraction(theCon->LB + linrhs[0]);
    }
    for(int i = 2; i <= NumberOfConstraints; i++) {
      theCon = InProb->GetConstraintLI(i);
      out << ",";
      if (i % MAXWIDTH == 0) {
	out << "\n";
      }
      if (theCon->LB + linrhs[i-1] <= -1e30) {
	coninfL.push_back(-1);
	out << float2fraction(-infrep);
      } else {
	coninfL.push_back(0);
	out << float2fraction(theCon->LB + linrhs[i - 1]);
      }
    }
    out << "];\n";
    vector<int> coninfU;
    out << "bU := [";
    theCon = InProb->GetConstraintLI(1);
    if (theCon->UB + linrhs[0] >= 1e30) { 
      coninfU.push_back(1);
      out << float2fraction(infrep);
    } else {
      coninfU.push_back(0);
      out << float2fraction(theCon->UB + linrhs[0]);
    }
    for(int i = 2; i <= NumberOfConstraints; i++) {
      theCon = InProb->GetConstraintLI(i);
      out << ",";
      if (i % MAXWIDTH == 0) {
	out << "\n";
      }
      if (theCon->UB + linrhs[i - 1] >= 1e30) { 
	coninfU.push_back(1);
	out << float2fraction(infrep);
      } else {
	coninfU.push_back(0);
	out << float2fraction(theCon->UB + linrhs[i - 1]);
      }
    }
    out << "];\n" << endl;
    out << "Cdir := [";
    if (coninfL[0] == -1 && coninfU[0] == 0) {
      out << -1;
    } else if (coninfL[0] == 0 && coninfU[0] == 1) {
      out << 1;
    } else if (coninfL[0] == 0 && coninfU[0] == 0) {
      out << 0;
    }
    for(int i = 1; i < NumberOfConstraints; i++) {
      out << ",";
      if (i % MAXWIDTH == 0) {
	out << "\n";
      }
      if (coninfL[i] == -1 && coninfU[i] == 0) {
	out << -1;
      } else if (coninfL[i] == 0 && coninfU[i] == 1) {
	out << 1;
      } else if (coninfL[i] == 0 && coninfU[i] == 0) {
	out << 0;
      }
    }
    out << "];\n" << endl;
  }
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
    cout << "SymmgroupSolver: best solution f* = " << OptimalObjVal << endl;
  }
  InProb->SetOptimalObjectiveValue(FIRSTOBJ, OptimalObjVal);
  InProb->SetSolved(true);    
  if (IsFeasible) {
    InProb->SetFeasible(1);
  } else {
    if (!Quiet) {
      cout << "SymmgroupSolver: this solution is infeasible by at least " 
           << discrepancy << endl;
    }
    InProb->SetFeasible(-1);
  }

  return ret;
}

void SymmgroupSolver::SetOptimizationDirection(int theoptdir) {
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

void SymmgroupSolver::GetSolution(map<int,double>& objfunval, 
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

void SymmgroupSolver::GetSolutionLI(vector<double>& objfunval, 
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

void SymmgroupSolver::GetSolutionLI(double& objfunval, 
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


void SymmgroupSolver::SetStartingPoint(int varID, double sp) {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

void SymmgroupSolver::SetStartingPointLI(int localindex, double sp) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  x[localindex - 1] = sp;
}

double SymmgroupSolver::GetStartingPoint(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

double SymmgroupSolver::GetStartingPointLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return x[localindex - 1];
}

void SymmgroupSolver::SetVariableLB(int varID, double LB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableLBLI(localindex, LB);
}

double SymmgroupSolver::GetVariableLB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void SymmgroupSolver::SetVariableUB(int varID, double UB) {
  int localindex = InProb->GetVarLocalIndex(varID);
  SetVariableUBLI(localindex, UB);
}

double SymmgroupSolver::GetVariableUB(int varID) const {
  int localindex = InProb->GetVarLocalIndex(varID);
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}

void SymmgroupSolver::SetVariableLBLI(int localindex, double LB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vlb[localindex - 1] = LB;
}

double SymmgroupSolver::GetVariableLBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vlb[localindex - 1];
}

void SymmgroupSolver::SetVariableUBLI(int localindex, double UB) {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  vub[localindex - 1] = UB;
}

double SymmgroupSolver::GetVariableUBLI(int localindex) const {
  assert(localindex >= 1 && localindex <= NumberOfVariables);
  return vub[localindex - 1];
}
