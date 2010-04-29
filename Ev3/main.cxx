/**********************************************************************
* Author:      Leo Liberti                                            *
* Name:        testfbbt.cxx                                           *
* Source:      GNU C++                                                *
* Purpose:     FBBT test main driver                                  *
* History:     091020 0.0 work started                                *
* License:     (C) Leo Liberti, all rights reserved. Code published under the
               Common Public License.
***********************************************************************/

#include <iostream>

#include "expression.h"
#include "parser.h"

using namespace std;

int main(int argc, char** argv) {

  ExpressionParser p;
  Expression e[2];

  int nerr;

  e[0] = p.Parse("2*x - y", nerr);
  e[1] = p.Parse("x - 2*y", nerr);

  map<int,double> xL;
  map<int,double> xU;
  double elb[2];
  double eub[2];
  xL[1] = 0;  xU[1] = 1;
  xL[2] = 0;  xU[2] = 1;
  elb[0] = -LARGE; eub[0] = 0;
  elb[1] = -LARGE; eub[1] = 0;

  int m = 2;
  int n = xL.size();
  map<int,double> XL(xL);
  map<int,double> XU(xU);
  map<int,double> ZL(xL);
  map<int,double> ZU(xU);

  cout << "## fbbt.mod instance .dat file ######################" << endl;
  cout << "## WARNING: read (x) as (1*x)" << endl;
  for (int i = 0; i < m; i++) {
    cout << "# " 
	 << elb[i] << " <= " << e[i]->ToString() << " <= " << eub[i] << endl;
  }
  for(int j = 1; j <= n; j++) {
    cout << "#  x" << j << " in [" << ZL[j] << "," << ZU[j] << "]" << endl;
  }
  cout << "###################################################" << endl;

  // print out AMPL .dat instance
  vector<std::pair<int,int> > edgeList;
  map<int,int> nodeType; 
  map<int,int> varindex;
  map<int,double> coeff;
  vector<int> nNodes;
  map<int,pair<double,double> > bnds;
  cout << "param n := " << n << ";" << endl;
  cout << "param m := " << m << ";" << endl;
  for(int i = 0; i < m; i++) {
    cout << "set V[" << i+1 << "] := ";
    nNodes.push_back
      (e[i]->TreeData(1, edgeList, nodeType, varindex, coeff, bnds));
    for(int v = 1; v <= nNodes[i] - 1; v++) {
      cout << v << ", ";
    }
    cout << nNodes[i] << ";" << endl;
    cout << "set A[" << i+1 << "] := ";
    int a;
    for(a = 0; a < (int) edgeList.size(); a++) {
      cout << "(" << edgeList[a].first << "," << edgeList[a].second << ")";
      if (a < (int) edgeList.size() - 1) {
	cout << ", ";
      }
    }
    cout << ";" << endl;
  }

  cout << "param lambda := " << endl;
  for(int i = 0; i < m; i++) {
    e[i]->TreeData(1, edgeList, nodeType, varindex, coeff, bnds);
    for(int v = 1; v <= nNodes[i]; v++) {
      cout << "  " << i+1 << " " << v << "  " << nodeType[v]+1 << endl;
    }
  }
  cout << ";" << endl;

  cout << "param index := " << endl;
  for(int i = 0; i < m; i++) {
    e[i]->TreeData(1, edgeList, nodeType, varindex, coeff, bnds);
    for(int v = 0; v < nNodes[i]; v++) {
      cout << "  " << i+1 << " " << v+1 << "  " << varindex[v+1] << endl;
    }
  }
  cout << ";" << endl;

  cout << "param a := " << endl;
  for(int i = 0; i < m; i++) {
    e[i]->TreeData(1, edgeList, nodeType, varindex, coeff, bnds);
    for(int v = 0; v < nNodes[i]; v++) {
      cout << "  " << i+1 << " " << v+1 << "  " << coeff[v] << endl;
    }
  }
  cout << ";" << endl;

  cout << "param r := " << endl;
  for(int i = 0; i < m; i++) {
    cout << "  " << i+1 << "  1" << endl;
  }
  cout << ";" << endl;
  
  cout << "param gL := " << endl;
  for(int i = 0; i < m; i++) {
    cout << "  " << i+1 << "  " << elb[i] << endl;
  }
  cout << ";" << endl;

  cout << "param gU := " << endl;
  for(int i = 0; i < m; i++) {
    cout << "  " << i+1 << "  " << eub[i] << endl;
  }
  cout << ";" << endl;

  cout << "param x0L := " << endl;
  for(int j = 1; j <= n; j++) {
    cout << "  " << j << "  " << xL[j] << endl;
  }
  cout << ";" << endl;

  cout << "param x0U := " << endl;
  for(int j = 1; j <= n; j++) {
    cout << "  " << j << "  " << xU[j] << endl;
  }
  cout << ";" << endl;

  cout << "param y00L := " << endl;
  for(int i = 0; i < m; i++) {
    e[i]->TreeData(1, edgeList, nodeType, varindex, coeff, bnds);
    for(int v = 0; v < nNodes[i]; v++) {
      cout << "  " << i+1 << " " << v+1 << "  " << bnds[v].first << endl;
    }
  }
  cout << ";" << endl;
  
  cout << "param y00U := " << endl;
  for(int i = 0; i < m; i++) {
    e[i]->TreeData(1, edgeList, nodeType, varindex, coeff, bnds);
    for(int v = 0; v < nNodes[i]; v++) {
      cout << "  " << i+1 << " " << v+1 << "  " << bnds[v].second << endl;
    }
  }
  cout << ";" << endl;

  cout << "## FBBT gives: ###################################" << endl;
  cout << "# FBBT start" << endl;
  for(int j = 1; j <= n; j++) {
    cout << "#  x" << j << " in [" << ZL[j] << "," << ZU[j] << "]" << endl;
  }

  // FBBT
  int index;
  int iteration = 1;
  map<int,double>::iterator mi;
  const double epsilon = 1e-6;
  double discrepancy = 2*epsilon;
  while(discrepancy > epsilon) {
    // FBBT Up/Down iteration
    for (int h = 0; h < m; h++) {
      e[h]->FBBTUpDown(XL, XU, elb[h], eub[h]);
    }
    // compute discrepancy
    discrepancy = 0;
    for(mi = XL.begin(); mi != XL.end(); mi++) {
      index = mi->first;
      discrepancy += fabs(mi->second - ZL[index]);
    } 
    for(mi = XU.begin(); mi != XU.end(); mi++) {
      index = mi->first;
      discrepancy += fabs(ZU[index] - mi->second);
    } 
    // print out current intervals
    if (discrepancy > epsilon) {
      cout << "# FBBT iteration " << iteration << ": discrepancy = " 
	   << discrepancy << endl;
      for(index = 1; index <= n; index++) {
	cout << "#  x" << index << " in [" << XL[index] << "," << XU[index] 
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

  cout << "###################################################" << endl;


  return 0;
}
