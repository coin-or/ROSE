/*
** Name:       newsolver.cxx
** Author:     Leo Liberti
** Source:     GNU C++
** Purpose:    Solver selection function
** License:    (C) Leo Liberti, all rights reserved. Code published under the
               Common Public License.
** History:    050429 work started
*/

#include <iostream>
#include <cassert>
#include "solver.h"
#include "newsolver.h"

// list of included solvers
// #include "solver_snopt6.h"
// #include "solver_ipopt.h"
// #include "solver_glpk.h"
#include "solver_tabu.h"
#include "solver_vns.h"
#include "solver_gomory.h"
#include "solver_limitedbranch.h"
#include "solver_localbranch.h"
#include "solver_analyser.h"
#include "solver_Rprodbincont.h"
#include "solver_Rprint.h"
#include "solver_Rsmith.h"
#include "solver_Rcopy.h"
#include "solver_Rdisaggr.h"
#include "solver_Roa.h"
#include "solver_Rsymmgroup.h"
#include "solver_Rprintmod.h"
#include "solver_Rvinci.h"
#include "solver_Rporta.h"
#include "solver_RQuarticConvex.h"
#include "solver_Rconvexifier.h"
#include "solver_Rconvexifiermod.h"
#include "solver_Rprintdat.h"
#include "solver_Rcdd.h"
//#include "solver_sbb.h"
#include "solver_Rrelaxation.h"
#include "solver_subgradient.h"
#include "solver_Rprintncvxdiscr.h"
#include "solver_Rfbbtfp.h"
#include "solver_null.h"
#include "solver_Rmilp2gph.h"

// add new solver header files here

Solver* NewSolver(std::string solvername) {
  Solver* ret = NULL;
  
  if (solvername == "null" || solvername == "NULL") {
    solvername = "null";
    ret = new NullSolver;
  } else if (solvername == "snopt6" || solvername == "snopt" ||
	     solvername == "SNOPT6" || solvername == "SNOPT" ||
	     solvername == "Snopt6" || solvername == "Snopt") {
    solvername = "snopt6";
#ifdef HAS_SNOPT
    ret = new Snopt6Solver;
#else
    cout << "Snopt disabled\n" << endl;
#endif 

    /*
  } else if (solvername == "glpk" || solvername == "GLPK") {
    solvername = "glpk";
    ret = new GLPKSolver;
    */

    /*
  } else if (solvername == "ipopt" || solvername == "Ipopt" ||
             solvername == "IpOpt" || solvername == "IPOPT") {
    solvername = "ipopt";
    ret = new IpoptSolver;
    */

  } else if (solvername == "Tabu" || solvername == "tabu" ||
	     solvername == "letsgo" || solvername == "LETSGO") {
    solvername = "tabu";
    ret = new TabuSolver;
  } else if (solvername == "vns" || solvername == "VNS") {
    solvername = "vns";
    ret = new VNSSolver;
  } else if (solvername == "gomory" || solvername == "Gomory") {
    solvername = "gomory";
    ret = new GomorySolver;
  } else if (solvername == "limitedbranch" || solvername == "LimitedBranch") {
    solvername = "limitedbranch";
    ret = new LimitedBranchSolver;
  } else if (solvername == "localbranch" || solvername == "LocalBranch") {
    solvername = "localbranch";
    ret = new LocalBranchSolver;
  } else if (solvername == "analyser" || solvername == "Analyser") {
    solvername = "analyser";
    ret = new AnalyserSolver;
  } else if (solvername == "prodbincont" || solvername == "Rprodbincont" ||
	     solvername == "ProdBinCont" || solvername == "RProdBinCont") {
    solvername = "prodbincont";
    ret = new ProdBinContSolver;
  } else if (solvername == "print" || solvername == "Rprint" ||
	     solvername == "Print" || solvername == "RPrint") {
    solvername = "print";
    ret = new PrintSolver;
  } else if (solvername == "smith" || solvername == "Rsmith" ||
	     solvername == "Smith" || solvername == "RSmith") {
    solvername = "smith";
    ret = new SmithSolver;
  } else if (solvername == "copy" || solvername == "Rcopy" ||
	     solvername == "Copy" || solvername == "RCopy") {
    solvername = "copy";
    ret = new CopySolver;
  } else if (solvername == "disaggr" || solvername == "Rdisaggr" ||
	     solvername == "Disaggr" || solvername == "RDisaggr") {
    solvername = "disaggr";
    ret = new DisaggrSolver;
  } else if (solvername == "oa" || solvername == "Roa" ||
	     solvername == "Oa" || solvername == "ROa") {
    solvername = "oa";
    ret = new OaSolver;
  } else if (solvername == "symmgroup" || solvername == "Rsymmgroup" ||
	     solvername == "SymmGroup" || solvername == "RSymmGroup") {
    solvername = "symmgroup";
    ret = new SymmgroupSolver;
  } else if (solvername == "printmod" || solvername == "Rprintmod" ||
	     solvername == "PrintMod" || solvername == "RPrintMod") {
    solvername = "printmod";
    ret = new PrintmodSolver;
  } else if (solvername == "rvinci" || solvername == "Vinci" ||
	     solvername == "Rvinci" || solvername == "VINCI") {
    solvername = "rvinci";
    ret = new RvinciSolver;
  } else if (solvername == "rporta" || solvername == "Porta" ||
	     solvername == "Rporta" || solvername == "PORTA") {
    solvername = "rporta";
    ret = new RportaSolver;
  } else if (solvername == "rquarticconvex" || solvername == "QuarticConvex" ||
	     solvername == "RQuarticConvex" || solvername == "QUARTICCONVEX") {
    solvername = "rquarticconvex";
    ret = new RQuarticConvexSolver;
  } else if (solvername == "rconvexifier" || solvername == "Convexifier" ||
	     solvername == "Rconvexifier" || solvername == "CONVEXIFIER") {
    solvername = "rconvexifier";
    ret = new RconvexifierSolver;
  } else if (solvername == "rconvexifiermod"||solvername == "Convexifiermod" ||
	     solvername == "Rconvexifiermod"||solvername == "CONVEXIFIERMOD") {
    solvername = "rconvexifiermod";
    ret = new RconvexifiermodSolver;
  } else if (solvername == "printdat" || solvername == "Rprintdat" ||
	     solvername == "PrintDat" || solvername == "RPrintDat") {
    solvername = "rprintdat";
    ret = new RprintdatSolver;
  } else if (solvername == "rcdd" || solvername == "Cdd" ||
	     solvername == "Rcdd" || solvername == "CDD") {
    solvername = "rcdd";
    ret = new RcddSolver;

    /*
  } else if (solvername == "sbb" || solvername == "sBB") {
    solvername = "sbb";
    ret = new SbbSolver; 
    */

  } else if (solvername == "rrelaxation" || solvername == "Relaxation" ||
	     solvername == "Rrelaxation" || solvername == "RELAXATION") {
    solvername = "rrelaxation";
    ret = new RrelaxationSolver;
  } else if (solvername == "subgradient" || solvername == "Subgradient" ||
	     solvername == "SubGradient" || solvername == "SUBGRADIENT") {
    solvername = "subgradient";
    ret = new SubgradientSolver;
  } else if (solvername == "rprintncvxdiscrepancy" || solvername == "Printncvxdiscrepancy" ||
	     solvername == "Rprintncvxdiscrepancy" || solvername == "PRINTNCVXDISCREPANCY") {
    solvername = "rprintncvxdiscrepancy";
    ret = new RprintncvxdiscrSolver;
  } else if (solvername == "fbbtfp" || solvername == "rfbbtfp" ||
	     solvername == "Rfbbtfp" || solvername == "RFBBTFP") {
    solvername = "rfbbtfp";
    ret = new RfbbtfpSolver;
  } else if (solvername == "milp2gph" || solvername == "rmilp2gph" ||
	     solvername == "Rmilp2gph" || solvername == "RMILP2GPH") {
    solvername = "rmilp2gph";
    ret = new Rmilp2gphSolver;
  }


  // add else if clauses here
  if (!ret) {
    cerr << "NewSolver(): error: could not find solver " << solvername << endl;
    assert(false);
  }
  return ret;
}

void DeleteSolver(Solver* s) {
  // can put some delete-checking mechanisms here
  if (s) {
    string solvername = s->GetName();
    if (solvername == "snopt6") {
      // LEO081115: why is the following line commented?
      //delete dynamic_cast<Snopt6Solver*>(s);

      /*
    } else if (solvername == "glpk") {
      delete dynamic_cast<GLPKSolver*>(s);
      */

      /*
    } else if (solvername == "ipopt") {
      delete dynamic_cast<IpoptSolver*>(s);
      */

    } else if (solvername == "tabu") {
      delete dynamic_cast<TabuSolver*>(s);
    } else if (solvername == "vns") {
      delete dynamic_cast<VNSSolver*>(s);
    } else if (solvername == "gomory") {
      delete dynamic_cast<GomorySolver*>(s);
    } else if (solvername == "limitedbranch") {
      delete dynamic_cast<LimitedBranchSolver*>(s);
    } else if (solvername == "localbranch") {
      delete dynamic_cast<LocalBranchSolver*>(s);
    } else if (solvername == "analyser") {
      delete dynamic_cast<AnalyserSolver*>(s);
    } else if (solvername == "prodbincont") {
      delete dynamic_cast<ProdBinContSolver*>(s);
    } else if (solvername == "print") {
      delete dynamic_cast<PrintSolver*>(s);
    } else if (solvername == "smith") {
      delete dynamic_cast<SmithSolver*>(s);
    } else if (solvername == "copy") {
      delete dynamic_cast<CopySolver*>(s);
    } else if (solvername == "disaggr") {
      delete dynamic_cast<DisaggrSolver*>(s);
    } else if (solvername == "oa") {
      delete dynamic_cast<OaSolver*>(s);
    } else if (solvername == "symmgroup") {
      delete dynamic_cast<SymmgroupSolver*>(s);
    } else if (solvername == "printmodgroup") {
      delete dynamic_cast<PrintmodSolver*>(s);
    } else if (solvername == "rvinci") {
      delete dynamic_cast<RvinciSolver*>(s);
    } else if (solvername == "rporta") {
      delete dynamic_cast<RportaSolver*>(s);
    } else if (solvername == "rquarticconvex") {
      delete dynamic_cast<RQuarticConvexSolver*>(s);
    } else if (solvername == "rconvexifier") {
      delete dynamic_cast<RconvexifierSolver*>(s);
    } else if (solvername == "rconvexifiermod") {
      delete dynamic_cast<RconvexifiermodSolver*>(s);
    } else if (solvername == "rprintdat") {
      delete dynamic_cast<RprintdatSolver*>(s);
    } else if (solvername == "rcdd") {
      delete dynamic_cast<RcddSolver*>(s);

      /*
    } else if (solvername == "sbb") {
	    delete dynamic_cast<SbbSolver*>(s); 
      */

    } else if (solvername == "rrelaxation") {
      delete dynamic_cast<RrelaxationSolver*>(s);
    } else if (solvername == "subgradient") {
      delete dynamic_cast<SubgradientSolver*>(s);
    } else if (solvername == "rprintncvxdiscrepancy") {
      delete dynamic_cast<RprintncvxdiscrSolver*>(s);
    } else if (solvername == "rmilp2gph") {
      delete dynamic_cast<Rmilp2gphSolver*>(s);
    }

  }
}
