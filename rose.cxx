/*
** Name:       rose.cxx
** Author:     Leo Liberti
** Source:     GNU C++
** Purpose:    main for ROSE
** License:    (C) Leo Liberti, all rights reserved. Code published under the
               Common Public License.
** History:    050215 work started
*/

#include "rose.h"
#include "version.h"
#include "solverlist.h"
#include <sys/time.h>
#include <cstring>
#include <unistd.h>

// AMPL stuff
#undef filename
#undef amplflag
#undef write_sol
#undef basename

#include <libgen.h>

#define _STRINGIFY(x) #x
#define STRINGIFY(x) _STRINGIFY(x)


// to add an AMPL option, also change problem.cxx, problem.h
// NOTE: keep all parameter names in ALPHABETICAL ORDER in all lists!!!
// this is important for the way AMPL works - it will fail otherwise
namespace Ampl {
  static int convexifieroutampl = 0;
  static int disaggrk = 5;
  static int glpksolutionmethod = 0;
  static double limitedbrancharbitrary = 0.2;
  static int limitedbranchlocalsolver = 0;
  static int limitedbranchmaxtime = 0;
  static int localbranchkmax = 10;
  static int localbranchwarmup = 1;
  static int mainsolver = 3;
  static int nosimplifier = 0;
  static int relaxinteger = 0;
  static int quiet = 0;
  static int rquarticconvextype = 0;
  static int snopt6majoriterations = 1000;
  static int snopt6minoriterations = 1000;
  static double vnsepsilon = 1e-4;
  static int vnskmax = 10;
  static int vnskmin = 0;
  static int vnslocalsolver = -1;
  static int vnsmaxtime = 0;
  static int vnsneighbourhood = 1;
  static int vnssamples = 1;
  static int vnswarmup = 1;
  static int symmgroupouttype = 0;
  static keyword keywds[] = { // must be alphabetical
    KW(CHR"convexifieroutampl", I_val, &convexifieroutampl,
       CHR"Whether Rconvexifier should produce AMPL output"),
    KW(CHR"disaggrk", I_val, &disaggrk,
       CHR"Disaggr approximation (k->oo => opt-reformulation)"),
    KW(CHR"glpksolutionmethod", I_val, &glpksolutionmethod,
       CHR"GLPK solution method (0=simplex, 1=interior)"),
    KW(CHR"limitedbrancharbitrary", D_val, &limitedbrancharbitrary,
       CHR"LimitedBranch arbitrary choice probability (in [0,1])"),
    KW(CHR"limitedbranchlocalsolver", I_val, &limitedbranchlocalsolver,
       CHR"LimitedBranch local solver"),
    KW(CHR"limitedbranchmaxtime", I_val, &limitedbranchmaxtime,
       CHR"LimitedBranch max running time (secs)"),
    KW(CHR"localbranchkmax", I_val, &localbranchkmax,
       CHR"LocalBranch max neighbourhood size"),
    KW(CHR"localbranchwarmup", I_val, &localbranchwarmup,
       CHR"LocalBranch local solver warm up"),
    KW(CHR"mainsolver", I_val, &mainsolver, CHR"set the Rose solver"),
    KW(CHR"nosimplifier", I_val, &nosimplifier,
       CHR"no automatic symbolic simplification"),
    KW(CHR"quiet", I_val, &quiet, CHR"suppress messages"),
    KW(CHR"relaxinteger", I_val, &relaxinteger, CHR"relax integer variables"),
    KW(CHR"rquarticconvextype", I_val, &rquarticconvextype,
       CHR"type of quartic convexification"),
    KW(CHR"snopt6majoriterations", I_val, &snopt6majoriterations,
       CHR"Snopt6 major iterations limit"),
    KW(CHR"snopt6minoriterations", I_val, &snopt6minoriterations,
       CHR"Snopt6 minor iterations limit"),
    KW(CHR"symmgroupouttype", I_val, &symmgroupouttype,
       CHR"output type for symmgroup reformulator"),
    KW(CHR"vnsepsilon", D_val, &vnsepsilon, CHR"VNS epsilon tolerance"),
    KW(CHR"vnskmax", I_val, &vnskmax, CHR"VNS max neighbourhood size"),
    KW(CHR"vnskmin", I_val, &vnskmax, CHR"VNS min neighbourhood size"),
    KW(CHR"vnslocalsolver", I_val, &vnslocalsolver, CHR"VNS local solver"),
    KW(CHR"vnsmaxtime", I_val, &vnsmaxtime, CHR"VNS max running time (secs)"),
    KW(CHR"vnsneighbourhood", I_val, &vnsneighbourhood,
       CHR"VNS neighbourhood (0=shells, 1=filled rectangles)"),
    KW(CHR"vnssamples", I_val, &vnssamples, CHR"VNS number of samples in nb."),
    KW(CHR"vnswarmup", I_val, &vnswarmup, CHR"VNS local solver warm up"),
  };
  // this holds the AMPL structure
  ASL* asl;
  // these data structures used to return information to AMPL
  static Option_Info Oinfo = { CHR"rose",
			       CHR"Rose -- Leo Liberti (c) 2005-2009",
			       CHR"roseamp_options", keywds, nkeywds, 0,
			       CHR"Rose for AMPL"};
  Option_Info* TheOInfo = &Oinfo;
  typedef struct { char *msg; int code, wantsol; } Sol_info;
  Sol_info *Si;
  static Sol_info solinfo[] = {
    { /* 0 */ "Optimal solution found", 0, 1 },
    { /* 1 */ "Infeasible problem" , 200, 1 },
    { /* 2 */ "Constraint infeasible" , 201, 1 },
    { /* 3 */ "Integer infeasible", 202, 1}
  };
  SufDecl suftab[] = {
    { "integer", 0, ASL_Sufkind_var | ASL_Sufkind_output },
    { "linear", 0, ASL_Sufkind_con | ASL_Sufkind_output },
    { "feasible", 0, ASL_Sufkind_con | ASL_Sufkind_output },
    { "evconvex", 0, ASL_Sufkind_con | ASL_Sufkind_output },
    { "evconcave", 0, ASL_Sufkind_con | ASL_Sufkind_output },
    { "amplsolverindex", 0, ASL_Sufkind_var | ASL_Sufkind_output },
    { "linear", 0, ASL_Sufkind_obj | ASL_Sufkind_output },
    { "infeasquant", 0, ASL_Sufkind_con | ASL_Sufkind_real | ASL_Sufkind_output}
  };

  int* TheConvexifierOutAmpl = &convexifieroutampl;
  int* TheDisaggrK = &disaggrk;
  int* TheGLPKSolutionMethod = &glpksolutionmethod;
  double* TheLimitedBranchArbitrary = &limitedbrancharbitrary;
  int* TheLimitedBranchLocalSolver = &limitedbranchlocalsolver;
  int* TheLimitedBranchMaxTime = &limitedbranchmaxtime;
  int* TheLocalBranchKmax = &localbranchkmax;
  int* TheLocalBranchWarmUp = &localbranchwarmup;
  int* TheNoSimplifier = &nosimplifier;
  int* TheMainSolver = &mainsolver;
  int* TheQuiet = &quiet;
  int* TheRelaxInteger = &relaxinteger;
  int* TheRQuarticConvexType = &rquarticconvextype;
  int* TheSnopt6MajorIterations = &snopt6majoriterations;
  int* TheSnopt6MinorIterations = &snopt6minoriterations;
  double* TheVNSEpsilon = &vnsepsilon;
  int* TheSymmgroupOutType = &symmgroupouttype;
  int* TheVNSKmax = &vnskmax;
  int* TheVNSKmin = &vnskmin;
  int* TheVNSLocalSolver = &vnslocalsolver;
  int* TheVNSMaxTime = &vnsmaxtime;
  int* TheVNSNeighbourhood = &vnsneighbourhood;
  int* TheVNSSamples = &vnssamples;
  int* TheVNSWarmUp = &vnswarmup;
};

#define MAXBUFSIZE 8192

#ifdef PROFILER
extern "C" int main(int argc, char** argv);
extern "C" void etext();
extern "C" int monstartup(int (*)(int, char**), void (*)());
//extern "C" void monitor(void *first, void *last = 0, short *buf = 0,
// int bufsize = 0, int nfuncs = 0);
extern "C" int moncontrol(int mode);
#endif

void Help(bool amplflag) {
  using namespace std;
  cout << "ROSE (Reformulation/Optimization Software Engine) - (c) L. Liberti 2005-2010" << endl
  #ifdef ROSEVERSION
       <<"     version " << STRINGIFY(ROSEVERSION)
  #endif
  #ifdef ROSEBUILDDATE
       <<"     build " << ROSEBUILDDATE
  #endif
         << endl;
  if (amplflag) {
    cout << "Syntax: roseamp is an AMPL (www.ampl.com) solver" << endl;
    cout << "  insert \"option solver roseamp\" in your AMPL directives" << endl;
    cout << "  after making sure that the executable \"roseamp\" "
	 << "is in the path" << endl;
    cout << "AMPL solver options (use \"option roseamp_options\" "
	 << "AMPL directive)" << endl;
    cout << "  convexifieroutampl (0=no, 1=yes)" << endl;
    //    cout << "  glpksolutionmethod (0=simplex, 1=interior)" << endl;
    cout << "  limitedbrancharbitrary (double in [0,1])" << endl;
    cout << "  limitedbranchlocalsolver (solver code)" << endl;
    cout << "  limitedbranchmaxtime (in seconds)" << endl;
    cout << "  localbranchwarmup (0 or 1)" << endl;
    cout << "  mainsolver (solver code)" << endl;
    cout << "  nosimplifier (0 or 1)" << endl;
    cout << "  quiet (0 or 1)" << endl;
    cout << "  relaxinteger (0 or 1)" << endl;
    //    cout << "  snopt6majoriterations (number of iterations)" << endl;
    //    cout << "  snopt6minoriterations (number of iterations)" << endl;
    cout << "  symmgroupouttype (0=nauty[default], 1=AMPL .dat [fpp.mod], 2=AMPL .dat [MILP])" << endl;
    cout << "  vnsepsilon (double)" << endl;
    cout << "  vnskmax (max VNS neighbourhood size (int))" << endl;
    cout << "  vnskmin (min VNS neighbourhood size (int))" << endl;
    cout << "  vnslocalsolver (solver code)" << endl;
    cout << "  vnsmaxtime (seconds)" << endl;
    cout << "  vnssamples (number of starting points in each neighbourhoods)"
	 << endl;
    cout << "  vnswarmup (0 or 1)" << endl;
    cout << "where the solver code is an integer defined as follows:" << endl;
    // cout << "   0: snopt" << endl;
    cout << "   1: tabu(letsgo) [don't use]" << endl;
    // cout << "   2: glpk" << endl;
    cout << "   3: vns [don't use]" << endl;
    cout << "   4: gomory [don't use]" << endl;
    cout << "   5: limitedbranch [don't use]" << endl;
    cout << "   6: localbranch [don't use]" << endl;
    cout << "   7: analyser [outputs problem info to AMPL suffixes]" << endl;
    cout << "   8: prodbincont [reformulates xy where x binary]" << endl;
    cout << "   9: print [prints out the formulation]" << endl;
    cout << "  10: smith [Smith standard form]" << endl;
    cout << "  11: copy [internal use]" << endl;
    cout << "  12: disaggr [disaggregates vars into sums thereof]" << endl;
    cout << "  13: oa [outer approximation (unfinished)]" << endl;
    cout << "  14: symmgroup [output coloured DAG representation of formulation]" << endl;
    cout << "  15: printmod [output AMPL flat formulation]" << endl;
    // cout << "  16: ipopt" << endl;
    cout << "  17: rvinci [output Vinci input form]" << endl;
    cout << "  18: rporta [output Porta input form]" << endl;
    cout << "  19: rquarticconvex [output cvx rel of quartic formulations]" << endl;
    cout << "  20: rconvexifier [output cvx relaxation from Smith std form]" << endl;
    cout << "  21: rprintdat [output LP in .mod/.dat form]" << endl;
    cout << "  22: rcdd [output cdd input form]" << endl;
    cout << "  23: rrelaxation [unfinished (empty)]" << endl;
    cout << "  24: subgradient [unfinished (empty)]" << endl;
    cout << "  25: rconvexifiermod [output cvx rel of polynomial programs]" << endl;
    cout << "  26: rprintncvxdiscrepancy [discrepancy between defvars and their cvxrel]" << endl;
    cout << "  27: rfbbtfp [perform FBBT and output FBBTLP formulation]" << endl;
    cout << "  28: rmilp2gph [output bipartite constr/var incidence graph as .gph]" << endl;
  } else {
    cout << "Syntax: rose [options] file\n";
    cout << "  where file specifies the optimization problems (.ros format),\n";
    cout << "  and options can be:\n";
    cout << "    -h for this help\n";
    cout << "    -p parameter_list (\"param val, param val, ...\", quote the list)\n";
    cout << "       (param names are capitalized, e.g. SymmgroupOutType)\n";
    cout << "    -n no automatic simplifier\n";
    cout << "    -q for quiet mode\n";
    cout << "    -s solvername (" << ROSESOLVERLIST << ")\n";
    cout << "    -t table output\n";
    cout << "Rose may also be used as an AMPL solver -- "
	 << "rename it as \"roseamp\" and\n";
    cout << "put it in the path, then call ampl with "
	 << "\"option solver roseamp\"\n";
  }
}

int main(int argc, char** argv) {

  using namespace std;

#ifdef PROFILER
  monstartup(main, etext);
#endif

  // initialize randomizer
  struct timeval theTV;
  struct timezone theTZ;
  gettimeofday(&theTV, &theTZ);
  srandom(theTV.tv_usec);

  // see if we were called from AMPL
  bool amplflag = false;
  char* tfn;
  if ((tfn = strrchr(argv[0], '/')) == NULL) {
    tfn = argv[0];
  } else {
    tfn++;
  }
  if (strncmp(tfn, "roseamp", 7) == 0 || strncmp(tfn, "moronamp", 8) == 0) {
    amplflag = true;
  }

  // create the problem
  Problem *p;
  int ret = 0;
  if (argc < 2) {
    Help(amplflag);
    exit(2);
  }

  bool Quiet = false;

  // parse the command line
  int c;
  char parms[MAXBUFSIZE];
  char parmname[MAXBUFSIZE];
  char parmvalue[MAXBUFSIZE];
  char* t1 = parms;
  bool boolval;
  int intval;
  double doubleval;
  string stringval;
  string pn;
  int paramtype;
  char filename[MAXBUFSIZE];
  string solvername = "none";
  bool TableOut = false;
  bool nosimp = false;
  bool output = false;

  // see if filename is a stub .nl file
  for(int i = 1; i < argc; i++) {
    if (strstr(argv[i], ".nl") != NULL) {
      amplflag = true;
      strncpy(filename, argv[i], MAXBUFSIZE);
      break;
    }
    if (strstr(argv[i], "-AMPL") != NULL) {
      amplflag = true;
      strncpy(filename, "FakeProblemName", MAXBUFSIZE);
      break;
    }
  }

  if (amplflag) {
    p = new Problem(nosimp);
    Ampl::asl = p->ParseAMPL(argv, argc);
    suf_declare_ASL(Ampl::asl, Ampl::suftab,
		    sizeof(Ampl::suftab)/sizeof(Ampl::SufDecl));
  } else {
    while(1) {
      c = getopt(argc, argv, "nhop:qs:t");
      switch(c) {
      case 'n':
	nosimp = true;
	break;
      case 'o':
	output = true;
	break;
      case 'p':
	strncpy(parms, optarg, MAXBUFSIZE);
	break;
      case 'q':
	Quiet = true;
	break;
      case 'h':
	Help(amplflag);
	exit(2);
      case 's':
	solvername = optarg;
	break;
      case 't':
	TableOut = true;
	Quiet = true;
	break;
      default:
	break;
      }
      if (optind < argc) {
	strncpy(filename, argv[optind], MAXBUFSIZE);
      } else {
	cerr << "rose: error: need filename on command line\n";
	exit(1);
      }
      if (c == -1) {
	break;
      }
    }
    // parse the problem
    p = new Problem(nosimp);
    p->Parse(&filename[0]);
  }

  // get problem name
  char* thebasec = strdup(filename);
  char* thebasename = basename(thebasec);
  char* thebt = strrchr(thebasename, '.');
  if (thebt != NULL) {
    *thebt = '\0';
  }
  string ProbName = thebasename;
  free(thebasec);
  p->SetName(ProbName);

  // get parameters from problem
  Parameters& blob = p->GetParamsRef();

  // merge in parameters from command line
  if (strlen(parms) > 0) {
    while(1) {
      // skip spaces
      while(*t1 == ' ') {
	t1++;
      }
      // find parameter name
      c = 0;
      while(*t1 != ' ') {
	parmname[c] = *t1;
	c++;
	t1++;
      }
      parmname[c] = '\0';
      pn = parmname;
      // skip spaces
      while(*t1 == ' ') {
	t1++;
      }
      // find parameter value
      c = 0;
      while(*t1 != ',' && *t1 != ' ' && *t1 != '\0') {
	parmvalue[c] = *t1;
	c++;
	t1++;
      }
      parmvalue[c] = '\0';
      if (pn.empty() || strlen(parmvalue) == 0) {
	cerr << "rose: parsing error in cmd line param " << parmname << endl;
	exit(3);
      }
      // found valid parameter, parse it
      paramtype = parmtype(parmvalue, boolval, intval, doubleval, stringval);
      if (paramtype == BoolType) {
	blob.SetBoolParameter(pn, boolval);
      } else if (paramtype == IntType) {
	blob.SetIntParameter(pn, intval);
      } else if (paramtype == DoubleType) {
	blob.SetDoubleParameter(pn, doubleval);
      } else if (paramtype == StringType) {
	blob.SetStringParameter(pn, stringval);
      }
      if (*t1 == ',') {
	t1++;
      }
      while(*t1 == ' ') {
	t1++;
      }
      if (*t1 == '\0') {
	break;
      }
    }
  }
  if (solvername != "none") {
    blob.SetStringParameter("MainSolver", solvername);
  }

  if (Quiet) {
    blob.SetBoolParameter("Quiet", true);
  }

  if (amplflag) {
    blob.SetBoolParameter("AmplSolver", true);
  } else {
    blob.SetBoolParameter("AmplSolver", false);
  }

  // get number of parameters
  int np = blob.GetNumberOfParameters();

  // find out solver
  if (solvername == "none") {
    solvername = "print";
  }
  for(int i = 1; i <= np; i++) {
    if (blob.GetParameterName(i) == "Solver" ||
	blob.GetParameterName(i) == "MainSolver" ||
	blob.GetParameterName(i) == "solver" ||
	blob.GetParameterName(i) == "mainsolver") {
      solvername = blob.GetParameterStringValue(i);
      blob.SetStringParameter("MainSolver", solvername);
    } else if (blob.GetParameterName(i) == "Quiet") {
      Quiet = blob.GetParameterBoolValue(i);
    }
  }
  if(!amplflag&&solvername=="analyser") {
     cout << "This reformulator can only be called in the AMPL version" << endl;
     exit(0);
  }

  Solver* s = NewSolver(solvername);
  s->SetProblem(p);

  if (output && !Quiet) {
    cout << *p;
  }

  if (!Quiet) {
    // display parameters
    cout << "Rose: called with parameters:\n";
    int pt;
    for(int i = 1; i <= np; i++) {
      cout << " param " << i << ": ";
      cout << blob.GetParameterName(i) << " = ";
      pt = blob.GetParameterType(i);
      if (pt == BoolType) {
	cout << blob.GetParameterBoolValue(i) << " (bool)" << endl;
      } else if (pt == IntType) {
	cout << blob.GetParameterIntValue(i) << " (int)" << endl;
      } else if (pt == DoubleType) {
	cout << blob.GetParameterDoubleValue(i) << " (double)" << endl;
      } else if (pt == StringType) {
	cout << blob.GetParameterStringValue(i) << " (string)" << endl;
      }
    }
  }

  // solve the problem
  double tstart[2], tend[2];
  double starttime, stoptime, usertimelap, totaltimelap;
  starttime = getcputime(tstart);
  if (!Quiet) {
    cout << "***** Solver (" << solvername << ") messages *****\n";
  }
  ret = s->Solve();
  if (!Quiet) {
    cout << "***** End Solver messages *****\n";
  }
  stoptime = getcputime(tend);
  usertimelap = tend[0] - tstart[0];
  totaltimelap = stoptime - starttime;

  // display the solution
  int nvar;
  double* x;
  double fval;
  map<int,double> objval;
  map<int,double> solution;
  s->GetSolution(objval, solution);
  fval = objval.begin()->second;
  map<int,double>::iterator mi;
  map<int,double>::iterator miend;
  miend = solution.end();

  if (amplflag) {
    // called from AMPL
    Ampl::Si = Ampl::solinfo + ret;
    Ampl::asl->p.solve_code_ = Ampl::Si->code;
    nvar = solution.size();
    x = new double [nvar];
    int i = 0;
    for(mi = solution.begin(); mi != miend; mi++) {
      x[i] = mi->second;
      i++;
    }
    char buf[512];
    Ampl::Sprintf(buf, "f* = %f", fval);
    Ampl::write_sol_ASL(Ampl::asl, buf, x, 0, Ampl::TheOInfo);
    delete [] x;

  } else {
    // called from cmd line
    if (TableOut) {
      // table output
      cout << "{\\tt " << ProbName << "}" << " & "
	   << p->GetNumberOfVariables() << " & "
	   << p->GetNumberOfIntegerVariables() << " & "
	   << p->GetNumberOfConstraints() << " & ";
      if (p->IsFeasible() == 1) {
	cout << fval;
      } else {
	cout << " $\\infty$ ";
      }
      cout << " & " << usertimelap << " \\\\ \\hline" << endl;
    } else {
      // normal output
      cout << "Rose: solving problem " << ProbName << " with "
	   << solvername << ":\n";
      cout << " Feasible: ";
      if (p->IsFeasible() == 1) {
	cout << "yes" << endl;
      } else {
	cout << "NO" << endl;
      }
      cout << " f* = ";
      if (s->GetOptimizationDirection() == Maximization) {
	cout << "max f = ";
      } else {
	cout << "min f = ";
      }
      cout << fval << endl;
      cout << " x* = (";
      miend--;
      for (mi = solution.begin(); mi != miend; mi++) {
	cout << mi->second << ",";
      }
      cout << miend->second << ")\n";
      cout << " return code: " << ret << endl;
      cout << " user CPU time: " << usertimelap << endl;
    }
  }

  // clean up
  DeleteSolver(s);
  //delete p;
  if (amplflag) {
    ASL_free(&Ampl::asl);
  }

#ifdef PROFILER
  moncontrol(0);
#endif

#ifdef MEMDEBUG
  cerr << "MEMCHECK: up to " << memcheckdebug.size()
       << " undeleted expression(s)" << endl;
  int mcdcount = 0;
  for(int i = 0; i < memcheckdebug.size(); i++) {
    if (memcheckdebug[i].second.second &&
	*memcheckdebug[i].second.second > 0) {
      cerr << "MEMCHECK: node = " << memcheckdebug[i].first << " ("
	   << memcheckdebug[i].second.first << "); *ncount = "
	   << *(memcheckdebug[i].second.second) << endl;
      mcdcount++;
    }
  }
  cerr << "MEMCHECK: " << mcdcount << " undeleted expression(s)" << endl;
#endif

  return 0;

}
