/**********************************************************************
* Author:      Leo Liberti                                            *
* Name:        exceptions.h                                           *
* Source:      GNU C++                                                *
* Purpose:     Expression v3 exceptions                               *
* History:     010622 0.0 work started                                *
* License:     (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
***********************************************************************/

#ifndef __EV3EXCEPTIONSH__
#define __EV3EXCEPTIONSH__

#define RCS8 "$Id: exceptions.h,v 1.5 2003/10/08 11:13:58 liberti Exp liberti $"

#define HELPURL "http://liberti.dhs.org/"
#define NONE "[none]"
#define STDACTION std::cerr << interface << "::" << scope << ": in [" << operation << "]: " << description << ", code = " << code << ", see " << moreinfo << std::endl

class ErrBase {
 public:
  unsigned long code;
  std::string interface;
  std::string scope;
  std::string operation;
  std::string description;
  std::string moreinfo;
  ErrBase() :
    code(0), interface(NONE), scope(NONE), operation(NONE), 
    description(NONE),  moreinfo(HELPURL) { STDACTION; }
  ErrBase(unsigned long mycode,
	  std::string myif,
	  std::string myscope,
	  std::string myop,
	  std::string mydesc,
	  std::string myinfo) :
    code(mycode), interface(myif), scope(myscope), operation(myop), 
    description(mydesc), moreinfo(myinfo) { STDACTION; }
};

class ErrUnknown : public ErrBase { 
 public:
  ErrUnknown(unsigned long mycode,
	     std::string myif,
	     std::string myscope,
	     std::string myop,
	     std::string mydesc,
	     std::string myinfo) : 
    ErrBase(mycode, myif, myscope, myop, mydesc, myinfo) { STDACTION; }
};

class ErrNotPermitted : public ErrBase { 
 public:
  ErrNotPermitted(unsigned long mycode,
		  std::string myif,
		  std::string myscope,
		  std::string myop,
		  std::string mydesc,
		  std::string myinfo) :
    ErrBase(mycode, myif, myscope, myop, mydesc, myinfo) { STDACTION; }
};

class ErrDivideByZero : public ErrBase { 
 public:
  std::string dividend;
  ErrDivideByZero(unsigned long mycode,
		  std::string myif,
		  std::string myscope,
		  std::string myop,
		  std::string mydesc,
		  std::string myinfo, 
		  std::string mydiv) : 
    ErrBase(mycode, myif, myscope, myop, mydesc, myinfo),
    dividend(mydiv) { STDACTION; }
};

#endif
