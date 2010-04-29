/*
** Name:       parameters.h
** Author:     Leo Liberti
** Source:     GNU C++
** Purpose:    Parameters class header file
** License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
** History:    050429 work started
*/

#ifndef _ROSEPARAMETERSH
#define _ROSEPARAMETERSH

#include <iostream>
#include <vector>
#include <string>

//#define NOTIMPL cerr << "method not implemented\n"; exit(11)

enum { IntType, BoolType, DoubleType, StringType };

class Parameters {

 public:
  // parameters
  // parameters: pi stands for 'paramindex' and ranges from 1
  int GetNumberOfParameters(void) const { return ParamName.size(); }
  std::string GetParameterName(int pi) const { return ParamName[pi - 1]; }
  int GetParameterType(int pi) const { return ParamType[pi - 1]; }
  std::string GetParameterStringValue(int pi) const { 
    return ParamStrVal[pi - 1]; 
  }
  double GetParameterDoubleValue(int pi) const { 
    return ParamDoubleVal[pi - 1]; 
  }
  int GetParameterIntValue(int pi) const { return ParamIntVal[pi - 1]; }
  bool GetParameterBoolValue(int pi) const { return ParamBoolVal[pi - 1]; }
  std::string GetStringParameter(std::string pname) {    
    std::string ret;
    for(int i = 0; i < (int) ParamName.size(); i++) {
      if (ParamName[i] == pname) {
	ret = ParamStrVal[i];
	break;
      }
    }
    return ret;
  }
  double GetDoubleParameter(std::string pname) {    
    double ret;
    for(int i = 0; i < (int) ParamName.size(); i++) {
      if (ParamName[i] == pname) {
	ret = ParamDoubleVal[i];
	break;
      }
    }
    return ret;
  }
  int GetIntParameter(std::string pname) {
    int ret;
    for(int i = 0; i < (int) ParamName.size(); i++) {
      if (ParamName[i] == pname) {
	ret = ParamIntVal[i];
	break;
      }
    }
    return ret;
  }
  bool GetBoolParameter(std::string pname) {
    bool ret;
    for(int i = 0; i < (int) ParamName.size(); i++) {
      if (ParamName[i] == pname) {
	ret = ParamBoolVal[i];
	break;
      }
    }
    return ret;
  }
  void SetStringParameter(std::string pn, std::string value) {
    int i;
    for(i = 0; i < (int) ParamName.size(); i++) {
      if (ParamName[i] == pn) {
        ParamType[i] = StringType;
        ParamStrVal[i] = value;
        ParamDoubleVal[i] = 0;
        ParamIntVal[i] = 0;
        ParamBoolVal[i] = false;
        return;
      } 
    }
    ParamName.push_back(pn);
    ParamType.push_back(StringType);
    ParamStrVal.push_back(value);
    ParamDoubleVal.push_back(0);
    ParamIntVal.push_back(0);
    ParamBoolVal.push_back(false);
  }
  void SetIntParameter(std::string pn, int value) {
    int i;
    for(i = 0; i < (int) ParamName.size(); i++) {
      if (ParamName[i] == pn) {
        ParamType[i] = IntType;
        ParamStrVal[i] = "";
        ParamDoubleVal[i] = 0;
        ParamIntVal[i] = value;
        ParamBoolVal[i] = false;
        return;
      } 
    }
    ParamName.push_back(pn);
    ParamType.push_back(IntType);
    ParamStrVal.push_back("");
    ParamDoubleVal.push_back(0);
    ParamIntVal.push_back(value);
    ParamBoolVal.push_back(false);
  }
  void SetBoolParameter(std::string pn, bool value) {
    int i;
    for(i = 0; i < (int) ParamName.size(); i++) {
      if (ParamName[i] == pn) {
        ParamType[i] = BoolType;
        ParamStrVal[i] = "";
        ParamDoubleVal[i] = 0;
        ParamIntVal[i] = 0;
        ParamBoolVal[i] = value;
        return;
      } 
    }
    ParamName.push_back(pn);
    ParamType.push_back(BoolType);
    ParamStrVal.push_back("");
    ParamDoubleVal.push_back(0);
    ParamIntVal.push_back(0);
    ParamBoolVal.push_back(value);
  }
  void SetDoubleParameter(std::string pn, double value) {
    int i;
    for(i = 0; i < (int) ParamName.size(); i++) {
      if (ParamName[i] == pn) {
        ParamType[i] = DoubleType;
        ParamStrVal[i] = "";
        ParamDoubleVal[i] = value;
        ParamIntVal[i] = 0;
        ParamBoolVal[i] = false;
        return;
      } 
    }
    ParamName.push_back(pn);
    ParamType.push_back(DoubleType);
    ParamStrVal.push_back("");
    ParamDoubleVal.push_back(value);
    ParamIntVal.push_back(0);
    ParamBoolVal.push_back(false);
  }


 protected:
  // parameters:
  std::vector<std::string> ParamName;
  std::vector<int> ParamType;
  std::vector<std::string> ParamStrVal;
  std::vector<double> ParamDoubleVal;
  std::vector<int> ParamIntVal;
  std::vector<bool> ParamBoolVal; 

};

#endif
