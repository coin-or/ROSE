/*
** Name:       newsolver.h
** Author:     Leo Liberti
** Source:     GNU C++
** Purpose:    Solver selection function (header file)
** License:    (C) Leo Liberti, all rights reserved. Code published under the 
               Common Public License.
** History:    050429 work started
*/

#ifndef NEWSOLVERH
#define NEWSOLVERH

#include <string>
#include "solver.h"

Solver* NewSolver(std::string solvername);
void DeleteSolver(Solver* s);

#endif
