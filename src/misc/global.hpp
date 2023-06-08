/***********************************************************************
  File        [ global.hpp ]
  System      [ OLSQ: optimal quantum layout synthesis tool]
  Package     [ misc ]
  Synopsis    [ Global header file for OLSQ project ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 22, Nov., 2022 ]
***********************************************************************/
#ifndef GLOBAL_HPP
#define GLOBAL_HPP

////////////////////////////////////////////////////////////////////////
///                          INCLUDES                                ///
////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#include <iostream>

////////////////////////////////////////////////////////////////////////
///                          NO DEBUG                                ///
////////////////////////////////////////////////////////////////////////

// #define NDEBUG

#pragma GCC diagnostic ignored "-Wstrict-overflow"
#pragma GCC diagnostic ignored "-Wunused-result"

////////////////////////////////////////////////////////////////////////
///                         NAMESPACES                               ///
////////////////////////////////////////////////////////////////////////

#define OLSQ_NAMESPACE __OLSQ__NAMESPACE__

#define OLSQ_NAMESPACE_HPP_START  namespace OLSQ_NAMESPACE { using namespace std;
#define OLSQ_NAMESPACE_HPP_END    }
#define OLSQ_NAMESPACE_CPP_START  namespace OLSQ_NAMESPACE { using namespace std;
#define OLSQ_NAMESPACE_CPP_END    }

OLSQ_NAMESPACE_HPP_START

////////////////////////////////////////////////////////////////////////
///                         BASIC TYPES                              ///
////////////////////////////////////////////////////////////////////////

constexpr unsigned  MAX_UNSIGNED =  numeric_limits<unsigned>::max()     / 3;
// constexpr int       MAX_INT      =  numeric_limits<int>::max()          / 3;
// constexpr int       MIN_INT      =  numeric_limits<int>::lowest()       / 3;
// constexpr size_t    MAX_SIZE_T   =  numeric_limits<size_t>::max()       / 3;
// constexpr double    MAX_DOUBLE   =  numeric_limits<double>::max()       / 3;
// constexpr double    MIN_DOUBLE   =  numeric_limits<double>::lowest()    / 3;
// constexpr long long MAX_LL       =  numeric_limits<long long>::max()    / 3;
// constexpr long long MIN_LL       =  numeric_limits<long long>::lowest() / 3;
// constexpr double    EPSILON      =  1e-8;

constexpr double    TIME_SCALE   = 1000000.0;
constexpr double    MEMORY_SCALE = 1024.0;

////////////////////////////////////////////////////////////////////////
///                         TYPEDEF                                  ///
////////////////////////////////////////////////////////////////////////

typedef unsigned   unsigned_t;
typedef int        int_t;
typedef double     double_t;
typedef float      float_t;

OLSQ_NAMESPACE_HPP_END

#endif // GLOBAL_HPP
