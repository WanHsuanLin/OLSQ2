/***********************************************************************

  File        [ timeUsage.hpp ]

  System      [ CRoute: Academic Router with Cell Movement ]

  Package     [ misc ]

  Synopsis    [ TimeUsage class header ]

  Author      [  ]

  Affiliation [ National Taiwan University ]

  Date        [ 7, Jun., 2021 ]

***********************************************************************/

#ifndef TIME_USAGE_HPP
#define TIME_USAGE_HPP

#include "global.hpp"

OLSQ_NAMESPACE_HPP_START

class TimeState {
public:
  TimeState(double r = 0, double u = 0, double s = 0)
    : realTime(r), userTime(u), sysTime(s) {}
  double realTime, userTime, sysTime;
};

class TimeUsage {
public:
  TimeUsage();

  enum TimeType { FULL, PARTIAL };

  void start(TimeType type);
  void showUsage(const char* comment, TimeType type);

  double fullRealTime() const;

private:
  TimeState diff(TimeState& start, TimeState& end);
  void checkUsage(TimeState& st) const;
  TimeState _fullStart, _partialStart; // total, period
};

OLSQ_NAMESPACE_HPP_END

#endif // TIME_USAGE_HPP