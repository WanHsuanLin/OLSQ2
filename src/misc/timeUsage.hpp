/***********************************************************************

  File        [ timeUsage.hpp ]

  System      [ CRoute: Academic OLSQ2 with Cell Movement ]

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


struct TimeState {
  TimeState(double_t r = 0, double_t u = 0, double_t s = 0)
    : realTime(r), userTime(u), sysTime(s){}
  double_t realTime, userTime, sysTime;
  void checkUsage();
  void printState(const char* comment){
    fprintf(stdout, "%-20s: Real: %fs; User: %fs; System: %fs\n\n", comment, realTime, userTime, sysTime);
  }
  void printStateRatio(const char* comment, TimeState & totalTime){
    fprintf(stdout, "%-20s: Real: %fs. (%2f%%)\n", comment, realTime, (realTime/(totalTime.realTime))*100);
  }
};

TimeState diff(TimeState& start, TimeState& end);



class TimeUsage {
public:
  TimeUsage(unsigned_t timeout = 1200);

  enum TimeType { FULL, PARTIAL };

  void start(TimeType type);
  void showUsage(const char* comment, TimeType type);
  bool isTimeout(){
    TimeState curSt;
    curSt.checkUsage();
    TimeState dur = diff(_fullStart, curSt);
    return (_timeout < dur.userTime);
  }
  void setTimeout(unsigned_t s) {
    _timeout = s;
  }
  double_t fullRealTime() const;

  TimeState _fullStart, _partialStart; // total, period
  unsigned_t _timeout;
};

OLSQ_NAMESPACE_HPP_END
#endif // TIME_USAGE_HPP
