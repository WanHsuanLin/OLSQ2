/***********************************************************************

  File        [ timeUsage.cpp ]

  System      [ CRoute: Academic Router with Cell Movement ]

  Package     [ misc ]

  Synopsis    [ TimeUsage class implementation ]

  Author      [  ]

  Affiliation [ National Taiwan University ]

  Date        [ 7, Jun., 2021 ]

***********************************************************************/

#include <sys/time.h>
#include <sys/resource.h>

#include "timeUsage.hpp"

OLSQ_NAMESPACE_CPP_START

TimeUsage::TimeUsage() {
  start(FULL);
  start(PARTIAL);
}

void TimeUsage::start(TimeType type) {
  (type == FULL) ? checkUsage(_fullStart) : checkUsage(_partialStart);
}

void TimeUsage::showUsage(const char* comment, TimeType type) {
  TimeState curSt;
  checkUsage(curSt);
  TimeState dur = (type == FULL) ? diff(_fullStart, curSt) : diff(_partialStart, curSt);
  if (type == FULL) {
    fprintf(stderr, "---------- %-20s total  time usage -----------\n", comment);
  }
  else {
    fprintf(stderr, "---------- %-20s period time usage -----------\n", comment);
  }
  fprintf(stderr, "Real: %fs; User: %fs; System: %fs\n\n", dur.realTime, dur.userTime, dur.sysTime);
}

TimeState TimeUsage::diff(TimeState& start, TimeState& end) {
  return TimeState(end.realTime - start.realTime,
                    end.userTime - start.userTime,
                    end.sysTime  - start.sysTime);
}

void TimeUsage::checkUsage(TimeState& st) const {
  rusage tUsg;
  getrusage(RUSAGE_SELF, &tUsg);
  timeval tReal;
  gettimeofday(&tReal, NULL);
  st.userTime = tUsg.ru_utime.tv_sec + tUsg.ru_utime.tv_usec / TIME_SCALE;
  st.sysTime  = tUsg.ru_stime.tv_sec + tUsg.ru_stime.tv_usec / TIME_SCALE;
  st.realTime = tReal.tv_sec + tReal.tv_usec / TIME_SCALE;
}

double TimeUsage::fullRealTime() const {
  timeval tReal;
  gettimeofday(&tReal, NULL);
  double curReal = tReal.tv_sec + tReal.tv_usec / TIME_SCALE;
  return curReal - _fullStart.realTime;
}

OLSQ_NAMESPACE_CPP_END