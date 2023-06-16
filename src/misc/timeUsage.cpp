/***********************************************************************

  File        [ timeUsage.cpp ]

  System      [ CRoute: Academic OLSQ2 with Cell Movement ]

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

TimeState diff(TimeState& start, TimeState& end) {
  return TimeState(end.realTime - start.realTime,
                    end.userTime - start.userTime,
                    end.sysTime  - start.sysTime);
}

void TimeState::checkUsage() {
  rusage tUsg;
  getrusage(RUSAGE_SELF, &tUsg);
  timeval tReal;
  gettimeofday(&tReal, NULL);
  userTime = tUsg.ru_utime.tv_sec + tUsg.ru_utime.tv_usec / TIME_SCALE;
  sysTime  = tUsg.ru_stime.tv_sec + tUsg.ru_stime.tv_usec / TIME_SCALE;
  realTime = tReal.tv_sec + tReal.tv_usec / TIME_SCALE;
}

TimeUsage::TimeUsage(unsigned_t timeout) {
  _timeout = timeout;
  start(FULL);
  start(PARTIAL);
}

void TimeUsage::start(TimeType type) {
  (type == FULL) ? _fullStart.checkUsage() : _partialStart.checkUsage();
}

void TimeUsage::showUsage(const char* comment, TimeType type) {
  TimeState curSt;
  curSt.checkUsage();
  TimeState dur = (type == FULL) ? diff(_fullStart, curSt) : diff(_partialStart, curSt);
  if (type == FULL) {
    fprintf(stdout, "---------- %-20s total  time usage -----------\n", comment);
  }
  else {
    fprintf(stdout, "---------- %-20s period time usage -----------\n", comment);
  }
  fprintf(stdout, "Real: %fs; User: %fs; System: %fs\n\n", dur.realTime, dur.userTime, dur.sysTime);
}

double_t TimeUsage::fullRealTime() const {
  timeval tReal;
  gettimeofday(&tReal, NULL);
  double_t curReal = tReal.tv_sec + tReal.tv_usec / TIME_SCALE;
  return curReal - _fullStart.realTime;
}


OLSQ_NAMESPACE_CPP_END