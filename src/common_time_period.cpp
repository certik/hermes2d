#include <time.h>
#include "common.h"
#include "common_time_period.h"

using namespace std;

#ifndef WIN32
# define H2D_TIME_SEC_FRAC 1000
#endif

TimePeriod::TimePeriod(const char *name) : period_name(name == NULL ? "unnamed" : name) {
  //initialization
#ifdef WIN32 //Windows
  LARGE_INTEGER freq;
  if (QueryPerformanceFrequency(&freq))
    frequency = (double)freq.QuadPart;
  else
    frequency = -1;
#endif
  tick_reset();
}

uint64_t TimePeriod::get_time() const {
#ifdef WIN32 //Windows
  if (frequency > 0) {
    LARGE_INTEGER ticks;
    QueryPerformanceCounter(&ticks);
    return ticks.QuadPart;
  }
  else {
    return clock();
  }
#else //Linux
  timespec tm;
  clock_gettime(CLOCK_REALTIME, &tm);
  return (uint64_t)tm.tv_sec * H2D_TIME_SEC_FRAC + (tm.tv_nsec / (1000000 / H2D_TIME_SEC_FRAC));
#endif
}

double TimePeriod::to_seconds(const uint64_t& period) const {
#ifdef WIN32 //Windows
  if (frequency > 0)
    return period / frequency;
  else
    return period / (double)CLOCKS_PER_SEC;
#else //Linux
  return period / (double)H2D_TIME_SEC_FRAC;
#endif
}

const TimePeriod& TimePeriod::tick(TimerPeriodTickType type) {
  uint64_t cur_time = get_time();
  double secs = to_seconds(cur_time - last_time);
  if (type == H2D_ACCUMULATE) {
    accum += secs;
    last_period = secs;
  }
  else {
    last_period = -1;
  }
  last_time = cur_time;
  return *this;
}

const TimePeriod& TimePeriod::tick_reset() {
  tick(H2D_SKIP);
  reset();
  return *this;
}

const TimePeriod& TimePeriod::reset() {
  accum = 0;
  last_period = -1;
  return *this;
}

string TimePeriod::to_string(double secs) const {
  if (secs < 0)
    return "NO TIME";
  else {
  int hours = (int) secs / (3600);
  int mins = (int) fmod(secs, 3600) / 60;
  secs = fmod(secs, 60);

  stringstream str;
  if (hours > 0)
    str << hours << "h ";
  if (hours > 0 || mins > 0)
    str << mins << "m ";
  str << secs << "s";

  return str.str();
  }
}

HERMES2D_API ostream& operator<<(ostream& stream, const TimePeriod& period) {
  stream << period.accumulated_str();
  return stream;
}
