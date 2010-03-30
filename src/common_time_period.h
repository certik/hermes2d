// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __HERMES2D_COMMON_TIME_PERIOD_H
#define __HERMES2D_COMMON_TIME_PERIOD_H

enum TimerPeriodTickType { ///< Tick type.
  H2D_ACCUMULATE, ///< Accumulate a period between ticks.
  H2D_SKIP ///< Skip period between ticks, i.e., do not accumulate it.
};

//TODO: Measure time that CPU spent on the task instead of a global time.
/// \brief A named time period measurement with accumulation.
/// An instance of the timer should not be used across threads. Timer is not thread-safe.
class HERMES2D_API TimePeriod {
public:
  TimePeriod(const char *name = NULL); ///< Constructs internal structures and starts measuring.

  const TimePeriod& reset(); ///< Resets accumulated time.
  const TimePeriod& tick_reset(); ///< Starts a new period and resets accumulated time.
  const TimePeriod& tick(TimerPeriodTickType type = H2D_ACCUMULATE); ///< Starts/ends a new period.

  const std::string& name() const { return period_name; }

  double accumulated() const { return accum; }; ///< Returns accumulated time (in seconds).
  std::string accumulated_str() const { return to_string(accum); }; ///< Returns accumulated time in human readable form.
  double last() const { return last_period; }; ///< Returns last measured period (in seconds). -1 if period was skipped.
  std::string last_str() const { return to_string(last_period); }; ///< Returns last measured period in human readable form.

private:
  const std::string period_name; ///< Name of the timer (can be empty)
  double last_period; ///< Time of the last measured period.
  uint64_t last_time; ///< Time when the timer was started/resumed (in platform-dependent units).
  double accum; ///< Time accumulator (in seconds).
#ifdef WIN32
  double frequency; ///< Frequency of the performance timer. If zero, no hi-res timer is supported. (Win32 only)
#endif

  uint64_t get_time() const; ///< Returns current time (in platform-dependent units).
  double to_seconds(const uint64_t& period) const; ///< Coverts time from platform-dependent units to seconds.
  std::string to_string(const double time) const; ///< Converts time from seconds to human readable form.

  friend HERMES2D_API std::ofstream& operator<<(std::ofstream& stream, const TimePeriod& period);
};
extern HERMES2D_API std::ostream& operator<<(std::ostream& stream, const TimePeriod& period);

#endif
