#include <sys/time.h>

/**
 * @brief The Timer class can be used as a stop watch to measure running times.
 * The resolution is milliseconds.
 */
class Timer {
 public:
  Timer() {
    mStart = 0;
    mEnd = 0;
  }

  /**
   * @brief Resets and starts the timer.
   */
  void start() {
    mStart = now();
    mEnd = 0;
  }

  /**
   * @brief Stops the timer.
   */
  void stop() { mEnd = now(); }

  /**
   * @brief Stops the timer if it is still running and returns the elapsed time.
   * @return Elapsed time in milliseconds.
   */
  uint64_t elapsed() {
    if (mEnd == 0) {
      stop();
    }
    return mEnd - mStart;
  }

 private:
  /**
   * @brief Helper function to retrieve a timestamp.
   * @return Number of milliseconds elapsed since 1970-01-01 00:00:00 +0000
   * (UTC).
   */
  uint64_t now() {
    timeval tv;
    gettimeofday(&tv, NULL);
    return 1000 * tv.tv_sec + tv.tv_usec / 1000;
  }

  uint64_t mStart;
  uint64_t mEnd;
};
