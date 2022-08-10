/*
 * timing.h
 *
 *  Created on: 31/lug/2010
 *      Author: Luca
 */

#ifndef DEBORAH_TIMING_H_
#define DEBORAH_TIMING_H_

#include <time.h>

#include <chrono>

#define READTSC(buf) __asm__ __volatile__ ( \
  "rdtsc\n\t" \
  "movl %%eax,0(%0)\n\t" \
  "movl %%edx,4(%0)" \
        : : "r" (buf) : "eax", "edx" )

// Coarse-precision timer based on clock() API
class BasicTimer
{
 public:
  BasicTimer()
  {
      t0 = t1 = clock();
  };

  double time()
  {
      return (double) (t1 - t0) / (double) CLOCKS_PER_SEC;
  };

  double mark(bool t)
  {
      if (t)
      {
          t0 = clock();
      }
      else
      {
          t1 = clock();
      }
      return time();
  };
 private:
  clock_t t0, t1;
};

// High-precision timer based on std::chrono
class NanoTimer
{
 public:
  NanoTimer()
  {
      t0 = t1 = std::chrono::high_resolution_clock::now();
  }

  double time()
  {
      std::chrono::duration<float> diff = t1 - t0;
      return diff.count();
  }

  double mark(bool reset)
  {
      if (reset)
      {
          t0 = std::chrono::high_resolution_clock::now();
      }
      else
      {
          t1 = std::chrono::high_resolution_clock::now();
      }

      return time();
  }

 private:
  std::chrono::time_point<std::chrono::high_resolution_clock> t0, t1;
};

#endif /* TIMING_H_ */
