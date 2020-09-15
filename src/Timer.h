#ifndef LF_CODEC_TIMER_H
#define LF_CODEC_TIMER_H

#include <ctime>
#include <chrono>
#include <iostream>


#if LFCODEC_USE_TIMER
#define TIMER_TIC(timer) (timer).tic()
#define TIMER_TOC(timer) (timer).toc()
#else
#define TIMER_TIC(timer)
#define TIMER_TOC(timer)
#endif


class Timer {
public:
  Timer();

  Timer(const Timer &orig);

    virtual ~Timer();

    void tic();

    void toc();

    double getTotalTime();

    long get_nCalls();

private:
    long n_calls{};
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> duration{};
};


#endif // LF_CODEC_TIMER_H
