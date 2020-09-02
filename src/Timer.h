#ifndef LF_CODEC_TIMER_H
#define LF_CODEC_TIMER_H

#include <ctime>
#include <chrono>
#include <iostream>

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
