#ifndef LF_CODEC_TIME_H
#define LF_CODEC_TIME_H

#include <ctime>
#include <chrono>
#include <iostream>

class Time {
public:
    Time();

    Time(const Time &orig);

    virtual ~Time();

    void tic();

    void toc();

    double getTotalTime();

    long get_nCalls();

private:
    long n_calls{};
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> duration{};
};


#endif //LF_CODEC_TIME_H
