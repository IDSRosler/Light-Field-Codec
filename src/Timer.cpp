#include "Timer.h"

Timer::Timer() {
    this->n_calls = 0;
    std::chrono::duration<double>::zero();
}

Timer::Timer(const Timer &orig) {
}

Timer::~Timer() = default;

void Timer::tic() {
    this->start = std::chrono::system_clock::now();
    ++this->n_calls;
}

void Timer::toc() {
    this->end = std::chrono::system_clock::now();
    this->duration += this->end - this->start;
}

double Timer::getTotalTime() {
    return this->duration.count();
}

long Timer::get_nCalls() {
    return this->n_calls;
}