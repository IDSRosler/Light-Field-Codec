#include "Time.h"

Time::Time() {
    this->n_calls = 0;
    std::chrono::duration<double>::zero();
}

Time::Time(const Time &orig) {
}

Time::~Time() = default;

void Time::tic() {
    this->start = std::chrono::system_clock::now();
    ++this->n_calls;
}

void Time::toc() {
    this->end = std::chrono::system_clock::now();
    this->duration += this->end - this->start;
}

double Time::getTotalTime() {
    return this->duration.count();
}

long Time::get_nCalls() {
    return this->n_calls;
}