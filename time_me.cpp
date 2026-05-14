#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include <functional>

using namespace std;
using namespace chrono;

template<typename F>
double time_me(F&& f) {
    auto start = high_resolution_clock::now();
    f();
    auto end = high_resolution_clock::now();
    return duration_cast<microseconds>(end - start).count() / 1000.0;
}

#endif
