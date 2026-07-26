#pragma once

#include <chrono>

namespace ybi
{

struct Timer
{
    std::chrono::steady_clock::time_point start;
};

inline Timer StartTimer()
{
    return {std::chrono::steady_clock::now()};
}

// Returns elapsed milliseconds.
inline double StopTimer(const Timer &timer)
{
    const auto end = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(end - timer.start).count();
}

} // namespace ybi
