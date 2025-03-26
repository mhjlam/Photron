#pragma once

#include <chrono>
#include <format>
#include <string>


class Timer
{
public:
    Timer() : rt0(std::chrono::system_clock::now()) {}

    void Reset()
    {
        rt0 = std::chrono::system_clock::now();
    }

    long long Punch()
    {
        auto now = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed = now - rt0;
        return std::chrono::duration_cast<std::chrono::seconds>(elapsed).count();
    }

    std::string HoursMinutesSeconds(long long seconds)
    {
        std::chrono::duration<long long> elapsed(seconds);
        return std::format("{:%H:%M:%S}", elapsed);
    }

private:
    std::chrono::system_clock::time_point rt0;
};
