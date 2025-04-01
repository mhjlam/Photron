#pragma once

#include <chrono>
#include <format>
#include <string>


class Timer
{
public:
    Timer() : rt0(std::chrono::system_clock::now()) {}

    void reset()
    {
        rt0 = std::chrono::system_clock::now();
    }

    long long punch() const
    {
        auto now = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed = now - rt0;
        return std::chrono::duration_cast<std::chrono::seconds>(elapsed).count();
    }

    std::string hoursMinSec(long long add_seconds = 0) const
    {
        auto seconds = punch() + add_seconds;
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::seconds(seconds));
        return std::format("{:%H:%M:%S}", duration);
    }

private:
    std::chrono::system_clock::time_point rt0;
};
