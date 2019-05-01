#include <iostream>
#include <chrono>
#include <iomanip>

int main() {
    auto currTime = std::chrono::system_clock::now();
    auto durr = currTime + std::chrono::minutes(60);
    auto newTime = currTime + std::chrono::hours(24);
    auto outCurrTime = std::chrono::system_clock::to_time_t(currTime);
    auto outNewTime = std::chrono::system_clock::to_time_t(newTime);
    std::cout << "Current Time: " << std::put_time(std::localtime(&outCurrTime), "%F %T") << '\n';
    std::cout << "New Time: " << std::put_time(std::localtime(&outNewTime), "%F %T") << '\n';
}
