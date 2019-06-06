#include <iostream>
#include <chrono>
#include <vector>


using Timer = std::chrono::steady_clock;
int main() {
    std::vector<int> test = std::vector<int>(10);
    for (auto elem : test)
        std::cout << ' ' << elem;
    std::cout << '\n';
    // int numInts { 1000000 };
    // Timer::time_point start = Timer::now();
    //
    // std::vector<int> vec1 { };
    // for (int i { 0 }; i < numInts; ++i) {
    //     vec1.push_back(i);
    // }
    // std::cout << "No reserve: " << std::chrono::duration<double>(Timer::now() - start).count() << '\n';
    //
    // Timer::time_point start2 = Timer::now();
    // std::vector<int> vec2 { };
    // vec2.reserve(numInts);
    // for (int i { 0 }; i < numInts; ++i) {
    //     vec2.push_back(i);
    // }
    // std::cout << "Reserve: " << std::chrono::duration<double>(Timer::now() - start2).count() << '\n';
    // std::cout << "Vec2 capacity before: " << vec2.capacity() << '\n';
    // vec2.clear();
    // std::cout << "Vec2 capacity after: " << vec2.capacity() << '\n';
}
