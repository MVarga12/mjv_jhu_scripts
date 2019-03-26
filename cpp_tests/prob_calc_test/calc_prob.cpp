#include "rand_gsl.hpp"

#include <chrono>
#include <iostream>
#include <random>

struct Parms {
    double usToS { 1E-6 };
    double rate { 100 };
    double dt { 0.1 };

    Parms() = default;
};

using Timer = std::chrono::steady_clock;

int main() {
    // Set up rng
    std::random_device rd {};
    unsigned seed { rd() };
    srand(seed);

    // Parameters, the values don't matter
    Parms params{};

    // probability
    const double probConst { 1 - exp(-params.rate * params.dt * params.usToS) };

    int i { 0 };
    
    // Set up timer
    Timer::time_point start { Timer::now() };
    while (++i < 10000) {
        for (unsigned itr { 0 }; itr < 10000; ++itr) {
            double probConst { 1 - exp(-params.rate * params.dt * params.usToS) };
            double rand { 1 * rand_gsl() };
        }
    }
    std::cout << "With prob. calculation every time: " << std::chrono::duration<double>(Timer::now() - start).count() << '\n';


    i = 0;
    //
    // Set up timer
    Timer::time_point start2 { Timer::now() };

    while (++i < 10000) {
        for (unsigned itr { 0 }; itr < 10000; ++itr) {
            double prob { probConst };
            double rand { 1 * rand_gsl() };
        }
    }
    std::cout << "With calculation once: " << std::chrono::duration<double>(Timer::now() - start2).count() << '\n';
}
