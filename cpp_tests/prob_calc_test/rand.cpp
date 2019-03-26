#include "rand_gsl.hpp"
#include <iostream>
#include <random>

int main() {
    std::random_device rd{};

    unsigned i { 0 };
    while (++i < 100) {
        std::cout << rd() << '\n';
    }
}
