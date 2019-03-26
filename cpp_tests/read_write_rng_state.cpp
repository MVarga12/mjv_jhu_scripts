// Simple test of writing and reading GSL RNG states for restarting FPR runs 
#include <iostream>
#include <fstream>
#include "gsl/gsl_rng.h"

static gsl_rng* the_generator = nullptr;

double rand_gsl()
{
    if (!the_generator)
        the_generator = gsl_rng_alloc(gsl_rng_taus);
    return gsl_rng_uniform(the_generator);
}

void srand_gsl(int num)
{
    if (!the_generator)
        the_generator = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(the_generator, num);
}

int main() {
    int seed { 288347862 };
    srand(seed);
    rand_gsl();

    auto stateOut = fopen("rng_state.dat", "w");

    for (int i { 0 }; i < 10; ++i)
        rand_gsl();

    gsl_rng_fwrite(stateOut, the_generator);
    fclose(stateOut);

    std::cout << "11th random number: " << rand_gsl() << '\n';
    std::cout << "12th random number: " << rand_gsl() << '\n';
    std::cout << "13th random number: " << rand_gsl() << '\n';
    std::cout << "14th random number: " << rand_gsl() << '\n';

    auto stateIn = fopen("rng_state.dat", "r");
    auto readState = gsl_rng_fread(stateIn, the_generator);
    if (readState == GSL_EFAILED) {
        std::cerr << "Failed to read RNG state file.\n";
        exit(1);
    }

    std::cout << "11th random number: " << rand_gsl() << '\n';
    std::cout << "12th random number: " << rand_gsl() << '\n';
    std::cout << "13th random number: " << rand_gsl() << '\n';
    std::cout << "14th random number: " << rand_gsl() << '\n';
}
