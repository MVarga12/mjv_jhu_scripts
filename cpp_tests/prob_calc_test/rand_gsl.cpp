#include "rand_gsl.hpp"
// #include "gsl/gsl_rng.h"
#include <iostream>
#include <cmath>

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

double GaussV() {
    double R {2.0};
    double V1{};

    while (R >= 1.0) {
        V1 = 2.0 * rand_gsl() - 1.0;
        double V2 = 2.0 * rand_gsl() - 1.0;
        R = (V1 * V1) + (V2 * V2);
    }
    return (V1 * sqrt(-2.0 * log(R) / R));
}
