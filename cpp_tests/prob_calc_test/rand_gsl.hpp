#pragma once

#include "gsl/gsl_rng.h"

/*!
 * \brief Uses the previously initialized random number generator to return a random number
 * \param[out] double Uniformly distributed random double.
 */
double rand_gsl();

/*!
 * \brief Initializes the GSL random number generator.
 */
void srand_gsl(int);

/*!
 * \brief Uses Box-Mueller method to greate Gaussian-distributed random numbers from a uniform random number generator.
 * \param[out] double Gaussian-distributed random double.
 */
double GaussV();
