/*
 * Copyright (c) 2022 Michael Weitzel <mich@elweitzel.de>
 *
 * This file is licensed under the terms of the standard MIT License;
 * see the file LICENSE for details.
 */
#include "rosenbrock.h"
#include <math.h>

static inline double
sqr(double v)
{
    return v * v;
}

double
rosenbrock_f(const double* x, void* param)
{
    const Rosenbrock* rb = (const Rosenbrock*)param;

    return rb->b * sqr(x[1] - sqr(x[0])) + sqr(rb->a - x[0]);
}

void
rosenbrock_Df(const double* x, double* Df, void* param)
{
    const Rosenbrock* rb = (const Rosenbrock*)param;

    Df[0] = -4. * rb->b * x[0] * (x[1] - sqr(x[0])) - 2. * (rb->a - x[0]);
    Df[1] = 2. * rb->b * (x[1] - sqr(x[0]));
}

// vim: fenc=utf-8 et:
