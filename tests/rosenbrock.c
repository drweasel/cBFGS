/*
 * Copyright (c) 2022 Michael Weitzel <mich@elweitzel.de>
 *
 * This file is licensed under the terms of the standard MIT License;
 * see the file LICENSE for details.
 */
#include <math.h>

static inline double sqr(double v) { return v*v; }

double rosenbrock_f(const double* x, void* unused)
{
	return 10.*sqr(x[1]-sqr(x[0])) + sqr(1 - x[0]);
}

void rosenbrock_Df(const double* x, double* Df, void* unused)
{
	Df[0] = -40.*x[0]*(x[1]-sqr(x[0]))-2.*(1.-x[0]);
	Df[1] = 20.*(x[1]-sqr(x[0]));
}

// vim: fenc=utf-8 noet:
