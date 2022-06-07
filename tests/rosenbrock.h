#pragma once
/*
 * Copyright (c) 2022 Michael Weitzel <mich@elweitzel.de>
 *
 * This file is licensed under the terms of the standard MIT License;
 * see the file LICENSE for details.
 */

double
rosenbrock_f(const double* x, void*);

void
rosenbrock_Df(const double* x, double* Df, void*);

// vim: fenc=utf-8 noet:
