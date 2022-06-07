#pragma once
/*
 * Copyright (c) 2022 Michael Weitzel <mich@elweitzel.de>
 *
 * This file is licensed under the terms of the standard MIT License;
 * see the file LICENSE for details.
 */

#include <stdbool.h>

/**
 * Quasi-Newton optimizer implementing the BFGS algorithm for
 * unconstraint minization.
 *
 * @param f objective function
 * @param Df gradient vector of the objective function
 * @param user_ptr a user pointer passed to the function handles
 * @param x parameter vector
 * @param n number of parameters
 * @param max_iter maximum number of iterations
 *
 * @author Michael Weitzel <mich@elweitzel.de>
 */
bool
bfgs(
  double (*f)(const double*, void*),
  void (*Df)(const double*, double*, void*),
  void* user_ptr,
  double* x,
  const int n,
  const int max_iter);

// vim: fenc=utf-8 noet:
