#pragma once
/*
 * Copyright (c) 2022 Michael Weitzel <mich@elweitzel.de>
 *
 * This file is licensed under the terms of the standard MIT License;
 * see the file LICENSE for details.
 */

#include <stdbool.h>

/** Computes the size of the required workspace, in bytes. */
#define bfgs_WORKSPACE_SIZE(n) \
    (((unsigned)(n) * ((unsigned)(n) + 1) + 6 * (unsigned)(n)) * sizeof(double))

typedef struct
{
    enum BFGS_State
    {
        bfgs_invalid_arguments = -3,
        bfgs_max_iter_exceeded = -2,
        bfgs_not_a_minimum = -1,
        bfgs_unknown = 0,
        bfgs_success_small_gradient = 1,
        bfgs_success_small_step = 2
    } state;
    unsigned int niters;
    double step_length;
} BFGS_Result;

/**
 * Tests the previous optimisation run for success.
 */
bool
is_bfgs_success(const BFGS_Result* result);

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
 * @param workspace pre-allocated memory (see bfgs_WORKSPACE)
 */
BFGS_Result
bfgs(
  double (*f)(const double*, void*),
  void (*Df)(const double*, double*, void*),
  void* user_ptr,
  double* x,
  const unsigned int n,
  const unsigned int max_iter,
  void* workspace);

/**
 * Convenience wrapper for bfgs which computes derivated based on the finite
 * difference method.
 *
 * @see bfgs
 */
BFGS_Result
bfgs_fd(
  double (*f)(const double*, void*),
  void* user_ptr,
  double* x,
  const unsigned int n,
  const unsigned int max_iter,
  void* workspace);

// vim: fenc=utf-8 et:
