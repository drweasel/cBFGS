/*
 * Copyright (c) 2022 Michael Weitzel <mich@elweitzel.de>
 *
 * This file is licensed under the terms of the standard MIT License;
 * see the file LICENSE for details.
 */

#include "finite_diff.h"
#include "stdlib.h"

#define MACHEPS 2.2204e-16

#if 0
// Halley iteration for approximating the function a^(1/order)
// up to machine precision. This can be evaluated
static inline double
orderth_root(const unsigned int _order, const double a)
{
    const double order = (const double)_order;
    double xj = a;

    for (;;)
    {
        // general Halley step for f(x) = x^k-a =(!) 0 => x = a^(1/k)
        double xj_k = xj;
        for (unsigned int l = 0; l < _order - 1; ++l)
            xj_k *= xj;

        double xj_n = order * a * xj / (xj_k + (order - 1.0) * a);
        if (xj_n == xj)
            return xj_n;
        xj = xj_n;
    }
}
#endif

static double
finite_diff(
  double (*f)(const double*, void*),
  void* user_ptr,
  double* x,
  const unsigned int k)
{
    const unsigned int L = 1;
    const unsigned int R = 1;
    const unsigned int order = L + R;
#if 0
    const double h = (order == 1) ? 3.4526698e-4 : orderth_root(order, MACHEPS);
#endif
    static const double all_h[] = { 0.,
                                    0.00034526698,
                                    1.490116119384766e-08,
                                    6.055454452393343e-06,
                                    0.0001220703125,
                                    0.0007400959797414051,
                                    0.002460783300575925,
                                    0.005804665191941207,
                                    0.01104854345603981,
                                    0.01822701624337682,
                                    0.02720470510300388 };
    const double h = all_h[order];

    static const double all_coeffs[][6][11] = {
        // clang-format off
		{ // L = 0
			{ 0 }, // R = 0
			{ -1, 1 }, // R = 1
			{ -3/2., 2, -1/2. }, // R = 2
			{ -11/6., 3, -3/2., 1/3. }, // R = 3
			{ -25/12., 4, -3, 4/3., -1/4. }, // R = 4
			{ -137/60., 5, -5, 10/3., -5/4., 1/5. } // R = 5
		},
		{ // L = 1
			{ -1, 1. }, // R = 0
			{ -1/2., 0, 1/2. }, // R = 1
			{ -1/3., -1/2., 1, -1/6. }, // R = 2
			{ -1/4., -5/6., 3/2., -1/2., 1/12. }, // R = 3
			{ -1/5., -13/12., 2, -1, 1/3., -1/20. }, // R = 4
			{ -1/6., -77/60., 5/2., -5/3., 5/6., -1/4., 1/30. } // R = 5
		},
		{ // L = 2
			{ 1/2., -2, 3/2. }, // R = 0
			{ 1/6., -1, 1/2., 1/3. }, // R = 1
			{ 1/12., -2/3., 0, 2/3., -1/12. }, // R = 2
			{ 1/20., -1/2., -1/3., 1, -1/4., 1/30. }, // R = 3
			{ 1/30., -2/5., -7/12., 4/3., -1/2., 2/15., -1/60. }, // R = 4
			{ 1/42., -1/3., -47/60., 5/3., -5/6., 1/3., -1/12., 1/105. } // R = 5
		},
		{ // L = 3
			{ -1/3., 3/2., -3, 11/6. }, // R = 0
			{ -1/12., 1/2., -3/2., 5/6., 1/4. }, // R = 1
			{ -1/30., 1/4., -1, 1/3., 1/2., -1/20. }, // R = 2
			{ -1/60., 3/20., -3/4., 0, 3/4., -3/20., 1/60. }, // R = 3
			{ -1/105., 1/10., -3/5., -1/4., 1, -3/10., 1/15., -1/140. }, // R = 4
			{ -1/168., 1/14., -1/2., -9/20., 5/4., -1/2., 1/6., -1/28., 1/280. } // R = 5
		},
		{ // L = 4
			{ 1/4., -4/3., 3, -4, 25/12. }, // R = 0
			{ 1/20., -1/3., 1, -2, 13/12., 1/5. }, // R = 1
			{ 1/60., -2/15., 1/2., -4/3., 7/12., 2/5., -1/30. }, // R = 2
			{ 1/140., -1/15., 3/10., -1, 1/4., 3/5., -1/10., 1/105. }, // R = 3
			{ 1/280., -4/105., 1/5., -4/5., 0, 4/5., -1/5., 4/105., -1/280. }, // R = 4
			{ 1/504., -1/42., 1/7., -2/3., -1/5., 1, -1/3., 2/21., -1/56., 1/630. } // R = 5
		},
		{ // L = 5
			{ -1/5., 5/4., -10/3., 5, -5, 137/60. }, // R = 0
			{ -1/30., 1/4., -5/6., 5/3., -5/2., 77/60., 1/6. }, // R = 1
			{ -1/105., 1/12., -1/3., 5/6., -5/3., 47/60., 1/3., -1/42. }, // R = 2
			{ -1/280., 1/28., -1/6., 1/2., -5/4., 9/20., 1/2., -1/14., 1/168. }, // R = 3
			{ -1/630., 1/56., -2/21., 1/3., -1, 1/5., 2/3., -1/7., 1/42., -1/504. }, // R = 4
			{ -1/1260., 5/504., -5/84., 5/21., -5/6., 0, 5/6., -5/21., 5/84., -5/504., 1/1260. } // R = 5
		}
        // clang-format on
    };

    const double* coeff = all_coeffs[L][R];
    const double v = x[k];
    double result = 0.;

    for (int j = -(int)L; j <= (int)R; ++j)
    {
        if (coeff[(unsigned int)j + L] != 0.)
        {
            x[k] = v + (double)j * h;
            result += f(x, user_ptr) * coeff[(unsigned int)j + L];
        }
    }
    x[k] = v;
    return result / h;
}

/**
 * Computes a function's gradient based on finite differences.
 *
 * @param f function pointer
 * @param user_ptr a user-supplied pointer passed to function f
 * @param x the current parameters of f
 * @param df the computed gradient of f at x
 * @param n the number of arguments in x
 */
static void
fd_gradient(
  double (*f)(const double*, void*),
  void* user_ptr,
  double* x,
  double* df,
  const unsigned int n)
{
    for (unsigned int k = 0; k < n; ++k)
        df[k] = finite_diff(f, user_ptr, x, k);
}

static double
fdwrapper_f_(const double* x, void* user_ptr)
{
    FDWrapper* fdw = user_ptr;
    return fdw->orig_f(x, fdw->user_ptr);
}

static void
fdwrapper_Df_(const double* x, double* df, void* user_ptr)
{
    FDWrapper* fdw = user_ptr;
    fd_gradient(fdw->orig_f, fdw->user_ptr, (double*)x /* ouch! */, df, fdw->n);
}

FDWrapper
make_fd_wrapper(
  double (*f)(const double*, void*),
  void* user_ptr,
  unsigned int n)
{
    FDWrapper fdw;
    fdw.user_ptr = user_ptr;
    fdw.orig_f = f;
    fdw.n = n;
    fdw.f = fdwrapper_f_;
    fdw.Df = fdwrapper_Df_;
    return fdw;
}

// vim: fenc=utf-8 et:
