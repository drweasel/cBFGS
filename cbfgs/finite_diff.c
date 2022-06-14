#include "finite_diff.h"
#include "stdlib.h"

#define MACHEPS 2.2204e-16

// constexpr Halley iteration for approximating the function a^(1/order)
// up to machine precision
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
    const double h = (order == 1) ? 3.4526698e-4 : orderth_root(order, MACHEPS);

    typedef double V;
    const double all_coeffs[][6][11] = {
        // clang-format off
		{ // L = 0
			{ 0 }, // R = 0
			{ -1, 1 }, // R = 1
			{ -3/(V)(2), 2, -1/(V)(2) }, // R = 2
			{ -11/(V)(6), 3, -3/(V)(2), 1/(V)(3) }, // R = 3
			{ -25/(V)(12), 4, -3, 4/(V)(3), -1/(V)(4) }, // R = 4
			{ -137/(V)(60), 5, -5, 10/(V)(3), -5/(V)(4), 1/(V)(5) } // R = 5
		},
		{ // L = 1
			{ -1, (V)(1) }, // R = 0
			{ -1/(V)(2), 0, 1/(V)(2) }, // R = 1
			{ -1/(V)(3), -1/(V)(2), 1, -1/(V)(6) }, // R = 2
			{ -1/(V)(4), -5/(V)(6), 3/(V)(2), -1/(V)(2), 1/(V)(12) }, // R = 3
			{ -1/(V)(5), -13/(V)(12), 2, -1, 1/(V)(3), -1/(V)(20) }, // R = 4
			{ -1/(V)(6), -77/(V)(60), 5/(V)(2), -5/(V)(3), 5/(V)(6), -1/(V)(4), 1/(V)(30) } // R = 5
		},
		{ // L = 2
			{ 1/(V)(2), -2, 3/(V)(2) }, // R = 0
			{ 1/(V)(6), -1, 1/(V)(2), 1/(V)(3) }, // R = 1
			{ 1/(V)(12), -2/(V)(3), 0, 2/(V)(3), -1/(V)(12) }, // R = 2
			{ 1/(V)(20), -1/(V)(2), -1/(V)(3), 1, -1/(V)(4), 1/(V)(30) }, // R = 3
			{ 1/(V)(30), -2/(V)(5), -7/(V)(12), 4/(V)(3), -1/(V)(2), 2/(V)(15), -1/(V)(60) }, // R = 4
			{ 1/(V)(42), -1/(V)(3), -47/(V)(60), 5/(V)(3), -5/(V)(6), 1/(V)(3), -1/(V)(12), 1/(V)(105) } // R = 5
		},
		{ // L = 3
			{ -1/(V)(3), 3/(V)(2), -3, 11/(V)(6) }, // R = 0
			{ -1/(V)(12), 1/(V)(2), -3/(V)(2), 5/(V)(6), 1/(V)(4) }, // R = 1
			{ -1/(V)(30), 1/(V)(4), -1, 1/(V)(3), 1/(V)(2), -1/(V)(20) }, // R = 2
			{ -1/(V)(60), 3/(V)(20), -3/(V)(4), 0, 3/(V)(4), -3/(V)(20), 1/(V)(60) }, // R = 3
			{ -1/(V)(105), 1/(V)(10), -3/(V)(5), -1/(V)(4), 1, -3/(V)(10), 1/(V)(15), -1/(V)(140) }, // R = 4
			{ -1/(V)(168), 1/(V)(14), -1/(V)(2), -9/(V)(20), 5/(V)(4), -1/(V)(2), 1/(V)(6), -1/(V)(28), 1/(V)(280) } // R = 5
		},
		{ // L = 4
			{ 1/(V)(4), -4/(V)(3), 3, -4, 25/(V)(12) }, // R = 0
			{ 1/(V)(20), -1/(V)(3), 1, -2, 13/(V)(12), 1/(V)(5) }, // R = 1
			{ 1/(V)(60), -2/(V)(15), 1/(V)(2), -4/(V)(3), 7/(V)(12), 2/(V)(5), -1/(V)(30) }, // R = 2
			{ 1/(V)(140), -1/(V)(15), 3/(V)(10), -1, 1/(V)(4), 3/(V)(5), -1/(V)(10), 1/(V)(105) }, // R = 3
			{ 1/(V)(280), -4/(V)(105), 1/(V)(5), -4/(V)(5), 0, 4/(V)(5), -1/(V)(5), 4/(V)(105), -1/(V)(280) }, // R = 4
			{ 1/(V)(504), -1/(V)(42), 1/(V)(7), -2/(V)(3), -1/(V)(5), 1, -1/(V)(3), 2/(V)(21), -1/(V)(56), 1/(V)(630) } // R = 5
		},
		{ // L = 5
			{ -1/(V)(5), 5/(V)(4), -10/(V)(3), 5, -5, 137/(V)(60) }, // R = 0
			{ -1/(V)(30), 1/(V)(4), -5/(V)(6), 5/(V)(3), -5/(V)(2), 77/(V)(60), 1/(V)(6) }, // R = 1
			{ -1/(V)(105), 1/(V)(12), -1/(V)(3), 5/(V)(6), -5/(V)(3), 47/(V)(60), 1/(V)(3), -1/(V)(42) }, // R = 2
			{ -1/(V)(280), 1/(V)(28), -1/(V)(6), 1/(V)(2), -5/(V)(4), 9/(V)(20), 1/(V)(2), -1/(V)(14), 1/(V)(168) }, // R = 3
			{ -1/(V)(630), 1/(V)(56), -2/(V)(21), 1/(V)(3), -1, 1/(V)(5), 2/(V)(3), -1/(V)(7), 1/(V)(42), -1/(V)(504) }, // R = 4
			{ -1/(V)(1260), 5/(V)(504), -5/(V)(84), 5/(V)(21), -5/(V)(6), 0, 5/(V)(6), -5/(V)(21), 5/(V)(84), -5/(V)(504), 1/(V)(1260) } // R = 5
		}
        // clang-format on
    };

    const double* coeff = all_coeffs[L][R];
    const double v = x[k];
    double result = 0.;

    for (int j = -(int)L; j <= (int)R; ++j)
    {
        if (coeff[(unsigned int)j] + L != 0.)
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
    FDWrapper* fdw = (FDWrapper*)user_ptr;
    return fdw->orig_f(x, fdw->user_ptr);
}

static void
fdwrapper_Df_(const double* x, double* df, void* user_ptr)
{
    FDWrapper* fdw = (FDWrapper*)user_ptr;
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
