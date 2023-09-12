/*
 * Copyright (c) 2022 Michael Weitzel <mich@elweitzel.de>
 *
 * This file is licensed under the terms of the standard MIT License;
 * see the file LICENSE for details.
 */
#include "cbfgs/bfgs.h"
#include "cbfgs/finite_diff.h"
#include "rosenbrock.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

bool
test_rosenbrock(void)
{
    unsigned char workspace[bfgs_WORKSPACE_SIZE(2)];

    double x[2] = { 3, -4. };
    Rosenbrock rb_param = { 1., 100. };

#if 0
    BFGS_Result result =
      bfgs(rosenbrock_f, rosenbrock_Df, &rb_param, x, 2, 100, workspace);
#elif 0
    FDWrapper fdw = make_fd_wrapper(rosenbrock_f, &rb_param, 2);
    BFGS_Result result =
      bfgs(fdw.f, fdw.Df, &fdw, x, 2, 100, workspace);
#else
    BFGS_Result result =
      bfgs_fd(rosenbrock_f, &rb_param, x, 2, 10000, workspace);
#endif

    printf(
      "result=%s; #iters=%i; step_length=%g, x = [ %g, %g ]\n",
      is_bfgs_success(&result) ? "true" : "false",
      result.niters,
	  result.step_length,
      x[0],
      x[1]);

    return is_bfgs_success(&result);
}

// vim: fenc=utf-8 et:
