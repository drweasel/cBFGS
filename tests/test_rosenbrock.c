/*
 * Copyright (c) 2022 Michael Weitzel <mich@elweitzel.de>
 *
 * This file is licensed under the terms of the standard MIT License;
 * see the file LICENSE for details.
 */
#include "cbfgs/bfgs.h"
#include "rosenbrock.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

bool
test_rosenbrock()
{
    unsigned char workspace[bfgs_WORKSPACE_SIZE(2)];

    double x[2] = { -3., -4. };
    Rosenbrock rb_param = { 1., 100. };

    BFGS_Result result =
      bfgs(rosenbrock_f, rosenbrock_Df, &rb_param, x, 2, 100, workspace);

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
