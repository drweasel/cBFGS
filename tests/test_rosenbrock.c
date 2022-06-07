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
	unsigned char workspace[bfgs_WORKSPACE(2)];
    double x[2] = { -14., 10. };
    bool result =
      bfgs(rosenbrock_f, rosenbrock_Df, NULL, x, 2, 1000, workspace);
    printf(
      "result=%s; x = [ %g, %g ]\n", result ? "true" : "false", x[0], x[1]);

    return result;
}

// vim: fenc=utf-8 noet:
