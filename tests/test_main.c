/*
 * Copyright (c) 2022 Michael Weitzel <mich@elweitzel.de>
 *
 * This file is licensed under the terms of the standard MIT License;
 * see the file LICENSE for details.
 */
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

bool
test_rosenbrock(void);

bool
test_rgb2rgbw(void);

int
main(void)
{
    int result = EXIT_SUCCESS;

    printf("\n=== Rosenbrock function ...\n");
    if (!test_rosenbrock())
    {
        result = EXIT_FAILURE;
        printf("=== FAILED!\n");
    }
    else
    {
        printf("=== ok.\n");
    }

    printf("\n=== RGB2RGBW conversion ...\n");
    if (!test_rgb2rgbw())
    {
        result = EXIT_FAILURE;
        printf("=== FAILED!\n");
    }
    else
    {
        printf("=== ok.\n");
    }

    return result;
}

// vim: fenc=utf-8 et:
