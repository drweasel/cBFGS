/*
 * Copyright (c) 2022 Michael Weitzel <mich@elweitzel.de>
 *
 * This file is licensed under the terms of the standard MIT License;
 * see the file LICENSE for details.
 */
#include <stdbool.h>
#include <stdio.h>

bool
test_rosenbrock();

bool
test_rgb2rgbw();

int
main()
{
    int result = 0;

    printf("\n=== Rosenbrock function ...\n");
    if (!test_rosenbrock())
    {
        result = 1;
        printf("=== FAILED!\n");
    }
    else
        printf("=== ok.\n");

    printf("\n=== RGB2RGBW conversion ...\n");
    if (!test_rgb2rgbw())
    {
        result = 1;
        printf("=== FAILED!\n");
    }
    else
        printf("=== ok.\n");

    return result;
}

// vim: fenc=utf-8 et:
