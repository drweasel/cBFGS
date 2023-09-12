/*
 * Copyright (c) 2022 Michael Weitzel <mich@elweitzel.de>
 *
 * This file is licensed under the terms of the standard MIT License;
 * see the file LICENSE for details.
 */
#include "cbfgs/bfgs.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

typedef union
{
    struct
    {
        double R, G, B;
    };
    double vec[3];
} RGB;

typedef union
{
    struct
    {
        double X; // chrominance, [0,1]
        double Y; // luminance, [0,1]
        double Z; // chrominance, [0,1]
    };
    double vec[3];
} XYZ;

typedef struct
{
    RGB obj_rgb;
    RGB* base;
    double* s;
    int dim;
} RGB_Decomp;

RGB
XYZtoRGB(XYZ xyz)
{
    double x = xyz.X;
    double y = xyz.Y;
    double z = xyz.Z;

    RGB rgb;
    rgb.R = 3.240479 * x - 1.537150 * y - 0.498535 * z;
    rgb.G = -0.969256 * x + 1.875992 * y + 0.041556 * z;
    rgb.B = 0.055648 * x - 0.204043 * y + 1.057311 * z;

    for (int i = 0; i < 3; ++i)
    {
        if (rgb.vec[i] < 0.)
        {
            rgb.vec[i] = 0.;
        }
        else if (rgb.vec[i] > 1.)
        {
            rgb.vec[i] = 1.;
        }
    }
    return rgb;
}

XYZ
KtoXYZ(double T, double Y)
{
    double x, y;
    double t1 = 1e3 / T;
    double t2 = 1e6 / (T * T);
    double t3 = 1e9 / (T * T * T);

    if (T > 25000.)
    {
        T = 25000.;
    }
    if (T < 1667.)
    {
        T = 1667.;
    }

    // chromacity
    if (T >= 1667. && T <= 4000.)
    {
        x = -0.2661239 * t3 - 0.2343589 * t2 + 0.8776956 * t1 + 0.179910;
    }
    else // if (T > 4000. && T <= 25000.) */
    {
        x = -3.0258469 * t3 + 2.1070379 * t2 + 0.2226347 * t1 + 0.240390;
    }

    double x2 = x * x;
    double x3 = x * x * x;

    if (T >= 1667. && T < 2222.)
    {
        y = -1.1063814 * x3 - 1.34811020 * x2 + 2.18555832 * x - 0.20219683;
    }
    else if (T >= 2222. && T < 4000)
    {
        y = -0.9549476 * x3 - 1.37418593 * x2 + 2.09137015 * x - 0.16748867;
    }
    else // if (T >= 4000 && T <= 25000.)
    {
        y = +3.0817580 * x3 - 5.87338670 * x2 + 3.75112997 * x - 0.37001483;
    }

    XYZ xyz;
    xyz.X = Y / y * x;
    xyz.Y = Y;
    xyz.Z = Y / y * (1. - x - y);
    return xyz;
}

RGB
KtoRGB(double T, double Y)
{
    return XYZtoRGB(KtoXYZ(T, Y));
}

RGB_Decomp
new_rgb_decomp(RGB obj_rgb, RGB* base, int dim)
{
    RGB_Decomp dcmp;
    for (int k = 0; k < 3; ++k)
    {
        dcmp.obj_rgb.vec[k] = obj_rgb.vec[k];
    }

    dcmp.base = (RGB*)malloc((size_t)dim * sizeof(RGB));
    for (int k = 0; k < dim; ++k)
    {
        dcmp.base[k] = base[k];
    }

    dcmp.s = (double*)malloc((size_t)dim * sizeof(double));
    for (int k = 0; k < dim; ++k)
    {
        dcmp.s[k] = 0.5;
    }

    dcmp.dim = dim;
    return dcmp;
}

void
free_rgb_decomp(RGB_Decomp* dcmp)
{
    free(dcmp->base);
    free(dcmp->s);
}

double
rgb_decomp_obj_f(const double* s, void* uptr)
{
    RGB_Decomp* dcmp = (RGB_Decomp*)uptr;
    double e = 0.;
    for (int i = 0; i < 3; ++i)
    {
        double e_i = -dcmp->obj_rgb.vec[i];
        for (int j = 0; j < dcmp->dim; ++j)
        {
            e_i += dcmp->base[j].vec[i] * s[j];
        }
        e += e_i * e_i;
    }
    return e;
}

void
rgb_decomp_obj_Df(const double* s, double* Df, void* uptr)
{
    const RGB_Decomp* dcmp = uptr;
    RGB e;
    e.R = e.G = e.B = 0.;
    for (int i = 0; i < 3; ++i)
    {
        double e_i = -dcmp->obj_rgb.vec[i];
        for (int j = 0; j < dcmp->dim; ++j)
        {
            e_i += dcmp->base[j].vec[i] * s[j];
        }
        e.vec[i] = e_i;
    }
    for (int i = 0; i < dcmp->dim; ++i)
    {
        Df[i] = 0.;
        for (int j = 0; j < 3; ++j)
        {
            Df[i] += dcmp->base[i].vec[j] * e.vec[j];
        }
        Df[i] *= 2.;
    }
}

bool
decompose_RGB(RGB_Decomp* dcmp)
{
    BFGS_Result result;
    if (dcmp == NULL || dcmp->dim <= 0)
        return false;

    for (int k = 0; k < dcmp->dim; ++k)
        dcmp->s[k] = 0.5;

    void* workspace = malloc(bfgs_WORKSPACE_SIZE(dcmp->dim));
    result = bfgs(
      rgb_decomp_obj_f,
      rgb_decomp_obj_Df,
      dcmp,
      dcmp->s,
      (unsigned int)dcmp->dim,
      1000,
      workspace);
    free(workspace);

    for (int k = 0; k < dcmp->dim; ++k)
    {
        if (dcmp->s[k] < 0.)
            dcmp->s[k] = 0.;
        else if (dcmp->s[k] > 1.)
            dcmp->s[k] = 1.;
    }

    return is_bfgs_success(&result);
}

bool
test_rgb2rgbw(void)
{
    RGB base[] = { { { 0.9, 0.1, 0.05 } },
                   { { 0.0, 0.9, 0.1 } },
                   { { 0.05, 0.01, 1 } },
                   { { 0.95, 0.9, 1 } } };
    RGB obj = { { 0.6, 0.3, 0.7 } };

    RGB_Decomp dcmp = new_rgb_decomp(obj, base, 4);

    double x[4] = { 0.5, 0.5, 0.5, 0.5 };

    unsigned char workspace[bfgs_WORKSPACE_SIZE(4)];
    BFGS_Result result =
      bfgs(rgb_decomp_obj_f, rgb_decomp_obj_Df, &dcmp, x, 4, 1000, workspace);

    printf("x = (%g,%g,%g,%g)\n", x[0], x[1], x[2], x[3]);

    RGB probe = { { 0., 0., 0. } };
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < dcmp.dim; ++j)
        {
            probe.vec[i] += dcmp.base[j].vec[i] * x[j];
        }
    }
    printf("probe: (%g,%g,%g)\n", probe.R, probe.G, probe.B);

#if 0
    if (decompose_RGB(&dcmp))
    {
        printf("success:");
        for (int i = 0; i < dcmp.dim; ++i)
            printf(" %g", dcmp.s[i]);
        printf("\n");

        RGB probe = { 0., 0., 0. };
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < dcmp.dim; ++j)
            {
                probe.vec[i] += dcmp.base[j].vec[i] * dcmp.s[j];
            }
        }
        printf("probe: (%g,%g,%g)\n", probe.R, probe.G, probe.B);
    }
#    if 0
	for (double K=2500.; K<8000.; K+=500.)
	{
		RGB bbr = KtoRGB(K,0.75);
		printf("%.0f K => (%.0f,%.0f,%.0f)\n", K, bbr.R*255., bbr.G*255., bbr.B*255.);
	}
#    endif
#endif
    free_rgb_decomp(&dcmp);
    return is_bfgs_success(&result);
}

// vim: fenc=utf-8 et:
