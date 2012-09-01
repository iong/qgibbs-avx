/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 * This file is part of GROMACS.
 * Copyright (c) 2012-  
 *
 * Written by the Gromacs development team under coordination of
 * David van der Spoel, Berk Hess, and Erik Lindahl.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 * 
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */
#ifndef _gmx_math_x86_sse4_1_single_h_
#define _gmx_math_x86_sse4_1_single_h_

#include <stdio.h>
#include <math.h>

//#include "gmx_x86_sse4_1.h"



#ifndef M_PI
#  define M_PI 3.14159265358979323846264338327950288
#endif




/************************
 *                      *
 * Simple math routines *
 *                      *
 ************************/

/* 1.0/sqrt(x) */
static __inline__ __m128
gmx_mm_invsqrt_ps(__m128 x)
{
    const __m128 half  = _mm_set_ps(0.5,0.5,0.5,0.5);
    const __m128 three = _mm_set_ps(3.0,3.0,3.0,3.0);

    __m128 lu = _mm_rsqrt_ps(x);

    return _mm_mul_ps(half,_mm_mul_ps(_mm_sub_ps(three,_mm_mul_ps(_mm_mul_ps(lu,lu),x)),lu));
}


static __inline__ __m128
gmx_mm_abs_ps(__m128 x)
{
    const __m128 signmask  = _mm_castsi128_ps( _mm_set1_epi32(0x7FFFFFFF) );

    return _mm_and_ps(x,signmask);
}


/* Exponential function. This could be calculated from 2^x as Exp(x)=2^(y), where y=log2(e)*x,
 * but there will then be a small rounding error since we lose some precision due to the
 * multiplication. This will then be magnified a lot by the exponential.
 *
 * Instead, we calculate the fractional part directly as a minimax approximation of
 * Exp(z) on [-0.5,0.5]. We use extended precision arithmetics to calculate the fraction
 * remaining after 2^y, which avoids the precision-loss.
 * The final result is correct to within 1 LSB over the entire argument range.
 */
static __m128
gmx_mm_exp_ps(__m128 x)
{
    const __m128  argscale      = _mm_set1_ps(1.44269504088896341f);
    /* Lower bound: Disallow numbers that would lead to an IEEE fp exponent reaching +-127. */
    const __m128  arglimit      = _mm_set1_ps(126.0f);
    const __m128i expbase       = _mm_set1_epi32(127);

    const __m128  invargscale0  = _mm_set1_ps(0.693359375f);
    const __m128  invargscale1  = _mm_set1_ps(-2.12194440e-4f);

    const __m128  CC5           = _mm_set1_ps(1.9875691500e-4f);
    const __m128  CC4           = _mm_set1_ps(1.3981999507e-3f);
    const __m128  CC3           = _mm_set1_ps(8.3334519073e-3f);
    const __m128  CC2           = _mm_set1_ps(4.1665795894e-2f);
    const __m128  CC1           = _mm_set1_ps(1.6666665459e-1f);
    const __m128  CC0           = _mm_set1_ps(5.0000001201e-1f);
    const __m128  one           = _mm_set1_ps(1.0f);

    __m128  y,x2;
    __m128  p0,p1;
    __m128  valuemask;
    __m128i iexppart;
    __m128  fexppart;
    __m128  intpart;

    y = _mm_mul_ps(x,argscale);

    iexppart  = _mm_cvtps_epi32(y);
    intpart   = _mm_round_ps(y,_MM_FROUND_TO_NEAREST_INT);

    iexppart  = _mm_slli_epi32(_mm_add_epi32(iexppart,expbase),23);
    valuemask = _mm_cmpge_ps(arglimit,gmx_mm_abs_ps(y));
    fexppart  = _mm_and_ps(valuemask,_mm_castsi128_ps(iexppart));

    /* Extended precision arithmetics */
    x         = _mm_sub_ps(x,_mm_mul_ps(invargscale0,intpart));
    x         = _mm_sub_ps(x,_mm_mul_ps(invargscale1,intpart));

    x2        = _mm_mul_ps(x,x);

    p1        = _mm_mul_ps(CC5,x2);
    p0        = _mm_mul_ps(CC4,x2);
    p1        = _mm_add_ps(p1,CC3);
    p0        = _mm_add_ps(p0,CC2);
    p1        = _mm_mul_ps(p1,x2);
    p0        = _mm_mul_ps(p0,x2);
    p1        = _mm_add_ps(p1,CC1);
    p0        = _mm_add_ps(p0,CC0);
    p1        = _mm_mul_ps(p1,x);
    p0        = _mm_add_ps(p0,p1);
    p0        = _mm_mul_ps(p0,x2);
    x         = _mm_add_ps(x,one);
    x         = _mm_add_ps(x,p0);

    x         = _mm_mul_ps(x,fexppart);

    return x;
}

#endif /* _gmx_math_x86_sse4_1_single_h_ */
