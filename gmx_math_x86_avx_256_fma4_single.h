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
 *
 * Truncated by Ionut Georgescu to the bare minimum.
 */
#ifndef _gmx_math_x86_avx_256_fma4_single_h_
#define _gmx_math_x86_avx_256_fma4_single_h_

//#include "gmx_x86_avx_256.h"


#ifndef M_PI
#  define M_PI 3.14159265358979323846264338327950288
#endif



/************************
 *                      *
 * Simple math routines *
 *                      *
 ************************/

/* 1.0/sqrt(x), 256-bit wide version */
static __inline__ __m256
gmx_mm256_invsqrt_ps(__m256 x)
{
    const __m256 half  = _mm256_set1_ps(0.5f);
    const __m256 three = _mm256_set1_ps(3.0f);

    __m256 lu = _mm256_rsqrt_ps(x);

    return _mm256_mul_ps(half,_mm256_mul_ps(_mm256_sub_ps(three,_mm256_mul_ps(_mm256_mul_ps(lu,lu),x)),lu));
}

static __inline__ __m256
gmx_mm256_abs_ps(__m256 x)
{
    const __m256 signmask  = _mm256_castsi256_ps( _mm256_set1_epi32(0x7FFFFFFF) );

    return _mm256_and_ps(x,signmask);
}


/* Exponential function, 256 bit wide. This could be calculated from 2^x as Exp(x)=2^(y), 
 * where y=log2(e)*x, but there will then be a small rounding error since we lose some 
 * precision due to the multiplication. This will then be magnified a lot by the exponential.
 *
 * Instead, we calculate the fractional part directly as a minimax approximation of
 * Exp(z) on [-0.5,0.5]. We use extended precision arithmetics to calculate the fraction
 * remaining after 2^y, which avoids the precision-loss.
 * The final result is correct to within 1 LSB over the entire argument range.
 */
static __m256
gmx_mm256_exp_ps(__m256 x)
{
    const __m256  argscale      = _mm256_set1_ps(1.44269504088896341f);
    /* Lower bound: Disallow numbers that would lead to an IEEE fp exponent reaching +-127. */
    const __m256  arglimit      = _mm256_set1_ps(126.0f);
    const __m128i expbase       = _mm_set1_epi32(127);

    const __m256  invargscale0  = _mm256_set1_ps(0.693359375f);
    const __m256  invargscale1  = _mm256_set1_ps(-2.12194440e-4f);

    const __m256  CE5           = _mm256_set1_ps(1.9875691500e-4f);
    const __m256  CE4           = _mm256_set1_ps(1.3981999507e-3f);
    const __m256  CE3           = _mm256_set1_ps(8.3334519073e-3f);
    const __m256  CE2           = _mm256_set1_ps(4.1665795894e-2f);
    const __m256  CE1           = _mm256_set1_ps(1.6666665459e-1f);
    const __m256  CE0           = _mm256_set1_ps(5.0000001201e-1f);
    const __m256  one           = _mm256_set1_ps(1.0f);

    __m256  x2,exp2arg;
    __m256  p0,p1;
    __m256  valuemask;
    __m256i iexppart;
    __m128i iexppart128a,iexppart128b;
    __m256  fexppart;
    __m256  intpart;

    exp2arg = _mm256_mul_ps(x,argscale);

    iexppart  = _mm256_cvtps_epi32(exp2arg);
    intpart   = _mm256_round_ps(exp2arg,_MM_FROUND_TO_NEAREST_INT);

    iexppart128b = _mm256_extractf128_si256(iexppart,0x1);
    iexppart128a = _mm256_castsi256_si128(iexppart);

    iexppart128a = _mm_slli_epi32(_mm_add_epi32(iexppart128a,expbase),23);
    iexppart128b = _mm_slli_epi32(_mm_add_epi32(iexppart128b,expbase),23);

    iexppart  = _mm256_castsi128_si256(iexppart128a);
    iexppart  = _mm256_insertf128_si256(iexppart,iexppart128b,0x1);
    valuemask = _mm256_cmp_ps(arglimit,gmx_mm256_abs_ps(exp2arg),_CMP_GE_OQ);
    fexppart  = _mm256_and_ps(valuemask,_mm256_castsi256_ps(iexppart));

    /* Extended precision arithmetics */
    x         = _mm256_nmacc_ps(invargscale0,intpart,x);
    x         = _mm256_nmacc_ps(invargscale1,intpart,x);

    x2        = _mm256_mul_ps(x,x);

    p1        = _mm256_macc_ps(CE5,x2,CE3);
    p0        = _mm256_macc_ps(CE4,x2,CE2);
    p1        = _mm256_macc_ps(p1,x2,CE1);
    p0        = _mm256_macc_ps(p0,x2,CE0);
    p0        = _mm256_macc_ps(p1,x,p0);
    p0        = _mm256_macc_ps(p0,x2,one);

    x         = _mm256_add_ps(x,p0);

    x         = _mm256_mul_ps(x,fexppart);

    return x;
}

#endif
