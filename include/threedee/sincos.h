#ifndef THREEDEE_SINCOS_H
#define THREEDEE_SINCOS_H

#include <threedee/vector.h>

#ifdef stdmathlib

#include <math.h>

static inline vec4 vsincos(vec4 x, vec4 * restrict s, vec4 * restrict c)
{
    float *xs = (float*)&x;
    *s = vec(sin(xs[0]), sin(xs[1]), sin(xs[2]), sin(xs[3]));
    *c = vec(cos(xs[0]), cos(xs[1]), cos(xs[2]), cos(xs[3]));
}

#else

#include <mipp/math/sse_mathfun.h>

static inline void vsincos(vec4 x, vec4 * restrict s, vec4 * restrict c)
{
    sincos_ps(x, s, c);
}

#endif

#endif
