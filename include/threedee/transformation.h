#ifndef THREEDEE_TRANSFORMATION_H
#define THREEDEE_TRANSFORMATION_H

#include <threedee/vector.h>
#include <threedee/matrix.h>

static inline mat4 mtranslate_zero(vec4 delta) __attribute__((always_inline));
static inline mat4 mtranslate_zero(vec4 delta)
{
    vec4 one = _mm_set_ss(1.0);
    vec4 row0 = vshuffle(one, delta, 0, 1, 3, 0);
    vec4 row1 = vshuffle(one, delta, 1, 0, 3, 1);
    vec4 row2 = vshuffle(one, vshuffle(one, delta, 0, 1, 2, 2), 1, 1, 0, 2);
    vec4 row3 = vshuffle(one, one, 1, 1, 1, 0);

    mat4 m = {{ row0, row1, row2, row3 }};
    return m;
}

static inline mat4 mtranslate(vec4 delta) __attribute__((always_inline));
static inline mat4 mtranslate(vec4 delta)
{
    return mtranslate_zero(vxyz(delta));
}

static inline mat4 mscale_zero(vec4 factor) __attribute__((always_inline));
static inline mat4 mscale_zero(vec4 factor)
{
    vec4 one = _mm_set_ss(1.0);
    vec4 row0 = vshuffle(factor, vzero(), 0, 3, 0, 0);
    vec4 row1 = vshuffle(factor, vzero(), 3, 1, 0, 0);
    vec4 row2 = vshuffle(vzero(), factor, 0, 0, 2, 3);
    vec4 row3 = vshuffle(one, one, 1, 1, 1, 0);
    mat4 m = {{ row0, row1, row2, row3 }};
    return m;
}

static inline mat4 mscale(vec4 factor) __attribute__((always_inline));
static inline mat4 mscale(vec4 factor)
{
    return mscale_zero(vxyz(factor));
}

#endif
