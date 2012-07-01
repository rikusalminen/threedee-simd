#ifndef THREEDEE_TRANSFORMATION_H
#define THREEDEE_TRANSFORMATION_H

#include <threedee/vector.h>
#include <threedee/matrix.h>

static inline mat4 mtranslate_zero(vec4 delta) __attribute__((always_inline));
static inline mat4 mtranslate_zero(vec4 delta)
{
    vec4 one = _mm_set_ss(1.0);
    vec4 row0 = vshuffle(one, one, 0, 1, 1, 3);
    vec4 row1 = vshuffle(one, one, 1, 0, 1, 3);
    vec4 row2 = vshuffle(one, one, 1, 1, 0, 3);
    vec4 row3 = vshuffle(delta, vshuffle(delta, one, 2, 2, 0, 1), 0, 1, 0, 2);

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
    vec4 col0 = vshuffle(factor, vzero(), 0, 3, 0, 0);
    vec4 col1 = vshuffle(factor, vzero(), 3, 1, 0, 0);
    vec4 col2 = vshuffle(vzero(), factor, 0, 0, 2, 3);
    vec4 col3 = vshuffle(one, one, 1, 1, 1, 0);
    mat4 m = {{ col0, col1, col2, col3 }};
    return m;
}

static inline mat4 mscale(vec4 factor) __attribute__((always_inline));
static inline mat4 mscale(vec4 factor)
{
    return mscale_zero(vxyz(factor));
}

#endif
