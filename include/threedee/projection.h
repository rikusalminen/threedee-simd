#ifndef THREEDEE_PROJECTION_H
#define THREEDEE_PROJECTION_H

#include <math.h>

#include <threedee/vector.h>
#include <threedee/matrix.h>

static inline mat4 mat_orthov_zero(vec4 mins, vec4 maxs)
{
    vec4 rcp = vrcp(maxs - mins);
    rcp = vshuffle(rcp, vshuffle(rcp, vzero(), 2, 3, 0, 0), 0, 1, 0, 3);

    vec4 t = -(maxs + mins) * rcp;
    vec4 s = vscalar(2.0) * rcp;
    vec4 one = _mm_set_ss(1.0);

    vec4 row0 = vshuffle(s, t, 0, 3, 3, 0);
    vec4 row1 = vshuffle(s, t, 3, 0, 3, 0);
    vec4 row2 = vshuffle(vzero(), vshuffle(s, t, 2, 2, 2, 2), 0, 0, 0, 2);
    vec4 row3 = vshuffle(one, one, 1, 1, 1, 0);

    mat4 mat = {{ row0, row1, row2, row3 }};
    return mat;
}

static inline mat4 mat_orthov(vec4 mins, vec4 maxs)
{
    return mat_orthov_zero(vxyz(mins), vxyz(maxs));
}

static inline mat4 mat_ortho(scalar left, scalar right, scalar bottom, scalar top, scalar near, scalar far)
{
    return mat_orthov_zero(vec(left, bottom, near, 0.0), vec(right, top, far, 0.0));
}

static inline mat4 mat_frustumv(vec4 mins, vec4 maxs)
{
    vec4 sum = mins + maxs;
    vec4 dif = vrcp(maxs - mins);

    vec4 div = sum * dif;
    vec4 dep = vscalar(2.0) * vsplat(mins, 2) * dif;

    vec4 shuf = vshuffle(div, dep, 0, 1, 0, 1);

    vec4 shuf0 = vshuffle(shuf, vzero(), 0, 2, 0, 0);
    vec4 row0 = vshuffle(shuf0, shuf0, 1, 3, 0, 3);

    vec4 shuf1 = vshuffle(shuf, vzero(), 1, 3, 0, 0);
    vec4 row1 = vshuffle(shuf1, shuf1, 3, 1, 0, 3);

    vec4 nearxfar = smul(vshuffle(maxs, maxs, 2, 0, 0, 0), vshuffle(dep, dep, 2, 0, 0, 0));
    vec4 row2 = _mm_xor_ps(vshuffle(vzero(), _mm_set_ss(-0.0), 0, 0, 0, 0),
            vshuffle(vzero(), vshuffle(div, nearxfar, 2, 2, 0, 0), 0, 0, 0, 2));

    vec4 row3 = vshuffle(vzero(), _mm_set_ss(-1.0), 0, 0, 0, 1);

    mat4 mat = {{ row0, row1, row2, row3 }};
    return mat;
}

static inline mat4 mat_frustum_scalar(scalar left, scalar right, scalar bottom, scalar top, scalar near, scalar far)
{
    mat4 result;

    float *matrix = (float*)&result;
    matrix[0] = 2.0 * near / (right - left);
    matrix[1] = 0.0;
    matrix[2] = (right + left) / (right - left);
    matrix[3] = 0.0;

    matrix[4] = 0.0;
    matrix[5] = 2.0 * near / (top - bottom);
    matrix[6] = (top + bottom) / (top - bottom);
    matrix[7] = 0.0;

    matrix[8] = 0.0;
    matrix[9] = 0.0;
    matrix[10] = -(far + near) / (far - near);
    matrix[11] = (-2.0 * near * far) / (far - near);

    matrix[12] = 0.0;
    matrix[13] = 0.0;
    matrix[14] = -1.0;
    matrix[15] = 0.0;

    return result;
}

static inline mat4 mat_frustum(scalar left, scalar right, scalar bottom, scalar top, scalar near, scalar far)
{
    return mat_frustumv(vec(left, bottom, near, 0.0), vec(right, top, far, 0.0));
}

static inline mat4 mat_frustum_inf_zv(vec4 mins, vec4 maxs)
{
    vec4 sum = mins + maxs;
    vec4 dif = vrcp(maxs - mins);

    vec4 div = sum * dif;
    vec4 near2 = vscalar(2.0) * vsplat(mins, 2);
    vec4 dep = near2 * dif;

    vec4 shuf = vshuffle(div, dep, 0, 1, 0, 1);

    vec4 shuf0 = vshuffle(shuf, vzero(), 0, 2, 0, 0);
    vec4 row0 = vshuffle(shuf0, shuf0, 1, 3, 0, 3);

    vec4 shuf1 = vshuffle(shuf, vzero(), 1, 3, 0, 0);
    vec4 row1 = vshuffle(shuf1, shuf1, 3, 1, 0, 3);

    vec4 row2 = _mm_xor_ps(vshuffle(vzero(), _mm_set_ss(-0.0), 0, 0, 0, 0),
            vshuffle(vzero(), vshuffle(_mm_set_ss(1.0), near2, 0, 0, 0, 0), 0, 0, 0, 2));

    vec4 row3 = vshuffle(vzero(), _mm_set_ss(-1.0), 0, 0, 0, 1);

    mat4 mat = {{ row0, row1, row2, row3 }};
    return mat;
}

static inline mat4 mat_frustum_inf_z(scalar left, scalar right, scalar bottom, scalar top, scalar near)
{
    return mat_frustum_inf_zv(vec(left, bottom, near, 0.0), vec(right, top, 0.0, 0.0));
}

static inline mat4 mat_perspective_fovy(scalar fovy, scalar aspect, scalar near, scalar far)
{
    float ymax = near * tanf(fovy / 2.0);
    float xmax = ymax * aspect;
    return mat_frustumv(vec(-xmax, -ymax, near, 0.0), vec(xmax, ymax, far, 0.0));
}

static inline mat4 mat_perspective_fovy_scalar(scalar fovy, scalar aspect, scalar near, scalar far)
{
    float ymax = near * tanf(fovy / 2.0);
    float xmax = ymax * aspect;
    return mat_frustum_scalar(-xmax, xmax, -ymax, ymax, near, far);
}

// static inline mat4 mat_perspective_fovx(scalar fovx, scalar aspect, scalar near, scalar far);

static inline mat4 mat_perspective_fovy_inf_z(scalar fovy, scalar aspect, scalar near)
{
    float ymax = near * tanf(fovy / 2.0);
    float xmax = ymax * aspect;
    return mat_frustum_inf_zv(vec(-xmax, -ymax, near, 0.0), vec(xmax, ymax, 0.0, 0.0));
}

// static inline mat4 mat_perspective_fovx_inf_z(scalar fovx, scalar aspect, scalar near, scalar far);

#endif
