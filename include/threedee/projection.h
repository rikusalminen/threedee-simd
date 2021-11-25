#ifndef THREEDEE_PROJECTION_H
#define THREEDEE_PROJECTION_H

#include <math.h>

#include <threedee/vector.h>
#include <threedee/matrix.h>

static inline mat4 mat_orthov(vec4 mins, vec4 maxs) __attribute__((always_inline));
static inline mat4 mat_orthov(vec4 mins, vec4 maxs)
{
    vec4 rcp = vrcp(maxs - mins);

    vec4 t = vxyz1(-(maxs + mins) * rcp);
    vec4 s = vxyz(vscalar(2.0) * rcp);

    vec4 col0 = vshuffle(s, s, 0, 3, 3, 3);
    vec4 col1 = vshuffle(s, s, 3, 1, 3, 3);
    vec4 col2 = -vshuffle(s, s, 3, 3, 2, 3);
    vec4 col3 = t;

    mat4 m = {{ col0, col1, col2, col3 }};
    return m;
}

static inline mat4 mat_ortho(scalar left, scalar right, scalar bottom, scalar top, scalar near, scalar far) __attribute__((always_inline));
static inline mat4 mat_ortho(scalar left, scalar right, scalar bottom, scalar top, scalar near, scalar far)
{
    return mat_orthov(vec(left, bottom, near, 0.0), vec(right, top, far, 0.0));
}

static inline mat4 mat_frustumv(vec4 mins, vec4 maxs) __attribute__((always_inline));
static inline mat4 mat_frustumv(vec4 mins, vec4 maxs)
{
    vec4 sum = mins + maxs;
    vec4 dif = vrcp(maxs - mins);

    vec4 div = vxyz1(sum * dif);
    vec4 dep = vxyz(vscalar(2.0) * vsplat(mins, 2) * dif);
    vec4 nearxfar = _mm_xor_ps(_mm_set_ss(-0.0),
            smul(vshuffle(dep, dep, 2, 3, 3, 3), vshuffle(maxs, maxs, 2, 0, 0, 0)));

    vec4 col0 = vshuffle(dep, dep, 0, 3, 3, 3);
    vec4 col1 = vshuffle(dep, dep, 3, 1, 3, 3);
    vec4 col2 = _mm_xor_ps(vshuffle(vzero(), _mm_set_ss(-0.0), 0, 0, 0, 0), div);
    vec4 col3 = vshuffle(nearxfar, nearxfar, 3, 3, 0, 3);

    mat4 mat = {{ col0, col1, col2, col3 }};
    return mat;
}

static inline mat4 mat_frustum_scalar(scalar left, scalar right, scalar bottom, scalar top, scalar near, scalar far) __attribute__((always_inline));
static inline mat4 mat_frustum_scalar(scalar left, scalar right, scalar bottom, scalar top, scalar near, scalar far)
{
    mat4 result;

    scalar *matrix = (scalar*)&result;
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

static inline mat4 mat_frustum(scalar left, scalar right, scalar bottom, scalar top, scalar near, scalar far) __attribute__((always_inline));
static inline mat4 mat_frustum(scalar left, scalar right, scalar bottom, scalar top, scalar near, scalar far)
{
    return mat_frustumv(vec(left, bottom, near, 0.0), vec(right, top, far, 0.0));
}

static inline mat4 mat_frustum_inf_zv(vec4 mins, vec4 maxs) __attribute__((always_inline));
static inline mat4 mat_frustum_inf_zv(vec4 mins, vec4 maxs)
{
    vec4 sum = mins + maxs;
    vec4 dif = vrcp(maxs - mins);

    vec4 div = sum * dif;
    vec4 near2 = vxyz(vscalar(2.0) * vsplat(mins, 2));
    vec4 dep = vxyz(near2 * dif);
    near2 = _mm_xor_ps(_mm_set_ss(-0.0), near2);

    vec4 col0 = vshuffle(dep, dep, 0, 3, 3, 3);
    vec4 col1 = vshuffle(dep, dep, 3, 1, 3, 3);
    vec4 col2 = vshuffle(div, _mm_set_ss(-1.0), 0, 1, 0, 0);
    vec4 col3 = vshuffle(near2, near2, 3, 3, 0, 3);

    mat4 mat = {{ col0, col1, col2, col3 }};
    return mat;
}

static inline mat4 mat_frustum_inf_z(scalar left, scalar right, scalar bottom, scalar top, scalar near) __attribute__((always_inline));
static inline mat4 mat_frustum_inf_z(scalar left, scalar right, scalar bottom, scalar top, scalar near)
{
    return mat_frustum_inf_zv(vec(left, bottom, near, 0.0), vec(right, top, 0.0, 0.0));
}

static inline mat4 mat_perspective_fovy(scalar fovy, scalar aspect, scalar near, scalar far) __attribute__((always_inline));
static inline mat4 mat_perspective_fovy(scalar fovy, scalar aspect, scalar near, scalar far)
{
    scalar ymax = near * tanf(fovy / 2.0);
    scalar xmax = ymax * aspect;
    return mat_frustumv(vec(-xmax, -ymax, near, 0.0), vec(xmax, ymax, far, 0.0));
}

static inline mat4 mat_perspective_fovy_scalar(scalar fovy, scalar aspect, scalar near, scalar far) __attribute__((always_inline));
static inline mat4 mat_perspective_fovy_scalar(scalar fovy, scalar aspect, scalar near, scalar far)
{
    scalar ymax = near * tanf(fovy / 2.0);
    scalar xmax = ymax * aspect;
    return mat_frustum_scalar(-xmax, xmax, -ymax, ymax, near, far);
}

// static inline mat4 mat_perspective_fovx(scalar fovx, scalar aspect, scalar near, scalar far);

static inline mat4 mat_perspective_fovy_inf_z(scalar fovy, scalar aspect, scalar near) __attribute__((always_inline));
static inline mat4 mat_perspective_fovy_inf_z(scalar fovy, scalar aspect, scalar near)
{
    scalar ymax = near * tanf(fovy / 2.0);
    scalar xmax = ymax * aspect;
    return mat_frustum_inf_zv(vec(-xmax, -ymax, near, 0.0), vec(xmax, ymax, 0.0, 0.0));
}

// static inline mat4 mat_perspective_fovx_inf_z(scalar fovx, scalar aspect, scalar near, scalar far);

#endif
