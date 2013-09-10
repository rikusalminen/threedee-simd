#ifndef THREEDEE_ROTATION_H
#define THREEDEE_ROTATION_H

#include <threedee/vector.h>
#include <threedee/matrix.h>
#include <threedee/sincos.h>

static inline mat4 mat_euler_scalar(vec4 angles) __attribute__((always_inline));
static inline mat4 mat_euler_scalar(vec4 angles)
{
    vec4 vsn, vcs;
    vsincos(angles, &vsn, &vcs);

    float *sn = (float*)&vsn, *cs = (float*)&vcs;

    mat4 mat = {{
        {  cs[1]*cs[2],                     cs[1]*sn[2],                        -sn[1],         0.0 },
        {  sn[0]*sn[1]*cs[2] - cs[0]*sn[2], sn[0]*sn[1]*sn[2] + cs[0]*cs[2],   sn[0]*cs[1],   0.0 },
        {  cs[0]*sn[1]*cs[2] + sn[0]*sn[2],  cs[0]*sn[1]*sn[2] - sn[0]*cs[2],    cs[0]*cs[1],   0.0 },
        { 0.0, 0.0, 0.0, 1.0 }
    }};

    return mat;
}


static inline mat4 mat_euler3(vec4 angles) __attribute__((always_inline));
static inline mat4 mat_euler3(vec4 angles)
{
    vec4 sn, cs;
    vsincos(angles, &sn, &cs);

    vec4 shuf = vshuffle(sn, cs, 0, 2, 0, 2);

    vec4 sam = vshuffle(shuf, cs, 0, 2, 1, 1) * vshuffle(shuf, shuf, 1, 3, 2, 3);

    vec4 neg_sn = -sn;
    vec4 mix = vshuffle(sn, neg_sn, 2, 0, 2, 0) * vshuffle(cs, cs, 1, 1, 0, 2);

    vec4 temp1 = vshuffle(sam, mix, 0, 1, 2, 3);
    vec4 tri =
        (vsplat(sn, 1) * vshuffle(shuf, shuf, 0, 0, 2, 2) * vshuffle(shuf, shuf, 3, 1, 3, 1))
        + vshuffle(temp1, temp1, 2, 1, 0, 3);

    mat4 m;
    m.cols[0] = vshuffle(vshuffle(sam, mix, 3, 3, 0, 0), neg_sn, 0, 2, 1, 1);
    m.cols[1] = vshuffle(tri, mix, 0, 1, 1, 0);
    m.cols[2] = vshuffle(tri, sam, 2, 3, 2, 0);

    return m;
}

static inline mat4 mat_euler(vec4 angles) __attribute__((always_inline));
static inline mat4 mat_euler(vec4 angles)
{
    return mat3_to_mat4(mat_euler3(angles));
}

// static inline mat4 mat_axisangle(vec4 axisangle); __attribute__((always_inline));
// static inline mat4 mat_axisangle(vec4 axisangle);

static inline mat4 quat_to_mat_mmul(vec4 quat) __attribute__((always_inline));
static inline mat4 quat_to_mat_mmul(vec4 quat)
{
    float *q = (float*)&quat;

    mat4 m1 = {{
        { q[3],  q[2], -q[1],  q[0] },
        {-q[2],  q[3],  q[0],  q[1] },
        { q[1], -q[0],  q[3],  q[2] },
        {-q[0], -q[1], -q[2],  q[3] },
    }};

    mat4 m2 = {{
        { q[3],  q[2], -q[1], -q[0] },
        {-q[2],  q[3],  q[0], -q[1] },
        { q[1], -q[0],  q[3], -q[2] },
        { q[0],  q[1],  q[2],  q[3] },
     }};

    return mmmul(m1, m2);
}

static inline mat4 quat_to_mat_scalar(vec4 quat) __attribute__((always_inline));
static inline mat4 quat_to_mat_scalar(vec4 quat)
{
    float *q = (float*)&quat;

    mat4 mat = {{
        { 1.0f - 2.0f*(q[1]*q[1] + q[2]*q[2]), 2.0f*(q[0]*q[1] + q[2]*q[3]), 2.0f*(q[0]*q[2] - q[1]*q[3]), 0.0 },
        { 2.0f*(q[0]*q[1] - q[2]*q[3]), 1.0f - 2.0f*(q[0]*q[0] + q[2]*q[2]), 2.0f*(q[1]*q[2] + q[0]*q[3]), 0.0 },
        { 2.0f*(q[0]*q[2] + q[1]*q[3]), 2.0f*(q[1]*q[2] - q[0]*q[3]), 1.0f - 2.0f*(q[0]*q[0] + q[1]*q[1]), 0.0 },
        { 0.0f, 0.0f, 0.0f, 1.0f }
    }};

    return mat;
}

static inline mat4 quat_to_mat3(vec4 quat) __attribute__((always_inline));
static inline mat4 quat_to_mat3(vec4 quat)
{
    vec4 xs = vsplat(quat, 0) * quat;
    vec4 ys = vsplat(quat, 1) * quat;
    vec4 ws = vsplat(quat, 3) * quat;
    vec4 sq = quat * quat;

    vec4 sum = vscalar(2.0) * (vshuffle(ys, xs, 0, 2, 2, 0) + vshuffle(ws, ws, 2, 0, 1, 0));
    vec4 dif = vscalar(2.0) * (vshuffle(xs, ys, 2, 1, 2, 0) - vshuffle(ws, ws, 1, 2, 0, 0));
    vec4 sqs = vnmadd(vshuffle(sq, sq, 1, 0, 0, 0) + vshuffle(sq, sq, 2, 2, 1, 0), vscalar(2.0), vscalar(1.0));

    mat4 result;
    result.cols[0] = vshuffle(vshuffle(sqs, dif, 0, 0, 1, 1), sum, 0, 2, 2, 2);
    result.cols[1] = vshuffle(vshuffle(sum, sqs, 0, 0, 1, 1), dif, 0, 2, 2, 2);
    result.cols[2] = vshuffle(vshuffle(dif, sum, 0, 0, 1, 1), sqs, 0, 2, 2, 2);

    return result;
}

static inline mat4 quat_to_mat(vec4 quat) __attribute__((always_inline));
static inline mat4 quat_to_mat(vec4 quat)
{
    return mat3_to_mat4(quat_to_mat3(quat));
}

#include <math.h>
static inline vec4 quat_euler_gems(vec4 angles) __attribute__((always_inline));
static inline vec4 quat_euler_gems(vec4 angles)
{
    const float *ptr = (const float *)&angles;
    const float fPitch = ptr[0], fYaw = ptr[1], fRoll = ptr[2];

    const float fSinPitch = sinf(fPitch*0.5F);
    const float fCosPitch = cosf(fPitch*0.5F);
    const float fSinYaw = sinf(fYaw*0.5F);
    const float fCosYaw = cosf(fYaw*0.5F);
    const float fSinRoll = sinf(fRoll*0.5F);
    const float fCosRoll = cosf(fRoll*0.5F);
    const float fCosPitchCosYaw = fCosPitch*fCosYaw;
    const float fSinPitchSinYaw = fSinPitch*fSinYaw;
    float X = fSinRoll * fCosPitchCosYaw     - fCosRoll * fSinPitchSinYaw;
    float Y = fCosRoll * fSinPitch * fCosYaw + fSinRoll * fCosPitch * fSinYaw;
    float Z = fCosRoll * fCosPitch * fSinYaw - fSinRoll * fSinPitch * fCosYaw;
    float W = fCosRoll * fCosPitchCosYaw     + fSinRoll * fSinPitchSinYaw;

    return vec(X, Y, Z, W);
}

static inline vec4 quat_euler_scalar(vec4 angles) __attribute__((always_inline));
static inline vec4 quat_euler_scalar(vec4 angles)
{
    vec4 vsn, vcs;

    vsincos(vscalar(0.5) * angles, &vsn, &vcs);

    float *sn = (float*)&vsn, *cs = (float*)&vcs;

    vec4 result = {
        sn[2] * cs[0] * cs[1] - cs[2] * sn[0] * sn[1],
        cs[2] * sn[0] * cs[1] + sn[2] * cs[0] * sn[1],
        cs[2] * cs[0] * sn[1] - sn[2] * sn[0] * cs[1],
        cs[2] * cs[0] * cs[1] + sn[2] * sn[0] * sn[1]
    };
    return result;
}

static inline vec4 quat_euler(vec4 angles) __attribute__((always_inline));
static inline vec4 quat_euler(vec4 angles)
{
    vec4 sn, cs;
    vsincos(vscalar(0.5) * angles, &sn, &cs);

    vec4 sc = vshuffle(sn, cs, 1, 2, 2, 1);

    vec4 l =
        vshuffle(sn, sc, 2, 0, 0, 2) *
        vshuffle(cs, cs, 0, 2, 2, 0) *
        vshuffle(cs, cs, 1, 1, 0, 1);

    vec4 r =
        vshuffle(cs, sc, 2, 0, 3, 1) *
        vshuffle(sn, sn, 0, 2, 2, 0) *
        vshuffle(sn, sn, 1, 1, 0, 1);

#ifdef xor_negate
    vec4 negative = _mm_set1_ps(-0.0);
    vec4 result = l + _mm_xor_ps(vshuffle(negative, negative, 0, 1, 0, 1), r);
#else
    vec4 sum = l + r, dif = l - r;
    vec4 sumdif = vshuffle(sum, dif, 0, 2, 1, 3);
    vec4 result = vshuffle(sumdif, sumdif, 0, 2, 1, 3);
#endif

    return result;
}

/**
 * convert an axis-angle vector to a quaternion
 * angle in radians
 *
 * [axis.x, axis.y, axis.z, angle] => quaternion
 */
static inline vec4 quat_axisangle(vec4 axisangle) __attribute__((always_inline));
static inline vec4 quat_axisangle(vec4 axisangle) {
    const float halfAngle = axisangle[3] * 0.5f;

    /* overwrite the last entry */
    vec4 result = axisangle * vscalar(sinf(halfAngle));
    result[3] = cosf(halfAngle);

    return result;
}

// static inline vec4 quat_to_axisangle(vec4 quat); __attribute__((always_inline));
// static inline vec4 quat_to_axisangle(vec4 quat);
// static inline vec4 mat_to_quat(mat4 mat); __attribute__((always_inline));
// static inline vec4 mat_to_quat(mat4 mat);

static inline vec4 qprod_scalar(vec4 x, vec4 y) __attribute__((always_inline));
static inline vec4 qprod_scalar(vec4 x, vec4 y)
{
    float *xs = (float*)&x, *ys = (float*)&y;

    vec4 result = {
        xs[3]*ys[0] + xs[0]*ys[3] + xs[1]*ys[2] - xs[2]*ys[1],
        xs[3]*ys[1] + xs[1]*ys[3] + xs[2]*ys[0] - xs[0]*ys[2],
        xs[3]*ys[2] + xs[2]*ys[3] + xs[0]*ys[1] - xs[1]*ys[0],
        xs[3]*ys[3] - xs[0]*ys[0] - xs[1]*ys[1] - xs[2]*ys[2],
    };
    return result;
}

static inline vec4 qprod(vec4 x, vec4 y) __attribute__((always_inline));
static inline vec4 qprod(vec4 x, vec4 y)
{
    vec4 negative = _mm_set_ss(-0.0);
    return
    vshuffle(x, x, 3, 3, 3, 3) * y +
        _mm_xor_ps(vshuffle(negative, negative, 1, 1, 1, 0),
            vshuffle(x, x, 0, 1, 2, 0) * vshuffle(y, y, 3, 3, 3, 0) +
            vshuffle(x, x, 1, 2, 0, 1) * vshuffle(y, y, 2, 0, 1, 1)) -
        vshuffle(x, x, 2, 0, 1, 2) * vshuffle(y, y, 1, 2, 0, 2);
}

static inline vec4 qprod_mad(vec4 x, vec4 y) __attribute__((always_inline));
static inline vec4 qprod_mad(vec4 x, vec4 y)
{
    vec4 negative = _mm_set_ss(-0.0);
    return vnmadd(
        vshuffle(x, x, 2, 0, 1, 2), vshuffle(y, y, 1, 2, 0, 2),
            vmadd(vshuffle(x, x, 3, 3, 3, 3), y,
            _mm_xor_ps(vshuffle(negative, negative, 1, 1, 1, 0),
                vmadd(vshuffle(x, x, 0, 1, 2, 0), vshuffle(y, y, 3, 3, 3, 0),
                vshuffle(x, x, 1, 2, 0, 1) * vshuffle(y, y, 2, 0, 1, 1)))));
}

static inline vec4 qconj(vec4 x) __attribute__((always_inline));
static inline vec4 qconj(vec4 x)
{
    vec4 negative = _mm_set_ss(-0.0);
    return _mm_xor_ps(vshuffle(negative, negative, 0, 0, 0, 1), x);
}

#endif
