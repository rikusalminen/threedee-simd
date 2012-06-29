#ifndef THREEDEE_MATRIX_H
#define THREEDEE_MATRIX_H

#include <threedee/types.h>
#include <threedee/vector.h>

static inline mat4 mtranspose(mat4 m)
{
    _MM_TRANSPOSE4_PS(m.rows[0], m.rows[1], m.rows[2], m.rows[3]);
    return m;
}

static inline mat4 mload(const mat4 *ptr)
{
    mat4 m = {{ vload((vec4*)ptr+0), vload((vec4*)ptr+1), vload((vec4*)ptr+2), vload((vec4*)ptr+3) }};
    return m;
}

static inline mat4 mloadu(const float *ptr)
{
    mat4 m = {{ vloadu((float*)ptr+0), vloadu((float*)ptr+4), vloadu((float*)ptr+8), vloadu((float*)ptr+12) }};
    return m;
}

static inline void mstore(mat4 *ptr, mat4 mat)
{
    vstore((vec4*)ptr + 0, mat.rows[0]);
    vstore((vec4*)ptr + 1, mat.rows[1]);
    vstore((vec4*)ptr + 2, mat.rows[2]);
    vstore((vec4*)ptr + 3, mat.rows[3]);
}

static inline void mstoreu(float *ptr, mat4 mat)
{
    vstoreu(ptr + 0, mat.rows[0]);
    vstoreu(ptr + 4, mat.rows[1]);
    vstoreu(ptr + 8, mat.rows[2]);
    vstoreu(ptr + 12, mat.rows[3]);
}

static inline void mstream(mat4 *ptr, mat4 mat)
{
    vstream((vec4*)ptr + 0, mat.rows[0]);
    vstream((vec4*)ptr + 1, mat.rows[1]);
    vstream((vec4*)ptr + 2, mat.rows[2]);
    vstream((vec4*)ptr + 3, mat.rows[3]);
}

static inline mat4 mloadt(const mat4 *ptr) { return mtranspose(mload(ptr)); }
static inline mat4 mloadut(const float *ptr) { return mtranspose(mloadu(ptr)); }
static inline void mstoret(mat4 *ptr, mat4 mat) { mstore(ptr, mtranspose(mat)); }
static inline void mstoreut(float *ptr, mat4 mat) { mstoreu(ptr, mtranspose(mat)); }
static inline void mstreamt(mat4 *ptr, mat4 mat) { mstream(ptr, mtranspose(mat)); }

static inline mat4 smmul(scalar s, mat4 m)
{
    mat4 result = {{ vscalar(s) * m.rows[0], vscalar(s) * m.rows[1], vscalar(s) * m.rows[2], vscalar(s) * m.rows[3] }};
    return result;
}

static inline mat4 msmul(mat4 m, scalar s) { return smmul(s, m); }

mat4 midentity()
{
    vec4 one = _mm_set_ss(1.0);
    vec4 row0 = vshuffle(one, one, 0, 1, 1, 1);
    vec4 row1 = vshuffle(one, one, 1, 0, 1, 1);
    vec4 row2 = vshuffle(one, one, 1, 1, 0, 1);
    vec4 row3 = vshuffle(one, one, 1, 1, 1, 0);
    mat4 identity = {{ row0, row1, row2, row3 }};
    return identity;
}

static inline mat4 mzero()
{
    mat4 zero = {{ vzero(), vzero(), vzero(), vzero() }};
    return zero;
}

static inline mat4 mat3_to_mat4(mat4 mat)
{
    vec4 one = _mm_set_ss(1.0);
    vec4 row3 = vshuffle(one, one, 1, 1, 1, 0);
    mat4 result = {{ vxyz(mat.rows[0]), vxyz(mat.rows[1]), vxyz(mat.rows[2]), row3 }};
    return result;
}

static inline vec4 mvmul_rows(vec4 x, vec4 y, vec4 z, vec4 w, vec4 v)
{
    return vshuffle(
        vshuffle(
            vdot(x, v),
            vdot(y, v),
            0, 0, 0, 0),
        vshuffle(
            vdot(z, v),
            vdot(w, v),
            0, 0, 0, 0),
        0, 2, 0, 2);
}

static inline vec4 mvmul(mat4 m, vec4 v)
{
    return mvmul_rows(m.rows[0], m.rows[1], m.rows[2], m.rows[3], v);
}

static inline mat4 mmmul(mat4 l, mat4 r)
{
    _MM_TRANSPOSE4_PS(r.rows[0], r.rows[1], r.rows[2], r.rows[3]);

    vec4 row0 = vshuffle(
        vshuffle(vdot(l.rows[0], r.rows[0]), vdot(l.rows[0], r.rows[1]), 0, 0, 0, 0),
        vshuffle(vdot(l.rows[0], r.rows[2]), vdot(l.rows[0], r.rows[3]), 0, 0, 0, 0),
        0, 2, 0, 2);
    vec4 row1 = vshuffle(
        vshuffle(vdot(l.rows[1], r.rows[0]), vdot(l.rows[1], r.rows[1]), 0, 0, 0, 0),
        vshuffle(vdot(l.rows[1], r.rows[2]), vdot(l.rows[1], r.rows[3]), 0, 0, 0, 0),
        0, 2, 0, 2);
    vec4 row2 = vshuffle(
        vshuffle(vdot(l.rows[2], r.rows[0]), vdot(l.rows[2], r.rows[1]), 0, 0, 0, 0),
        vshuffle(vdot(l.rows[2], r.rows[2]), vdot(l.rows[2], r.rows[3]), 0, 0, 0, 0),
        0, 2, 0, 2);
    vec4 row3 = vshuffle(
        vshuffle(vdot(l.rows[3], r.rows[0]), vdot(l.rows[3], r.rows[1]), 0, 0, 0, 0),
        vshuffle(vdot(l.rows[3], r.rows[2]), vdot(l.rows[3], r.rows[3]), 0, 0, 0, 0),
        0, 2, 0, 2);

    mat4 result = { { row0, row1, row2, row3 } };
    return result;
}

static inline mat4 mmmul_freevec(mat4 l, mat4 r)
{
    vec4 row0 =
        vmadd(vsplat(l.rows[0], 3), r.rows[0],
        vmadd(vsplat(l.rows[0], 2), r.rows[0],
        vmadd(vsplat(l.rows[0], 1), r.rows[0],
        (vsplat(l.rows[0], 0) * r.rows[0]))));
    vec4 row1 =
        vmadd(vsplat(l.rows[1], 3), r.rows[1],
        vmadd(vsplat(l.rows[1], 2), r.rows[1],
        vmadd(vsplat(l.rows[1], 1), r.rows[1],
        (vsplat(l.rows[1], 0) * r.rows[1]))));
    vec4 row2 =
        vmadd(vsplat(l.rows[2], 3), r.rows[2],
        vmadd(vsplat(l.rows[2], 2), r.rows[2],
        vmadd(vsplat(l.rows[2], 1), r.rows[2],
        (vsplat(l.rows[2], 0) * r.rows[2]))));
    vec4 row3 =
        vmadd(vsplat(l.rows[3], 3), r.rows[3],
        vmadd(vsplat(l.rows[3], 2), r.rows[3],
        vmadd(vsplat(l.rows[3], 1), r.rows[3],
        (vsplat(l.rows[3], 0) * r.rows[3]))));

    mat4 result = { { row0, row1, row2, row3 } };
    return result;
}

static inline mat4 mmmul_add(mat4 l, mat4 r)
{
    vec4 row0 =
        vsplat(l.rows[0], 3) * r.rows[0] +
        vsplat(l.rows[0], 2) * r.rows[0] +
        vsplat(l.rows[0], 1) * r.rows[0] +
        vsplat(l.rows[0], 0) * r.rows[0];
    vec4 row1 =
        vsplat(l.rows[1], 3) * r.rows[1] +
        vsplat(l.rows[1], 2) * r.rows[1] +
        vsplat(l.rows[1], 1) * r.rows[1] +
        vsplat(l.rows[1], 0) * r.rows[1];
    vec4 row2 =
        vsplat(l.rows[2], 3) * r.rows[2] +
        vsplat(l.rows[2], 2) * r.rows[2] +
        vsplat(l.rows[2], 1) * r.rows[2] +
        vsplat(l.rows[2], 0) * r.rows[2];
    vec4 row3 =
        vsplat(l.rows[3], 3) * r.rows[3] +
        vsplat(l.rows[3], 2) * r.rows[3] +
        vsplat(l.rows[3], 1) * r.rows[3] +
        vsplat(l.rows[3], 0) * r.rows[3];

    mat4 result = { { row0, row1, row2, row3 } };
    return result;
}

static inline mat4 minverse_transpose_rows(vec4 row0, vec4 row1, vec4 row2, vec4 row3)
{
    vec4 minor0, minor1, minor2, minor3;
    vec4 det, tmp1;

    row1 = vshuffle(row1, row1, 2, 3, 0, 1);
    row3 = vshuffle(row3, row3, 2, 3, 0, 1);

    tmp1 = row2 * row3;
    tmp1 = vshuffle(tmp1, tmp1, 1, 0, 3, 2);
    minor0 = row1 * tmp1;
    minor1 = row0 * tmp1;
    tmp1 = vshuffle(tmp1, tmp1, 2, 3, 0, 1);
    minor0 = vmsub(row1, tmp1, minor0);
    minor1 = vmsub(row0, tmp1, minor1);
    minor1 = vshuffle(minor1, minor1, 2, 3, 0, 1);

    tmp1 = row1 * row2;
    tmp1 = vshuffle(tmp1, tmp1, 1, 0, 3, 2);
    minor0 = vmadd(row3, tmp1, minor0);
    minor3 = row0 * tmp1;
    tmp1 = vshuffle(tmp1, tmp1, 2, 3, 0, 1);
    minor0 = vnmadd(row3, tmp1, minor0);
    minor3 = vmsub(row0, tmp1, minor3);
    minor3 = vshuffle(minor3, minor3, 2, 3, 0, 1);

    tmp1 = vshuffle(row1, row1, 2, 3, 0, 1) * row3;
    tmp1 = vshuffle(tmp1, tmp1, 1, 0, 3, 2);
    row2 = vshuffle(row2, row2, 2, 3, 0, 1);
    minor0 = vmadd(row2, tmp1, minor0);
    minor2 = row0 * tmp1;
    tmp1 = vshuffle(tmp1, tmp1, 2, 3, 0, 1);
    minor0 = vnmadd(row2, tmp1, minor0);
    minor2 = vmsub(row0, tmp1, minor2);
    minor2 = vshuffle(minor2, minor2, 2, 3, 0, 1);

    tmp1 = row0 * row1;
    tmp1 = vshuffle(tmp1, tmp1, 1, 0, 3, 2);
    minor2 = vmadd(row3, tmp1, minor2);
    minor3 = vmsub(row2, tmp1, minor3);
    tmp1 = vshuffle(tmp1, tmp1, 2, 3, 0, 1);
    minor2 = vmsub(row3, tmp1, minor2);
    minor3 = vnmadd(row2, tmp1, minor3);

    tmp1 = row0 * row3;
    tmp1 = vshuffle(tmp1, tmp1, 1, 0, 3, 2);
    minor1 = vnmadd(row2, tmp1, minor1);
    minor2 = vmadd(row1, tmp1, minor2);
    tmp1 = vshuffle(tmp1, tmp1, 2, 3, 0, 1);
    minor1 = vmadd(row2, tmp1, minor1);
    minor2 = vnmadd(row1, tmp1, minor2);

    tmp1 = row0 * row2;
    tmp1 = vshuffle(tmp1, tmp1, 1, 0, 3, 2);
    minor1 = vmadd(row3, tmp1, minor1);
    minor3 = vnmadd(row1, tmp1, minor3);
    tmp1 = vshuffle(tmp1, tmp1, 2, 3, 0, 1);
    minor1 = vnmadd(row3, tmp1, minor1);
    minor3 = vmadd(row1, tmp1, minor3);

    det = row0 * minor0;
    det = vshuffle(det, det, 2, 3, 0, 1) + det;
    det = sadd(vshuffle(det, det, 1, 0, 3, 2), det);
    tmp1 = srcp(det);
    det = ssub(sadd(tmp1, tmp1), smul(det, smul(tmp1, tmp1)));
    det = vshuffle(det, det, 0, 0, 0, 0);

    minor0 = det * minor0;
    minor1 = det * minor1;
    minor2 = det * minor2;
    minor3 = det * minor3;

    mat4 result = { { minor0, minor1, minor2, minor3 } };
    return result;
}

static inline mat4 minverse(mat4 m)
{
    return mtranspose(minverse_transpose_rows(m.rows[0], m.rows[1], m.rows[2], m.rows[3]));
}

static inline mat4 minverse_scalar(mat4 mat)
{
    mat4 result;
    float *m = (float*)&mat, *inv = (float*)&result;

    inv[0] =   m[5]*m[10]*m[15] - m[5]*m[11]*m[14] - m[9]*m[6]*m[15]
             + m[9]*m[7]*m[14] + m[13]*m[6]*m[11] - m[13]*m[7]*m[10];
    inv[4] =  -m[4]*m[10]*m[15] + m[4]*m[11]*m[14] + m[8]*m[6]*m[15]
             - m[8]*m[7]*m[14] - m[12]*m[6]*m[11] + m[12]*m[7]*m[10];
    inv[8] =   m[4]*m[9]*m[15] - m[4]*m[11]*m[13] - m[8]*m[5]*m[15]
             + m[8]*m[7]*m[13] + m[12]*m[5]*m[11] - m[12]*m[7]*m[9];
    inv[12] = -m[4]*m[9]*m[14] + m[4]*m[10]*m[13] + m[8]*m[5]*m[14]
             - m[8]*m[6]*m[13] - m[12]*m[5]*m[10] + m[12]*m[6]*m[9];
    inv[1] =  -m[1]*m[10]*m[15] + m[1]*m[11]*m[14] + m[9]*m[2]*m[15]
             - m[9]*m[3]*m[14] - m[13]*m[2]*m[11] + m[13]*m[3]*m[10];
    inv[5] =   m[0]*m[10]*m[15] - m[0]*m[11]*m[14] - m[8]*m[2]*m[15]
             + m[8]*m[3]*m[14] + m[12]*m[2]*m[11] - m[12]*m[3]*m[10];
    inv[9] =  -m[0]*m[9]*m[15] + m[0]*m[11]*m[13] + m[8]*m[1]*m[15]
             - m[8]*m[3]*m[13] - m[12]*m[1]*m[11] + m[12]*m[3]*m[9];
    inv[13] =  m[0]*m[9]*m[14] - m[0]*m[10]*m[13] - m[8]*m[1]*m[14]
             + m[8]*m[2]*m[13] + m[12]*m[1]*m[10] - m[12]*m[2]*m[9];
    inv[2] =   m[1]*m[6]*m[15] - m[1]*m[7]*m[14] - m[5]*m[2]*m[15]
             + m[5]*m[3]*m[14] + m[13]*m[2]*m[7] - m[13]*m[3]*m[6];
    inv[6] =  -m[0]*m[6]*m[15] + m[0]*m[7]*m[14] + m[4]*m[2]*m[15]
             - m[4]*m[3]*m[14] - m[12]*m[2]*m[7] + m[12]*m[3]*m[6];
    inv[10] =  m[0]*m[5]*m[15] - m[0]*m[7]*m[13] - m[4]*m[1]*m[15]
             + m[4]*m[3]*m[13] + m[12]*m[1]*m[7] - m[12]*m[3]*m[5];
    inv[14] = -m[0]*m[5]*m[14] + m[0]*m[6]*m[13] + m[4]*m[1]*m[14]
             - m[4]*m[2]*m[13] - m[12]*m[1]*m[6] + m[12]*m[2]*m[5];
    inv[3] =  -m[1]*m[6]*m[11] + m[1]*m[7]*m[10] + m[5]*m[2]*m[11]
             - m[5]*m[3]*m[10] - m[9]*m[2]*m[7] + m[9]*m[3]*m[6];
    inv[7] =   m[0]*m[6]*m[11] - m[0]*m[7]*m[10] - m[4]*m[2]*m[11]
             + m[4]*m[3]*m[10] + m[8]*m[2]*m[7] - m[8]*m[3]*m[6];
    inv[11] = -m[0]*m[5]*m[11] + m[0]*m[7]*m[9] + m[4]*m[1]*m[11]
             - m[4]*m[3]*m[9] - m[8]*m[1]*m[7] + m[8]*m[3]*m[5];
    inv[15] =  m[0]*m[5]*m[10] - m[0]*m[6]*m[9] - m[4]*m[1]*m[10]
             + m[4]*m[2]*m[9] + m[8]*m[1]*m[6] - m[8]*m[2]*m[5];

    float det = m[0]*inv[0] + m[1]*inv[4] + m[2]*inv[8] + m[3]*inv[12];

    for(int i = 0; i < 16; i++)
        inv[i] = inv[i] * (1.0 / det);

    return result;
}

static inline mat4 minverse3_scalar(mat4 mat)
{
    mat4 result;
    float *m = (float*)&mat, *inv = (float*)&result;

    float inv_det = 1.0 / (
            m[0] * (m[5]*m[10] - m[6]*m[9])
            - m[4] * (m[1]*m[10] - m[2]*m[9])
            + m[8] * (m[1]*m[6] - m[2]*m[5]));

    inv[0] = inv_det * (m[5]*m[10] - m[9]*m[6]);
    inv[1] = inv_det * -(m[1]*m[10] - m[9]*m[2]);
    inv[2] = inv_det * (m[1]*m[6] - m[2]*m[5]);

    inv[4] = inv_det * -(m[4]*m[10] - m[6]*m[8]);
    inv[5] = inv_det * (m[0]*m[10] - m[2]*m[8]);
    inv[6] = inv_det * -(m[4]*m[6] - m[2]*m[4]);

    inv[8] = inv_det * (m[4]*m[9] - m[5]*m[8]);
    inv[9] = inv_det * -(m[0]*m[9] - m[1]*m[8]);
    inv[10] = inv_det * (m[0]*m[5] - m[4]*m[1]);

    inv[3] = inv[7] = inv[11] = 0.0; // 4th row
    inv[12] = inv[13] = inv[14] = 0.0; // 4th column
    inv[15] = 1.0;

    return result;
}

static inline mat4 minverse3(mat4 mat)
{
    return minverse(mat3_to_mat4(mat));
}

#endif
