#ifndef THREEDEE_MATRIX_H
#define THREEDEE_MATRIX_H

#include <threedee/types.h>
#include <threedee/vector.h>

/* Transpose the 4x4 matrix composed of row[0-3].  */
#define _MM_TRANSPOSE4_PD(row0, row1, row2, row3)			\
do {									\
  __v4sd __r0 = (row0), __r1 = (row1), __r2 = (row2), __r3 = (row3);	\
  __v4sd __t0 = __builtin_ia32_unpcklpd (__r0, __r1);			\
  __v4sd __t1 = __builtin_ia32_unpcklpd (__r2, __r3);			\
  __v4sd __t2 = __builtin_ia32_unpckhpd (__r0, __r1);			\
  __v4sd __t3 = __builtin_ia32_unpckhpd (__r2, __r3);			\
  (row0) = __builtin_ia32_movlhpd (__t0, __t1);				\
  (row1) = __builtin_ia32_movhlpd (__t1, __t0);				\
  (row2) = __builtin_ia32_movlhpd (__t2, __t3);				\
  (row3) = __builtin_ia32_movhlpd (__t3, __t2);				\
} while (0)

#define _MM_TRANSPOSE4_PI(row0, row1, row2, row3)			\
do {									\
  __v4si __r0 = (row0), __r1 = (row1), __r2 = (row2), __r3 = (row3);	\
  __v4si __t0 = __builtin_ia32_unpcklpi (__r0, __r1);			\
  __v4si __t1 = __builtin_ia32_unpcklpi (__r2, __r3);			\
  __v4si __t2 = __builtin_ia32_unpckhpi (__r0, __r1);			\
  __v4si __t3 = __builtin_ia32_unpckhpi (__r2, __r3);			\
  (row0) = __builtin_ia32_movlhpi (__t0, __t1);				\
  (row1) = __builtin_ia32_movhlpi (__t1, __t0);				\
  (row2) = __builtin_ia32_movlhpi (__t2, __t3);				\
  (row3) = __builtin_ia32_movhlpi (__t3, __t2);				\
} while (0)

static inline mat4 mtranspose(mat4 m) __attribute__((always_inline));
static inline mat4 mtranspose(mat4 m)
{
    _MM_TRANSPOSE4_PS(m.cols[0], m.cols[1], m.cols[2], m.cols[3]);
    return m;
}

static inline mat4 mload(const mat4 *ptr) __attribute__((always_inline));
static inline mat4 mload(const mat4 *ptr)
{
    mat4 m = {{ vload((vec4*)ptr+0), vload((vec4*)ptr+1), vload((vec4*)ptr+2), vload((vec4*)ptr+3) }};
    return m;
}

static inline mat4 mloadu(const scalar *ptr) __attribute__((always_inline));
static inline mat4 mloadu(const scalar *ptr)
{
    mat4 m = {{ vloadu((scalar*)ptr+0), vloadu((scalar*)ptr+4), vloadu((scalar*)ptr+8), vloadu((scalar*)ptr+12) }};
    return m;
}

static inline void mstore(mat4 *ptr, mat4 mat) __attribute__((always_inline));
static inline void mstore(mat4 *ptr, mat4 mat)
{
    vstore((vec4*)ptr + 0, mat.cols[0]);
    vstore((vec4*)ptr + 1, mat.cols[1]);
    vstore((vec4*)ptr + 2, mat.cols[2]);
    vstore((vec4*)ptr + 3, mat.cols[3]);
}

static inline void mstoreu(scalar *ptr, mat4 mat) __attribute__((always_inline));
static inline void mstoreu(scalar *ptr, mat4 mat)
{
    vstoreu(ptr + 0, mat.cols[0]);
    vstoreu(ptr + 4, mat.cols[1]);
    vstoreu(ptr + 8, mat.cols[2]);
    vstoreu(ptr + 12, mat.cols[3]);
}

static inline void mstream(mat4 *ptr, mat4 mat) __attribute__((always_inline));
static inline void mstream(mat4 *ptr, mat4 mat)
{
    vstream((vec4*)ptr + 0, mat.cols[0]);
    vstream((vec4*)ptr + 1, mat.cols[1]);
    vstream((vec4*)ptr + 2, mat.cols[2]);
    vstream((vec4*)ptr + 3, mat.cols[3]);
}

static inline mat4 mloadt(const mat4 *ptr) __attribute__((always_inline));
static inline mat4 mloadt(const mat4 *ptr) { return mtranspose(mload(ptr)); }
static inline mat4 mloadut(const scalar *ptr) __attribute__((always_inline));
static inline mat4 mloadut(const scalar *ptr) { return mtranspose(mloadu(ptr)); }
static inline void mstoret(mat4 *ptr, mat4 mat) __attribute__((always_inline));
static inline void mstoret(mat4 *ptr, mat4 mat) { mstore(ptr, mtranspose(mat)); }
static inline void mstoreut(scalar *ptr, mat4 mat) __attribute__((always_inline));
static inline void mstoreut(scalar *ptr, mat4 mat) { mstoreu(ptr, mtranspose(mat)); }
static inline void mstreamt(mat4 *ptr, mat4 mat) __attribute__((always_inline));
static inline void mstreamt(mat4 *ptr, mat4 mat) { mstream(ptr, mtranspose(mat)); }

static inline mat4 smmul(scalar s, mat4 m) __attribute__((always_inline));
static inline mat4 smmul(scalar s, mat4 m)
{
    mat4 result = {{ vscalar(s) * m.cols[0], vscalar(s) * m.cols[1], vscalar(s) * m.cols[2], vscalar(s) * m.cols[3] }};
    return result;
}

static inline mat4 msmul(mat4 m, scalar s) __attribute__((always_inline));
static inline mat4 msmul(mat4 m, scalar s) { return smmul(s, m); }

static inline mat4 midentity() __attribute__((always_inline));
static inline mat4 midentity()
{
    vec4 one = _mm_set_ss(1.0);
    vec4 row0 = vshuffle(one, one, 0, 1, 1, 1);
    vec4 row1 = vshuffle(one, one, 1, 0, 1, 1);
    vec4 row2 = vshuffle(one, one, 1, 1, 0, 1);
    vec4 row3 = vshuffle(one, one, 1, 1, 1, 0);
    mat4 identity = {{ row0, row1, row2, row3 }};
    return identity;
}

static inline mat4 mzero() __attribute__((always_inline));
static inline mat4 mzero()
{
    mat4 zero = {{ vzero(), vzero(), vzero(), vzero() }};
    return zero;
}

static inline mat4 mat3_to_mat4(mat4 mat) __attribute__((always_inline));
static inline mat4 mat3_to_mat4(mat4 mat)
{
    vec4 one = _mm_set_ss(1.0);
    vec4 row3 = vshuffle(one, one, 1, 1, 1, 0);
    mat4 result = {{ vxyz(mat.cols[0]), vxyz(mat.cols[1]), vxyz(mat.cols[2]), row3 }};
    return result;
}

static inline vec4 mvmul_add_cols(vec4 col0, vec4 col1, vec4 col2, vec4 col3, vec4 v) __attribute__((always_inline));
static inline vec4 mvmul_add_cols(vec4 col0, vec4 col1, vec4 col2, vec4 col3, vec4 v)
{
    return col0 * vsplat(v, 0) + col1 * vsplat(v, 1) + col2 * vsplat(v, 2) + col3 * vsplat(v, 3);
}

static inline vec4 mvmul_madd_cols(vec4 col0, vec4 col1, vec4 col2, vec4 col3, vec4 v)
{
    return vmadd(col0, vsplat(v, 0),
            vmadd(col1, vsplat(v, 1),
            vmadd(col2, vsplat(v, 2),
            vmul(col3, vsplat(v, 3)))));
}

static inline vec4 mvmul_dot_rows(vec4 x, vec4 y, vec4 z, vec4 w, vec4 v) __attribute__((always_inline));
static inline vec4 mvmul_dot_rows(vec4 x, vec4 y, vec4 z, vec4 w, vec4 v)
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

static inline vec4 mvmul(mat4 m, vec4 v) __attribute__((always_inline));
static inline vec4 mvmul(mat4 m, vec4 v)
{
    mtranspose(m);
    return mvmul_dot_rows(m.cols[0], m.cols[1], m.cols[2], m.cols[3], v);
}

static inline mat4 mmmul_dot(mat4 l, mat4 r) __attribute__((always_inline));
static inline mat4 mmmul_dot(mat4 l, mat4 r)
{
    mtranspose(l);

    vec4 row0 = vshuffle(
        vshuffle(vdot(l.cols[0], r.cols[0]), vdot(l.cols[0], r.cols[1]), 0, 0, 0, 0),
        vshuffle(vdot(l.cols[0], r.cols[2]), vdot(l.cols[0], r.cols[3]), 0, 0, 0, 0),
        0, 2, 0, 2);
    vec4 row1 = vshuffle(
        vshuffle(vdot(l.cols[1], r.cols[0]), vdot(l.cols[1], r.cols[1]), 0, 0, 0, 0),
        vshuffle(vdot(l.cols[1], r.cols[2]), vdot(l.cols[1], r.cols[3]), 0, 0, 0, 0),
        0, 2, 0, 2);
    vec4 row2 = vshuffle(
        vshuffle(vdot(l.cols[2], r.cols[0]), vdot(l.cols[2], r.cols[1]), 0, 0, 0, 0),
        vshuffle(vdot(l.cols[2], r.cols[2]), vdot(l.cols[2], r.cols[3]), 0, 0, 0, 0),
        0, 2, 0, 2);
    vec4 row3 = vshuffle(
        vshuffle(vdot(l.cols[3], r.cols[0]), vdot(l.cols[3], r.cols[1]), 0, 0, 0, 0),
        vshuffle(vdot(l.cols[3], r.cols[2]), vdot(l.cols[3], r.cols[3]), 0, 0, 0, 0),
        0, 2, 0, 2);

    mat4 result = { { row0, row1, row2, row3 } };
    return result;
}

static inline mat4 mmmul_madd(mat4 l, mat4 r) __attribute__((always_inline));
static inline mat4 mmmul_madd(mat4 l, mat4 r)
{
    vec4 col0 =
        vmadd(vsplat(r.cols[0], 3), l.cols[3],
        vmadd(vsplat(r.cols[0], 2), l.cols[2],
        vmadd(vsplat(r.cols[0], 1), l.cols[1],
        vmadd(vsplat(r.cols[0], 0), l.cols[0],
            vzero()))));
    vec4 col1 =
        vmadd(vsplat(r.cols[1], 3), l.cols[3],
        vmadd(vsplat(r.cols[1], 2), l.cols[2],
        vmadd(vsplat(r.cols[1], 1), l.cols[1],
        vmadd(vsplat(r.cols[1], 0), l.cols[0],
            vzero()))));
    vec4 col2 =
        vmadd(vsplat(r.cols[2], 3), l.cols[3],
        vmadd(vsplat(r.cols[2], 2), l.cols[2],
        vmadd(vsplat(r.cols[2], 1), l.cols[1],
        vmadd(vsplat(r.cols[2], 0), l.cols[0],
            vzero()))));
    vec4 col3 =
        vmadd(vsplat(r.cols[3], 3), l.cols[3],
        vmadd(vsplat(r.cols[3], 2), l.cols[2],
        vmadd(vsplat(r.cols[3], 1), l.cols[1],
        vmadd(vsplat(r.cols[3], 0), l.cols[0],
            vzero()))));

    mat4 result = { { col0, col1, col2, col3 } };
    return result;
}

static inline mat4 mmmul_add(mat4 l, mat4 r) __attribute__((always_inline));
static inline mat4 mmmul_add(mat4 l, mat4 r)
{
    vec4 col0 =
        vsplat(r.cols[0], 3) * l.cols[3] +
        vsplat(r.cols[0], 2) * l.cols[2] +
        vsplat(r.cols[0], 1) * l.cols[1] +
        vsplat(r.cols[0], 0) * l.cols[0];
    vec4 col1 =
        vsplat(r.cols[1], 3) * l.cols[3] +
        vsplat(r.cols[1], 2) * l.cols[2] +
        vsplat(r.cols[1], 1) * l.cols[1] +
        vsplat(r.cols[1], 0) * l.cols[0];
    vec4 col2 =
        vsplat(r.cols[2], 3) * l.cols[3] +
        vsplat(r.cols[2], 2) * l.cols[2] +
        vsplat(r.cols[2], 1) * l.cols[1] +
        vsplat(r.cols[2], 0) * l.cols[0];
    vec4 col3 =
        vsplat(r.cols[3], 3) * l.cols[3] +
        vsplat(r.cols[3], 2) * l.cols[2] +
        vsplat(r.cols[3], 1) * l.cols[1] +
        vsplat(r.cols[3], 0) * l.cols[0];

    mat4 result = { { col0, col1, col2, col3 } };
    return result;
}

static inline mat4 mmmul(mat4 a, mat4 b) { return mmmul_add(a, b); }

static inline mat4 minverse_transpose_cols(vec4 col0, vec4 col1, vec4 col2, vec4 col3) __attribute__((always_inline));
static inline mat4 minverse_transpose_cols(vec4 col0, vec4 col1, vec4 col2, vec4 col3)
{
    vec4 minor0, minor1, minor2, minor3;
    vec4 det, tmp1;

    col1 = vshuffle(col1, col1, 2, 3, 0, 1);
    col3 = vshuffle(col3, col3, 2, 3, 0, 1);

    tmp1 = col2 * col3;
    tmp1 = vshuffle(tmp1, tmp1, 1, 0, 3, 2);
    minor0 = col1 * tmp1;
    minor1 = col0 * tmp1;
    tmp1 = vshuffle(tmp1, tmp1, 2, 3, 0, 1);
    minor0 = vmsub(col1, tmp1, minor0);
    minor1 = vmsub(col0, tmp1, minor1);
    minor1 = vshuffle(minor1, minor1, 2, 3, 0, 1);

    tmp1 = col1 * col2;
    tmp1 = vshuffle(tmp1, tmp1, 1, 0, 3, 2);
    minor0 = vmadd(col3, tmp1, minor0);
    minor3 = col0 * tmp1;
    tmp1 = vshuffle(tmp1, tmp1, 2, 3, 0, 1);
    minor0 = vnmadd(col3, tmp1, minor0);
    minor3 = vmsub(col0, tmp1, minor3);
    minor3 = vshuffle(minor3, minor3, 2, 3, 0, 1);

    tmp1 = vshuffle(col1, col1, 2, 3, 0, 1) * col3;
    tmp1 = vshuffle(tmp1, tmp1, 1, 0, 3, 2);
    col2 = vshuffle(col2, col2, 2, 3, 0, 1);
    minor0 = vmadd(col2, tmp1, minor0);
    minor2 = col0 * tmp1;
    tmp1 = vshuffle(tmp1, tmp1, 2, 3, 0, 1);
    minor0 = vnmadd(col2, tmp1, minor0);
    minor2 = vmsub(col0, tmp1, minor2);
    minor2 = vshuffle(minor2, minor2, 2, 3, 0, 1);

    tmp1 = col0 * col1;
    tmp1 = vshuffle(tmp1, tmp1, 1, 0, 3, 2);
    minor2 = vmadd(col3, tmp1, minor2);
    minor3 = vmsub(col2, tmp1, minor3);
    tmp1 = vshuffle(tmp1, tmp1, 2, 3, 0, 1);
    minor2 = vmsub(col3, tmp1, minor2);
    minor3 = vnmadd(col2, tmp1, minor3);

    tmp1 = col0 * col3;
    tmp1 = vshuffle(tmp1, tmp1, 1, 0, 3, 2);
    minor1 = vnmadd(col2, tmp1, minor1);
    minor2 = vmadd(col1, tmp1, minor2);
    tmp1 = vshuffle(tmp1, tmp1, 2, 3, 0, 1);
    minor1 = vmadd(col2, tmp1, minor1);
    minor2 = vnmadd(col1, tmp1, minor2);

    tmp1 = col0 * col2;
    tmp1 = vshuffle(tmp1, tmp1, 1, 0, 3, 2);
    minor1 = vmadd(col3, tmp1, minor1);
    minor3 = vnmadd(col1, tmp1, minor3);
    tmp1 = vshuffle(tmp1, tmp1, 2, 3, 0, 1);
    minor1 = vnmadd(col3, tmp1, minor1);
    minor3 = vmadd(col1, tmp1, minor3);

    det = col0 * minor0;
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

static inline mat4 minverse_transpose(mat4 m) __attribute__((always_inline));
static inline mat4 minverse_transpose(mat4 m)
{
    return minverse_transpose_cols(m.cols[0], m.cols[1], m.cols[2], m.cols[3]);
}

static inline mat4 minverse(mat4 m) __attribute__((always_inline));
static inline mat4 minverse(mat4 m)
{
    return mtranspose(minverse_transpose_cols(m.cols[0], m.cols[1], m.cols[2], m.cols[3]));
}

static inline mat4 minverse_scalar(mat4 mat) __attribute__((always_inline));
static inline mat4 minverse_scalar(mat4 mat)
{
    mat4 result;
    scalar *m = (scalar*)&mat, *inv = (scalar*)&result;

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

    scalar det = m[0]*inv[0] + m[1]*inv[4] + m[2]*inv[8] + m[3]*inv[12];

    for(int i = 0; i < 16; i++)
        inv[i] = inv[i] * (1.0 / det);

    return result;
}

static inline mat4 minverse3_scalar(mat4 mat)
{
    mat4 result;
    scalar *m = (scalar*)&mat, *inv = (scalar*)&result;

    scalar inv_det = 1.0 / (
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
