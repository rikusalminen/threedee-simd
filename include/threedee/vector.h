#ifndef THREEDEE_VECTOR_H
#define THREEDEE_VECTOR_H

#include <stdint.h>

#include <threedee/types.h>

#include <x86intrin.h>

static inline vec4 vec(scalar x, scalar y, scalar z, scalar w) { return _mm_setr_ps(x, y, z, w); }
static inline vec4 vscalar(scalar x) __attribute__((always_inline));
static inline vec4 vscalar(scalar x) { return _mm_set1_ps(x); }
static inline vec4 vscalari(int32_t x) __attribute__((always_inline));
static inline vec4 vscalari(int32_t x) { union { float f; int32_t i; } u = { .i = x }; return _mm_set1_ps(u.f); }
static inline vec4 vscalaru(uint32_t x) __attribute__((always_inline));
static inline vec4 vscalaru(uint32_t x) { union { float f; uint32_t ui; } u = { .ui = x }; return _mm_set1_ps(u.f); }
static inline vec4 vzero(void) __attribute__((always_inline));
static inline vec4 vzero(void) { return _mm_setzero_ps(); }

static inline vec4 vadd(vec4 x, vec4 y) __attribute__((always_inline));
static inline vec4 vadd(vec4 x, vec4 y) { return _mm_add_ps(x, y); }
static inline vec4 vsub(vec4 x, vec4 y) __attribute__((always_inline));
static inline vec4 vsub(vec4 x, vec4 y) { return _mm_sub_ps(x, y); }
static inline vec4 vmul(vec4 x, vec4 y) __attribute__((always_inline));
static inline vec4 vmul(vec4 x, vec4 y) { return _mm_mul_ps(x, y); }
static inline vec4 vdiv(vec4 x, vec4 y) __attribute__((always_inline));
static inline vec4 vdiv(vec4 x, vec4 y) { return _mm_div_ps(x, y); }
static inline vec4 vmin(vec4 x, vec4 y) __attribute__((always_inline));
static inline vec4 vmin(vec4 x, vec4 y) { return _mm_min_ps(x, y); }
static inline vec4 vmax(vec4 x, vec4 y) __attribute__((always_inline));
static inline vec4 vmax(vec4 x, vec4 y) { return _mm_max_ps(x, y); }

static inline vec4 vsqrt(vec4 x) __attribute__((always_inline));
static inline vec4 vsqrt(vec4 x) { return _mm_sqrt_ps(x); }
static inline vec4 vrcp(vec4 x) __attribute__((always_inline));
static inline vec4 vrcp(vec4 x) { return _mm_rcp_ps(x); }
static inline vec4 vrsqrt(vec4 x) __attribute__((always_inline));
static inline vec4 vrsqrt(vec4 x) { return _mm_rsqrt_ps(x); }

static inline vec4 sadd(vec4 x, vec4 y) __attribute__((always_inline));
static inline vec4 sadd(vec4 x, vec4 y) { return _mm_add_ss(x, y); }
static inline vec4 ssub(vec4 x, vec4 y) __attribute__((always_inline));
static inline vec4 ssub(vec4 x, vec4 y) { return _mm_sub_ss(x, y); }
static inline vec4 smul(vec4 x, vec4 y) __attribute__((always_inline));
static inline vec4 smul(vec4 x, vec4 y) { return _mm_mul_ss(x, y); }
static inline vec4 sdiv(vec4 x, vec4 y) __attribute__((always_inline));
static inline vec4 sdiv(vec4 x, vec4 y) { return _mm_div_ss(x, y); }
static inline vec4 smin(vec4 x, vec4 y) __attribute__((always_inline));
static inline vec4 smin(vec4 x, vec4 y) { return _mm_min_ss(x, y); }
static inline vec4 smax(vec4 x, vec4 y) __attribute__((always_inline));
static inline vec4 smax(vec4 x, vec4 y) { return _mm_max_ss(x, y); }

static inline vec4 ssqrt(vec4 x) __attribute__((always_inline));
static inline vec4 ssqrt(vec4 x) { return _mm_sqrt_ss(x); }
static inline vec4 srcp(vec4 x) __attribute__((always_inline));
static inline vec4 srcp(vec4 x) { return _mm_rcp_ss(x); }
static inline vec4 srsqrt(vec4 x) __attribute__((always_inline));
static inline vec4 srsqrt(vec4 x) { return _mm_rsqrt_ss(x); }

static inline vec4 vload(const vec4 *ptr) __attribute__((always_inline));
static inline vec4 vload(const vec4 *ptr) { return _mm_load_ps((float*)ptr); }
static inline vec4 vloadu(const float *ptr) __attribute__((always_inline));
static inline vec4 vloadu(const float *ptr) { return _mm_loadu_ps((float*)ptr); }

static inline void vstore(vec4 *ptr, vec4 v) __attribute__((always_inline));
static inline void vstore(vec4 *ptr, vec4 v) { _mm_store_ps((float*)ptr, v); }
static inline void vstoreu(float *ptr, vec4 v) __attribute__((always_inline));
static inline void vstoreu(float *ptr, vec4 v) { _mm_storeu_ps(ptr, v); }

static inline void vstream(vec4 *ptr, vec4 v) __attribute__((always_inline));
static inline void vstream(vec4 *ptr, vec4 v) { _mm_stream_ps((float*)ptr, v); }

#ifdef __clang__
#define vshuffle(x, y, a, b, c, d) (__builtin_shufflevector((x), (y), (a), (b), (4+c), (4+d)))
#else
#define vshuffle_mask(a, b, c, d) (((a) << 0) | ((b) << 2) | ((c) << 4) | ((d) << 6))
#define vshuffle(x, y, a, b, c, d) (__builtin_ia32_shufps((x), (y), vshuffle_mask((a),(b),(c),(d))))
#endif

#define vsplat(v, x) vshuffle((v), (v), (x), (x),  (x),  (x))

static inline vec4 vxyz(vec4 x) __attribute__((always_inline));
static inline vec4 vxyz(vec4 x) { return vshuffle(x, vshuffle(x, vzero(), 2, 3, 0, 0), 0, 1, 0, 3); }

static inline vec4 vxyz1(vec4 x) __attribute__((always_inline));
static inline vec4 vxyz1(vec4 x) { return vshuffle(x, vshuffle(x, _mm_set_ss(1.0), 2, 2, 0, 1), 0, 1, 0, 2); }

#ifdef __FMA4__

static inline vec4 vmadd(vec4 a, vec4 b, vec4 c) __attribute__((always_inline));
static inline vec4 vmadd(vec4 a, vec4 b, vec4 c)
{
    return _mm_macc_ps(a, b, c);
}

static inline vec4 vnmadd(vec4 a, vec4 b, vec4 c) __attribute__((always_inline));
static inline vec4 vnmadd(vec4 a, vec4 b, vec4 c)
{
return _mm_nmacc_ps(a, b, c);
}

static inline vec4 vmsub(vec4 a, vec4 b, vec4 c) __attribute__((always_inline));
static inline vec4 vmsub(vec4 a, vec4 b, vec4 c)
{
    return _mm_msub_ps(a, b, c);
}

static inline vec4 vnmsub(vec4 a, vec4 b, vec4 c) __attribute__((always_inline));
static inline vec4 vnmsub(vec4 a, vec4 b, vec4 c)
{
    return _mm_nmsub_ps(a, b, c);
}

#else

static inline vec4 vmadd(vec4 a, vec4 b, vec4 c) __attribute__((always_inline));
static inline vec4 vmadd(vec4 a, vec4 b, vec4 c)
{
    return a * b + c;
}

static inline vec4 vnmadd(vec4 a, vec4 b, vec4 c) __attribute__((always_inline));
static inline vec4 vnmadd(vec4 a, vec4 b, vec4 c)
{
    return -(a * b) + c;
}

static inline vec4 vmsub(vec4 a, vec4 b, vec4 c) __attribute__((always_inline));
static inline vec4 vmsub(vec4 a, vec4 b, vec4 c)
{
    return a * b - c;
}

static inline vec4 vnmsub(vec4 a, vec4 b, vec4 c) __attribute__((always_inline));
static inline vec4 vnmsub(vec4 a, vec4 b, vec4 c)
{
    return -(a * b) - c;
}

#endif


#ifdef __SSE4_1__

static inline vec4 vblend(vec4 a, vec4 b, vec4 mask) __attribute__((always_inline));
static inline vec4 vblend(vec4 a, vec4 b, vec4 mask) { return _mm_blendv_ps(a, b, mask); }

static inline vec4 vdot(vec4 x, vec4 y) __attribute__((always_inline));
static inline vec4 vdot(vec4 x, vec4 y)
{
    return _mm_dp_ps(x, y, 0xff);
}

static inline vec4 vdot3(vec4 x, vec4 y) __attribute__((always_inline));
static inline vec4 vdot3(vec4 x, vec4 y)
{
    return _mm_dp_ps(x, y, 0x7f);
}

#else

// static inline vec4 vblend(vec4 a, vec4 b, vec4 mask); __attribute__((always_inline));
// static inline vec4 vblend(vec4 a, vec4 b, vec4 mask);
// static inline vec4 vblend_mask(vec4 a, vec4 b, const int mask); __attribute__((always_inline));
// static inline vec4 vblend_mask(vec4 a, vec4 b, const int mask);

#ifdef __SSE3__

static inline vec4 vdot(vec4 x, vec4 y) __attribute__((always_inline));
static inline vec4 vdot(vec4 x, vec4 y)
{
    vec4 prod = x * y;
    vec4 sum1 = _mm_hadd_ps(prod, prod);
    vec4 sum2 = _mm_hadd_ps(sum1, sum1);
    return sum2;
}

#else

static inline vec4 vdot(vec4 x, vec4 y) __attribute__((always_inline));
static inline vec4 vdot(vec4 x, vec4 y)
{
    vec4 prod = x * y;
    vec4 sum1 = prod + vshuffle(prod, prod, 1, 0, 3, 2);
    vec4 sum2 = sum1 + vshuffle(sum1, sum1, 2, 2, 0, 0);
    return sum2;
}

#endif

static inline vec4 vdot3(vec4 x, vec4 y) __attribute__((always_inline));
static inline vec4 vdot3(vec4 x, vec4 y)
{
    return vdot(vxyz(x), vxyz(y));

}
#endif

static inline vec4 vcross_scalar(vec4 x, vec4 y) __attribute__((always_inline));
static inline vec4 vcross_scalar(vec4 x, vec4 y)
{
    float *xs = (float*)&x, *ys = (float*)&y;

    vec4 result = {
        xs[1]*ys[2] - xs[2]*ys[1],
        xs[2]*ys[0] - xs[0]*ys[2],
        xs[0]*ys[1] - xs[1]*ys[0],
        0.0 };
    return result;
}

static inline vec4 vcross(vec4 x, vec4 y) __attribute__((always_inline));
static inline vec4 vcross(vec4 x, vec4 y)
{
    return vshuffle(x, x, 1, 2, 0, 3) * vshuffle(y, y, 2, 0, 1, 3)
        - vshuffle(x, x, 2, 0, 1, 3) * vshuffle(y, y, 1, 2, 0, 3);
}

static inline vec4 vmag(vec4 x) __attribute__((always_inline));
static inline vec4 vmag(vec4 x) { return vsqrt(vdot(x, x)); }
static inline vec4 vmags(vec4 x) __attribute__((always_inline));
static inline vec4 vmags(vec4 x) { return ssqrt(vdot(x, x)); }
static inline vec4 vmag3(vec4 x) __attribute__((always_inline));
static inline vec4 vmag3(vec4 x) { return vsqrt(vdot3(x, x)); }
static inline vec4 vmags3(vec4 x) __attribute__((always_inline));
static inline vec4 vmags3(vec4 x) { return ssqrt(vdot3(x, x)); }

static inline vec4 vunit(vec4 x) __attribute__((always_inline));
static inline vec4 vunit(vec4 x) { return x * vrsqrt(vdot(x, x)); }
static inline vec4 vunit3(vec4 x) __attribute__((always_inline));
static inline vec4 vunit3(vec4 x) { return x * vrsqrt(vdot3(x, x)); }

#endif
