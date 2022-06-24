#ifndef THREEDEE_TYPES_H
#define THREEDEE_TYPES_H

#include <stdint.h>
#include <stddef.h>

/* scalar type */
#if defined(__EXTENDED_PRECSION_SCALAR__)
typedef long double scalar;
#elif defined(__DOUBLE_PRECISION_SCALAR__)
typedef double scalar;
#else
typedef float scalar;
#endif

/* alignment of vector type */
#if ! ((defined __sun__ || defined __VXWORKS__) && defined __i386__)
#define REQUIRED_ALIGNMENT (uint8_t)(4 * sizeof(scalar) < __alignof__(max_align_t) ? 4 * sizeof(scalar) : __alignof__(max_align_t))
#else
#define REQUIRED_ALIGNMENT (uint8_t)(4 * sizeof(scalar))
#endif
#define SCALAR_SUFFIX s

/* vector types */
typedef scalar vec2 __attribute__((vector_size( sizeof(scalar) * 2 )));
typedef scalar vec4 __attribute__((vector_size( sizeof(scalar) * 4 )));
typedef float vec2f __attribute__((vector_size( sizeof(float) * 2 )));
typedef float vec4f __attribute__((vector_size( sizeof(float) * 4 )));
typedef double vec2d __attribute__((vector_size( sizeof(double) * 2 )));
typedef double vec4d __attribute__((vector_size( sizeof(double) * 4 )));

/* matrix type */
struct mat4_t
{
    vec4 cols[4];
} __attribute__((aligned( REQUIRED_ALIGNMENT )));

typedef struct mat4_t mat4;

#endif
