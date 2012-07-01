#ifndef THREEDEE_TYPES_H
#define THREEDEE_TYPES_H

typedef float vec4 __attribute__((vector_size(16)));
typedef float scalar;

struct mat4_t
{
    vec4 cols[4];
} __attribute__((aligned(16)));

typedef struct mat4_t mat4;

#endif
