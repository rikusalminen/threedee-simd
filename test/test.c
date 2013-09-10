/**
 * Copyright (c) 2013 Nicolas Hillegeer
 * Copyright (c) 2013 Riku Salminen
 *
 * All matrices used by this library are in column-major format
 *
 * The reference values for cos/sin/tan et cetera are calculated
 * by the standard C library, as they are slow but accurate
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <threedee/threedee.h>

/* default tolerance */
#define TOL 0.001f

#define T_PI 3.14159265f
#define T_PI2 (T_PI * 0.50f)
#define T_PI3 (T_PI / 3.00f)
#define T_PI4 (T_PI * 0.25f)
#define T_PI6 (T_PI / 6.00f)

static const char *printv(vec4 vec) {
    static char buffer[256];

    const float *x = (const float*)&vec;
    snprintf(buffer, 256, "%2.3f\t%2.3f\t%2.3f\t%2.3f", x[0], x[1], x[2], x[3]);

    return buffer;
}

static const char *printvi(vec4 vec) {
    static char buffer[256];

    const uint32_t *x = (const uint32_t*)&vec;
    snprintf(buffer, 256, "%8x\t%8x\t%8x\t%8x\n", x[0], x[1], x[2], x[3]);

    return buffer;
}

static const char *printm(mat4 mat) {
    static char buffer[256];
    int accum = 0;

    const float *m = (const float*)&mat;

    accum += snprintf(buffer + accum, 256 - accum, "\n");
    for (int i = 0; i < 4; ++i) {
        accum += snprintf(buffer + accum, 256 - accum, "%2.3f\t%2.3f\t%2.3f\t%2.3f\n", m[i + 0], m[i + 4], m[i + 8], m[i + 12]);
    }
    accum += snprintf(buffer + accum, 256 - accum, "\n");

    return buffer;
}

const char *printm3(mat4 mat) {
    static char buffer[256];
    int accum = 0;

    const float *m = (const float*)&mat;

    accum += snprintf(buffer + accum, 256 - accum, "\n");
    for (int i = 0; i < 3; ++i) {
        accum += snprintf(buffer + accum, 256 - accum, "%2.3f\t%2.3f\t%2.3f\n", m[i + 0], m[i + 4], m[i + 8]);
    }
    accum += snprintf(buffer + accum, 256 - accum, "\n");

    return buffer;
}

/**
 * test if two matrices are equal to each other, within a certain tolerance
 *
 * return 1 if the matrices are (approximately) equal, 0 if not
 */
static int mat_equal_tol(mat4 am, mat4 bm, float tol) {
    const float *a = (const float*)&am;
    const float *b = (const float*)&bm;

    for (int i = 0; i < 16; ++i) {
        if (fabsf(b[i] - a[i]) > tol) {
            return 0;
        }
    }

    return 1;
}

static void compare_mat_equal(mat4 a, mat4 b, const char *name) {
    if (!mat_equal_tol(a, b, TOL)) {
        fprintf(stderr, "FAILED: %s\n", name);
        fprintf(stderr, "%s\n", printm(a));
        fprintf(stderr, "!=\n");
        fprintf(stderr, "%s\n", printm(b));
    }
    else {
        fprintf(stdout, "PASSED: %s\n", name);
    }
}

static void test_quat_to_mat() {
    const float cos45 = cosf(T_PI4);
    const float sin45 = sinf(T_PI4);

    /* rotate 45 degrees around the y-axis */
    {
        mat4 rot_y_45_ref = {{
            { cos45,    0,  -sin45,  0 },
            { 0,        1,  0,      0 },
            { sin45,    0,  cos45,  0 },
            { 0,        0,  0,      1 }
        }};

        {
            const char *name  = "quaternion to matrix (default)";
            vec4 axisangle    = vec(0.0f, 1.0f, 0.0f, T_PI4);
            vec4 qua          = quat_axisangle(axisangle);
            mat4 mat_from_qua = quat_to_mat(qua);

            compare_mat_equal(rot_y_45_ref, mat_from_qua, name);
        }

        {
            const char *name  = "quaternion to matrix (scalar)";
            vec4 axisangle    = vec(0.0f, 1.0f, 0.0f, T_PI4);
            vec4 qua          = quat_axisangle(axisangle);
            mat4 mat_from_qua = quat_to_mat_scalar(qua);

            compare_mat_equal(rot_y_45_ref, mat_from_qua, name);
        }
    }
}

int main() {
    test_quat_to_mat();

    return 0;
}
