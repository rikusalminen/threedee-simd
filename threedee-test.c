#include <stdio.h>

#include <threedee/types.h>
#include <threedee/vector.h>
#include <threedee/matrix.h>
#include <threedee/rotation.h>

void printv(vec4 vec)
{
    float x[4];
    _mm_storeu_ps(x, vec);
    printf("%2.3f\t%2.3f\t%2.3f\t%2.3f\n", x[0], x[1], x[2], x[3]);
}

void printvi(vec4 vec)
{
    uint32_t x[4];
    _mm_storeu_ps((float*)x, vec);
    printf("%8x\t%8x\t%8x\t%8x\n", x[0], x[1], x[2], x[3]);
}

void printm(mat4 mat)
{
    printf("\n");
    for(int i = 0; i < 4; ++i)
        printv(mat.rows[i]);
    printf("\n");
}

int main(int argc, char *argv[])
{
    (void)argc; (void)argv;

    printv(vec(1, 2, 3, 4));
    printv(vzero());

    printm(mzero());
    printm(midentity());

    vec4 a = { 1, 2, 3, 4 }, b = { 5, 6, 7, 8 };

    //vec4 blendmask = { 1, -1, 1, -1 };
    //printv(vblend(x, y, blendmask));

    vec4 x = { 1, 0, 0, 0 }, y = { 0, 1, 0, 0 }, z = { 0, 0, 1, 0 }, w = { 0, 0, 0, 1 };

    printf("\ncross products:\n");

    printv(vcross(x, y));
    printv(vcross(y, x));

    printv(vcross_scalar(x, y));
    printv(vcross_scalar(y, x));

    printf("\nquaternion products:\n");

    printv(qprod(x, y));
    printv(qprod(y, x));

    printv(qprod_mad(x, y));
    printv(qprod_mad(y, x));

    printv(qprod_scalar(x, y));
    printv(qprod_scalar(y, x));

    printf("\nquaternion conjugates:\n");

    printv(qconj(x));
    printv(qconj(y));
    printv(qconj(z));
    printv(qconj(w));

    printf("\nmat from quat:\n");
    printm(quat_to_mat(w));
    printm(quat_to_mat_mmul(w));
    printm(quat_to_mat_scalar(w));

    vec4 angles = { 0.0, 0.0, 0.0, 0.0 };
    printf("\neuler to quaternion:\n");
    printv(quat_euler(angles));
    printv(quat_euler_scalar(angles));
    printv(quat_euler_gems(angles));

    printf("\neuler to matrix:\n");
    printm(mat_euler(angles));
    printm(mat_euler_scalar(angles));
    printm(quat_to_mat(quat_euler(angles)));

    return 0;
}
