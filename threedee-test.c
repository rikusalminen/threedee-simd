#include <stdio.h>

#include <threedee/threedee.h>


void printv(vec4 vec)
{
    const float *x = (const float*)&vec;
    printf("%2.3f\t%2.3f\t%2.3f\t%2.3f\n", x[0], x[1], x[2], x[3]);
}

void printvi(vec4 vec)
{
    const uint32_t *x = (const uint32_t*)&vec;
    printf("%8x\t%8x\t%8x\t%8x\n", x[0], x[1], x[2], x[3]);
}

void printm(mat4 mat)
{
    const float *m = (const float*)&mat;
    printf("\n");
    for(int i = 0; i < 4; ++i)
        printf("%2.3f\t%2.3f\t%2.3f\t%2.3f\n", m[i + 0], m[i + 4], m[i + 8], m[i + 12]);
    printf("\n");
}

int main(int argc, char *argv[])
{
    (void)argc; (void)argv;
#ifdef __SSE__
    printf("SSE ");
#endif
#ifdef __SSE2__
    printf("SSE2 ");
#endif
#ifdef __SSE3__
    printf("SSE3 ");
#endif
#ifdef __SSE4__
    printf("SSE4 ");
#endif
#ifdef __SSE4_1__
    printf("SSE4.1 ");
#endif
#ifdef __SSE4_2__
    printf("SSE4.2 ");
#endif
#ifdef __AVX__
    printf("AVX ");
#endif
#ifdef __FMA4__
    printf("FMA4 ");
#endif
    printf("\n");

    printv(vec(1, 2, 3, 4));
    printv(vzero());

    printm(mzero());
    printm(midentity());

    vec4 a = { 1, 2, 3, 4 }, b = { 5, 6, 7, 8 };

    printv(a);
    printv(b);

    printf("\nshuffles:\n");
    printv(vshuffle(a, a, 0, 1, 2, 3));
    printv(vshuffle(a, a, 3, 2, 1, 0));
    printv(vshuffle(a, b, 0, 1, 0, 1));
    printv(vshuffle(a, b, 2, 3, 2, 3));

    printf("\ndot products:\n");

    printv(vdot(a, b));
    printv(vdot(b, a));

    printv(vdot3(a, b));
    printv(vdot3(b, a));

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

    printf("\nperspective matrix:\n");
    printm(mat_perspective_fovy((float)M_PI/4.0f, 16.0f/9.0f, 0.1f, 100.0f));
    printm(mat_perspective_fovy_inf_z((float)M_PI/4.0f, 16.0f/9.0f, 0.1f));
    printm(mat_perspective_fovy_scalar((float)M_PI/4.0f, 16.0f/9.0f, 0.1f, 100.0f));

    printf("\northogonal matrix:\n");
    printm(mat_ortho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0));
    printm(mat_ortho(-1.0, 2.0, -1.0, 2.0, -1.0, 2.0));

    printf("\ntranslate matrix:\n");
    printm(mtranslate(a));

    printf("\nscale matrix:\n");
    printm(mscale(a));

    return 0;
}
