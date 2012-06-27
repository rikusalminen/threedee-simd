#ifndef THREEDEE_PROJECTION_H
#define THREEDEE_PROJECTION_H

static inline mat4 mat_ortho(scalar left, scalar right, scalar bottom, scalar top, scalar near, scalar far);

static inline mat4 mat_frustum(scalar left, scalar right, scalar bottom, scalar top, scalar near, scalar far);
static inline mat4 mat_frustum_inf_z(scalar left, scalar right, scalar bottom, scalar top, scalar near);

static inline mat4 mat_projection_fovy(scalar fovy, scalar aspect, scalar near, scalar far);
static inline mat4 mat_projection_fovx(scalar fovx, scalar aspect, scalar near, scalar far);

static inline mat4 mat_projection_fovy_inf_z(scalar fovy, scalar aspect, scalar near, scalar far);
static inline mat4 mat_projection_fovx_inf_z(scalar fovx, scalar aspect, scalar near, scalar far);

#endif
