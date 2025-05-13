/*
 * Copyright Yinan Liao. All rights reserved.
 *
 * MIT: https://github.com/neluca/madgwick_ahrs/blob/main/LICENSE
 */

#ifndef SRC_QUATERNION_H
#define SRC_QUATERNION_H

#include <math.h>

typedef float fp_t;
#define q_inv_sqrt(a)  (1 / sqrtf(a))
#define q_atan2(a, b)  atan2f(a, b)
#define q_asin(a)      asinf(a)

// q = l * r
static inline void quaternion_prod(fp_t *q, const fp_t *l, const fp_t *r) {
    q[0] = l[0] * r[0] - l[1] * r[1] - l[2] * r[2] - l[3] * r[3];
    q[1] = l[0] * r[1] + l[1] * r[0] + l[2] * r[3] - l[3] * r[2];
    q[2] = l[0] * r[2] - l[1] * r[3] + l[2] * r[0] + l[3] * r[1];
    q[3] = l[0] * r[3] + l[1] * r[2] - l[2] * r[1] + l[3] * r[0];
}

// q += r
static inline void quaternion_add(fp_t *q, const fp_t *r) {
    q[0] += r[0];
    q[1] += r[1];
    q[2] += r[2];
    q[3] += r[3];
}

// q -= r
static inline void quaternion_sub(fp_t *q, const fp_t *r) {
    q[0] -= r[0];
    q[1] -= r[1];
    q[2] -= r[2];
    q[3] -= r[3];
}

// scalar multiplication : q *= n
static inline void quaternion_scalar(fp_t *q, fp_t n) {
    q[0] *= n;
    q[1] *= n;
    q[2] *= n;
    q[3] *= n;
}

// q* = [s, -v] if q = [s, v]
static inline void quaternion_conj(fp_t *q_hat, const fp_t *q) {
    q_hat[0] = q[0];
    q_hat[1] = -q[1];
    q_hat[2] = -q[2];
    q_hat[3] = -q[3];
}

// normalize q
static inline void quaternion_norm(fp_t *q) {
    float norm = q_inv_sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
    q[0] *= norm;
    q[1] *= norm;
    q[2] *= norm;
    q[3] *= norm;
}

// quaternion to euler angles
// x: roll  y: pitch  z: yaw
static inline void quaternion2euler(fp_t *x, fp_t *y, fp_t *z, const fp_t *q) {
    *z = q_atan2((2 * q[1] * q[2] - 2 * q[0] * q[3]), (2 * q[0] * q[0] + 2 * q[1] * q[1] - 1));
    *y = -q_asin(2 * q[1] * q[3] + 2 * q[0] * q[2]);
    *x = q_atan2((2 * q[2] * q[3] - 2 * q[0] * q[1]), (2 * q[0] * q[0] + 2 * q[3] * q[3] - 1));
}

#define Q_(sign_name) fp_t sign_name[4]

#endif  /* SRC_QUATERNION_H */
