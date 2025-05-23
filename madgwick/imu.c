/*
 * Copyright Yinan Liao. and other contributors. All rights reserved.
 *
 * Reference: https://github.com/neluca/madgwick_ahrs/blob/main/madgwick_internal_report.pdf
 * MIT: https://github.com/neluca/madgwick_ahrs/blob/main/LICENSE
 */

#include "imu.h"
#include "quaternion.h"

float imu_beta = 1.0f;

#ifdef IMU_USE_ZETA
float imu_zeta = 0.0f;
#endif

void imu_update(float *q, float ax, float ay, float az, float gx, float gy, float gz, float dt) {
    Q_(q_a) = {0, ax, ay, az};
    quaternion_norm(q_a);

    // Compute the objective function for gravity,
    // simplified to equation due to the 0's in the acceleration reference quaternion
    float f[3];
    f[0] = 2.0f * (q[1] * q[3] - q[0] * q[2]) - q_a[1];
    f[1] = 2.0f * (q[0] * q[1] - q[2] * q[3]) - q_a[2];
    f[2] = 2.0f * (0.5f - q[1] * q[1] - q[2] * q[2]) - q_a[3];

    // Compute the Jacobian matrix for gravity
    float j[3][4];
    j[0][0] = -2.0f * q[2];
    j[0][1] = 2.0f * q[3];
    j[0][2] = -2.0f * q[0];
    j[0][3] = 2.0f * q[1];

    j[1][0] = 2.0f * q[1];
    j[1][1] = 2.0f * q[0];
    j[1][2] = 2.0f * q[3];
    j[1][3] = 2.0f * q[2];

    j[2][0] = 0.0f;
    j[2][1] = -4.0f * q[1];
    j[2][2] = -4.0f * q[2];
    j[2][3] = 0.0f;

    // Gradient decent algorithm corrective step, step = j' * f
    Q_(step);
    step[0] = j[0][0] * f[0] + j[1][0] * f[1] + j[2][0] * f[2];
    step[1] = j[0][1] * f[0] + j[1][1] * f[1] + j[2][1] * f[2];
    step[2] = j[0][2] * f[0] + j[1][2] * f[1] + j[2][2] * f[2];
    step[3] = j[0][3] * f[0] + j[1][3] * f[1] + j[2][3] * f[2];

    quaternion_norm(step);
    quaternion_scale(step, imu_beta);

    // Compute rate of change of quaternion
    Q_(q_g) = {0, gx, gy, gz};
    Q_(q_dot);
    quaternion_prod(q_dot, q, q_g);
    quaternion_scale(q_dot, 0.5f);
    quaternion_sub(q_dot, step);

    // Integrate to yield quaternion
    quaternion_scale(q_dot, dt);
    quaternion_add(q, q_dot);
    quaternion_norm(q);
}


void imu_ahrs_update(float *q,
                     float ax, float ay, float az, float gx, float gy, float gz, float mx, float my, float mz,
                     float dt
) {
    Q_(q_a) = {0, ax, ay, az};
    quaternion_norm(q_a);

    Q_(q_m) = {0, mx, my, mz};
    quaternion_norm(q_m);

    // Reference direction of Earth's magnetic field
    Q_(h);
    quaternion_prod(h, q, q_m);
    Q_(q_hat);
    quaternion_conj(q_hat, q);
    quaternion_prod(h, h, q_hat);

    float b[4];
    b[0] = 0.0f;
    b[1] = sqrtf(h[1] * h[1] + h[2] * h[2]);
    b[2] = 0.0f;
    b[3] = h[3];

    // Compute the objective function for gravity and magnetic field
    float f[6];
    f[0] = 2.0f * (q[1] * q[3] - q[0] * q[2]) - q_a[1];
    f[1] = 2.0f * (q[0] * q[1] - q[2] * q[3]) - q_a[2];
    f[2] = 2.0f * (0.5f - q[1] * q[1] - q[2] * q[2]) - q_a[3];
    f[3] = 2.0f * b[1] * (0.5f - q[2] * q[2] - q[3] * q[3]) + 2.0f * b[3] * (q[1] * q[3] - q[0] * q[2]) - q_m[1];
    f[4] = 2.0f * b[1] * (q[1] * q[2] - q[0] * q[3]) + 2.0f * b[3] * (q[0] * q[1] + q[2] * q[3]) - q_m[2];
    f[5] = 2.0f * b[1] * (q[0] * q[2] + q[1] * q[3]) + 2.0f * b[3] * (0.5f - q[1] * q[1] - q[2] * q[2]) - q_m[3];

    // Compute the Jacobian matrix for gravity and magnetic field
    float j[6][4];
    j[0][0] = -2.0f * q[2];
    j[0][1] = 2.0f * q[3];
    j[0][2] = -2.0f * q[0];
    j[0][3] = 2.0f * q[1];

    j[1][0] = 2.0f * q[1];
    j[1][1] = 2.0f * q[0];
    j[1][2] = 2.0f * q[3];
    j[1][3] = 2.0f * q[2];

    j[2][0] = 0.0f;
    j[2][1] = -4.0f * q[1];
    j[2][2] = -4.0f * q[2];
    j[2][3] = 0.0f;

    j[3][0] = -2.0f * b[3] * q[2];
    j[3][1] = 2.0f * b[3] * q[3];
    j[3][2] = -4.0f * b[1] * q[2] - 2 * b[3] * q[0];
    j[3][3] = -4.0f * b[1] * q[3] + 2 * b[3] * q[1];

    j[4][0] = -2.0f * b[1] * q[3] + 2 * b[3] * q[1];
    j[4][1] = 2.0f * b[1] * q[2] + 2 * b[3] * q[0];
    j[4][2] = 2.0f * b[1] * q[1] + 2 * b[3] * q[3];
    j[4][3] = -2.0f * b[1] * q[0] + 2 * b[3] * q[2];

    j[5][0] = 2.0f * b[1] * q[2];
    j[5][1] = 2.0f * b[1] * q[3] - 4 * b[3] * q[1];
    j[5][2] = 2.0f * b[1] * q[0] - 4 * b[3] * q[2];
    j[5][3] = 2.0f * b[1] * q[1];

    // Gradient decent algorithm corrective step, step = j' * f
    Q_(step);
    step[0] = j[0][0] * f[0] + j[1][0] * f[1] + j[2][0] * f[2] + j[3][0] * f[3] + j[4][0] * f[4] + j[5][0] * f[5];
    step[1] = j[0][1] * f[0] + j[1][1] * f[1] + j[2][1] * f[2] + j[3][1] * f[3] + j[4][1] * f[4] + j[5][1] * f[5];
    step[2] = j[0][2] * f[0] + j[1][2] * f[1] + j[2][2] * f[2] + j[3][2] * f[3] + j[4][2] * f[4] + j[5][2] * f[5];
    step[3] = j[0][3] * f[0] + j[1][3] * f[1] + j[2][3] * f[2] + j[3][3] * f[3] + j[4][3] * f[4] + j[5][3] * f[5];

    quaternion_norm(step);
    quaternion_scale(step, imu_beta);

    Q_(q_g) = {0, gx, gy, gz};
#ifdef IMU_USE_ZETA
    // Gyroscope compensation drift
    Q_(q_tmp);
    quaternion_prod(q_tmp, q_hat, step);
    quaternion_scale(q_tmp, 2.0f);
    quaternion_scale(q_tmp, dt);
    quaternion_scale(q_tmp, imu_zeta);
    quaternion_sub(q_g, q_tmp);
#endif

    // Compute rate of change of quaternion
    Q_(q_dot);
    quaternion_prod(q_dot, q, q_g);
    quaternion_scale(q_dot, 0.5f);
    quaternion_sub(q_dot, step);

    // Integrate to yield quaternion
    quaternion_scale(q_dot, dt);
    quaternion_add(q, q_dot);
    quaternion_norm(q);
}
