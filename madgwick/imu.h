/*
 * Copyright Yinan Liao. and other contributors. All rights reserved.
 *
 * Reference: https://github.com/neluca/madgwick_ahrs/blob/main/madgwick_internal_report.pdf
 * MIT: https://github.com/neluca/madgwick_ahrs/blob/main/LICENSE
 */

#ifndef SRC_IMU_H
#define SRC_IMU_H

#define IMU_USE_ZETA

extern float imu_beta;

#ifdef IMU_USE_ZETA
extern float imu_zeta;
#endif

void imu_update(float *q,
                float ax, float ay, float az,
                float gx, float gy, float gz,
                float dt);

void imu_ahrs_update(float *q,
                     float ax, float ay, float az,
                     float gx, float gy, float gz,
                     float mx, float my, float mz,
                     float dt);

#endif  /* SRC_IMU_H */
