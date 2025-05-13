/*
 * Copyright Yinan Liao. and other contributors. All rights reserved.
 */

#ifndef SRC_IMU_H
#define SRC_IMU_H

extern float imu_beta;
extern float imu_zeta;

void imu_update(float *q, float ax, float ay, float az, float gx, float gy, float gz, float dt);

void imu_ahrs_update(float *q,
                     float ax, float ay, float az,
                     float gx, float gy, float gz,
                     float mx, float my, float mz,
                     float dt);

#endif  /* SRC_IMU_H */
