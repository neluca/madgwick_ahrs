/*
 * Copyright Yinan Liao. and other contributors. All rights reserved.
 *
 * Reference: https://github.com/neluca/madgwick_ahrs/blob/main/madgwick_internal_report.pdf
 * MIT: https://github.com/neluca/madgwick_ahrs/blob/main/LICENSE
 */

#ifndef SRC_IMU_H
#define SRC_IMU_H

#define IMU_USE_ZETA

extern float imu_beta;  // algorithm gain beta

#ifdef IMU_USE_ZETA
extern float imu_zeta;  // algorithm gain zeta
#endif


void imu_update(float *q,
                float ax, float ay, float az, // Accelerometer (m/s2)
                float gx, float gy, float gz, // Gyroscope (rad/s)
                float dt);                    // Sampling time (s)


void imu_ahrs_update(float *q,
                     float ax, float ay, float az, // Accelerometer (m/s2)
                     float gx, float gy, float gz, // Gyroscope (rad/s)
                     float mx, float my, float mz, // Magnetometer (uT)
                     float dt);                    // Sampling time (s)

#endif  /* SRC_IMU_H */
