#include <stdio.h>
#include "imu.h"
#include "quaternion.h"

#define SIMPLING_TIME 0.01f

int main(void) {
    float roll, pitch, yaw;
    Q_(q) = {1.0f, 0.0f, 0.0f, 0.0f};
    imu_beta = 3.1415926f * (5.0f / 180.0f) * sqrtf(3.0f / 4.0f);

    printf("Upright to solve gradient descent problem\n");
    for (int i = 0; i < 1000; i++) {
        imu_ahrs_update(q, 0.05f, 0.05f, 9.4f, 0, 0, 0, 0, 10, 0, SIMPLING_TIME);
        quaternion2euler(&roll, &pitch, &yaw, q);
        printf("Time (s): %.3f, roll: %f, pitch: %f, yaw: %f\n", ((float) i * SIMPLING_TIME), roll, pitch, yaw);
    }

    return 0;
}