cmake_minimum_required(VERSION 3.10)

project(madgwick_ahrs C)

set(CMAKE_C_STANDARD 11)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/madgwick)

add_executable(imu
        main.c
        madgwick/imu.c
        madgwick/imu.h)

add_executable(imu_ahrs
        main_ahrs.c
        madgwick/imu.c
        madgwick/imu.h)

add_executable(imu_test
        dataset_test.c
        madgwick/imu.c
        madgwick/imu.h)