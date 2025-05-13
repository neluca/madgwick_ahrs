cmake_minimum_required(VERSION 3.10)

project(madgwick_ahrs C)

set(CMAKE_C_STANDARD 11)

include_directories(madgwick)

add_executable(imu main.c madgwick/imu.c madgwick/imu.h)

add_executable(imu_ahrs main_ahrs.c madgwick/imu.c madgwick/imu.h)