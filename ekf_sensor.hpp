#ifndef EKF_SENSER_HPP
#define EKF_SENSER_HPP
#include <Eigen/Dense>
#include <stdio.h>
#include<cmath>
#include<random>
#include<cmath>
#include <iostream>
#include "pico/stdlib.h"
#include <string.h>
#include "pico/binary_info.h"
#include "hardware/spi.h"
#include "lsm9ds1_reg.h"

#define SENSOR_BUS spi0
#define BOOT_TIME 20 //ms
#define PIN_CSAG  1
#define PIN_MISO  4
#define PIN_CSM   5
#define PIN_SCK   6
#define PIN_MOSI  7
#define GRAV (9.80665)
#define MN (-315.0)
#define MD (440.0)
#define PI (3.14159)

using Eigen::MatrixXd;
using Eigen::MatrixXf;
using Eigen::Matrix4f;
using Eigen::Matrix;
using Eigen::PartialPivLU;
using namespace Eigen;

extern  int16_t data_raw_acceleration[3];
extern  int16_t data_raw_angular_rate[3];
extern  int16_t data_raw_magnetic_field[3];
extern  float acceleration_mg[3];
extern  float angular_rate_mdps[3];
extern  float magnetic_field_mgauss[3];
extern  lsm9ds1_id_t whoamI;
extern  lsm9ds1_status_t reg;
extern  uint8_t rst;
extern  uint8_t tx_buffer_imu[1000];
extern  uint8_t tx_buffer_mag[1000];
extern float Acceleration_mg[3];
extern float Angular_rate_mdps[3];
extern float Magnetic_field_mgauss[3];
extern sensbus_t Ins_bus;
extern sensbus_t Mag_bus;
extern stmdev_ctx_t Imu_h;
extern stmdev_ctx_t Mag_h;

//std::random_device rnd;
//std::mt19937 mt(rnd());  
//std::normal_distribution<> norm(0.0, 1.0);

uint8_t observation_equation(Matrix<float, 7, 1>x, Matrix<float, 6, 1>&z, float g, float mn, float md);
void imu_mag_data_read(void);
void imu_mag_init(void);
float Psi(Matrix<float, 7, 1>x);
float Theta(Matrix<float, 7, 1>x);
float Phl(Matrix<float, 7, 1>x);
uint8_t ekf( Matrix<float, 7, 1> &xe,
             Matrix<float, 7, 1> &xp,
             Matrix<float, 7, 7> &P,
             Matrix<float, 6, 1> z,
             Matrix<float, 3, 1> omega,
             Matrix<float, 3, 3> Q, 
             Matrix<float, 6, 6> R, 
             Matrix<float, 7, 3> G,
             Matrix<float, 3, 1> beta,

             float dt);
#endif
