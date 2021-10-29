#include "pico/stdlib.h"
#include "hardware/pwm.h"
#include "ekf_sensor.hpp"

#define DUTYMIN 1330
#define DUTYMAX 1520
int16_t data_raw_acceleration[3];
int16_t data_raw_angular_rate[3];
int16_t data_raw_magnetic_field[3];
float acceleration_mg[3];
float angular_rate_mdps[3];
float magnetic_field_mgauss[3];
lsm9ds1_id_t whoamI;
lsm9ds1_status_t reg;
uint8_t rst;
uint8_t tx_buffer_imu[1000];
uint8_t tx_buffer_mag[1000];
float Acceleration_mg[3];
float Angular_rate_mdps[3];
float Magnetic_field_mgauss[3];
sensbus_t Ins_bus={spi0, PIN_CSAG};
sensbus_t Mag_bus={spi0, PIN_CSM};
stmdev_ctx_t Imu_h;
stmdev_ctx_t Mag_h;
uint8_t slice_num[2];


uint8_t pwm_settei(){ 
     //pwmの設定
    gpio_set_function(2, GPIO_FUNC_PWM);
    gpio_set_function(3, GPIO_FUNC_PWM);
   
    gpio_set_function(4, GPIO_FUNC_PWM);
    gpio_set_function(5, GPIO_FUNC_PWM);
    // Find out which PWM slice is connected to GPIO 0 (it's slice 0)
    slice_num[0] = pwm_gpio_to_slice_num(3);
    slice_num[1] = pwm_gpio_to_slice_num(4);
    pwm_set_wrap(slice_num[0], 3124);
    pwm_set_clkdiv(slice_num[0], 100.0);
    pwm_set_wrap(slice_num[1], 3124);
    pwm_set_clkdiv(slice_num[1], 100.0);
    pwm_set_chan_level(slice_num[0], PWM_CHAN_A, DUTYMAX);
    pwm_set_chan_level(slice_num[0], PWM_CHAN_B, DUTYMAX);
    pwm_set_chan_level(slice_num[1], PWM_CHAN_A, DUTYMAX);
    pwm_set_chan_level(slice_num[1], PWM_CHAN_B, DUTYMAX);
    pwm_set_enabled(slice_num[0], true);
    pwm_set_enabled(slice_num[1], true);
    sleep_ms(2000);
    pwm_set_chan_level(slice_num[0], PWM_CHAN_A, DUTYMIN);
    pwm_set_chan_level(slice_num[0], PWM_CHAN_B, DUTYMIN);
    pwm_set_chan_level(slice_num[1], PWM_CHAN_A, DUTYMIN);
    pwm_set_chan_level(slice_num[1], PWM_CHAN_B, DUTYMIN);
    return 0;
}

//Runge-Kutta Method 
uint8_t rk4(uint8_t (*func)(float t, 
                            Matrix<float, 7, 1> x, 
                            Matrix<float, 3, 1> omega, 
                            Matrix<float, 3, 1> beta, 
                            Matrix<float, 7, 1> &k),
         		    float t, 
          		    float h, 
      		 	    Matrix<float, 7 ,1> &x, 
            		    Matrix<float, 3, 1> omega,
        		    Matrix<float, 3, 1> beta)
{
  Matrix<float, 7, 1> k1,k2,k3,k4;

  func(t,       x,            omega, beta, k1);
  func(t+0.5*h, x + 0.5*h*k1, omega, beta, k2);
  func(t+0.5*h, x + 0.5*h*k2, omega, beta, k3);
  func(t+h,     x +     h*k3, omega, beta, k4);
  x = x + h*(k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
  
  return 0;
}
#if 0
//Continuous Time State Equation for simulation
uint8_t xdot( float t, 
              Matrix<float, 7, 1> x, 
              Matrix<float, 3, 1> omega , 
              Matrix<float, 3, 1> beta, 
              Matrix<float, 7, 1> &k)
{
  float q0 = x(0,0);
  float q1 = x(1,0);
  float q2 = x(2,0);
  float q3 = x(3,0);
  float dp = x(4,0);
  float dq = x(5,0);
  float dr = x(6,0);
  float p = omega(0,0);
  float q = omega(1,0);
  float r = omega(2,0);
  float betax = beta(0,0);
  float betay = beta(1,0);
  float betaz = beta(2,0);
 
  k(0,0) = 0.5*(-p*q1 - q*q2 - r*q3);
  k(1,0) = 0.5*( p*q0 + r*q2 - q*q3);
  k(2,0) = 0.5*( q*q0 - r*q1 + p*q3);
  k(3,0) = 0.5*( r*q0 + q*q1 - p*q2);
  k(4,0) =-betax*dp + norm(mt);
  k(5,0) =-betay*dq + norm(mt);
  k(6,0) =-betaz*dr + norm(mt);

  return 0;
}
#endif
//Discrite Time State Equation
uint8_t state_equation( Matrix<float, 7, 1> &xe, 
                        Matrix<float, 3, 1> omega_m, 
                        Matrix<float, 3, 1> beta, 
                        float dt,
                        Matrix<float, 7, 1> &xp)
{
  float q0=xe(0, 0);
  float q1=xe(1, 0);
  float q2=xe(2, 0);
  float q3=xe(3, 0);
  float dp=xe(4, 0);
  float dq=xe(5, 0);
  float dr=xe(6, 0);
  float pm=omega_m(0, 0);
  float qm=omega_m(1, 0);
  float rm=omega_m(2, 0);
  float betax=beta(0, 0);
  float betay=beta(1, 0);
  float betaz=beta(2, 0);

  xp(0, 0) = q0 + 0.5*(-(pm-dp)*q1 -(qm-dq)*q2 -(rm-dr)*q3)*dt;
  xp(1, 0) = q1 + 0.5*( (pm-dp)*q0 +(rm-dr)*q2 -(qm-dq)*q3)*dt;
  xp(2, 0) = q2 + 0.5*( (qm-dq)*q0 -(rm-dr)*q1 +(pm-dp)*q3)*dt;
  xp(3, 0) = q3 + 0.5*( (rm-dr)*q0 +(qm-dq)*q1 -(pm-dp)*q2)*dt;
  xp(4, 0) = dp -betax*dp*dt;
  xp(5, 0) = dq -betay*dq*dt;
  xp(6, 0) = dr -betaz*dr*dt;
  
  return 0;
}

//Observation Equation
uint8_t observation_equation(Matrix<float, 7, 1>x, Matrix<float, 6, 1>&z, float g, float mn, float md){
  float q0 = x(0, 0);
  float q1 = x(1, 0);
  float q2 = x(2, 0);
  float q3 = x(3, 0);

  z(0, 0) = 2.0*(q1*q3 - q0*q2)*g;
  z(1, 0) = 2.0*(q2*q3 + q0*q1)*g;
  z(2, 0) = (q0*q0 - q1*q1 - q2*q2 + q3*q3)*g ;
  z(3, 0) = (q0*q0 + q1*q1 - q2*q2 - q3*q3)*mn + 2.0*(q1*q3 - q0*q2)*md;
  z(4, 0) = 2.0*(q1*q2 - q0*q3)*mn + 2.0*(q2*q3 + q0*q1)*md;
  z(5, 0) = 2.0*(q1*q3 + q0*q2)*mn + (q0*q0 - q1*q1 - q2*q2 + q3*q3)*md;
 // printf("%f %f %f %f %f %f\n",z(0,0),z(1,0),z(2,0),z(3,0),z(4,0),z(5,0));
  return 0;
}

//Make Jacobian matrix F
uint8_t F_jacobian( Matrix<float, 7, 7>&F, 
                    Matrix<float, 7, 1> x, 
                    Matrix<float, 3, 1> omega, 
                    Matrix<float, 3, 1> beta, 
                    float dt)
{
  float q0=x(0, 0);
  float q1=x(1, 0);
  float q2=x(2, 0);
  float q3=x(3, 0);
  float dp=x(4, 0);
  float dq=x(5, 0);
  float dr=x(6, 0);
  float pm=omega(0, 0);
  float qm=omega(1, 0);
  float rm=omega(2, 0);
  float betax=beta(0, 0);
  float betay=beta(1, 0);
  float betaz=beta(2, 0);

  //x(0, 0) = q0 + 0.5*(-(pm-dp)*q1 -(qm-dq)*q2 -(rm-dr)*q3)*dt;
  F(0, 0)= 1.0;
  F(0, 1)=-0.5*(pm -dp)*dt;
  F(0, 2)=-0.5*(qm -dq)*dt;
  F(0, 3)=-0.5*(rm -dr)*dt;
  F(0, 4)= 0.5*q1*dt;
  F(0, 5)= 0.5*q2*dt;
  F(0, 6)= 0.5*q3*dt;

  //x(1, 0) = q1 + 0.5*( (pm-dp)*q0 +(rm-dr)*q2 -(qm-dq)*q3)*dt;
  F(1, 0)= 0.5*(pm -dp)*dt;
  F(1, 1)= 1.0;
  F(1, 2)= 0.5*(rm -dr)*dt;
  F(1, 3)=-0.5*(qm -dq)*dt;
  F(1, 4)=-0.5*q0*dt;
  F(1, 5)= 0.5*q3*dt;
  F(1, 6)=-0.5*q2*dt;

  //x(2, 0) = q2 + 0.5*( (qm-dq)*q0 -(rm-dr)*q1 +(pm-dp)*q3)*dt; <-miss!
  F(2, 0)= 0.5*(qm -dq)*dt;
  F(2, 1)=-0.5*(rm -dr)*dt;
  F(2, 2)= 1.0;
  F(2, 3)= 0.5*(pm -dp)*dt;
  F(2, 4)=-0.5*q3*dt;
  F(2, 5)=-0.5*q0*dt;
  F(2, 6)= 0.5*q1*dt;
  
  //x(3, 0) = q3 + 0.5*( (rm-dr)*q0 +(qm-dq)*q1 -(pm-dp)*q2)*dt;
  F(3, 0)= 0.5*(rm -dr)*dt;
  F(3, 1)= 0.5*(qm -dq)*dt;
  F(3, 2)=-0.5*(pm -dp)*dt;
  F(3, 3)= 1.0;
  F(3, 4)= 0.5*q2*dt;
  F(3, 5)=-0.5*q1*dt;
  F(3, 6)=-0.5*q0*dt;
  
  //x(4, 0) = dp -betax*dp*dt;
  F(4, 0)= 0.0;
  F(4, 1)= 0.0;
  F(4, 2)= 0.0;
  F(4, 3)= 0.0;
  F(4, 4)= 1.0 - betax*dt;
  F(4, 5)= 0.0;
  F(4, 6)= 0.0;
  
  //x(5, 0) = dq - betay*dq*dt;
  F(5, 0)= 0.0;
  F(5, 1)= 0.0;
  F(5, 2)= 0.0;
  F(5, 3)= 0.0;
  F(5, 4)= 0.0;
  F(5, 5)= 1.0 - betay*dt;
  F(5, 6)= 0.0;
  
  //x(6, 0) = dr -betaz*dr*dt;
  F(6, 0)= 0.0;
  F(6, 1)= 0.0;
  F(6, 2)= 0.0;
  F(6, 3)= 0.0;
  F(6, 4)= 0.0;
  F(6, 5)= 0.0;
  F(6, 6)= 1.0 - betaz*dt;
  
  return 0;
}

//Make Jacobian matrix H
uint8_t H_jacobian( Matrix<float, 6, 7> &H, 
                    Matrix<float, 7, 1> x, 
                    float g, 
                    float mn, 
                    float md)
{
  float q0 = x(0, 0);
  float q1 = x(1, 0);
  float q2 = x(2, 0);
  float q3 = x(3, 0);

  //z(0, 0) = 2.0*(q1*q3 - q0*q2)*g;
  H(0, 0) =-2.0*q2*g;
  H(0, 1) = 2.0*q3*g;
  H(0, 2) =-2.0*q0*g;
  H(0, 3) = 2.0*q1*g;
  H(0, 4) = 0.0;
  H(0, 5) = 0.0;
  H(0, 6) = 0.0;
  
  //z(1, 0) = 2.0*(q2*q3 + q0*q1)*g;
  H(1, 0) = 2.0*q1*g;
  H(1, 1) = 2.0*q0*g;
  H(1, 2) = 2.0*q3*g;
  H(1, 3) = 2.0*q2*g;
  H(1, 4) = 0.0;
  H(1, 5) = 0.0;
  H(1, 6) = 0.0;
  
  //z(2, 0) = (q0*q0 - q1*q1 - q2*q2 + q3*q3)*g ;
  H(2, 0) = 2.0*q0*g;
  H(2, 1) =-2.0*q1*g;
  H(2, 2) =-2.0*q2*g;
  H(2, 3) = 2.0*q3*g;
  H(2, 4) = 0.0;
  H(2, 5) = 0.0;
  H(2, 6) = 0.0;
  
  //z(3, 0) = (q0*q0 + q1*q1 - q2*q2 - q3*q3)*mn + 2.0*(q1*q3 - q0*q2)*md;
  H(3, 0) = 2.0*( q0*mn - q2*md);
 H(3, 1) = 2.0*( q1*mn + q3*md);
  H(3, 2) = 2.0*(-q2*mn - q0*md);
  H(3, 3) = 2.0*(-q3*mn + q1*md);
  H(3, 4) = 0.0;
  H(3, 5) = 0.0;
  H(3, 6) = 0.0;
  
  //z(4, 0) = 2.0*(q1*q2 - q0*q3)*mn + 2.0*(q2*q3 + q0*q1)*md;
  H(4, 0) = 2.0*(-q3*mn + q1*md);
  H(4, 1) = 2.0*( q2*mn + q0*md);
  H(4, 2) = 2.0*( q1*mn + q3*md);
  H(4, 3) = 2.0*(-q0*mn + q2*md);
  H(4, 4) = 0.0;
  H(4, 5) = 0.0;
  H(4, 6) = 0.0;
  
  //z(5, 0) = 2.0*(q1*q3 + q0*q2)*mn + (q0*q0 - q1*q1 - q2*q2 + q3*q3)*md;
  H(5, 0) = 2.0*( q2*mn + q0*md);
  H(5, 1) = 2.0*( q3*mn - q1*md);
  H(5, 2) = 2.0*( q0*mn - q2*md);
  H(5, 3) = 2.0*( q1*mn + q3*md);
  H(5, 4) = 0.0;
  H(5, 5) = 0.0;
  H(5, 6) = 0.0;
  
  return 0;
}

//Extended Kalman Filter
uint8_t ekf( Matrix<float, 7, 1> &xe,
             Matrix<float, 7, 1> &xp,
             Matrix<float, 7, 7> &P,
             Matrix<float, 6, 1> z,
             Matrix<float, 3, 1> omega,
             Matrix<float, 3, 3> Q, 
             Matrix<float, 6, 6> R, 
             Matrix<float, 7, 3> G,
             Matrix<float, 3, 1> beta,
             float dt)
{
  Matrix<float, 7, 7> F;
  Matrix<float, 6, 7> H;
  Matrix<float, 6, 6> Den;
  Matrix<float, 6, 6> I6=MatrixXf::Identity(6,6);
  Matrix<float, 7, 6> K;
  Matrix<float, 6, 1> zbar;
  
  //Update
  H_jacobian(H, xp, GRAV, MN, MD);
  Den = H * P * H.transpose() + R;
  //PartialPivLU< Matrix<float, 6, 6> > dec(Den);
  //Den = dec.solve(I6);
  K = P * H.transpose() * Den.inverse();
  observation_equation(xp, zbar,GRAV, MN, MD);
  xe = xp + K*(z - zbar);
  P = P - K*H*P;

  //Predict
  state_equation(xe, omega, beta, dt, xp);
  F_jacobian(F, xe, omega, beta, dt);
  P = F*P*F.transpose() + G*Q*G.transpose();

  return 0;
}


float Phl(Matrix<float, 7, 1>x){
        double e23,e33;
        double phl=0;
        float qzero= x(0, 0);
        float qone = x(1, 0);
        float qtwo = x(2, 0);
        float qthree = x(3, 0);
        e23=2*(qtwo*qthree+qzero*qone);
        e33=qzero*qzero-qone*qone-qtwo*qtwo+qthree*qthree;
        phl=atan2(e23,e33);
        return phl;
}
float Theta(Matrix<float, 7, 1>x){
        double e13,e23,e33;
        double theta=0;
        float qzero= x(0, 0);
        float qone = x(1, 0);
        float qtwo = x(2, 0);
        float qthree = x(3, 0);
        e13=2*(qone*qthree-qzero*qtwo);
        e23=2*(qtwo*qthree+qzero*qone);
        e33=qzero*qzero-qone*qone-qtwo*qtwo+qthree*qthree;
        theta=atan2(-e13,sqrt(e23*e23+e33*e33));
        return theta;
}
float Psi(Matrix<float, 7, 1>x){
        float qzero= x(0, 0);
        float qone = x(1, 0);
        float qtwo = x(2, 0);
        float qthree = x(3, 0);
        double e12,e11;
        double psi=0;
        e11= (qzero*qzero+qone*qone-qtwo*qtwo-qthree*qthree);
        e12= 2*(qone*qtwo+qzero*qthree);
        psi=atan2(e12,e11);
        return psi;
}





void imu_mag_init(void)
{
  /* Initialize platform specific hardware */
  platform_init(  &Ins_bus,/* INS spi bus & cs pin */
                  &Mag_bus,/* MAG spi bus & cs pin */
                  &Imu_h,  /* IMU Handle */
                  &Mag_h,  /* MAG Handle */
                  1000*1000, /* SPI Frequency */
                  PIN_MISO,/* MISO Pin number */
                  PIN_SCK, /* SCK  Pin number */
                  PIN_MOSI,/* MOSI Pin number */
                  PIN_CSAG,/* CSAG Pin number */
                  PIN_CSM  /* CSM  Pin number */
  );
  /* Wait sensor boot time */
  platform_delay(BOOT_TIME);
  sleep_ms(1000);
  /* Check device ID */
  lsm9ds1_dev_id_get(&Mag_h, &Imu_h, &whoamI);
  if (whoamI.imu != LSM9DS1_IMU_ID || whoamI.mag != LSM9DS1_MAG_ID)     {
    while (1) {
      /* manage here device not found */
      printf("Device not found !\n");      sleep_ms(1000);
   }
}
 /* Restore default configuration */
 lsm9ds1_dev_reset_set(&Mag_h, &Imu_h, PROPERTY_ENABLE);
 do {
  lsm9ds1_dev_reset_get(&Mag_h, &Imu_h, &rst);
 } while (rst);
 /* Enable Block Data Update */
 lsm9ds1_block_data_update_set(&Mag_h, &Imu_h,
                               PROPERTY_ENABLE);
 /* Set full scale */
 lsm9ds1_xl_full_scale_set(&Imu_h, LSM9DS1_4g);
 lsm9ds1_gy_full_scale_set(&Imu_h, LSM9DS1_2000dps);
 lsm9ds1_mag_full_scale_set(&Mag_h, LSM9DS1_16Ga);
  /* Configure filtering chain - See datasheet for filtering chain     details */
  /* Accelerometer filtering chain */
  lsm9ds1_xl_filter_aalias_bandwidth_set(&Imu_h, LSM9DS1_AUTO);
  lsm9ds1_xl_filter_lp_bandwidth_set(&Imu_h,
                                     LSM9DS1_LP_ODR_DIV_50);
  lsm9ds1_xl_filter_out_path_set(&Imu_h, LSM9DS1_LP_OUT);
  /* Gyroscope filtering chain */
  lsm9ds1_gy_filter_lp_bandwidth_set(&Imu_h,
                                     LSM9DS1_LP_ULTRA_LIGHT);
  lsm9ds1_gy_filter_hp_bandwidth_set(&Imu_h, LSM9DS1_HP_MEDIUM);
  lsm9ds1_gy_filter_out_path_set(&Imu_h,
                                 LSM9DS1_LPF1_HPF_LPF2_OUT);
  /* Set Output Data Rate / Power mode */
  lsm9ds1_imu_data_rate_set(&Imu_h, LSM9DS1_IMU_476Hz);
  lsm9ds1_mag_data_rate_set(&Mag_h, LSM9DS1_MAG_MP_560Hz);
}
void imu_mag_data_read(void)
{
  /* Read device status register */
  lsm9ds1_dev_status_get(&Mag_h, &Imu_h, &reg);
  if ( reg.status_imu.xlda && reg.status_imu.gda ) {
    /* Read imu data */
    memset(data_raw_acceleration, 0x00, 3 * sizeof(int16_t));
    memset(data_raw_angular_rate, 0x00, 3 * sizeof(int16_t));
   lsm9ds1_acceleration_raw_get(&Imu_h,
                                 data_raw_acceleration);
   lsm9ds1_angular_rate_raw_get(&Imu_h,
                                data_raw_angular_rate);
   acceleration_mg[0] = lsm9ds1_from_fs4g_to_mg(
                          data_raw_acceleration[0]);
   acceleration_mg[1] = lsm9ds1_from_fs4g_to_mg(
                          data_raw_acceleration[1]);
   acceleration_mg[2] = lsm9ds1_from_fs4g_to_mg(
                          data_raw_acceleration[2]);
   angular_rate_mdps[0] = lsm9ds1_from_fs2000dps_to_mdps(
                           data_raw_angular_rate[0]);
   angular_rate_mdps[1] = lsm9ds1_from_fs2000dps_to_mdps(
                           data_raw_angular_rate[1]);
   angular_rate_mdps[2] = lsm9ds1_from_fs2000dps_to_mdps(
                            data_raw_angular_rate[2]);
  }
  if ( reg.status_mag.zyxda ) {
    /* Read magnetometer data */
    memset(data_raw_magnetic_field, 0x00, 3 * sizeof(int16_t));
    lsm9ds1_magnetic_raw_get(&Mag_h, data_raw_magnetic_field);
    magnetic_field_mgauss[0] = lsm9ds1_from_fs16gauss_to_mG(
                                 data_raw_magnetic_field[0]);
    magnetic_field_mgauss[1] = lsm9ds1_from_fs16gauss_to_mG(
                                 data_raw_magnetic_field[1]);
    magnetic_field_mgauss[2] = lsm9ds1_from_fs16gauss_to_mG(
                                 data_raw_magnetic_field[2]);
  }

}












#if 0
int main(void)
{
  Matrix<float, 7 ,1> xp = MatrixXf::Zero(7,1);
  Matrix<float, 7 ,1> xe = MatrixXf::Zero(7,1);
  Matrix<float, 7 ,1> x_sim = MatrixXf::Zero(7,1);
  Matrix<float, 7 ,7> P = MatrixXf::Identity(7,7);
  Matrix<float, 6 ,1> z = MatrixXf::Zero(6,1);
  Matrix<float, 6 ,1> z_sim = MatrixXf::Zero(6,1);
  Matrix<float, 6 ,1> z_noise = MatrixXf::Zero(6,1);
  Matrix<float, 3, 1> omega_m = MatrixXf::Zero(3, 1);
  Matrix<float, 3, 1> omega_sim;
  Matrix<float, 3, 1> domega;
  Matrix<float, 3, 1> domega_sim;
  Matrix<float, 3, 3> Q = MatrixXf::Identity(3, 3)*1;
  Matrix<float, 6, 6> R = MatrixXf::Identity(6, 6)*1;
  Matrix<float, 7 ,3> G;
  Matrix<float, 3 ,1> beta;
  float t=0.0,dt=0.01;
  uint64_t s_time=0,e_time=0,d_time=0;
  double phl,theta,psi;
  short i,waittime=5;
  float p_com[10]={0.1*PI, 0, 0.01*PI, 0, -0.01*PI, 0, 0.02*PI, 0, -0.05*PI, 0};
  float q_com[10]={0.1*PI, 0, 0      , 0, -0.01*PI, 0, 0.02*PI, 0, -0.02*PI, 0};
  float r_com[10]={0.1*PI, 0,-0.02*PI, 0,      -PI, 0, 0.02*PI, 0,  0.02*PI, 0};
  float endtime=10000.0;
  float control_period=5.0;
  float control_time=5.0;
  int counter=0;
  int sample=1;
  int control_counter=0;
  int control_counter_max=0;
  double pi=3.14159265358;
  float ax,ay,az,wp,wq,wr,mx,my,mz; 
  float p11=-0.60526553,p12=0.79021444,p13=-0.09599364,p21=0.78892428,p22=0.61155945,p23=0.05994594,p31=-0.10607597,p32=0.0394485,p33=0.99357521;
  float r1=-4.96008102e-6,r2=-4.60371715e-6,r3=-4.17649591e-6;
  float w;
  float dmx,dmy,dmz,dxx,dxy,dxz;
  float ddxx, ddxy,ddxz;
  float dddxx,dddxy,dddxz;
  float ddddxx,ddddxy,ddddxz;
  //Variable Initalize
  xe << 1.0, 0.0, 0.0, 0.0, -0.078, 0.0016, 0.00063;
  xp =xe;
  x_sim << 1.0, 0.0, 0.0, 0.0, 0.01, 0.02, 0.03;
  observation_equation(x_sim, z_sim, GRAV, MN, MD);

  G <<  0.0,0.0,0.0, 
        0.0,0.0,0.0, 
        0.0,0.0,0.0, 
        0.0,0.0,0.0, 
        1.0,0.0,0.0, 
        0.0,1.0,0.0, 
        0.0,0.0,1.0;

  beta << 0.003, 0.003, 0.003;

  P <<  1,0,0,0,0,0,0,  
        0,1,0,0,0,0,0,
        0,0,1,0,0,0,0,  
        0,0,0,1,0,0,0, 
        0,0,0,0,1,0,0,  
        0,0,0,0,0,1,0,  
        0,0,0,0,0,0,1;

  Q << 7.34944e-6,0,0,
       0,6.861e-6,0,
       0,0,5.195e-6;

  R << 3.608e-6,0,0,0,0,0,
       0,6.261e-6,0,0,0,0,
       0,0,1.889e-5,0,0,0,
       0,0,0,1.0,0,0,
       0,0,0,0,1.0,0,
       0,0,0,0,0,1.0;
  
  //Initilize Console Input&Output
  stdio_init_all();
  imu_mag_init();

  

#if 1
  //Start up wait for Pico
  for (i=0;i<waittime;i++)
  {
    printf("#Please wait %d[s] ! Random number test:%f\n",waittime-i, norm(mt) );
    sleep_ms(1000);
  }
  printf("#Start Kalman Filter\n");
#endif    
  while(t<endtime)
  {
    s_time=time_us_64();
    //Control
    if(t>control_time)
    {
      control_time = control_time + control_period;
      control_counter++;
      if(control_counter>control_counter_max)control_counter=0;
    }
    ax=   -acceleration_mg[0]*1/1000*GRAV;
    ay=   -acceleration_mg[1]*1/1000*GRAV;
    az=    acceleration_mg[2]*1/1000*GRAV;
    wp=    angular_rate_mdps[0]/1000*pi/180;
    wq=    angular_rate_mdps[1]/1000*pi/180;
    wr=   -angular_rate_mdps[2]/1000*pi/180;
    dmx=  -(magnetic_field_mgauss[0]-310);
    dmy=   (magnetic_field_mgauss[1]-10);
    dmz=  -(magnetic_field_mgauss[2]);




    //校正作業
    dxx=p11*dmx+p21*dmy+p31*dmz;
    dxy=p12*dmx+p22*dmy+p32*dmz;
    dxz=p13*dmx+p23*dmy+p33*dmz;

    ddxx=dxx-22.760831749415342;
    ddxy=dxy-19.734355196006327;
    ddxz=dxz-141.33565570453044;


    w=-1.087745370038146;
  
    dddxx=ddxx*0.0020572671658147883;
    dddxy=ddxy*0.0021354074993493823;
    dddxz=ddxz*0.0019594870993397107;

    mx=p11*dddxx+p22*dddxy+p13*dddxz;
    my=p21*dddxx+p22*dddxy+p23*dddxz;
    mz=p31*dddxx+p32*dddxy+p33*dddxz;
 
    omega_m <<wp,wq,wr;
    z       <<ax,ay,az,mx,my,mz;//ここに入れる
    //--Begin Extended Kalman Filter--
    ekf(xp, xe, P, z, omega_m, Q, R, G*dt, beta, dt);
   // e_time=time_us_64();
    //--End   Extended Kalman Filter--

    //Result output
    if(counter%sample==0)
    {
		imu_mag_data_read();
                phl=Phl(xe);
                theta=Theta(xe);
                psi=Psi(xe);
              //  printf("%9.6f %9.6f %9.6f %9.6f  ",t,phl,theta,psi);
      printf("%9.2f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f\n", 
               t, xe(0,0), xe(1,0), xe(2,0),xe(3,0), xe(4,0), xe(5,0),xe(6,0)
     // printf("%9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f\n",ax,ay,az,wp,wq,wr,mx,my,mz);       
	 x_sim(0,0), x_sim(1,0), x_sim(2,0), x_sim(3,0),
	        x_sim(4,0), x_sim(5,0), x_sim(6,0),
                p_com[control_counter], q_com[control_counter], r_com[control_counter],
                e_time-s_time);  
      printf( "%9.6fIMU-[mg]:\t%9.6f\t%9.6f\t%9.6f\t[mdps]:\t%9.6f\t%9.6f\t%9.6f\t",t
              ,acceleration_mg[0]*1/1000*GRAV, acceleration_mg[1]*1/1000*GRAV, acceleration_mg[2]*1/1000*GRAV*(-1),
                angular_rate_mdps[0]/1000*pi/180, angular_rate_mdps[1]/1000*pi/180, angular_rate_mdps[2]/1000*pi/180*(-1));
      printf( "MAG-[mG]:\t%9.6f\t%9.6f\t%9.6f\r\n"
               ,(magnetic_field_mgauss[0]-310)*(-1), magnetic_field_mgauss[1]-10,magnetic_field_mgauss[2]*(-1));
                //Control
               // domega<<xe(4,0), xe(5,0), xe(6,0);
              	//printf("%f,%f,%f\n",angular_rate_mdps[0], angular_rate_mdps[1], angular_rate_mdps[2]*(-1));*/
                //omega=omega_sim+domega;
                //Simulation
               // observation_equation(xe, z, GRAV, MN, MD);
                //z=z_sim;
                //rk4(quatdot, t, dt, quat_sim, omega_sim);
                //x_sim << quat_sim(0,0), quat_sim(1,0), quat_sim(2,0), quat_sim(3,0), 0,0,0;
                //t=t+dt;
                //Begin Extended Kalman Filter
                // s_time=time_us_64();
                //ekf(xe, P, z, omega_m, beta, dt);
      printf("%9.2f %9.6f %9.6f %9.6f \n",t,mx,my,mz);

     }  
    counter++;
    
    
    t=t+dt;
    while (time_us_64()-s_time<10000);

  }
  return 0;
}
#endif
