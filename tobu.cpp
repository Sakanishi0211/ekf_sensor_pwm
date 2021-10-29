#include "ekf_sensor.hpp"
#include <stdio.h>
#include "pico/stdlib.h"
#include "hardware/pwm.h"
#include "hardware/uart.h"
#include "hardware/irq.h"
#include "pico/binary_info.h"
#include "hardware/spi.h"
#include "lsm9ds1_reg.h"

#define UART_ID uart0
#define BAUD_RATE 100000      //115200
#define DATA_BITS 8
#define STOP_BITS 2             //1
#define PARITY    UART_PARITY_EVEN
#define DUTYMIN 1330
#define DUTYMAX 1520
#define CH1MAX 1748
#define CH1MIN 332
#define CH2MAX 1720
#define CH2MIN 324
#define CH3MAX 1724
#define CH3MIN 323
#define CH4MAX 1715
#define CH4MIN 340
#define CH5MAX 1904
#define CH5MIN 144
#define CH6MAX 1904
#define CH6MIN 144
//０番と1番ピンに接続
#define UART_TX_PIN 0
#define UART_RX_PIN 1
/* Private macro -------------------------------------------------------------*/

extern uint8_t slice_num[2];
extern uint8_t pwm_settei();

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
 //pwm_uart
  float duty_rr,duty_rl,duty_fl,duty_fr;
  float seigyo=0.5;  
  pwm_settei();
 // elevator();
  sleep_ms(2000);
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
    printf("#Please wait %d[s] ! \n",waittime-i);
    sleep_ms(1000);
  }
  printf("#Start Kalman Filter\n");
#endif    
  while(t<endtime)
  {
    s_time=time_us_64();
    sleep_ms(10);
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
     /* printf("%9.2f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f\n", 
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
               ,(magnetic_field_mgauss[0]-310)*(-1), magnetic_field_mgauss[1]-10,magnetic_field_mgauss[2]*(-1));*/
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
