#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

// defino mis funciones 

double f_1(double t, double x, double y, double vx, double vy){
    double G = 4.*M_PI*M_PI;
    double M_Sol = 1.;
    return -G*M_Sol*(x/sqrt((x*x + y*y)))/(x*x + y*y);
}

double f_2(double t, double x, double y, double vx, double vy){
    double G = 4.*M_PI*M_PI;
    double M_Sol = 1.;
    return -G*M_Sol*(y/sqrt((x*x + y*y)))/(x*x + y*y);
}

double g_1(double t, double x, double y, double vx, double vy)
{return vx;}
double g_2(double t, double x, double y, double vx, double vy)
{return vy;}

/// guia metodos https://github.com/ComputoCienciasUniandes/MetodosComputacionales/blob/master/secciones/10.ODEs/ODE.ipynb

///// METODO DE EULER //////

// para la velocidad
double euler_vx(double t_n, double x_n, double y_n, double vx_n, double vy_n, double h){
    return vx_n + h*f_1(t_n,x_n,y_n,vx_n,vy_n);
}
double euler_vy(double t_n, double x_n, double y_n, double vx_n, double vy_n, double h){
    return vy_n + h*f_2(t_n,x_n,y_n,vx_n,vy_n);
}

/// para la posicion
double euler_x(double t_n, double x_n, double y_n, double vx_n, double vy_n, double h){
    return x_n + h*g_1(t_n,x_n,y_n,vx_n,vy_n);
}
double euler_y(double t_n, double x_n, double y_n, double vx_n, double vy_n, double h){
    return y_n + h*g_2(t_n,x_n,y_n,vx_n,vy_n);
    
//// LEAP FROG /////
    
double leap_vx(double t_n, double x_n, double y_n, double vx_n, double vy_n, double vx_n_2,double h){
    return vx_n_2 + 2*h*f_1(t_n,x_n,y_n,vx_n,vy_n);
}
double leap_vy(double t_n, double x_n, double y_n, double vx_n, double vy_n, double vy_n_2,double h){
    return vy_n_2 + 2*h*f_2(t_n,x_n,y_n,vx_n,vy_n);
}
double leap_x(double t_n, double x_n, double y_n, double vx_n, double vy_n, double x_n_2,double h){
    return x_n_2 + 2*h*g_1(t_n,x_n,y_n,vx_n,vy_n);
}
double leap_y(double t_n, double x_n, double y_n, double vx_n, double vy_n,  double y_n_2,double h){
    return y_n_2 + 2*h*g_2(t_n,x_n,y_n,vx_n,vy_n);
}
////////////RUNGE- KUTTA 4 ////////
    
    
double runge_kutta_vx(double t_n, double x_n, double y_n, double vx_n, double vy_n, double h){
    double k11 = h*g_1(t_n,x_n,y_n,vx_n,vy_n);
    double k12 = h*g_2(t_n,x_n,y_n,vx_n,vy_n);
    double k13 = h*f_1(t_n,x_n,y_n,vx_n,vy_n);
    double k14 = h*f_2(t_n,x_n,y_n,vx_n,vy_n);

    double k21 = h*g_1(t_n + h*0.5,x_n + k11*0.5, y_n + k12*0.5, vx_n + k13*0.5, vy_n + k14*0.5);
    double k22 = h*g_2(t_n + h*0.5,x_n + k11*0.5, y_n + k12*0.5, vx_n + k13*0.5, vy_n + k14*0.5); 
    double k23 = h*f_1(t_n + h*0.5,x_n + k11*0.5, y_n + k12*0.5, vx_n + k13*0.5, vy_n + k14*0.5); 
    double k24 = h*f_2(t_n + h*0.5,x_n + k11*0.5, y_n + k12*0.5, vx_n + k13*0.5, vy_n + k14*0.5); 

    double k31 = h*g_1(t_n + h*0.5,x_n + k21*0.5, y_n + k22*0.5, vx_n + k23*0.5, vy_n + k24*0.5);
    double k32 = h*g_2(t_n + h*0.5,x_n + k21*0.5, y_n + k22*0.5, vx_n + k23*0.5, vy_n + k24*0.5); 
    double k33 = h*f_1(t_n + h*0.5,x_n + k21*0.5, y_n + k22*0.5, vx_n + k23*0.5, vy_n + k24*0.5); 
    double k34 = h*f_2(t_n + h*0.5,x_n + k21*0.5, y_n + k22*0.5, vx_n + k23*0.5, vy_n + k24*0.5); 

    double k41 = h*g_1(t_n + h*0.5,x_n + k31, y_n + k32, vx_n + k33, vy_n + k34);
    double k42 = h*g_2(t_n + h*0.5,x_n + k31, y_n + k32, vx_n + k33, vy_n + k34); 
    double k43 = h*f_1(t_n + h*0.5,x_n + k31, y_n + k32, vx_n + k33, vy_n + k34); 
    double k44 = h*f_2(t_n + h*0.5,x_n + k31, y_n + k32, vx_n + k33, vy_n + k34); 

    return vx_n + (1./6.)*(k13 + 2.*k23 + 2.*k33 + k43);
}
double runge_kutta_vy(double t_n, double x_n, double y_n, double vx_n, double vy_n, double h){
    double k11 = h*g_1(t_n,x_n,y_n,vx_n,vy_n);
    double k12 = h*g_2(t_n,x_n,y_n,vx_n,vy_n);
    double k13 = h*f_1(t_n,x_n,y_n,vx_n,vy_n);
    double k14 = h*f_2(t_n,x_n,y_n,vx_n,vy_n);

    double k21 = h*g_1(t_n + h*0.5,x_n + k11*0.5, y_n + k12*0.5, vx_n + k13*0.5, vy_n + k14*0.5);
    double k22 = h*g_2(t_n + h*0.5,x_n + k11*0.5, y_n + k12*0.5, vx_n + k13*0.5, vy_n + k14*0.5); 
    double k23 = h*f_1(t_n + h*0.5,x_n + k11*0.5, y_n + k12*0.5, vx_n + k13*0.5, vy_n + k14*0.5); 
    double k24 = h*f_2(t_n + h*0.5,x_n + k11*0.5, y_n + k12*0.5, vx_n + k13*0.5, vy_n + k14*0.5); 

    double k31 = h*g_1(t_n + h*0.5,x_n + k21*0.5, y_n + k22*0.5, vx_n + k23*0.5, vy_n + k24*0.5);
    double k32 = h*g_2(t_n + h*0.5,x_n + k21*0.5, y_n + k22*0.5, vx_n + k23*0.5, vy_n + k24*0.5); 
    double k33 = h*f_1(t_n + h*0.5,x_n + k21*0.5, y_n + k22*0.5, vx_n + k23*0.5, vy_n + k24*0.5); 
    double k34 = h*f_2(t_n + h*0.5,x_n + k21*0.5, y_n + k22*0.5, vx_n + k23*0.5, vy_n + k24*0.5); 

    double k41 = h*g_1(t_n + h*0.5,x_n + k31, y_n + k32, vx_n + k33, vy_n + k34);
    double k42 = h*g_2(t_n + h*0.5,x_n + k31, y_n + k32, vx_n + k33, vy_n + k34); 
    double k43 = h*f_1(t_n + h*0.5,x_n + k31, y_n + k32, vx_n + k33, vy_n + k34); 
    double k44 = h*f_2(t_n + h*0.5,x_n + k31, y_n + k32, vx_n + k33, vy_n + k34); 

    return vy_n + (1./6.)*(k14 + 2.*k24 + 2.*k34 + k44);
}
double runge_kutta_x(double t_n, double x_n, double y_n, double vx_n, double vy_n, double h){
    double k11 = h*g_1(t_n,x_n,y_n,vx_n,vy_n);
    double k12 = h*g_2(t_n,x_n,y_n,vx_n,vy_n);
    double k13 = h*f_1(t_n,x_n,y_n,vx_n,vy_n);
    double k14 = h*f_2(t_n,x_n,y_n,vx_n,vy_n);

    double k21 = h*g_1(t_n + h*0.5,x_n + k11*0.5, y_n + k12*0.5, vx_n + k13*0.5, vy_n + k14*0.5);
    double k22 = h*g_2(t_n + h*0.5,x_n + k11*0.5, y_n + k12*0.5, vx_n + k13*0.5, vy_n + k14*0.5); 
    double k23 = h*f_1(t_n + h*0.5,x_n + k11*0.5, y_n + k12*0.5, vx_n + k13*0.5, vy_n + k14*0.5); 
    double k24 = h*f_2(t_n + h*0.5,x_n + k11*0.5, y_n + k12*0.5, vx_n + k13*0.5, vy_n + k14*0.5); 

    double k31 = h*g_1(t_n + h*0.5,x_n + k21*0.5, y_n + k22*0.5, vx_n + k23*0.5, vy_n + k24*0.5);
    double k32 = h*g_2(t_n + h*0.5,x_n + k21*0.5, y_n + k22*0.5, vx_n + k23*0.5, vy_n + k24*0.5); 
    double k33 = h*f_1(t_n + h*0.5,x_n + k21*0.5, y_n + k22*0.5, vx_n + k23*0.5, vy_n + k24*0.5); 
    double k34 = h*f_2(t_n + h*0.5,x_n + k21*0.5, y_n + k22*0.5, vx_n + k23*0.5, vy_n + k24*0.5); 

    double k41 = h*g_1(t_n + h*0.5,x_n + k31, y_n + k32, vx_n + k33, vy_n + k34);
    double k42 = h*g_2(t_n + h*0.5,x_n + k31, y_n + k32, vx_n + k33, vy_n + k34); 
    double k43 = h*f_1(t_n + h*0.5,x_n + k31, y_n + k32, vx_n + k33, vy_n + k34); 
    double k44 = h*f_2(t_n + h*0.5,x_n + k31, y_n + k32, vx_n + k33, vy_n + k34); 

    return x_n + (1./6.)*(k11 + 2.*k21 + 2.*k31 + k41);
}
double runge_kutta_y(double t_n, double x_n, double y_n, double vx_n, double vy_n, double h){
    double k11 = h*g_1(t_n,x_n,y_n,vx_n,vy_n);
    double k12 = h*g_2(t_n,x_n,y_n,vx_n,vy_n);
    double k13 = h*f_1(t_n,x_n,y_n,vx_n,vy_n);
    double k14 = h*f_2(t_n,x_n,y_n,vx_n,vy_n);

    double k21 = h*g_1(t_n + h*0.5,x_n + k11*0.5, y_n + k12*0.5, vx_n + k13*0.5, vy_n + k14*0.5);
    double k22 = h*g_2(t_n + h*0.5,x_n + k11*0.5, y_n + k12*0.5, vx_n + k13*0.5, vy_n + k14*0.5); 
    double k23 = h*f_1(t_n + h*0.5,x_n + k11*0.5, y_n + k12*0.5, vx_n + k13*0.5, vy_n + k14*0.5); 
    double k24 = h*f_2(t_n + h*0.5,x_n + k11*0.5, y_n + k12*0.5, vx_n + k13*0.5, vy_n + k14*0.5); 

    double k31 = h*g_1(t_n + h*0.5,x_n + k21*0.5, y_n + k22*0.5, vx_n + k23*0.5, vy_n + k24*0.5);
    double k32 = h*g_2(t_n + h*0.5,x_n + k21*0.5, y_n + k22*0.5, vx_n + k23*0.5, vy_n + k24*0.5); 
    double k33 = h*f_1(t_n + h*0.5,x_n + k21*0.5, y_n + k22*0.5, vx_n + k23*0.5, vy_n + k24*0.5); 
    double k34 = h*f_2(t_n + h*0.5,x_n + k21*0.5, y_n + k22*0.5, vx_n + k23*0.5, vy_n + k24*0.5); 

    double k41 = h*g_1(t_n + h*0.5,x_n + k31, y_n + k32, vx_n + k33, vy_n + k34);
    double k42 = h*g_2(t_n + h*0.5,x_n + k31, y_n + k32, vx_n + k33, vy_n + k34); 
    double k43 = h*f_1(t_n + h*0.5,x_n + k31, y_n + k32, vx_n + k33, vy_n + k34); 
    double k44 = h*f_2(t_n + h*0.5,x_n + k31, y_n + k32, vx_n + k33, vy_n + k34); 

    return y_n + (1./6.)*(k12 + 2.*k22 + 2.*k32 + k42);
}

int main(void){


    int N = 20*1000;

    double t0 = 0;
    double tf = 1;

    double dt;

    double x_euler[N];
    double y_euler[N];
    double vx_euler[N];
    double vy_euler[N];

    double x_leap[N];
    double y_leap[N];
    double vx_leap[N];
    double vy_leap[N];

    double x_runge[N];
    double y_runge[N];
    double vx_runge[N];
    double vy_runge[N];

    x_euler[0] = 0.1163;
    y_euler[0] = 0.9772;
    vx_euler[0] = -6.35;
    vy_euler[0] = 0.606;

    x_leap[0] = 0.1163;
    y_leap[0] = 0.9772;
    vx_leap[0] = -6.35;
    vy_leap[0] = 0.606;
    
    x_runge[0] = 0.1163;
    y_runge[0] = 0.9772;
    vx_runge[0] = -6.35;
    vy_runge[0] = 0.606;

    for(int k = 1; k < 4; k++){

        dt = k*1e-3;

        x_leap[1] = x_leap[0] + dt*g_1(t0,x_leap[0],y_leap[0],vx_leap[0],vy_leap[0]);
        y_leap[1] = y_leap[0] + dt*g_2(t0,x_leap[0],y_leap[0],vx_leap[0],vy_leap[0]);
        vx_leap[1] = vx_leap[0] + dt*f_1(t0,x_leap[0],y_leap[0],vx_leap[0],vy_leap[0]);
        vy_leap[1] = vy_leap[0] + dt*f_2(t0,x_leap[0],y_leap[0],vx_leap[0],vy_leap[0]);

        for(int i = 1; i < N; i++){
            vx_euler[i] = euler_vx(t0+i*dt,x_euler[i-1],y_euler[i-1],vx_euler[i-1],vy_euler[i-1],dt);
            vy_euler[i] = euler_vy(t0+i*dt,x_euler[i-1],y_euler[i-1],vx_euler[i-1],vy_euler[i-1],dt);
            x_euler[i] = euler_x(t0+i*dt,x_euler[i-1],y_euler[i-1],vx_euler[i-1],vy_euler[i-1],dt);
            y_euler[i] = euler_y(t0+i*dt,x_euler[i-1],y_euler[i-1],vx_euler[i-1],vy_euler[i-1],dt);
        }

        for(int i = 2; i < N; i++){
            vx_leap[i] = leap_vx(t0+i*dt,x_leap[i-1],y_leap[i-1],vx_leap[i-1],vy_leap[i-1],vx_leap[i-2],dt);
            vy_leap[i] = leap_vy(t0+i*dt,x_leap[i-1],y_leap[i-1],vx_leap[i-1],vy_leap[i-1],vy_leap[i-2],dt);
            x_leap[i] = leap_x(t0+i*dt,x_leap[i-1],y_leap[i-1],vx_leap[i-1],vy_leap[i-1],x_leap[i-2],dt);
            y_leap[i] = leap_y(t0+i*dt,x_leap[i-1],y_leap[i-1],vx_leap[i-1],vy_leap[i-1],y_leap[i-2],dt);
        }

        for(int i = 1; i < N; i++){
            vx_runge[i] = runge_kutta_vx(t0+i*dt,x_runge[i-1],y_runge[i-1],vx_runge[i-1],vy_runge[i-1],dt);
            vy_runge[i] = runge_kutta_vy(t0+i*dt,x_runge[i-1],y_runge[i-1],vx_runge[i-1],vy_runge[i-1],dt);
            x_runge[i] = runge_kutta_x(t0+i*dt,x_runge[i-1],y_runge[i-1],vx_runge[i-1],vy_runge[i-1],dt);
            y_runge[i] = runge_kutta_y(t0+i*dt,x_runge[i-1],y_runge[i-1],vx_runge[i-1],vy_runge[i-1],dt);
        }

        ofstream data_euler("data_euler_dt_"+to_string(k)+"e-3.dat");    
        for(int i = 0; i < N; i++){
            data_euler << t0+i*dt << '\t' << x_euler[i] << '\t' << y_euler[i] << '\t' << vx_euler[i] << '\t' << vy_euler[i] << endl;
        }

        ofstream data_leap("data_leap_dt_"+to_string(k)+"e-3.dat");    
        for(int i = 0; i < N; i++){
            data_leap << t0+i*dt << '\t' << x_leap[i] << '\t' << y_leap[i] << '\t' << vx_euler[i] << '\t' << vy_euler[i] << endl;
        }

        ofstream data_runge("data_runge_dt_"+to_string(k)+"e-3.dat");    
        for(int i = 0; i < N; i++){
            data_runge << t0+i*dt << '\t' << x_runge[i] << '\t' << y_runge[i] << '\t' << vx_euler[i] << '\t' << vy_euler[i] << endl;
        }
    }



   
    return 0;
}