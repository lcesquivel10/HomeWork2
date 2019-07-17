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

//Euler
double euler_vx(double t_n, double x_n, double y_n, double vx_n, double vy_n, double h){
    return vx_n + h*f_1(t_n,x_n,y_n,vx_n,vy_n);
}
double euler_vy(double t_n, double x_n, double y_n, double vx_n, double vy_n, double h){
    return vy_n + h*f_2(t_n,x_n,y_n,vx_n,vy_n);
}
double euler_x(double t_n, double x_n, double y_n, double vx_n, double vy_n, double h){
    return x_n + h*g_1(t_n,x_n,y_n,vx_n,vy_n);
}
double euler_y(double t_n, double x_n, double y_n, double vx_n, double vy_n, double h){
    return y_n + h*g_2(t_n,x_n,y_n,vx_n,vy_n);