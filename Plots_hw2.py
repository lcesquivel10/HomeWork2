import numpy as np
import matplotlib.pyplot as plt

for i in range(1,4): ## hacer un recorrido por los datos obtenidos
    data_euler = np.loadtxt('data_euler_dt_'+str(i)+'e-3.dat') ## cargar los datos y asignarles nombre con el dt
    data_leap = np.loadtxt('data_leap_dt_'+str(i)+'e-3.dat')
    data_runge = np.loadtxt('data_runge_dt_'+str(i)+'e-3.dat')

    t = data_euler[:,0]  # primera columna 

    if i == 1:
        t_max = t[len(t)-1]

    max_index=0
    for j in range(len(t)):
        if t[j] >= t_max:
            max_index = j
            break

    x_euler = data_euler[:,1] ## asignar a las columnas su valor para cada metodo 
    y_euler = data_euler[:,2]
    vx_euler = data_euler[:,3]
    vy_euler = data_euler[:,4]

    x_leap = data_leap[:,1]
    y_leap = data_leap[:,2]
    vx_leap = data_leap[:,3]
    vy_leap = data_leap[:,4]

    x_runge = data_runge[:,1]
    y_runge = data_runge[:,2]
    vx_runge = data_runge[:,3]
    vy_runge = data_runge[:,4]