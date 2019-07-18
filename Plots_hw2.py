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
    
    energy_euler = 0.5*(vx_euler**2. + vy_euler**2.) ### formula para sacar energia 
    energy_leap = 0.5*(vx_leap**2. + vy_leap**2.)
    energy_runge = 0.5*(vx_runge**2. + vy_runge**2.)

    angular_euler = x_euler*vy_euler - y_euler*vx_euler ## para sacar el momentum 
    angular_leap = x_leap*vy_leap - y_leap*vx_leap
    angular_runge = x_runge*vy_runge - y_runge*vx_runge
    
######################################33 ORBITAS X vs Y
    plt.figure(figsize = (5,5))
    plt.plot(x_euler[:max_index],y_euler[:max_index], label = "Euler")
    plt.xlabel("X (UA)")
    plt.ylabel("Y (UA)")
    plt.title("Orbita")
    plt.legend()
    plt.savefig("Orbita_euler_dt_"+str(i)+"e-3.png")

    plt.figure(figsize = (5,5))
    plt.plot(x_leap[:max_index],y_leap[:max_index], label= "Leap Frog")
    plt.xlabel("X (UA)")
    plt.ylabel("Y (UA)")
    plt.title("Orbita")
    plt.legend()
    plt.savefig("Orbita_leap_dt_"+str(i)+"e-3.png")

    plt.figure(figsize = (5,5))
    plt.plot(x_runge[:max_index],y_runge[:max_index], label = "Runge Kutta")
    plt.xlabel("X (UA)")
    plt.ylabel("Y (UA)")
    plt.title("Orbita")
    plt.legend()
    plt.savefig("Orbita_runge_dt_"+str(i)+"e-3.png")



    
    