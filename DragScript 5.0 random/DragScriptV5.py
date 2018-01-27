from numpy import *
from matplotlib.pyplot import *
from astropy.io import ascii
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d

#style.use('classic')

def num_solver(thrust_profile, rocket_mass, motor_mass, propellant_mass, time_res, temp, burn_time, plots):
    
    
    rocket_time = ascii.read(thrust_profile)['time']
    thrust_curve = ascii.read(thrust_profile)['thrust']

    thrust_function = interp1d(rocket_time, thrust_curve)
    
    #constants
    g = 9.81 #m/s^2
    
    temps = [0., 20., 40., 60.]
    densities = [1.293, 1.205, 1.127, 1.067]
    
    temp_func = interp1d(temps, densities, kind = 'quadratic')
    rho = temp_func(temp) # kg/m^3
    

    #initial conditions
    h0 = 0. #m  (height)
    v0 = 0. #m/s (velocity)

    #rocket characteristics

    Wet_mass = motor_mass + rocket_mass + propellant_mass  #kg
    Dry_mass = motor_mass + rocket_mass #kg
    
    
    # 1D trajectory plot via Forward Euler Integration

    time_step = time_res #sec
    
    #Get the mass flow rate from the thrust curve
    t_dummy = 0     #dummy time variable
    total_thrust = 0
    while t_dummy<burn_time+1:
        total_thrust = total_thrust+thrust_function(t_dummy)
        t_dummy = t_dummy+time_step
        
    
    def mass_flow_rate(t):
        if t<=burn_time:
            return thrust_function(t)*(Wet_mass-Dry_mass)/total_thrust
        else:
            return 0
        
    def height_step(hi, vi):
        return  hi + vi * time_step 
    
    def velocity_step(hi, vi, mi, ti):
        
        area = 0.006032066874
        coeff_drag=0.44
        
        return vi + (-g - 0.5*air_pressure(hi, rho)*area*coeff_drag*vi*vi/mi + thrust_function(ti)/mi)*time_step
    
    def mass_step(mi, t): 
        return mi + (-mass_flow_rate(t))
    
    def air_pressure(h, p_i):
        R=8.31447       #Universal gas constant
        M=0.0289644     #Molar mass of air kg/mol
        T=temp+273.15   #temperature of the air at launch altitude in Kelvin
        return p_i*e**(-1*(g*M*h)/(R*T))

    #arrays to push results to

    h_arr = [h0]
    v_arr = [v0]
    a_arr = [0]
    m_arr = [Wet_mass]
    time_array = [0]

    #calculation
    i = 0
    
    
    #ascent
   
    while (time_array[i]<burn_time): 
    
        ti = (time_step)*(i)
        time_array.append(ti)
    
        hi = height_step(h_arr[i], v_arr[i])
        h_arr.append(hi)
    
        vi = velocity_step(h_arr[i], v_arr[i], m_arr[i], time_array[i])
        v_arr.append(vi)
    
        mi = mass_step(m_arr[i], ti)
        m_arr.append(mi)
    
        i = i+1
    
    # coast
    c = i #mark transition to coast

    while (v_arr[i]>0): 
    
        ti = (time_step)*(i)
        time_array.append(ti)
    
        hi = height_step(h_arr[i], v_arr[i])
        h_arr.append(hi)
    
        vi = velocity_step(h_arr[i], v_arr[i], m_arr[i], time_array[i])
        v_arr.append(vi)
    
        mi = Dry_mass
        m_arr.append(mi)
    
        i = i+1
       
    #descent

    d = i #mark transition to descent
    
    t_arr = linspace(0,time_step*(len(h_arr)+1),len(h_arr))  #check

    if plots:
        
        #############a function that differentiates the data##############3
        
        def derive_accel(t, V):
            accel_array = zeros(len(V)-1)
            for i in range(0, len(accel_array)):
                accel_array[i] = (V[i+1] - V[i])/(t[i+1]-t[i])
            return accel_array
   
        figure(figsize=(10, 15))
        subplot(3,1,1)
        plot(t_arr,array(h_arr), '--', color = 'black', label = 'Trajectory')
        grid()
        title('Trajectory')
        xlabel('Time [s]')
        ylabel('Height [m]')

        plot(t_arr[c],h_arr[c], 'o', color = 'black', label = 'Main Engine Cutoff')
        print 'MECO at', round(t_arr[c], 2), 'sec, at' , round(h_arr[c], 2), 'm'

        plot(t_arr[d],h_arr[d], 'o', color = 'tomato', label = 'Apogee')
        print 'Apogee at', round(t_arr[d], 2), 'sec, at' , round(h_arr[d], 2), 'm'
        legend(loc = 'best')
        xlim(0, time_array[-1]+0.1)
        ylim(-10, max(array(h_arr))+50)

        subplot(3,1,2)
        plot(t_arr, array(v_arr), '--', color = 'black', label = 'Velocity')
        grid()
        title('Velocity')
        xlabel('Time [s]')
        ylabel('Velocity [m/s]')

        #plot(t_arr[c],v_arr[c], 'o', color = 'black', label = 'Main Engine Cutoff')
        print 'MECO at', round(v_arr[c], 2), 'm/s'

        plot(t_arr[d], v_arr[d], 'o', color = 'tomato', label = 'Apogee')
        xlim(0, time_array[-1]+0.1)
        ylim(min(v_arr)-5, v_arr[c] + 40)
        legend(loc = 'best')
        
        
        
        subplot(3,1,3)
        accel_array = derive_accel(t_arr, array(v_arr))
        plot(t_arr[0:-1], accel_array, '-', color = 'black', label = 'Acceleration')
        #plot(t_arr[c], accel_array[c], 'o', color = 'black', label = 'Main Engine Cutoff')
        #plot(t_arr[d], accel_array[d], 'o', color = 'tomato', label = 'Apogee')
        xlabel('Time [s]')
        ylabel('Acceleration [m/sec^2]')
        title('Acceleration')
        grid()
        legend(loc = 'best')
        xlim(0, time_array[-1]+0.1)
    
    return t_arr, v_arr, accel_array, h_arr, m_arr, round(h_arr[-1], 2), round(v_arr[c], 2)


