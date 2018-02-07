from numpy import *
from matplotlib.pyplot import *
from astropy.io import ascii
from scipy.interpolate import interp1d

def num_solver(thrust_profile, rocket_mass, motor_mass, propellant_mass, time_res, temp, burn_time, max_deploy, t_start, t_deploy, plots, drag_f):
    
    gauss_steepness= 0.005  #this measures the proportion that the gaussian starts off, so the bigger it is the bigger the jump but the faster ADAS deploys
    sigma_squared = ((-1)*(t_deploy)**2)/(2*log(gauss_steepness))
    t_start = t_start+burn_time
    ##t_deploy = t_deploy + burn_time - t_start
    
    
    rocket_time = ascii.read(thrust_profile)['time']
    thrust_curve = ascii.read(thrust_profile)['thrust']

    thrust_function = interp1d(rocket_time, thrust_curve)

    drag_curve_arr = []
    ADAS_deployment_arr = []
    
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
    

    time_step = time_res #sec

    #Get the mass flow rate from the thrust curve
    t_dummy = 0     #dummy time variable
    total_thrust = 0
    while t_dummy<burn_time+1:
        total_thrust = total_thrust+thrust_function(t_dummy)
        t_dummy = t_dummy+time_step
        
    ADAS_drag_times_area = [0.002397941193, 0.002398131584, 0.002403933158, 0.002431870269, 0.002475983869, 0.002525946346, 0.002618169871, 0.002683394918, 0.00274357945, 0.002788191485, 0.002887541676, 0.002993059848, 0.003095213036, 0.003197404042, 0.003298209003]
    
    ADAS_deployments = [0.000000, 5.555556, 11.111111, 16.666667, 22.222222, 27.777778, 33.333333, 38.888889, 44.444444, 50.000000, 55.555556, 66.666667, 77.777778, 88.888889, 100.000000]
    
    
    drag_function = interp1d(ADAS_deployments, ADAS_drag_times_area, kind='quadratic')
    
    
    def mass_flow_rate(t):
        if t<=burn_time:
            return thrust_function(t)*(Wet_mass-Dry_mass)/total_thrust
        else:
            return 0
    
    def drag_curve(hi, vi, t): #
        '''
        if t<=t_start:
            deployment = 0
        else:
            if t_start<t<t_deploy+t_start:
                deployment = (max_deploy)*e**(-1*(t-t_deploy-t_start)**2/(2*sigma_squared))
            else:
                deployment = max_deploy
        '''
                
        #stop at 20, 40, 60, 80 and 100 percent
        deployment = 0
        
        if t<=t_start:
            deployment = 0
        for j in range (0,10):
            if t_start+t_deploy*(j-1)<t<=t_start+t_deploy*(j):
                deployment = 25*j
        
        drag = drag_function(deployment)
        drag_curve_arr.append(0.5*air_pressure(hi)*drag*vi*vi)
        ADAS_deployment_arr.append(deployment)
        return drag
        
    def height_step(hi, vi):
        return  hi + vi * time_step 
    
    def velocity_step(hi, vi, mi, ti):
        return vi +(-g -drag_f * 0.5*air_pressure(hi) *vi*vi* drag_curve(hi, vi, ti)/mi + thrust_function(ti)/mi)*time_step
    
    def mass_step(mi, t): 
        return mi + (-mass_flow_rate(t))
    
    def air_pressure(h):
        R=8.31447       #Universal gas constant
        M=0.0289644     #Molar mass of air kg/mol
        T=temp+273.15   #temperature of the air at launch altitude in Kelvin
        return rho*e**(-1*(g*M*h)/(R*T))

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
    """
    while (h_arr[i]>0):
    
        ti = (time_step)*(i)
        time_array.append(ti)
        
        hi = height_step(h_arr[i], v_arr[i])
        h_arr.append(hi)
    
        vi = velocity_step(h_arr[i], v_arr[i], m_arr[i], time_array[i])
        v_arr.append(vi)
    
        mi = Dry_mass
        m_arr.append(mi)
    
        i = i+1
    """
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
        plot(t_arr[0:-1], drag_curve_arr, '--', color = 'tomato', label = 'Trajectory')
        grid()
        xlim(0, time_array[-1]+0.1)
        xlabel('Time [sec]')
        ylabel('Drag [N]')
   
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
        
        figure(figsize=(10, 15))
        subplot(3,1,1)
        plot(t_arr[0:-2], ADAS_deployment_arr[0:-1], '-', color = 'black', label = 'ADAS Deployment')
        xlabel('Time [s]')
        ylabel('ADAS Dployment[deg]')
        title('ADAS Deployment')
        grid()
        legend(loc = 'best')
        xlim(0, time_array[-1]+0.1)
    
    return t_arr, v_arr, accel_array, h_arr, m_arr, ADAS_deployment_arr, round(h_arr[-1], 2), round(v_arr[c], 2)


