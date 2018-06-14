from numpy import *
from matplotlib.pyplot import *
from astropy.io import ascii
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d


def num_solver(thrust_profile, rocket_mass, motor_mass, propellant_mass, time_res, temp, burn_time, max_deploy, t_start, t_deploy, plots, drag_f):
    
    gauss_steepness= 0.005  #this measures the proportion that the gaussian starts off, so the bigger it is the bigger the jump but the faster ADAS deploys
    sigma_squared = ((-1)*(t_deploy)**2)/(2*log(gauss_steepness))
    t_start = t_start+burn_time
    ##t_deploy = t_deploy + burn_time - t_start
    
    
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

    time_step = time_res #rename
    
    #Get the mass flow rate from the thrust curve
    t_dummy = 0     #dummy time variable
    total_thrust = 0
    while t_dummy<burn_time+1:
        total_thrust = total_thrust+thrust_function(t_dummy)
        t_dummy = t_dummy+time_step
        
        
    drag_array = [
[0, 0.59182,2.2296,4.83375,8.55161,13.4116,20.0177,28.9221,38.1124,49.3355,62.0703,76.0614,92.0093,100.165],
[0, 0.591264,2.2295,4.83807,8.55265,13.4071,20.0157,28.9142,38.2055,49.2593,62.0156,76.2053,92.08,100.07],
[0, 0.592957,2.23896,4.85798,8.58222,13.4421,20.1047,28.965,38.327,49.34,62.133,76.276,92.111,100.161],
[0, 0.6,2.268,4.925,8.696,13.609,20.378,29.306,38.839,49.893,62.725,77.065,92.937,101.081],
[0, 0.611581,2.31152,5.02217,8.86343,13.8963,20.7639,29.8256,39.5542,50.7838,63.7972,78.3321,94.4829,102.622],
[0, 0.624584,2.36323,5.14082,9.07321,14.2113,21.2211,30.4146,40.2816,51.7134,65.021,79.7737,96.0817,104.417],
[0, 0.646593,2.45084,5.34285,9.43213,14.7778,22.0341,31.5141,41.7977,53.6209,67.1877,82.4599,99.4136,108.013],
[0, 0.661885,2.51603,5.49391,9.7072,15.2164,22.668,32.3353,42.7805,54.8745,68.7875,84.2702,101.464,110.237],
[0, 0.676409,2.57771,5.61966,9.95661,15.6125,23.1988,33.0944,43.7218,56.0576,70.2006,86.0549,103.6,112.28],
[0, 0.686241,2.61922,5.74393,10.1423,15.8945,23.5704,33.5766,44.4382,56.8404,71.3742,87.3685,105.032,114.097],
[0, 0.706858,2.7187,5.96159,10.5506,16.5219,24.4577,34.7788,46.0257,58.8808,73.8345,90.2789,108.543,117.836],
[0, 0.755909,2.91253,6.39603,11.3087,17.7266,26.2551,37.2794,49.436,63.0928,78.9964,96.7176,116.246,126.51],
[0, 0.796638,3.09445,6.83227,12.1446,19.0946,28.2083,39.7859,52.7009,67.1758,84.0781,102.727,123.357,134.042]]#DUNCAN GOT HIS SHIT TOGETHER :D
    
    ADAS_vel_array = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 250]
    ADAS_deploy_array = [0, 5.5555, 11.1111, 16.6667, 22.2222, 27.7778, 33.3333, 38.8888, 44.4444, 50, 55.5555, 77.7778, 100]
    
    drag_function = interp2d(ADAS_vel_array, ADAS_deploy_array, drag_array, kind='cubic')
    
    
    def mass_flow_rate(t):
        if t<=burn_time:
            return thrust_function(t)*(Wet_mass-Dry_mass)/total_thrust
        else:
            return 0
        
    def air_drag(hi, vi, f): #f is the angle
        return drag_function(vi, f) * air_pressure(hi, rho)/1.225
    
    def drag_curve(hi, vi, t):
        
        if t<=t_start:
            deployment = 0
        else:
            if t_start<t<t_deploy+t_start:
                deployment = (max_deploy)*e**(-1*(t-t_deploy-t_start)**2/(2*sigma_squared))
            else:
                deployment = max_deploy
        drag = air_drag(hi, vi, deployment)
        
        drag_curve_arr.append(drag)
        ADAS_deployment_arr.append(deployment)
        return drag
        
    def height_step(hi, vi):
        return  hi + vi * time_step 
    
    def velocity_step(vi, ai):
        return vi + ai*time_step
    
    def acc_step(hi, vi, mi, ti):
        return -g - drag_f * drag_curve(hi, vi, ti)/mi + thrust_function(ti)/mi
    
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
    drag_curve_arr = [0]
    ADAS_deployment_arr = [0]

    #calculation
    i = 1
    
    
    #ascent
   
    while (time_array[i-1]<burn_time): 
    
        ti = (time_step)*(i)
        time_array.append(ti)
    
        hi = height_step(h_arr[i-1], v_arr[i-1])
        h_arr.append(hi)
    
        vi = velocity_step(v_arr[i-1], a_arr[i-1])
        v_arr.append(vi)
        
        ai = acc_step(h_arr[i-1], v_arr[i-1], m_arr[i-1], time_array[i-1])
        a_arr.append(ai)
    
        mi = mass_step(m_arr[i-1], ti)
        m_arr.append(mi)
    
        i = i+1
    
    # coast
    c = i-1 #mark transition to coast

    while (v_arr[i-1]>0): 
    
        ti = (time_step)*(i)
        time_array.append(ti)
    
        hi = height_step(h_arr[i-1], v_arr[i-1])
        h_arr.append(hi)
    
        vi = velocity_step(v_arr[i-1], a_arr[i-1])
        v_arr.append(vi)
        
        ai = acc_step(h_arr[i-1], v_arr[i-1], m_arr[i-1], time_array[i-1])
        a_arr.append(ai)
    
        mi = Dry_mass
        m_arr.append(mi)
    
        i = i+1
       
    #descent

    d = i-1 #mark transition to descent
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
    
        figure(figsize=(10, 15))
        subplot(3,1,1)
        plot(t_arr[0:-1], drag_curve_arr[0:-1], '--', color = 'tomato', label = 'Trajectory')
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
        plot(t_arr[0:-1], a_arr[0:-1], '-', color = 'black', label = 'Acceleration')
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
        plot(t_arr[0:-1], ADAS_deployment_arr[0:-1], '-', color = 'black', label = 'ADAS Deployment')
        xlabel('Time [s]')
        ylabel('ADAS Dployment[deg]')
        title('ADAS Deployment')
        grid()
        legend(loc = 'best')
        xlim(0, time_array[-1]+0.1)
    
    return t_arr, v_arr, a_arr, h_arr, m_arr, ADAS_deployment_arr, round(h_arr[-1], 2), round(v_arr[c], 2)

#################################################################################################################################################################################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################################################################################################################################################################################

def PID(thrust_profile, rocket_mass, motor_mass, propellant_mass, time_res, temp, burn_time, max_deploy, t_start, t_deploy, plots, drag_f, k_p, k_d, k_i, PID_scale, ADAS_max_speed, update_interval, mixing_factor):
    
    nom_t, nom_v, nom_a, nom_h, nom_m, nom_ADAS, nom_apogee, nom_v_MECO = num_solver(thrust_profile, rocket_mass, motor_mass, propellant_mass, time_res, temp, burn_time, max_deploy, t_start, t_deploy, 0, drag_f)
    #calculate the nominal flight path
    
    
    
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

    time_step = time_res #rename
    
    #Get the mass flow rate from the thrust curve
    t_dummy = 0     #dummy time variable
    total_thrust = 0
    while t_dummy<burn_time+1:
        total_thrust = total_thrust+thrust_function(t_dummy)
        t_dummy = t_dummy+time_step
        
        
    drag_array = [
[0, 0.59182,2.2296,4.83375,8.55161,13.4116,20.0177,28.9221,38.1124,49.3355,62.0703,76.0614,92.0093,100.165],
[0, 0.591264,2.2295,4.83807,8.55265,13.4071,20.0157,28.9142,38.2055,49.2593,62.0156,76.2053,92.08,100.07],
[0, 0.592957,2.23896,4.85798,8.58222,13.4421,20.1047,28.965,38.327,49.34,62.133,76.276,92.111,100.161],
[0, 0.6,2.268,4.925,8.696,13.609,20.378,29.306,38.839,49.893,62.725,77.065,92.937,101.081],
[0, 0.611581,2.31152,5.02217,8.86343,13.8963,20.7639,29.8256,39.5542,50.7838,63.7972,78.3321,94.4829,102.622],
[0, 0.624584,2.36323,5.14082,9.07321,14.2113,21.2211,30.4146,40.2816,51.7134,65.021,79.7737,96.0817,104.417],
[0, 0.646593,2.45084,5.34285,9.43213,14.7778,22.0341,31.5141,41.7977,53.6209,67.1877,82.4599,99.4136,108.013],
[0, 0.661885,2.51603,5.49391,9.7072,15.2164,22.668,32.3353,42.7805,54.8745,68.7875,84.2702,101.464,110.237],
[0, 0.676409,2.57771,5.61966,9.95661,15.6125,23.1988,33.0944,43.7218,56.0576,70.2006,86.0549,103.6,112.28],
[0, 0.686241,2.61922,5.74393,10.1423,15.8945,23.5704,33.5766,44.4382,56.8404,71.3742,87.3685,105.032,114.097],
[0, 0.706858,2.7187,5.96159,10.5506,16.5219,24.4577,34.7788,46.0257,58.8808,73.8345,90.2789,108.543,117.836],
[0, 0.755909,2.91253,6.39603,11.3087,17.7266,26.2551,37.2794,49.436,63.0928,78.9964,96.7176,116.246,126.51],
[0, 0.796638,3.09445,6.83227,12.1446,19.0946,28.2083,39.7859,52.7009,67.1758,84.0781,102.727,123.357,134.042]]#DUNCAN GOT HIS SHIT TOGETHER :D
    
     ##############PID things#########
    
    start_heights = [298, 396, 489, 577, 660, 738, 812, 881, 946, 1008, 1066, 1120, 1172, 1220, 1265, 1306, 1345, 1381, 1415, 1437, 1454, 1470, 1486, 1501, 1514, 1527, 1539, 1549, 1559, 1568, 1576, 1584, 1590, 1595, 1600, 1604, 1606, 1608, 1609]

    function_index = 1 #keep track of which function we are currently calculating with
    prev_signal = 0  #the cur_signal from the previous loop
    deriv_signal = 0
    cur_signal = 0

    velocity_functions = [[ -0.000455874984234 ,  0.216170183163 ,  176.710523811 ],
                          [ 1.04868657812e-05 ,  -0.124391243726 ,  238.721409985 ],
                          [ 7.92643542567e-06 ,  -0.121867443237 ,  238.09878419 ],
                          [ -2.9669386665e-05 ,  -0.0770086057444 ,  224.7201655 ],
                          [ -5.85983840182e-07 ,  -0.117759718173 ,  238.93759641 ],
                          [ -9.63650725967e-07 ,  -0.117073698461 ,  238.637791425 ],
                          [ -2.45244048555e-06 ,  -0.114671767 ,  237.668761566 ],
                          [ -4.85400385168e-06 ,  -0.110416808308 ,  235.783785956 ],
                          [ -9.94098238419e-06 ,  -0.100765719385 ,  231.205588252 ],
                          [ -1.76104026014e-05 ,  -0.0852958397462 ,  223.403682762 ],
                          [ -2.62789842423e-05 ,  -0.0668150640905 ,  213.5528557 ],
                          [ -3.59707836192e-05 ,  -0.0450937677849 ,  201.381490807 ],
                          [ -4.71445236453e-05 ,  -0.0189039780727 ,  186.034211748 ],
                          [ -6.05064680359e-05 ,  0.0136992809133 ,  166.145195871 ],
                          [ -7.68529722073e-05 ,  0.0550534303274 ,  139.989260124 ],
                          [ -9.7330441225e-05 ,  0.108575152886 ,  105.015874164 ],
                          [ -0.000124134067604 ,  0.18072143409 ,  56.4662923145 ],
                          [ -0.000160151635125 ,  0.280264534494 ,  -12.3128839095 ],
                          [ -0.000201442981411 ,  0.396962285167 ,  -94.7667852488 ],
                          [ -0.000242138394989 ,  0.513843709934 ,  -178.691181158 ],
                          [ -0.000288541421794 ,  0.64879502734 ,  -276.809459232 ],
                          [ -0.000346447983581 ,  0.819110120969 ,  -402.042870368 ],
                          [ -0.000419664626094 ,  1.03672808983 ,  -563.747025569 ],
                          [ -0.00051360798885 ,  1.3186884767 ,  -775.315703504 ],
                          [ -0.000636158067973 ,  1.68984862901 ,  -1056.344088 ],
                          [ -0.000799086511593 ,  2.18743601704 ,  -1436.25585816 ],
                          [ -0.00102072532098 ,  2.86954097586 ,  -1961.06016126 ],
                          [ -0.00133010049267 ,  3.82836671438 ,  -2703.96735091 ],
                          [ -0.00177503542906 ,  5.21615510456 ,  -3786.12498453 ],
                          [ -0.00243830432254 ,  7.29690176513 ,  -5418.00872656 ],
                          [ -0.0034714283374 ,  10.5546767442 ,  -7986.21589242 ],
                          [ -0.0051716052346 ,  15.9403897623 ,  -12251.3488201 ],
                          [ -0.00817502186666 ,  25.492319445 ,  -19845.9828478 ],
                          [ -0.0140087260584 ,  44.1087432237 ,  -34698.0960943 ],
                          [ -0.026961899123 ,  85.5618903791 ,  -67863.015371 ],
                          [ -0.0622133407046 ,  198.630881534 ,  -158530.21714 ],
                          [ -0.198323557315 ,  635.951816446 ,  -509807.320707 ],
                          [ -1.32577931543 ,  4262.63229352 ,  -3426287.93913 ],
                          [ -53.2846829669 ,  171488.316639 ,  -137976998.018 ]]

    def calc_vel(height):
        function_index = 0
        while (height > start_heights[function_index]):
            function_index = function_index + 1
            if function_index == len(start_heights):
                break
        if (function_index == len(start_heights)):
            return 0; #retract fins
        order = len(velocity_functions[function_index])-1
  
        function_value = 0
        index = 0
        while (index <= order):
            function_value = function_value + velocity_functions[function_index][index]*height**(order-index)
            index = index + 1
        return function_value
    
    
    
    
    ADAS_vel_array = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 250]
    ADAS_deploy_array = [0, 5.5555, 11.1111, 16.6667, 22.2222, 27.7778, 33.3333, 38.8888, 44.4444, 50, 55.5555, 77.7778, 100]
    
    drag_function = interp2d(ADAS_vel_array, ADAS_deploy_array, drag_array, kind='cubic')
    
    prop_signal_arr = []
    deriv_signal_arr = []
    total_signal_arr = []
    
    max_index = len(nom_t)-1
    
    
    motor_commands_arr = [0] #stores the list of commands sent to the motor
    wanted_deployments = [0] #stores every calculation done by PID
    last_seized_index = [0]  #stores the index of where the data has already been drawn from from wanted_deployments
    
    def mass_flow_rate(t):
        if t<=burn_time:
            return thrust_function(t)*(Wet_mass-Dry_mass)/total_thrust
        else:
            return 0
        
    def air_drag(hi, vi, f): #f is the angle
        return drag_function(vi, f) * air_pressure(hi, rho)/1.225
    
    def update_PID(hi, vi, ai, t):
        prev_deployment = ADAS_deployment_arr[-1]
        if len(prop_signal_arr) != 0:
            prev_signal = prop_signal_arr[-1]
        if t<=burn_time:
            cur_signal = 0
            deriv_signal = 0
        else:
            my_height = hi
            wanted_vel = calc_vel(hi)
            cur_signal = (vi-wanted_vel)
            deriv_signal = (cur_signal-prev_signal)/(time_step)
        
        total_signal = PID_scale*(k_p*cur_signal + k_d*deriv_signal)
        prop_signal_arr.append(cur_signal)
        deriv_signal_arr.append(deriv_signal)
        total_signal_arr.append(total_signal)
        
        deployment = total_signal*mixing_factor + prev_deployment
        
        
        #check that it isn't trying to deploy too fast
        if abs(deployment-prev_deployment) > ADAS_max_speed*update_interval:
            sign = abs(deployment-prev_deployment)/(deployment-prev_deployment)
            deployment = prev_deployment + sign*ADAS_max_speed*update_interval
        
        #check that the signal isn't above 100 or below 0
        if deployment > 100:
            deployment = 100
        if deployment < 0:
            deployment = 0
            
        wanted_deployments.append(deployment)
        
    
    def move_ADAS(t):
        '''
        if (t/time_step)%(update_interval/time_step) == 0:    #if it is time to update the motor configuration
            av_depl = average(wanted_deployments[last_seized_index[-1]:-1])
            last_seized_index.append(int(t/time_step))              #update for the next iteration
            motor_commands_arr.append(av_depl)
            
        specific_wanted_deployment = motor_commands_arr[-1]
            
        current_ADAS_config = ADAS_deployment_arr[-1]
        if specific_wanted_deployment == current_ADAS_config:
            ADAS_deployment_arr.append(current_ADAS_config)
            return
        
        sign = abs(specific_wanted_deployment-current_ADAS_config)/(specific_wanted_deployment-current_ADAS_config) #whether to deploy more or less
        
        new_deployment = current_ADAS_config + sign*ADAS_max_speed*time_step
        '''
        new_deployment = wanted_deployments[-1]
        ADAS_deployment_arr.append(new_deployment)
    
    def drag_curve(hi, vi, t):
        drag = air_drag(hi, vi, ADAS_deployment_arr[-1])
        drag_curve_arr.append(drag)
        
        return drag
        
    def height_step(hi, vi):
        return  hi + vi * time_step 
    
    def velocity_step(vi, ai):
        return vi + ai*time_step
    
    def acc_step(hi, vi, mi, ti):
        return -g - drag_f * drag_curve(hi, vi, ti)/mi + thrust_function(ti)/mi
    
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
    drag_curve_arr = [0]
    ADAS_deployment_arr = [0]

    #calculation
    i = 1
    
    
    #ascent
   
    while (time_array[i-1] < burn_time): 
    
        ti = (time_step)*(i)
        time_array.append(ti)
    
        hi = height_step(h_arr[i-1], v_arr[i-1])
        h_arr.append(hi)
    
        vi = velocity_step(v_arr[i-1], a_arr[i-1])
        v_arr.append(vi)
        
        ai = acc_step(h_arr[i-1], v_arr[i-1], m_arr[i-1], time_array[i-1])
        a_arr.append(ai)
    
        mi = mass_step(m_arr[i-1], ti)
        m_arr.append(mi)
        
        #update ADAS calculation
        update_PID(h_arr[i-1], v_arr[i-1], a_arr[i-1], time_array[i-1])
        #move motor
        move_ADAS(ti)
    
        i = i+1
    
    # coast
    c = i-1 #mark transition to coast

    while (v_arr[i-1]>0): 
    
        ti = (time_step)*(i)
        time_array.append(ti)
    
        hi = height_step(h_arr[i-1], v_arr[i-1])
        h_arr.append(hi)
    
        vi = velocity_step(v_arr[i-1], a_arr[i-1])
        v_arr.append(vi)
        
        ai = acc_step(h_arr[i-1], v_arr[i-1], m_arr[i-1], time_array[i-1])
        a_arr.append(ai)
    
        mi = Dry_mass
        m_arr.append(mi)
        
        #update ADAS calculation
        update_PID(h_arr[i-1], v_arr[i-1], a_arr[i-1], time_array[i-1])
        #move motor
        move_ADAS(ti)
    
        i = i+1
       
    #descent

    d = i-1 #mark transition to descent
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
        '''
        figure(figsize=(10, 15))
        subplot(3,1,1)
        plot(t_arr[0:-1], drag_curve_arr, '--', color = 'tomato', label = 'Trajectory')
        grid()
        xlim(0, time_array[-1]+0.1)
        xlabel('Time [sec]')
        ylabel('Drag [N]')
        '''
        
        figure(figsize=(10, 15))
        subplot(3,1,1)
        plot(t_arr[0:-1], ADAS_deployment_arr[0:-1], '-', color = 'black', label = 'ADAS Deployment')
        plot(nom_t[0:-1], nom_ADAS[0:-1], '--', color = 'tomato', label = 'Nominal ADAS Deployment')
        xlabel('Time [s]')
        ylabel('ADAS Dployment[%]')
        title('ADAS Deployment')
        grid()
        legend(loc = 'best')
        xlim(0, time_array[-1]+0.1)
        
        figure(figsize=(10, 15))
        subplot(3,1,1)
        plot(t_arr[0:-2], prop_signal_arr[0:-1], '-', color = 'black', label = 'Proportional Signal')
        plot(t_arr[0:-2], deriv_signal_arr[0:-1], '-', color = 'tomato', label = 'Derivative Signal')
        plot(t_arr[0:-2], total_signal_arr[0:-1], '-', color = 'blue', label = 'Total Signal')
        xlabel('Time [s]')
        ylabel('Signal')
        title('PID Signals')
        grid()
        legend(loc = 'best')
        xlim(0, time_array[-1]+0.1)
        ylim(-10,10)
        
        figure(figsize=(10, 15))
        subplot(3,1,1)
        plot(t_arr,array(h_arr), '-', color = 'black', label = 'Trajectory')
        plot(nom_t,array(nom_h), '--', color = 'tomato', label = 'Nominal Trajectory')
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
        plot(t_arr, array(v_arr), '-', color = 'black', label = 'Velocity')
        plot(nom_t, nom_v, '--', color = 'tomato', label = 'Nominal Velocity')
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
        plot(t_arr[0:-1], a_arr[0:-1], '-', color = 'black', label = 'Acceleration')
        plot(nom_t, nom_a, '--', color = 'tomato', label = 'Nominal Acceleration')
        xlabel('Time [s]')
        ylabel('Acceleration [m/sec^2]')
        title('Acceleration')
        grid()
        legend(loc = 'best')
        xlim(0, time_array[-1]+0.1)
        
        
    
    return t_arr, v_arr, a_arr, h_arr, m_arr, ADAS_deployment_arr, round(h_arr[-1], 2), round(v_arr[c], 2)