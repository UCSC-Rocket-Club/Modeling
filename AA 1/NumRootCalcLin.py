from numpy import *

alpha = 0.1025  #ie 1.225/(2*5)
g = 9.8
height_diff = 206.0  #x_t-x_0
v_0 = 65.0


'''
#for debugging, saves the iterative guesses
x_values = []
y_values = []
'''

def updateParam(newMass, newAirDensity, newHeightDiff, newVel0):
    #newMass [kg]
    #newAirDensity [kg m^-3]
    #newHeightDiff [m]
    #newVel0 [m/s]
    global alpha
    alpha = newAirDensity/(2*newMass)
    global height_diff
    height_diff = newHeightDiff
    global v_0
    v_0 = newVel0
    
#returns the value of f(DELTA), the function as in the document (equation _)
def function(x):#returns the value of the function
    return exp(x*alpha*height_diff)-sqrt(1+alpha*v_0**2*x/g)
    
#returns the value of f'(DELTA) as in the document (equation _)
def grad(x):
    return alpha*height_diff*exp(x*alpha*height_diff)-(alpha*v_0**2)/(2*g*sqrt(1+alpha*v_0**2*x/g))

def calculateRoots(guess_x, threshold, max_iterations):
    
    y = function(guess_x)
    counter = 0
    while abs(y)>threshold and counter<max_iterations:
        #calculate the new x guess
        y = float(function(guess_x))
        m = float(grad(guess_x))
        
        #there is always a root at x=0 so we want to avoid convergence towards that
        if m < 0:
            guess_x = guess_x*2
        else:
            guess_x = guess_x-y/m
        counter = counter+1
        
        '''
        #for debugging, saves the iterative guesses
        x_values.append(guess_x)
        y_values.append(y)
        '''
        
    if counter == max_iterations:
        print('Root could not be converged on within the iterations specified')
    return guess_x