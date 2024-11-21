using SNOW
using OrdinaryDiffEq
using Plots

#global variables
const barrel_radius = 1 #cm
const wire_radius = 0.06438 / 2 #cm
const unit_resistance = 16.14 * (1/1000) * (1/12) * (1/2.54) #ohms/cm 
const EMF_constant = 0.1
const mass = 0.01 #kg
const friction_coefficient = 0.03 
const max_voltage = 12 #volts
const diode_resistance = 5 #ohms
const forward_voltage = 0.51 #volts

"""
coil(min_radius, coil_length, wire_radius, wire_length, unit_resistance)
This function takes in some of the constraints for building a coil and
    estimates the resistance (in ohms) and inductance (in mH)

min_radius is the radius of the object the wire is being wrapped around (in cm)
coil_length is the length of the coiled wire (as opposed to the length of the wire) (in cm)
wire_radius is the radius of the wire (in cm)
wire_length is the length of the wire before it is wound into a coil (in cm)
unit_resistance is the resistance per unit length of the wire (ohms/cm)

this function returns a tuple with resistance(ohms) and inductance (Henrys)

This function is based off of "Multilayer air-core coil" on the wikipedia page on inductors
This function treats the the coils like helixes where each progressive layer sits in the gap 
    made by two coils on a lower layer.
"""
function coil(min_radius, coil_length, wire_radius, wire_length, unit_resistance)
    #calculate the resistance
    resistance = unit_resistance * wire_length
    #initialize some variables
    turns = 0 #this is N on wikipedia
    layers = 0 
    r_add = 0.0 #this is to help me calculate the average radius
    wire_left = wire_length #each layer will use up some of the wire
    layer_turns = round(coil_length / (2*wire_radius)) # we can make this many turns on the first layer
    layer_radius = min_radius + wire_radius #The radius from the center of the core to the center of the wire
    layer_length = layer_turns * sqrt((layer_radius^2) + (wire_radius^2)) #how much wire does it take to make the first layer
    while wire_left >= layer_length && layer_turns > 0 #while we still have enough wire to wind the next layer, do it
        layers = layers + 1 #I don't think I need this variable
        r_add = r_add + (layer_radius * layer_turns) #the summation part of an average radius of all the turns
        turns = turns + layer_turns #add up how many turns we have so far
        wire_left = wire_left - layer_length #how much wire do we have left after making this layer?
        layer_radius = layer_radius + (wire_radius*sqrt(3)) #what is the next layer's radius (think eqilateral triangle)
        layer_turns = layer_turns - 1 #the next layer will have one fewer turns because those turns will be nested inbetween this layer's
        layer_length = layer_turns * sqrt((layer_radius^2) + (layer_length^2)) #how much wire will the next layer take?
    end
    average_radius = r_add / turns #finish calculating the average radius of all the turns
    depth = (layer_radius - (wire_radius*sqrt(3))) - (min_radius + wire_radius)#go back to the radius of the last layer, and get the depth
    inductance = ((average_radius^2)*(turns^2)) / ((19*average_radius)+(29*coil_length)+(32*depth)) #calculate the inductance
    inductance = inductance / (10^6) #convert uH to H
    return resistance, inductance
end

"""
single_stage(du, u, p, t)
This function has the EOM for a simple, single-stage coil gun. It is used in the ODE solver

du is a vector with the derivatives of the state variables
u is a vector of the state variables
p is a vector of parameters:
    p = [inductance, resistance, EMF_constant, mass, friction_coefficient, voltage]
t is the time
"""
function single_stage(du, u, p, t)
    #pull the variables out of the vector p, and give them more legible names
    #pull out the values from p
    L_a = p[1]
    R = p[2]
    K_e = p[3]#assuming K_e and K_f are the same
    m = p[4]
    mu_k = p[5]
    g = 9.81
    V_a = p[6]

    #this is the juicy part... the EOM (du is the derivative and u is the state vector)
    du[1] = -(R/L_a)*u[1] - (K_e/L_a)*u[3] + (V_a/L_a)
    du[2] = u[3]
    du[3] = (K_e/m)*u[1] + mu_k*g
end

"""
get_solution(x)
This function takes the design variables and runs the solver on them. It also 
    contains a lot of conditions for when things within the equations change. 
    Some examples are when the mass 

x is a vector of the design variables:
    x[1] = cutoff time
    x[2] = length of the coil
    x[3] = length of wire
    x[4] = end of barrel (this one is constrained to only be the one value)

This function returns a solution object which carries the information from the
     differential equation solver
"""
function get_solution(x)
    #set up variables
    cutoff_time = x[1] # s
    coil_length = x[2] #cm
    wire_length = x[3] #cm
    
    #calculate the resistance and inductance of the coil
    resistance, inductance = coil(barrel_radius, coil_length, wire_radius, wire_length, unit_resistance)
    
    #initialize the parameters (see the values at the top of this file)
    p = [inductance, resistance, EMF_constant, mass, friction_coefficient, max_voltage] #put the variables into the vector
    #time
    t = (0.0, 1.0) #The solver will evaluate the differential equations for 1 second unless something stops it sooner
    #set the initial conditions
    u0 = [0.0, -0.05, 0.0] #[current (amps), position(meters), velocity(m/s)]
    
    #objective (velocity when x = pos)
    pos = x[4] #pos is representative of the end of the barrel, and we want to stop simulating when the mass crosses that point
        
    #set up end condition for simulation (end the simulation when the mass reaches pos)
    condition(u,t,integrator) = u[2] - pos #the condition is triggered when u[2] - pos == 0, AKA: when the mass reaches pos
    affect!(integrator) = terminate!(integrator) #the affect on the solver is that it ends the simulation
    cb = ContinuousCallback(condition, affect!) #pack the condition and affect into a single variable
    
    #set up voltage changes
    #voltage change when t = cutoff_time
    tstop = [cutoff_time]
    function condition2(u, t, integrator) #condition2 is satisfied when t == cutoff_time
        t in tstop
    end
    function affect2!(integrator) #affect2 changes the voltage to the negative of the forward_voltage of the diode
        integrator.p[6] = -forward_voltage 
    end
    cb2 = DiscreteCallback(condition2, affect2!) #might need to add something to save the position <---I don't think we do
    #voltage change when diode_resistance*current < forward_voltage
    function condition3(u,t,integrator) #condition3 is satisfied when the voltage has 1) already been set to the negative of the forward_voltage, 2) the diode_resistance * current == forward_voltage of the diode
        if integrator.p[6] == -forward_voltage
            return diode_resistance*u[1] - forward_voltage
        else
            return 1
        end
    end
    function affect3!(integrator)#affect3 is that the voltage becomes the diode_resistance * the current at the time condition3 is satisfied. I couldn't find out how to make it a function of current, but it shouldn't have too large an affect
        integrator.p[6] = diode_resistance * u[3]#notice that this puts the voltage at a fixed voltage
    end
    cb3 = ContinuousCallback(condition3, affect3!)
    
    #change sign of k_e whenever x == 0
    condition4(u,t,integrator) = u[2] #condition4 is satisfied whenever the mass is at position 0
    function affect4!(integrator) #affect4 is to change the sign of k_e
        integrator.p[3] = -1 * integrator.p[3]
    end
    cb4 = ContinuousCallback(condition, affect!)

    #change sign of mu_k whenever xdot == 0
    condition5(u,t,integrator) = u[3] # condition5 is satisfied whenever the velocity of the mass is 0
    function affect5!(integrator) #affect5 is to change the sign of the friction_coefficient
        integrator.p[5] = -1 * integrator.p[5]
    end
    cb5 = ContinuousCallback(condition, affect!)
        
    ##combine the callbacks
    cbs = CallbackSet(cb, cb2, cb3, cb4, cb5)
    
    #simulate
    prob = ODEProblem(single_stage, u0, t, p) #set up the problem
    sol = solve(prob, callback = cbs); #solve the problem

    return sol
end

"""
function ob!(g, x)
This function is used in the optimizer. 
The vector g is a vector of constraint values (which we don't really use in this problem) and x is a vector of 
    variables the optimizer can vary to find the optimum

    x[1] = cutoff time
    x[2] = length of the coil
    x[3] = length of wire
    x[4] = end of barrel (this one is constrained to only be the one value)

"""
#objective function to optimize coil gun
#let's vary: cutoff_time
function ob!(g, x)
    #run the simulation
    sol = get_solution(x)

    # get the velocity when the position is pos
    f = -sol[length(sol)][3] #this is negative because we are minimizing in the optimizer

    #constraints
    g[1] = x[1] #the voltage can't be turned off at a negative value nor one larger than 1s
    
    #return
    return f
end

"""
optimize_coil_gun()
This function runs, without any inputs, the optimizer for the single stage coil gun.
    x[1] = cutoff time
    x[2] = length of the coil
    x[3] = length of wire
    x[4] = end of barrel (this one is constrained to only be the one value)

To read the results, look for "xstar = _____" This is the solution the optimizer found.
    "fstar = _____" is the value of the objective. (The velocity of the projectile at the end of the barrel)
"""
function optimize_coil_gun()
    max_wire_length = 200 #cm
    barrel_length = 0.1 #m -- we can change this variable to match the length of the actual barrel
    shortest_coil = 0.08 #cm -- here I am trying to make it so that we have at least one turn
    longest_coil = (2*0.06438*max_wire_length) / (sqrt((1.5^2)+(0.06438^2))) #the case where we use all of the wire on a single layer
    x0 = [0.0, 0.5, 100, barrel_length]  # starting point
    lx = [0.0, shortest_coil, 5, barrel_length]  # lower bounds on x
    ux = [1.0, longest_coil, 200, barrel_length]  # upper bounds on x
    ng = 1  # number of constraints
    lg = [0.0]  # lower bounds on g
    ug = [1.0]  # upper bounds on g
    options = Options(solver=IPOPT())  # choosing IPOPT solver

    xopt, fopt, info = minimize(ob!, x0, ng, lx, ux, lg, ug, options)

    println("xstar = ", xopt)
    println("fstar = ", fopt)
    println("info = ", info)

    return xopt
end

"""
plot_a_sol(solution)
This function lets you plot the results of a simulation

solution is a vector of design variables:
    solution[1] = cutoff time
    solution[2] = length of the coil
    solution[3] = length of wire
    solution[4] = end of barrel (this one is constrained to only be the one value)
"""
function plot_a_sol(solution)
    #calculate the resistance and inductance of the coil
    resistance, inductance = coil(barrel_radius, solution[2], wire_radius, solution[3], unit_resistance)

    println("Resistance: " * string(resistance) * " Ohms")
    println("Inductance: " * string(inductance) * " Henrys")

    sol = get_solution(solution)

    plot(sol, label = ["current (A)" "position(m)" "velocity(m/s)"])
    plot!(display = true)
end

"""
plot_vary_cutoff(coil_length, wire_length, barrel_length)
This function takes the variables that define a coil (and a barrel length) and plots
    the end velocities of the mass for a range of different cutoff times.
"""
function plot_vary_cutoff(coil_length, wire_length, barrel_length)
    #make a vector of vectors with different cutoff times
    inputs = []
    times = range(start = 0.001, stop = 0.1, step = 0.0001)
    for t in times
        push!(inputs, [t, coil_length, wire_length, barrel_length])
    end
    #run a simulation for each cutofftime and store the end velocities
    velocities = []
    for i in inputs
        sol = get_solution(i)
        push!(velocities, sol[length(sol)][3])
    end
    plot(times, velocities, xlabel = "time(s)", ylabel = "velocity(m/s)")
    plot!(display = true)
end

"""
plot_vary_wire_length(cutoff_time, coil_length, barrel_length)
This function takes the cutoff time and coil length (and a barrel length) and plots
    the end velocities of the mass for a range of different wire lengths.
"""
function plot_vary_wire_length(cutoff_time, coil_length, barrel_length)
    #make a vector of vectors with different cutoff times
    inputs = []
    lengths = range(start = 10, stop = 200, step = 1)
    for l in lengths
        push!(inputs, [cutoff_time, coil_length, l, barrel_length])
    end
    #run a simulation for each cutofftime and store the end velocities
    velocities = []
    for i in inputs
        sol = get_solution(i)
        push!(velocities, sol[length(sol)][3])
    end
    plot(lengths, velocities, xlabel = "wire length(cm)", ylabel = "velocity(m/s)")
    plot!(display = true)
end

"""
plot_vary_coil_length(cutoff_time, coil_length, barrel_length)
This function takes the cutoff time and coil length (and a barrel length) and plots
    the end velocities of the mass for a range of different wire lengths.
"""
function plot_vary_coil_length(cutoff_time, wire_length, barrel_length)
    #make a vector of vectors with different cutoff times
    inputs = []
    barrel_length = 0.1 #m -- we can change this variable to match the length of the actual barrel
    shortest_coil = 0.08 #cm -- here I am trying to make it so that we have at least one turn
    longest_coil = (2*0.06438*wire_length) / (sqrt((1.5^2)+(0.06438^2))) #the case where we use all of the wire on a single layer
    lengths = range(start = shortest_coil, stop = longest_coil, step = 0.1)
    for l in lengths
        push!(inputs, [cutoff_time, l, wire_length, barrel_length])
    end
    #run a simulation for each cutofftime and store the end velocities
    velocities = []
    for i in inputs
        sol = get_solution(i)
        push!(velocities, sol[length(sol)][3])
    end
    plot(lengths, velocities, xlabel = "coil length(cm)", ylabel = "velocity(m/s)")
    plot!(display = true)
end

#this is the code that actually runs
solution = optimize_coil_gun() #This finds the optimal solution

#For some reason on my machine, I can only create one plot per run. so I comment out everything but the one I want
plot_a_sol(solution) #this plots what the state space for that solution looks like
#plot_vary_cutoff(solution[2], solution[3], solution[4]) #this plots a what-if we varyied the cutoff time
#plot_vary_coil_length(solution[1], solution[3], solution[4]) #this plots a what-if we varyied the coil length
#plot_vary_wire_length(solution[1], solution[2], solution[4])#this plots a what-if we varyied the wire length