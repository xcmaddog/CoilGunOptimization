using SNOW
using OrdinaryDiffEq
using Plots
using Colors

#global variables
const EMF_constant = 0.1
const mass = 0.001472 #kg
const friction_coefficient = 0.28 
const max_voltage = 12 #volts
const diode_resistance = 5 #ohms
const forward_voltage = 0.51 #volts
const selected_L = 51e-6 #Henrys
const selected_R = 0.2 #ohms
const length_of_wire = 800 #cm
const length_of_coil = 4.5 #cm
const barrel_length = 0.1 #m

"""
single_stage(du, u, p, t)
This function has the EOM for a simple, single-stage coil gun. It is used in the ODE solver

du is a vector with the derivatives of the state variables
u is a vector of the state variables
p is a vector of parameters:
    p = [inductance, resistance, EMF_constant, mass, friction_coefficient, voltage, against_diode]
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

    #this bit allows the callbacks to change the voltage between voltage sources vs diode resistance
    
    if p[7] == 0.0 
        V_a = p[6]
    else 
        V_a = -diode_resistance*u[1]
    end

    #this is the juicy part... the EOM (du is the derivative and u is the state vector)
    du[1] = -(R/L_a)*u[1] - (K_e/L_a)*u[3] + (V_a/L_a)
    du[2] = u[3]
    du[3] = (K_e/m)*u[1] + mu_k*g
end

"""
get_solution(x)
This function takes the design variables and runs the solver on them. It also 
    contains a conditions for when things within the equations change. 

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

    #grab the global variables
    inductance = selected_L
    resistance = selected_R
    
    #initialize the parameters (see the values at the top of this file)
    p = [inductance, resistance, EMF_constant, mass, friction_coefficient, max_voltage, 0.0] #put the variables into the vector
    #time
    t = (0.0, 0.2) #The solver will evaluate the differential equations for 0.2 seconds unless something stops it sooner
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
    condition2(u, t, integrator) = abs(t - tstop[1]) < 1e-3 #condition2 is satisfied when t == cutoff_time
    function affect2!(integrator) #affect2 changes the voltage to the negative of the forward_voltage of the diode
        integrator.p[6] = -forward_voltage 
    end
    cb2 = DiscreteCallback(condition2, affect2!) #might need to add something to save the position <---I don't think we do
    
    #voltage change when diode_resistance*current < forward_voltage
    function condition3(u,t,integrator) #condition3 is satisfied when the voltage has 1) already been set to the negative of the forward_voltage, 2) the diode_resistance * current == forward_voltage of the diode
        if integrator.p[6] == -forward_voltage
            return diode_resistance*u[1] - forward_voltage
        else
            return 1 #hopefully this doesn't mess with the finder too hard, but will still keep it from hitting 0 before the cutoff time
        end
    end
    function affect3!(integrator)#affect3 is that the voltage becomes the diode_resistance * the current at the time condition3 is satisfied. I couldn't find out how to make it a function of current, but it shouldn't have too large an affect
        integrator.p[7] = 1.0 #change this so that the ODE will swap to using a voltage that varies with current
    end
    cb3 = ContinuousCallback(condition3, affect3!)
    
    #change sign of k_e whenever x == 0
    condition4(u,t,integrator) = u[2] #condition4 is satisfied whenever the mass is at position 0
    function affect4!(integrator) #affect4 is to change the sign of k_e
        integrator.p[3] = -1 * integrator.p[3]
    end
    cb4 = ContinuousCallback(condition4, affect4!)

    #change sign of mu_k whenever xdot == 0
    condition5(u,t,integrator) = u[3] # condition5 is satisfied whenever the velocity of the mass is 0
    function affect5!(integrator) #affect5 is to change the sign of the friction_coefficient
        integrator.p[5] = -1 * integrator.p[5]
    end
    cb5 = ContinuousCallback(condition5, affect5!)

    #this would be good code to implement in a real situation, but it works so well that it doesn't give us a 
        #chance to explain the physics of the system.
    #set up end condition for simulation (end the simulation when the mass exits going backwards)
    #function condition6(u,t,integrator)#the condition is triggered when u[2] == 0 and u[3] is negative
    #    if u[3] < 0
    #        return u[2]
    #    else
    #        return 1
    #    end
    #end
    #affect6!(integrator) = terminate!(integrator) #the affect on the solver is that it ends the simulation
    #cb6 = ContinuousCallback(condition6, affect6!) #pack the condition and affect into a single variable
        
    ##combine the callbacks
    cbs = CallbackSet(cb, cb2, cb3, cb4, cb5)
    
    #simulate
    prob = ODEProblem(single_stage, u0, t, p) #set up the problem
    sol = solve(prob, callback = cbs, adaptive = false, dt = 1e-6); #solve the problem
    #sol = solve(prob, callback = cbs)

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
    g[1] = 0.5 #this is a filler constraint. It is always satisfied

    #return
    return f
end

"""
function ob_big!(g,x)
This function lets me quickly simplify the much more complicated ob! function. 
    It depends on the global variables to fill in the rest of the values that ob! requires.

It would be nice to get rid of this function and be able to optimize for the length of the coil and the length of the wire too.
"""
function ob_big!(g,x)
    xPrime = [x[1] length_of_coil length_of_wire barrel_length]
    return ob!(g, xPrime)
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

    x0 = [0.01] #starting cutoff time
    lx = [0.0] #min cutoff time
    ux = [0.2] #max cutoff time

    ng = 1  # number of constraints
    lg = [0.0]  # lower bounds on g
    ug = [1.0]  # upper bounds on g

    options = Options(solver=IPOPT())  # choosing IPOPT solver

    xopt, fopt, info = minimize(ob_big!, x0, ng, lx, ux, lg, ug, options)

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
"""
function plot_a_sol(solution)
    cut_time = solution[1]

    sol = get_solution([cut_time, length_of_coil, length_of_wire, barrel_length])

    #get time range
    tspan = range(start = 0, stop = sol.t[length(sol.t)], step = 0.0001)
    data = zeros(length(tspan), 3)
    i = 1
    for t in tspan
        s = sol(t)
        data[i, 1] = s[1]
        data[i, 2] = s[2]
        data[i, 3] = s[3]
        i = i + 1
    end
    # First Plot: Left Y-Axis
    p = plot(tspan, data[:, 1],
        label = "", 
        ylabel = "Current (A)",
        xlabel = "Time (s)",
        linewidth=2,
        color = RGB(0.337, 0.631, 0.749),
        ytickfontcolor = RGB(0.337, 0.631, 0.749))

    # Second Plot: First Right Y-Axis
    p2 = twinx(p)  # Attach the first twin axis to `p`
        plot!(p2, tspan, data[:, 2], 
        label = "",
        ylabel = "\n"*"Position (m)",
        linewidth=2,
        color = RGB(0.2, 0.468, 0.2),
        foreground_color_guide = RGB(0.2, 0.468, 0.2),
        ytickfontcolor = RGB(0.2, 0.468, 0.2))

    # Third Plot: Second Right Y-Axis
    p3 = twinx(p)  # Attach another twin axis to `p`
        plot!(p3, tspan, data[:, 3], 
        label = "",
        ylabel = "Velocity (m/s)",
        linewidth=2,
        color = RGB(0.949, 0.352, 0.219),
        foreground_color_guide = RGB(0.949, 0.352, 0.219),
        ytickfontcolor = RGB(0.949, 0.352, 0.219))

    #save plot
    #savefig("C:\\Users\\xcmad\\Documents\\Homework\\ME EN 335\\CoilGunProject\\plot_of_a_shot")
end

"""
plot_vary_cutoff(;coil_length = length_of_coil, wire_length = length_of_wire, barrel_length = barrel_length)
This function optionally takes the variables that define a coil (and a barrel length) and plots
    the end velocities of the mass for a range of different cutoff times.
If variables aren't provided, it uses the global variables defined at the top of this file.
"""
function plot_vary_cutoff(;coil_length = length_of_coil, wire_length = length_of_wire, barrel_length = barrel_length)
    #make a vector of vectors with different cutoff times
    inputs = []
    times = range(start = 0.001, stop = 0.1, step = 0.001)
    for t in times
        push!(inputs, [t, coil_length, wire_length, barrel_length])
    end
    #run a simulation for each cutofftime and store the end velocities
    velocities = []
    #currents = []
    positions = []
    for i in inputs
        sol = get_solution(i)
        push!(velocities, sol[length(sol)][3])
        #push!(currents, sol[length(sol)][1])
        push!(positions, sol[length(sol)][2])

        #println(i)
    end
    p = plot(times, velocities,
            xlabel = "Cutoff Time (s)",
            ylabel = "Velocity (m/s)",
            label = "",
            linewidth=2,
            color = RGB(0.337, 0.631, 0.749),
            ytickfontcolor = RGB(0.337, 0.631, 0.749))
    p2 = twinx(p)
        plot!(p2, times, positions,
            ylabel = "Position (m)",
            label = "",
            linewidth=2,
            color = RGB(0.2, 0.468, 0.2),
            ytickfontcolor = RGB(0.2, 0.468, 0.2))
    
    #save plot
    #savefig("C:\\Users\\xcmad\\Documents\\Homework\\ME EN 335\\CoilGunProject\\varying_cuttof_time_show_position")

    plot!(display = true)
end


#this is the code that actually runs
#solution = optimize_coil_gun() #This finds the optimal solution

#For some reason on my machine, I can only create one plot per run. so I comment out everything but the one I want

#plot_a_sol(solution) #this plots what the state space for that solution looks like

plot_vary_cutoff()#this plots a what-if we varyied the cutoff time

#these two lines let us plot what would happen if we had a cutoff time at selected_x
#selected_x = [0.1]
#plot_a_sol(selected_x)