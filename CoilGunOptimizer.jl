using SNOW
using OrdinaryDiffEq

#returns the voltage as if it was on until it was turned off at a variable time
function voltage(max_v, cutoff_time, curr_time)
    if curr_time >= cutoff_time
        return 0.0
    else
        return max_v
    end
end

#ODE for EOM of a single stage coil gun
function single_stage(du, u, p, t)
    L_a = p[1]
    R = p[2]
    K_e = p[3]#assuming K_e and K_f are the same
    m = p[4]
    mu_k = p[5]
    g = 9.81
    max_v = p[6]
    cutoff_time = p[7]
    V_a = voltage(max_v, cutoff_time, t)

    du[1] = -(R/L_a)*u[1] - (K_e/L_a)*u[3] + (V_a/L_a)
    du[2] = u[3]
    du[3] = (K_e/m)*u[1] + mu_k*g
end

#objective function to optimize coil gun
#let's vary: cutoff_time
function ob!(g, x)
    #set up variables
    p = [1.0, 1.5, 0.1, 0.01, 0.03, 1500.0, x[1]] #important variables
    t = (0.0, 2.0) #time
    u0 = [0.0, -0.1, 0.0] #initial conditions
    
    #objective (velocity when x = pos)
    pos = 0.1
    
    #set up end condition for simulation
    condition(u,t,integrator) = u[2] - 0.1 # Is zero when u[2] = 0.1
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(condition, affect!)
    #simulate
    prob = ODEProblem(single_stage, u0, t, p) #set up the problem
    sol = solve(prob, callback = cb); #solve the problem

    # get the velocity when the position is pos
    f = -sol[length(sol)][3]

    #constraints
    g[1] = x[1] #the voltage can't be turned off at a negative value nor one larger than 2s
    
    #return
    return f
end

function optimize_coil_gun()
    x0 = [0.0]  # starting point
    lx = [0.0]  # lower bounds on x
    ux = [2.0]  # upper bounds on x
    ng = 1  # number of constraints
    lg = [0.0]  # lower bounds on g
    ug = [2.0]  # upper bounds on g
    options = Options(solver=IPOPT())  # choosing IPOPT solver

    xopt, fopt, info = minimize(ob!, x0, ng, lx, ux, lg, ug, options)

    println("xstar = ", xopt)
    println("fstar = ", fopt)
    println("info = ", info)
end

optimize_coil_gun()