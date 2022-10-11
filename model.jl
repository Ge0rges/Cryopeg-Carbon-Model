using DifferentialEquations
using DiffEqSensitivity
using ForwardDiff
using GlobalSensitivity
using Statistics


# Define model domain. Reject any values below 0.
function is_invalid_domain(u,p,t)
    return u[1] < 0 || u[2] < 0 || u[3] < 0 || u[4] < 0
end

# Defines a problem object and solves it.
function solve_model(p, u0)
    # Set the time timespan
    tspan = (0.0, last(u0))

    # Build a callback to introduce a carbon addition as it is a discontinuity
    stops = []
    carbon = []
    punctual_organic_carbon_addition = p[end] == -1 ? [] : p[end]
    for (time_to_add, (pOC_to_add, dOC_to_add)) in punctual_organic_carbon_addition
        push!(stops, time_to_add)
        push!(carbon, (convert(Float64, pOC_to_add), convert(Float64, dOC_to_add)))
    end

    additions = Dict(Pair.(stops, carbon))

    addition_condition(u,t,integrator) = any(x->t==x, stops)  # If we hit any of the stops return true
    function addition!(integrator)
        integrator.u[1] += additions[integrator.t][1]
        integrator.u[2] += additions[integrator.t][2]

        if integrator.u[3] < 1
            integrator.u[3] = 1  # Add a viable cell if none exist
        end
    end
    carbon_add_cb = DiscreteCallback(addition_condition, addition!)

    # Convert p and remove puncutal_punctual_organic_carbon_addition
    p = convert(Array{Float64}, p[1:end-1])

    # Build the ODE Problem and Solve
    prob = ODEProblem(model, u0[1:end-1], tspan, p)
    sol = solve(prob, Rosenbrock23(), callback=carbon_add_cb, tstops=stops, isoutofdomain=is_invalid_domain, maxiters=1e6)

    return sol
end


# Runs the model using solve_model and packages results nicely.
function run_model(p, u0)
    u0 = convert(Array{Float64}, u0)

    sol = solve_model(p, u0)

    # Return the solution array - [pOC, dOC, IC, Cells, t]
    return [[x[1] for x in sol.u], [x[2] for x in sol.u], [x[3] for x in sol.u], [x[4] for x in sol.u], sol.t]
end


# Runs the sensitivity analysis using Sobol method.
function run_sensitivity_analysis(p, p_bounds, u0)
    u0 = convert(Array{Float64}, u0)

    p_bounds = [p_bounds[i, :] for i in 1:size(p_bounds, 1)]

    # Define a function that remakes the problem and gets its result
    f1 = function (p)
        sol = solve_model(p, u0)
        [last(sol[1,:] * 1e200)]  # This hacky solution to handle tiny (10^-300) values sucks.
    end

    # Run GSA
    sobol_result = GlobalSensitivity.gsa(f1, Sobol(), p_bounds, Ei_estimator=:Sobol2007, samples=2^12)
    return (sobol_result.ST[1,:], sobol_result.S1[1,:])
end


# Implements the differential equations that define the model
function model(du,u,p,t)
    # Paramaters
    mu_max = p[1]
    maintenance_per_cell = p[2]
    dOC_per_cell = p[3]

    carrying_capacity = p[4]
    pOC_input_rate = p[5]
    dOC_input_rate = p[6]
    inorganic_carbon_input_rate = p[7]
    inorganic_carbon_fixing_rate = p[8]
    inorganic_carbon_per_cell = p[9]
    eea_rate = p[10]
    Ks = p[11]

    # Load state conditions
    pOC_content = u[1]
    dOC_content = u[2]
    inorganic_carbon_content = u[3]
    cell_count = u[4]

    ## CELL COUNT
    # Growth
    growth = mu_max * (dOC_content / (Ks + dOC_content)) * cell_count * (1 - cell_count / carrying_capacity)

    # Organic carbon requirement
    required_dOC_per_cell = maintenance_per_cell
    required_dOC = required_dOC_per_cell * cell_count

    # Starvation Deaths
    dOC_missing = max(required_dOC - dOC_content, 0)
    starvation_deaths = dOC_missing == 0 ? 0 : dOC_missing / required_dOC_per_cell

    # Total Deaths
    deaths =  starvation_deaths

    ## CARBON
    dOC_consumption = required_dOC_per_cell * (cell_count - deaths)
    fixed_carbon = inorganic_carbon_fixing_rate * inorganic_carbon_content

    # Particulate Organic carbon
    du[1] = pOC_input_rate - eea_rate*pOC_content*cell_count
    # Do i need to multiply the rate by t?
    #TODO

    # Dissolved Organic carbon
    du[2] = dOC_input_rate + dOC_per_cell * (deaths - growth) - dOC_consumption + eea_rate*pOC_content*cell_count + fixed_carbon

    # Inorganic carbon
    du[3] = inorganic_carbon_input_rate + inorganic_carbon_per_cell * (deaths - growth) + required_dOC - fixed_carbon

    # Net cell count change
    du[4] = growth - deaths
end
