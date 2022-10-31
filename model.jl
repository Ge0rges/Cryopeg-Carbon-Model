using DifferentialEquations
using DiffEqSensitivity
using ForwardDiff
using GlobalSensitivity
using Statistics
using Printf


# Defines a problem object and solves it.
function solve_model(p, u0)
    # Build a callback to introduce a carbon addition as it is a discontinuity
    stops = []
    carbon = []
    punctual_organic_carbon_addition = p[end] == -1 ? [] : p[end]
    for (time_to_add, (pOC_to_add, dOC_to_add)) in punctual_organic_carbon_addition
        push!(stops, time_to_add)
        push!(carbon, (convert(Float64, pOC_to_add), convert(Float64, dOC_to_add)))
    end

    # Convert p and remove puncutal_punctual_organic_carbon_addition
    p = convert(Array{Float64}, p[1:end-1])

    # Callback for punctual addition - causes a discontinuity
    additions = Dict(Pair.(stops, carbon))
    function addition!(integrator)
        integrator.u[1] += additions[integrator.t][1]
        integrator.u[2] += additions[integrator.t][2]

        if integrator.u[3] < 1
            integrator.u[3] = 1  # Add a viable cell if none exist
        end
    end
    carbon_add_cb = PresetTimeCallback(stops, addition!)

    # Callback for the max() function in the model - causes a discontinuity.
    # If max equation changes, this condition will have to change.
    max_condition(u, t, integrator) = p[2] * u[4] - u[2]
    do_nothing(integrator) = nothing
    max_cb = ContinuousCallback(max_condition, do_nothing)

    # Min condition for EEA rate
    min_condition(u, t, integrator) = p[10] * u[4] - u[1]
    min_cb = ContinuousCallback(min_condition, do_nothing)

    # Set things to 0 if they are less than 1e-100
    zero_condition(u, t, integrator) = u[1] < 1e-100 || u[2] < 1e-100 || u[3] < 1e-100 || u[4] < 1e-100
    function zero_out!(integrator)
        for i in 1:4
            if integrator.u[i] < 1e-100
                integrator.u[i] = 0
            end
        end
    end
    zero_cb = DiscreteCallback(zero_condition, zero_out!)

    # Callback list
    cbs = CallbackSet(carbon_add_cb, max_cb, min_cb, zero_cb)

    # Out of domain function
    is_invalid_domain(u,p,t) = u[1] < 0 || u[2] < 0 || u[3] < 0 || u[4] < 0

    # Build the ODE Problem and Solve
    prob = ODEProblem(model, u0[1:end-1], (0.0, last(u0)), p)
    sol = solve(prob, Rosenbrock23(), callback=cbs, isoutofdomain=is_invalid_domain, maxiters=1e6)

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
function run_sensitivity_analysis(p_bounds, u0, carbon_output)
    u0 = convert(Array{Float64}, u0)

    p_bounds = [p_bounds[i, :] for i in 1:size(p_bounds, 1)]

    # Define a function that remakes the problem and gets its result. Called for each sample.
    f1 = function (p)
        sol = solve_model(p, u0)
        sol[:, end]
    end

    # Run GSA
    sobol_result = GlobalSensitivity.gsa(f1, eFAST(), p_bounds, samples=2^10)
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
    Kd = p[11]

    # Load state conditions
    pOC_content = u[1]
    dOC_content = u[2]
    inorganic_carbon_content = u[3]
    cell_count = u[4]

    ## CELL COUNT
    # Growth
    growth = mu_max * (dOC_content / (Kd + dOC_content)) * cell_count * (1 - cell_count / carrying_capacity)

    # Organic carbon requirement
    required_dOC_per_cell = maintenance_per_cell
    required_dOC = required_dOC_per_cell * cell_count

    # Starvation Deaths
    dOC_missing = max(required_dOC - dOC_content, 0)
    starvation_deaths = dOC_missing / required_dOC_per_cell

    # Total Deaths
    deaths =  starvation_deaths

    ## CARBON
    dOC_consumption = required_dOC_per_cell * (cell_count - deaths)
    fixed_carbon = inorganic_carbon_fixing_rate * inorganic_carbon_content

    # EEA rate
    eea_rate = min(eea_rate*cell_count, pOC_content)

    # Particulate Organic carbon
    du[1] = pOC_input_rate - eea_rate*cell_count

    # Dissolved Organic carbon
    du[2] = dOC_input_rate + dOC_per_cell * (deaths - growth) - dOC_consumption + eea_rate*cell_count + fixed_carbon

    # Inorganic carbon
    du[3] = inorganic_carbon_input_rate + inorganic_carbon_per_cell * (deaths - growth) + required_dOC - fixed_carbon

    # Net cell count change
    du[4] = growth - deaths
end
