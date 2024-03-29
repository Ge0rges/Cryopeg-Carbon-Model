using DifferentialEquations
using DiffEqSensitivity
using ForwardDiff
using GlobalSensitivity
using Statistics
using Printf


# Runs the model using solve_model and packages results nicely.
function run_model(p, u0)
    cbs = make_callbacks(sensitivity_analysis=false, p=p)

    # Convert p and u0, and remove puncutal_punctual_organic_carbon_addition from p
    p = convert(Array{Float64}, p[1:end-1])
    u0 = convert(Array{Float64}, u0)

    prob = ODEProblem(model, u0[1:end-1], (0.0, last(u0)), p)
    sol = solve(prob, Rosenbrock23(), callback=cbs, maxiters=1e6)

    # Return the solution array - [pOC, dOC, DIC, Cells, t]
    return [[x[1] for x in sol.u], [x[2] for x in sol.u], [x[3] for x in sol.u], [x[4] for x in sol.u], sol.t]
end


# Runs the sensitivity analysis using Sobol method.
function run_sensitivity_analysis(p_bounds, u0_length)
    # In Sensitivity analysis, parameters include the initial conditions.
    p_bounds = [p_bounds[i, :] for i in 1:size(p_bounds, 1)] # Fix array shape
    p_bounds = convert(Array{Array{Float64}}, p_bounds[1:end])

    lb = [p_bounds[i][1] for i in 1:size(p_bounds, 1)]
    ub = [p_bounds[i][2] for i in 1:size(p_bounds, 1)]

    # Timespan must stay constant
    @assert ub[end] == lb[end] "Timespan must not be subjected to SA"

    # Get cbs
    cbs = make_callbacks(sensitivity_analysis=true)

    # Define a function that remakes the problem and gets its result. Called for each sample.
    f1 = function (p)
        # u0 is the last n indices
        prob = ODEProblem(model, p[end-u0_length+1:end-1], (0.0, p[end]), p[1:end-u0_length])
        sol = solve(prob, Rosenbrock23(); callback=cbs, maxiters=1e6, saveat=ub[end])
        sol[:, end]
    end

    # Run GSA
    sobol_result = GlobalSensitivity.gsa(f1, Sobol(order=[0, 1], nboot=5, conf_level=0.95), p_bounds, samples=2^20)

    return (sobol_result.ST, sobol_result.S1, sobol_result.ST_Conf_Int, sobol_result.S1_Conf_Int)
end


# Defines a problem object and solves it.
function make_callbacks(;sensitivity_analysis=false, p=nothing)
    # Build a callback to introduce a carbon addition as it is a discontinuity
    carbon_add_cb = nothing
    if !sensitivity_analysis
        stops = []
        carbon = []
        punctual_organic_carbon_addition = p[end] == -1 ? [] : p[end]
        for (time_to_add, (pOC_to_add, dOC_to_add)) in punctual_organic_carbon_addition
            push!(stops, time_to_add)
            push!(carbon, (convert(Float64, pOC_to_add), convert(Float64, dOC_to_add)))
        end

        # Callback for punctual addition - causes a discontinuity
        additions = Dict(Pair.(stops, carbon))
        function addition!(integrator)
            integrator.u[1] += additions[integrator.t][1]
            integrator.u[2] += additions[integrator.t][2]

            if integrator.u[4] < 1
                integrator.u[4] = 1  # Add a viable cell if none exist
            end
        end
        carbon_add_cb = PresetTimeCallback(stops, addition!)
    end

    # Callback for the max() function in the model - causes a discontinuity.
    # If max equation changes, this condition will have to change.
    max_condition(u, t, integrator) = integrator.p[2] * u[4] - u[2]
    do_nothing(integrator) = nothing
    max_cb = ContinuousCallback(max_condition, do_nothing)

    # Min condition for EEA rate
    eea_min_condition(u, t, integrator) = integrator.p[10] * u[4] - u[1]
    eea_min_cb = ContinuousCallback(eea_min_condition, do_nothing)

    # Min condition for DIC fixation rate
    ic_min_condition(u, t, integrator) = integrator.p[8] - u[3]
    ic_min_cb = ContinuousCallback(ic_min_condition, do_nothing)

    # Set things to 0 if they are less than 1e-100
    zero_condition(u, t, integrator) = any(x -> x < 1e-20, u)
    function zero_out!(integrator)
        for i in 1:4
            if integrator.u[i] < 1e-20
                integrator.u[i] = 0
            end
        end
    end
    zero_cb = DiscreteCallback(zero_condition, zero_out!)

    # Callback list
    return CallbackSet(carbon_add_cb, max_cb, eea_min_cb, ic_min_cb, zero_cb, PositiveDomain())
end


# Implements the differential equations that define the model
function model(du, u, p, t)
    # Parameters
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
    deaths = starvation_deaths

    ## CARBON
    dOC_consumption = required_dOC_per_cell * (cell_count - deaths)
    fixed_carbon = min(inorganic_carbon_fixing_rate,  inorganic_carbon_content)

    # EEA rate
    eea_removal = min(eea_rate*cell_count, pOC_content)

    # Particulate Organic carbon
    du[1] = pOC_input_rate - eea_removal

    # Dissolved Organic carbon
    du[2] = dOC_input_rate + dOC_per_cell * (deaths - growth) - dOC_consumption + eea_removal + fixed_carbon

    # Inorganic carbon
    du[3] = inorganic_carbon_input_rate + inorganic_carbon_per_cell * (deaths - growth) + required_dOC - fixed_carbon

    # Net cell count change
    du[4] = growth - deaths
end
