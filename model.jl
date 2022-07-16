using DifferentialEquations
using DiffEqSensitivity
using ForwardDiff
using GlobalSensitivity
using Statistics
using Printf

function run_model(p, u0)
    u0 = convert(Array{Float64}, u0)

    # Check len of p and u0 to be 10 and 3
    @assert size(p)[1] == 11 "Wrong paramaters passed"
    @assert size(u0)[1] == 4 "Wrong ivp passed"

    # Build the ODE Problem and Solve
    tspan = (0.0, last(u0))
    prob = ODEProblem(model, u0[1:3], tspan, p)
    sol = solve(prob, reltol=1e-9, abstol=1e-9)

    # Return the solution array - OC, IC, Cells, t
    return [[x[1] for x in sol.u], [x[2] for x in sol.u], [x[3] for x in sol.u], sol.t]
end


function run_sensitivity_analysis(p_bounds, u0)
    u0 = convert(Array{Float64}, u0)

    # Check len of p and u0 to be 10 and 3
    @assert size(p)[1] == 11 "Wrong paramaters passed"
    @assert size(u0)[1] == 4 "Wrong ivp passed"

    # Build the initial ODE Problem
    tspan = (0.0, last(u0))
    prob = ODEProblem(model, u0[1:3], tspan, p)

    # Define a function that remakes the problem and gets its result
    f1 = function (p)
      prob1 = remake(prob;p=p)
      sol = solve(prob1, Rosenbrock23(); saveat=collect(range(0, stop=duration, length=200)), maxiters=Int(1e6))
      [sol[1,:]]
    end

    # Run GSA
    sobol_result = GlobalSensitivity.gsa(f1, Sobol(), p_bounds, N=1000)
    return (sobol_result.ST[1,:], sobol_result.S1[1,:])
end


function model(du,u,p,t)
    # PARAMATERS
    # Load paramaters
    carrying_capacity = p[1]
    organic_carbon_input = p[2]
    inorganic_carbon_input = p[3]
    natural_death_fraction = p[4]
    inorganic_carbon_fixing_factor = p[5]
    ks = p[6]
    organic_carbon_content_per_cell = p[7]
    inorganic_carbon_content_per_cell = p[8]
    mu_max = p[9]
    base_maintenance_per_cell = p[10]
    m_prime = p[11]
    @assert m_prime == "empty" "m_prime not implemented"

    # Load state conditions
    organic_carbon_content = u[1]
    inorganic_carbon_content = u[2]
    cell_count = u[3]

    ## CELL COUNT
    # Growth
    growth = mu_max * (organic_carbon_content / (ks + organic_carbon_content)) * cell_count * (1 - cell_count / carrying_capacity)

#     # Specific growth rate
#     next_cell_count = cell_count + growth
#     specific_growth_rate = 0
#     if next_cell_count > 0 && cell_count > 0
#         specific_growth_rate = max((np.log(next_cell_count) - np.log(cell_count)), 0)
#     end

    # Organic carbon requirement
    required_organic_carbon_per_cell = base_maintenance_per_cell  #+ m_prime * (1 - specific_growth_rate/mu_max)
    required_organic_carbon = required_organic_carbon_per_cell * cell_count

    # Starvation Deaths
    organic_carbon_missing = max(required_organic_carbon - organic_carbon_content, 0)
    starvation_deaths = organic_carbon_missing == 0 ? 0 : organic_carbon_missing / required_organic_carbon_per_cell

    # Natural Death
    natural_deaths = natural_death_fraction * cell_count - natural_death_fraction * starvation_deaths

    # Total Deaths
    deaths = natural_deaths + starvation_deaths

    # Net cell count change
    du[3] = growth - deaths

    ## CARBON
    carbon_consumption = required_organic_carbon_per_cell * (cell_count - deaths)
    fixed_carbon = inorganic_carbon_fixing_factor * inorganic_carbon_content

    # Inorganic carbon
    du[2] = inorganic_carbon_input + inorganic_carbon_content_per_cell * (deaths - growth) + carbon_consumption - fixed_carbon

    # Organic carbon
    du[1] = organic_carbon_input + organic_carbon_content_per_cell * (deaths - growth) - carbon_consumption + fixed_carbon
end

