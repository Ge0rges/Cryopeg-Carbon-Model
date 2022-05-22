def sympy_model(paramaters, ivp):
    sympy.init_printing()

    N, S, I = symbols("N S I", cls=Function)
    mu, t, m_base, m_p, u_max, YG, delta, KS, N_max, d = symbols("mu t m_base m_p mu_max Y_G delta K_s N_max d")
    Sin, alphaS, Ifr, alphaI, Iin = symbols("S_in alpha_S I_fr I_in alpha_I")

    # m_base = paramaters["base_maintenance_per_cell"]
    # m_p = paramaters["m_prime"]
    # u_max = paramaters["mu_max"]
    # YG = paramaters["maximum_growth_yield"]
    # delta = paramaters["growth_penalty"]
    # KS = paramaters["ks"]
    # N_max = paramaters["carrying_capacity"]
    # d = paramaters["natural_death_fraction"]
    # Sin = paramaters["organic_carbon_input"]
    # alphaS = paramaters["organic_carbon_content_per_cell"]
    # Ifr = paramaters["inorganic_carbon_fixing_factor"]
    # alphaI = paramaters["inorganic_carbon_content_per_cell"]
    # Iin = paramaters["inorganic_carbon_input"]

    m = m_base  # + m_p*(1 - mu/u_max) + mu/YG
    G = (1 - delta) * u_max * t * S(t) / (KS + S(t)) * N(t) * (1 - N(t) / N_max)
    D = d * (N(t) - Max(m * N(t) - S(t), 0) / m) + Max(m * N(t) - S(t), 0) / m

    dndt = G - D
    dsdt = Sin + alphaS * (D - G) - m * (N(t) - D) + Ifr * I(t)
    didt = Iin + alphaI * (D - G) + m * (N(t) - D) - Ifr * I(t)

    system = [Eq(N(t).diff(t), dndt), Eq(S(t).diff(t), dsdt), Eq(I(t).diff(t), didt)]
    ics = {N(0): ivp["initial_cell_count"], S(0): ivp["initial_carbon_content"], I(0): ivp["initial_inorganic_carbon_content"]}

    system = [dndt, dsdt, didt]
    for eq in system:
        sympy.pprint(diff(eq, m))
        print("END##############")

    # dsolve_system(system, [N(t), S(t), I(t)], ics=ics, t=t)
