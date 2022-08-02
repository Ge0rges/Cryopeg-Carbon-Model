'''
This file contains master functions for the different cryopeg brine scenarios.
Run the method that corresponds to the scenario you'd like to simulate.
All relevant paramaters will be calculated for that scenario including:
    - Maintenance energy high and low bounds
    - Model output given the calculated maintenance energy
    - Brine expansion factor
The main function of this file executes all scenarios and the sensitivity analytsis.
Call run_all_analysis() with custom paramaters to create your own scenario.
'''

from analysis import *
from plots import *
from utils import *
import numpy as np


def cb1_common_scenario(do_sensitivity_analysis=False, use_me_growth_rate=False, use_me_lower_bound=True):
    dry_sediment_density = 2.625 * 10**6  # ug/ml average density density of kaolinite and sand
    expansion_factor = 0.0905  # Assumed expansion factor of porewater from liquid to solid. (French & Shur 2010)
    beo_volumetric_ice_content = 0.731  # L solid porewater/L permafrost from Go Iwahana (unpublished, but see Meyer et al. 2010)
    beo_poc = 0.0232  # ug C/ug dry sed (our data)
    beo_doc = 51.2323  # mg C/L thawed porewater (our data)
    cb1_18_poc = 1.24*10**4  # in uM from Cooper et al. 2019
    cb1_18_doc = 1.02*10**5  # in uM from Cooper et al. 2019

    # Conversions
    # Converting from ug C/ug dry sediment to fg C/ml permafrost
    # L dry sediment / L permafrost = (1 - L solid porewater / L permafrost)
    # ug C/ug dry sed * ug dry sediment/ml dry sediment * ml/dry sediment/ml permafrost
    beo_poc *= dry_sediment_density * (1 - beo_volumetric_ice_content) * 10**9

    # Converting from mg C/L thawed porewater * L thawed porewater/L permafrost to fg C/ml
    beo_doc *= beo_volumetric_ice_content * (1-expansion_factor) * 10**9

    # Converting from micromolar C to fg C/ml using carbon molar mass of 12.011
    cb1_18_poc = convert_micromolar_carbon_to_fgC_per_ml(cb1_18_poc)
    cb1_18_doc = convert_micromolar_carbon_to_fgC_per_ml(cb1_18_doc)

    # PARAMATERS
    start_carbon = beo_poc + beo_doc
    end_carbon = cb1_18_poc + cb1_18_doc

    print("CB1 start carbon in micromolar: " + str(convert_fgC_per_ml_to_micromolar_carbon(start_carbon)))
    print("CB1 end carbon in micromolar: " + str(convert_fgC_per_ml_to_micromolar_carbon(end_carbon)))

    start_cell = 10**5  # Average cell/ml order of magnitue for sea-ice (Cooper et al. 2019)
    observed_end_cell = 5.7 * 10**6  # Average cell/ml order of magnitude (Cooper et al. 2019)
    timespan = 40000  # Age in years based on carbon dating (Iwanaha et al. 2021)
    growth_rate = 0.06  # Marinobacter aerobic growth rate in-situ based on lab experiments (unpublished), 50% anaerobic penalty
    organic_carbon_per_cell = 15.7  # Litterature based value (Nguyen & Maranger 2011)
    inorganic_carbon_per_cell = 0  # IC not taken into account yet
    inorganic_carbon_fixation_factor = 0  # IC not taken into account yet

    return run_analysis(start_carbon, end_carbon, start_cell, observed_end_cell, timespan, growth_rate,
                        organic_carbon_per_cell, inorganic_carbon_per_cell, inorganic_carbon_fixation_factor,
                        0, do_sensitivity_analysis, use_me_growth_rate, use_me_lower_bound)


def cb4_scenario(do_sensitivity_analysis=False, use_me_growth_rate=False, use_me_lower_bound=True):
    dry_sediment_density = 2.625 * 10**6  # ug/ml average density density of kaolinite and sand
    expansion_factor = 0.0905  # Assumed expansion factor of porewater from liquid to solid. (French & Shur 2010)
    cb4_volumetric_ice_content = 0.527  # L solid porewater/L permafrost from Go Iwahana (unpublished, but see Meyer et al. 2010)
    cb4_poc = 0.0136  # ug C/ug dry sed (our data)
    cb4_doc = 1286.81  # mg C/L thawed porewater (our data)
    cb4_brine_poc = 4.14 * 10**3  # in uM (Cooper et al. 2019)
    cb4_brine_doc = 8.50 * 10**4  # in uM (Cooper et al. 2019)

    # Conversions
    # Converting from ug C/ug dry sediment to fg C/ml permafrost
    # L dry sediment / L permafrost = (1 - L solid porewater / L permafrost)
    # ug C/ug dry sed * ug dry sediment/ml dry sediment * ml/dry sediment/ml permafrost
    cb4_poc *= dry_sediment_density * (1 - cb4_volumetric_ice_content) * 10**9

    # Converting from mg C/L thawed porewater * L thawed porewater/L permafrost  to fg C/ml
    cb4_doc *= cb4_volumetric_ice_content * (1-expansion_factor) * 10**9

    # Converting from micromolar C to fg C/ml using carbon molar mass of 12.011
    cb4_brine_poc = convert_micromolar_carbon_to_fgC_per_ml(cb4_brine_poc)
    cb4_brine_doc = convert_micromolar_carbon_to_fgC_per_ml(cb4_brine_doc)

    # PARAMATERS
    start_carbon = cb4_poc + cb4_doc
    end_carbon = cb4_brine_poc + cb4_brine_doc

    print("CB4 start carbon in micromolar: " + str(convert_fgC_per_ml_to_micromolar_carbon(start_carbon)))
    print("CB4 end carbon in micromolar: " + str(convert_fgC_per_ml_to_micromolar_carbon(end_carbon)))

    start_cell = 10**5  # Average cell/ml order of magnitue for sea-ice (Cooper et al. 2019)
    observed_end_cell = 1.14 * 10**7  # cell/ml for CB4_18 (Cooper et al. 2019)
    timespan = 40000  # Age in years based on carbon dating (Iwanaha et al. 2021)
    growth_rate = 0.016  # Psychrobacter cryohalolentis growth rate in-situ based on lab experiments (Bakermans et al. 2003)
    organic_carbon_per_cell = 54.04  # Took average P. cryohalolentis size: 0.365014 um3 (Bakermans et al. 2006), and carbon conversion factor as 148 fg C/um3 (Kirchman et al. 2009)
    inorganic_carbon_per_cell = 0  # IC not taken into account yet
    inorganic_carbon_fixation_factor = 0  # IC not taken into account yet

    return run_analysis(start_carbon, end_carbon, start_cell, observed_end_cell, timespan, growth_rate,
                        organic_carbon_per_cell, inorganic_carbon_per_cell, inorganic_carbon_fixation_factor,
                        0, do_sensitivity_analysis, use_me_growth_rate, use_me_lower_bound)


def cbiw_scenario(do_sensitivity_analysis=False, use_me_growth_rate=False, use_me_lower_bound=True):
    """
    CBIW is a intra-ice cryopeg brine. It is thought to have migrated from sediment to massive ice at around 11000 years
    ago. Therefore, this scenario is unique in including an injection of carbon at 11 000 years. We use the lower bound
    estimate of how much carbon was injected to be the organic carbon measured in the massive ice around the brine today.
    It is likely that more carbon was added than our estimate, however given our current understanding it is impossible
    to quantify precisely how much carbon was added into the brine during its migration.
    """
    dry_sediment_density = 2.625 * 10**6  # ug/ml average density density of kaolinite and sand
    expansion_factor_porewater = 0.0905  # Assumed expansion factor of porewater from liquid to solid. (French & Shur 2010)
    beo_volumetric_ice_content = 0.731  # L solid porewater/L permafrost from Go Iwahana (unpublished, but see Meyer et al. 2010)
    expansion_factor_ice = 0.08042  # % Expansion factor of pure water (massive ice) from liquid to solid
    beo_poc = 0.0232  # ug C/ug dry sed (our data)
    beo_doc = 51.2323  # mg C/L thawed porewater (our data)
    massive_ice_poc = 4.5  # ug C/ml thawed massive ice (Collangelo-Lillis 2016)
    massive_ice_doc = 0.65  # μg C/ml thawed massive ice. EPS as a proxy for DOC. (Colangelo-Lillis et al. 2016)
    cbiw_brine_poc = 1.98 * 10**3  # in uM (Cooper et al. 2019)
    cbiw_brine_doc = 3 * 10**4  # in uM (Cooper et al. 2019)

    # Conversions
    # Converting from ug C/ug dry sediment to fg C/ml permafrost
    # L dry sediment / L permafrost = (1 - L solid porewater / L permafrost)
    # ug C/ug dry sed * ug dry sediment/ml dry sediment * ml/dry sediment/ml permafrost
    beo_poc *= dry_sediment_density * (1 - beo_volumetric_ice_content) * 10**9

    # Converting from mg C/L thawed porewater * L thawed porewater/L permafrost to fg C/ml
    beo_doc *= beo_volumetric_ice_content * (1 - expansion_factor_porewater) * 10 ** 9

    # Converting from ug C/ml thawed massive ice to fg C/ml massive ice
    massive_ice_poc *= massive_ice_poc * (1-expansion_factor_ice) * 10**9
    massive_ice_doc *= massive_ice_doc * (1-expansion_factor_ice) * 10**9

    # Converting from micromolar C to fg C/ml using carbon molar mass of 12.011
    cbiw_brine_poc = convert_micromolar_carbon_to_fgC_per_ml(cbiw_brine_poc)
    cbiw_brine_doc = convert_micromolar_carbon_to_fgC_per_ml(cbiw_brine_doc)

    # PARAMATERS
    start_carbon = beo_poc + beo_doc
    punctual_organic_carbon_addition = [(massive_ice_poc + massive_ice_doc, 11000)]  # Add massive ice carbon at 11000y.
    end_carbon = cbiw_brine_poc + cbiw_brine_doc

    print("CBIW start carbon in micromolar: " + str(convert_fgC_per_ml_to_micromolar_carbon(start_carbon)))
    print("CBIW end carbon in micromolar: " + str(convert_fgC_per_ml_to_micromolar_carbon(end_carbon)))

    start_cell = 10 ** 5  # Average cell/ml order of magnitue for sea-ice  (Cooper et al. 2019)
    observed_end_cell = 1.30 * 10 ** 8  # average cell/ml of CBIW (Cooper et al. 2019)
    timespan = 40000  # Age in years based on carbon dating (Iwanaha et al. 2021)
    growth_rate = 0.06  # Marinobacter aerobic growth rate in-situ based on lab experiments (unpublished), 50% anaerobic penalty
    organic_carbon_per_cell = 15.7  # Litterature based value (Nguyen & Maranger 2011)
    inorganic_carbon_per_cell = 0  # IC not taken into account yet
    inorganic_carbon_fixation_factor = 0  # IC not taken into account yet

    return run_analysis(start_carbon, end_carbon, start_cell, observed_end_cell, timespan, growth_rate,
                        organic_carbon_per_cell, inorganic_carbon_per_cell, inorganic_carbon_fixation_factor,
                        punctual_organic_carbon_addition, do_sensitivity_analysis, use_me_growth_rate, use_me_lower_bound)


def run_analysis(start_carbon, end_carbon, start_cell, observed_end_cell, timespan, growth_rate, org_carbon_per_cell, inorg_carbon_per_cell, inorg_carbon_fixation_factor, punctual_organic_carbon_addition=0, do_sensitivity_analysis=False, use_me_growth_rate=False, use_me_lower_bound=True):

    # Convert timespawn from years to days
    timespan *= 365.25

    # When estimating maintenance energy we include all added carbon.
    if punctual_organic_carbon_addition != 0:
        # noinspection PyTypeChecker
        for i, (carbon_addded, time_to_add) in enumerate(punctual_organic_carbon_addition):
            start_carbon += carbon_addded
            punctual_organic_carbon_addition[i] = (carbon_addded, timespan - time_to_add*365.25)  # Convert from years to days, and go from end not from start

    # Estimate bounds for maintenance energy
    me_lower = estimate_me_no_growth(start_carbon, end_carbon, observed_end_cell, timespan)
    me_upper_growth_rate, doubling_time, me_upper = estimate_me_exp_growth(start_carbon, end_carbon, start_cell,
                                                                           observed_end_cell, org_carbon_per_cell,
                                                                           timespan)

    # Build paramater, ivp arrays (ORDER MATTERS)                                      IVP:
    #  carrying_capacity = p[1]                organic_carbon_input = p[2]               organic_carbon_content = u[1]
    #  cell_count = u[3]                       natural_death_fraction = p[4]             inorganic_carbon_content = u[2]
    #  inorganic_carbon_fixing_factor = p[5]   ks = p[6]                                 start_cell_count = u[3]
    #  organic_carbon_content_per_cell = p[7]  inorganic_carbon_content_per_cell = p[8]  timespan = u[4]
    #  mu_max = p[9]                           base_maintenance_per_cell = p[10]
    #  punctual_carbon_added = p[11]           m_prime = p[12]
    growth_rate = me_upper_growth_rate if use_me_growth_rate else growth_rate  # Switch out growth rate to calcualted one
    me = me_lower if use_me_lower_bound else me_upper  # Switch out maintenance energy based on paramater
    death_rate = 0 if use_me_growth_rate else 0.001  # Switch out death rate to 0 for net calculated growth rate

    p = [10**9, 0, 0, death_rate, inorg_carbon_fixation_factor, 882000000, org_carbon_per_cell,
         inorg_carbon_per_cell, growth_rate, me, punctual_organic_carbon_addition, 0]
    ivp = [start_carbon, 0, start_cell, timespan]

    # Run the model through PyJulia using the lower ME
    S, I, N, t = run_model(p, ivp)

    # Calculate brine expansion
    _, _, ratio_dimensions, _, _ = calculate_brine_expansion(start_carbon, (S[-1] - S[0])/timespan)

    # Run the sensitivity analysis
    ST, S1 = None, None
    if do_sensitivity_analysis:
        ST, S1 = run_sensitivity_analysis(p, default_paramater_bounds, ivp)

    return (me_lower, me_upper, me_upper_growth_rate, doubling_time), (S, I, N, t), ratio_dimensions, (ST, S1)


def log_results(results, savefig=True):
    # Store for collective plot
    S_array = []
    I_array = []
    N_array = []
    t_array = []
    labels = []

    for result, label in results:
        (me_lower, me_upper, me_growth_rate, doubling_time), (S, I, N, t), ratio_dimensions, (ST, S1) = result

        # Store for collective plot
        S_array.append(S)
        I_array.append(I)
        N_array.append(N)
        t_array.append(t)
        labels.append(label)

        # Log endpoints and ME calculations

        print(label + " end organic carbon is: {}".format(np.format_float_scientific(S[-1], precision=2)))
        print(label + " end cell density is: {}".format(np.format_float_scientific(N[-1], precision=2)))

        if "- Calculated Growth Rate" not in label:  # Same results, don't log twice
            print(label + " Maintenance energy (low, high, low growth rate, doubling_time) in fg C/cell day and /day: ({}, {}. {}, {})"
                  .format(np.format_float_scientific(me_lower, precision=4),
                          np.format_float_scientific(me_upper, precision=4),
                          np.format_float_scientific(me_growth_rate, precision=4),
                          doubling_time/365.25
                          )
                  )
        # print(label + " Brine expansion ratio for linear dimensions is: {.2f} % along each axis per year.".format(ratio_dimensions))

        # Individual plots
        model_fig = plot_model(S, I, N, t, label)
        # model_fig.show()

        if savefig:
            model_fig.savefig("Plots/" + label + "_model.pdf", format="pdf", dpi=500, bbox_inches='tight')

        # Check if sensitivity analysis was done
        if ST is not None:
            sa_fig = plot_sensitivity(ST, S1, default_paramater_names, label)
            # sa_fig.show()

            if savefig:
                sa_fig.savefig("Plots/" + label + "_sa.pdf", format="pdf", dpi=500, bbox_inches='tight')

    # Make collective plot
    colors = ["darkgreen", "red", "violet", "pink", "peru", "black"]
    scenarios_fig = plot_multiple_scenarios(S_array, I_array, N_array, t_array, labels, "Model outputs of all considered scenarios", colors)
    # scenarios_fig.show()

    if savefig:
        scenarios_fig.savefig("Plots/all_scenarios_model.pdf", format="pdf", dpi=500, bbox_inches='tight')


if __name__ == "__main__":
    # Generates all figures and data points presented.
    all_results = []

    all_results.append((cb1_common_scenario(do_sensitivity_analysis=True, use_me_growth_rate=False, use_me_lower_bound=True), "Intra-sediment"))
    all_results.append((cb4_scenario(do_sensitivity_analysis=False, use_me_growth_rate=False, use_me_lower_bound=True), "CB4"))
    all_results.append((cbiw_scenario(do_sensitivity_analysis=False, use_me_growth_rate=False, use_me_lower_bound=True), "Intra-ice"))

    log_results(all_results)

    hg_fig = hypothetical_growth_scenarios()
    hg_fig.savefig("Plots/" + "hypothetical_growth.pdf", format="pdf", dpi=500, bbox_inches='tight')
