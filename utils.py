"""
Some helper functions. Includes classes to encapsulate analysis results.
"""

from decimal import Decimal

import pandas
import numpy as np


def convert_micromolar_carbon_to_fgC_per_ml(micromolar):
    """Converts a micromolar of carbon value to femtograms of carbon per ml."""
    return (micromolar * 12.011 * 10**-6 * 10**15)/1000


def convert_fgC_per_ml_to_micromolar_carbon(fgC_per_ml):
    """Converts a femtograms of carbon per ml value to micromolar of carbon."""
    return (fgC_per_ml * 1000)/(12.011 * 10**-6 * 10**15)


class MaintenanceResult:
    """
    Class that encapsulates the result of a maintenance energy estimation.
    """
    lower_bound_me: float = None
    upper_bound_me: float = None
    minimum_growth_rate: float = None
    minimum_doubling_time: float = None


class SensitivityResult:
    """
    Class that encapsulates the results of a sensitivity analysis.
    """
    total_sobol_indices: [float] = None
    first_order_sobol_indices: [float] = None
    total_conf_int: [[float]] = None
    first_order_conf_int: [[float]] = None
    variables: [str] = ["P", "D", "I", "N", "Mean"]
    assert len(variables) == len(total_sobol_indices) + 1

    def get_dataframe(self, scenario):
        """
        Returns model results as a dataframe. Removes values that are exactly 0 (SA wasn't run).
        """
        p_names = scenario._paramater_names
        ST = self.total_sobol_indices
        S1 = self.first_order_sobol_indices
        st_err = self.total_conf_int
        s1_err = self.first_order_conf_int
       
        df = []
        for row in range(len(p_names[0:-len(scenario.get_julia_ordered_ivp())])):
            if all(ST[i, row] == 0 for i in range(4)):
                continue

            for i, output in enumerate(self.variables):
                if output == "Mean":
                    print(ST[:, row])
                    df.append([p_names[row], output, np.mean(ST[:, row]), np.mean(S1[:, row]), np.mean(st_err[:, row]), np.mean(s1_err[:, row])])

                else:

        df = pandas.DataFrame(data=df, columns=["Parameter", "Output", "Total-effect", "First-order",
        return df.astype({"Parameter": str, "Output": str, "Total-effect": float, "First-order": float,
                          "Total Error": float, "First Error": float})


class ModelResult:
    """
    Class that encapsulates the results of a model run.
    """
    pOC: [float] = None
    dOC: [float] = None
    IC: [float] = None
    cells: [float] = None
    t: [float] = None

    def get_dataframe(self, scenario, variable_title):
        """
        Returns model results as a dataframe. Removes duplicate values (e.g. discontinuities).
        """
        data = [self.t/365.25, self.pOC, self.dOC, self.IC, self.cells, [scenario]*len(self.t), [variable_title]*len(self.t)]
        df = pandas.DataFrame(data=np.column_stack(data),
                                columns=["Years from start", "POC", "DOC", "IC", "Cells", "Scenario", "Analysis type"])

        df.drop_duplicates(subset=["Years from start", "Scenario", "Analysis type"], keep="last", inplace=True)

        return df.astype({"Years from start": float, "POC": float, "DOC": float, "IC": float, "Cells": float,
                           "Scenario": str, "Analysis type": str})


class ExpansionResult:
    """
    Class that encapsulates the results of a brine expansion calculation.
    """
    ratio_volume: Decimal = None  # Ratio of expansion at 3D
    ratio_area: Decimal = None  # Ratio of expansion at 2D
    ratio_dimensions: Decimal = None  # Ratio of expansion at 1D
    expansion_a: Decimal = None  # Expansion along the A-axis of the prolate spheroid.
    expansion_b: Decimal = None  # Expansion along the B-axis of the prolate spheroid.


class EEAResult:
    """
    Class that encapsulates a EEA rate calculation result.
    """
    eea_upper: float = None
    eea_lower: float = None
    predicted_timespan: float = None
