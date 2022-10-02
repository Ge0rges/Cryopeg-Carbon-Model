"""
Some helper functions. Includes classes to encapsulate analysis results.
"""

from decimal import Decimal


def convert_micromolar_carbon_to_fgC_per_ml(micromolar):
    """Converts a micromolar of carbon value to femtograms of carbon per ml."""
    return (micromolar * 12.011 * 10**-6 * 10**15)/1000


def convert_fgC_per_ml_to_micromolar_carbon(fgC_per_ml):
    """Converts a femtograms of carbon per ml value to micromolar of carbon."""
    return (fgC_per_ml * 1000)/(12.011 * 10**-6 * 10**15)


class MaintenanceResult:
    """
    Class to encapsulate the maintenance energy results.
    """
    lower_bound_me: float = None
    upper_bound_me: float = None
    minimum_growth_rate: float = None
    minimum_doubling_time: float = None


class SensitivityResult:
    """
    Class to encapsulate results of a sensitivity analysis.
    """
    total_sobol_indices: [float] = None
    first_order_sobol_indices: [float] = None


class ModelResult:
    """
    Class to encapsulate results of a model run.
    """
    pOC: [float] = None
    dOC: [float] = None
    IC: [float] = None
    cells: [float] = None
    t: [float] = None


class ExpansionResult:
    """
    Class to encapsulate results of a brine expansion calculation.
    """
    ratio_volume: Decimal = None  # Ratio of expansion at 3D
    ratio_area: Decimal = None  # Ratio of expansion at 2D
    ratio_dimensions: Decimal = None  # Ratio of expansion at 1D
    expansion_a: Decimal = None  # Expansion along the A-axis of the prolate spheroid.
    expansion_b: Decimal = None  # Expansion along the B-axis of the prolate spheroid.
