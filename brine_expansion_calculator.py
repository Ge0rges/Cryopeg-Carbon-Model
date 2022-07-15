from decimal import Decimal, getcontext, ROUND_CEILING

getcontext().rounding = ROUND_CEILING

## BRINE VOLUME PARAMATERS
# Prolate spheroid dimensions
a = Decimal('2.185099')  # cm - width
b = Decimal('4')  # cm - length
brine_volume = Decimal('4') / Decimal('3') * Decimal('3.14159') * (a / Decimal('2')) ** Decimal('2') * (b / Decimal('2'))  # Volume of prolate spheroid: V = 4/3 * pi * (a/2) ** 2 * (b/2)


def calculate_brine_expansion(carbon_density_in_permafrost, carbon_required_per_year):
    # Calculate volume needed to expand
    volume_needed = carbon_required_per_year / carbon_density_in_permafrost  # The volume that must be covered to get the required amount of C

    # Expansion ratio of volume
    ratio = volume_needed / brine_volume

    # Expansion ratio of area
    ratio_area = (1 + volume_needed / brine_volume) ** (Decimal('2/3'))

    # Expansion ratio linearly
    ratio_dimensions = (1 + volume_needed / brine_volume) ** (Decimal('1/3'))

    expansion_a = a * ratio_dimensions
    expansion_b = b * ratio_dimensions

    return ratio, ratio_area, ratio_dimensions, expansion_a, expansion_b

    # print("% expansion of volume based on brine volume per year is: " + "{:.0f}".format(ratio * 100) + "%")
    #
    # print("% expansion of area based on brine volume per year is: " + "{:.0f}".format((ratio_area - 1) * 100) + "%")
    #
    # print("Brine will expand by (a, b) cm per year: (" + "{:.2f}".format
    # (expansion_a) + ", " + "{:.2f}".format(expansion_b) + ") = " + "{:.0f}".format((ratio_dimensions - 1) * 100) + "%")
