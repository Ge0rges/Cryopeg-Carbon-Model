from decimal import Decimal, getcontext, ROUND_CEILING

getcontext().rounding = ROUND_CEILING

## BRINE VOLUME PARAMATERS
# Prolate spheroid dimensions
a = Decimal('2.185099')  # cm - width
b = Decimal('4')  # cm - length
brine_volume = Decimal('4') / Decimal('3') * Decimal('3.14159') * (a / Decimal('2')) ** Decimal('2') * (b / Decimal('2'))  # Volume of prolate spheroid: V = 4/3 * pi * (a/2) ** 2 * (b/2)


def calculate_brine_expansion():
    # Lower and upper bounds for how much migration we need. From Go excel file.
    carbon_density_lower = Decimal(51.23) / Decimal('1e6')  # g/ml - Carbon surrounding the brine (lower bound)
    carbon_density_upper = Decimal(1286.81) / Decimal('1e6')  # g/ml - Carbon surrounding the brine (upper bound)

    # Calculate system organic carbon requirements
    total_carbon_required = (I[-1] - initial_carbon_content - organic_carbon_content_per_cell * initial_cell_count) / one_gram_in_femto * brine_volume  # grams
    total_carbon_required_per_year = total_carbon_required / (duration / Decimal('365.25'))  # (grams) /year

    # Calculate volume needed to expand
    volume_needed_lower = total_carbon_required_per_year / carbon_density_lower  # The volume that must be covered to get the required amount of C
    volume_needed_upper = total_carbon_required_per_year / carbon_density_upper  # The volume that must be covered to get the required amount of C

    # Expansion ratio of volime
    ratio_lower = volume_needed_lower / brine_volume
    ratio_upper = volume_needed_upper / brine_volume

    # Expansion ratio of area
    ratio_area_lower = (1 + volume_needed_lower / brine_volume) ** (Decimal('2/3'))
    ratio_area_upper = (1 + volume_needed_upper / brine_volume) ** (Decimal('2/3'))

    # Expansion ratio linearly
    ratio_dimensions_lower = (1 + volume_needed_lower / brine_volume) ** (Decimal('1/3'))
    ratio_dimensions_upper = (1 + volume_needed_upper / brine_volume) ** (Decimal('1/3'))

    a_low = a * ratio_dimensions_lower
    b_low = b * ratio_dimensions_lower

    a_high = a * ratio_dimensions_upper
    b_high = b * ratio_dimensions_upper

    print("Total carbon required over timespan is: " + "{:.3f}".format(total_carbon_required) + "g")
    print("Total carbon required on average per year is: " + "{:.3f}".format(total_carbon_required_per_year) + "g")

    print("Lower bound for % expansion of volume based on brine volume per year is: " + "{:.0f}".format(
        ratio_lower * 100) + "%")
    print("Upper bound for % expansion of volume based on brine volume per year is: " + "{:.0f}".format(
        ratio_upper * 100) + "%")

    print("Lower bound for % expansion of area based on brine volume per year is: " + "{:.0f}".format(
        (ratio_area_lower - 1) * 100) + "%")
    print("Upper bound for % expansion of area based on brine volume per year is: " + "{:.0f}".format(
        (ratio_area_upper - 1) * 100) + "%")

    print("Lower bound of carbon, brine will expand by (a, b) cm per year: (" + "{:.2f}".format(
        a_low) + ", " + "{:.2f}".format(b_low) + ") = " + "{:.0f}".format((ratio_dimensions_lower - 1) * 100) + "%")
    print("Upper bound of carbon, brine will expand by (a, b) cm per year: (" + "{:.2f}".format(
        a_high) + ", " + "{:.2f}".format(b_high) + ") = " + "{:.0f}".format((ratio_dimensions_upper - 1) * 100) + "%")
