"""
This function determines the chemical potential (μ) by integrating the Gibbs-Helmholtz
equation, using provided thermochemical data (ΔH, ΔG) and a user-defined
polynomial for the heat capacity (Cp). It handles multiple components in a single
call by looking up their respective heat capacity coefficients by name.

Args:
    T (float): The target temperature for the calculation, in Kelvin (K).
    components (List[Component]): A list of 'Component' objects, where each
        object contains properties like name, 'deltaHf' (standard enthalpy
        of formation), and 'deltaGf' (standard Gibbs free energy of formation).
    cp_polynomial (CpFactory): A factory function that accepts heat capacity
        coefficients (e.g., a, b, c, d) and returns a callable function,
        `cp(T)`, which calculates the heat capacity at a given temperature.
    cp_coefficients (Dict[str, Dict[str, float]]): A dictionary where keys are
        component names (matching those in the `components` list) and values
        are dictionaries containing the coefficients for the `cp_polynomial`
        factory.

Returns:
    List[dict]: A list of dictionaries, where each dictionary contains the
        name of the component and its calculated chemical potential ('mu_i')
        in Joules per mole (J/mol).
"""

from typing import Callable, List, Dict
from scipy.integrate import quad
from tes.utils import Component

CpFunction = Callable[[float], float]
CpFactory = Callable[..., CpFunction]

def gibbs_pad(T: float,                                                              # Default value in Kelvin but the user can indicate any temperature, it will be converted.
              components: List[Component],                                           # List of components.
              cp_polynomial: CpFactory,                                              # Cp polynomial.
              cp_coefficients: Dict[str, Dict[str, float]]) -> List[float]:          # List of coefficients for the Cp polynomial.
    
    from tes.utils import setup_logger
    logger = setup_logger()
    T0 = 298.15  # Reference temperature in Kelvin
    
    results = []

    for comp in components:
        props = comp.get_properties()

        deltaH = props['deltaHf']
        deltaG = props['deltaGf']
        name = props['name']
        coeffs = cp_coefficients.get(name)

        if not coeffs:
            logger.warning(f"Component '{name}' not found in cp_coefficients. Skipping.")
            continue

        cp = cp_polynomial(**coeffs)

        def inner_integral(T_prime):
            h_integral, _ = quad(cp, T0, T_prime)
            return (deltaH + h_integral) / T_prime ** 2
        integral_value, _ = quad(inner_integral, T0, T)
        mu_i = T * ((deltaG / T0) - integral_value)
        results.append({'name': name, 'mu_i': mu_i})

    return results