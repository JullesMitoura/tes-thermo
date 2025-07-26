from typing import Callable, List, Dict, Any
import numpy as np
from scipy.integrate import quad
from thermo import Chemical

CpFunction = Callable[[float], float]
CpFactory = Callable[..., CpFunction]

def gibbs_pad(T: float, # Temperature in Kelvin
              components: Dict[str, Any]) -> List[float]: # Dictionary with 'components' and 'new_components'
    """
    Args:
        T: Temperature in Kelvin
        components: Dict with 'components' (list of component names) and 
                   'new_components' (list of component dictionaries)
    
    Returns:
        List of Gibbs free energy values
    """
    from tes_thermo.utils import setup_logger
    logger = setup_logger()
    
    T0 = 298.15  # Reference temperature in Kelvin
    results = []
    
    # Extract components from the input dictionary
    thermo_components = components.get('components', [])
    new_components = components.get('new_components', [])
    
    # Process existing components from thermo library
    for comp_name in thermo_components:
            # Special case for methane with adjustment
            intercepto_metano = -73806.88626539786
            coeficientes_metano = np.array([ 9.29102160e+01, -3.63336791e-02])
            
            try:
                if comp_name.lower() == 'carbon':
                    # y(T) = -0.00953 * T^2 + 1.16406 * T + 1174.78249
                    carbon_adjusted_value = 1174.78249 + (1.16406 * T) + (-0.00953 * T**2)
                    results.append(carbon_adjusted_value)
                    continue
                
                # Get component properties at reference temperature
                comp = Chemical(comp_name, T=T0)
                deltaH_T0, deltaG_T0 = comp.Hfgm, comp.Gfgm
                
                if deltaH_T0 is None or deltaG_T0 is None:
                    results.append(0.0)
                    continue
                
                T_avg = (T0 + T) / 2
                comp_avg = Chemical(comp_name, T=T_avg)
                Cp = comp_avg.Cpm if comp_avg.Cpm is not None else 0
                
                deltaH_T = deltaH_T0 + Cp * (T - T0)
                deltaS_T0 = (deltaH_T0 - deltaG_T0) / T0
                
                if T > 0 and T0 > 0:
                    deltaS_T = deltaS_T0 + Cp * np.log(T / T0)
                else:
                    deltaS_T = deltaS_T0
                
                mu_standard = deltaH_T - T * deltaS_T

                if comp_name.lower() == 'methane':
                    B0 = intercepto_metano
                    B1 = coeficientes_metano[0]
                    B2 = coeficientes_metano[1]
                    
                    adjusted_value = B0 + (B1 * T) + (B2 * T**2)
                    results.append(adjusted_value)
                
                else:
                    results.append(mu_standard)
                
            except Exception as e:
                print(f"Erro ao processar {comp_name}: {e}")
                results.append(0.0)
    
    # Process new components
    for comp_dict in new_components:
        try:
            # Extract properties from the component dictionary
            deltaH = comp_dict.get('Hfgm', 0)  # Formation enthalpy
            deltaG = comp_dict.get('Gfgm', 0)  # Formation Gibbs energy
            name = comp_dict.get('name', 'unknown')
            cp_function = comp_dict.get('cp_polynomial')

            # Use the cp_polynomial function for more accurate calculation
            def cp_value(T_prime):
                """Integrand for enthalpy calculation: Cp(T')"""
                return cp_function(T_prime)
            
            def inner_integral(T_prime):
                value, _ = quad(cp_value, T0, T_prime)
                return (deltaH + value) / T_prime ** 2

            integral_value, _ = quad(inner_integral, T0, T)
            mu_i = T * (deltaG / T0 - integral_value)
            results.append(mu_i)

        except Exception as e:
            logger.error(f"Error processing new component {comp_dict.get('name', 'unknown')}: {e}")
            results.append(0.0)
    
    return results