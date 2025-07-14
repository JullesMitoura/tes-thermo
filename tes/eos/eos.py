"""
This function supports ideal gas and virial equations of state (EoS) to determine 
the fugacity coefficients for components in the gas phase. It handles solid 
components by assigning them a fugacity coefficient of 1.

Args:
    T (float): System temperature in Kelvin (K).
    P (float): System pressure in bar.
    eq (str): The name of the equation of state to use for calculations.
        Supported options: 'Ideal Gas', 'Virial'.
    n (list): A list of mole numbers for each component in the mixture.
        The order must match the order of components in the `components` dictionary.
    components (dict): A dictionary containing the thermodynamic data for each
        component. The keys are component names, and the values are dictionaries
        with properties like 'Phase', 'Tc', 'Pc', 'omega', 'Vc', 'Zc'.

Returns:
    list: A list of fugacity coefficients (φ) for each component, in the same
    order as the input `components` dictionary. Returns np.nan for components
    where calculation is not possible.
"""


import numpy as np
from tes.utils import setup_logger

logger = setup_logger()

def fug(T: float,
        P: float,
        eq: str,
        n: list,
        components: dict) -> list:

    R = 8.314462    # Universal gas constant in J/(mol*K) or Pa*m^3/(mol*K)
    P_pa = P * 1e5  # Convert pressure from bar to Pa

    comp_names = list(components.keys())
    total_n = sum(n)
    
    if total_n == 0:
        return [np.nan] * len(comp_names)
    
    if not comp_names:
        return []

    mole_fractions = {name: n_i / total_n for name, n_i in zip(comp_names, n)}
    results_list = [0.0] * len(comp_names)

    # Separate components by phase
    gas_components = {name: data for name, data in components.items() if data.get('phase', 'g').lower() != 's'}
    solid_components = {name: data for name, data in components.items() if data.get('phase', 'g').lower() == 's'}

    # Assign fugacity coefficient of 1 to solid components
    for name in solid_components:
        idx = comp_names.index(name)
        results_list[idx] = 1.0
        
    if not gas_components:
        return results_list

    if eq == 'Ideal Gas':
        for name in gas_components:
            idx = comp_names.index(name)
            results_list[idx] = 1.0
        return results_list

    if eq == 'Virial':
        gas_comp_names = list(gas_components.keys())
        y = np.array([mole_fractions[name] for name in gas_comp_names])
        
        # Virial Equation (Truncated at the 2nd Coefficient)
        Tc = np.array([gas_components[name]['Tc'] for name in gas_comp_names])
        omega = np.array([gas_components[name]['omega'] for name in gas_comp_names])
        Zc = np.array([gas_components[name]['Zc'] for name in gas_comp_names])
        Vc_cm3_mol = np.array([gas_components[name]['Vc'] for name in gas_comp_names])
        
        # Convert Vc from cm^3/mol to m^3/mol
        Vc = Vc_cm3_mol / 1e6

        num_comps = len(gas_comp_names)
        B_matrix = np.zeros((num_comps, num_comps))

        for i in range(num_comps):
            for j in range(num_comps):
                # Calculate kij using the formula: kij = 1 - 8*(Vc_i*Vc_j)^0.5 / (Vc_i^(1/3) + Vc_j^(1/3))^3
                if i == j:
                    kij = 0.0  # kii = 0 for pure components
                else:
                    numerator = 8 * (Vc[i] * Vc[j])**0.5
                    denominator = (Vc[i]**(1/3) + Vc[j]**(1/3))**3
                    kij = 1 - numerator / denominator
                
                Tcij = np.sqrt(Tc[i] * Tc[j]) * (1 - kij)
                wij = (omega[i] + omega[j]) / 2
                Vcij = ((Vc[i]**(1/3) + Vc[j]**(1/3)) / 2)**3
                Zcij = (Zc[i] + Zc[j]) / 2
                Pcij_pa = Zcij * R * Tcij / Vcij
                
                Tr_ij = T / Tcij
                B0 = 0.083 - 0.422 / (Tr_ij**1.6)
                B1 = 0.139 - 0.172 / (Tr_ij**4.2)
                B_matrix[i, j] = (R * Tcij / Pcij_pa) * (B0 + wij * B1)
        
        B_mix = y.T @ B_matrix @ y
        sum_yB = B_matrix @ y
        ln_phi_k = (2 * sum_yB - B_mix) * P_pa / (R * T)
        
        phi_k = np.exp(ln_phi_k)

        for i, name in enumerate(gas_comp_names):
            idx = comp_names.index(name)
            results_list[idx] = phi_k[i]

        return results_list
    
    else:
        raise ValueError(f"Equation of state '{eq}' is not supported. Only 'Ideal Gas' and 'Virial' are available.")