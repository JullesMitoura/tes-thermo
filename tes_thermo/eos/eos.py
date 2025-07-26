"""
This function supports ideal gas, virial, Peng-Robinson, and Soave-Redlich-Kwong equations 
of state (EoS) to determine the fugacity coefficients for components in the gas phase. 
It handles solid components by assigning them a fugacity coefficient of 1.

Args:
    T (float): System temperature in Kelvin (K).
    P (float): System pressure in bar.
    eq (str): The name of the equation of state to use for calculations.
        Supported options: 'Ideal Gas', 'Virial', 'Peng-Robinson', 'Soave-Redlich-Kwong'.
    n (list): A list of mole numbers for each component in the mixture.
        The order must match the order of components in the `components` list.
    components (list): A list of dictionaries containing the thermodynamic data for each
        component. Each dictionary must have keys like 'name', 'phase', 'Tc', 'Pc', 
        'omega', 'Vc', 'Zc'.
    kij (array-like, optional): Matrix with binary interaction parameters (kij).
        Only used for Peng-Robinson and Soave-Redlich-Kwong equations.

Returns:
    list: A list of fugacity coefficients (φ) for each component, in the same
    order as the input `components` list. Returns np.nan for components
    where calculation is not possible.
"""


import numpy as np
from tes_thermo.utils import setup_logger

logger = setup_logger()

def fug(T: float,
        P: float,
        eq: str,
        n: list,
        components: list,
        kij=None) -> list:

    R = 8.314462    # Universal gas constant in J/(mol*K) or Pa*m^3/(mol*K)
    P_pa = P

    num_components = len(components)
    total_n = sum(n)
    
    if total_n == 0:
        return [np.nan] * num_components
    
    if not components:
        return []

    # Create mole fractions list
    mole_fractions = [n_i / total_n for n_i in n]
    results_list = [0.0] * num_components

    # Separate components by phase
    gas_indices = []
    solid_indices = []
    
    for i, comp in enumerate(components):
        phase = comp.get('phase', 'g').lower()
        if phase == 's':
            solid_indices.append(i)
        else:
            gas_indices.append(i)

    # Assign fugacity coefficient of 1 to solid components
    for idx in solid_indices:
        results_list[idx] = 1.0
        
    if not gas_indices:
        return results_list

    if eq == 'Ideal Gas':
        for idx in gas_indices:
            results_list[idx] = 1.0
        return results_list

    elif eq == 'Virial':
        # Get gas components data
        gas_components = [components[i] for i in gas_indices]
        y = np.array([mole_fractions[i] for i in gas_indices])
        
        # Virial Equation (Truncated at the 2nd Coefficient)
        Tc = np.array([comp['Tc'] for comp in gas_components])
        omega = np.array([comp['omega'] for comp in gas_components])
        Zc = np.array([comp['Zc'] for comp in gas_components])
        Vc_input = np.array([comp['Vc'] for comp in gas_components])
        
        # Check if Vc is in cm^3/mol or m^3/mol and convert to m^3/mol
        if np.mean(Vc_input) > 1e-3:  # Likely in cm^3/mol
            Vc = Vc_input / 1e6  # Convert from cm^3/mol to m^3/mol
        else:  # Already in m^3/mol
            Vc = Vc_input

        num_gas_comps = len(gas_components)
        B_matrix = np.zeros((num_gas_comps, num_gas_comps))

        for i in range(num_gas_comps):
            for j in range(num_gas_comps):
                # Calculate kij using the formula: kij = 1 - 8*(Vc_i*Vc_j)^0.5 / (Vc_i^(1/3) + Vc_j^(1/3))^3
                if i == j:
                    kij_val = 0.0  # kii = 0 for pure components
                else:
                    numerator = 8 * (Vc[i] * Vc[j])**0.5
                    denominator = (Vc[i]**(1/3) + Vc[j]**(1/3))**3
                    kij_val = 1 - numerator / denominator
                
                Tcij = np.sqrt(Tc[i] * Tc[j]) * (1 - kij_val)
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

        # Assign fugacity coefficients to gas components
        for i, gas_idx in enumerate(gas_indices):
            results_list[gas_idx] = phi_k[i]

        return results_list

    elif eq in ['Peng-Robinson', 'Soave-Redlich-Kwong']:
        # Parameters for different equations of state
        eos_params = {
            'Peng-Robinson': {
                'Omega_a': 0.45724, 'Omega_b': 0.07780,
                'm_func': lambda w: 0.37464 + 1.54226 * w - 0.26992 * w**2,
                'alpha_func': lambda Tr, m: (1 + m * (1 - np.sqrt(Tr)))**2,
                'Z_coeffs': lambda A, B: [1, B - 1, A - 2*B - 3*B**2, -A*B + B**2 + B**3],
                'ln_phi_term': lambda Z, B: (1 / (2 * np.sqrt(2))) * np.log((Z + (1 + np.sqrt(2)) * B) / (Z + (1 - np.sqrt(2)) * B))
            },
            'Soave-Redlich-Kwong': {
                'Omega_a': 0.42748, 'Omega_b': 0.08664,
                'm_func': lambda w: 0.480 + 1.574 * w - 0.176 * w**2,
                'alpha_func': lambda Tr, m: (1 + m * (1 - np.sqrt(Tr)))**2,
                'Z_coeffs': lambda A, B: [1, -1, A - B - B**2, -A*B],
                'ln_phi_term': lambda Z, B: np.log(1 + B/Z)
            }
        }
        
        params = eos_params[eq]
        
        # Get gas components data
        gas_components = [components[i] for i in gas_indices]
        n_gas = np.array([n[i] for i in gas_indices])
        y = n_gas / np.sum(n_gas)  # Mole fractions for gas phase only
        
        # Extract thermodynamic properties for gas components
        Tcs = np.array([comp['Tc'] for comp in gas_components])
        Pcs = np.array([comp['Pc'] for comp in gas_components])
        omegas = np.array([comp['omega'] for comp in gas_components])
        
        # Set up kij matrix for gas components
        if kij is not None:
            # Extract kij for gas components only
            kij_gas = np.array([[kij[gas_indices[i], gas_indices[j]] for j in range(len(gas_indices))] 
                               for i in range(len(gas_indices))])
        else:
            kij_gas = np.zeros((len(gas_indices), len(gas_indices)))
        
        # Calculate EOS parameters
        m = params['m_func'](omegas)
        Tr = T / Tcs
        alpha = params['alpha_func'](Tr, m)
        
        a_i = params['Omega_a'] * (R**2 * Tcs**2 / Pcs) * alpha
        b_i = params['Omega_b'] * (R * Tcs / Pcs)
        
        # Mixing rules
        a_ij = np.sqrt(np.outer(a_i, a_i)) * (1 - kij_gas)
        
        a_mix = np.sum(np.outer(y, y) * a_ij)
        b_mix = np.sum(y * b_i)
        
        A = a_mix * P_pa / (R**2 * T**2)
        B = b_mix * P_pa / (R * T)
        
        # Solve cubic equation for compressibility factor
        coeffs = params['Z_coeffs'](A, B)
        Z_roots = np.roots(coeffs)
        
        real_roots = Z_roots[np.isreal(Z_roots)].real
        positive_real_roots = real_roots[real_roots > 0]
        
        if len(positive_real_roots) == 0:
            # If no valid roots, assign NaN to all gas components
            for gas_idx in gas_indices:
                results_list[gas_idx] = np.nan
            return results_list
            
        Z = positive_real_roots.max()  # Use largest positive real root (vapor phase)
        
        if Z <= B:
            # If Z is not physically meaningful
            for gas_idx in gas_indices:
                results_list[gas_idx] = np.nan
            return results_list
        
        # Calculate fugacity coefficients
        term1 = b_i / b_mix * (Z - 1)
        term2 = -np.log(Z - B)
        sum_y_a_ij = np.dot(y, a_ij)
        
        term3_dyn = (2 * sum_y_a_ij / a_mix) - (b_i / b_mix)
        term3_log = params['ln_phi_term'](Z, B)
        ln_phi_i = term1 + term2 - (A / B) * term3_dyn * term3_log
        
        phi_i = np.exp(ln_phi_i)
        
        # Assign fugacity coefficients to gas components
        for i, gas_idx in enumerate(gas_indices):
            results_list[gas_idx] = phi_i[i]

        return results_list
    
    else:
        raise ValueError(f"Equation of state '{eq}' is not supported. Available options: 'Ideal Gas', 'Virial', 'Peng-Robinson', 'Soave-Redlich-Kwong'.")