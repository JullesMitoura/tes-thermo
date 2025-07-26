from thermo import Chemical
import pyomo.environ as pyo
import pandas as pd
import numpy as np
from thermo import Chemical, ChemicalConstantsPackage
from thermo.interaction_parameters import IPDB

def gibbs_pad(T,                        # temperatura em K
              component_list):          # lista com nomes dos componentes: ['methane', 'water', 'carbon monoxide', ...]
    T0 = 298.15
    results = []

    intercepto_metano = -73806.88626539786
    coeficientes_metano = np.array([ 9.29102160e+01, -3.63336791e-02])
    
    for comp_name in component_list:
        try:
            if comp_name.lower() == 'carbon':
                # y(T) = -0.00953 * T^2 + 1.16406 * T + 1174.78249
                carbon_adjusted_value = 1174.78249 + (1.16406 * T) + (-0.00953 * T**2)
                results.append(carbon_adjusted_value)
                continue

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
    
    return results

# Function to calculate gas fugacity

def fug(T,                          # Temperatura (K)
        P,                          # Pressão (bar)
        eq,                         # Nome da equação de estado ('Peng-Robinson', etc.)
        n,                          # Lista com o número de mols dos componentes
        components,                 # Lista de objetos com dados termodinâmicos (Tc, Pc, omega)
        kij):                       # Matriz ou DataFrame com parâmetros de interação binária (kij)

    R = 8.314462                    # Constante dos gases em J/(mol*K) ou Pa*m^3/(mol*K)
    P_pa = P * 1e5                  # Converte pressão de bar para Pa

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
    
    if eq not in eos_params:
        raise ValueError(f"Equação de estado '{eq}' não suportada.")
    
    params = eos_params[eq]
    n_array = np.array(n, dtype=float)
    y = n_array / np.sum(n_array)
    
    Tcs = np.array([c.Tc for c in components])
    Pcs = np.array([c.Pc for c in components])
    omegas = np.array([c.omega for c in components])
    
    kij_array = np.array(kij)

    m = params['m_func'](omegas)
    Tr = T / Tcs
    alpha = params['alpha_func'](Tr, m)
    
    a_i = params['Omega_a'] * (R**2 * Tcs**2 / Pcs) * alpha
    b_i = params['Omega_b'] * (R * Tcs / Pcs)
    
    a_ij = np.sqrt(np.outer(a_i, a_i)) * (1 - kij_array)
    
    a_mix = np.sum(np.outer(y, y) * a_ij)
    b_mix = np.sum(y * b_i)
    
    A = a_mix * P_pa / (R**2 * T**2)
    B = b_mix * P_pa / (R * T)
    
    coeffs = params['Z_coeffs'](A, B)
    Z_roots = np.roots(coeffs)
    
    real_roots = Z_roots[np.isreal(Z_roots)].real
    positive_real_roots = real_roots[real_roots > 0]

    if len(positive_real_roots) == 0:
        return np.full_like(y, np.nan)

    Z = positive_real_roots.max() 
    if Z <= B:
        return np.full_like(y, np.nan)

    term1 = b_i / b_mix * (Z - 1)
    term2 = -np.log(Z - B)
    sum_y_a_ij = np.dot(y, a_ij)
    
    term3_dyn = (2 * sum_y_a_ij / a_mix) - (b_i / b_mix)
    term3_log = params['ln_phi_term'](Z, B)
    ln_phi_i = term1 + term2 - (A / B) * term3_dyn * term3_log
    
    phi_i = np.exp(ln_phi_i)
        
    return phi_i

class Gibbs:
    def __init__(self, components, equation='Ideal Gas'):
        self.components = components
        self.equation = equation
        self.total_components = len(components)
        
        try:
            self.components_chemical = [Chemical(ID) for ID in self.components]
            self.constants, _ = ChemicalConstantsPackage.from_IDs(self.components)
            if self.equation == "Peng-Robinson":
                db_name = 'ChemSep PR'
            else:
                db_name = None

            if db_name:
                self.kijs = IPDB.get_ip_asymmetric_matrix(db_name, self.constants.CASs, 'kij')
            else:
                self.kijs = np.zeros((self.total_components, self.total_components))
        except Exception as e:
            print(f"Alerta: Falha ao carregar dados dos componentes ou kijs. {e}")
            self.components_chemical = []
            self.constants = None
            self.kijs = np.zeros((self.total_components, self.total_components))

        element_set = set()
        for chemical in self.components_chemical:
            if chemical.atoms:
                element_set.update(chemical.atoms.keys())
        
        self.species = sorted(list(element_set))
        self.total_species = len(self.species)
        self.A = np.array([
            [chemical.atoms.get(element, 0) for element in self.species]
            for chemical in self.components_chemical
        ])
        
        self.gas_indices = []
        self.solid_indices = []
        self.liquid_indices = []

        for i, component_name in enumerate(self.components):
            if component_name.lower() == 'carbon':
                self.solid_indices.append(i)
            else:
                self.gas_indices.append(i)

    def bnds_values(self, initial):
        max_species = np.dot(initial, self.A)
        epsilon = 1e-5
        bnds_aux = []
        for i in range(self.total_components):
            divisores = np.where(self.A[i] > 1e-12, self.A[i], np.inf)
            a = max_species / divisores
            positive_a = a[a > 1e-12]
            upper_bound = np.min(positive_a) if positive_a.size > 0 else epsilon
            bnds_aux.append((1e-8, max(upper_bound, epsilon)))
        return tuple(bnds_aux)

    def solve_gibbs(self, initial, T, P):
        initial = np.array(initial, dtype=float)
        initial[initial == 0] = 1e-10
        bnds = self.bnds_values(initial)

        model = pyo.ConcreteModel()
        model.n = pyo.Var(range(self.total_components), 
                          domain=pyo.NonNegativeReals, 
                          bounds=bnds)
        
        for i in range(self.total_components):
            model.n[i].value = max(bnds[i][0], min(initial[i], bnds[i][1]))

        def gibbs_rule(model):
            R = 8.314462
            P_bar = P  # P já está em bar
            
            n_values = np.array([pyo.value(model.n[i]) for i in range(self.total_components)])
            
            df_pad = gibbs_pad(T, self.components)
            
            if self.equation != 'Ideal Gas':
                phi_g_list = fug(T=T, 
                                 P=P_bar, 
                                 eq=self.equation, 
                                 n=n_values, 
                                 components=self.components_chemical,
                                 kij=self.kijs)
            else:
                phi_g_list = np.ones(self.total_components)
            
            total_gibbs = 0
            
            if self.gas_indices:
                gas_moles = [model.n[i] for i in self.gas_indices]
                total_gas_moles = sum(gas_moles)

                mu_gas = [
                    df_pad[i] + R * T * (
                        pyo.log(phi_g_list[i]) + 
                        pyo.log(P_bar) +
                        pyo.log(model.n[i] / total_gas_moles)
                    ) for i in self.gas_indices
                ]
                
                total_gibbs += sum(gas_moles[k] * mu_gas[k] for k in range(len(gas_moles)))
            
            if self.solid_indices:
                mu_solids = [df_pad[i] for i in self.solid_indices]
                solid_moles = [model.n[i] for i in self.solid_indices]
            
                total_gibbs += sum(solid_moles[k] * mu_solids[k] for k in range(len(solid_moles)))

            return total_gibbs

        model.obj = pyo.Objective(rule=gibbs_rule, sense=pyo.minimize)
        
        model.element_balance = pyo.ConstraintList()
        initial_species_moles = np.dot(initial, self.A)
        for i in range(self.total_species):
            lhs = sum(self.A[j, i] * model.n[j] for j in range(self.total_components))
            rhs = initial_species_moles[i]
            model.element_balance.add(lhs == rhs)

        solver = pyo.SolverFactory('ipopt')
        solver.options['tol'] = 1e-8
        solver.options['max_iter'] = 300
        
        results = solver.solve(model, tee=False)

        if results.solver.termination_condition in [pyo.TerminationCondition.optimal, pyo.TerminationCondition.locallyOptimal]:
            return [pyo.value(model.n[i]) for i in range(self.total_components)]
        else:
            print(f"Otimização falhou ou não convergiu. Condição: {results.solver.termination_condition}")
            return initial.tolist()
        
    def run(self,
                initial,
                Tmin,
                Tmax,
                Pmin,
                Pmax,
                nT,
                nP):
            
            TRange = np.linspace(Tmin, Tmax, nT)
            PRange = np.linspace(Pmin, Pmax, nP)

            all_results = []

            for T in TRange:
                for P in PRange:
                    equilibrium_moles = self.solve_gibbs(initial=initial,
                                                        T=T,
                                                        P=P)
                    
                    row_data = {'T': T, 'P': P}
                    component_data = dict(zip(self.components, equilibrium_moles))
                    
                    row_data.update(component_data)
                    all_results.append(row_data)

            return pd.DataFrame(all_results)