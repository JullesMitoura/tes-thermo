import numpy as np
import pandas as pd
import pyomo.environ as pyo
from thermo import Chemical, ChemicalConstantsPackage
from thermo.interaction_parameters import IPDB
from tes_thermo.utils import setup_logger, UnitConverter, get_solver, gibbs_pad
from tes_thermo.eos import fug

logger = setup_logger(__name__)

class Gibbs():
    def __init__(self,
                 components: dict = None,
                 equation: str = 'Ideal Gas',
                 inhibited_component: str = None,
                 solver_path: str = "tes/solver/bin/ipopt.exe"):
        
        self.component_objects = components
        self.thermo_components = components['components']
        self.components_chemical = [Chemical(ID) for ID in self.thermo_components]
        self.new_components = components['new_components']
        self.new_components_kij = self.new_components
        self.inhibited_component = inhibited_component
        self.equation = equation
        self.solver_path = solver_path
        self.components_info = self.define_components()
        self.component_names = [comp['name'] for comp in self.components_info]
        self.total_components = len(self.component_names)
        self.total_species = 0  # Will be set after element matrix creation

        logger.info(f"Initializing Gibbs class...")
        logger.info(f"Components created as Chemical from the thermo library: {self.thermo_components}")
        logger.info(f"Manually defined components: {self.new_components}")
        logger.info(f"Equation of state: {equation}")

        if self.equation == "Peng-Robinson":
            db_name = 'ChemSep PR'
            constants, _ = ChemicalConstantsPackage.from_IDs(self.thermo_components)
            try:
                self.kijs = IPDB.get_ip_asymmetric_matrix(db_name, constants.CASs, 'kij')
            except:
                logger.warning("Could not retrieve kijs from database. Using zeros.")
                self.kijs = np.zeros((len(self.thermo_components), len(self.thermo_components)))
        else:
            self.kijs = np.zeros((len(self.thermo_components), len(self.thermo_components)))

        self.kijs = np.array(self.kijs)
        logger.debug(f"Initial kijs matrix shape: {self.kijs.shape}")


        if self.new_components and len(self.new_components) > 0:
            for comp in self.new_components:
                if 'kijs' not in comp or comp['kijs'] is None:
                    comp['kijs'] = self._calculate_kijs_for_new_component(comp)
                    logger.info(f"Calculated kijs for {comp['name']}: {comp['kijs']}")
                else:
                    expected_size = self.kijs.shape[0] + 1
                    actual_size = len(comp['kijs'])
                    logger.debug(f"Component {comp['name']}: expected kijs size {expected_size}, actual size {actual_size}")
                    
                    if actual_size != expected_size:
                        logger.warning(f"kijs size mismatch for {comp['name']}. Expected {expected_size}, got {actual_size}. Recalculating...")
                        comp['kijs'] = self._calculate_kijs_for_new_component(comp)
                    else:
                        logger.info(f"Using provided kijs for {comp['name']}: {comp['kijs']}")
                
                self.kijs = self.add_component_kij(self.kijs, comp['kijs'])
                logger.debug(f"Updated kijs matrix shape after adding {comp['name']}: {self.kijs.shape}")

        element_set = set()
        for chemical in self.components_chemical:
            if chemical.atoms:
                element_set.update(chemical.atoms.keys())

        for comp in self.new_components:
            if 'structure' in comp:
                element_set.update(comp['structure'].keys())

        self.species = sorted(list(element_set))
        self.total_species = len(self.species)  # Set total_species here
        
        matrix_rows = []
        for chemical in self.components_chemical:
            row = [chemical.atoms.get(element, 0) for element in self.species]
            matrix_rows.append(row)

        for comp in self.new_components:
            structure = comp.get('structure', {})
            row = [structure.get(element, 0) for element in self.species]
            matrix_rows.append(row)

        self.A = np.array(matrix_rows)

    def _calculate_kijs_for_new_component(self, new_comp):
        """
        Calcula os parâmetros de interação kij para um novo componente usando a fórmula:
        kij = 1 - 8*(Vc_i*Vc_j)^0.5 / (Vc_i^(1/3) + Vc_j^(1/3))^3
        """
        # Extrair Vc do novo componente
        Vc_new = new_comp.get('Vc')
        if Vc_new is None:
            logger.warning(f"Vc not found for {new_comp['name']}. Setting all kijs to 0.")
            return [0.0] * (len(self.components_chemical) + 1)
        
        # Converter para m³/mol se necessário (assumindo que está em cm³/mol)
        if Vc_new > 1e-3:  # Likely in cm³/mol
            Vc_new = Vc_new / 1e6  # Convert to m³/mol
        
        kijs = []
        
        # Calcular kij com cada componente existente
        for chem in self.components_chemical:
            Vc_existing = getattr(chem, 'Vc', None)  # Usar getattr para evitar AttributeError
            if Vc_existing is None:
                kij = 0.0
                logger.warning(f"Vc not found for {chem.name}. Setting kij to 0.")
            else:
                # Converter para m³/mol se necessário
                if Vc_existing > 1e-3:  # Likely in cm³/mol
                    Vc_existing = Vc_existing / 1e6  # Convert to m³/mol
                
                try:
                    # Aplicar a fórmula: kij = 1 - 8*(Vc_i*Vc_j)^0.5 / (Vc_i^(1/3) + Vc_j^(1/3))^3
                    numerator = 8 * (Vc_new * Vc_existing)**0.5
                    denominator = (Vc_new**(1/3) + Vc_existing**(1/3))**3
                    
                    if denominator == 0:
                        kij = 0.0
                        logger.warning(f"Division by zero when calculating kij between {new_comp['name']} and {chem.name}. Setting kij to 0.")
                    else:
                        kij = 1 - numerator / denominator
                        
                except Exception as e:
                    kij = 0.0
                    logger.warning(f"Error calculating kij between {new_comp['name']} and {chem.name}: {e}. Setting kij to 0.")
            
            kijs.append(kij)
        
        # kii = 0 (componente consigo mesmo)
        kijs.append(0.0)
        
        return kijs

    def define_components(self):
        components_info = []
        for comp in self.components_chemical:
            props = {
                "name": comp.name,
                "Tc": comp.Tc,
                "Pc": comp.Pc,
                "omega": comp.omega,
                "Vc": comp.Vc,
                "Zc": comp.Zc,
                "Hfgm": comp.Hfgm,
                "Gfgm": comp.Gfgm,
                "phase": 's' if comp.phase == 's' else 'g'
            }
            components_info.append(props)

        for comp in self.new_components:
            props = {
                "name": comp.get('name', 'unknown'),
                "Tc": comp.get('Tc', None),
                "Pc": comp.get('Pc', None),
                "omega": comp.get('omega', None),
                "Vc": comp.get('Vc', None),
                "Zc": comp.get('Zc', None),
                "Hfgm": comp.get('Hfgm', None),
                "Gfgm": comp.get('Gfgm', None),
                "phase": comp.get('phase', 'g'),
            }
            components_info.append(props)
        return components_info

    @staticmethod
    def add_component_kij(kij_matrix, new_kij_row):
        kij_matrix = np.array(kij_matrix, dtype=float)
        new_kij_row = np.array(new_kij_row, dtype=float)

        if len(new_kij_row) != kij_matrix.shape[0] + 1:
            raise ValueError("Length of new_kij_row must be one more than the size of kij_matrix")

        new_matrix = np.zeros((kij_matrix.shape[0] + 1, kij_matrix.shape[1] + 1))
        new_matrix[:-1, :-1] = kij_matrix

        new_matrix[-1, :] = new_kij_row
        new_matrix[:, -1] = new_kij_row

        return new_matrix
    
    def identify_phases(self, phase_type):
        phases = [i for i, comp in enumerate(self.components_info) if comp.get("phase") == phase_type]
        logger.debug(f"Components identified for phase '{phase_type}': {phases}")
        return phases
    
    def _get_bounds(self, initial_moles: np.ndarray) -> tuple:
        """Calculates the upper and lower bounds for the mole number of each component.
            This allows the user to inhibit the formation of a component, so we define its maximum value as a very small value.
        """
        logger.debug("Calculating bounds for mole numbers of each component")
        
        max_species_moles = np.dot(initial_moles, self.A)
        epsilon = 1e-5
        bounds_list = []

        inhibited_idx = -1
        if self.inhibited_component and self.inhibited_component in self.component_names:
            inhibited_idx = self.component_names.index(self.inhibited_component)
            logger.info(f"Inhibited component found: {self.inhibited_component} (index: {inhibited_idx})")

        for i in range(len(self.component_names)):
            if i == inhibited_idx:
                bounds_list.append((1e-8, epsilon))
                logger.debug(f"Bounds for inhibited component {self.component_names[i]}: (1e-8, {epsilon})")
            else:
                with np.errstate(divide='ignore'):
                    a = np.multiply(1 / np.where(self.A[i] != 0, self.A[i], np.inf), max_species_moles)
                
                positive_limits = a[a > 0]
                upper_bound = np.min(positive_limits) if positive_limits.size > 0 else epsilon
                bounds_list.append((1e-8, max(upper_bound, epsilon)))
                logger.debug(f"Bounds for {self.component_names[i]}: (1e-8, {max(upper_bound, epsilon)})")

        logger.debug(f"Bounds calculated for all components: {len(bounds_list)} bounds")
        return tuple(bounds_list)
    
    def solve_gibbs(self, initial, T, P, T_unit, P_unit, progress_callback=None):

        T = UnitConverter.convert_temperature(T, T_unit)
        P = UnitConverter.convert_pressure(P, P_unit)
        logger.info(f"Starting Gibbs minimization problem resolution")
        logger.info(f"Temperature: {T} K, Pressure: {P} Pa.")
        logger.debug(f"Initial moles: {initial}")
        
        initial[initial == 0] = 0.00001
        logger.debug("Zero values in initial moles replaced with 0.00001")
        
        bnds = self._get_bounds(initial)
        solids = self.identify_phases('s')
        gases = self.identify_phases('g')
 
        logger.info(f"Phases identified - Solids: {len(solids)}, Gases: {len(gases)}")

        logger.debug("Creating Pyomo model")
        model = pyo.ConcreteModel()
        model.n = pyo.Var(range(len(self.component_names)), domain=pyo.NonNegativeReals, bounds=lambda m, i: bnds[i])
        
        def gibbs_rule(model):
            logger.debug("Calculating objective function (Gibbs energy)")
            R = 8.314  # J/mol·K
            
            df_pad = gibbs_pad(components=self.component_objects,T=T)
            logger.debug("Thermodynamic data obtained via gibbs_pad")
            print(df_pad)
            # Calculate fugacity coefficients using current mole values
            # For the initial calculation, we need to convert model.n to a list of values
            current_moles = [pyo.value(model.n[i]) if hasattr(model.n[i], 'value') and model.n[i].value is not None else initial[i] for i in range(self.total_components)]
            
            phii = fug(T=T, P=P, eq=self.equation, n=current_moles, components=self.components_info)
            logger.debug(f"Fugacity coefficients calculated using equation: {self.equation}")
            
            # Ensure phii is a list and convert any numpy types to Python floats
            if not isinstance(phii, list):
                phii = [float(phii)] * self.total_components
            else:
                phii = [float(phi) for phi in phii]
            
            logger.debug(f"Fugacity coefficients: {phii}")

            # df_pad now returns a list of mu_i values in the same order as components
            # Create a dictionary mapping component names to mu_i values
            mu_dict = {self.component_names[i]: df_pad[i] for i in range(len(self.component_names))}
            
            # Calculate chemical potentials for gas phase components
            mi_gas = []
            for i in gases:
                mu_i = mu_dict[self.component_names[i]]
                phi_i = phii[i]
                
                gas_moles = [model.n[i] for i in gases]
                total_gas_moles = sum(gas_moles)
                
                epsilon = 1e-18
                # Build the chemical potential expression
                chemical_potential = (
                    mu_i + 
                    R * T * (
                        pyo.log(phi_i) + 
                        pyo.log(model.n[i]/total_gas_moles + epsilon) + 
                        pyo.log(P*1e-5)
                    )
                )
                mi_gas.append(chemical_potential)

            # Chemical potentials for solid phase components
            mi_solids = [mu_dict[self.component_names[i]] for i in solids]

            # Calculate total Gibbs energy
            regularization_term = 1e-6
            total_gibbs = (
                sum(mi_gas[idx] * model.n[gases[idx]] for idx in range(len(gases))) + 
                sum(mi_solids[idx] * model.n[solids[idx]] for idx in range(len(solids))) + 
                regularization_term
            )
            
            logger.debug("Objective function calculated successfully")
            return total_gibbs
        
        model.obj = pyo.Objective(rule=gibbs_rule, sense=pyo.minimize)
        logger.debug("Objective function defined in model")
        
        model.element_balance = pyo.ConstraintList()
        for i in range(self.total_species):
            tolerance = 1e-8
            lhs = sum(self.A[j, i] * model.n[j] for j in range(self.total_components))
            rhs = sum(self.A[j, i] * initial[j] for j in range(self.total_components))
            model.element_balance.add(pyo.inequality(-tolerance, lhs - rhs, tolerance))

        logger.info(f"Element balance constraints added for {self.total_species} species")

        logger.debug("Getting solver")
        solver = get_solver(self.solver_path)

        solver.options['tol'] = 1e-8
        solver.options['max_iter'] = 5000
        logger.debug("Solver options configured: tol=1e-8, max_iter=5000")

        logger.info("Starting optimization problem resolution")
        results = solver.solve(model, tee=False)

        if results.solver.termination_condition == pyo.TerminationCondition.optimal:
            logger.info("Optimal solution found successfully")
            
            solution = {
                    "Temperature (K)": T,
                    "Pressure (bar)": P*1e-5
                }
            solution.update({
                            name.capitalize().replace("_", " "): pyo.value(model.n[i])
                            for i, name in enumerate(self.component_names)
                        })
            
            logger.debug("Solution compiled:")
            for key, value in solution.items():
                if key not in ["Temperature (K)", "Pressure (bar)"]:
                    logger.debug(f"  {key}: {value:.6e} mol")
            
            return solution
        else:
            error_msg = f"Optimal solution not found. Termination condition: {results.solver.termination_condition}"
            logger.error(error_msg)
            raise Exception("Optimal solution not found.")
        
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
                    component_data = dict(zip(self.component_names, equilibrium_moles))
                    
                    row_data.update(component_data)
                    all_results.append(row_data)

            return pd.DataFrame(all_results)