"""
    Calculates chemical equilibrium by minimizing the total Gibbs energy of a system.

    This class sets up and solves a constrained optimization problem using Pyomo.
    It processes component data, including thermodynamic properties and atomic
    structure, to enforce elemental balance constraints.
"""


import pyomo.environ as pyo
import pandas as pd
import numpy as np
from tes.utils import gibbs_pad, get_solver, setup_logger, UnitConverter
from tes.eos import fug

logger = setup_logger()

class Gibbs():
    def __init__(self,
                 components: list,
                 cp_coefficients: dict,
                 cp_polynomial_factory: callable,
                 kij: pd.DataFrame = None,
                 equation: str = 'Ideal Gas',
                 inhibited_component: str = None,
                 solver_path: str = "tes/solver/bin/ipopt.exe"):
        """
        Initializes the Gibbs energy minimizer.

        Args:
            components (list): A list of pre-configured Component objects.
            cp_coefficients (dict): Dictionary with heat capacity coefficients.
            cp_polynomial_factory (callable): A factory function that returns a Cp(T) function.
            kij (pd.DataFrame): DataFrame with binary interaction parameters.
            equation (str): Equation of state to use ('Ideal Gas', 'Peng-Robinson', etc.).
            inhibited_component (str, optional): Name of a component whose formation
                                                 is to be suppressed. Defaults to None.
            solver_path (str, optional): Absolute path to the solver executable (e.g., ipopt).
                                         Defaults to None, assuming it's in the system's PATH.
        """

        logger.info("Initializing Gibbs class")
        logger.debug(f"Equation of state: {equation}")
        logger.debug(f"Inhibited component: {inhibited_component}")
        logger.debug(f"Solver path: {solver_path}")

        self.component_objects = components
        self.cp_coefficients = cp_coefficients
        self.cp_polynomial_factory = cp_polynomial_factory
        self.inhibited_component = inhibited_component
        self.equation = equation
        self.kij = kij
        self.solver_path = solver_path
        

        self.components_data = {comp.get_properties()['name']: comp.get_properties() for comp in self.component_objects}
        self.component_names = list(self.components_data.keys())
        self.total_components = len(self.component_names)

        logger.info(f"Total number of components: {self.total_components}")
        logger.debug(f"Components: {self.component_names}")

        species_set = set()
        for props in self.components_data.values():
            species_set.update(props.get('structure', {}).keys())
        self.species = sorted(list(species_set))
        self.total_species = len(self.species)
        
        logger.info(f"Total number of species: {self.total_species}")
        logger.debug(f"Species: {self.species}")

        self.A = np.array([
            [self.components_data[name].get('structure', {}).get(spec, 0) for spec in self.species]
            for name in self.component_names
        ])

        logger.debug(f"Stoichiometric matrix A created with shape: {self.A.shape}")

    def identify_phases(self, phase_type):
        """
        Identifies the components that belong to the given phase type ('s' for solids, 'g' for gases).
        """
        phases = [i for i, comp in enumerate(self.components_data) if self.components_data[comp].get("phase") == phase_type]
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

        for i in range(self.total_components):
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
        logger.info(f"Temperature: {T} K, Pressure: {P} bar")
        logger.debug(f"Initial moles: {initial}")
        
        initial[initial == 0] = 0.00001
        logger.debug("Zero values in initial moles replaced with 0.00001")
        
        bnds = self._get_bounds(initial)
        solids = self.identify_phases('s')
        gases = self.identify_phases('g')

        logger.info(f"Phases identified - Solids: {len(solids)}, Gases: {len(gases)}")

        logger.debug("Creating Pyomo model")
        model = pyo.ConcreteModel()
        model.n = pyo.Var(range(self.total_components), domain=pyo.NonNegativeReals, bounds=lambda m, i: bnds[i])
        
        def gibbs_rule(model):
            logger.debug("Calculating objective function (Gibbs energy)")
            R = 8.314  # J/mol·K
            
            df_pad = gibbs_pad(T = T,
                               components = self.component_objects,
                               cp_polynomial = self.cp_polynomial_factory,
                               cp_coefficients = self.cp_coefficients)
            
            logger.debug("Thermodynamic data obtained via gibbs_pad")
            
            phii = fug(T=T, P=P, eq=self.equation, n=model.n, components=self.components_data)
            logger.debug(f"Fugacity coefficients calculated using equation: {self.equation}")

            if isinstance(phii, (int, float)):  
                phii = [phii] * self.total_components
                logger.debug("Fugacity coefficient converted to list")

            mu_dict = {item['name']: item['mu_i'] for item in df_pad}
            
            mi_gas = [
                    mu_dict[self.component_names[i]] + R * T * (
                        pyo.log(phii[i]) + 
                        pyo.log(model.n[i] / sum(model.n[j] for j in range(self.total_components))) + 
                        pyo.log(P)
                    ) for i in gases
                ]

            mi_solids = [
                            mu_dict[self.component_names[i]] for i in solids
                        ]

            regularization_term = 1e-6
            total_gibbs = sum(mi_gas[i] * model.n[gases[i]] for i in range(len(mi_gas))) + \
                        sum(mi_solids[i] * model.n[solids[i]] for i in range(len(mi_solids))) + \
                        regularization_term
            
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
                    "Pressure (bar)": P
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