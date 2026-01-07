"""
Thermodynamic simulation function for OpenAI function calling.
Replaces LangChain BaseTool with native OpenAI function format.
"""
from typing import Dict, Any, Optional
import numpy as np
import pandas as pd
from tes_thermo.utils.prompts import Prompts
from tes_thermo.utils.units import convert_pressure_to_bar, convert_temperature_to_K
from tes_thermo.gibbs import Gibbs
from tes_thermo.utils import Component


def min_g_function(Tmin: float = 600.0,
                   Tmax: float = 1200.0,
                   Tunit: str = "K",
                   Pmin: float = 1.0,
                   Pmax: float = 10.0,
                   Punit: str = "bar",
                   Equation: str = "Peng-Robinson",
                   SelectedComponents: Dict[str, float] = None) -> str:
    """
    Simulate an isothermal reactor using Gibbs energy minimization.
    
    Args:
        Tmin: Minimum temperature value
        Tmax: Maximum temperature value
        Tunit: Unit of measurement for temperature (K, F, C)
        Pmin: Minimum pressure value
        Pmax: Maximum pressure value
        Punit: Unit of measurement for pressure (bar, Pa, MPa)
        Equation: Equation of state (Peng-Robinson, Ideal Gas)
        SelectedComponents: Dictionary with component formulas as keys and quantities as values
    
    Returns:
        JSON string with simulation results as DataFrame
    """
    if SelectedComponents is None or len(SelectedComponents) == 0:
        return "Error: SelectedComponents must be provided with at least one component."
    
    try:
        # Extract components and compositions in order (Python 3.7+ maintains dict order)
        components = list(SelectedComponents.keys())
        compositions = list(SelectedComponents.values())
        
        # Convert to numpy array and ensure proper types (fix for composition bug)
        compositions = np.array(compositions, dtype=float)
        
        # Create component structure using Component class
        comp_obj = Component(components=components, new_component={})
        components_dict = comp_obj.get_components()
        
        # Initialize Gibbs with the correct structure
        gibbs = Gibbs(components=components_dict, equation=Equation)
        
        Tmin_K = convert_temperature_to_K(Tmin, Tunit)
        Tmax_K = convert_temperature_to_K(Tmax, Tunit)
        Pmin_bar = convert_pressure_to_bar(Pmin, Punit)
        Pmax_bar = convert_pressure_to_bar(Pmax, Punit)
        
        if Tmin_K != Tmax_K:
            TRange = np.linspace(Tmin_K, Tmax_K, 10)
        else:
            TRange = np.array([Tmin_K])
        if Pmin_bar != Pmax_bar:
            PRange = np.linspace(Pmin_bar, Pmax_bar, 10)
        else:
            PRange = np.array([Pmin_bar])
        
        all_results = []
        for T in TRange:
            for P in PRange:
                # solve_gibbs now requires T_unit and P_unit, and returns a dict
                solution = gibbs.solve_gibbs(
                    initial=compositions,
                    T=T,
                    P=P,
                    T_unit='K',
                    P_unit='bar'
                )
                
                # Convert solution dict to row_data format
                row_data = {'T': T, 'P': P}
                component_data = {}
                
                # Use gibbs.component_names to get the correct order and format names
                if hasattr(gibbs, 'component_names'):
                    # Map each input component name to its converted and formatted name in solution
                    for i, comp_name in enumerate(components):
                        if i < len(gibbs.component_names):
                            # Get the component name from Gibbs (converted name, e.g., "methane" not "CH4")
                            gibbs_comp_name = gibbs.component_names[i]
                            # Format it the same way Gibbs does in solve_gibbs
                            formatted_name = gibbs_comp_name.capitalize().replace("_", " ")
                            # Get the value from solution dict using the formatted name
                            if formatted_name in solution:
                                # Map the original input name (e.g., "CH4") to the value from solution
                                component_data[comp_name] = float(solution[formatted_name])
                            else:
                                component_data[comp_name] = 0.0
                        else:
                            component_data[comp_name] = 0.0
                else:
                    # Fallback if component_names not available
                    for comp_name in components:
                        component_data[comp_name] = 0.0
                
                row_data.update(component_data)
                all_results.append(row_data)
        
        df = pd.DataFrame(all_results)
        # Convert DataFrame to JSON string for OpenAI response
        return df.to_json(orient='records', indent=2)
        
    except Exception as e:
        return f"Error during simulation: {str(e)}"


def get_min_g_function_schema() -> Dict[str, Any]:
    """
    Returns OpenAI function calling schema for min_g_function.
    """
    return {
        "type": "function",
        "function": {
            "name": "min_g_calc",
            "description": Prompts.ming(),
            "parameters": {
                "type": "object",
                "properties": {
                    "Tmin": {
                        "type": "number",
                        "description": "Minimum temperature value (e.g., 600). Default: 600.0"
                    },
                    "Tmax": {
                        "type": "number",
                        "description": "Maximum temperature value (e.g., 1200). Default: 1200.0"
                    },
                    "Tunit": {
                        "type": "string",
                        "description": "Unit of measurement for temperature (e.g., K, F, C). Default: K",
                        "enum": ["K", "F", "C"]
                    },
                    "Pmin": {
                        "type": "number",
                        "description": "Minimum pressure value (e.g., 1). Default: 1.0"
                    },
                    "Pmax": {
                        "type": "number",
                        "description": "Maximum pressure value (e.g., 10). Default: 10.0"
                    },
                    "Punit": {
                        "type": "string",
                        "description": "Unit of measurement for pressure (e.g., bar, Pa, MPa). Default: bar",
                        "enum": ["bar", "Pa", "MPa"]
                    },
                    "Equation": {
                        "type": "string",
                        "description": "Equation of state. Default: Peng-Robinson",
                        "enum": ["Peng-Robinson", "Ideal Gas"]
                    },
                    "SelectedComponents": {
                        "type": "object",
                        "description": "REQUIRED: Dictionary with ALL chemical components involved in the reaction. Keys must be exact chemical formulas (e.g., 'CH4', 'H2O', 'CO', 'CO2', 'H2'). Values are initial mole amounts (numbers). You MUST include ALL components mentioned by the user, even if some have zero initial amounts. Example for methane steam reforming with 0.5 mole CH4 and 1 mole H2O: {\"CH4\": 0.5, \"H2O\": 1.0, \"CO\": 0.0, \"CO2\": 0.0, \"H2\": 0.0}",
                        "additionalProperties": {
                            "type": "number"
                        }
                    }
                },
                "required": ["SelectedComponents"]
            }
        }
    }

