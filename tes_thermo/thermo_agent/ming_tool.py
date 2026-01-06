from langchain.tools import BaseTool
from tes_thermo.utils.prompts import Prompts
from typing import Optional, Dict
import numpy as np
import pandas as pd
from pydantic import BaseModel, Field
from tes_thermo.utils.units import convert_pressure_to_bar, convert_temperature_to_K
from tes_thermo.gibbs import Gibbs
from tes_thermo.utils import Component

class MinGInputs(BaseModel):
    Tmin: Optional[float] = Field(default=600.0, 
                                  description="Minimum temperature value (e.g., 600).")
    Tmax: Optional[float] = Field(default=1200.0, 
                                  description="Maximum temperature value (e.g., 1200).")
    Tunit: Optional[str] = Field(default="K", 
                                 description="Unit of measurement for temperature (e.g., K, F, C).")
    Pmin: Optional[float] = Field(default=1.0, 
                                  description="Minimum pressure value (e.g., 1).")
    Pmax: Optional[float] = Field(default=10.0, 
                                  description="Maximum pressure value (e.g., 10).")
    Punit: Optional[str] = Field(default="bar", 
                                 description="Unit of measurement for pressure (e.g., bar, Pa, MPa).")
    Equation: Optional[str] = Field(default="Peng-Robinson", 
                                    description="This parameter defines the equation of state specified by the user. If not selected, consider 'Peng-Robinson' as the default.")
    SelectedComponents: Optional[Dict[str, float]] = Field(
        default_factory=dict,
        description="Dictionary with selected components as keys and their quantities as values. The chemical components must be presented as formula: CH4, CO2 ..."
    )

class MinG(BaseTool):
    name: str = "minG"
    description: str = Prompts.ming()
    args_schema = MinGInputs

    def _run(self, 
             Tmin: float = 600, 
             Tmax: float = 1200.0, 
             Tunit: str = "K", 
             Pmin: float = 1, 
             Pmax: float = 10, 
             Punit: str = "bar",
             Equation:str = "Peng-Robinson",
             SelectedComponents:list = None) -> str:

 
        components = [k for k, v in SelectedComponents.items()]
        compositions = [v for v in SelectedComponents.values()]
        
        # Create component structure using Component class
        comp_obj = Component(components=components, new_component={})
        components_dict = comp_obj.get_components()
        
        # Initialize Gibbs with the correct structure
        gibbs = Gibbs(components=components_dict,
                      equation=Equation)
        
        Tmin_K = convert_temperature_to_K(Tmin, Tunit)
        Tmax_K = convert_temperature_to_K(Tmax, Tunit)
        Pmin_bar = convert_pressure_to_bar(Pmin, Punit)
        Pmax_bar = convert_pressure_to_bar(Pmax, Punit)

        if Tmin_K != Tmax_K:
            TRange = np.linspace(Tmin_K, Tmax_K, 10)
        else:
            TRange = np.linspace(Tmin_K, Tmax_K, 1)
        if Pmin_bar != Pmax_bar:
            PRange = np.linspace(Pmin_bar, Pmax_bar, 10)
        else:
            PRange = np.linspace(Pmin_bar, Pmax_bar, 1)

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
                # The solution dict has component names as keys (capitalized and formatted)
                # The component_names in Gibbs are in the same order as the input components
                # So we can map by index
                solution_keys = [k for k in solution.keys() if k not in ["Temperature (K)", "Pressure (bar)"]]
                component_data = {}
                
                # Map by index - solution_keys should be in the same order as gibbs.component_names
                # which should match the order of input components
                for i, comp_name in enumerate(components):
                    if i < len(solution_keys):
                        # Use the solution key at the same index
                        component_data[comp_name] = solution[solution_keys[i]]
                    else:
                        # Fallback: use original composition if index is out of range
                        component_data[comp_name] = compositions[i] if i < len(compositions) else 0.0
                
                row_data.update(component_data)
                all_results.append(row_data)
        return pd.DataFrame(all_results)