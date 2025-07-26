from langchain.tools import BaseTool
from tes_thermo.utils.prompts import Prompts
from typing import Optional, Dict
import numpy as np
import pandas as pd
from pydantic import BaseModel, Field
from tes_thermo.utils.units import convert_pressure_to_bar, convert_temperature_to_K
from tes_thermo.thermo_agent.ming import Gibbs

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
        gibbs = Gibbs(components=components,
                      equation=Equation)
        
        Tmin = convert_temperature_to_K(Tmin, Tunit)
        Tmax = convert_temperature_to_K(Tmax, Tunit)
        Pmin = convert_pressure_to_bar(Pmin, Punit)
        Pmax = convert_pressure_to_bar(Pmax, Punit)

        if Tmin != Tmax:
            TRange = np.linspace(Tmin, Tmax, 10)
        else:
            TRange = np.linspace(Tmin, Tmax, 1)
        if Pmin != Pmax:
            PRange = np.linspace(Pmin, Pmax, 10)
        else:
            PRange = np.linspace(Pmin, Pmax, 1)

        all_results = []
        for T in TRange:
            for P in PRange:
                equilibrium_moles = gibbs.solve_gibbs(initial=compositions,T=T,P=P)
                
                row_data = {'T': T, 
                            'P': P}
                component_data = dict(zip(components, 
                                          equilibrium_moles))
                row_data.update(component_data)
                all_results.append(row_data)
        return pd.DataFrame(all_results)