from tes_thermo.utils import UnitConverter
from types import SimpleNamespace
from tes_thermo.utils.logging import setup_logger
import numpy as np

logger = setup_logger(__name__)

class Component:
    """
    Class representing a chemical component.
    This class is used to create and manage chemical components with their properties.
    
    Initialization parameters:
        components: list[str] - list with component names. In this case, the properties are taken from thermo package.
        new_component: dict - dictionary with properties of the new component to be added.

    The format of the new_component dictionary should be:
    {
        "component_name": {
            "name": "methane",
            "Tc": 190.6, "Tc_unit": "K",
            "Pc": 45.99, "Pc_unit": "bar",
            "omega": 0.012,
            "Vc": 98.6, "Vc_unit": "cm³/mol",
            "Zc": 0.286,
            "Hfgm": -74873 / 1000, "Hfgm_unit": "kJ/mol",
            "Gfgm": -50870 / 1000, "Gfgm_unit": "kJ/mol",
            "structure": {"C": 1, "H": 4},
            "phase": "g",
            "kijs": [0, 0, 0, 0], # Interaction parameters with other components (OPCIONAL)
            "cp_polynomial": lambda T: 2.211 + 12.216e-3 * T - 3.450e-6 * T**2
        }
    }
    """
    
    def __init__(self, 
                 components: list[str] = None,      # list with component names
                 new_component: dict = None,        # dictionary with properties of new components
                 **kwargs):                         # for individual component properties
        
        self.components = components or []
        self.new_component = new_component or {}
        
        # For individual component creation (used by create method)
        for key, value in kwargs.items():
            setattr(self, key, value)

    @classmethod
    def create(cls, components_data):
        """
        Args:
            components_data (dict): A dictionary with component data.

        Returns:
            list: A list of fully configured Component objects.
        """
        components = []

        property_mapping = {
            'Tc': UnitConverter.convert_temperature,
            'Pc': UnitConverter.convert_pressure,
            'omega': None,         # Dimensionless
            'Vc': UnitConverter.convert_volume,
            'Zc': None,            # Dimensionless
            'Hfgm': UnitConverter.convert_energy,
            'Gfgm': UnitConverter.convert_energy,
            'phase': None,         # Optional, no conversion
            'structure': None,     # Optional, no conversion
            'kijs': None,          # Optional, no conversion - REMOVIDO DA OBRIGATORIEDADE
            'cp_coeffs': None,     # Optional, no conversion
            'cp_polynomial': None  # Optional, no conversion
        }

        default_values = {
            'phase': None,
            'structure': None,
            'kijs': None,          # Agora None por padrão
            'cp_coeffs': None,
            'cp_polynomial': None
        }

        for name, properties in components_data.items():
            processed_props = {'name': name}
            processed_props.update(default_values)

            for prop_name, converter in property_mapping.items():
                if prop_name in properties:
                    value = properties[prop_name]

                    if converter:
                        unit_key = f"{prop_name}_unit"
                        if unit_key in properties:
                            unit = properties[unit_key]
                            try:
                                processed_props[prop_name] = converter(value, unit)
                            except ValueError as e:
                                print(f"Error converting {prop_name} for {name}: {e}")
                        else:
                            print(f"Warning: Unit not specified for {prop_name} in {name}. Assuming base units.")
                            processed_props[prop_name] = value
                    else:
                        processed_props[prop_name] = value

            component = cls(**processed_props)
            components.append(component)

        return components
    
    def get_properties(self):
        properties = {}

        expected_properties = [
            'name', 'Tc', 'Pc', 'omega', 'Vc', 'Zc', 'Hfgm', 'Gfgm', 
            'phase', 'structure', 'kijs', 'cp_coeffs', 'cp_polynomial'
        ]

        for attr_name in expected_properties:
            if hasattr(self, attr_name):
                properties[attr_name] = getattr(self, attr_name)
            else:
                if attr_name == 'kijs':
                    properties[attr_name] = None  # Explicitamente None se não definido
                elif attr_name in ['phase', 'structure', 'cp_coeffs', 'cp_polynomial']:
                    properties[attr_name] = None
        
        return SimpleNamespace(**properties)

    def get_components(self):
        """
        Returns:
            dict: Dictionary with 'components' and 'new_components' keys
        """
        result = {
            "components": self.components,
            "new_components": []
        }
        
        if self.new_component:
            new_components_list = self.create(self.new_component)
            for comp in new_components_list:
                properties = comp.get_properties()
                properties_dict = vars(properties)
                result["new_components"].append(properties_dict)
        
        return result