"""
    Represents a chemical component and its physical properties.
    
    This class acts as a data container for properties like critical temperature,
    pressure, acentric factor, and more.
    
    Instances are best created using the `create_from_dict` class method,
    which can process a dictionary containing data for multiple components,
    handle unit conversions automatically, and return a list of Component objects.

    Example Usage:

    from tes.utils import Component

    data = {
        'methane': {
            'Tc': 190.56, 'Tc_unit': 'K',
            'Pc': 45.99, 'Pc_unit': 'bar',
            'omega': 0.011,
            'Vc': 98.6, 'Vc_unit': 'cm³/mol',
            'Zc': 0.286,
            'deltaHf': -74.81, 'deltaHf_unit': 'kJ/mol',
            'deltaGf': -50.72, 'deltaGf_unit': 'kJ/mol'
        },
        'ethane': {
            'Tc': 305.32, 'Tc_unit': 'K',
            'Pc': 48.72, 'Pc_unit': 'bar',
            'omega': 0.100
        }
    }

    my_components = Component.create(data)
    for comp in comps:
    print(comp.get_properties())

    Response:
    {'name': 'methane', 'Tc': 190.56, 'Pc': 4599000.0, 'omega': 0.011, 'Vc': 9.859999999999998e-05, 'Zc': 0.286, 'deltaHf': -74810.0, 'deltaGf': -50720.0}
    {'name': 'ethane', 'Tc': 305.32, 'Pc': 4872000.0, 'omega': 0.099, 'Vc': 0.00014549999999999999, 'Zc': 0.279, 'deltaHf': -84000.0, 'deltaGf': -32000.0}
"""

from tes.utils import UnitConverter

class Component:
    
    def __init__(self, 
                 name=None, 
                 Tc=None, 
                 Pc=None, 
                 omega=None, 
                 Vc=None, 
                 Zc=None, 
                 deltaHf=None, 
                 deltaGf=None,
                 phase=None,
                 structure=None):
        self.name = name
        self.Tc = Tc                # Critical Temperature (K)
        self.Pc = Pc                # Critical Pressure (Pa)
        self.omega = omega          # Acentric Factor (-)
        self.Vc = Vc                # Critical Volume (m³/mol)
        self.Zc = Zc                # Critical Compressibility Factor (-)
        self.deltaHf = deltaHf      # Standard Enthalpy of Formation (J/mol)
        self.deltaGf = deltaGf      # Standard Gibbs Energy of Formation (J/mol)
        self.phase = phase          # Phase of the component (e.g., 'g' for gas, 's' for solid)
        self.structure = structure  # Chemical structure as a dictionary (e.g., {"C":

    @classmethod
    def create(cls, components_data):
        """
        Factory method that processes a dictionary where each key is a component
        name and its value is a dictionary of properties. It automatically
        handles unit conversions to the standard base units (K, Pa, m³/mol, J/mol),
        and processes optional fields like phase and structure.

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
            'deltaHf': UnitConverter.convert_energy,
            'deltaGf': UnitConverter.convert_energy,
            'phase': None,         # Optional, no conversion
            'structure': None      # Optional, no conversion
        }

        for name, properties in components_data.items():
            processed_props = {'name': name}

            for prop_name, converter in property_mapping.items():
                if prop_name in properties:
                    value = properties[prop_name]

                    if converter:
                        # Unit conversion required
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
                        # No conversion needed (dimensionless or optional)
                        processed_props[prop_name] = value

            component = cls(**processed_props)
            components.append(component)

        return components

    def get_properties(self):
        return {
            'name': self.name,
            'Tc': self.Tc,
            'Pc': self.Pc,
            'omega': self.omega,
            'Vc': self.Vc,
            'Zc': self.Zc,
            'deltaHf': self.deltaHf,
            'deltaGf': self.deltaGf,
            'phase': self.phase,
            'structure': self.structure
        }