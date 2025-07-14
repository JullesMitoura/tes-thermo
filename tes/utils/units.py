"""
A comprehensive unit conversion utility for scientific and engineering calculations.

This class provides a centralized and straightforward way to convert between various 
units of temperature, pressure, molar volume, and molar energy. All conversion 
methods are static, allowing them to be called directly from the class without 
needing to create an instance. The base units for conversion are Kelvin (K) for 
temperature, Pascal (Pa) for pressure, cubic meters per mole (m³/mol) for molar 
volume, and Joules per mole (J/mol) for molar energy.

The class is designed for ease of use and extensibility, with conversion factors 
and formulas stored in dictionaries. This makes it simple to add new units or 
modify existing ones.
"""
class UnitConverter:
    TEMPERATURE_CONVERSIONS = {
        'K': 1.0,
        'C': lambda x: x + 273.15,
        '°C': lambda x: x + 273.15,
        'F': lambda x: (x - 32) * 5/9 + 273.15,
        '°F': lambda x: (x - 32) * 5/9 + 273.15,
        'R': lambda x: x * 5/9  # Rankine to Kelvin
    }
    PRESSURE_CONVERSIONS = {
    'Pa': 1 / 100000.0,
    'kPa': 1 / 100.0,
    'MPa': 10.0,
    'bar': 1.0,
    'atm': 1.01325,
    'psi': 0.0689476,
    'mmHg': 0.00133322,
    'torr': 0.00133322
}
    VOLUME_CONVERSIONS = {
        'm³/mol': 1.0,
        'L/mol': 0.001,
        'cm³/mol': 1e-6,
        'mL/mol': 1e-6
    }
    ENERGY_CONVERSIONS = {
        'J/mol': 1.0,
        'kJ/mol': 1000.0,
        'cal/mol': 4.184,
        'kcal/mol': 4184.0,
        'BTU/mol': 1055.06
    }
    
    @staticmethod
    def convert_temperature(value, from_unit):
        """Converts temperature to Kelvin"""
        if from_unit in UnitConverter.TEMPERATURE_CONVERSIONS:
            converter = UnitConverter.TEMPERATURE_CONVERSIONS[from_unit]
            if callable(converter):
                return converter(value)
            else:
                return value * converter
        else:
            raise ValueError(f"Temperature unit '{from_unit}' not supported")
    
    @staticmethod
    def convert_pressure(value, from_unit):
        """Converts pressure to Pascal"""
        if from_unit in UnitConverter.PRESSURE_CONVERSIONS:
            return value * UnitConverter.PRESSURE_CONVERSIONS[from_unit]
        else:
            raise ValueError(f"Pressure unit '{from_unit}' not supported")
    
    @staticmethod
    def convert_volume(value, from_unit):
        """Converts molar volume to m³/mol"""
        if from_unit in UnitConverter.VOLUME_CONVERSIONS:
            return value * UnitConverter.VOLUME_CONVERSIONS[from_unit]
        else:
            raise ValueError(f"Volume unit '{from_unit}' not supported")
    
    @staticmethod
    def convert_energy(value, from_unit):
        """Converts energy to J/mol"""
        if from_unit in UnitConverter.ENERGY_CONVERSIONS:
            return value * UnitConverter.ENERGY_CONVERSIONS[from_unit]
        else:
            raise ValueError(f"Energy unit '{from_unit}' not supported")