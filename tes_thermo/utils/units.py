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

def convert_temperature_to_K(value: float, unit: str) -> float:
    if value:
        unit = unit.lower().strip()
        if unit == "k":
            return value
        elif unit == "c" or unit == "°c" or unit == "celsius":
            return value + 273.15
        elif unit == "f" or unit == "°f" or unit == "fahrenheit":
            return (value - 32) * 5/9 + 273.15
        else:
            raise ValueError(f"Unidade de temperatura '{unit}' não reconhecida. Use 'K', 'C' ou 'F'.")
    
def convert_pressure_to_bar(value: float, unit: str) -> float:
    if value:
        unit = unit.lower().strip()
        if unit == "bar":
            return value
        elif unit in ["pa", "pascal"]:
            return value / 1e5
        elif unit in ["kpa"]:
            return value / 100
        elif unit in ["mpa"]:
            return value * 10
        elif unit in ["atm"]:
            return value * 1.01325
        elif unit in ["psi"]:
            return value * 0.0689476
        else:
            raise ValueError(f"Unidade de pressão '{unit}' não reconhecida.")
        
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
    'Pa': 1.0,           
    'kPa': 1e3,           
    'MPa': 1e6,           
    'bar': 1e5,            
    'atm': 101325.0,     
    'psi': 6894.76,        
    'mmHg': 133.322,        
    'torr': 133.322         
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