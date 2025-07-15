from tes.utils import Component
from tes.gibbs import Gibbs
import numpy as np

components_data = {
    'methane': {
        'Tc': 190.60, 'Tc_unit': 'K',
        'Pc': 45.99, 'Pc_unit': 'bar',
        'omega': 0.01,
        'Vc': 98.60, 'Vc_unit': 'cm³/mol',
        'Zc': 0.29,
        'deltaHf': -74528 / 1000, 'deltaHf_unit': 'kJ/mol', # Converted from J/mol
        'deltaGf': -50460 / 1000, 'deltaGf_unit': 'kJ/mol', # Converted from J/mol
        'structure': {"C": 1, "H": 4},
        'phase': 'g'
    },
    'water': {
        'Tc': 647.10, 'Tc_unit': 'K',
        'Pc': 220.55, 'Pc_unit': 'bar',
        'omega': 0.35,
        'Vc': 55.90, 'Vc_unit': 'cm³/mol',
        'Zc': 0.23,
        'deltaHf': -241818 / 1000, 'deltaHf_unit': 'kJ/mol',
        'deltaGf': -228572 / 1000, 'deltaGf_unit': 'kJ/mol',
        'structure': {"H": 2, "O": 1},
        'phase': 'g'
    },
    'carbon_monoxide': {
        'Tc': 132.90, 'Tc_unit': 'K',
        'Pc': 34.99, 'Pc_unit': 'bar',
        'omega': 0.05,
        'Vc': 93.40, 'Vc_unit': 'cm³/mol',
        'Zc': 0.30,
        'deltaHf': -110525 / 1000, 'deltaHf_unit': 'kJ/mol',
        'deltaGf': -137169 / 1000, 'deltaGf_unit': 'kJ/mol',
        'structure': {"C": 1, "O": 1},
        'phase': 'g'
    },
    'carbon_dioxide': {
        'Tc': 304.20, 'Tc_unit': 'K',
        'Pc': 73.83, 'Pc_unit': 'bar',
        'omega': 0.22,
        'Vc': 94.0, 'Vc_unit': 'cm³/mol',
        'Zc': 0.27,
        'deltaHf': -393509 / 1000, 'deltaHf_unit': 'kJ/mol',
        'deltaGf': -394359 / 1000, 'deltaGf_unit': 'kJ/mol',
        'structure': {"C": 1, "O": 2},
        'phase': 'g'
    },
    'hydrogen': {
        'Tc': 33.19, 'Tc_unit': 'K',
        'Pc': 13.13, 'Pc_unit': 'bar',
        'omega': -0.22,
        'Vc': 64.10, 'Vc_unit': 'cm³/mol',
        'Zc': 0.31,
        'deltaHf': 0.0, 'deltaHf_unit': 'kJ/mol',
        'deltaGf': 0.0, 'deltaGf_unit': 'kJ/mol',
        'structure': {"H": 2},
        'phase': 'g'
    },
    'carbon': {
        'Tc': 0, 'Tc_unit': 'K',
        'Pc': 0, 'Pc_unit': 'bar',
        'omega': 0,
        'Vc': 0, 'Vc_unit': 'cm³/mol',
        'Zc': 0,
        'deltaHf': 0.0, 'deltaHf_unit': 'kJ/mol',
        'deltaGf': 0.0, 'deltaGf_unit': 'kJ/mol',
        'structure': {"C": 1},
        'phase': 's'
    }
}

comps = Component.create(components_data)

cp_coeffs = {
    'methane':          {'a': 1.70, 'b': 0.01,   'c': 0.00,    'd': 0},
    'water':            {'a': 3.47, 'b': 0.00,   'c': 0,       'd': 12100},
    'carbon_monoxide':  {'a': 3.38, 'b': 0.00,   'c': 0,       'd': -3100},
    'carbon_dioxide':   {'a': 5.46, 'b': 0.00,   'c': 0,       'd': -115700},
    'hydrogen':         {'a': 3.25, 'b': 0.00,   'c': 0,       'd': 8300},
    'carbon':           {'a': 1.77, 'b': 0.00,   'c': 0,       'd': -86700}
}


def cp(a, b, c, d):
    R = 8.314  # Ideal gas constant in J/(mol*K)
    def cp_function(T):
        return R * (a + b * T + c * T**2 + d / T**2)
    return cp_function


gibbs = Gibbs(components = comps,
              cp_polynomial_factory=cp,
              cp_coefficients=cp_coeffs,
              equation="Virial")

print(gibbs.solve_gibbs(initial=np.array([1, 1, 0, 0, 0, 0]),
                  T=1200, P=1, T_unit='K', P_unit='bar',))
