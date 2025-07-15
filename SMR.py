import matplotlib.pyplot as plt
from tes_thermo.utils import Component
from tes_thermo.gibbs import Gibbs
import numpy as np

components_data = {
    'methane': {
        'Tc': 190.6, 'Tc_unit': 'K',
        'Pc': 45.99, 'Pc_unit': 'bar',
        'omega': 0.012,
        'Vc': 98.6, 'Vc_unit': 'cm³/mol',
        'Zc': 0.286, 
        'deltaHf': -74520 / 1000, 'deltaHf_unit': 'kJ/mol',
        'deltaGf': -50460 / 1000, 'deltaGf_unit': 'kJ/mol',
        'structure': {"C": 1, "H": 4},
        'phase': 'g'
    },
    'water': {
        'Tc': 647.1, 'Tc_unit': 'K',
        'Pc': 220.55, 'Pc_unit': 'bar',
        'omega': 0.345,
        'Vc': 55.9, 'Vc_unit': 'cm³/mol',
        'Zc': 0.229, 
        'deltaHf': -241818 / 1000, 'deltaHf_unit': 'kJ/mol',
        'deltaGf': -228572 / 1000, 'deltaGf_unit': 'kJ/mol',
        'structure': {"H": 2, "O": 1},
        'phase': 'g'
    },
    'carbon_monoxide': {
        'Tc': 132.9, 'Tc_unit': 'K',
        'Pc': 34.99, 'Pc_unit': 'bar',
        'omega': 0.048,
        'Vc': 93.4, 'Vc_unit': 'cm³/mol',
        'Zc': 0.294,
        'deltaHf': -110525 / 1000, 'deltaHf_unit': 'kJ/mol',
        'deltaGf': -137169 / 1000, 'deltaGf_unit': 'kJ/mol',
        'structure': {"C": 1, "O": 1},
        'phase': 'g'
    },
    'carbon_dioxide': {
        'Tc': 304.2, 'Tc_unit': 'K',
        'Pc': 73.83, 'Pc_unit': 'bar',
        'omega': 0.224,
        'Vc': 94.0, 'Vc_unit': 'cm³/mol',
        'Zc': 0.274,
        'deltaHf': -393509 / 1000, 'deltaHf_unit': 'kJ/mol',
        'deltaGf': -394359 / 1000, 'deltaGf_unit': 'kJ/mol',
        'structure': {"C": 1, "O": 2},
        'phase': 'g'
    },
    'hydrogen': {
        'Tc': 33.19, 'Tc_unit': 'K',
        'Pc': 13.13, 'Pc_unit': 'bar',
        'omega': -0.216,
        'Vc': 64.1, 'Vc_unit': 'cm³/mol',
        'Zc': 0.305,
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
    },
    'methanol': {
        'Tc': 512.6, 'Tc_unit': 'K',
        'Pc': 80.97, 'Pc_unit': 'bar',
        'omega': 0.564,
        'Vc': 118.0, 'Vc_unit': 'cm³/mol', 
        'Zc': 0.224,
        'deltaHf': -200660 / 1000, 'deltaHf_unit': 'kJ/mol',
        'deltaGf': -161960 / 1000, 'deltaGf_unit': 'kJ/mol',
        'structure': {"C": 1, "H": 4, "O": 1},
        'phase': 'g'
    }
}

comps = Component.create(components_data)

cp_coeffs = {
    'methane':          {'a': 1.702, 'b': 9.081e-3, 'c': -2.164e-6, 'd': 0},
    'water':            {'a': 3.470, 'b': 1.450e-3, 'c': 0,         'd': 12100},
    'carbon_monoxide':  {'a': 3.376, 'b': 0.557e-3, 'c': 0,         'd': -3100},
    'carbon_dioxide':   {'a': 5.457, 'b': 1.045e-3, 'c': 0,         'd': -115700},
    'hydrogen':         {'a': 3.249, 'b': 0.422e-3, 'c': 0,         'd': 8300},
    'carbon':           {'a': 1.77,  'b': 0.771e-3,     'c': 0,         'd': -86700},
    'methanol':         {'a': 2.211, 'b': 12.216e-3,'c': -3.450e-6, 'd': 0}
}


def cp(a, b, c, d):
    R = 8.314  # Ideal gas constant in J/(mol*K)
    def cp_function(T):
        return R * (a + b * T + c * T**2 + d / T**2)
    return cp_function


gibbs = Gibbs(components = comps,
              cp_polynomial_factory=cp,
              cp_coefficients=cp_coeffs,
              equation="Ideal Gas",)

res = []
T_range = np.linspace(773.15, 1273.15, 10)  # Temperature range from 600K to 1200K
for t in T_range:
    res.append(gibbs.solve_gibbs(initial=np.array([1, 1, 0, 0, 0, 0, 0]),
                  T=t, P=1, T_unit='K', P_unit='bar',))
    
temperatures = [r["Temperature (K)"] for r in res]
components = ["Methane", "Water", "Carbon monoxide", "Carbon dioxide", "Hydrogen", "Carbon", "Methanol"]
data = {comp: [r[comp] for r in res] for comp in components}

plt.figure(figsize=(10, 6))
for comp in components:
    plt.plot(temperatures, data[comp], label=comp)

plt.xlabel("Temperature (K)")
plt.ylabel("Mols")
plt.title("Species composition vs Temperature")
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend()
plt.tight_layout()
plt.show()