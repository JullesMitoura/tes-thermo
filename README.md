# TeS - Themodynamic Equilibrium Simulation

TeS - Thermodynamic Equilibrium Simulation is an open-source software designed to optimize studies in thermodynamic equilibrium and related subjects. TeS is recommended for initial analyses of reactional systems. The current version contains the following simulation module:

### 1. Gibbs Energy Minimization (minG):

This module allows the user to simulate an isothermal reactor using the Gibbs energy minimization approach. References on the mathematical development can be found in previous work reported by Mitoura and Mariano (2024).

As stated, the objective is to minimize the Gibbs energy, which is formulated as a non-linear programming problem, as shown in the equation below:

$$min G = \sum_{i=1}^{NC} \sum_{j=1}^{NF} n_i^j \mu_i^j$$

The next step is the calculation of the Gibbs energy. The equation below shows the relationship between enthalpy and heat capacity.

$$\frac{\partial \bar{H}_i^g}{\partial T} = Cp_i^g \text{  para } i=1,\ldots,NC$$

Knowing the relationship between enthalpy and temperature, the next step is to calculate the chemical potential. The equation below presents the correlation for calculating chemical potentials.

$$\frac{\partial}{\partial T} \left( \frac{\mu_i^g}{RT} \right) = -\frac{\bar{H}_i^g}{RT^2} \quad \text{para } i=1,\ldots,NC$$

We then have the calculation of the chemical potential for component i:

$$
\mu_i^0 = \frac {T}{T^0} \Delta G_f^{298.15 K} - T \int_{T_0}^{T} \frac {\Delta H_f^{298.15 K} + \int_{T_0}^{T} (CPA + CPB \cdot T + CPC \cdot T^2 + \frac{CPD}{T^2}) \, dT}{T^2} \, dT
$$

With the chemical potentials known, we can define the objective function:

$$\min G = \sum_{i=1}^{NC} n_i^g \mu_i^g $$

Where:

$$\mu _i^g = \mu _i^0 + R.T.(ln(\phi_i)+ln(P)+ln(y_i)) $$

For the calculation of fugacity coefficients, we will have two possibilities:

1. Ideal Gas:

$$\phi = 1 $$

2. Non-ideal Gas:
For non-ideal gases, the calculation of fugacity coefficients is based on the Virial equation of state, as detailed in section 1.1.

The space of possible solutions must be restricted by two conditions:
1. Non-negativity of moles:

$$ n_i^j \geq 0 $$

2. Conservation of atoms:

$$
\sum_{i=1}^{NC} a_{mi} \left(\sum_{j=1}^{NF} n_{i}^{j}\right) = \sum_{i=1}^{NC} a_{mi} n_{i}^{0}
$$

References:

Mitoura, Julles.; Mariano, A.P. Gasification of Lignocellulosic Waste in Supercritical Water: Study of Thermodynamic Equilibrium as a Nonlinear Programming Problem. Eng 2024, 5, 1096-1111. https://doi.org/10.3390/eng5020060

### 1.1 Fugacity Coefficient Calculation:

### Virial Equation (2nd Term)

The Virial equation truncated at the second term relates the compressibility factor to pressure:

$$Z = 1 + \frac{B_{mix} P}{RT}$$

The second Virial coefficient for the mixture is calculated using the following mixing rule:

$$B_{mix} = \sum_{i=1}^{NC} \sum_{j=1}^{NC} y_i y_j B_{ij}$$

The logarithm of the fugacity coefficient for each component i in the mixture is given by:

$$\ln \phi_i = \left[ 2 \sum_{j=1}^{NC} y_j B_{ij} - B_{mix} \right] \frac{P}{RT}$$

Finally, for any of the models:

$$\phi_i = \exp(\ln \phi_i)$$


---
### Usage Example:
#### Methane Steam Reforming Process

First, you need to install tes-thermo:

```python
pip intsall -qU tes-thermo
```
Now you have access to tes-thermo code. With this, you just need to import:

```python
from tes_thermo.utils import Component
from tes_thermo.gibbs import Gibbs
import numpy as np
```

Define the componentes propriets:

```python
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
```

With this informations, you can create your components instances using `Component`:

```python
comps = Component.create(components_data)
```

The next step is define a polynomial to calculte `Cp`:

```python
def cp(a, b, c, d):
    R = 8.314  # Ideal gas constant in J/(mol*K)
    def cp_function(T):
        return R * (a + b * T + c * T**2 + d / T**2)
    return cp_function
```

And define the values of coefficientes:

```python

cp_coeffs = {
    'methane':          {'a': 1.702, 'b': 9.081e-3, 'c': -2.164e-6, 'd': 0},
    'water':            {'a': 3.470, 'b': 1.450e-3, 'c': 0,         'd': 12100},
    'carbon_monoxide':  {'a': 3.376, 'b': 0.557e-3, 'c': 0,         'd': -3100},
    'carbon_dioxide':   {'a': 5.457, 'b': 1.045e-3, 'c': 0,         'd': -115700},
    'hydrogen':         {'a': 3.249, 'b': 0.422e-3, 'c': 0,         'd': 8300},
    'carbon':           {'a': 1.77,  'b': 0.771e-3,     'c': 0,         'd': -86700},
    'methanol':         {'a': 2.211, 'b': 12.216e-3,'c': -3.450e-6, 'd': 0}
}

```

Now, you have all you need to simulate:

```python

gibbs = Gibbs(components = comps,
              cp_polynomial_factory=cp,
              cp_coefficients=cp_coeffs,
              equation="Ideal Gas",)

res = gibbs.solve_gibbs(initial=np.array([1, 1, 0, 0, 0, 0, 0]),
                  T=1200, P=1, T_unit='K', P_unit='bar',)
                  
res
```
As a result:

```python
{'Temperature (K)': 1200.0, 'Pressure (bar)': 220.0, 'Methane': 0.17607201808452386, 'Water': 1.0024576179445397, 'Carbon monoxide': 0.15036090499183952, 'Carbon dioxide': 0.17360440825679002, 'Hydrogen': 0.6454229503976403, 'Carbon': 2.7718349964766464e-09, 'Methanol': 2.685497942223023e-06}
```


---

### Third-Party Dependencies and Licenses

This project uses the Ipopt solver, which is made available under the Eclipse Public License v1.0 (EPL-1.0). A full copy of the Ipopt license can be verified here: https://github.com/coin-or/Ipopt/blob/stable/3.14/LICENSE

---
