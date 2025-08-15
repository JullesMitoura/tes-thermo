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

### Peng-Robinson and Soave-Redlich-Kwong:

Based on the chosen equation of state, the following parameters are defined:
$$\Omega_a$$
$$\Omega_b$$

#### Temperature-Adjusted Attraction Parameter

$$
m_i =
\begin{cases}
0.37464 + 1.54226 \cdot \omega_i - 0.26992 \cdot \omega_i^2 & \text{(Peng-Robinson)} \\
0.480 + 1.574 \cdot \omega_i - 0.176 \cdot \omega_i^2 & \text{(SRK)}
\end{cases}
$$

$$
\alpha_i = \left(1 + m_i(1 - \sqrt{T/T_{c,i}})\right)^2 \quad \text{(PR and SRK)}
$$

$$
a_i = \Omega_a \cdot \left( \frac{R^2 T_{c,i}^2}{P_{c,i}} \right) \cdot \alpha_i
\quad ; \quad
b_i = \Omega_b \cdot \left( \frac{R T_{c,i}}{P_{c,i}} \right)
$$

#### Binary Interaction Parameter
$$k_{ij}$$

$$a_{ij} = (1 - k_{ij}) \cdot \sqrt{a_i \cdot a_j}$$

$$
a_{\text{mix}} = \sum_i \sum_j y_i y_j a_{ij}
\quad ; \quad
b_{\text{mix}} = \sum_i y_i b_i
$$

$$
A = \frac{a_{\text{mix}} P}{R^2 T^2}
\quad ; \quad
B = \frac{b_{\text{mix}} P}{R T}
$$

The cubic equation is written as:
$$Z^3 + c_2 Z^2 + c_1 Z + c_0 = 0$$

The coefficients depend on the EOS:

### Peng-Robinson (PR):
$$Z^3 + (B - 1)Z^2 + (A - 2B - 3B^2)Z + (-AB + B^2 + B^3) = 0$$

### SRK:
$$Z^3 - Z^2 + (A - B - B^2)Z - AB = 0$$

#### Solution
Select the largest positive real root (Z) that represents the gas phase.

#### Fugacity Coefficient
For each component ($i$):
$$\ln \phi_i = \frac{b_i}{b_{\text{mix}}}(Z - 1) - \ln(Z - B) - \frac{A}{B} \cdot \left( \frac{2 \sum_j y_j a_{ij}}{a_{\text{mix}}} - \frac{b_i}{b_{\text{mix}}} \right) \cdot f(Z, B)$$

Where:

#### For PR:
$$f(Z, B) = \frac{1}{2\sqrt{2}} \cdot \ln\left( \frac{Z + (1 + \sqrt{2})B}{Z + (1 - \sqrt{2})B} \right)$$

#### For SRK:
$$f(Z, B) = \ln\left(1 + \frac{B}{Z} \right)$$

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
To define componentes:
```python
new_components = {
        "methane": {
            "name": "methane",
            "Tc": 190.6, "Tc_unit": "K",
            "Pc": 45.99, "Pc_unit": "bar",
            "omega": 0.012,
            "Vc": 98.6, "Vc_unit": "cm³/mol",
            "Zc": 0.286,
            "Hfgm": -74520, "Hfgm_unit": "J/mol",
            "Gfgm": -50460, "Gfgm_unit": "J/mol",
            "structure": {"C": 1, "H": 4},
            "phase": "g",
#            "kijs": [0, 0, 0, 0, 0, 0]
            "cp_polynomial": lambda T: 8.314 * (1.702 + 0.009081* T -0.000002164*T**2),
        }
    }

components = ['water','carbon monoxide', 'carbon dioxide', 'hydrogen', 'methanol']
```
In the example above, `new_components` refers to the components to be added. For these, the user must specify all thermodynamic properties as well as the polynomial to be used to calculate `Cp`. For this example, the following polynomial was used:

$$C_p(T) = R \times \left( 1.702 + 0.009081T - 0.000002164T^2 \right)$$

where $T$ is the temperature in Kelvin and $ C_p $ is the heat capacity in J/(mol·K).

``components`` refers to the components that will be queried using the thermo library, so it is not necessary to indicate thermodynamic properties.

Note that when adding a new component, adding ``kij`` values is optional. If not specified, they will be estimated using the critical volume values.

The next step is to instantiate the components using the ``Component`` class.

```python
comps = Component(components, new_components)
comps = comps.get_components()
gibbs = Gibbs(components=comps,equation='Peng-Robinson')
res = gibbs.solve_gibbs(T=800, T_unit= 'K',
                         P=60, P_unit='bar',
                         initial=np.array([0, 1, 0, 0, 1, 0]))
```

After defining the components, the ``Gibbs`` class is used, and the parameters must be the components and the equation of state to be used. The Gibbs class has the ``solve_gibbs`` method, where the user must specify the parameters to be considered in the simulation.

The results are as follows:
```python
{'Temperature (K)': 800.0,
 'Pressure (bar)': 60.00000000000001,
 'Water': 0.04337941510960611,
 'Carbon monoxide': 0.11003626275025695,
 'Carbon dioxide': 0.923292149916935,
 'Hydrogen': 0.023277493427407873,
 'Methanol': 5.225802852767064e-08,
 'Methane': 0.9666715055286526}
````
The current version of `tes-thermo` also includes `thermo-agent`. To verify its use, use the example shown in:

```
tes-thermo
├─ notebooks
│  ├─ smr.ipynb
│  └─ thermo_agent.ipynb <- This example!
```
Repository link: https://github.com/JullesMitoura/tes-thermo

---

### Third-Party Dependencies and Licenses

This project uses the Ipopt solver, which is made available under the Eclipse Public License v1.0 (EPL-1.0). A full copy of the Ipopt license can be verified here: https://github.com/coin-or/Ipopt/blob/stable/3.14/LICENSE

---