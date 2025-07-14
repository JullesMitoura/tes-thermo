## Third-Party Dependencies and Licenses

This project uses the Ipopt solver, which is made available under the Eclipse Public License v1.0 (EPL-1.0). A full copy of the Ipopt license can be verified here: https://github.com/coin-or/Ipopt/blob/stable/3.14/LICENSE

## TeS - Themodynamic Equilibrium Simulation

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

The second Virial coefficient for the mixture ($B_{mix}$) is calculated using the following mixing rule:

$$B_{mix} = \sum_{i=1}^{NC} \sum_{j=1}^{NC} y_i y_j B_{ij}$$

Where $B_{ii}$ is the coefficient for the pure component and $B_{ii}$ is the cross-coefficient for the i-j pair. These coefficients are temperature-dependent and are generally obtained from empirical correlations based on critical properties.

The logarithm of the fugacity coefficient for each component i in the mixture is given by:

$$\ln \phi_i = \left[ 2 \sum_{j=1}^{NC} y_j B_{ij} - B_{mix} \right] \frac{P}{RT}$$

Finally, for any of the models:

$$\phi_i = \exp(\ln \phi_i)$$

For solid components, it is assumed that ($\phi_i = 1.0$).

IPOPT Solver:

In all simulation modules of this project, the core of the problem lies in finding the optimal solution for a system described by complex, non-linear equations. For this task, the IPOPT solver was chosen.

IPOPT (Interior Point Optimizer) is a high-performance, open-source software package designed specifically for solving large-scale Nonlinear Programming (NLP) problems.

The choice of IPOPT for TeS v3 was not accidental and is based on several key technical reasons:

Nature of the Thermodynamic Problem: Chemical and phase equilibrium problems are inherently non-linear. The objective functions (like Gibbs energy) and constraints depend on logarithms of mole fractions (ln(y_i)), complex equations of state, and other relationships that cannot be solved with linear programming methods. IPOPT is specialized for this class of problems.

Robust Constraint Handling: IPOPT uses an interior-point method (or barrier method), which is extremely effective for handling problems with both equality (e.g., atomic balance, energy balance) and inequality constraints (e.g., non-negativity of moles, n_i
ge0). It approaches the optimal solution from within the feasible region, making it robust for the types of constraints found in thermodynamics.

Convergence and Stability: Thermodynamic models can be numerically sensitive. IPOPT is recognized in the academic and industrial communities for its excellent performance and ability to converge to an optimal (or locally optimal) solution in a stable and efficient manner, even for complex and ill-conditioned functions.

Integration with Pyomo: The code utilizes the Pyomo modeling framework. IPOPT is one of the standard and best-integrated solvers with Pyomo, allowing the mathematical formulation of the problem to be translated directly and reliably into a format the solver understands. This synergy simplifies code development and maintenance.

Open-Source Philosophy: Like TeS, IPOPT is an open-source project (part of the COIN-OR initiative). Its use ensures that the software is completely free, accessible, and transparent, aligning with the project's goals of fostering study and research in thermodynamics.

In summary, IPOPT was chosen because it is a powerful, validated, open-source tool that is technically suited for the class of non-linear optimization problems found in thermodynamic equilibrium simulation, ensuring accurate and reliable results for TeS users.

The solver can be downloaded from this address:

https://github.com/coin-or/Ipopt/releases

---