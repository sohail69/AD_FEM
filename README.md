# AD FEM Overview
An attempt at making a fully functional AD framework in MFEM, for problems which
can be described as constrianed minimals of energy functionals. The general problem
setup is as follows, suppose you want to minimise an energy function subject to
constraints(A):

```math
\displaylines{Min[e(u_1, u_2,...,u_n)] \\
              A.u = c}
```
In finite element problems we want this functional to be stationary with regards to
some defined DOFs, this means that the first derivative of the energy functional with
regards to the DOF's of interest must be zero, this is the residual form:

```math
\displaylines{R_m(u) = \frac{\partial e(u)}{\partial u_m} = 0}
```

# Picard-Iteration and Newtons Method
To find a solution (u) such that the possibly (generally) non-linear residual form is zero
a number of methods can be used to solve this problem the most popular being Picard iteration
and Newtons-Method. In the Case of Picard Iteration only the residual is needed to solve the
problem, however in Newtons method the Jacobian is also needed.

Picards Iteration is given by:
```math
\displaylines{u_{m}^{i+1} = u_{m}^{i} +  R_m(u^{i})}
```
Newtons Method is given by:
```math
\displaylines{u_{m}^{i+1} = u_{m}^{i} +  J_{mn}^{-1}(u^{i}) R_n(u^{i})}
```
The Jacobian is given by:
```math
\displaylines{J_mn(u) = \frac{\partial e(u)}{\partial u_m \partial u_n}}
```

# Dual Numbers and Forward Auto-Diff
Dual number are numbers that have the property:
.. math:: \epsilon \neq 0
   :label: HyperRealP1

.. math:: \epsilon * \epsilon = 0
   :label: HyperRealP2

```math
\displaylines{  \epsilon \neq 0  \\ \epsilon * \epsilon = 0}
```
