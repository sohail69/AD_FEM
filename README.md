# AD FEM Overview
An attempt at making a fully functional AD framework in MFEM, for problems which
can be described as constrained minimals of energy functionals. The general problem
setup is as follows, suppose you want to minimise an energy function subject to
constraints(A):

```math
\displaylines{Min[e(u_1, u_2,...,u_n)] \\
              A.u = c}
```
Some example functionals of integrated/sampled continuous variables include:
```math
\displaylines{ e_{LE}(u) = \int_{\Omega} \left( \frac{1}{2} \vec{u} \cdot K \vec{u} 
                         + \vec{F} \cdot \vec{u} \right) d\Omega \quad Linear-elasticity \\
               e_{HE}(T) = \int_{\Omega} \left( T \cdot \frac{\partial T}{\partial t} 
                         + \nabla T \cdot \nabla T \right) d\Omega \quad Heat-equation \\
               e_{HY}(u,p) = \int_{\Omega} \left( \delta_{ij} C_{ij}(u) + P ln(J(u)) 
                           - \frac{P^{2}}{2\lambda} \right) d\Omega \quad Mixed-NeoHookean-equation}
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
\displaylines{J_{mn}(u) = \frac{\partial e(u)}{\partial u_m \partial u_n}}
```

# Dual Numbers and Forward Auto-Diff
Dual number are numbers that have the property:
```math
\displaylines{  \epsilon \neq 0  \\
                \epsilon * \epsilon = 0}
```
A dual number can be added and multiplied by any other type of number
without issue giving another dual number:
```math
\displaylines{  \epsilon + \epsilon =  2 \epsilon  \\
                5 * \epsilon = 5 \epsilon          \\
                5 + \epsilon = 5 + \epsilon        \\
                (5 + \epsilon)^{2} = 25 + 10 \epsilon}
```
Dual numbers in conjunction with templated functions can be used to
exactly calculate derivatives. Given a function f(x) using limits the derivative
of that function is defined as:
```math
\displaylines{ f'(x) = Lim_{h \rightarrow 0} \frac{ f(x+h) - f(x) }{h} }
```
In the context above (h) can be replaced by a dual number, such that the equation
becomes transformed into:
```math
\displaylines{   f'(x)*\epsilon = f(x+\epsilon) - f(x) }
```
The key difference between this and the equation above is that the limits don't need
to be taken and neither do manual derivatives. It may seem like a simple multiplication
but rather the Non-dual part of the number cancels out and the remainder left is the dual
part multiplied by a coefficient, the coefficient is the exact evaluated derivative of the
function.

The C++ dual number template is as follows:
```c++
val_t  val;
grad_t grad;
dualNumber<val_t,grad_t> x(val,grad);
```
To calculate the first derivative (residual) only the dual number of the type of the
precision of float is needed setting the grad component to one and evaluating the function, e.g. 
```c++
double val=10.0 /*some arbitrary value for function eval*/
dualNumber<double,double> x(val,1.0);
double df = f(x).grad;
```

For the second derivative (Jacobian) a slightly more sophisticated route is needed. you need a dual
number of dual numbers
```c++
double val=10.0 /*some arbitrary value for function eval*/
dualNumber<double,double> a(val,1.0), b(1.0,0.0);
dualNumber<dualNumber<double,double>,<dualNumber<double,double>> x(a,b);
double d2f = f(x).grad.grad;
```
For Vector input functionals (scalars with vector inputs) this makes the direct evaulation of the
residual and Jacobian quite natural on the element level without having to manually calculate 
the derivatives of the functional. One just perturbs each component of the element vector 
(one by one) and returns the appropriate dual component.

# Finite element continuous function sampling
Finite element weak forms evaluate integrals, these integrals are often (except in special cases) 
approximated by a weighted sampling rule, where the discrete DOF variables are sampled as (piecewise) continuous
variables at the integration points. The approximation of the integral can be defined as follows:
```math
\displaylines{  e(u) = \int_{\Omega} f(u) d\Omega \approx \sum^{N_{ip}}_{ip=1} f(u_{ip}) \cdot det(J) \cdot w_{ip} }
```
The variable (u_ip) is the sampled continuous variable which is given by:
```math
\displaylines{u^{ip} = H^{ip}_{m} \tilde{u}_{m} }
```
Where (H) is the discrete field interpolator sampled at the integration point and ($$\tilde{u}$$) is the
discrete variable/DOF that is being interpolated/solved-for. The interpolation is often a linear 
matrix/function which weights. For example if the user wanted to get the sampled gradient of a field
equivalently :
```math
\displaylines{\frac{\partial u^{ip}}{\partial x_{j}} = \frac{\partial H^{ip}_{m}}{\partial x_{j}} \tilde{u}_{m} }
```
For the general case you may need a number of different derivatives and linear transforms of derivatives of
a particular variables e.g.
```math
\displaylines{ \vec{u}, \nabla(\vec{u}), \nabla \times (\vec{u}) }
```
These can all be described as a set of linear transforms of the discrete data multiplied by derivatives of 
the basis functions:
```math
\displaylines{ \left (\vec{u}, \nabla \vec{u}, \nabla \times \vec{u} \right) = \mathbf{Q} \tilde{u} }
```
Where (Q) is the lumped interpolation operator which combines all the basis (and derivatives) into a single
flat tensor output. Each tensor has an (standard structured) iterator which can be used to recover the
original structure of the variable. Some basis functions and integral types require additional information
to transform from the global to the local element frame besides simply element restriction, for example 
H(div) and H(curl) elements.
