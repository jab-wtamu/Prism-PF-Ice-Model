---
title: Phase-field evolution (Demange ice-growth model)

---

# Phase-field evolution (Demange ice-growth model)
## 1. Strong form (kinetics)
The phase-field equation governing ice growth is given by (Demange et al.):
**(Eq.1)**
$$
\begin{aligned}
A(\mathbf{n})^2 \frac{\partial \phi}{\partial t} = f'(\phi) + \lambda B(\mathbf{n}) g'(\phi)\, u + \frac{1}{2}\nabla_\Gamma \cdot \left(
|\nabla \phi|^2
\frac{\partial \left(A(\mathbf{n})^2\right)}{\partial (\nabla \phi)}
\right) + \nabla_\Gamma \cdot \left(
A(\mathbf{n})^2 \nabla_\Gamma \phi
\right)
\end{aligned} 
$$
where:
* $\phi$ is the phase-field variable,
* $u$ is the supersaturation,
* n is $- \nabla \phi /|\nabla \phi|,$
* $\nabla_\Gamma=(\partial_x,\; \partial_y,\; \Gamma\,\partial_z)$
* A(n) and B(n) encode anisotropy

**Equation (2): Supersaturation evolution**
$$
\partial_t u=
\tilde{D}\,\nabla_\Gamma \cdot \left( q(\phi)\,\nabla_\Gamma u \right)-
\frac{L_{\mathrm{sat}}}{2}\,B(\mathbf{n})\,\partial_t \phi.
$$
where $u$ is the supersaturation field,
$q(\phi)=1-\phi,$
$\mathbf{n}=-\frac{\nabla \phi}{|\nabla \phi|},$
$\nabla_\Gamma = (\partial_x,\partial_y,\Gamma\,\partial_z).$



## Time discretization (Forward Euler)
Using explicit Forward Euler time stepping:
(Eq.1)


where:
$$F_1^n
=\frac{1}{2}\,|\nabla \phi^n|^2
\frac{\partial\!\left((A^n)^2\right)}{\partial (\nabla \phi^n)},      
F_2^n=(A^n)^2 \nabla_\Gamma \phi^n.
$$
Solving for $\phi^{n+1}$:
$$
\phi^{n+1}=\phi^n + \frac{\Delta t}{(A^n)^2}[f'(\phi^n)+ \lambda B^n g'(\phi^n)\,u^n + \nabla_\Gamma \cdot F_1^n + \nabla_\Gamma \cdot F_2^n]
$$

**(Eq.2)**

$$
u^{n+1}=
u^n+
\Delta t\,\tilde{D}\,\nabla_\Gamma \cdot \left(
q(\phi^n)\,\nabla_\Gamma u^n
\right)-\Delta t\
\frac{L_{\mathrm{sat}}}{2}\,B^n\,(\phi^{n+1}-\phi^n).
$$
where:
$$
q(\phi^n)=1-\phi^n,
\qquad
B^n = B(\mathbf{n}^n),
\qquad
\mathbf{n}^n=-\frac{\nabla \phi^n}{|\nabla \phi^n|}.
$$






## Weak Formulation
**(Eq.1)**

$$
\int_\Omega w \phi^{n+1} \, dV=\int_\Omega w (\phi^n+\frac{\Delta t}{(A^n)^2}(f'(\phi^n)+ \lambda B^n g'(\phi^n) u^n))dV-\int_\Omega
\frac{\Delta t}{(A^n)^2}
\nabla_\Gamma w \cdot F_1^n \, dV - \int_\Omega
\frac{\Delta t}{(A^n)^2}
\nabla_\Gamma w \cdot F_2^n \, dV.
$$
Weak form becomes:
$$
\int_\Omega w\,\phi^{n+1}\,dV=\int_\Omega [
w\,r_{\phi}+ \nabla_\Gamma w \cdot r_{\phi x}]
 dV.
$$
Residual form:
$$
r_\phi=\phi^n+\frac{\Delta t}{(A^n)^2}
[f'(\phi^n)+ \lambda B^n g'(\phi^n)\,u^n
]
$$
$$r_{\phi x}=-\frac{\Delta t}{(A^n)^2}
(F_1^n+ F_2^n
).
$$

**(Eq.2)**
$$
\int_\Omega w\,u^{n+1}\,dV=
\int_\Omega w\,u^n\,dV+
\int_\Omega w\,\Delta t\,\tilde{D}\,\nabla_\Gamma\!\cdot
\left(q(\phi^n)\nabla_\Gamma u^n\right)\,dV-
\int_\Omega w\,\Delta t\
\frac{L_{\mathrm{sat}}}{2}\,B^n(\phi^{n+1}-\phi^n)\,dV.
$$
Applying integration by parts to the diffusion term and assuming periodic
or no-flux boundary conditions,
$$
\int_\Omega w\,\nabla_\Gamma\cdot \mathbf{H}\,dV=
\int_{\partial \Omega} w\, (\mathbf{H} \cdot \mathbf{n}_{\Gamma})\, dS
-\int_\Omega \nabla_\Gamma w \cdot \mathbf{H}\,dV,
\qquad
\mathbf{H}=q(\phi^n)\nabla_\Gamma u^n.
$$
Substituting into the weak form gives
$$
\int_{\Omega} w\,u^{n+1}\,dV=\int_{\Omega}[w\,u^{n}-
\Delta t\,\tilde{D}\,\nabla_{\Gamma} w \cdot
( q(\phi^{n})\,\nabla_{\Gamma} u^{n})-
w\,
\Delta t\
\frac{L_{\mathrm{sat}}}{2}\,B^{n}\,
( \phi^{n+1} - \phi^{n}] dV.
$$

final weak form:
$$
\int_\Omega w\,u^{n+1}\,dV=
\int_\Omega \left(
w\,r_u+
\nabla_\Gamma w \cdot r_{ux}
\right)\,dV.
$$
Residual form:
$$
r_u=
u^n-
\frac{L_{\mathrm{sat}}}{2}\,B^n(\phi^{n+1}-\phi^n),
\qquad
r_{ux}=
-\Delta t\,\tilde{D}\,q(\phi^n)\nabla_\Gamma u^n.
$$














