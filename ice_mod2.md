# Implementation steps for the Demange ice-growth model

## 1. New app plumbing with only $u$ and $\phi$

### Goal

Start a new application with only two primary scalar fields and confirm PRISMS-PF runs and outputs them.

### What changed

In `loadVariableAttributes()`:

* Variable $0$ renamed to $u$
* Variable $1$ renamed to $\phi$
* Both set as `SCALAR` and `EXPLICIT_TIME_DEPENDENT`
* No physics added yet; this step only makes the app recognize the fields

---

## 2. First part of Eq. (2): diffusion term for $u$

### Goal

Implement only the diffusion-like part of the supersaturation equation:

* only the diffusion term,
* use the standard gradient $\nabla$ instead of $\nabla_\Gamma$,
* no coupling term yet.

### Physics added

$$
\partial_t u=\widetilde{D}\,\nabla\cdot\left(q(\phi)\nabla u\right)
$$

with

$$
q(\phi)=1-\phi.
$$

---

## 3. First part of Eq. (1): local source term for $\phi$ with anisotropy off

### Goal

Add the simplest local part of the phase-field equation:

* only the local driving/source term,
* anisotropy off, so $A(\mathbf{n})=1$ and $B(\mathbf{n})=1$,
* no divergence anisotropy terms yet.

### Physics added

$$
\partial_t \phi=f'(\phi)+\lambda g'(\phi)\,u.
$$

with typical choices

$$
f'(\phi)=-\phi+\phi^3
$$

and

$$
g'(\phi)=(1-\phi^2)^2.
$$

---

## 4. Add Laplacian smoothing term for $\phi$

### Physics added

$$
\nabla\cdot(\nabla\phi).
$$

In weak form, this contributes a flux proportional to

$$
-\nabla\phi.
$$

### Verification

* $\phi$ began to show smooth interface evolution
* this was the first step where the phase field behaved like a diffuse interface instead of only reacting locally

---

## 5. Add $u$--$\phi$ coupling back into $u$: latent heat term

### Physics added

$$
u^{n+1}=u^n-\frac{L_{\mathrm{sat}}}{2}\left(\phi^{n+1}-\phi^n\right).
$$

---

## 6. Add anisotropy scaffolding: define $\nabla_\Gamma$, compute $\mathbf{n}$, $A(\mathbf{n})$, and $B(\mathbf{n})$

### Definitions added

The anisotropic gradient operator:

$$
\nabla_\Gamma=(\partial_x,\partial_y,\Gamma\,\partial_z).
$$

The anisotropy factor $B(\mathbf{n})$:

$$
B(\mathbf{n})=\sqrt{n_x^2+n_y^2+\Gamma^2 n_z^2}.
$$

The anisotropy factor $A(\mathbf{n})$ in the paper form:

$$
A(\mathbf{n})=1+\epsilon_{xy}\cos(6\theta)+\epsilon_z\cos(2\psi).
$$

### What changed

Added helper code to compute:

* $\nabla_\Gamma$
* interface normal $\mathbf{n}$
* $A(\mathbf{n})$
* $B(\mathbf{n})$

Added parameters/constants in `customPDE.h` and `parameters.prm`:

* $\Gamma$
* $\epsilon_{xy}$
* $\epsilon_z$

---

## 7. Turn on $B(\mathbf{n})$ in the $\phi$ source term and apply $1/A(\mathbf{n})^2$ scaling

### Goal

Activate the paper's kinetic anisotropy in the local source term:

* multiply the coupling by $B(\mathbf{n})$,
* scale the time derivative term by $1/A(\mathbf{n})^2$.

### Physics added

Starting from

$$
A(\mathbf{n})^2\partial_t\phi=f'(\phi)+\lambda B(\mathbf{n})g'(\phi)\,u,
$$

solve for $\partial_t\phi$ to obtain

$$
\partial_t\phi=
\frac{f'(\phi)+\lambda B(\mathbf{n})g'(\phi)\,u}{A(\mathbf{n})^2}.
$$

---

## 8. Implement $\Gamma$-consistent weak-form fluxes and run 3D tests

### Goal

Ensure that $\nabla_\Gamma$ is implemented correctly in weak form.

For an operator of the form

$$
\nabla_\Gamma\cdot\left(q\nabla_\Gamma u\right),
$$

the weak form implies a $\Gamma^2$ scaling in the $z$ component of the flux.

### Flux form

In 3D, the correct flux behaves as

$$
q(\phi)\left(u_x,u_y,\Gamma^2 u_z\right).
$$

### Verification

* ran 3D nucleation tests with two seeds
* in ParaView, $\phi$ isosurfaces showed two nuclei as expected
* the $u$ field smoothed and spread as expected for diffusion

---

## 9. Fix the paper normal definition

### Goal

Match the paper's definition of the interface normal used inside $A(\mathbf{n})$ and $B(\mathbf{n})$.

### Physics corrected

$$
\mathbf{n}=-\frac{\nabla\phi}{|\nabla\phi|}.
$$

This replaced using $\nabla_\Gamma\phi$ inside the normal definition.

---

## 10. Turn on full Eq. (2) coupling anisotropy in $u$

### Goal

Match the orientation dependence in Eq. (2) of the paper.

### Physics updated

From

$$
u^{n+1}=u^n-\frac{L_{\mathrm{sat}}}{2}\left(\phi^{n+1}-\phi^n\right)
$$

to

$$
u^{n+1}=u^n-\frac{L_{\mathrm{sat}}}{2}B(\mathbf{n})\left(\phi^{n+1}-\phi^n\right).
$$

---

## 11. Rewrite the $\phi$ Laplacian/$F_2$ term in the paper's explicit-update form

### Paper structure used

$$
\phi^{n+1}=\phi^n+\frac{\Delta t}{A(\mathbf{n})^2}
\left[
\text{local}
+\nabla_\Gamma\cdot(F_1+F_2)
\right].
$$

with

$$
F_2=A(\mathbf{n})^2\nabla_\Gamma\phi.
$$

This puts the $F_2$ contribution into the same explicit-update layout as the paper.

---

## 12. Add the hard-anisotropy or stiffness term $F_1$

### Goal

Add the missing anisotropic stiffness-like divergence term responsible for stronger faceting and branching behavior.

### Physics added

$$
\frac{1}{2}\nabla_\Gamma\cdot\left(
|\nabla\phi|^2
\frac{\partial\left(A(\mathbf{n})^2\right)}{\partial(\nabla\phi)}
\right).
$$

This is the $F_1$-type contribution from the Demange formulation.
