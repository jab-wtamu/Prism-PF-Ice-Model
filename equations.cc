// =============================================================================================
// equations.cc  (ICE APP: u, phi)
// =============================================================================================
//
// STATUS / INCREMENTAL PLAN
//
// This week’s incremental updates (STEP 3):
//
// 1) u equation (index 0):
//    - Kept diffusion term with cutoff:
//        ∂t u = Dtilde * ∇·( (1 - phi) ∇u )
//    - Added the Demange “minus/coupling” term using the explicit phi update:
//        u^{n+1} = u^n - (Lsat/2) * (phi^{n+1} - phi^n)
//
// 2) phi equation (index 1):
//    - Kept local source term (anisotropy OFF):
//        ∂t phi = f'(phi) + lambda * g'(phi) * u
//    - Added spatial smoothing term (normal Laplacian):
//        + ∇·(∇phi)  
//
// =============================================================================================


// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================

void
customAttributeLoader::loadVariableAttributes()
{
  // ---------------------------------------------------------------------------
  // Variable 0: u (reduced supersaturation)
  // ---------------------------------------------------------------------------
  set_variable_name(0, "u");
  set_variable_type(0, SCALAR);
  set_variable_equation_type(0, EXPLICIT_TIME_DEPENDENT);

  // STEP 1 deps (kept): u equation uses u, phi, and grad(u)
  set_dependencies_value_term_RHS(0, "u,phi");
  set_dependencies_gradient_term_RHS(0, "grad(u)");

  // ---------------------------------------------------------------------------
  // Variable 1: phi (phase field)
  // ---------------------------------------------------------------------------
  set_variable_name(1, "phi");
  set_variable_type(1, SCALAR);
  set_variable_equation_type(1, EXPLICIT_TIME_DEPENDENT);

  // STEP 2 deps (UPDATED):
  // phi equation now depends on phi and u values (still no gradients for this step)
  set_dependencies_value_term_RHS(1, "phi,u");
  
// STEP 3 (NEW): phi now includes the normal Laplacian term ∇·(∇phi),
// which appears in weak form as a gradient-test term
set_dependencies_gradient_term_RHS(1, "grad(phi)");
}


// =============================================================================================
// explicitEquationRHS
// =============================================================================================

template <int dim, int degree>
void
customPDE<dim, degree>::explicitEquationRHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{
  // ---------------------------------------------------------------------------
  // Get the values/gradients of the primary fields
  // ---------------------------------------------------------------------------

  // u (variable 0)
  scalarvalueType u  = variable_list.get_scalar_value(0);
  scalargradType  ux = variable_list.get_scalar_gradient(0); //  gradient ∇u

  // phi (variable 1)
  scalarvalueType phi = variable_list.get_scalar_value(1);
  
  // STEP 3 (NEW): normal gradient of phi, needed for ∇·(∇phi) term
scalargradType phix = variable_list.get_scalar_gradient(1); //  gradient ∇phi



  // ---------------------------------------------------------------------------
  // Common: build a ZERO gradient tensor (needed for any "no-gradient" terms)
  // ---------------------------------------------------------------------------
  scalargradType zero_grad;
  for (unsigned int d = 0; d < dim; ++d)
    zero_grad[d] = constV(0.0);

  // ===========================================================================================
  // STEP 1 (kept): u diffusion only
  //
  //   ∂t u = Dtilde * ∇ · ( q(phi) ∇u )
  //   q(phi) = 1 - phi
  //
  // Weak-form gradient residual / flux term:
  //   r_ux = -dt * Dtilde * q(phi) * ∇u
  // ===========================================================================================

  // Cutoff function q(phi) = 1 - phi
  scalarvalueType q = constV(1.0) - phi;

  // scale = -(dt * Dtilde) * q(phi)
  scalarvalueType scale_u = constV(-userInputs.dtValue * Dtilde) * q;

  scalargradType eqx_u;
  for (unsigned int d = 0; d < dim; ++d){
    eqx_u[d] = ux[d] * scale_u;
    }

  scalarvalueType eq_u = u;  // placeholder for now, waiting for phi update

  
  // ===========================================================================================
  // STEP 2 (NEW): phi local source term only (anisotropy OFF: A=1, B=1)
  //
  // Target partial equation:
  //   A(n)^2 * ∂t phi = f'(phi) + lambda * B(n) * g'(phi) * u
  //
  // With anisotropy OFF:
  //   A = 1, B = 1  =>  ∂t phi = f'(phi) + lambda * g'(phi) * u
  //
  // And since we are NOT implementing divergence terms (F1/F2) yet,
  // phi has NO gradient RHS term in this step.
  //
  // Explicit update:
  //   phi^{n+1} = phi^n + dt * ( f'(phi^n) + lambda * g'(phi^n) * u^n )
  // ===========================================================================================

  // f'(phi) for f(phi) = -phi^2/2 + phi^4/4  =>  f'(phi) = -phi + phi^3
  scalarvalueType fprime = -phi + phi * phi * phi;

  // g'(phi) used by Demange: g'(phi) = (1 - phi^2)^2
  scalarvalueType one_minus_phi2 = constV(1.0) - phi * phi;
  scalarvalueType gprime         = one_minus_phi2 * one_minus_phi2;

  // Local RHS for phi (anisotropy OFF: B=1)
  scalarvalueType rhs_phi = fprime + constV(lambda) * gprime * u;

  // Explicit value update for phi (scalar coefficient for w)
scalarvalueType eq_phi = phi + constV(userInputs.dtValue) * rhs_phi;

  // ===========================================================================================

// STEP 3 (NEW): Add the  Laplacian term ∇·(∇phi) in weak form.
// Weak form contribution:  -∫ dt * ∇w · ∇phi dV
// So the gradient-term coefficient is:  r_phix = -dt * ∇phi (gradient term coeff)


scalargradType eqx_phi;                         //(vector coefficient for grad(w))
for (unsigned int d = 0; d < dim; ++d){
  eqx_phi[d] = phix[d] * constV(-userInputs.dtValue);
}

// STEP 3 (NEW): u–phi coupling (Demange Eq.2 minus term, anisotropy OFF) ---
scalarvalueType delta_phi = eq_phi - phi;               // eq_phi(new phi at n+1) // phi (old)
eq_u = u - constV(0.5 * Lsat) * delta_phi;


  // ---------------------------------------------------------------------------
  // Submit RHS terms
  // ---------------------------------------------------------------------------

  // u
  variable_list.set_scalar_value_term_RHS(0, eq_u);
  variable_list.set_scalar_gradient_term_RHS(0, eqx_u);

  // phi
  variable_list.set_scalar_value_term_RHS(1, eq_phi);
  variable_list.set_scalar_gradient_term_RHS(1, eqx_phi);
}


// =============================================================================================
// nonExplicitEquationRHS
// =============================================================================================

template <int dim, int degree>
void
customPDE<dim, degree>::nonExplicitEquationRHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{
  // No non-explicit equations in this incremental step.
}


// =============================================================================================
// equationLHS
// =============================================================================================

template <int dim, int degree>
void
customPDE<dim, degree>::equationLHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{
  // No time-independent equations in this incremental step.
}
