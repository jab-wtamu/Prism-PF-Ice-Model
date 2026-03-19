// =============================================================================================
// equations.cc  (ICE APP: u, phi)
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
set_dependencies_gradient_term_RHS(0, "grad(u),grad(phi)");
  // ---------------------------------------------------------------------------
  // Variable 1: phi (phase field)
  // ---------------------------------------------------------------------------
  set_variable_name(1, "phi");
  set_variable_type(1, SCALAR);
  set_variable_equation_type(1, EXPLICIT_TIME_DEPENDENT);

  // STEP 2 deps 
  // phi equation now depends on phi and u values (still no gradients for this step)
  set_dependencies_value_term_RHS(1, "phi,u");
  
// STEP 3 : phi now includes the normal Laplacian term ∇·(∇phi),
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
  
  // STEP 3 : normal gradient of phi, needed for ∇·(∇phi) term
scalargradType phix = variable_list.get_scalar_gradient(1); //  gradient ∇phi




 // ===========================================================================================
 // GRADIENT GAMMA VARIABLE (if used)
 // ===========================================================================================
 
// grad_gamma  (David's Motiel's suggestion)
// In 2D: returns grad unchanged (no z component)
  
  
 [[maybe_unused]] auto grad_gamma = [&](const scalargradType &g, const double Gamma_val) -> scalargradType {
    scalargradType out = g;  //gradient vector
    if constexpr (dim == 3)
      {
        out[2] = g[2] * constV(Gamma_val);      //scale the z component of the gradient by Γ
      }
    return out;
  };
  
  // STEP 6 NEW: Γ-gradient of phi (paper uses ∇_Γ phi in n, A(n), B(n))
[[maybe_unused]] const scalargradType phix_gamma = grad_gamma(phix, Gamma);
  
  
// ===========================================================================================
  // Compute unit normal vector of phi n = -∇phi / |∇phi|  
  // Uses eps to avoid division by zero.
// ===========================================================================================  


  auto compute_normal = [&](const scalargradType &grad_phi) -> scalargradType {
    scalarvalueType mag2 = constV(0.0);
    for (unsigned int d = 0; d < dim; ++d)
      mag2 = mag2 + grad_phi[d] * grad_phi[d];

    // small epsilon (vectorized) to avoid division by zero
    const scalarvalueType eps = constV(1e-12);
    scalarvalueType inv_mag;
    for (unsigned int v = 0; v < mag2.size(); ++v)
      inv_mag[v] = 1.0 / std::sqrt(mag2[v] + eps[v]);           //∣∇ϕ∣≈∣∇ϕ∣2+ϵ

    scalargradType n;
    for (unsigned int d = 0; d < dim; ++d)
      n[d] = -grad_phi[d] * inv_mag;            //−∇ϕ/∣∇ϕ∣

    return n;
  };
  
  
// =================================================================================================
  // B(n) = sqrt(nx^2 + ny^2 + Gamma^2 nz^2)  BðnÞ¼ n2 x þ n2 y þ Γ2n2
// ===============================================================================================


  auto compute_B = [&](const scalargradType &n, const double Gamma_val) -> scalarvalueType {
    scalarvalueType B2 = n[0] * n[0] + n[1] * n[1];
    if constexpr (dim == 3)
      {
        B2 = B2 + constV(Gamma_val * Gamma_val) * n[2] * n[2];
      }

    scalarvalueType B;
    for (unsigned int v = 0; v < B2.size(); ++v)
      B[v] = std::sqrt(B2[v]);                                  //B2=nx2​+ny2​+Γ2nz2

    return B;
  };
// ==============================================================================================
  // A(n) = 1 + eps_xy cos(6θ) + eps_z cos(2ψ)
  // θ=atan2(ny,nx)
  // ψ=atan2(sqrt(nx^2+ny^2), nz)   (3D only)
  //
  // In 2D: ψ not defined; we implement A using only the 6-fold in-plane term:
  //   A = 1 + eps_xy cos(6θ)
// ==============================================================================================  


  auto compute_A = [&](const scalargradType &n,
                       const double eps_xy_val,
                       const double eps_z_val) -> scalarvalueType {
    scalarvalueType A = constV(1.0);

    // theta = atan2(ny, nx)
    scalarvalueType theta;      //n[1][v] = ny ,        //n[0][v] = nx 
    for (unsigned int v = 0; v < theta.size(); ++v)
      theta[v] = std::atan2(n[1][v], n[0][v]);

    // 6-fold term
    scalarvalueType cos6t;
    for (unsigned int v = 0; v < cos6t.size(); ++v)
      cos6t[v] = std::cos(6.0 * theta[v]);

    A = A + constV(eps_xy_val) * cos6t;

    if constexpr (dim == 3)
      {
        // psi = atan2( sqrt(nx^2+ny^2), nz )
        scalarvalueType r;
        for (unsigned int v = 0; v < r.size(); ++v)
          r[v] = std::sqrt(n[0][v] * n[0][v] + n[1][v] * n[1][v]);

        scalarvalueType psi;
        for (unsigned int v = 0; v < psi.size(); ++v)
          psi[v] = std::atan2(r[v], n[2][v]);

        scalarvalueType cos2p;
        for (unsigned int v = 0; v < cos2p.size(); ++v)
          cos2p[v] = std::cos(2.0 * psi[v]);

        A = A + constV(eps_z_val) * cos2p;
      }

    return A;
  };
// STEP 6
  // Compute n, A(n), B(n) using Γ-gradient (paper-consistent)
const scalargradType nvec = compute_normal(phix);
const scalarvalueType Bn  = compute_B(nvec, Gamma);
[[maybe_unused]] const scalarvalueType An = compute_A(nvec, eps_xy, eps_z);

 

 // ===========================================================================================
// STEP 1 : u diffusion only, with Γ-anisotropy in 3D via flux tensor KΓ = diag(1,1,Γ^2)
//
//   ∂t u = Dtilde * ∇ · ( q(phi) KΓ ∇u )
//   q(phi) = 1 - phi
//
// Weak-form gradient residual / flux term:
//   r_ux = -dt * Dtilde * q(phi) * (KΓ ∇u)
// ===========================================================================================

// Cutoff function q(phi) = 1 - phi
scalarvalueType q = constV(1.0) - phi;

// scale = -(dt * Dtilde) * q(phi)
scalarvalueType scale_u = constV(-userInputs.dtValue * Dtilde) * q;

// Build flux eqx_u = scale_u * (KΓ ∇u)
scalargradType eqx_u;
for (unsigned int d = 0; d < dim; ++d)
  eqx_u[d] = ux[d] * scale_u;

// Only in 3D: z-flux gets Γ^2
if constexpr (dim == 3)
  eqx_u[2] = ux[2] * (scale_u * constV(Gamma * Gamma));

  // ===========================================================================================


  scalarvalueType eq_u = u;  // placeholder for now, waiting for phi update

  
  // ===========================================================================================
  // STEP 2: 
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
  
  // ===========================================================================================

                        // STEP 6

  // ===========================================================================================
//Local RHS for phi (paper form): includes B(n)
// CHANGE: multiply coupling by Bn to activate kinetic anisotropy
scalarvalueType rhs_phi = (-fprime) + constV(lambda)*Bn*gprime*u;

// STEP NEXT: paper form is  A(n)^2 * dphi/dt = rhs_phi
//          => dphi/dt = rhs_phi / A(n)^2
// We apply this by scaling the explicit update RHS by 1/(A^2).
scalarvalueType An2 = (An * An);

// Avoid divide-by-zero (shouldn't happen for small eps, but safe)
scalarvalueType inv_An2;
for (unsigned int v = 0; v < An2.size(); ++v)
  inv_An2[v] = 1.0 / (An2[v] + 1e-12);

// Scaled RHS used in explicit update
scalarvalueType rhs_phi_scaled = rhs_phi * inv_An2;

  // Explicit value update for phi (scalar coefficient for w)
scalarvalueType eq_phi = phi + constV(userInputs.dtValue) * rhs_phi_scaled;

// ===========================================================================================
// STEP 4 (paper structure for F2 term):
//
// Paper update: phi^{n+1} = phi^n + (dt / A(n)^2) * [ ... + ∇_Γ · (F2) + ... ]
// where F2 = A(n)^2 * (KΓ ∇phi)
//
// ===========================================================================================

// Build F2 = (1/2) * A(n)^2 * (KΓ ∇phi)
scalargradType F2;
for (unsigned int d = 0; d < dim; ++d)
  F2[d] = phix[d] * (An2 * constV(0.5));

if constexpr (dim == 3)
  F2[2] = phix[2] * (An2 * constV(0.5 * Gamma * Gamma));

// Now apply the global (dt / A^2) scaling to the divergence term:
// eqx_phi = -(dt/A^2) * F2
scalargradType eqx_phi;
for (unsigned int d = 0; d < dim; ++d)
  eqx_phi[d] = F2[d] * (constV(-userInputs.dtValue) * inv_An2);
  
// ===========================================================================================
// STEP 8: NEW (stiffness / hard anisotropy term)
//
// Add:  + (1/2) ∇_Γ · ( |∇_Γ φ|^2 * ∂(A(n)^2)/∂(∇_Γ φ) )
//
// In PRISMS flux form (since gradient term is -∫ ∇w · flux dV):
//   eqx_phi +=  -dt * (1/2) * ( |g|^2 * dA2_dg )
//
// Here we implement dA2_dg numerically by central differences:
//   d(A^2)/dg_i ≈ (A^2(g+δ e_i) - A^2(g-δ e_i)) / (2δ)
//
// where:
//   g = ∇_Γ φ  (we use phix_gamma)
// ===========================================================================================

// Toggle on/off while debugging
const bool enable_stiffness = true;

if (enable_stiffness)
{
  // ---- helper: compute A^2 given a gradient g (g treated as ∇_Γ φ) ----
  auto A2_from_grad = [&](const scalargradType &g) -> scalarvalueType {
    const scalargradType n_tmp  = compute_normal(g);                 // n = -g/|g|
    const scalarvalueType A_tmp = compute_A(n_tmp, eps_xy, eps_z);    // A(n)
    return A_tmp * A_tmp;                                            // A^2
  };

  const scalargradType g0 = phix;                              // g = ∇ φ

  scalarvalueType g0_mag2 = constV(0.0);
  for (unsigned int d = 0; d < dim; ++d)
    g0_mag2 = g0_mag2 + g0[d] * g0[d];

  // Skip stiffness away from the interface (reduces noise + cost)
  const double g2_threshold = 1e-10;
  bool any_lane_active = false;
  for (unsigned int v = 0; v < g0_mag2.size(); ++v)
    if (g0_mag2[v] > g2_threshold) { any_lane_active = true; break; }

  if (!any_lane_active)
  {
    // do nothing: keep eqx_phi as-is
  }
  else
  {
    // Numerical derivative step (tune if noisy)
    const double delta = 1e-6;

    // d(A^2)/dg as a vector
    scalargradType dA2_dg;

    for (unsigned int d = 0; d < dim; ++d)
    {
      scalargradType gp = g0;
      scalargradType gm = g0;

      gp[d] = gp[d] + constV(delta);
      gm[d] = gm[d] - constV(delta);

      const scalarvalueType A2p = A2_from_grad(gp);
      const scalarvalueType A2m = A2_from_grad(gm);

      // central difference: (A2(g+δ) - A2(g-δ)) / (2δ)
      dA2_dg[d] = (A2p - A2m) * constV(1.0 / (2.0 * delta));
    }

    // Build stiffness flux:
//   F_stiff = -dt * 0.5 * ( |g|^2 * dA2_dg )
scalargradType F_stiff;
for (unsigned int d = 0; d < dim; ++d)
{
  F_stiff[d] = dA2_dg[d] * (g0_mag2 * constV(-0.5 * userInputs.dtValue));
}

// Apply KΓ to stiffness flux the same way as F2/u-diffusion (ONLY if g0 = phix)
if constexpr (dim == 3)
{
  F_stiff[2] = F_stiff[2] * constV(Gamma * Gamma);
}

// Add to existing phi gradient flux (NO extra dt; dt is already inside F_stiff)
for (unsigned int d = 0; d < dim; ++d)
{
  eqx_phi[d] = eqx_phi[d] + F_stiff[d] * inv_An2;
}
} // closes: else (any_lane_active)
}   // closes: if (enable_stiffness)

// ===========================================================================================

// STEP (u–phi coupling): use the SAME phi increment we actually computed here (local part)
// NOTE: PRISMS adds divergence terms via eqx_phi later, so eq_phi is not full phi^{n+1}.
scalarvalueType delta_phi_local = constV(userInputs.dtValue) * rhs_phi_scaled;  // dt * dphi/dt_local

eq_u = u - constV(0.5 * Lsat) * Bn * delta_phi_local;

// ===========================================================================================
//lambda

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
