// =============================================================================================
// equations.cc  (ICE APP: u, phi, xi1)
// =============================================================================================

void
customAttributeLoader::loadVariableAttributes()
{
  // Variable 0: u
  set_variable_name(0, "u");
  set_variable_type(0, SCALAR);
  set_variable_equation_type(0, EXPLICIT_TIME_DEPENDENT);
  set_dependencies_value_term_RHS(0, "u,phi,xi1");
  set_dependencies_gradient_term_RHS(0, "grad(u),grad(phi)");

  // Variable 1: phi
  set_variable_name(1, "phi");
  set_variable_type(1, SCALAR);
  set_variable_equation_type(1, EXPLICIT_TIME_DEPENDENT);
  set_dependencies_value_term_RHS(1, "phi,xi1");
  set_dependencies_gradient_term_RHS(1, "grad(phi)");

  // Variable 2: xi1
  set_variable_name(2, "xi1");
  set_variable_type(2, SCALAR);
  set_variable_equation_type(2, EXPLICIT_TIME_DEPENDENT);
  set_dependencies_value_term_RHS(2, "phi,u");
  set_dependencies_gradient_term_RHS(2, "grad(phi)");
}

template <int dim, int degree>
void
customPDE<dim, degree>::explicitEquationRHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double>                             element_volume) const
{
  const scalarvalueType u   = variable_list.get_scalar_value(0);
  const scalargradType  ux  = variable_list.get_scalar_gradient(0);
  const scalarvalueType phi = variable_list.get_scalar_value(1);
  const scalargradType  phix = variable_list.get_scalar_gradient(1);
  const scalarvalueType xi1 = variable_list.get_scalar_value(2);

  auto apply_Gamma = [&](const scalargradType &g) -> scalargradType {
    scalargradType out = g;
    if constexpr (dim == 3)
      out[2] = g[2] * constV(Gamma);
    return out;
  };

  auto apply_GammaT = [&](const scalargradType &g) -> scalargradType {
    scalargradType out = g;
    if constexpr (dim == 3)
      out[2] = g[2] * constV(Gamma);
    return out;
  };

  auto compute_normal = [&](const scalargradType &grad_phi) -> scalargradType {
    scalarvalueType mag2 = constV(0.0);
    for (unsigned int d = 0; d < dim; ++d)
      mag2 = mag2 + grad_phi[d] * grad_phi[d];

    const scalarvalueType eps = constV(1e-12);
    scalarvalueType inv_mag;
    for (unsigned int v = 0; v < mag2.size(); ++v)
      inv_mag[v] = 1.0 / std::sqrt(mag2[v] + eps[v]);

    scalargradType n;
    for (unsigned int d = 0; d < dim; ++d)
      n[d] = -grad_phi[d] * inv_mag;

    return n;
  };

  auto compute_B = [&](const scalargradType &n) -> scalarvalueType {
    scalarvalueType B2 = n[0] * n[0] + n[1] * n[1];
    if constexpr (dim == 3)
      B2 = B2 + constV(Gamma * Gamma) * n[2] * n[2];

    scalarvalueType B;
    for (unsigned int v = 0; v < B2.size(); ++v)
      B[v] = std::sqrt(B2[v]);

    return B;
  };

  auto compute_A = [&](const scalargradType &n) -> scalarvalueType {
    scalarvalueType A = constV(1.0);

    scalarvalueType theta;
    for (unsigned int v = 0; v < theta.size(); ++v)
      theta[v] = std::atan2(n[1][v], n[0][v]);

    scalarvalueType cos6t;
    for (unsigned int v = 0; v < cos6t.size(); ++v)
      cos6t[v] = std::cos(6.0 * theta[v]);

    A = A + constV(eps_xy) * cos6t;

    if constexpr (dim == 3)
      {
        scalarvalueType r;
        for (unsigned int v = 0; v < r.size(); ++v)
          r[v] = std::sqrt(n[0][v] * n[0][v] + n[1][v] * n[1][v]);

        scalarvalueType psi;
        for (unsigned int v = 0; v < psi.size(); ++v)
          psi[v] = std::atan2(r[v], n[2][v]);

        scalarvalueType cos2p;
        for (unsigned int v = 0; v < cos2p.size(); ++v)
          cos2p[v] = std::cos(2.0 * psi[v]);

        A = A + constV(eps_z) * cos2p;
      }

    return A;
  };

  const scalargradType nvec = compute_normal(phix);
  const scalarvalueType Bn = compute_B(nvec);
  const scalarvalueType An = compute_A(nvec);
  const scalarvalueType An2 = An * An;

  scalarvalueType inv_An2;
  for (unsigned int v = 0; v < An2.size(); ++v)
    inv_An2[v] = 1.0 / (An2[v] + 1e-12);

  const scalarvalueType fprime = -phi + phi * phi * phi;
  const scalarvalueType one_minus_phi2 = constV(1.0) - phi * phi;
  const scalarvalueType gprime = one_minus_phi2 * one_minus_phi2;

  // f1(phi, grad(phi), u) = -f'(phi) + lambda * B(n) * g'(phi) * u
  const scalarvalueType f1 = (-fprime) + constV(lambda) * Bn * gprime * u;

  auto A2_from_grad = [&](const scalargradType &g) -> scalarvalueType {
    const scalargradType n_tmp = compute_normal(g);
    const scalarvalueType A_tmp = compute_A(n_tmp);
    return A_tmp * A_tmp;
  };

  // d(A^2)/d(grad(phi)) via central differences
  const double delta = 1e-6;
  scalargradType dA2_dg;
  for (unsigned int d = 0; d < dim; ++d)
    {
      scalargradType gp = phix;
      scalargradType gm = phix;
      gp[d] = gp[d] + constV(delta);
      gm[d] = gm[d] - constV(delta);

      const scalarvalueType A2p = A2_from_grad(gp);
      const scalarvalueType A2m = A2_from_grad(gm);
      dA2_dg[d] = (A2p - A2m) * constV(1.0 / (2.0 * delta));
    }

  scalarvalueType gradphi_mag2 = constV(0.0);
  for (unsigned int d = 0; d < dim; ++d)
    gradphi_mag2 = gradphi_mag2 + phix[d] * phix[d];

  const scalargradType gamma_grad_phi = apply_Gamma(phix);

  scalargradType F1_inner;
  for (unsigned int d = 0; d < dim; ++d)
    F1_inner[d] = gradphi_mag2 * dA2_dg[d] + An2 * gamma_grad_phi[d];

  scalargradType F1 = apply_GammaT(F1_inner);
  for (unsigned int d = 0; d < dim; ++d)
    F1[d] = F1[d] * constV(0.5);

  // F2(phi, grad(u)) = D_tilde * Gamma^T * q(phi) * Gamma * grad(u)
  const scalarvalueType q = constV(1.0) - phi;
  const scalargradType gamma_grad_u = apply_Gamma(ux);
  scalargradType F2;
  for (unsigned int d = 0; d < dim; ++d)
    F2[d] = gamma_grad_u[d] * (constV(Dtilde) * q);
  F2 = apply_GammaT(F2);

  // Weak-form residual targets requested in prompt
  const scalarvalueType dt = constV(userInputs.dtValue);

  const scalarvalueType r_xi = f1;
  scalargradType r_xi_x;
  for (unsigned int d = 0; d < dim; ++d)
    r_xi_x[d] = -F1[d];

  const scalarvalueType r_phi = phi + dt * inv_An2 * xi1;

  const scalarvalueType r_u =
    u - dt * constV(0.5 * Lsat) * (Bn * inv_An2) * xi1;

  scalargradType r_u_x;
  for (unsigned int d = 0; d < dim; ++d)
    r_u_x[d] = -dt * F2[d];

  variable_list.set_scalar_value_term_RHS(0, r_u);
  variable_list.set_scalar_gradient_term_RHS(0, r_u_x);

  variable_list.set_scalar_value_term_RHS(1, r_phi);

  variable_list.set_scalar_value_term_RHS(2, r_xi);
  variable_list.set_scalar_gradient_term_RHS(2, r_xi_x);
}

template <int dim, int degree>
void
customPDE<dim, degree>::nonExplicitEquationRHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double>                             element_volume) const
{}

template <int dim, int degree>
void
customPDE<dim, degree>::equationLHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double>                             element_volume) const
{}
