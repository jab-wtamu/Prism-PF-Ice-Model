// Demange snow-crystal growth model -- PRISMS-PF 4.0.0, 3D.
//
// The modular phase-field formulation provides configurable conventions:
//   half_convention: variational | demange_literal
//   q_normalization: one_sided | demange_literal
//   gamma_structure: demange_literal | stretched_variational
//   anisotropy: pole_regular | supplement_literal
//   pole_weight: sin6 | sin4 | sin8 | quintic
//
// Required Model constants are validated at initialization. Missing and unused
// constants are rejected to keep configurations reproducible.

#pragma once

#include "demange_anisotropy.h"

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/pde_operator_base.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <initializer_list>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class CustomPDE : public PDEOperatorBase<dim, degree, number>
{
  static_assert(dim == 3, "This snow crystal application is 3D only.");

public:
  using ScalarValue = dealii::VectorizedArray<number>;
  using ScalarGrad  = dealii::Tensor<1, dim, ScalarValue>;
  using ScalarHess  = dealii::Tensor<2, dim, ScalarValue>;
  using VectorValue = dealii::Tensor<1, dim, ScalarValue>;
  using VectorGrad  = dealii::Tensor<2, dim, ScalarValue>;
  using VectorHess  = dealii::Tensor<3, dim, ScalarValue>;
  using PDEOperatorBase<dim, degree, number>::get_user_inputs;
  using PDEOperatorBase<dim, degree, number>::get_pf_tools;

  /**
   * @brief Constructor. Reads and validates every Model constant, then prints a model summary.
   */
  CustomPDE(const UserInputParameters<dim> &_user_inputs, PhaseFieldTools<dim> &_pf_tools)
    : PDEOperatorBase<dim, degree, number>(_user_inputs, _pf_tools)
    , u0(read_double(_user_inputs, "u0"))
    , eps_xy(read_double(_user_inputs, "eps_xy"))
    , eps_z(read_double(_user_inputs, "eps_z"))
    , Gamma(read_double(_user_inputs, "Gamma"))
    , lambda(read_double(_user_inputs, "lambda"))
    , D_tilde(read_double(_user_inputs, "D_tilde"))
    , Lsat(read_double(_user_inputs, "Lsat"))
    , half_convention(
        read_choice(_user_inputs, "half_convention", {"variational", "demange_literal"}))
    , q_normalization(
        read_choice(_user_inputs, "q_normalization", {"one_sided", "demange_literal"}))
    , gamma_structure(read_choice(_user_inputs,
                                  "gamma_structure",
                                  {"demange_literal", "stretched_variational"}))
    , anisotropy_model(
        read_choice(_user_inputs, "anisotropy", {"pole_regular", "supplement_literal"}))
    , pole_weight(
        read_choice(_user_inputs, "pole_weight", {"sin6", "sin4", "sin8", "quintic"}))
    , seed(read_choice(_user_inputs, "seed", {"hollow_prism", "disc"}))
    , seed_hex_orientation(seed == "hollow_prism" ? read_choice(_user_inputs,
                                                                "seed_hex_orientation",
                                                                {"model_consistent", "pf3"}) :
                                                    std::string("none"))
    , seed_center {{read_double(_user_inputs, "seed_center_x"),
                    read_double(_user_inputs, "seed_center_y"),
                    read_double(_user_inputs, "seed_center_z")}}
    , ic_interface_width(read_double(_user_inputs, "ic_interface_width"))
    , prism_outer_radius(
        read_seed_double(_user_inputs, seed == "hollow_prism", "prism_outer_radius"))
    , prism_outer_half_height(
        read_seed_double(_user_inputs, seed == "hollow_prism", "prism_outer_half_height"))
    , prism_inner_radius(
        read_seed_double(_user_inputs, seed == "hollow_prism", "prism_inner_radius"))
    , prism_inner_half_height(
        read_seed_double(_user_inputs, seed == "hollow_prism", "prism_inner_half_height"))
    , disc_radius(read_seed_double(_user_inputs, seed == "disc", "disc_radius"))
    , disc_half_thickness(read_seed_double(_user_inputs, seed == "disc", "disc_half_thickness"))
    , aniso(anisotropy_model == "pole_regular" ? demange::Model::pole_regular :
                                                 demange::Model::supplement_literal,
            static_cast<double>(eps_xy),
            static_cast<double>(eps_z),
            demange::parse_pole_weight(pole_weight))
  {
    check_model_constant_names(_user_inputs);
    flux_options.variational_half = (half_convention == "variational");
    flux_options.stretched        = (gamma_structure == "stretched_variational");
    flux_options.Gamma            = static_cast<double>(Gamma);
    // Literal [S] is singular at n = +-z; reproduce PF3's angle floor there (not a physical choice).
    flux_options.pf3_pole_floor = (anisotropy_model == "supplement_literal");
    // Section Q: q(phi) = q_factor (1 - phi); q_factor = 1/2 makes D_tilde the vapour diffusivity.
    q_factor       = (q_normalization == "one_sided") ? 0.5 : 1.0;
    seed_is_disc   = (seed == "disc");
    hex_swap_xy    = (seed_hex_orientation == "model_consistent");
    validate_and_report();
  }

private:
  // ================================================================================================
  // Parameter reading
  // ================================================================================================
  static double
  read_double(const UserInputParameters<dim> &ui, const std::string &name)
  {
    try
      {
        return ui.user_constants.get_double(name);
      }
    catch (const std::exception &e)
      {
        throw std::runtime_error("parameters: Model constant '" + name +
                                 "' (DOUBLE) is missing or has the wrong type (" + e.what() + ")");
      }
  }

  static std::string
  read_choice(const UserInputParameters<dim>     &ui,
              const std::string                  &name,
              std::initializer_list<const char *> allowed)
  {
    std::string value;
    try
      {
        value = ui.user_constants.get_string(name);
      }
    catch (const std::exception &e)
      {
        throw std::runtime_error("parameters: Model constant '" + name +
                                 "' (STRING) is missing or has the wrong type (" + e.what() + ")");
      }
    for (const char *a : allowed)
      if (value == a)
        return value;
    std::string msg = "parameters: Model constant '" + name + "' = '" + value + "' must be one of:";
    for (const char *a : allowed)
      msg += std::string(" ") + a;
    throw std::runtime_error(msg);
  }

  // Seed geometry is read only for the selected seed; the other seed's constants must be absent.
  static double
  read_seed_double(const UserInputParameters<dim> &ui, const bool used, const std::string &name)
  {
    return used ? read_double(ui, name) : 0.0;
  }

  // A parameters file must list exactly the Model constants this model reads: missing ones are
  // reported by read_double/read_choice, every other name is rejected here. PRISMS-PF itself accepts
  // any Model constant, so without this check a constant of the removed pf3_frozen code path
  // (implementation, regval) or of the other seed would be silently ignored. The names come from the
  // parser PRISMS-PF uses to declare the constants, so no absent constant is ever looked up.
  void
  check_model_constant_names(const UserInputParameters<dim> &ui) const
  {
    std::set<std::string> used = {"u0",
                                  "eps_xy",
                                  "eps_z",
                                  "Gamma",
                                  "lambda",
                                  "D_tilde",
                                  "Lsat",
                                  "half_convention",
                                  "q_normalization",
                                  "gamma_structure",
                                  "anisotropy",
                                  "pole_weight",
                                  "seed",
                                  "seed_center_x",
                                  "seed_center_y",
                                  "seed_center_z",
                                  "ic_interface_width"};
    if (seed == "hollow_prism")
      used.insert({"seed_hex_orientation",
                   "prism_outer_radius",
                   "prism_outer_half_height",
                   "prism_inner_radius",
                   "prism_inner_half_height"});
    else
      used.insert({"disc_radius", "disc_half_thickness"});
    for (const std::string &name : UserConstants<dim>::get_names(ui.user_constants.file_name))
      if (used.count(name) == 0)
        throw std::runtime_error(
          "parameters: Model constant '" + name + "' is not used by this model (" +
          (name == "implementation" || name == "regval" ?
             std::string("it belonged to the removed pf3_frozen code path") :
             "seed = " + seed) +
          "); delete it from the parameters file");
  }

  void
  validate_and_report() const
  {
    auto require = [](const bool ok, const std::string &what) {
      if (!ok)
        throw std::runtime_error("parameters: " + what);
    };
    require(std::isfinite(static_cast<double>(u0)), "u0 must be finite");
    require(Gamma > 0.0, "Gamma must be > 0");
    require(lambda >= 0.0, "lambda must be >= 0");
    require(D_tilde > 0.0, "D_tilde must be > 0");
    require(Lsat >= 0.0, "Lsat must be >= 0");
    require(ic_interface_width > 0.0, "ic_interface_width must be > 0");
    if (seed_is_disc)
      require(disc_radius > 0.0 && disc_half_thickness > 0.0,
              "disc_radius and disc_half_thickness must be > 0");
    else
      {
        require(prism_outer_radius > prism_inner_radius && prism_inner_radius > 0.0,
                "need 0 < prism_inner_radius < prism_outer_radius");
        require(prism_outer_half_height > 0.0 && prism_inner_half_height > 0.0,
                "prism half heights must be > 0");
      }

    // A > 0 is guaranteed analytically for 0 <= eps_xy < 1 - eps_z (checked by demange::Anisotropy);
    // the sampled minimum is reported so the log documents the value (pole-regular: 1 - eps_z - eps_xy).
    const double amin = demange::sampled_min_A(aniso);
    require(amin > 0.0, "anisotropy A(n) is not positive");

    auto &out = ConditionalOStreams::pout_base();
    out << "\n==== Demange snow-crystal model (IceGrowthDemange) ====\n"
        << "  implementation       = modular (the only code path; pf3_frozen removed)\n"
        << "  half_convention      = " << half_convention << "\n"
        << "  q_normalization      = " << q_normalization << "  (q = " << q_factor
        << " (1 - phi))\n"
        << "  gamma_structure      = " << gamma_structure << "\n"
        << "  anisotropy           = " << anisotropy_model;
    if (anisotropy_model == "pole_regular")
      out << "  (pole_weight = " << pole_weight << ", psi_c = " << aniso.psi_c()
          << (aniso.has_basal_cap() ? ", basal cap B0 = " : ", no basal cap")
          << (aniso.has_basal_cap() ? std::to_string(aniso.B0()) : std::string()) << ")";
    out << "\n"
        << "  derivatives of A     = analytic (exact)\n"
        << "  u0 = " << u0 << ", eps_xy = " << eps_xy << ", eps_z = " << eps_z
        << ", Gamma = " << Gamma << ", lambda = " << lambda << ", D_tilde = " << D_tilde
        << ", Lsat = " << Lsat << "\n"
        << "  min A(n) (sampled)   = " << amin << "\n"
        << "  initial condition    = "
        << seed + (seed_is_disc ? std::string() : ", hex " + seed_hex_orientation)
        << ", center (" << seed_center[0] << ", " << seed_center[1] << ", " << seed_center[2]
        << "), ic_interface_width = " << ic_interface_width << "\n";
    if (seed_is_disc)
      out << "  seed geometry        = disc_radius " << disc_radius << ", disc_half_thickness "
          << disc_half_thickness << "\n";
    else
      out << "  seed geometry        = prism outer radius " << prism_outer_radius
          << ", outer half height " << prism_outer_half_height << ", inner radius "
          << prism_inner_radius << ", inner half height " << prism_inner_half_height << "\n";
    out << "==============================================================\n" << std::flush;
  }

  // ================================================================================================
  // Initial conditions
  // ================================================================================================
  void
  set_initial_condition([[maybe_unused]] const unsigned int       &index,
                        [[maybe_unused]] const unsigned int       &component,
                        [[maybe_unused]] const dealii::Point<dim> &point,
                        [[maybe_unused]] number                   &scalar_value,
                        [[maybe_unused]] number &vector_component_value) const override
  {
    set_initial_condition_modular(index, point, scalar_value);
  }

  // u = u0; phi = -tanh(d / (sqrt(2) w)) with d a signed-distance-like function (< 0 in ice) and
  // w = ic_interface_width. w = 1 is the planar equilibrium profile of the variational operator
  // (-tanh(x/sqrt(2))); use w = 1/sqrt(2) with half_convention = demange_literal (-tanh(x)).
  // xi1 (index 2) has no IC: it is computed from (u, phi) at increment 0.
  void
  set_initial_condition_modular(const unsigned int       &index,
                                const dealii::Point<dim> &point,
                                number                   &scalar_value) const
  {
    if (index == 0)
      {
        scalar_value = u0;
        return;
      }
    if (index != 1)
      return;
    const double dx = point[0] - seed_center[0];
    const double dy = point[1] - seed_center[1];
    const double dz = point[2] - seed_center[2];
    double       d  = 0.0;
    if (seed_is_disc)
      {
        // [P] Methods: circular-disc germ of radius R (thickness not stated: disc_half_thickness)
        d = std::max(std::hypot(dx, dy) - disc_radius, std::abs(dz) - disc_half_thickness);
      }
    else
      {
        // Hollow hexagonal prism (PF3 geometry). hex_radius <= R is a hexagon with apothem R.
        // PF3 orientation: facet normals at theta = 0, 60 deg (maxima of A). model_consistent swaps
        // x and y: facet normals at 30, 90 deg (minima of A, where the Wulff shape has its facets),
        // vertices along theta = 0.
        double ax = std::abs(dx);
        double ay = std::abs(dy);
        if (hex_swap_xy)
          std::swap(ax, ay);
        const double hex_radius = std::max(ax, 0.5 * ax + 0.5 * std::sqrt(3.0) * ay);
        const double outer =
          std::max(hex_radius - prism_outer_radius, std::abs(dz) - prism_outer_half_height);
        const double inner =
          std::max(hex_radius - prism_inner_radius, std::abs(dz) - prism_inner_half_height);
        d = std::max(outer, -inner);
      }
    scalar_value = -std::tanh(d / (std::sqrt(2.0) * ic_interface_width));
  }

  // ================================================================================================
  // Right-hand sides
  // ================================================================================================
  void
  compute_rhs(FieldContainer<dim, degree, number> &variable_list,
              const SimulationTimer               &sim_timer,
              unsigned int                         solve_block_id) const override
  {
    compute_rhs_modular(variable_list, sim_timer, solve_block_id);
  }

  // Explicit time-integration scheme for the phase-field formulation:
  //   block 0 (old values):  phi^{n+1} = phi^n + dt xi1^n / A^2
  //                          u^{n+1}   = u^n + dt div F2 - (Lsat/2) B (phi^{n+1} - phi^n)
  //                          F2 = D_tilde q(phi) G^T G grad u,  q = q_factor (1 - phi)
  //   block 1 (new values):  xi1 = phi - phi^3 + lambda B g'(phi) u + div F1,  g' = (1 - phi^2)^2
  //                          F1 = G^T (1/2 |g|^2 dA^2/dg + c_i A^2 G grad phi)   (demange_anisotropy.h)
  // Weak form: value terms r, gradient terms r_x = -F (natural BC F.nu = 0).
  void
  compute_rhs_modular(FieldContainer<dim, degree, number> &variable_list,
                      const SimulationTimer               &sim_timer,
                      unsigned int                         solve_block_id) const
  {
    const number dt = sim_timer.get_timestep();

    if (solve_block_id == 0)
      {
        const ScalarValue u    = variable_list.template get_value<Scalar, OldOne>(0);
        const ScalarGrad  ux   = variable_list.template get_gradient<Scalar, OldOne>(0);
        const ScalarValue phi  = variable_list.template get_value<Scalar, OldOne>(1);
        const ScalarGrad  phix = variable_list.template get_gradient<Scalar, OldOne>(1);
        const ScalarValue xi1  = variable_list.template get_value<Scalar, OldOne>(2);

        // A(n_g)^2 (kinetic coefficient) and B(n) at each SIMD lane. A >= min A > 0 (checked at
        // start-up), so no denominator floor is needed.
        ScalarValue A2;
        ScalarValue B_n;
        for (unsigned int i = 0; i < ScalarValue::size(); ++i)
          {
            const double m[3] = {static_cast<double>(phix[0][i]),
                                 static_cast<double>(phix[1][i]),
                                 static_cast<double>(phix[2][i])};
            const demange::PointValues pv = demange::point_values(aniso, flux_options, m);
            A2[i]  = static_cast<number>(pv.A * pv.A);
            B_n[i] = static_cast<number>(pv.B);
          }

        const ScalarValue q_phi = q_factor * (1.0 - phi);
        ScalarGrad        F2;
        F2[0] = D_tilde * q_phi * ux[0];
        F2[1] = D_tilde * q_phi * ux[1];
        F2[2] = D_tilde * q_phi * Gamma * Gamma * ux[2];

        const ScalarValue dphi   = dt * xi1 / A2; // phi^{n+1} - phi^n
        const ScalarValue eq_phi = phi + dphi;
        const ScalarValue eq_u   = u - (Lsat / static_cast<number>(2.0)) * B_n * dphi;
        const ScalarGrad  eqx_u  = -dt * F2;

        variable_list.set_value_term(0, eq_u);
        variable_list.set_gradient_term(0, eqx_u);
        variable_list.set_value_term(1, eq_phi);
      }
    else if (solve_block_id == 1)
      {
        const ScalarValue u    = variable_list.template get_value<Scalar, Current>(0);
        const ScalarValue phi  = variable_list.template get_value<Scalar, Current>(1);
        const ScalarGrad  phix = variable_list.template get_gradient<Scalar, Current>(1);

        ScalarValue B_n;
        ScalarGrad  F1;
        for (unsigned int i = 0; i < ScalarValue::size(); ++i)
          {
            const double m[3] = {static_cast<double>(phix[0][i]),
                                 static_cast<double>(phix[1][i]),
                                 static_cast<double>(phix[2][i])};
            const demange::PointValues pv = demange::point_values(aniso, flux_options, m);
            B_n[i]   = static_cast<number>(pv.B);
            F1[0][i] = static_cast<number>(pv.F[0]);
            F1[1][i] = static_cast<number>(pv.F[1]);
            F1[2][i] = static_cast<number>(pv.F[2]);
          }

        const ScalarValue minus_fprime = phi - phi * phi * phi;               // -f'(phi)
        const ScalarValue gprime       = (1.0 - phi * phi) * (1.0 - phi * phi); // g'(phi)
        const ScalarValue eq_xi1       = minus_fprime + lambda * B_n * gprime * u;
        const ScalarGrad  eqx_xi1      = -F1;

        variable_list.set_value_term(2, eq_xi1);
        variable_list.set_gradient_term(2, eqx_xi1);
      }
  }

  // ================================================================================================
  // Members
  // ================================================================================================
  // Snow / Demange model constants (same names and meaning as PF3)
  number u0;
  number eps_xy;
  number eps_z;
  number Gamma;
  number lambda;
  number D_tilde;
  number Lsat;

  // Model formulation options
  std::string half_convention;
  std::string q_normalization;
  std::string gamma_structure;
  std::string anisotropy_model;
  std::string pole_weight;

  // Initial-condition parameters
  std::string           seed;
  std::string           seed_hex_orientation;
  std::array<double, 3> seed_center;
  double                ic_interface_width;
  double                prism_outer_radius;
  double                prism_outer_half_height;
  double                prism_inner_radius;
  double                prism_inner_half_height;
  double                disc_radius;
  double                disc_half_thickness;

  // Derived
  demange::Anisotropy  aniso;
  demange::FluxOptions flux_options;
  double               q_factor     = 0.5;
  bool                 seed_is_disc = false;
  bool                 hex_swap_xy  = true;
};

PRISMS_PF_END_NAMESPACE
