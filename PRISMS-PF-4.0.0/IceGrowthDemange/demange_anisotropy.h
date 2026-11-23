// Demange snow-crystal anisotropy A(n), its missing-orientation regularization, exact analytic
// derivatives, and the Cartesian phase-field flux F1 -- PPF 4.0 Physics Update.
//
// Header-only, no deal.II/PRISMS-PF dependency (unit-tested stand-alone:
// verification/cpp/test_demange_anisotropy.cc). Function-for-function mirror of the verified
// reference implementation verification/python/demange_aniso.py.
// Mathematics and evidence: docs/DERIVATIONS.md (sections H, P, R, D, G).
//
// Angle conventions ([P]; folded as in the PF3 code):
//   theta = atan2(n_y, n_x), psi = atan2(sqrt(n_x^2 + n_y^2), n_z)
//   th = |remainder(theta, pi/3)| in [0, pi/6]  (6-fold rotation + mirror theta -> -theta)
//   ps = atan2(rho, |n_z|)       in [0, pi/2]  (A is even in n_z)
// Every model returns A, dA/dth, dA/dps and A_th_over_sin = (dA/dth)/sin(ps) in folded angles.

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace demange
{
  inline constexpr double pi = 3.14159265358979323846;

  enum class Model
  {
    supplement_literal, // [P] A(n) + [S] Eq. 1-4 regularization (|theta| reading), exact roots
    pole_regular        // damp-then-regularize, smooth at n = +-z (DERIVATIONS.md section P)
  };

  enum class PoleWeight
  {
    sin4,
    sin6,
    sin8,
    quintic
  };

  enum Region : int
  {
    region_raw         = 0,
    region_theta_facet = 1,
    region_psi_cap     = 2, // [S] meridian-wise cap A2(theta) + B2(theta) cos(psi)
    region_blend       = 3, // [S] alpha blend
    region_basal_cap   = 4  // pole-regular axisymmetric cap B0 cos(psi)
  };

  struct AnisoEval
  {
    double A;
    double A_th;
    double A_ps;
    double A_th_over_sin;
    int    region;
  };

  // ------------------------------------------------------------------------------------------
  // [S] Eq. 2 (first line) in x = cos(2 theta_m):  E (20x^3 + 24x^2 - 3x - 6) = a
  // ------------------------------------------------------------------------------------------
  inline double
  P_cubic(const double x)
  {
    return 20.0 * x * x * x + 24.0 * x * x - 3.0 * x - 6.0;
  }

  inline double
  P_cubic_prime(const double x)
  {
    return 60.0 * x * x + 48.0 * x - 3.0;
  }

  // Root of P(x) = ratio on the admissible branch x in [1/2, 1] (theta_m in (0, pi/6]).
  // P(1/2) = 1, P(1) = 35 and P' > 0 there, so the root is unique for 1 <= ratio <= 35.
  // The caller only asks for ratio = a/E < 35 (theta-missing test E > a/35); ratio >= 1 is
  // guaranteed by the parameter check eps_xy < 1 - eps_z <= a.
  inline double
  solve_x_theta_m(double ratio)
  {
    if (!(ratio >= 1.0))
      throw std::domain_error("demange::solve_x_theta_m: a/E < 1 (A would not stay positive)");
    ratio      = std::min(ratio, 35.0);
    double lo  = 0.5;
    double hi  = 1.0;
    double x   = 0.75;
    const double eps = std::numeric_limits<double>::epsilon();
    for (int it = 0; it < 200; ++it)
      {
        const double f = P_cubic(x) - ratio;
        if (f > 0.0)
          hi = x;
        else
          lo = x;
        double xn = x - f / P_cubic_prime(x);
        if (!(xn > lo && xn < hi))
          xn = 0.5 * (lo + hi); // safeguarded Newton: bisection if Newton leaves the bracket
        if (std::abs(xn - x) <= 4.0 * eps * std::max(1.0, std::abs(x)))
          {
            x = xn;
            break;
          }
        x = xn;
      }
    return x;
  }

  // sin(6t)/sin(t) = 2 (cos t + cos 3t + cos 5t): exact identity, no 0/0 as t -> 0.
  inline double
  g_ratio(const double t)
  {
    return 2.0 * (std::cos(t) + std::cos(3.0 * t) + std::cos(5.0 * t));
  }

  struct ThetaFacet
  {
    double tm;  // theta_m
    double dtm; // d theta_m / d ps (only for the [S] alpha blend; diverges where theta_m -> 0)
    double A1;  // == 0 up to round-off (tangency)
    double dA1; // == 0 exactly
    double B1;
    double dB1;
  };

  // [S] theta-construction at one latitude for a*(1) + E cos(6 theta) (a, E and their ps-derivatives).
  // [S] Eq. 4: B1 = 6E sin(6tm)/sin(tm), A1 = a + E cos(6tm) - B1 cos(tm).
  // A1 == 0 for every latitude, and B1 sin(tm) = 6E sin(6tm), so differentiating A1 == 0 gives the
  // regular closed form dB1/dps = (da + dE cos 6tm)/cos tm (no dtheta_m/dps needed).
  inline ThetaFacet
  theta_facet(const double E, const double dE, const double a, const double da)
  {
    ThetaFacet r;
    const double x = solve_x_theta_m(a / E);
    r.tm           = 0.5 * std::acos(x);
    r.B1           = 6.0 * E * g_ratio(r.tm);
    r.A1           = a + E * std::cos(6.0 * r.tm) - r.B1 * std::cos(r.tm);
    r.dB1          = (da + dE * std::cos(6.0 * r.tm)) / std::cos(r.tm);
    r.dA1          = 0.0;
    const double dratio = (da * E - a * dE) / (E * E);
    const double dx     = dratio / P_cubic_prime(x);
    const double s2     = std::sqrt(std::max(0.0, 1.0 - x * x)); // = sin(2 tm)
    r.dtm = (s2 > 0.0) ? -dx / (2.0 * s2) :
                         std::copysign(std::numeric_limits<double>::infinity(), -dx);
    return r;
  }

  inline double
  ipow(const double x, const int p)
  {
    double r = 1.0;
    for (int k = 0; k < p; ++k)
      r *= x;
    return r;
  }

  // Latitude weight S(t), S'(t), S''(t) on t in [0, 1]; S(0) = 0, S(1) = 1, S'(1) = 0.
  struct Weight3
  {
    double S, S1, S2;
  };

  inline Weight3
  weight_shape(const PoleWeight kind, const double t)
  {
    if (t <= 0.0)
      return {0.0, 0.0, 0.0};
    if (t >= 1.0)
      return {1.0, 0.0, 0.0};
    if (kind == PoleWeight::quintic)
      {
        // unique quintic with S(0)=0, S(1)=1, S'=S''=0 at both ends
        return {t * t * t * (10.0 - 15.0 * t + 6.0 * t * t),
                30.0 * t * t * (1.0 - t) * (1.0 - t),
                60.0 * t * (1.0 - t) * (1.0 - 2.0 * t)};
      }
    // sin^p(pi t / 2); with psi_c = 0 this is sin^p(psi) (p = 6: harmonic sin^6 psi cos 6 theta)
    const int p = (kind == PoleWeight::sin4) ? 4 : (kind == PoleWeight::sin6) ? 6 : 8;
    const double k = 0.5 * pi;
    const double s = std::sin(k * t);
    const double c = std::cos(k * t);
    return {ipow(s, p),
            p * ipow(s, p - 1) * c * k,
            p * k * k * ((p - 1) * ipow(s, p - 2) * c * c - ipow(s, p))};
  }

  inline PoleWeight
  parse_pole_weight(const std::string &s)
  {
    if (s == "sin6")
      return PoleWeight::sin6;
    if (s == "sin4")
      return PoleWeight::sin4;
    if (s == "sin8")
      return PoleWeight::sin8;
    if (s == "quintic")
      return PoleWeight::quintic;
    throw std::invalid_argument("pole_weight must be sin6, sin4, sin8 or quintic, got '" + s + "'");
  }

  // ------------------------------------------------------------------------------------------
  // The anisotropy function
  // ------------------------------------------------------------------------------------------
  class Anisotropy
  {
  public:
    Anisotropy(const Model model_kind, const double eps_xy, const double eps_z,
               const PoleWeight pole_weight = PoleWeight::sin6)
      : model_(model_kind)
      , exy_(eps_xy)
      , ez_(eps_z)
      , weight_(pole_weight)
    {
      // A > 0 everywhere and theta_m well defined need 0 <= eps_xy < 1 - eps_z, 0 <= eps_z < 1.
      if (!(exy_ >= 0.0 && ez_ >= 0.0 && ez_ < 1.0 && exy_ < 1.0 - ez_))
        throw std::invalid_argument(
          "demange::Anisotropy: need 0 <= eps_xy < 1 - eps_z and 0 <= eps_z < 1");
      if (ez_ > 1.0 / 3.0)
        {
          // axisymmetric Eggleston cap: cos(2 psi_c) = 1/eps_z - 2, B0 = 4 eps_z cos(psi_c)
          psi_c_   = 0.5 * std::acos(1.0 / ez_ - 2.0);
          B0_      = std::sqrt(8.0 * ez_ * (1.0 - ez_));
          has_cap_ = true;
        }
      L_ = 0.5 * pi - psi_c_;
    }

    Model
    model() const
    {
      return model_;
    }

    double
    psi_c() const
    {
      return psi_c_;
    }

    double
    B0() const
    {
      return B0_;
    }

    bool
    has_basal_cap() const
    {
      return has_cap_;
    }

    AnisoEval
    eval(const double th, const double ps) const
    {
      return (model_ == Model::pole_regular) ? eval_pole_regular(th, ps) :
                                               eval_supplement(th, ps);
    }

    // Latitude weight w(ps) and its first two derivatives (pole-regular model).
    Weight3
    weight(const double ps) const
    {
      const Weight3 s = weight_shape(weight_, (ps - psi_c_) / L_);
      return {s.S, s.S1 / L_, s.S2 / (L_ * L_)};
    }

  private:
    // [P] A(n) with the [S] missing-orientation regularization (literal; singular at the poles).
    AnisoEval
    eval_supplement(const double th, const double ps) const
    {
      const double c6 = std::cos(6.0 * th);
      const double s6 = std::sin(6.0 * th);
      const double a  = 1.0 + ez_ * std::cos(2.0 * ps);
      const double da = -2.0 * ez_ * std::sin(2.0 * ps);

      AnisoEval raw{a + exy_ * c6, -6.0 * exy_ * s6, da, 0.0, region_raw};

      const bool tmiss = exy_ > a / 35.0;                   // [S] Eq. 1
      const bool pmiss = ez_ > (1.0 + exy_ * c6) / 3.0;     // [S] Eq. 1

      ThetaFacet tf{};
      AnisoEval  At{};
      if (tmiss)
        {
          tf = theta_facet(exy_, 0.0, a, da);
          At = {tf.A1 + tf.B1 * std::cos(th), -tf.B1 * std::sin(th),
                tf.dA1 + tf.dB1 * std::cos(th), 0.0, region_theta_facet};
        }
      double    pm = 0.0, dpm = 0.0;
      AnisoEval Ap{};
      if (pmiss)
        {
          const double C = std::clamp((1.0 + exy_ * c6) / ez_ - 2.0, -1.0, 1.0); // cos(2 psi_m)
          pm             = 0.5 * std::acos(C);
          const double dC  = -6.0 * exy_ * s6 / ez_;
          const double s2p = std::sin(2.0 * pm);
          dpm = (s2p > 0.0) ? -dC / (2.0 * s2p) :
                              std::copysign(std::numeric_limits<double>::infinity(), -dC);
          const double B2 = 4.0 * ez_ * std::cos(pm); // = 2 eps_z sin(2 pm)/sin(pm)
          const double A2 = 1.0 + exy_ * c6 + ez_ * std::cos(2.0 * pm) - B2 * std::cos(pm);
          // A2 == 0 for every theta (tangency) => dA2 = 0 and dB2/dth = -6 eps_xy sin(6th)/cos(pm)
          const double dB2 = -6.0 * exy_ * s6 / std::cos(pm);
          Ap = {A2 + B2 * std::cos(ps), dB2 * std::cos(ps), -B2 * std::sin(ps), 0.0,
                region_psi_cap};
        }

      AnisoEval r = raw;
      if (tmiss && pmiss)
        {
          if (th < tf.tm && ps >= pm)
            r = At;
          else if (th >= tf.tm && ps < pm)
            r = Ap;
          else if (th < tf.tm && ps < pm)
            {
              // [S] Eq. 4 alpha blend, alpha = |th - theta_m| / r (|theta| reading)
              const double d1 = th - tf.tm; // < 0 inside the blend region
              const double d2 = ps - pm;    // < 0
              const double rr = std::hypot(d1, d2);
              const double al = -d1 / rr;
              const double r3 = rr * rr * rr;
              // d alpha = -d2 (d2 dd1 - d1 dd2)/r^3 with dd1 = (1, -dtm), dd2 = (-dpm, 1)
              const double dal_th = -d2 * (d2 * 1.0 - d1 * (-dpm)) / r3;
              const double dal_ps = -d2 * (d2 * (-tf.dtm) - d1 * 1.0) / r3;
              r.A      = al * At.A + (1.0 - al) * Ap.A;
              r.A_th   = dal_th * (At.A - Ap.A) + al * At.A_th + (1.0 - al) * Ap.A_th;
              r.A_ps   = dal_ps * (At.A - Ap.A) + al * At.A_ps + (1.0 - al) * Ap.A_ps;
              r.region = region_blend;
            }
        }
      else if (tmiss)
        r = (th < tf.tm) ? At : raw;
      else if (pmiss)
        r = (ps < pm) ? Ap : raw;

      const double s = std::sin(ps);
      r.A_th_over_sin =
        (s > 0.0) ? r.A_th / s :
                    ((r.A_th == 0.0) ? 0.0 :
                                       std::copysign(std::numeric_limits<double>::infinity(), r.A_th));
      return r;
    }

    // Pole-regular (damp-then-regularize) model, DERIVATIONS.md section P.
    AnisoEval
    eval_pole_regular(const double th, const double ps) const
    {
      if (has_cap_ && ps < psi_c_)
        return {B0_ * std::cos(ps), 0.0, -B0_ * std::sin(ps), 0.0, region_basal_cap};

      const Weight3 w  = weight(ps);
      const double  E  = exy_ * w.S;
      const double  dE = exy_ * w.S1;
      const double  a  = 1.0 + ez_ * std::cos(2.0 * ps);
      const double  da = -2.0 * ez_ * std::sin(2.0 * ps);
      if (E > a / 35.0)
        {
          const ThetaFacet tf = theta_facet(E, dE, a, da);
          if (th < tf.tm)
            {
              const double A_th = -tf.B1 * std::sin(th);
              // sin(ps) > 0 here: E > a/35 > 0 needs w > 0, i.e. ps > psi_c >= 0
              return {tf.A1 + tf.B1 * std::cos(th), A_th, tf.dA1 + tf.dB1 * std::cos(th),
                      A_th / std::sin(ps), region_theta_facet};
            }
        }
      const double c6 = std::cos(6.0 * th);
      const double s6 = std::sin(6.0 * th);
      // A_th / sin(ps) = -6 eps_xy (w / sin ps) sin(6 th); w vanishes at least like sin^3 at the pole
      const double s          = std::sin(ps);
      const double w_over_sin = (s > 0.0) ? w.S / s : 0.0;
      return {a + E * c6, -6.0 * E * s6, da + dE * c6, -6.0 * exy_ * w_over_sin * s6, region_raw};
    }

    Model      model_;
    double     exy_;
    double     ez_;
    PoleWeight weight_;
    bool       has_cap_ = false;
    double     psi_c_   = 0.0;
    double     B0_      = 0.0;
    double     L_       = 0.5 * pi;
  };

  // ------------------------------------------------------------------------------------------
  // A(n) and its tangential (surface) gradient in Cartesian components for a unit normal n.
  //   grad_s A = A_psi e_psi + (A_theta / sin psi) e_theta,
  //   e_psi = (n_z cos theta, n_z sin theta, -rho), e_theta = (-sin theta, cos theta, 0).
  // pf3_floor_rel2 > 0 is used ONLY for the literal [S] model, to reproduce PF3's regularization of
  // the singular theta part: 1/sin(psi) -> sin(psi)/(sin^2(psi) + pf3_floor_rel2) with
  // pf3_floor_rel2 = (1e-8/|grad phi|)^2 (PF3 angle_floor); it has no physical justification.
  // ------------------------------------------------------------------------------------------
  inline double
  surface_gradient(const Anisotropy &aniso, const double n[3], double gs[3],
                   const double pf3_floor_rel2 = 0.0)
  {
    const double rho   = std::hypot(n[0], n[1]);
    const double ps    = std::atan2(rho, std::abs(n[2]));
    const double theta = std::atan2(n[1], n[0]);
    const double t     = std::remainder(theta, pi / 3.0);
    const double s_th  = (t >= 0.0) ? 1.0 : -1.0;
    const double s_z   = (n[2] >= 0.0) ? 1.0 : -1.0;
    const AnisoEval e  = aniso.eval(std::abs(t), ps);
    const double c     = std::cos(theta);
    const double s     = std::sin(theta);
    const double dpsi  = s_z * e.A_ps;
    double       dth_over_sin;
    if (aniso.model() == Model::supplement_literal && pf3_floor_rel2 > 0.0)
      dth_over_sin = s_th * e.A_th * rho / (rho * rho + pf3_floor_rel2);
    else
      dth_over_sin = s_th * e.A_th_over_sin;
    gs[0] = dpsi * n[2] * c - dth_over_sin * s;
    gs[1] = dpsi * n[2] * s + dth_over_sin * c;
    gs[2] = -dpsi * rho;
    return e.A;
  }

  // ------------------------------------------------------------------------------------------
  // Phase-field flux F1 (the phi equation contains + div F1) and the orientation values.
  // ------------------------------------------------------------------------------------------
  struct FluxOptions
  {
    bool   variational_half  = true;  // c_i = 1 (Karma-Rappel / d e/dm); false: [P] 1/2 on both
    bool   stretched         = false; // g = G m (variational in X = G^{-1} x); false: [P] g = m
    double Gamma             = 1.0;
    bool   pf3_pole_floor    = false; // literal [S] only: PF3-equivalent floor (see above)
  };

  struct PointValues
  {
    double A;    // A(n_g), used in the kinetic coefficient A^2 and in F1
    double B;    // B(n) = |G n| with the physical normal n ([P])
    double F[3]; // F1
  };

  // m = grad(phi). For m = 0 the normal is undefined; the flux is 0 (its exact limit) and A, B are
  // evaluated at n = +z. In exact arithmetic those values only multiply vanishing terms (bulk phases).
  inline PointValues
  point_values(const Anisotropy &aniso, const FluxOptions &o, const double m[3])
  {
    PointValues pv{};
    const double Gm[3] = {m[0], m[1], o.Gamma * m[2]};
    const double *g    = o.stretched ? Gm : m;
    const double  gn   = std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
    const double  mn   = std::sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);

    double nz_phys[3] = {0.0, 0.0, 1.0};
    if (mn > 0.0)
      {
        nz_phys[0] = -m[0] / mn;
        nz_phys[1] = -m[1] / mn;
        nz_phys[2] = -m[2] / mn;
      }
    pv.B = std::sqrt(nz_phys[0] * nz_phys[0] + nz_phys[1] * nz_phys[1] +
                     o.Gamma * o.Gamma * nz_phys[2] * nz_phys[2]);

    if (!(gn > 0.0))
      {
        const double up[3] = {0.0, 0.0, 1.0};
        double       gs[3];
        pv.A    = surface_gradient(aniso, up, gs);
        pv.F[0] = pv.F[1] = pv.F[2] = 0.0;
        return pv;
      }
    const double n[3] = {-g[0] / gn, -g[1] / gn, -g[2] / gn};
    double       gs[3];
    const double floor2 = (o.pf3_pole_floor) ? (1.0e-8 / gn) * (1.0e-8 / gn) : 0.0;
    const double A      = surface_gradient(aniso, n, gs, floor2);
    const double ci     = o.variational_half ? 1.0 : 0.5;
    // 1/2 |g|^2 dA^2/dg = -|g| A grad_s A ; the isotropic part is c_i A^2 (G m)
    double inner[3];
    for (int k = 0; k < 3; ++k)
      inner[k] = -gn * A * gs[k] + ci * A * A * Gm[k];
    pv.A    = A;
    pv.F[0] = inner[0];
    pv.F[1] = inner[1];
    pv.F[2] = o.Gamma * inner[2];
    return pv;
  }

  // Minimum of A over a dense sample of the folded orientation domain (start-up positivity check).
  inline double
  sampled_min_A(const Anisotropy &aniso, const int nt = 181, const int np = 361)
  {
    double amin = std::numeric_limits<double>::infinity();
    for (int i = 0; i < nt; ++i)
      for (int j = 0; j < np; ++j)
        amin = std::min(amin, aniso.eval((pi / 6.0) * i / (nt - 1), (0.5 * pi) * j / (np - 1)).A);
    return amin;
  }
} // namespace demange
