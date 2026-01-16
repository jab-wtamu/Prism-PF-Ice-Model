// =====================================================================================
// customPDE.h  (ICE APP: u, phi) — Step 1: u diffusion only needs Dtilde
// =====================================================================================
//
// WHAT CHANGED (and WHY)
//
// - Added model constant:
//     Dtilde  (reduced diffusion coefficient for u equation)
//
// WHY:
// - equations.cc now computes:  ∂t u = Dtilde ∇·( q(phi) ∇u )
// - So Dtilde must be provided in parameters.prm and read here.
//
// Everything else remains minimal to keep the app stable.
// =====================================================================================

#include <core/matrixFreePDE.h>

using namespace dealii;

template <int dim, int degree>
class customPDE : public MatrixFreePDE<dim, degree>
{
public:
  customPDE(userInputParameters<dim> _userInputs)
    : MatrixFreePDE<dim, degree>(_userInputs)
    , userInputs(_userInputs)
  {}

  void
  setInitialCondition([[maybe_unused]] const Point<dim>  &p,
                      [[maybe_unused]] const unsigned int index,
                      [[maybe_unused]] double            &scalar_IC,
                      [[maybe_unused]] Vector<double>    &vector_IC) override;

  void
  setNonUniformDirichletBCs([[maybe_unused]] const Point<dim>  &p,
                            [[maybe_unused]] const unsigned int index,
                            [[maybe_unused]] const unsigned int direction,
                            [[maybe_unused]] const double       time,
                            [[maybe_unused]] double            &scalar_BC,
                            [[maybe_unused]] Vector<double>    &vector_BC) override;

private:
#include <core/typeDefs.h>

  const userInputParameters<dim> userInputs;

  void
  explicitEquationRHS(
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
    [[maybe_unused]] const VectorizedArray<double> element_volume) const override;

  void
  nonExplicitEquationRHS(
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
    [[maybe_unused]] const VectorizedArray<double> element_volume) const override;

  void
  equationLHS(
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
    [[maybe_unused]] const VectorizedArray<double> element_volume) const override;

#ifdef POSTPROCESS_FILE_EXISTS
  void
  postProcessedFields(
    [[maybe_unused]] const variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>       &pp_variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>>                     q_point_loc,
    [[maybe_unused]] const VectorizedArray<double> element_volume) const override;
#endif

#ifdef NUCLEATION_FILE_EXISTS
  double
  getNucleationProbability([[maybe_unused]] variableValueContainer variable_value,
                           [[maybe_unused]] double                 dV) const override;
#endif

  // ================================================================
  // Model constants needed by IC scaffold (phi seeds)
  // ================================================================
  Tensor<1, dim> center1 = userInputs.get_model_constant_rank_1_tensor("center1");
  Tensor<1, dim> center2 = userInputs.get_model_constant_rank_1_tensor("center2");

  double radius1 = userInputs.get_model_constant_double("radius1");
  double radius2 = userInputs.get_model_constant_double("radius2");

  // ================================================================
  // NEW: diffusion coefficient for u equation
  // ================================================================
  double Dtilde = userInputs.get_model_constant_double("Dtilde");
  
  // Coupling constant for phi equation (Demange Eq.1 local term, anisotropy off)
double lambda = userInputs.get_model_constant_double("lambda");

};
