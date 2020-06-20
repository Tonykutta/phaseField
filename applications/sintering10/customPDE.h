#include "../../include/matrixFreePDE.h"

template <int dim, int degree>
class customPDE: public MatrixFreePDE<dim,degree>
{
public:
  // Constructor
  customPDE(userInputParameters<dim> _userInputs): MatrixFreePDE<dim,degree>(_userInputs) , userInputs(_userInputs) {};

  // Function to set the initial conditions (in ICs_and_BCs.h)
  void setInitialCondition(const dealii::Point<dim> &p, const unsigned int index, double & scalar_IC, dealii::Vector<double> & vector_IC);

  // Function to set the non-uniform Dirichlet boundary conditions (in ICs_and_BCs.h)
  void setNonUniformDirichletBCs(const dealii::Point<dim> &p, const unsigned int index, const unsigned int direction, const double time, double & scalar_BC, dealii::Vector<double> & vector_BC);

private:
  #include "../../include/typeDefs.h"

  const userInputParameters<dim> userInputs;

  // Function to set the RHS of the governing equations for explicit time dependent equations (in equations.h)
  void explicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
    dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

    // Function to set the RHS of the governing equations for all other equations (in equations.h)
    void nonExplicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
      dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

      // Function to set the LHS of the governing equations (in equations.h)
      void equationLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
        dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

        // Function to set postprocessing expressions (in postprocess.h)
        #ifdef POSTPROCESS_FILE_EXISTS
        void postProcessedFields(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
          variableContainer<dim,degree,dealii::VectorizedArray<double> > & pp_variable_list,
          const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;
          #endif

          // Function to set the nucleation probability (in nucleation.h)
          #ifdef NUCLEATION_FILE_EXISTS
          double getNucleationProbability(variableValueContainer variable_value, double dV) const;
          #endif

          // ================================================================
          // Methods specific to this subclass
          // ================================================================
          void anisotropy(const dealii::Tensor<1, dim, dealii::VectorizedArray<double> > &normal1,
            dealii::Tensor<1, dim, dealii::VectorizedArray<double> > &normal2,
            dealii::VectorizedArray<double> &gamma) const;

            void normtens(const dealii::Tensor<1, dim, dealii::VectorizedArray<double> > &normal1,
              dealii::Tensor<1, dim, dealii::VectorizedArray<double> > &normal2,
              dealii::Tensor<1, dim, dealii::VectorizedArray<double> > &normalc,
              dealii::Tensor<2, dim, dealii::VectorizedArray<double> > &normaltensn,
              dealii::Tensor<2, dim, dealii::VectorizedArray<double> > &normaltensc,
              dealii::Tensor<2, dim, dealii::VectorizedArray<double> > &unitens
            ) const;


            // ================================================================
            // Model constants specific to this subclass
            // ================================================================


            double MnV = userInputs.get_model_constant_double("MnV");
            double gamma_s = userInputs.get_model_constant_double("gamma_s");
            double delta = userInputs.get_model_constant_double("delta");
            double gamma_gb_koeff = userInputs.get_model_constant_double("gamma_gb_koeff");
            double D_surf = userInputs.get_model_constant_double("D_surf");
            double D_vap = userInputs.get_model_constant_double("D_vap");
            double D_vol = userInputs.get_model_constant_double("D_vol");
            double D_gb = userInputs.get_model_constant_double("D_gb");
            bool tensor = userInputs.get_model_constant_bool("tensor");
            dealii::Tensor<1,dim> center1 = userInputs.get_model_constant_rank_1_tensor("center1");
            dealii::Tensor<1,dim> center2 = userInputs.get_model_constant_rank_1_tensor("center2");

            double radius1 = userInputs.get_model_constant_double("radius1");
            double radius2 = userInputs.get_model_constant_double("radius2");

            double matrix_concentration = userInputs.get_model_constant_double("matrix_concentration");



            // ================================================================

          };