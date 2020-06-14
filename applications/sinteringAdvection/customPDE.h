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

  mutable double IntegrationTimer= 0.0;
  mutable double InitTimer= this->currentTime;
  mutable double mVolume1=0.0;
  mutable double mVolume2=0.0;
  mutable double mVolume3=0.0;

  mutable double mForce1_x=0.0;
  mutable double mForce1_y=0.0;
  mutable double mForce2_x=0.0;
  mutable double mForce2_y=0.0;
  mutable double mForce3_x=0.0;
  mutable double mForce3_y=0.0;

  mutable double mTorque1_z=0.0;
  mutable double mTorque2_z=0.0;
  mutable double mTorque3_z=0.0;
  mutable double mTorque1_y=0.0;
  mutable double mTorque2_y=0.0;
  mutable double mTorque3_y=0.0;
  mutable double mTorque1_x=0.0;
  mutable double mTorque2_x=0.0;
  mutable double mTorque3_x=0.0;

  mutable double mCenter1_x=0.0;
  mutable double mCenter1_y=0.0;
  mutable double mCenter2_x=0.0;
  mutable double mCenter2_y=0.0;
  mutable double mCenter3_x=0.0;
  mutable double mCenter3_y=0.0;


  mutable double mAngleX1=0.0;
  mutable double mAngleY1=0.0;
  mutable double mAngleZ1=0.0;
  mutable double mAngleX2=0.0;
  mutable double mAngleY2=0.0;
  mutable double mAngleZ2=0.0;
  mutable double mAngleX3=0.0;
  mutable double mAngleY3=0.0;
  mutable double mAngleZ3=0.0;
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
    public:
      void  computeIntegralUser(const std::vector<vectorType*> variableSet) const;

      void UpdateAngle() const;
      void rotation(std::vector<double> &Angle, std::vector<double> &RotationVector, std::vector<double> &Tensor, std::vector<double> &TensorRot) const;
      void rotation1(std::vector<double> &Angle, std::vector<double> &RotationVector) const;

      //
      // ================================================================
      // Model constants specific to this subclass
      // ================================================================

      double D_surf = userInputs.get_model_constant_double("D_surf");
      double D_vap = userInputs.get_model_constant_double("D_vap");
      double D_vol = userInputs.get_model_constant_double("D_vol");
      double D_gb = userInputs.get_model_constant_double("D_gb");
      double Au = userInputs.get_model_constant_double("Au");
      double Bu = userInputs.get_model_constant_double("Bu");
      double k_nu = userInputs.get_model_constant_double("k_nu");
      double k_cu = userInputs.get_model_constant_double("k_cu");//	double McV = userInputs.get_model_constant_double("McV");
      double MnV = userInputs.get_model_constant_double("MnV");
      double gamma_s = userInputs.get_model_constant_double("gamma_s");
      double delta = userInputs.get_model_constant_double("delta");
      double gamma_gb_koeff = userInputs.get_model_constant_double("gamma_gb_koeff");
      bool Advection = userInputs.get_model_constant_bool("Advection");

      double K = userInputs.get_model_constant_double("K");
      double mR = userInputs.get_model_constant_double("mR");
      double mT = userInputs.get_model_constant_double("mT");
      dealii::Tensor<1,dim> center1 = userInputs.get_model_constant_rank_1_tensor("center1");
      dealii::Tensor<1,dim> center2 = userInputs.get_model_constant_rank_1_tensor("center2");
      dealii::Tensor<1,dim> center3 = userInputs.get_model_constant_rank_1_tensor("center3");
      double radius1 = userInputs.get_model_constant_double("radius1");
      double radius2 = userInputs.get_model_constant_double("radius2");
      double radius3 = userInputs.get_model_constant_double("radius3");

      double matrix_concentration = userInputs.get_model_constant_double("matrix_concentration");



};
