#include "../../include/matrixFreePDE.h"


 std::vector<double> mAngle{0.0,45.0,30.0,5.0,10.0,40.0,10.0,24.0,12.0,40.0,40.0,10.0,5.0,30.0,23.0};
 std::vector<double> mForce_x(15);
 std::vector<double> mForce_y(15);
 std::vector<double> mTorque_z(15);
 std::vector<double> mVolume(15);
 std::vector<double> mCenter_x(15);
  std::vector<double> mCenter_y(15);
double mGammagb=0.0;
    std::vector<double> mvec(30);
std::vector<std::vector<double>> mgammaij(30,mvec);

template <int dim, int degree>
class customPDE: public MatrixFreePDE<dim,degree>
{
public:

  // Construct
  customPDE(userInputParameters<dim> _userInputs): MatrixFreePDE<dim,degree>(_userInputs) , userInputs(_userInputs) {};

  // Function to set the initial conditions (in ICs_and_BCs.h)
  void setInitialCondition(const dealii::Point<dim> &p, const unsigned int index, double & scalar_IC, dealii::Vector<double> & vector_IC);

  // Function to set the non-uniform Dirichlet boundary conditions (in ICs_and_BCs.h)
  void setNonUniformDirichletBCs(const dealii::Point<dim> &p, const unsigned int index, const unsigned int direction, const double time, double & scalar_BC, dealii::Vector<double> & vector_BC);


private:
  #include "../../include/typeDefs.h"




  mutable double IntegrationTimer= 0.0;
  mutable double InitTimer= this->currentTime;


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
      void  computeIntegralUser(const std::vector<vectorType*> variableSet,  std::vector<double> &Force_x,std::vector<double> &Force_y,std::vector<double> &Volume,std::vector<double> &Center_x,std::vector<double> &Center_y,std::vector<double> &Torque_z) const;

      void UpdateAngle(std::vector<double> &Torque_z, std::vector<double> &Angle,std::vector<double> &Volume) const;
      void rotation(std::vector<double> &Angle, std::vector<double> &RotationVector, std::vector<double> &Tensor, std::vector<double> &TensorRot) const;
      void rotation1(std::vector<double> &Angle, std::vector<double> &RotationVector) const;
      void normtens(const dealii::Tensor<1, dim, dealii::VectorizedArray<double> > &normal1,
      dealii::Tensor<1, dim, dealii::VectorizedArray<double> > &normal2,
       dealii::Tensor<1, dim, dealii::VectorizedArray<double> > &normalc,
       dealii::Tensor<2, dim, dealii::VectorizedArray<double> > &normaltensn,
       dealii::Tensor<2, dim, dealii::VectorizedArray<double> > &normaltensc,
       dealii::Tensor<2, dim, dealii::VectorizedArray<double> > &unitens
     ) const;


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
      bool tensor = userInputs.get_model_constant_bool("tensor");
      unsigned int NumIndex = userInputs.get_model_constant_int("NumIndex");
      double K = userInputs.get_model_constant_double("K");
      double mR = userInputs.get_model_constant_double("mR");
      double mT = userInputs.get_model_constant_double("mT");
      dealii::Tensor<1,dim> center1 = userInputs.get_model_constant_rank_1_tensor("center1");
      dealii::Tensor<1,dim> center2 = userInputs.get_model_constant_rank_1_tensor("center2");



      double radius1 = userInputs.get_model_constant_double("radius1");
      double radius2 = userInputs.get_model_constant_double("radius2");




      double matrix_concentration = userInputs.get_model_constant_double("matrix_concentration");



};
