// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
// This function sets attributes for each variable/equation in the app. The
// attributes are set via standardized function calls. The first parameter for each
// function call is the variable index (starting at zero). The first set of
// variable/equation attributes are the variable name (any string), the variable
// type (SCALAR/VECTOR), and the equation type (EXPLICIT_TIME_DEPENDENT/
// TIME_INDEPENDENT/AUXILIARY). The next set of attributes describe the
// dependencies for the governing equation on the values and derivatives of the
// other variables for the value term and gradient term of the RHS and the LHS.
// The final pair of attributes determine whether a variable represents a field
// that can nucleate and whether the value of the field is needed for nucleation
// rate calculations.

#include "UpdateAngle.h"
#include "Advection.h"
#include "../../include/matrixFreePDE.h"
#include <iostream>
#include <cmath>
#include <math.h>
void variableAttributeLoader::loadVariableAttributes(){
  // Variable 0
  set_variable_name				(0,"c");
  set_variable_type				(0,SCALAR);
  set_variable_equation_type		(0,EXPLICIT_TIME_DEPENDENT);
  set_dependencies_value_term_RHS(0, "n2,n1,c,grad(n2),grad(n1),grad(c),grad(mu),mu,n3,grad(n3)");
  set_dependencies_gradient_term_RHS(0, "grad(n2),grad(n1),c,grad(c),n2,n1,grad(mu),mu,n3,grad(n3)");

  // Variable 1
  set_variable_name				(1,"mu");
  set_variable_type				(1,SCALAR);
  set_variable_equation_type		(1,AUXILIARY);
  set_dependencies_value_term_RHS(1, "n2,n1,c,grad(n2),grad(n1),grad(c),mu,grad(mu),n3,grad(n3)");
  set_dependencies_gradient_term_RHS(1, "grad(n2),grad(n1),c,grad(c),n2,n1,n3,grad(n3),mu,grad(mu)");

  // Variable 2
  set_variable_name				(2,"n1");
  set_variable_type				(2,SCALAR);
  set_variable_equation_type		(2,EXPLICIT_TIME_DEPENDENT);
  set_dependencies_value_term_RHS(2, "n2,n1,c,grad(n2),grad(n1),grad(c),n3,grad(n3)");
  set_dependencies_gradient_term_RHS(2, "grad(n2),grad(n1),c,grad(c),n2,n1,n3,grad(n3)");
  // Variable 3
  set_variable_name				(3,"n2");
  set_variable_type				(3,SCALAR);
  set_variable_equation_type		(3,EXPLICIT_TIME_DEPENDENT);
  set_dependencies_value_term_RHS(3, "n2,n1,c,grad(n2),grad(n1),grad(c),n3,grad(n3)");
  set_dependencies_gradient_term_RHS(3, "grad(n2),grad(n1),c,grad(c),n2,n1,n3,grad(n3)");

  set_variable_name				(4,"n3");
  set_variable_type				(4,SCALAR);
  set_variable_equation_type		(4,EXPLICIT_TIME_DEPENDENT);
  set_dependencies_value_term_RHS(4, "n2,n1,c,grad(n2),grad(n1),grad(c),n3,grad(n3),mu,grad(mu)");
  set_dependencies_gradient_term_RHS(4, "grad(n2),grad(n1),c,grad(c),n2,n1,n3,grad(n3),mu,grad(mu)");


}

// =============================================================================================
// explicitEquationRHS (needed only if one or more equation is explict time dependent)
// =============================================================================================
// This function calculates the right-hand-side of the explicit time-dependent
// equations for each variable. It takes "variable_list" as an input, which is a list
// of the value and derivatives of each of the variables at a specific quadrature
// point. The (x,y,z) location of that quadrature point is given by "q_point_loc".
// The function outputs two terms to variable_list -- one proportional to the test
// function and one proportional to the gradient of the test function. The index for
// each variable in this list corresponds to the index given at the top of this file.

template <int dim, int degree>
void customPDE<dim,degree>::explicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
  dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

    scalarvalueType c = variable_list.get_scalar_value(0);
    scalargradType cx = variable_list.get_scalar_gradient(0);
    scalargradType mux = variable_list.get_scalar_gradient(1);
    scalarvalueType n2 = variable_list.get_scalar_value(3);
    scalargradType nx2 = variable_list.get_scalar_gradient(3);
    scalarvalueType n1 = variable_list.get_scalar_value(2);
    scalargradType nx1 = variable_list.get_scalar_gradient(2);
    scalarvalueType n3 = variable_list.get_scalar_value(4);
    scalargradType nx3 = variable_list.get_scalar_gradient(4);

    scalarvalueType advectionterm_n1;
    scalarvalueType advectionterm_n2;
    scalarvalueType advectionterm_n3;
    scalarvalueType advectionterm_c;




    scalarvalueType pos_x;
    pos_x[0]=q_point_loc[0][0];
    scalarvalueType pos_y;
    pos_y[0]=q_point_loc[1][0];
    scalarvalueType vt1_x;
    vt1_x[0]=-mR/(customPDE::mVolume1)*(customPDE::mTorque1_z)*(pos_y[0]-customPDE::mCenter1_y);
    scalarvalueType vt1_y;
    vt1_y[0]=mR/(customPDE::mVolume1)*(customPDE::mTorque1_z)*(pos_x[0]-customPDE::mCenter1_x);
    scalarvalueType vt2_x;
    vt2_x[0]=-mR/(customPDE::mVolume2)*(customPDE::mTorque2_z)*(pos_y[0]-customPDE::mCenter2_y);
    scalarvalueType vt2_y;
    vt2_y[0]=mR/(customPDE::mVolume2)*(customPDE::mTorque2_z)*(pos_x[0]-customPDE::mCenter2_x);
    scalarvalueType vt3_x;
    vt3_x[0]=-mR/(customPDE::mVolume3)*(customPDE::mTorque3_z)*(pos_y[0]-customPDE::mCenter3_y);
    scalarvalueType vt3_y;
    vt3_y[0]=mR/(customPDE::mVolume3)*(customPDE::mTorque3_z)*(pos_x[0]-customPDE::mCenter3_x);
    scalarvalueType vvv1_x;
    vvv1_x[0]=mT*customPDE::mForce1_x/(customPDE::mVolume1);
    scalarvalueType vvv1_y;
    vvv1_y[0]=mT*customPDE::mForce1_y/(customPDE::mVolume1);
    scalarvalueType vvv2_x;
    vvv2_x[0]=mT*customPDE::mForce2_x/(customPDE::mVolume2);
    scalarvalueType vvv2_y;
    vvv2_y[0]=mT*customPDE::mForce2_y/(customPDE::mVolume2);
    scalarvalueType vvv3_x;
    vvv3_x[0]=mT*customPDE::mForce3_x/(customPDE::mVolume3);
    scalarvalueType vvv3_y;
    vvv3_y[0]=mT*customPDE::mForce3_y/(customPDE::mVolume3);
    scalarvalueType  vtv1_x=(vvv1_x+vt1_x);
    scalarvalueType  vtv1_y=(vvv1_y+vt1_y);
    scalarvalueType  vtv2_x= (vvv2_x+vt2_x);
    scalarvalueType  vtv2_y=(vvv2_y+vt2_y);
    scalarvalueType  vtv3_x= (vvv3_x+vt3_x);
    scalarvalueType  vtv3_y=(vvv3_y+vt3_y);

    vectorvalueType v1;
    v1[0][0]=vtv1_x[0];
    v1[1][0]=vtv1_y[0];

    vectorvalueType v2;
    v2[0][0]=vtv2_x[0];
    v2[1][0]=vtv2_y[0];
    vectorvalueType v3;
    v3[0][0]=vtv3_x[0];
    v3[1][0]=vtv3_y[0];

    if (Advection==true){
      advectionterm_n1=(2.0*n1*nx1*v1);
      advectionterm_n2=(2.0*n2*nx2*v2);
      advectionterm_n3=(2.0*n3*nx3*v3);
      advectionterm_c=c*(nx1*v1+nx2*v2+nx3*v3)+cx*(n1*v1+n2*v2+n3*v3);

    }else{
      advectionterm_n1=constV(0.0);
      advectionterm_n2=constV(0.0);
      advectionterm_n3=constV(0.0);
      advectionterm_c=constV(0.0);
    }






    scalarvalueType fnV1 = (12.0*n1*(c*(n1-1.0)+n1*n1-2.0*n1+n2*n2+1.0+n3*n3))*constV(Bu);
    scalarvalueType fnV2 = (12.0*n2*(c*(n2-1.0)+n2*n2-2.0*n2+n1*n1+1.0+n3*n3))*constV(Bu);
    scalarvalueType fnV3 = (12.0*n3*(c*(n3-1.0)+n3*n3-2.0*n3+n2*n2+1.0+n1*n1))*constV(Bu);
    scalargradType diffmu=constV(D_surf)*mux;
    scalarvalueType eq_n2 = n2-constV(userInputs.dtValue)*(constV(MnV)*fnV2+advectionterm_n2);
    scalargradType eqx_n2 = (-constV(userInputs.dtValue*MnV)*k_nu*nx2);
    scalarvalueType eq_n1 = n1-constV(userInputs.dtValue)*(constV(MnV)*fnV1+advectionterm_n1);
    scalargradType eqx_n1 = (-constV(userInputs.dtValue*MnV)*k_nu*nx1);
    scalarvalueType eq_n3 = n3-constV(userInputs.dtValue)*(constV(MnV)*fnV3+advectionterm_n3);
    scalargradType eqx_n3 = (-constV(userInputs.dtValue*MnV)*constV(k_nu)*nx3);
    scalarvalueType eq_c = c-constV(userInputs.dtValue)*advectionterm_c;
    scalargradType eqx_c = -constV(userInputs.dtValue)*diffmu;

    variable_list.set_scalar_value_term_RHS(2,eq_n1);
    variable_list.set_scalar_gradient_term_RHS(2,eqx_n1);
    variable_list.set_scalar_value_term_RHS(3,eq_n2);
    variable_list.set_scalar_gradient_term_RHS(3,eqx_n2);
    variable_list.set_scalar_value_term_RHS(4,eq_n3);
    variable_list.set_scalar_gradient_term_RHS(4,eqx_n3);
    variable_list.set_scalar_value_term_RHS(0,eq_c);
    variable_list.set_scalar_gradient_term_RHS(0,eqx_c);


  }

  // =============================================================================================
  // nonExplicitEquationRHS (needed only if one or more equation is time independent or auxiliary)
  // =============================================================================================
  // This function calculates the right-hand-side of all of the equations that are not
  // explicit time-dependent equations. It takes "variable_list" as an input, which is
  // a list of the value and derivatives of each of the variables at a specific
  // quadrature point. The (x,y,z) location of that quadrature point is given by
  // "q_point_loc". The function outputs two terms to variable_list -- one proportional
  // to the test function and one proportional to the gradient of the test function. The
  // index for each variable in this list corresponds to the index given at the top of
  // this file.

  template <int dim, int degree>
  void customPDE<dim,degree>::nonExplicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
    dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

      scalarvalueType c = variable_list.get_scalar_value(0);
      scalargradType cx = variable_list.get_scalar_gradient(0);
      scalarvalueType n1 = variable_list.get_scalar_value(2);
      scalarvalueType n2 = variable_list.get_scalar_value(3);
      scalargradType nx2 = variable_list.get_scalar_gradient(3);
      scalargradType nx1 = variable_list.get_scalar_gradient(2);
      scalarvalueType n3 = variable_list.get_scalar_value(4);
      scalargradType nx3 = variable_list.get_scalar_gradient(4);

      // The derivative of the local free energy
      scalarvalueType fcV = 2.0*(constV(Au)*c*(2.0*c*c-3.0*c+1.0)+constV(Bu)*(c+2.0*n1*n1*n1-3.0*n1*n1+2.0*n2*n2*n2-3.0*n2*n2+2.0*n3*n3*n3-3.0*n3*n3));

      std::vector<vectorType*> variableSet=this->solutionSet;
      double t= this->currentTime;
      if(this->currentTime >customPDE::IntegrationTimer ||customPDE::InitTimer<userInputs.dtValue){
        computeIntegralUser(variableSet);
        customPDE::InitTimer=10.0*userInputs.dtValue;
        customPDE::IntegrationTimer=t;
        UpdateAngle();
      }

      scalarvalueType eq_mu = fcV;
      scalargradType eqx_mu = constV(k_cu)*cx;
      variable_list.set_scalar_value_term_RHS(1,eq_mu);
      variable_list.set_scalar_gradient_term_RHS(1,eqx_mu);

    }
    // =============================================================================================
    // equationLHS (needed only if at least one equation is time independent)
    // =============================================================================================
    // This function calculates the left-hand-side of time-independent equations. It
    // takes "variable_list" as an input, which is a list of the value and derivatives of
    // each of the variables at a specific quadrature point. The (x,y,z) location of that
    // quadrature point is given by "q_point_loc". The function outputs two terms to
    // variable_list -- one proportional to the test function and one proportional to the
    // gradient of the test function -- for the left-hand-side of the equation. The index
    // for each variable in this list corresponds to the index given at the top of this
    // file. If there are multiple elliptic equations, conditional statements should be
    // sed to ensure that the correct residual is being submitted. The index of the field
    // being solved can be accessed by "this->currentFieldIndex".

    template <int dim, int degree>
    void customPDE<dim,degree>::equationLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
      dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {


      }
