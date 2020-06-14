//#include "UpdateAngle.h"
//#include "try.h"
#include <iostream>
#include <cmath>
#include <math.h>

// =============================================================================================
// loadPostProcessorVariableAttributes: Set the attributes of the postprocessing variables
// =============================================================================================
// This function is analogous to 'loadVariableAttributes' in 'equations.h', but for
// the postprocessing expressions. It sets the attributes for each postprocessing
// expression, including its name, whether it is a vector or scalar (only scalars are
// supported at present), its dependencies on other variables and their derivatives,
// and whether to calculate an integral of the postprocessed quantity over the entire
// domain. Note: this function is not a member of customPDE.

void variableAttributeLoader::loadPostProcessorVariableAttributes(){

	// Variable 0
	set_variable_name				(0,"f_tot");
	set_variable_type				(0,SCALAR);
	set_dependencies_value_term_RHS(0, "c, grad(c),grad(n1),grad(n2),grad(mu),grad(c),n2,mu,n1,n3,grad(n3)");
	set_dependencies_gradient_term_RHS(0, "");
	set_output_integral         	(0,true);

	set_variable_name				(1,"vt1_x");
	set_variable_type				(1,SCALAR);
	set_dependencies_value_term_RHS(1, "c, grad(c),grad(n1),grad(n2),grad(mu),grad(c),n2,mu,n1,n3,grad(n3)");
	set_dependencies_gradient_term_RHS(1, "");
	set_output_integral         	(1,true);

	set_variable_name				(2,"vt1_y");
	set_variable_type				(2,SCALAR);
	set_dependencies_value_term_RHS(2, "c, grad(c),grad(n1),grad(n2),grad(mu),grad(c),n2,mu,n1,n3,grad(n3)");
	set_dependencies_gradient_term_RHS(2, "");
	set_output_integral         	(2,true);

	set_variable_name				(3,"vr1_x");
	set_variable_type				(3,SCALAR);
	set_dependencies_value_term_RHS(3, "c, grad(c),grad(n1),grad(n2),grad(mu),grad(c),n2,mu,n1,n3,grad(n3)");
	set_dependencies_gradient_term_RHS(3, "");
	set_output_integral         	(3,true);

	set_variable_name				(4,"vr1_y");
	set_variable_type				(4,SCALAR);
	set_dependencies_value_term_RHS(4, "c, grad(c),grad(n1),grad(n2),grad(mu),grad(c),n2,mu,n1,n3,grad(n3)");
	set_dependencies_gradient_term_RHS(4, "");
	set_output_integral         	(4,true);

	set_variable_name				(5,"center1_y");
	set_variable_type				(5,SCALAR);
	set_dependencies_value_term_RHS(5, "c, grad(c),grad(n1),grad(n2),grad(mu),grad(c),n2,mu,n1,n3,grad(n3)");
	set_dependencies_gradient_term_RHS(5, "");
	set_output_integral         	(5,true);

	set_variable_name				(6,"vt2_x");
	set_variable_type				(6,SCALAR);
	set_dependencies_value_term_RHS(6, "c, grad(c),grad(n1),grad(n2),grad(mu),grad(c),n2,mu,n1,n3,grad(n3)");
	set_dependencies_gradient_term_RHS(6, "");
	set_output_integral         	(6,true);

	set_variable_name				(7,"vt2_y");
	set_variable_type				(7,SCALAR);
	set_dependencies_value_term_RHS(7, "c, grad(c),grad(n1),grad(n2),grad(mu),grad(c),n2,mu,n1,n3,grad(n3)");
	set_dependencies_gradient_term_RHS(7, "");
	set_output_integral         	(7,true);

	set_variable_name				(8,"vr2_x");
	set_variable_type				(8,SCALAR);
	set_dependencies_value_term_RHS(8, "c, grad(c),grad(n1),grad(n2),grad(mu),grad(c),n2,mu,n1,n3,grad(n3)");
	set_dependencies_gradient_term_RHS(8, "");
	set_output_integral         	(8,true);

	set_variable_name				(9,"vr2_y");
	set_variable_type				(9,SCALAR);
	set_dependencies_value_term_RHS(9, "c, grad(c),grad(n1),grad(n2),grad(mu),grad(c),n2,mu,n1,n3,grad(n3)");
	set_dependencies_gradient_term_RHS(9, "");
	set_output_integral         	(9,true);


	set_variable_name				(10,"vt3_x");
	set_variable_type				(10,SCALAR);
	set_dependencies_value_term_RHS(10, "c, grad(c),grad(n1),grad(n2),grad(mu),grad(c),n2,mu,n1,n3,grad(n3)");
	set_dependencies_gradient_term_RHS(10, "");
	set_output_integral         	(10,true);

	set_variable_name				(11,"vt3_y");
	set_variable_type				(11,SCALAR);
	set_dependencies_value_term_RHS(11, "c, grad(c),grad(n1),grad(n2),grad(mu),grad(c),n2,mu,n1,n3,grad(n3)");
	set_dependencies_gradient_term_RHS(11, "");
	set_output_integral         	(11,true);

	set_variable_name				(12,"vr3_x");
	set_variable_type				(12,SCALAR);
	set_dependencies_value_term_RHS(12, "c, grad(c),grad(n1),grad(n2),grad(mu),grad(c),n2,mu,n1,n3,grad(n3)");
	set_dependencies_gradient_term_RHS(12, "");
	set_output_integral         	(12,true);

	set_variable_name				(13,"vr3_y");
	set_variable_type				(13,SCALAR);
	set_dependencies_value_term_RHS(13, "c, grad(c),grad(n1),grad(n2),grad(mu),grad(c),n2,mu,n1,n3,grad(n3)");
	set_dependencies_gradient_term_RHS(13, "");
	set_output_integral         	(13,true);

	set_variable_name				(14,"center1_x");
	set_variable_type				(14,SCALAR);
	set_dependencies_value_term_RHS(14, "c, grad(c),grad(n1),grad(n2),grad(mu),grad(c),n2,mu,n1,n3,grad(n3)");
	set_dependencies_gradient_term_RHS(14, "");
	set_output_integral         	(14,true);

	set_variable_name				(15,"center2_x");
	set_variable_type				(15,SCALAR);
	set_dependencies_value_term_RHS(15, "c, grad(c),grad(n1),grad(n2),grad(mu),grad(c),n2,mu,n1,n3,grad(n3)");
	set_dependencies_gradient_term_RHS(15, "");
	set_output_integral         	(15,true);

	set_variable_name				(16,"center2_y");
	set_variable_type				(16,SCALAR);
	set_dependencies_value_term_RHS(16, "c, grad(c),grad(n1),grad(n2),grad(mu),grad(c),n2,mu,n1,n3,grad(n3)");
	set_dependencies_gradient_term_RHS(16, "");
	set_output_integral         	(16,true);


	set_variable_name				(17,"center3_x");
	set_variable_type				(17,SCALAR);
	set_dependencies_value_term_RHS(17, "c, grad(c),grad(n1),grad(n2),grad(mu),grad(c),n2,mu,n1,n3,grad(n3)");
	set_dependencies_gradient_term_RHS(17, "");
	set_output_integral         	(17,true);

	set_variable_name				(18,"center3_y");
	set_variable_type				(18,SCALAR);
	set_dependencies_value_term_RHS(18, "c, grad(c),grad(n1),grad(n2),grad(mu),grad(c),n2,mu,n1,n3,grad(n3)");
	set_dependencies_gradient_term_RHS(18, "");
	set_output_integral         	(18,true);

	set_variable_name				(19,"angle1");
	set_variable_type				(19,SCALAR);
	set_dependencies_value_term_RHS(19, "c, grad(c),grad(n1),grad(n2),grad(mu),grad(c),n2,mu,n1,n3,grad(n3)");
	set_dependencies_gradient_term_RHS(19, "");
	set_output_integral         	(19,true);

	set_variable_name				(20,"angle2");
	set_variable_type				(20,SCALAR);
	set_dependencies_value_term_RHS(20, "c, grad(c),grad(n1),grad(n2),grad(mu),grad(c),n2,mu,n1,n3,grad(n3)");
	set_dependencies_gradient_term_RHS(20, "");
	set_output_integral         	(20,true);

	set_variable_name				(21,"angle3");
	set_variable_type				(21,SCALAR);
	set_dependencies_value_term_RHS(21, "c, grad(c),grad(n1),grad(n2),grad(mu),grad(c),n2,mu,n1,n3,grad(n3)");
	set_dependencies_gradient_term_RHS(21, "");
	set_output_integral         	(21,true);

}

// =============================================================================================
// postProcessedFields: Set the postprocessing expressions
// =============================================================================================
// This function is analogous to 'explicitEquationRHS' and 'nonExplicitEquationRHS' in
// equations.h. It takes in "variable_list" and "q_point_loc" as inputs and outputs two terms in
// the expression for the postprocessing variable -- one proportional to the test
// function and one proportional to the gradient of the test function. The index for
// each variable in this list corresponds to the index given at the top of this file (for
// submitting the terms) and the index in 'equations.h' for assigning the values/derivatives of
// the primary variables.

template <int dim,int degree>
void customPDE<dim,degree>::postProcessedFields(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
	variableContainer<dim,degree,dealii::VectorizedArray<double> > & pp_variable_list,
	const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

		scalarvalueType c = variable_list.get_scalar_value(0);
		scalargradType cx = variable_list.get_scalar_gradient(0);
		scalarvalueType n1 = variable_list.get_scalar_value(2);
		scalarvalueType n2 = variable_list.get_scalar_value(3);
		scalarvalueType n3 = variable_list.get_scalar_value(4);
		scalargradType nx1 = variable_list.get_scalar_gradient(2);
		scalargradType nx2 = variable_list.get_scalar_gradient(3);
		scalargradType nx3 = variable_list.get_scalar_gradient(4);
		scalargradType mux= variable_list.get_scalar_gradient(1);
		scalarvalueType mu = variable_list.get_scalar_value(1);

		scalarvalueType f_tot = constV(0.0);

		scalarvalueType f_chem = c*c*c*c - 2.0*c*c*c + c*c;

		scalarvalueType f_grad = constV(0.0);

		f_tot = f_chem + f_grad;


		pp_variable_list.set_scalar_value_term_RHS(0, f_tot);



		scalarvalueType pos_x;
		pos_x[0]=q_point_loc[0][0];
		scalarvalueType pos_y;
		pos_y[0]=q_point_loc[1][0];
		scalarvalueType vt1_x;
		vt1_x[0]=-mR/(customPDE::mVolume1)*(customPDE::mTorque1_z)*(pos_y[0]-customPDE::mCenter1_y)*n1[0];
		scalarvalueType vt1_y;
		vt1_y[0]=mR/(customPDE::mVolume1)*(customPDE::mTorque1_z)*(pos_x[0]-customPDE::mCenter1_x)*n1[0];
		scalarvalueType vt2_x;
		vt2_x[0]=-mR/(customPDE::mVolume2)*(customPDE::mTorque2_z)*(pos_y[0]-customPDE::mCenter2_y)*n2[0];
		scalarvalueType vt2_y;
		vt2_y[0]=mR/(customPDE::mVolume2)*(customPDE::mTorque2_z)*(pos_x[0]-customPDE::mCenter2_x)*n2[0];

		scalarvalueType vt3_x;
		vt3_x[0]=-mR/(customPDE::mVolume3)*(customPDE::mTorque3_z)*(pos_y[0]-customPDE::mCenter3_y)*n3[0];
		scalarvalueType vt3_y;
		vt3_y[0]=mR/(customPDE::mVolume3)*(customPDE::mTorque3_z)*(pos_x[0]-customPDE::mCenter3_x)*n3[0];

		scalarvalueType vvv1_x;
		vvv1_x[0]=mT*customPDE::mForce1_x/(customPDE::mVolume1)*n1[0];
		scalarvalueType vvv1_y;
		vvv1_y[0]=mT*customPDE::mForce1_y/(customPDE::mVolume1)*n1[0];
		scalarvalueType vvv2_x;
		vvv2_x[0]=mT*customPDE::mForce2_x/(customPDE::mVolume2)*n2[0];
		scalarvalueType vvv2_y;
		vvv2_y[0]=mT*customPDE::mForce2_y/(customPDE::mVolume2)*n2[0];
		scalarvalueType vvv3_x;
		vvv3_x[0]=mT*customPDE::mForce3_x/(customPDE::mVolume3)*n3[0];
		scalarvalueType vvv3_y;
		vvv3_y[0]=mT*customPDE::mForce3_y/(customPDE::mVolume3)*n3[0];





		scalarvalueType cent1_x;
		scalarvalueType cent1_y;
		scalarvalueType cent2_x;
		scalarvalueType cent2_y;
		scalarvalueType cent3_x;
		scalarvalueType cent3_y;

		scalarvalueType ang1;
		ang1[0]=customPDE::mAngleX1*n1[0];
		scalarvalueType ang2;
		ang2[0]=customPDE::mAngleX2*n2[0];
		scalarvalueType ang3;
		ang3[0]=customPDE::mAngleX3*n3[0];

		cent1_y[0]=customPDE::mCenter1_y*n1[0];
		cent1_x[0]=customPDE::mCenter1_x*n1[0];

		cent2_y[0]=customPDE::mCenter2_y*n2[0];
		cent2_x[0]=customPDE::mCenter2_x*n2[0];

		cent3_y[0]=customPDE::mCenter3_y*n3[0];
		cent3_x[0]=customPDE::mCenter3_x*n3[0];

		pp_variable_list.set_scalar_value_term_RHS(1, vvv1_x);
		pp_variable_list.set_scalar_value_term_RHS(2, vvv1_y);
		pp_variable_list.set_scalar_value_term_RHS(3, vt1_x);
		pp_variable_list.set_scalar_value_term_RHS(4, vt1_y);
		pp_variable_list.set_scalar_value_term_RHS(6, vvv2_x);
		pp_variable_list.set_scalar_value_term_RHS(7, vvv2_y);
		pp_variable_list.set_scalar_value_term_RHS(8, vt2_x);
		pp_variable_list.set_scalar_value_term_RHS(9, vt2_y);
		pp_variable_list.set_scalar_value_term_RHS(10, vvv3_x);
		pp_variable_list.set_scalar_value_term_RHS(11, vvv3_y);
		pp_variable_list.set_scalar_value_term_RHS(12, vt3_x);
		pp_variable_list.set_scalar_value_term_RHS(13, vt3_y);
		pp_variable_list.set_scalar_value_term_RHS(5, cent1_y);
		pp_variable_list.set_scalar_value_term_RHS(14, cent1_x);
		pp_variable_list.set_scalar_value_term_RHS(16, cent2_y);
		pp_variable_list.set_scalar_value_term_RHS(15, cent2_x);
		pp_variable_list.set_scalar_value_term_RHS(18, cent3_y);
		pp_variable_list.set_scalar_value_term_RHS(17, cent3_x);
		pp_variable_list.set_scalar_value_term_RHS(19, ang1);
		pp_variable_list.set_scalar_value_term_RHS(20, ang2);
		pp_variable_list.set_scalar_value_term_RHS(21, ang3);
	}
