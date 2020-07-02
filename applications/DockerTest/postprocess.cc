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
	set_dependencies_value_term_RHS(0, "c, grad(c),grad(n1),grad(n2),grad(mu),grad(c),n2,mu,n1");
	set_dependencies_gradient_term_RHS(0, "");
	set_output_integral         	(0,true);

int numin=2;
unsigned int var_index=0;
unsigned int var_indexa=0;
unsigned int var_indexb=0;
unsigned int var_indexc=0;
unsigned int var_indexd=0;
	for ( var_index=1; var_index<=numin; var_index++){

	  // For the input file 'parameters_large_2D.in'
	  //for (unsigned int var_index=0; var_index<12; var_index++){
	  std::string var_name = "vrot_x";
	  var_name.append(std::to_string(var_index));

	  set_variable_name				(var_index,var_name);
	  set_variable_type				(var_index,SCALAR);
	 	set_dependencies_value_term_RHS(var_index, "c, grad(c),grad(n1),grad(n2),grad(mu),grad(c),n2,mu,n1");

	  // For the input file 'parameters.in'

	  set_dependencies_gradient_term_RHS(var_index, "");
	//		set_output_integral         	(var_index,true);

}

for ( var_indexa=1; var_indexa<=numin; var_indexa++){

	// For the input file 'parameters_large_2D.in'
	//for (unsigned int var_index=0; var_index<12; var_index++){
	std::string var_name = "vrot_y";
	var_name.append(std::to_string(var_indexa));

	set_variable_name				(var_indexa+numin,var_name);
	set_variable_type				(var_indexa+numin,SCALAR);
	set_dependencies_value_term_RHS(var_indexa+numin, "c, grad(c),grad(n1),grad(n2),grad(mu),grad(c),n2,mu,n1");

	// For the input file 'parameters.in'

	set_dependencies_gradient_term_RHS(var_indexa+numin, "");
	//	set_output_integral         	(var_indexa+var_index,true);

}

for ( var_indexb=1; var_indexb<=numin; var_indexb++){

	std::string var_name = "vtr_x";
	var_name.append(std::to_string(var_indexb));
	set_variable_name				(var_indexb+numin*2,var_name);
	set_variable_type				(var_indexb+numin*2,SCALAR);
	set_dependencies_value_term_RHS(var_indexb+numin*2, "c, grad(c),grad(n1),grad(n2),grad(mu),grad(c),n2,mu,n1");
	set_dependencies_gradient_term_RHS(var_indexb+numin*2, "");


}

for ( var_indexc=1; var_indexc<=numin; var_indexc++){

	// For the input file 'parameters_large_2D.in'
	//for (unsigned int var_index=0; var_index<12; var_index++){
	std::string var_name = "vtr_y";
	var_name.append(std::to_string(var_indexc));

	set_variable_name				(var_indexc+numin*3,var_name);
	set_variable_type				(var_indexc+numin*3,SCALAR);
	set_dependencies_value_term_RHS(var_indexc+numin*3, "c, grad(c),grad(n1),grad(n2),grad(mu),grad(c),n2,mu,n1");
  set_dependencies_gradient_term_RHS(var_indexc+numin*3, "");

}

for ( var_indexd=1; var_indexd<=numin; var_indexd++){

	std::string var_name = "ang";
	var_name.append(std::to_string(var_indexd));

	set_variable_name			         	(var_indexd+numin*4,var_name);
	set_variable_type				       (var_indexd+numin*4,SCALAR);
	set_dependencies_value_term_RHS(var_indexd+numin*4, "c, grad(c),grad(n1),grad(n2),grad(mu),grad(c),n2,mu,n1");
	set_dependencies_gradient_term_RHS(var_indexd+numin*4, "");


}

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

		scalarvalueType f_tot = constV(0.0);

		scalarvalueType f_chem = c*c*c*c - 2.0*c*c*c + c*c;

		scalarvalueType f_grad = constV(0.0);

		f_tot = f_chem + f_grad;


		pp_variable_list.set_scalar_value_term_RHS(0, f_tot);


		scalarvalueType pos_x;
	  pos_x[0]=q_point_loc[0][0];
	  scalarvalueType pos_y;
	  pos_y[0]=q_point_loc[1][0];



	  for(int j=0;j<NumIndex;j++){
	    scalarvalueType nj = variable_list.get_scalar_value(j+2);
	    scalargradType nxj = variable_list.get_scalar_gradient(j+2);

	    scalarvalueType vr_x;
	    vr_x[0]=-mR/(mVolume[j]+1.0e-16)*(mTorque_z[j])*(pos_y[0]-mCenter_y[j])*nj[0];
	    scalarvalueType vr_y;
	    vr_y[0]=mR/(mVolume[j]+1.0e-16)*(mTorque_z[j])*(pos_x[0]-mCenter_x[j])*nj[0];

	    scalarvalueType vt_x;
	    vt_x[0]=mT*mForce_x[j]/(mVolume[j] +1.0e-16)*nj[0];
	    scalarvalueType vt_y;
	    vt_y[0]=mT*mForce_y[j]/(mVolume[j]+1.0e-16)*nj[0];
			scalarvalueType an;
			an[0]=mAngle[j]*nj[0];

			pp_variable_list.set_scalar_value_term_RHS(j+1, vr_x);
      pp_variable_list.set_scalar_value_term_RHS(j+1+NumIndex, vr_y);
			pp_variable_list.set_scalar_value_term_RHS(j+1+NumIndex+NumIndex, vt_x);
			pp_variable_list.set_scalar_value_term_RHS(j+1+NumIndex+NumIndex+NumIndex, vt_y);
			pp_variable_list.set_scalar_value_term_RHS(j+1+NumIndex+NumIndex+NumIndex+NumIndex, an);

		}




	}
