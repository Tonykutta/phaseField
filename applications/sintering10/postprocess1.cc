// =================================================================================
// Set the attributes of the postprocessing variables
// =================================================================================
// This function is analogous to 'loadVariableAttributes' in 'equations.h', but for
// the postprocessing expressions. It sets the attributes for each postprocessing
// expression, including its name, whether it is a vector or scalar (only scalars are
// supported at present), its dependencies on other variables and their derivatives,
// and whether to calculate an integral of the postprocessed quantity over the entire
// domain.

void variableAttributeLoader::loadPostProcessorVariableAttributes(){

	// Variable 0
	set_variable_name				(0,"f_tot");
	set_variable_type				(0,SCALAR);

    set_dependencies_value_term_RHS(0, "n1,grad(n1),n2,grad(n2),c,grad(c)");
    set_dependencies_gradient_term_RHS(0, "");

    set_output_integral         	(0,true);

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

// --- Getting the values and derivatives of the model variables ---

// The order parameter and its derivatives
scalarvalueType c = variable_list.get_scalar_value(0);
scalargradType cx = variable_list.get_scalar_gradient(0);
scalarvalueType n1 = variable_list.get_scalar_value(2);
scalargradType nx1 = variable_list.get_scalar_gradient(2);
scalarvalueType n2 = variable_list.get_scalar_value(3);
scalargradType nx2 = variable_list.get_scalar_gradient(3);

// --- Setting the expressions for the terms in the postprocessing expressions ---

scalarvalueType f_tot = constV(0.0);
scalarvalueType f_grad = constV(0.0);
scalarvalueType f_chem = constV(0.0);

scalarvalueType gamma;
anisotropy(nx1,nx2,gamma);


scalarvalueType A= (12.0*constV(gamma_s)-7.0*constV(gamma_gb_koeff)*gamma)/constV(delta);
scalarvalueType B= constV(gamma_gb_koeff)*gamma/constV(delta);
scalarvalueType k_n =3.0/4.0*constV(delta)*constV(gamma_gb_koeff)*gamma;
scalarvalueType k_c=3.0/4.0*constV(delta)*(2.0*constV(gamma_s)-constV(gamma_gb_koeff)*gamma);

// The homogenous free energy
f_chem = A*c*c*(constV(1.0)-c)*(constV(1.0)-c)+B*(c*c+constV(6.0)*(constV(1.0)-c)*(n1*n1+n2*n2)-constV(4.0)*(constV(2.0)-c)*(n1*n1*n1+n2*n2*n2)+constV(3.0)*(n1*n1+n2*n2)*(n1*n1+n2*n2));

// The gradient free energy
f_grad = constV(0.5)*k_n*(nx1*nx1+nx2*nx2)+constV(0.5)*k_c*cx*cx;


// The total free energy
f_tot = f_chem + f_grad;

// --- Submitting the terms for the postprocessing expressions ---


pp_variable_list.set_scalar_value_term_RHS(0, f_tot);

}
