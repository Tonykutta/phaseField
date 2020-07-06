#include "anisotropy_facet.h"

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

void variableAttributeLoader::loadVariableAttributes(){
	// Variable 0
	set_variable_name				(0,"c");
	set_variable_type				(0,SCALAR);
	set_variable_equation_type		(0,EXPLICIT_TIME_DEPENDENT);

	set_dependencies_value_term_RHS(0, "c,biharm,grad(biharm)");
	set_dependencies_gradient_term_RHS(0, "c, grad(c), n, grad(n),grad(mu),biharm,grad(biharm)");

	set_variable_name				(1,"mu");
	set_variable_type				(1,SCALAR);
	set_variable_equation_type		(1,AUXILIARY);

	set_dependencies_value_term_RHS(1, "c,n,biharm,grad(biharm)");
	set_dependencies_gradient_term_RHS(1, "grad(c),grad(biharm),n,grad(n)");
	// Variable 1
	set_variable_name				(2,"n");
	set_variable_type				(2,SCALAR);
	set_variable_equation_type		(2,EXPLICIT_TIME_DEPENDENT);

	set_dependencies_value_term_RHS(2, "c, n");
	set_dependencies_gradient_term_RHS(2, "");

	// Variable 2
	set_variable_name				(3,"biharm");
	set_variable_type				(3,SCALAR);
	set_variable_equation_type		(3,AUXILIARY);

	set_dependencies_value_term_RHS(3, "");
	set_dependencies_gradient_term_RHS(3, "grad(c),mu,grad(mu)");


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

		// --- Getting the values and derivatives of the model variables ---

		// Concentration
		scalarvalueType c = variable_list.get_scalar_value(0);
		scalargradType cx = variable_list.get_scalar_gradient(0);

		// Order parameter
		scalarvalueType n = variable_list.get_scalar_value(2);
		scalargradType nx = variable_list.get_scalar_gradient(2);
		scalargradType mux = variable_list.get_scalar_gradient(1);
		// Field for split formulation of the biharmonic term


		// --- Setting the expressions for the terms in the governing equations ---

		// Bulk terms
		scalarvalueType faV = 0.5*c*c/16.0;
		scalarvalueType facV = 0.5*c/8.0;
		scalarvalueType faccV = constV(0.5/8.0);
		scalarvalueType fbV = 0.5*(c-1.0)*(c-1.0)/16.0;
		scalarvalueType fbcV = 0.5*(c-1.0)/8.0;
		scalarvalueType fbccV = constV(0.5/8.0);
		scalarvalueType hV = 3.0*n*n-2.0*n*n*n;
		scalarvalueType hnV = 6.0*n-6.0*n*n;


		scalarvalueType normgradc = std::sqrt(cx.norm_square());
		scalargradType  normal = cx/(normgradc+constV(1.0e-16));

		scalarvalueType gamma;
		scalargradType dgammadnormal;
		// Anisotropy function calculates gamma and dgamma/dn
		anisotropy(normal, gamma, dgammadnormal);

		// Product of projection matrix and dgammadnorm vector
		scalargradType aniso;
		for (unsigned int i=0; i<dim; ++i){
			for (unsigned int j=0; j<dim; ++j){
				aniso[i] += -normal[i]*normal[j]*dgammadnormal[j];
				if (i==j) aniso[i] +=dgammadnormal[j];
			}
		}
		// Anisotropic gradient term
		aniso = gamma*(aniso*normgradc+gamma*cx);



		// --- Setting the expressions for the terms in the governing equations ---
		scalarvalueType eq_c = c;
		scalargradType eqx_c = constV(-McV*userInputs.dtValue)*mux;

		scalarvalueType eq_n = n-constV(userInputs.dtValue*MnV)*(fbV-faV)*hnV*gammal;
		//scalargradType eqx_n = constV(userInputs.dtValue*MnV)*(-aniso+constV(delta2)*biharmx);

		// --- Submitting the terms for the governing equations ---
		variable_list.set_scalar_value_term_RHS(0,eq_c);
		variable_list.set_scalar_gradient_term_RHS(0,eqx_c);

		variable_list.set_scalar_value_term_RHS(2,eq_n);
		//variable_list.set_scalar_gradient_term_RHS(1,eqx_n);

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


			// --- Getting the values and derivatives of the model variables ---
			scalarvalueType c = variable_list.get_scalar_value(0);
			scalargradType cx = variable_list.get_scalar_gradient(0);
			scalargradType nx = variable_list.get_scalar_gradient(2);
			scalarvalueType n = variable_list.get_scalar_value(2);
			scalarvalueType faV = 0.5*c*c/16.0;
			scalarvalueType facV = 0.5*c/8.0;
			scalarvalueType faccV = constV(0.5/8.0);
			scalarvalueType fbV = 0.5*(c-1.0)*(c-1.0)/16.0;
			scalarvalueType fbcV = 0.5*(c-1.0)/8.0;
			scalarvalueType fbccV = constV(0.5/8.0);
			scalarvalueType hV = 3.0*n*n-2.0*n*n*n;
			scalarvalueType hnV = 6.0*n-6.0*n*n;
			scalargradType biharmx = variable_list.get_scalar_gradient(3);
			// --- Setting the expressions for the terms in the governing equations ---


			// --- Submitting the terms for the governing equations ---


			scalarvalueType normgradc = std::sqrt(cx.norm_square());
			scalargradType  normal = cx/(normgradc+constV(1.0e-16));

			scalarvalueType gamma;
			scalargradType dgammadnormal;
			// Anisotropy function calculates gamma and dgamma/dn
			anisotropy(normal, gamma, dgammadnormal);

			// Product of projection matrix and dgammadnorm vector
			scalargradType aniso;
			for (unsigned int i=0; i<dim; ++i){
				for (unsigned int j=0; j<dim; ++j){
					aniso[i] += -normal[i]*normal[j]*dgammadnormal[j];
					if (i==j) aniso[i] +=dgammadnormal[j];
				}
			}
			// Anisotropic gradient term
			aniso = gamma*(aniso*normgradc+gamma*cx);

			scalarvalueType fcV = gamma*(facV*(1.0-hV)+fbcV*hV);

			scalarvalueType eq_mu = fcV;
			scalargradType eqx_mu = (aniso-constV(delta2)*biharmx);

			// --- Submitting the terms for the governing equations ---

			variable_list.set_scalar_value_term_RHS(1,eq_mu);
			variable_list.set_scalar_gradient_term_RHS(1,eqx_mu);

			variable_list.set_scalar_gradient_term_RHS(3,-cx);

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