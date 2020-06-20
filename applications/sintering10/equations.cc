
#include "aniso.h"
#include "Diff.h"

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

	set_dependencies_value_term_RHS(0, "c");
	set_dependencies_gradient_term_RHS(0, "grad(mu),grad(c),grad(n1),grad(n2),n1,n2");

	// Variable 1
	set_variable_name				(1,"mu");
	set_variable_type				(1,SCALAR);
	set_variable_equation_type		(1,AUXILIARY);

	set_dependencies_value_term_RHS(1, "c,n1,n2");
	set_dependencies_gradient_term_RHS(1, "grad(c),grad(n1),grad(n2)");

	// Variable 2
	set_variable_name				(2,"n1");
	set_variable_type				(2,SCALAR);
	set_variable_equation_type		(2,EXPLICIT_TIME_DEPENDENT);

	set_dependencies_value_term_RHS(2, "n1,n2,c");
	set_dependencies_gradient_term_RHS(2, "grad(n1),grad(n2)");
	// Variable 3
	set_variable_name				(3,"n2");
	set_variable_type				(3,SCALAR);
	set_variable_equation_type		(3,EXPLICIT_TIME_DEPENDENT);

	set_dependencies_value_term_RHS(3, "n2,n1,c");
	set_dependencies_gradient_term_RHS(3, "grad(n2),grad(n1)");


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
		scalarvalueType c = variable_list.get_scalar_value(0);
		scalargradType cx = variable_list.get_scalar_gradient(0);

		scalargradType mux = variable_list.get_scalar_gradient(1);
		scalarvalueType n1 = variable_list.get_scalar_value(2);
		scalargradType nx1 = variable_list.get_scalar_gradient(2);
		scalarvalueType n2 = variable_list.get_scalar_value(3);
		scalargradType nx2 = variable_list.get_scalar_gradient(3);
		scalarvalueType gamma;
		anisotropy(nx1,nx2,gamma);

		vectorgradType tensn;
		vectorgradType tensc;
		vectorgradType unity;
		scalarvalueType phi;
		phi=c*c*c*(constV(10.0)-constV(15.0)*c+constV(6.0)*c*c);

		vectorgradType DIFF;
		scalarvalueType MOB;
		scalargradType diffmu;

		scalarvalueType A= (12.0*constV(gamma_s)-7.0*constV(gamma_gb_koeff)*gamma)/constV(delta);
		scalarvalueType B= constV(gamma_gb_koeff)*gamma/constV(delta);
		scalarvalueType k_n =3.0/4.0*constV(delta)*constV(gamma_gb_koeff)*gamma;
		if (tensor==true){
			normtens(nx1,nx2,cx,tensn,tensc,unity);
			DIFF=(constV(D_vol)*phi+constV(D_vap)*(constV(1.0)-phi))*unity+constV(D_surf)*c*(constV(1.0)-c)*tensc+constV(D_gb)*n1*n2*tensn;
			diffmu=DIFF*mux;
		}else
		{
			MOB=(constV(D_vol)*phi+constV(D_vap)*(constV(1.0)-phi))+constV(D_surf)*c*(constV(1.0)-c)+constV(D_gb)*n1*n2;
			diffmu=MOB*mux;
		}





		scalarvalueType fnV1 = (12.0*n1*(c*(n1-1.0)+n1*n1-2.0*n1+n2*n2+1.0))*B;
		scalarvalueType eq_n1 = (n1-constV(userInputs.dtValue*MnV)*k_n*fnV1);
		scalargradType eqx_n1 = (-constV(userInputs.dtValue*MnV)*nx1);
		scalarvalueType fnV2 = (12.0*n2*(c*(n2-1.0)+n2*n2-2.0*n2+n1*n1+1.0))*B;
		scalarvalueType eq_n2 = (n2-constV(userInputs.dtValue*MnV)*k_n*fnV2);
		scalargradType eqx_n2 = (-constV(userInputs.dtValue*MnV)*nx2);


		variable_list.set_scalar_value_term_RHS(2,eq_n1);
		variable_list.set_scalar_gradient_term_RHS(2,eqx_n1);

		variable_list.set_scalar_value_term_RHS(3,eq_n2);
		variable_list.set_scalar_gradient_term_RHS(3,eqx_n2);
		// --- Setting the expressions for the terms in the governing equations ---
		scalarvalueType eq_c = c;
		scalargradType eqx_c = -constV(userInputs.dtValue)*diffmu;

		// --- Submitting the terms for the governing equations ---
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

			// --- Getting the values and derivatives of the model variables ---

			scalarvalueType c = variable_list.get_scalar_value(0);
			scalargradType cx = variable_list.get_scalar_gradient(0);
			scalarvalueType n1 = variable_list.get_scalar_value(2);
			scalarvalueType n2 = variable_list.get_scalar_value(3);
			scalargradType nx1 = variable_list.get_scalar_gradient(2);
			scalargradType nx2 = variable_list.get_scalar_gradient(3);

			// --- Setting the expressions for the terms in the governing equations ---

			scalarvalueType gamma;
			anisotropy(nx1,nx2,gamma);
			scalarvalueType A= (12.0*constV(gamma_s)-7.0*constV(gamma_gb_koeff)*gamma)/constV(delta);
			scalarvalueType B= constV(gamma_gb_koeff)*gamma/constV(delta);
			// The derivative of the local free energy
			scalarvalueType fcV = 2.0*(A*c*(2.0*c*c-3.0*c+1.0)+B*(c+2.0*n1*n1*n1-3.0*n1*n1+2.0*n2*n2*n2-3.0*n2*n2));

			// The terms for the governing equations

			scalarvalueType k_c=3.0/4.0*constV(delta)*(2.0*constV(gamma_s)-constV(gamma_gb_koeff)*gamma);

			scalarvalueType eq_mu = fcV;
			scalargradType eqx_mu = k_c*cx;

			// --- Submitting the terms for the governing equations ---

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
