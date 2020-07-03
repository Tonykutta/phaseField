#include "anisotropy_facet.h"
#include "anisotropy_gb.h"
#define _USE_MATH_DEFINES

#include <cmath>
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
  set_dependencies_gradient_term_RHS(0, "c, grad(c), n1, grad(n1),grad(mu),biharm,grad(biharm),n2,grad(n2),grad(biharmn1),grad(biharmn2),biharmn1,biharmn2");

    // Variable 1
	set_variable_name				(1,"n1");
	set_variable_type				(1,SCALAR);
	set_variable_equation_type		(1,EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(1, "c, n1,n2");
  set_dependencies_gradient_term_RHS(1, "grad(n1),grad(n2),n1,n2,c,grad(c),grad(biharmn1)");//

	set_variable_name				(2,"n2");
	set_variable_type				(2,SCALAR);
	set_variable_equation_type		(2,EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(2, "c, n2,n1");
  set_dependencies_gradient_term_RHS(2, "grad(n1),grad(n2),n1,n2,c,grad(c),grad(biharmn2)"); //

	set_variable_name				(3,"biharm");
	set_variable_type				(3,SCALAR);
	set_variable_equation_type		(3,AUXILIARY);

  set_dependencies_value_term_RHS(3, "");
  set_dependencies_gradient_term_RHS(3, "grad(c),mu,grad(mu),biharmn1,biharmn2,grad(biharmn1),grad(biharmn2)");

	set_variable_name				(4,"mu");
	set_variable_type				(4,SCALAR);
	set_variable_equation_type		(4,AUXILIARY);

	set_dependencies_value_term_RHS(4, "c,n1,n2,grad(c),grad(biharm),n1,grad(n1),n2,grad(n2),c");
	set_dependencies_gradient_term_RHS(4, "grad(c),grad(biharm),n1,grad(n1),n2,grad(n2),c");//

	set_variable_name				(5,"biharmn1");
	set_variable_type				(5,SCALAR);
	set_variable_equation_type		(5,AUXILIARY);

	set_dependencies_value_term_RHS(5, "");
	set_dependencies_gradient_term_RHS(5, "grad(n1),mu,grad(n1),biharmn2,grad(biharmn2)");


		set_variable_name				(6,"biharmn2");
		set_variable_type				(6,SCALAR);
		set_variable_equation_type		(6,AUXILIARY);

		set_dependencies_value_term_RHS(6, "");
		set_dependencies_gradient_term_RHS(6, "grad(n2),mu,grad(mu),biharmn1,n2,grad(biharmn1)");
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
scalarvalueType n1 = variable_list.get_scalar_value(1);
scalarvalueType n2 = variable_list.get_scalar_value(2);
scalargradType nx1 = variable_list.get_scalar_gradient(1);
scalargradType nx2 = variable_list.get_scalar_gradient(2);
scalargradType mux = variable_list.get_scalar_gradient(4);
scalargradType biharmn1x = variable_list.get_scalar_gradient(5);
scalargradType biharmn2x = variable_list.get_scalar_gradient(6);


/*scalarvalueType faV = 0.5*c*c/16.0;
scalarvalueType facV = 0.5*c/8.0;
scalarvalueType faccV = constV(0.5/8.0);
scalarvalueType fbV = 0.5*(c-1.0)*(c-1.0)/16.0;
scalarvalueType fbcV = 0.5*(c-1.0)/8.0;
scalarvalueType fbccV = constV(0.5/8.0);
scalarvalueType n=n1;
scalarvalueType hV = 3.0*n*n-2.0*n*n*n;
scalarvalueType hnV = 6.0*n-6.0*n*n;
scalarvalueType quap=(fbV-faV)*hnV;*/
////////////////////////////////////////////////////////////////////////////////
scalarvalueType normgradn1 = std::sqrt(nx1.norm_square());
scalargradType  normaln1 = nx1/(normgradn1+constV(1.0e-16));
scalarvalueType gamman1;
scalargradType dgammadnormaln1;
anisotropy(normaln1, gamman1, dgammadnormaln1);
scalargradType anison1;
for (unsigned int i=0; i<dim; ++i){
		for (unsigned int j=0; j<dim; ++j){
				anison1[i] += -normaln1[i]*normaln1[j]*dgammadnormaln1[j];
				if (i==j) anison1[i] +=dgammadnormaln1[j];
		}
}
scalargradType anison1x;
if(normgradn1[0]<0.1){
	anison1x=constV(0.0);
}else{
anison1x=anison1;
}


anison1 = (anison1*normgradn1+gamman1*nx1);


scalarvalueType normgradn2 = std::sqrt(nx2.norm_square());

scalargradType  normaln2 = nx2/(normgradn2+constV(1.0e-16));
scalarvalueType gamman2;
scalargradType dgammadnormaln2;
anisotropy(normaln2, gamman2, dgammadnormaln2);
scalargradType anison2;
for (unsigned int i=0; i<dim; ++i){
		for (unsigned int j=0; j<dim; ++j){
				anison2[i] += -normaln2[i]*normaln2[j]*dgammadnormaln2[j];
				if (i==j) anison2[i] +=dgammadnormaln2[j];
		}
}

scalargradType anison2x;
if(normgradn2[0]<0.1){
	anison2x=constV(0.0);
}else{
anison2x=anison2;
}
anison2 = (anison2*normgradn2+gamman2*nx2);

////////////////////////////////////////////////////////////////////


scalarvalueType normgradgb = std::sqrt((nx1-nx2).norm_square());
scalargradType  normalgb = (nx1-nx2)/(normgradgb+constV(1.0e-16));
scalarvalueType gammagb;
scalargradType dgammadnormalgb;
scalarvalueType gamma_0;
scalargradType  normalnz;
scalargradType  anisogb;
scalargradType  anisogb1;
scalargradType  anisogb2;
double  alpha=M_PI/4.0;
double beta=0.0;
double a=0.2;
if(dim ==3){
normalnz[0][0]=std::cos(alpha)*std::cos(beta)*normalgb[0][0]-std::sin(alpha)*normalgb[1][0]+std::cos(alpha)*std::sin(beta)*normalgb[0][2];
normalnz[1][0]=std::sin(alpha)*std::cos(beta)*normalgb[0][0]+std::cos(alpha)*normalgb[1][0]+std::sin(alpha)*std::sin(beta)*normalgb[0][2];
normalnz[2][0]=-std::sin(beta)*normalgb[0][0]+std::cos(beta)*normalgb[0][2];
gammagb=constV(1.0)*(constV(1.0)-constV(a)*((4.0)*(normalnz[0]*normalnz[0]*normalnz[0]*normalnz[0]+normalnz[1]*normalnz[1]*normalnz[1]*normalnz[1]+normalnz[2]*normalnz[2]*normalnz[2]*normalnz[2])-constV(3.0)));
dgammadnormalgb[0][0]=a*16.0*(normalnz[0][0]*normalnz[0][0]*normalnz[0][0]*std::cos(alpha)*std::cos(beta)+normalnz[1][0]*normalnz[1][0]*normalnz[1][0]*std::sin(alpha)*std::cos(beta)-normalnz[2][0]*normalnz[2][0]*normalnz[2][0]*std::sin(beta));
dgammadnormalgb[1][0]=a*16.0*(-normalnz[0][0]*normalnz[0][0]*normalnz[0][0]*std::sin(alpha)+normalnz[1][0]*normalnz[1][0]*normalnz[1][0]*std::cos(alpha));
dgammadnormalgb[2][0]=a*16.0*(normalnz[0][0]*normalnz[0][0]*normalnz[0][0]*std::cos(alpha)*std::sin(beta)+normalnz[1][0]*normalnz[1][0]*normalnz[1][0]*std::sin(alpha)*std::sin(beta)+normalnz[2][0]*normalnz[2][0]*normalnz[2][0]*std::cos(beta));

}
else{
	normalnz[0][0]=std::cos(alpha)*normalgb[0][0]-std::sin(alpha)*normalgb[1][0];
	normalnz[1][0]=std::sin(alpha)*normalgb[0][0]+std::cos(alpha)*normalgb[1][0];
	gammagb=constV(1.0)*(constV(1.0)-constV(a)*((4.0)*(normalnz[0]*normalnz[0]*normalnz[0]*normalnz[0]+normalnz[1]*normalnz[1]*normalnz[1]*normalnz[1])-constV(3.0)));
	dgammadnormalgb[0][0]=a*16.0*(normalnz[0][0]*normalnz[0][0]*normalnz[0][0]*std::cos(alpha)+normalnz[1][0]*normalnz[1][0]*normalnz[1][0]*std::sin(alpha));
	dgammadnormalgb[1][0]=a*16.0*(-normalnz[0][0]*normalnz[0][0]*normalnz[0][0]*std::sin(alpha)+normalnz[1][0]*normalnz[1][0]*normalnz[1][0]*std::cos(alpha));

}
scalarvalueType phi;
phi[0]=std::atan(normalgb[0][0]/normalgb[1][0])+M_PI/4.0;

anisotropy_gb(normalnz, gammagb, dgammadnormalgb,gamma_0);

scalargradType anison;
for (unsigned int i=0; i<dim; ++i){
		for (unsigned int j=0; j<dim; ++j){
				anisogb[i] += -normalgb[i]*normalgb[j]*dgammadnormalgb[j];
				if (i==j) anisogb[i] +=dgammadnormalgb[j];
		}
}
scalargradType anisoc1= constV(3.0/4.0*0.5)*(2.0*anison1-anisogb)*cx*cx;
scalargradType anisoc2= constV(3.0/4.0*0.5)*(2.0*anison2-anisogb)*cx*cx;
anisogb1 = constV(3.0/4.0)*(anisogb*normgradgb+gammagb*nx1);
anisogb2 = constV(3.0/4.0)*(anisogb*normgradgb+gammagb*nx2);
gammagb=constV(0.5);
scalarvalueType fnV1 = 12.0*n1*(c*(n1-1.0)+n1*n1-2.0*n1+n2*n2+1.0)*gammagb/constV(del);
scalarvalueType fnV2 = 12.0*n2*(c*(n2-1.0)+n2*n2-2.0*n2+n1*n1+1.0)*gammagb/constV(del);
scalarvalueType eq_c = c;
scalargradType eqx_c = constV(-McV*userInputs.dtValue)*mux;
scalarvalueType eq_n1 = n1-constV(userInputs.dtValue*MnV)*(fnV1);
scalarvalueType eq_n2 = n2-constV(userInputs.dtValue*MnV)*(fnV2);

scalargradType der_n1= (anison1x)*(constV(12.0/del)*(c*c*((constV(1.0)-c)*(constV(1.0)-c)))+constV(3.0/4.0*del)*cx*cx);  //+anisogb*(c*c+constV(6.0)*(constV(1.0)-c)*(n1*n1+n2*n2)-constV(4.0)*(constV(2.0)-c)*(n1*n1*n1+n2*n2*n2)+constV(3.0)*(n1*n1+n2*n2)*(n1*n1+n2*n2));
scalargradType der_n2= (anison2x)*(constV(12.0/del)*(c*c*((constV(1.0)-c)*(constV(1.0)-c)))+constV(3.0/4.0*del)*cx*cx);  //  +anisogb*(c*c+constV(6.0)*(constV(1.0)-c)*(n1*n1+n2*n2)-constV(4.0)*(constV(2.0)-c)*(n1*n1*n1+n2*n2*n2)+constV(3.0)*(n1*n1+n2*n2)*(n1*n1+n2*n2));
//normgradn2
scalargradType eqx_n1 = -constV(userInputs.dtValue*MnV)*(constV(3.0/4.0*del)*gammagb*nx1-constV(delta1)*biharmn1x+der_n1);//-anison1  //-anisogb1+   der_n1
scalargradType eqx_n2 = -constV(userInputs.dtValue*MnV)*(constV(3.0/4.0*del)*gammagb*nx2-constV(delta1)*biharmn2x+der_n2);//-anison2  //-anisogb2+   der_n2
variable_list.set_scalar_value_term_RHS(0,eq_c);
variable_list.set_scalar_gradient_term_RHS(0,eqx_c);

variable_list.set_scalar_value_term_RHS(1,eq_n1);
variable_list.set_scalar_value_term_RHS(2,eq_n2);
variable_list.set_scalar_gradient_term_RHS(1,eqx_n1);
variable_list.set_scalar_gradient_term_RHS(2,eqx_n2);
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
 scalargradType nx1 = variable_list.get_scalar_gradient(1);
 scalargradType nx2 = variable_list.get_scalar_gradient(2);
 scalarvalueType n1 = variable_list.get_scalar_value(1);
 scalarvalueType n2 = variable_list.get_scalar_value(2);
 scalargradType biharmx = variable_list.get_scalar_gradient(3);
 variable_list.set_scalar_gradient_term_RHS(3,-cx),
 variable_list.set_scalar_gradient_term_RHS(5,-nx1);
 variable_list.set_scalar_gradient_term_RHS(6,-nx2);
 scalarvalueType faV = 0.5*c*c/16.0;
 scalarvalueType facV = 0.5*c/8.0;
 scalarvalueType faccV = constV(0.5/8.0);
 scalarvalueType fbV = 0.5*(c-1.0)*(c-1.0)/16.0;
 scalarvalueType fbcV = 0.5*(c-1.0)/8.0;
 scalarvalueType fbccV = constV(0.5/8.0);
 scalarvalueType n=n1;
 scalarvalueType hV = 3.0*n*n-2.0*n*n*n;
 scalarvalueType hnV = 6.0*n-6.0*n*n;
 scalarvalueType quip= facV*(1.0-hV)+fbcV*hV;



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
 scalargradType anisox;
 if(normgradc[0]<0.1){
 anisox=constV(0.0);
}else{
anisox=aniso;
}
 //aniso = (aniso*normgradc+(gamma)*cx);//-(gamma_0+gamman)

////////////////////////////////////////////////////////////////////

 /*scalarvalueType normgradgb = std::sqrt((nx1).norm_square());
 scalargradType  normalgb = (nx1)/(normgradgb+constV(1.0e-16));
 scalarvalueType gammagb;
 scalargradType dgammadnormalgb;
 scalarvalueType gamma_0;
 scalargradType  normalnz;
 normalnz[0][0]=std::cos(M_PI/4.0)*normalgb[0][0]-std::sin(M_PI/4.0)*normalgb[1][0];
 normalnz[1][0]=std::sin(M_PI/4.0)*normalgb[0][0]+std::cos(M_PI/4.0)*normalgb[1][0];
 scalarvalueType phi;
 phi[0]=std::atan(normalgb[0][0]/normalgb[1][0])+M_PI/4.0;

 anisotropy_gb(normalnz, gammagb, dgammadnormalgb,gamma_0);
 gammagb[0]=1.0-0.8*std::cos(phi[0]);


 gammagb=constV(1.0)*(constV(1.0)-0.9*((4.0)*(normalnz[0]*normalnz[0]*normalnz[0]*normalnz[0]+normalnz[1]*normalnz[1]*normalnz[1]*normalnz[1])-constV(3.0)));
*/
////////////////////////////////////////////////////////////////////////////

//scalarvalueType gammagb=constV(0.5);




scalarvalueType normgradgb = std::sqrt((nx1-nx2).norm_square());
scalargradType  normalgb = (nx1-nx2)/(normgradgb+constV(1.0e-16));
scalarvalueType gammagb;
scalargradType dgammadnormalgb;
scalarvalueType gamma_0;
scalargradType  normalnz;
scalargradType  anisogb;
scalargradType  anisogb1;
scalargradType  anisogb2;
double  alpha=M_PI/4.0;
double beta=0.0;
double a=0.2;
if(dim ==3){
normalnz[0][0]=std::cos(alpha)*std::cos(beta)*normalgb[0][0]-std::sin(alpha)*normalgb[1][0]+std::cos(alpha)*std::sin(beta)*normalgb[0][2];
normalnz[1][0]=std::sin(alpha)*std::cos(beta)*normalgb[0][0]+std::cos(alpha)*normalgb[1][0]+std::sin(alpha)*std::sin(beta)*normalgb[0][2];
normalnz[2][0]=-std::sin(beta)*normalgb[0][0]+std::cos(beta)*normalgb[0][2];
gammagb=constV(1.0)*(constV(1.0)-constV(a)*((4.0)*(normalnz[0]*normalnz[0]*normalnz[0]*normalnz[0]+normalnz[1]*normalnz[1]*normalnz[1]*normalnz[1]+normalnz[2]*normalnz[2]*normalnz[2]*normalnz[2])-constV(3.0)));
dgammadnormalgb[0][0]=a*16.0*(normalnz[0][0]*normalnz[0][0]*normalnz[0][0]*std::cos(alpha)*std::cos(beta)+normalnz[1][0]*normalnz[1][0]*normalnz[1][0]*std::sin(alpha)*std::cos(beta)-normalnz[2][0]*normalnz[2][0]*normalnz[2][0]*std::sin(beta));
dgammadnormalgb[1][0]=a*16.0*(-normalnz[0][0]*normalnz[0][0]*normalnz[0][0]*std::sin(alpha)+normalnz[1][0]*normalnz[1][0]*normalnz[1][0]*std::cos(alpha));
dgammadnormalgb[2][0]=a*16.0*(normalnz[0][0]*normalnz[0][0]*normalnz[0][0]*std::cos(alpha)*std::sin(beta)+normalnz[1][0]*normalnz[1][0]*normalnz[1][0]*std::sin(alpha)*std::sin(beta)+normalnz[2][0]*normalnz[2][0]*normalnz[2][0]*std::cos(beta));

}
else{
	normalnz[0][0]=std::cos(alpha)*normalgb[0][0]-std::sin(alpha)*normalgb[1][0];
	normalnz[1][0]=std::sin(alpha)*normalgb[0][0]+std::cos(alpha)*normalgb[1][0];
	gammagb=constV(1.0)*(constV(1.0)-constV(a)*((4.0)*(normalnz[0]*normalnz[0]*normalnz[0]*normalnz[0]+normalnz[1]*normalnz[1]*normalnz[1]*normalnz[1])-constV(3.0)));
	dgammadnormalgb[0][0]=a*16.0*(normalnz[0][0]*normalnz[0][0]*normalnz[0][0]*std::cos(alpha)+normalnz[1][0]*normalnz[1][0]*normalnz[1][0]*std::sin(alpha));
	dgammadnormalgb[1][0]=a*16.0*(-normalnz[0][0]*normalnz[0][0]*normalnz[0][0]*std::sin(alpha)+normalnz[1][0]*normalnz[1][0]*normalnz[1][0]*std::cos(alpha));

}
scalarvalueType phi;
phi[0]=std::atan(normalgb[0][0]/normalgb[1][0])+M_PI/4.0;

anisotropy_gb(normalnz, gammagb, dgammadnormalgb,gamma_0);
//gamagb[0]=(1.0-0.8*std::cos(phi[0]));
scalargradType anison;
for (unsigned int i=0; i<dim; ++i){
		for (unsigned int j=0; j<dim; ++j){
				anisogb[i] += -normalgb[i]*normalgb[j]*dgammadnormalgb[j];
				if (i==j) anisogb[i] +=dgammadnormalgb[j];
		}
}













gammagb=constV(0.5);



 scalarvalueType fcV = constV(2.0/del)*((constV(12.0)*gamma-constV(7.0)*gammagb)*c*(2.0*c*c-3.0*c+1.0)+gammagb*(c+2.0*n1*n1*n1-3.0*n1*n1+2.0*n2*n2*n2-3.0*n2*n2));

	scalarvalueType eq_mu =fcV;
scalargradType eqx_mu = (constV(3.0/4.0*del)*(2.0*gamma*cx-gammagb*cx)-1.0*constV(delta2)*biharmx+(constV(12.0)*c*c*(constV(1.0)-c)*(constV(1.0)-c)+(constV(3.0/4.0)*cx*cx))*(anisox/constV(del)));//aniso-constV(delta2)*biharmx;//-(gammagb)*cx normgradc+c
//+constV(0.5*3.0/4.0)*cx*cx
variable_list.set_scalar_value_term_RHS(4,eq_mu);
variable_list.set_scalar_gradient_term_RHS(4,eqx_mu);

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
