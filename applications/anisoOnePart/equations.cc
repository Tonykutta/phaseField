#include "anisotropy_facet.h"
#include "anisotropy_gb.h"
#define _USE_MATH_DEFINES

#include <cmath>
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
  set_dependencies_gradient_term_RHS(0, "c, grad(c), n1, grad(n1),grad(mu),biharm,grad(biharm),grad(biharmn1),biharmn1");

    // Variable 1
	set_variable_name				(1,"n1");
	set_variable_type				(1,SCALAR);
	set_variable_equation_type		(1,EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(1, "c, n1");
  set_dependencies_gradient_term_RHS(1, "grad(n1),n1,c,grad(c),grad(biharmn1)");//



	set_variable_name				(2,"biharm");
	set_variable_type				(2,SCALAR);
	set_variable_equation_type		(2,AUXILIARY);

  set_dependencies_value_term_RHS(2, "");
  set_dependencies_gradient_term_RHS(2, "grad(c),mu,grad(mu),biharmn1,grad(biharmn1)");

	set_variable_name				(3,"mu");
	set_variable_type				(3,SCALAR);
	set_variable_equation_type		(3,AUXILIARY);

	set_dependencies_value_term_RHS(3, "c,n1,grad(c),grad(biharm),n1,grad(n1),c");
	set_dependencies_gradient_term_RHS(3, "grad(c),grad(biharm),n1,grad(n1),c");//

		set_variable_name				(4,"biharmn1");
		set_variable_type				(4,SCALAR);
		set_variable_equation_type		(4,AUXILIARY);

		set_dependencies_value_term_RHS(4, "");
		set_dependencies_gradient_term_RHS(4, "grad(n1),mu,grad(mu),biharmn1,n1,grad(biharmn1)");
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
scalargradType nx1 = variable_list.get_scalar_gradient(1);
scalarvalueType  ni=constV(0.0);
scalarvalueType  nj=constV(0.0);
scalargradType nxi;
nxi=nxi*constV(0.0);
scalargradType nxj;
nxj=nxj*constV(0.0);

scalargradType mux = variable_list.get_scalar_gradient(3);
scalargradType biharmn1x = variable_list.get_scalar_gradient(4);



/*unsigned int n=1;

vectorgradType tensn;
vectorgradType tensc;
vectorgradType unity;
scalarvalueType phi;
scalarvalueType ds;
vectorgradType DIFF;
scalarvalueType MOB=constV(0.0);
scalargradType diffmu;
vectorgradType DS;
DS=DS*constV(0.0);
ds=constV(0.0);
phi=c*c*c*(constV(10.0)-constV(15.0)*c+constV(6.0)*c*c);
if (tensor==true){

normtens(nx1,nx1,cx,tensn,tensc,unity);

	DIFF=(constV(D_vol)*phi+constV(D_vap)*(constV(1.0)-phi))*unity+constV(D_surf)*c*(constV(1.0)-c)*tensc;
	diffmu=DIFF*mux;
}else
{

	MOB=(constV(D_vol)*phi+constV(D_vap)*(constV(1.0)-phi))+constV(D_surf)*c*(constV(1.0)-c);
	diffmu=MOB*mux;
}
*/



////////////////////////////////////////////////////////////////////////////////
scalarvalueType normgradn1 = std::sqrt(nx1.norm_square());
scalargradType  normaln1 = nx1/(normgradn1+constV(1.0e-16));
scalarvalueType gamman1;
scalarvalueType gammagb;
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
anison1x=anison1/(normgradn1+constV(1.0e-16));
}


anison1 = (anison1*normgradn1+gamman1*nx1);



////////////////////////////////////////////////////////////////////


/*scalarvalueType normgradgb = std::sqrt((nx1-nx2).norm_square());
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
*/

gammagb=constV(0.5);
scalarvalueType fnV1 = 12.0*n1*(c*(n1-1.0)+n1*n1-2.0*n1+1.0)*gammagb/constV(del);
scalarvalueType eq_c = c;
scalargradType eqx_c = constV(-userInputs.dtValue*McV)*mux;
scalarvalueType eq_n1 = n1-constV(userInputs.dtValue*MnV)*(fnV1);

scalargradType der_n1= (anison1x)*(constV(12.0/del)*(c*c*((constV(1.0)-c)*(constV(1.0)-c)))+constV(3.0/4.0*del)*cx*cx);  //+anisogb*(c*c+constV(6.0)*(constV(1.0)-c)*(n1*n1+n2*n2)-constV(4.0)*(constV(2.0)-c)*(n1*n1*n1+n2*n2*n2)+constV(3.0)*(n1*n1+n2*n2)*(n1*n1+n2*n2));
//normgradn2
scalargradType eqx_n1 = -constV(userInputs.dtValue*MnV)*(constV(3.0/4.0*del)*gammagb*nx1-constV(delta1)*biharmn1x+der_n1);//-anison1  //-anisogb1+   der_n1
variable_list.set_scalar_value_term_RHS(0,eq_c);
variable_list.set_scalar_gradient_term_RHS(0,eqx_c);

variable_list.set_scalar_value_term_RHS(1,eq_n1);
variable_list.set_scalar_gradient_term_RHS(1,eqx_n1);
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
 scalarvalueType n1 = variable_list.get_scalar_value(1);
 scalargradType biharmx = variable_list.get_scalar_gradient(2);
 variable_list.set_scalar_gradient_term_RHS(2,-cx),
 variable_list.set_scalar_gradient_term_RHS(4,-nx1);




 scalarvalueType normgradc = std::sqrt(cx.norm_square());
 scalargradType  normal = cx/(normgradc+constV(1.0e-16));

scalarvalueType gammagb;
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
anisox=aniso/(normgradc+constV(1.0e-16));
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




/*scalarvalueType normgradgb = std::sqrt((nx1-nx2).norm_square());
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


*/










gammagb=constV(0.5);



 scalarvalueType fcV = constV(2.0/del)*((constV(12.0)*gamma-constV(7.0)*gammagb)*c*(2.0*c*c-3.0*c+1.0)+gammagb*(c+2.0*n1*n1*n1-3.0*n1*n1));

	scalarvalueType eq_mu =fcV;
scalargradType eqx_mu = (constV(3.0/4.0*del)*(2.0*gamma*cx-gammagb*cx)-1.0*constV(delta2)*biharmx+(constV(12.0)*c*c*(constV(1.0)-c)*(constV(1.0)-c)+(constV(3.0/4.0)*cx*cx))*(anisox/constV(del)));//aniso-constV(delta2)*biharmx;//-(gammagb)*cx normgradc+c
//+constV(0.5*3.0/4.0)*cx*cx
variable_list.set_scalar_value_term_RHS(3,eq_mu);
variable_list.set_scalar_gradient_term_RHS(3,eqx_mu);

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
