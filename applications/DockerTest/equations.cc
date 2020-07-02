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
#include "Diff.h"
//#include "bulatov.h"

void variableAttributeLoader::loadVariableAttributes(){
  // Variable 0
  set_variable_name				(0,"c");
  set_variable_type				(0,SCALAR);
  set_variable_equation_type		(0,EXPLICIT_TIME_DEPENDENT);
  set_dependencies_value_term_RHS(0, "n2,n1,c,grad(n2),grad(n1),grad(c),grad(mu),mu");
  set_dependencies_gradient_term_RHS(0, "grad(n2),grad(n1),c,grad(c),n2,n1,grad(mu),mu");

  // Variable 1
  set_variable_name				(1,"mu");
  set_variable_type				(1,SCALAR);
  set_variable_equation_type		(1,AUXILIARY);
  set_dependencies_value_term_RHS(1, "n2,n1,c,grad(n2),grad(n1),grad(c),mu,grad(mu)");
  set_dependencies_gradient_term_RHS(1, "grad(n2),grad(n1),c,grad(c),n2,n1");

  for (unsigned int var_index=1; var_index<=2; var_index++){


    std::string var_name = "n";
    var_name.append(std::to_string(var_index));

    set_variable_name				(var_index+1,var_name);
    set_variable_type				(var_index+1,SCALAR);
    set_variable_equation_type		(var_index+1,EXPLICIT_TIME_DEPENDENT);


    set_dependencies_value_term_RHS(var_index+1, " n1, n2,c, grad(c),mu,grad(mu),grad(n1),grad(n2)");
    set_dependencies_gradient_term_RHS(var_index+1, "n1,n2,grad(n1),grad(n2),c,grad(c)");


  }


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
    unsigned int n;
    n= NumIndex;
    scalarvalueType c = variable_list.get_scalar_value(0);
    scalargradType cx = variable_list.get_scalar_gradient(0);
    scalargradType mux = variable_list.get_scalar_gradient(1);
    scalarvalueType  ni=constV(0.0);
    scalarvalueType  nj=constV(0.0);
    scalarvalueType  nk=constV(0.0);
    scalargradType nxi;
    nxi=nxi*constV(0.0);
    scalargradType nxj;
    nxj=nxj*constV(0.0);
    scalargradType nxk;
    nxk=nxk*constV(0.0);
    scalarvalueType A;
    scalarvalueType B;
    scalarvalueType k_n;
    scalarvalueType gammagb=constV(0.0);
    scalarvalueType gammagb_ex=constV(0.0);
    scalarvalueType advectionterm_c=constV(0.0);
    scalarvalueType advectionterm_n=constV(0.0);
    scalarvalueType fv=constV(0.0);
    scalarvalueType add1=constV(0.0);
    scalarvalueType add2=constV(0.0);
    scalarvalueType add3=constV(0.0);
    scalarvalueType eq_nj;
    scalargradType eqx_nj;
    std::vector<double>  Volume(NumIndex);
    std::vector<double>  Center_x(NumIndex);
    std::vector<double>  Center_y(NumIndex);
    std::vector<double>  Force_x(NumIndex);
    std::vector<double>  Force_y(NumIndex);
    std::vector<double>  Torque_z(NumIndex);
    std::vector<double>  Angle(NumIndex);
    double gamma;
    std::vector<vectorType*> variableSet=this->solutionSet;
    double sumeta=0.0;
    double t= this->currentTime;
    if(this->currentTime >customPDE::IntegrationTimer ||customPDE::InitTimer<userInputs.dtValue){

      if (Advection==true){
      computeIntegralUser(variableSet,Force_x,Force_y,Volume,Center_x,Center_y,Torque_z);
    }

      customPDE::IntegrationTimer=t;

     UpdateAngle(Torque_z, Angle,Volume);
      for (int i=0;i<NumIndex;i++){
        mTorque_z[i]=Torque_z[i];
        mForce_x[i]=Force_x[i];
        mForce_y[i]=Force_y[i];
        mAngle[i]=Angle[i];
        mCenter_x[i]=Center_x[i];
        mCenter_y[i]=Center_y[i];
        mVolume[i]=Volume[i];
      }
      ////////////////////////////////////////////////////////////////////////////
      if(customPDE::InitTimer<userInputs.dtValue){
      // bulatov(mAngle[0],mAngle[1],gamma);
       mGammagb=gamma;
     }
      /////////////////////////////////////////////////////////////////////////////7
      customPDE::InitTimer=10.0*userInputs.dtValue;

  /*    for (int k=0; k<NumIndex; k++){
        for (int j=0; j<NumIndex; j++) {
          if (j>k){
            bulatov(mAngle[k],mAngle[j],gamma);
            mgammaij[k][j]=gamma;
          }
        }
      }
    }
    for (int i=0; i<NumIndex; i++){
      ni=variable_list.get_scalar_value(i+2);
      for (int j=0; j<NumIndex; j++) {
        if (j>i){
          nj=variable_list.get_scalar_value(j+2);
          gammagb_ex[0]+=mgammaij[i][j]*ni[0]*ni[0]*nj[0]*nj[0];
          sumeta+=ni[0]*ni[0]*nj[0]*nj[0];
          mgammaij[j][i]=mgammaij[i][j];
        }
      }*/
    }
//    gammagb[0]=gammagb_ex[0]/(sumeta+1.0e-16);

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
      for (unsigned int i=0; i<n;i++){
        ni=variable_list.get_scalar_value(i+2);
        nxi=variable_list.get_scalar_gradient(i+2);
        for (unsigned int j=0; j<n;j++){
          nj=variable_list.get_scalar_value(j+2);
          nxj=variable_list.get_scalar_gradient(j+2);
          if(j>i){
            normtens(nxi,nxj,cx,tensn,tensc,unity);
            DS+=ni*nj*tensn;
          }
        }
      }
      DIFF=(constV(D_vol)*phi+constV(D_vap)*(constV(1.0)-phi))*unity+constV(D_surf)*c*(constV(1.0)-c)*tensc+constV(D_gb)*DS;
      diffmu=DIFF*mux;
    }else
    {
      for (unsigned int i=0; i<n;i++){
        ni=variable_list.get_scalar_value(i+2);
        for (unsigned int j=0; j<n;j++){
          nj=variable_list.get_scalar_value(j+2);
          if(j>i){
            ds+=ni*nj;
          }
        }
      }
      MOB=(constV(D_vol)*phi+constV(D_vap)*(constV(1.0)-phi))+constV(D_surf)*c*(constV(1.0)-c)+constV(D_gb)*ds;
      diffmu=MOB*mux;
    }

    vectorvalueType v;
      gammagb[0]=0.5;
    B=gammagb;
    k_n=constV(3.0/4.0)*gammagb;
    scalarvalueType dgammadeta=constV(0.0);
    scalarvalueType pos_x;
    pos_x[0]=q_point_loc[0][0];
    scalarvalueType pos_y;
    pos_y[0]=q_point_loc[1][0];

  /*  for (int k=0; k<NumIndex; k++){
      nk=variable_list.get_scalar_value(k+2);
      nxk=variable_list.get_scalar_gradient(k+2);
      add1+=nk*nk;
      add2+=nk*nk*nk;
      add3+=nxk*nxk;
    }*/
scalarvalueType vr_x;
scalarvalueType vr_y;
scalarvalueType vt_x;
scalarvalueType vt_y;
scalarvalueType  v_x;
scalarvalueType  v_y;
scalarvalueType fnV;
scalarvalueType add;
scalargradType s;
std::vector<dealii::VectorizedArray<double>> etas(n);
std::vector<dealii::Tensor<1, dim, dealii::VectorizedArray<double> >  >  n_grad(n,s);

  for(unsigned int j=0;j<n;j++){
etas[j]=  variable_list.get_scalar_value(j+2);
n_grad[j]= variable_list.get_scalar_gradient(j+2);

}



    for(unsigned int j=0;j<n;j++){


      nj = etas[j];
      nxj =n_grad[j];


      if (Advection==true){

      vr_x[0]=-mR/(mVolume[j]+1.0e-16)*(mTorque_z[j])*(pos_y[0]-mCenter_y[j]);

      vr_y[0]=mR/(mVolume[j]+1.0e-16)*(mTorque_z[j])*(pos_x[0]-mCenter_x[j]);

      vt_x[0]=mT*mForce_x[j]/(mVolume[j] +1.0e-16);

      vt_y[0]=mT*mForce_y[j]/(mVolume[j]+1.0e-16);

       v_x=(vr_x+vt_x);
       v_y=(vr_y+vt_y);
      v[0][0]=v_x[0];
      v[1][0]=v_y[0];

        advectionterm_n=(2.0*nj*nxj*v);
        advectionterm_c+=c*(nxj*v)+cx*(nj*v);

      }  else
        {
        advectionterm_n=constV(0.0);
        advectionterm_c=constV(0.0);
      }
        /*if (t<100.0*userInputs.dtValue){
        advectionterm_n=constV(0.0);
        advectionterm_c=constV(0.0);
      }*/


      fnV=constV(0.0);
      add=constV(0.0);
      for (unsigned int i=0; i<n; i++){
        if(i!=j){
          ni=etas[i];
          add+=ni*ni;

        }
      }

      /////////////////////////////////////////////////d_gamma/d_eta//////////////////////////////////////////////////////////////////////
    /*  scalarvalueType one=constV(0.0);
      scalarvalueType two=constV(0.0);
      for( int k=0; k<NumIndex; k++){
        if(k!=j){
          nk=variable_list.get_scalar_value(k+2);
          one+= 2.0*nj*mgammaij[j][k]*nk*nk*sumeta;
          two+=2.0*nj*nk*nk*gammagb_ex;
        }
      }
      dgammadeta= (one-two)/(sumeta*sumeta+constV(1.0e-16));*/
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////7

      //fv=dgammadeta*(-7.0*c*c*(1.0-c)*(1.0-c)+(c*c+6.0*(1.0-c)*add1-4.0*(2.0-c)*add2+3.0*add1*add1)+constV(3.0/4.0*0.5)*add3-constV(3.0/4.0*0.5)*cx*cx);
      fnV=(constV(12.0)*nj*(c*(nj-constV(1.0))+nj*nj-constV(2.0)*nj+constV(1.0)+add))*B;

     eq_nj = nj-constV(userInputs.dtValue)*(constV(MnV)*(fnV)+advectionterm_n);//+fv
      eqx_nj = -constV(userInputs.dtValue*MnV)*k_n*nxj;
      variable_list.set_scalar_value_term_RHS(j+2,eq_nj);
      variable_list.set_scalar_gradient_term_RHS(j+2,eqx_nj);

  }





    scalarvalueType eq_c = c-constV(userInputs.dtValue)*advectionterm_c;
    scalargradType eqx_c = -constV(userInputs.dtValue)*diffmu;
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
      scalarvalueType gammagb=constV(0.0);
      unsigned int n;
      n= NumIndex;
      scalarvalueType c = variable_list.get_scalar_value(0);
      scalargradType cx = variable_list.get_scalar_gradient(0);
      scalarvalueType add=constV(0.0);
      scalarvalueType  ni;
      scalarvalueType  nj;
      for(unsigned int j=0; j<n; j++){
        nj = variable_list.get_scalar_value(j+2);
        add+=constV(2.0)*nj*nj*nj-constV(3.0)*nj*nj;

      }

      scalarvalueType A;
      scalarvalueType B;
      scalarvalueType k_c;
      double sum=0.0;
    /*  for (int i=0; i<NumIndex; i++){
        ni=variable_list.get_scalar_value(i+2);
        for (int j=0; j<NumIndex; j++) {
          if (j>i){
            nj=variable_list.get_scalar_value(j+2);
            gammagb[0]+=mgammaij[i][j]*ni[0]*ni[0]*nj[0]*nj[0];
            sum+=ni[0]*ni[0]*nj[0]*nj[0];
          }
        }
      }
      gammagb[0]=gammagb[0]/(sum+1.0e-16);*/
        gammagb[0]=0.5;
      A=constV(12.0)*constV(gamma_s)-constV(7.0)*gammagb;
      B=gammagb;
      k_c=constV(3.0/4.0)*(constV(2.0)*constV(gamma_s)-gammagb);
      scalarvalueType fcV = constV(2.0)*(A*c*(constV(2.0)*c*c-constV(3.0)*c+constV(1.0))+B*(c+add));


      scalarvalueType eq_mu = fcV;
      scalargradType eqx_mu = k_c*cx;
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
