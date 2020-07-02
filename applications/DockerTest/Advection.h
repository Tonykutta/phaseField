
#include <iostream>
#include <cmath>
#include <math.h>
#include "../../include/matrixFreePDE.h"
#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>

template <int dim, int degree>
void  customPDE<dim,degree>::computeIntegralUser(const std::vector<vectorType*> variableSet,  std::vector<double> &Force_x,std::vector<double> &Force_y,std::vector<double> &Volume,std::vector<double> &Center_x,std::vector<double> &Center_y,std::vector<double> &Torque_z) const{

  QGauss<dim>  quadrature_formula(degree+1);
  FE_Q<dim> FE (QGaussLobatto<1>(degree+1));
  FEValues<dim> fe_values (FE, quadrature_formula, update_values | update_JxW_values | update_quadrature_points | update_gradients );
  const unsigned int   n_q_points    = quadrature_formula.size();
  int Ind= Force_x.size();
  std::vector<double> cVal_c(n_q_points);
  std::vector<double> n_val(n_q_points);
  std::vector<double> v_x(n_q_points);
  std::vector<double> v_y(n_q_points);
  std::vector<std::vector<double>> cVal_n(Ind,v_x);
  std::vector<std::vector<double>> cVal_n_x(Ind,v_x);
  std::vector<std::vector<double>> cVal_n_y(Ind,v_x);
  //std::vector<dealii::Tensor<1,dim,double> >  n_grad(n_q_points);

  std::vector<dealii::Tensor<1,dim,double> >  n_grad(n_q_points);
  // constraintsDirichletSet[index]->distribute(*variableSet[index]);
  // constraintsOtherSet[index]->distribute(*variableSet[index]);
  // variableSet[index]->update_ghost_values();

  typename DoFHandler<dim>::active_cell_iterator cell= this->dofHandlersSet[0]->begin_active(), endc = this->dofHandlersSet[0]->end();
  typename DoFHandler<dim>::active_cell_iterator cell1= this->dofHandlersSet[0]->begin_active(), endc1 = this->dofHandlersSet[0]->end();

  std::vector<double> Volume_x(Ind);
  std::vector<double> Volume_y(Ind);


  double pos_x=0.0;
  double pos_y=0.0;
  double k=K;
  double c_0=0.9816;
  double c_gb=0.14;

  double cond;
std::vector<double> vec(Ind);
  //  double sum_grad_eta_z=0.0;

  for (; cell!=endc; ++cell) {
    if (cell->is_locally_owned()){
      fe_values.reinit (cell);

      fe_values.get_function_values(*variableSet[0], cVal_c);
      for (int i=0; i<Ind; i++){
        fe_values.get_function_values(*variableSet[i+2], n_val);



        fe_values.get_function_gradients(*variableSet[i+2], n_grad);


        for (unsigned int p=0;p<n_q_points;p++){
          cVal_n_x[i][p]=n_grad[p][0];
          cVal_n_y[i][p]=n_grad[p][1];
          cVal_n[i][p]=n_val[p];
        }

      }


        for (unsigned int q=0; q<n_q_points; ++q){
          pos_x=fe_values.quadrature_point(q)[0];
          pos_y=fe_values.quadrature_point(q)[1];

        for(int j=0;j<Ind;j++){
          Volume[j]+=(cVal_n[j][q])*fe_values.JxW(q);

          Volume_x[j]+=(cVal_n[j][q])*fe_values.JxW(q)*pos_x;
          Volume_y[j]+=(cVal_n[j][q])*fe_values.JxW(q)*pos_y;

          for (int index=0; index<Ind; index++){
            if (index!=j){
              if(cVal_n[j][q]*cVal_n[index][q]>c_gb){
                cond=1.0;

              }else{

                cond=0.0;
              }

            Force_x[j]+=((cVal_n_x[j][q]-cVal_n_x[index][q])*cond*k*(cVal_c[q]-c_0))*fe_values.JxW(q);

            Force_y[j]+=((cVal_n_y[j][q]-cVal_n_y[index][q])*cond*k*(cVal_c[q]-c_0))*fe_values.JxW(q);
             }
          }

        }
      }
    }

  }



  for (int j=0; j<Ind; j++){

    Center_x[j]=Volume_x[j]/(Volume[j]+1.0e-16);
    Center_y[j]=Volume_y[j]/(Volume[j]+1.0e-16);

  }




  double df_x=0.0;
  double df_y=0.0;




  //////////////////////////////////////////////////




  for (; cell1!=endc1; ++cell1) {
    if (cell1->is_locally_owned()){
      fe_values.reinit (cell1);

      fe_values.get_function_values(*variableSet[0], cVal_c);
      for (int i=0; i<Ind; i++){
        fe_values.get_function_values(*variableSet[i+2], n_val);



        fe_values.get_function_gradients(*variableSet[i+2], n_grad);



        for (unsigned int p=0;p<n_q_points;p++){
          cVal_n_x[i][p]=n_grad[p][0];
          cVal_n_y[i][p]=n_grad[p][1];
          cVal_n[i][p]=n_val[p];
        }


      }

        for (unsigned int q=0; q<n_q_points; ++q){

          pos_x=fe_values.quadrature_point(q)[0];
          pos_y=fe_values.quadrature_point(q)[1];
        for(int j=0;j<Ind;j++){
          Volume[j]+=(cVal_n[j][q])*fe_values.JxW(q);

          Volume_x[j]+=(cVal_n[j][q])*fe_values.JxW(q)*pos_x;
          Volume_y[j]+=(cVal_n[j][q])*fe_values.JxW(q)*pos_y;

          for (int index=0; index<Ind; index++){
            if (index!=j){
              if(cVal_n[j][q]*cVal_n[index][q]>c_gb){
                cond=1.0;

              }else{

                cond=0.0;
              }

            df_x=(cVal_n_x[j][q]-cVal_n_x[index][q])*cond*k*(cVal_c[q]-c_0);

            df_y=(cVal_n_y[j][q]-cVal_n_y[index][q])*cond*k*(cVal_c[q]-c_0);
            Torque_z[j]+=((pos_x-Center_x[j])*df_y-(pos_y-Center_y[j])*df_x)*fe_values.JxW(q);
             }
          }

        }
      }
    }

  }


}
