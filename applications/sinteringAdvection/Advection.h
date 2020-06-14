
#include <iostream>
#include <cmath>
#include <math.h>
#include "../../include/matrixFreePDE.h"

template <int dim, int degree>
void  customPDE<dim,degree>::computeIntegralUser(const std::vector<vectorType*> variableSet) const{


  //  std::cout<<variableSet[0][0].size()<<std::endl;
  QGauss<dim>  quadrature_formula(degree+1);
  FE_Q<dim> FE (QGaussLobatto<1>(degree+1));
  FEValues<dim> fe_values (FE, quadrature_formula, update_values | update_JxW_values | update_quadrature_points | update_gradients );
  const unsigned int   n_q_points    = quadrature_formula.size();
  std::vector<double> cVal_c(n_q_points);
  std::vector<double> cVal_n1(n_q_points);
  std::vector<double> cVal_n2(n_q_points);
  std::vector<double> cVal_n3(n_q_points);
  std::vector<dealii::Tensor<1,dim,double> >  cVal_n1_grad(n_q_points);
  std::vector<dealii::Tensor<1,dim,double> >  cVal_n2_grad(n_q_points);
  std::vector<dealii::Tensor<1,dim,double> >  cVal_n3_grad(n_q_points);
  // constraintsDirichletSet[index]->distribute(*variableSet[index]);
  // constraintsOtherSet[index]->distribute(*variableSet[index]);
  // variableSet[index]->update_ghost_values();

  typename DoFHandler<dim>::active_cell_iterator cell= this->dofHandlersSet[0]->begin_active(), endc = this->dofHandlersSet[0]->end();
  typename DoFHandler<dim>::active_cell_iterator cell1= this->dofHandlersSet[0]->begin_active(), endc1 = this->dofHandlersSet[0]->end();
  double Volume1 = 0.0;
  double Volume1_x = 0.0;
  double Volume1_y = 0.0;
  double Volume2 = 0.0;
  double Volume2_x = 0.0;
  double Volume2_y = 0.0;
  double Volume3 = 0.0;
  double Volume3_x = 0.0;
  double Volume3_y = 0.0;
  double pos_x=0.0;
  double pos_y=0.0;
  double k=K;
  double c_0=0.9816;
  double c_gb=0.14;
  double force1_x=0.0;
  double force1_y=0.0;
  double force2_x=0.0;
  double force2_y=0.0;
  double force3_x=0.0;
  double force3_y=0.0;
  double cond;
  double cond1;
  double cond2;
  double sum_grad_eta_x=0.0;
  double sum_grad_eta_y=0.0;

  //  double sum_grad_eta_z=0.0;

  for (; cell!=endc; ++cell) {
    if (cell->is_locally_owned()){
      fe_values.reinit (cell);

      fe_values.get_function_values(*variableSet[0], cVal_c);
      fe_values.get_function_values(*variableSet[2], cVal_n1);
      fe_values.get_function_values(*variableSet[3], cVal_n2);
      fe_values.get_function_values(*variableSet[4], cVal_n3);
      fe_values.get_function_gradients(*variableSet[2], cVal_n1_grad);
      fe_values.get_function_gradients(*variableSet[3], cVal_n2_grad);
      fe_values.get_function_gradients(*variableSet[4], cVal_n3_grad);
      for (unsigned int q=0; q<n_q_points; ++q){

        pos_x=fe_values.quadrature_point(q)[0];
        pos_y=fe_values.quadrature_point(q)[1];
        //  pos_z=fe_values.quadrature_point(q)[2];
        Volume1+=(cVal_n1[q])*fe_values.JxW(q);
        Volume1_x+=(cVal_n1[q])*fe_values.JxW(q)*pos_x;
        Volume1_y+=(cVal_n1[q])*fe_values.JxW(q)*pos_y;
        //  Volume1_z+=(cVal_n1[q])*fe_values.JxW(q)*pos_z;
        Volume2+=(cVal_n2[q])*fe_values.JxW(q);
        Volume2_x+=(cVal_n2[q])*fe_values.JxW(q)*pos_x;
        Volume2_y+=(cVal_n2[q])*fe_values.JxW(q)*pos_y;
        //  Volume2_z+=(cVal_n2[q])*fe_values.JxW(q)*pos_z;
        Volume3+=(cVal_n3[q])*fe_values.JxW(q);
        Volume3_x+=(cVal_n3[q])*fe_values.JxW(q)*pos_x;
        Volume3_y+=(cVal_n3[q])*fe_values.JxW(q)*pos_y;

        if(cVal_n1[q]*cVal_n2[q]>c_gb){
          cond=1.0;
        }
        else{
          cond=0.0;
        }
        if(cVal_n1[q]*cVal_n3[q]>c_gb){
          cond1=1.0;
        }
        else{
          cond1=0.0;
        }
        if(cVal_n2[q]*cVal_n3[q]>c_gb){
          cond2=1.0;
        }
        else{
          cond2=0.0;
        }


        //std::cout<<"2--:"<<fe_values.JxW(q)<<'\n';
        force1_x+=((cVal_n1_grad[q][0]-cVal_n2_grad[q][0])*cond*k*(cVal_c[q]-c_0)+(cVal_n1_grad[q][0]-cVal_n3_grad[q][0])*cond1*k*(cVal_c[q]-c_0))*fe_values.JxW(q);
        force1_y+=((cVal_n1_grad[q][1]-cVal_n2_grad[q][1])*cond*k*(cVal_c[q]-c_0)+(cVal_n1_grad[q][1]-cVal_n3_grad[q][1])*cond1*k*(cVal_c[q]-c_0))*fe_values.JxW(q);

        force2_x+=((cVal_n2_grad[q][0]-cVal_n1_grad[q][0])*cond*k*(cVal_c[q]-c_0)+(cVal_n2_grad[q][0]-cVal_n3_grad[q][0])*cond2*k*(cVal_c[q]-c_0))*fe_values.JxW(q);
        force2_y+=((cVal_n2_grad[q][1]-cVal_n1_grad[q][1])*cond*k*(cVal_c[q]-c_0)+(cVal_n2_grad[q][1]-cVal_n3_grad[q][1])*cond2*k*(cVal_c[q]-c_0))*fe_values.JxW(q);

        force3_x+=((cVal_n3_grad[q][0]-cVal_n1_grad[q][0])*cond1*k*(cVal_c[q]-c_0)+(cVal_n3_grad[q][0]-cVal_n2_grad[q][0])*cond2*k*(cVal_c[q]-c_0))*fe_values.JxW(q);
        force3_y+=((cVal_n3_grad[q][1]-cVal_n1_grad[q][1])*cond1*k*(cVal_c[q]-c_0)+(cVal_n3_grad[q][1]-cVal_n2_grad[q][1])*cond2*k*(cVal_c[q]-c_0))*fe_values.JxW(q);



      }
    }
  }


  force1_x=Utilities::MPI::sum(force1_x, MPI_COMM_WORLD);
  force1_y=Utilities::MPI::sum(force1_y, MPI_COMM_WORLD);

  force2_x=Utilities::MPI::sum(force2_x, MPI_COMM_WORLD);
  force2_y=Utilities::MPI::sum(force2_y, MPI_COMM_WORLD);

  double center1_x=Volume1_x/(Volume1+1.0e-16);
  double center1_y=Volume1_y/(Volume1+1.0e-16);
  double center2_x=Volume2_x/(Volume2+1.0e-16);
  double center2_y=Volume2_y/(Volume2+1.0e-16);
  double center3_x=Volume3_x/(Volume3+1.0e-16);
  double center3_y=Volume3_y/(Volume3+1.0e-16);

  center1_x=Utilities::MPI::sum(center1_x, MPI_COMM_WORLD);
  center1_y=Utilities::MPI::sum(center1_y, MPI_COMM_WORLD);

  center2_x=Utilities::MPI::sum(center2_x, MPI_COMM_WORLD);
  center2_y=Utilities::MPI::sum(center2_y, MPI_COMM_WORLD);

  center3_x=Utilities::MPI::sum(center3_x, MPI_COMM_WORLD);
  center3_y=Utilities::MPI::sum(center3_y, MPI_COMM_WORLD);

  double torque1_x=0.0;
  double torque1_y=0.0;
  double torque2_x=0.0;
  double torque2_y=0.0;
  double torque1_z=0.0;
  double torque2_z=0.0;

  double torque3_x=0.0;
  double torque3_y=0.0;
  double torque3_z=0.0;
  double df1_x=0.0;
  double df1_y=0.0;
  double df2_x=0.0;
  double df2_y=0.0;

  double df3_x=0.0;
  double df3_y=0.0;

  for (; cell1!=endc1; ++cell1) {
    if (cell1->is_locally_owned()){
      fe_values.reinit (cell1);

      fe_values.get_function_values(*variableSet[0], cVal_c);
      fe_values.get_function_values(*variableSet[2], cVal_n1);
      fe_values.get_function_values(*variableSet[3], cVal_n2);
      fe_values.get_function_values(*variableSet[4], cVal_n3);
      fe_values.get_function_gradients(*variableSet[2], cVal_n1_grad);
      fe_values.get_function_gradients(*variableSet[3], cVal_n2_grad);
      fe_values.get_function_gradients(*variableSet[4], cVal_n3_grad);
      for (unsigned int q=0; q<n_q_points; ++q){
        pos_x=fe_values.quadrature_point(q)[0];
        pos_y=fe_values.quadrature_point(q)[1];
        if(cVal_n1[q]*cVal_n2[q]>c_gb){
          cond=1.0;
        }
        else{
          cond=0.0;
        }
        if(cVal_n1[q]*cVal_n3[q]>c_gb){
          cond1=1.0;
        }
        else{
          cond1=0.0;
        }
        if(cVal_n2[q]*cVal_n3[q]>c_gb){
          cond2=1.0;
        }
        else{
          cond2=0.0;
        }
        df1_x=(cVal_n1_grad[q][0]-cVal_n2_grad[q][0])*cond*k*(cVal_c[q]-c_0)+(cVal_n1_grad[q][0]-cVal_n3_grad[q][0])*cond1*k*(cVal_c[q]-c_0);
        df1_y=(cVal_n1_grad[q][1]-cVal_n2_grad[q][1])*cond*k*(cVal_c[q]-c_0)+(cVal_n1_grad[q][1]-cVal_n3_grad[q][1])*cond1*k*(cVal_c[q]-c_0);

        df2_x=(cVal_n2_grad[q][0]-cVal_n1_grad[q][0])*cond*k*(cVal_c[q]-c_0)+(cVal_n2_grad[q][0]-cVal_n3_grad[q][0])*cond2*k*(cVal_c[q]-c_0);
        df2_y=(cVal_n2_grad[q][1]-cVal_n1_grad[q][1])*cond*k*(cVal_c[q]-c_0)+(cVal_n2_grad[q][1]-cVal_n3_grad[q][1])*cond2*k*(cVal_c[q]-c_0);

        df3_x=(cVal_n3_grad[q][0]-cVal_n1_grad[q][0])*cond1*k*(cVal_c[q]-c_0)+(cVal_n3_grad[q][0]-cVal_n2_grad[q][0])*cond2*k*(cVal_c[q]-c_0);
        df3_y=(cVal_n3_grad[q][1]-cVal_n1_grad[q][1])*cond1*k*(cVal_c[q]-c_0)+(cVal_n3_grad[q][1]-cVal_n2_grad[q][1])*cond2*k*(cVal_c[q]-c_0);
        torque1_z+=((pos_x-center1_x)*df1_y-(pos_y-center1_y)*df1_x)*fe_values.JxW(q);
        torque2_z+=((pos_x-center2_x)*df2_y-(pos_y-center2_y)*df2_x)*fe_values.JxW(q);
        torque3_z+=((pos_x-center3_x)*df3_y-(pos_y-center3_y)*df3_x)*fe_values.JxW(q);

      }
    }
  }

  torque1_z=Utilities::MPI::sum(torque1_z, MPI_COMM_WORLD);
  torque2_z=Utilities::MPI::sum(torque2_z, MPI_COMM_WORLD);
  torque3_z=Utilities::MPI::sum(torque3_z, MPI_COMM_WORLD);

  customPDE::mForce1_x=force1_x;
  customPDE::mForce2_y=force2_y;
  customPDE::mForce1_y=force1_y;
  customPDE::mForce2_x=force2_x;
  customPDE::mForce3_x=force3_x;
  customPDE::mForce3_y=force3_y;
  customPDE::mVolume1=Volume1;
  customPDE::mVolume2=Volume2;
  customPDE::mVolume3=Volume3;
  customPDE::mTorque1_z=torque1_z;
  customPDE::mTorque2_z=torque2_z;
  customPDE::mTorque3_z=torque3_z;
  customPDE::mCenter1_x=center1_x;
  customPDE::mCenter1_y=center1_y;
  customPDE::mCenter2_x=center2_x;
  customPDE::mCenter2_y=center2_y;
  customPDE::mCenter3_x=center3_x;
  customPDE::mCenter3_y=center3_y;


}
