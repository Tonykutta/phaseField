// Number of orientation vectors used to generate anisotropy

#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include "bulatov.h"
template <int dim,int degree>
void customPDE<dim,degree>::anisotropy(const dealii::Tensor<1, dim, dealii::VectorizedArray<double> > &normal1,
dealii::Tensor<1, dim, dealii::VectorizedArray<double> > &normal2,
dealii::VectorizedArray<double> &gamma
) const {
double gammagb;
std::vector<double> v(dim);
double angle=0.0;
//std::vector<double> v_norm(dim);
for(int i=0;i<2;i++){
v[i]= normal1[i][0]-normal2[i][0];

}
//double quad=0.0;
//for(int j=0;j<dim;j++){
//    quad=v[i]*v[i];
//    }

//for(int i=0;i<dim;i++){
//    v_norm[i]=v[i]/quad;
//}
//angle=std::atan(v[1]/v[0]);
//double angle0=M_PI/4.0;
//angle=angle-angle0;
//gamma[0]=std::abs(std::cos(angle))+std::abs(std::sin(angle));
//gamma[0]=(1.0-0.24*std::cos(angle));

/*if(std::abs(v[0])>0.0001 ||std::abs(v[1])>0.0001  ){
bulatov(v,gammagb);
}
else{
  gammagb=0.0;
}*/
gamma[0]=0.6103; //10 degrees

}
