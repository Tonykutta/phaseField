// Number of orientation vectors used to generate anisotropy

#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>

template <int dim,int degree>
void customPDE<dim,degree>::normtens(const dealii::Tensor<1, dim, dealii::VectorizedArray<double> > &normal1,
dealii::Tensor<1, dim, dealii::VectorizedArray<double> > &normal2,
 dealii::Tensor<1, dim, dealii::VectorizedArray<double> > &normalc,
 dealii::Tensor<2, dim, dealii::VectorizedArray<double> > &normaltensn,
 dealii::Tensor<2, dim, dealii::VectorizedArray<double> > &normaltensc,
 dealii::Tensor<2, dim, dealii::VectorizedArray<double> > &unitens
) const {

std::vector<double> v(dim);
std::vector<double> vc(dim);
std::vector<double> v_norm(dim);
std::vector<double> vc_norm(dim);

double quad_sq=0.0;
double quadc_sq=0.0;


for(int i=0;i<dim;i++){
    unitens[i][i]=constV(1.0);
}

for(int i=0;i<dim;i++){
v[i]= normal1[i][0]-normal2[i][0];
vc[i]=normalc[i][0];
}
double quad=0.0;
double quadc=0.0;
for(int j=0;j<dim;j++){
    quad+=v[j]*v[j];
    quadc+=vc[j]*vc[j];
    }
quad_sq= std::sqrt(quad);
quadc_sq= std::sqrt(quadc);
for(int i=0;i<dim;i++){
    v_norm[i]=v[i]/(quad_sq+1.0e-16);
    vc_norm[i]=vc[i]/(quadc_sq+1.0e-16);
}

for (int i=0; i<dim; i++){
  for (int j=0; j<dim; j++){

  normaltensn[i][j][0]=v_norm[i]*v_norm[j];
  normaltensc[i][j][0]=vc_norm[i]*vc_norm[j];

  }
}

normaltensn = unitens - normaltensn;
normaltensc = unitens - normaltensc;


}
