
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <math.h>
 #include "../../include/matrixFreePDE.h"



template <int dim,int degree>
void customPDE<dim,degree>::rotation(std::vector<double>  &Angle, std::vector<double> &RotationVector, std::vector<double> &Tensor, std::vector<double> &TensorRot) const {




  std::vector<double> Angle_rad(3);

  for ( int i=0;i<2;i++ ){
  Angle_rad[i]=Angle[i]*M_PI/180.0;

  }

double phi_1= Angle_rad[0];
double Phi= Angle_rad[1];
double phi_2= Angle_rad[2];

double c1 = std::cos(phi_1);
double c2 = std::cos(Phi);
double c3 = std::cos(phi_2);

double s1 = std::sin(phi_1);
double s2 = std::sin(Phi);
double s3 = std::sin(phi_2);

RotationVector[0] = c1 * c3 - c2 * s1 * s3;  // R11
RotationVector[3] = -c1 * s3 - c2 * c3 * s1; // R12
RotationVector[6] = s1 * s2;                 // R12

RotationVector[1] = c3 * s1 + c1 * c2 * s3; // R21
RotationVector[4] = c1 * c2 * c3 - s1 * s3; // R22
RotationVector[7] = -c1 * s2;               // R23

RotationVector[2] = s2 * s3; // R31
RotationVector[5] = c3 * s2; // R32
RotationVector[8] = c2;      // R33


TensorRot[0]=RotationVector[0]*Tensor[0]+RotationVector[3]*Tensor[1]+RotationVector[6]*Tensor[2];
TensorRot[1]=RotationVector[1]*Tensor[0]+RotationVector[4]*Tensor[1]+RotationVector[7]*Tensor[2];
TensorRot[2]=RotationVector[2]*Tensor[0]+RotationVector[5]*Tensor[1]+RotationVector[8]*Tensor[2];

//std::cout<<"here1:"<<RotationVector[0]<<'\n';

}





template <int dim,int degree>
void customPDE<dim,degree>::rotation1(std::vector<double>  &Angle, std::vector<double> &RotationVector) const {




  std::vector<double> Angle_rad(3);

  for ( int i=0;i<2;i++ ){
  Angle_rad[i]=Angle[i]*M_PI/180.0;

  }

double phi_1= Angle_rad[0];
double Phi= Angle_rad[1];
double phi_2= Angle_rad[2];

double c1 = std::cos(phi_1);
double c2 = std::cos(Phi);
double c3 = std::cos(phi_2);

double s1 = std::sin(phi_1);
double s2 = std::sin(Phi);
double s3 = std::sin(phi_2);

RotationVector[0] = c1 * c3 - c2 * s1 * s3;  // R11
RotationVector[3] = -c1 * s3 - c2 * c3 * s1; // R12
RotationVector[6] = s1 * s2;                 // R13

RotationVector[1] = c3 * s1 + c1 * c2 * s3; // R21
RotationVector[4] = c1 * c2 * c3 - s1 * s3; // R22
RotationVector[7] = -c1 * s2;               // R23

RotationVector[2] = s2 * s3; // R31
RotationVector[5] = c3 * s2; // R32
RotationVector[8] = c2;      // R33




}
