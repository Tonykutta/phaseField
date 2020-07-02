// Number of orientation vectors used to generate anisotropy

#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <math.h>
 #include "../../include/matrixFreePDE.h"
#include "rotation.h"


template <int dim,int degree>
void customPDE<dim,degree>::UpdateAngle(std::vector<double> &Torque_z, std::vector<double> &Angle, std::vector<double> &Volume) const {


int Ind= Torque_z.size();
/*std::vector<double> T1={ customPDE::mTorque1_x, customPDE::mTorque1_y, customPDE::mTorque1_z};
std::vector<double> T2={customPDE::mTorque2_x, customPDE::mTorque2_y,customPDE::mTorque2_z};
std::vector<double> T3={ customPDE::mTorque3_x, customPDE::mTorque3_y, customPDE::mTorque3_z};
//std::vector<double> T4={ customPDE::mTorque4_x, customPDE::mTorque4_y, customPDE::mTorque4_z};

std::vector<double> AngleOld1={ customPDE::mAngleX1, customPDE::mAngleY1,customPDE::mAngleZ1};
std::vector<double> AngleOld2={ customPDE::mAngleX2, customPDE::mAngleY2,customPDE::mAngleZ3};
std::vector<double> AngleOld3={customPDE::mAngleX3, customPDE::mAngleY3, customPDE::mAngleZ3};
//std::vector<double> AngleOld4={customPDE::mAngleX4, customPDE::mAngleY4, customPDE::mAngleZ4};*/
/*std::vector<double> omega(3);
std::vector<double>  angle_change(3);
std::vector<double> Angle(3);
//std::cout<<"geg"<<AngleOld[2]<<'\n';
/*
std::vector<double> RotationVector(9);
//std::vector<double> RotationVector2(9);

std::vector<double> T1_rot(2);

//std::vector<double> T2_rot(2);

rotation(AngleOld, RotationVector,T1,T1_rot);
//rotation(Angle2, RotationVector2,T2,T2_rot);


std::vector<double> omega(2);



for (int i=0;i<3;i++){
omega[i]= mR/customPDE::mVolume1*T1_rot[i];
//omega2[i]= mR/customPDE::mVolume2*T2_rot[i];
}

/*std::cout<<omega[0]<<'\n';
std::cout<<omega[1]<<'\n';
std::cout<<omega[2]<<'\n';*/

//std::cout<<"1--:"<< omega[0]<<'\n';
//std::cout<<"omega2--:"<< omega[2]<<'\n';
/*std::vector<double>  angle_change(3);

angle_change[0] = omega[2] * userInputs.dtValue;
angle_change[1] = (omega[0] * std::cos(angle_change[0]) + omega[1] * std::sin(angle_change[0])) *userInputs.dtValue;
angle_change[2] = ((omega[0] * std::sin(angle_change[0]) * std::sin(angle_change[1]) -
                        omega[1] * std::cos(angle_change[0]) * std::sin(angle_change[1]) +
                        omega[2] * std::cos(angle_change[1]))) *userInputs.dtValue;


for (int i=0; i<3;i++) {
     angle_change[i] = angle_change[i]*(180.0 / M_PI);
}

//
//std::cout<<"1s--:"<<  M_PI<<'\n';
//std::cout<<"2s--:"<<  angle_change[1]<<'\n';
//std::cout<<"3s--:"<<  angle_change[2]<<'\n';

std::vector<double>  RotationVectorChange(9);
std::vector<double> NewRot(9);
rotation1(angle_change, RotationVectorChange);
//R11
NewRot[0]= RotationVectorChange[0]*RotationVector[0]+RotationVectorChange[3]*RotationVector[1]+RotationVectorChange[6]*RotationVector[2];
//R12
NewRot[1]= RotationVectorChange[0]*RotationVector[3]+RotationVectorChange[3]*RotationVector[4]+RotationVectorChange[6]*RotationVector[5];
//R13
NewRot[2]= RotationVectorChange[0]*RotationVector[6]+RotationVectorChange[3]*RotationVector[7]+RotationVectorChange[6]*RotationVector[8];
//R21
NewRot[3]= RotationVectorChange[1]*RotationVector[0]+RotationVectorChange[4]*RotationVector[1]+RotationVectorChange[7]*RotationVector[2];
//R22
NewRot[4]= RotationVectorChange[1]*RotationVector[3]+RotationVectorChange[4]*RotationVector[4]+RotationVectorChange[7]*RotationVector[5];
//R23
NewRot[5]= RotationVectorChange[1]*RotationVector[6]+RotationVectorChange[4]*RotationVector[7]+RotationVectorChange[7]*RotationVector[8];
//R31
NewRot[6]= RotationVectorChange[2]*RotationVector[0]+RotationVectorChange[5]*RotationVector[1]+RotationVectorChange[8]*RotationVector[2];
//R32
NewRot[7]= RotationVectorChange[2]*RotationVector[3]+RotationVectorChange[5]*RotationVector[4]+RotationVectorChange[8]*RotationVector[5];
//R33
NewRot[8]= RotationVectorChange[2]*RotationVector[6]+RotationVectorChange[5]*RotationVector[7]+RotationVectorChange[8]*RotationVector[8];

//std::cout<<NewRot[1]<<'\n';
//

//std::cout<<NewRot[4]<<'\n';

std::vector<double> Angle(3);

/*std::cout<<"1--:"<< RotationVectorChange[0]<<'\n';
std::cout<<"2--:"<< RotationVector[2]<<'\n';
std::cout<<"3--:"<< RotationVectorChange[0]*RotationVector[2]<<'\n';

std::cout<<"4--:"<< RotationVectorChange[1]<<'\n';
std::cout<<"5--:"<< RotationVector[5]<<'\n';
std::cout<<"6--:"<< RotationVectorChange[1]*RotationVector[5]<<'\n';

std::cout<<"7--:"<< RotationVectorChange[2]<<'\n';
std::cout<<"8--:"<< RotationVector[8]<<'\n';
std::cout<<"9--:"<< RotationVectorChange[2]*RotationVector[8]<<'\n';
std::cout<<"10--:"<< RotationVectorChange[2]*RotationVector[8]+RotationVectorChange[1]*RotationVector[5]+RotationVectorChange[0]*RotationVector[2]<<'\n';
std::cout<<"11--:"<< NewRot[2]<<'\n';
std::cout<<"12--:"<< NewRot[5]<<'\n';*/
  /*   if (NewRot[8] != 1.0 && NewRot[8] != -1.0) // checks if cos(Phi) = 1 or -1
      {
         Angle[0] = std::atan2(NewRot[6], -NewRot[7]) * (180.0 /M_PI);
         Angle[1] = std::acos(NewRot[8]) * (180.0 /M_PI);
         Angle[2]= std::atan2(NewRot[2], NewRot[5]) * (180.0 / M_PI);
std::cout<<"here"<<'\n';
       }
       else if (NewRot[8] == 1.0) // special case for Phi = 0.0
         {
            if (RotationVector[8] == 1.0)
              // when Phi_old = 0.0; all the rotations are about z axis and angles accumulates after each
              // rotation
              Angle[0] = AngleOld[0] + AngleOld[2] + angle_change[0];
            else
              Angle[0] = angle_change[0];
            // Comply with bunge euler angle definitions, 0.0 <= phi1 <= 360.0
            if (std::abs(Angle[0]) > 360.0)
            {
              int laps = Angle[0] / 360.0;
              Angle[0] -= laps * 360.0;
           }
            Angle[1] = 0.0;
            Angle[2] = -Angle[0] + std::atan2(NewRot[1], NewRot[4]) * (180.0 / M_PI);

          }
         else
          {
            if (RotationVector[8] == 1.0)
              Angle[0] = AngleOld[0] + Angle[2] + angle_change[0];
            else
              Angle[0] = angle_change[0];
            // Comply with bunge euler angle definitions, 0.0 <= phi1 <= 360.0
            if (std::abs(Angle[0]) > 360.0)
            {
              int laps = Angle[0] / 360.0;
              Angle[0] -= laps * 360.0;
            }
           Angle[1] = 180.0;
            Angle[2] = Angle[0] - std::atan2(-NewRot[1], -NewRot[4]) * (180.0 / M_PI);

          }

      if (Angle[0] < 0.0){
        Angle[0] += 360.0;
      }
     if (Angle[2] < 0.0){
        Angle[2]+= 360.0;
      }
      if (Angle[1] < 0.0){
      std::cout<<"error"<<'\n';
    }*/
    double omega;
    double angle;
    double angleOld;

    for (int j=0; j<Ind; j++){
      angleOld=mAngle[j];
      omega= mR/(Volume[j]+1e-16)*Torque_z[j]*180.0/M_PI;

      angle=omega*userInputs.dtValue+angleOld;
      if(angle<0.0){
        angle=angle+360.0;
      }
      Angle[j]=angle;

    }




/*for (int i=0;i<3;i++){
omega[i]= mR/(customPDE::mVolume4*T4[i]+1e-16)*180/M_PI;
//omega2[i]= mR/customPDE::mVolume2*T2_rot[i];
}
Angle[0]=omega[2]*userInputs.dtValue+AngleOld4[0];
customPDE::mAngleX4=Angle[0];
customPDE::mAngleY4=Angle[1];
customPDE::mAngleZ4=Angle[2];




//customPDE::mAngleX1=Utilities::MPI::sum(customPDE::mAngleX1, MPI_COMM_WORLD);
//customPDE::mAngleY1=Utilities::MPI::sum(customPDE::mAngleY1, MPI_COMM_WORLD);
//customPDE::mAngleZ1=Utilities::MPI::sum(customPDE::mAngleZ1, MPI_COMM_WORLD);
/*std::cout<<"T1:"<< T1[2]<<'\n';
std::cout<<"omega:"<< omega[2]<<'\n';
std::cout<<"Change:"<< omega[2]*userInputs.dtValue<<'\n';


std::cout<<"2--:"<< Angle[0]<<'\n';*/

}
