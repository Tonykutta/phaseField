#include <iostream>
#include <Eigen/Dense>
#include <string>
#include <math.h>
#include <Eigen/StdVector>
#include <Eigen/Core>
#include <vector>
#define _USE_MATH_DEFINES
using namespace Eigen;
#include<algorithm>



void   sortGeom(std::vector< std::vector<double> > & geom) {
  //sort values by distance, ksi, eta, and phi
  bool change = true; while (change) {
  change = false;
  for (unsigned int i = 0; i < (geom[0].size() - 1); ++i) {
    if (geom[0][i+1] < geom[0][i]) {
      double temp = geom[0][i];
      geom[0][i] = geom[0][i+1];
      geom[0][i+1] = temp;
      temp = geom[1][i];
      geom[1][i] = geom[1][i+1];
      geom[1][i+1] = temp;
      temp = geom[2][i];
      geom[2][i] = geom[2][i+1];
      geom[2][i+1] = temp;
      temp = geom[3][i]; geom[3][i] = geom[3][i+1];
      geom[3][i+1] = temp;
      change = true;
    }
    else if (geom[0][i+1] == geom[0][i])
    {
      if (geom[1][i+1] < geom[1][i])
      {
        double temp = geom[0][i];
        geom[0][i] = geom[0][i+1];
        geom[0][i+1] = temp;
        temp = geom[1][i];
        geom[1][i] = geom[1][i+1];
        geom[1][i+1] = temp;
        temp = geom[2][i];
        geom[2][i] = geom[2][i+1];
        geom[2][i+1] = temp;
        temp = geom[3][i];
        geom[3][i] = geom[3][i+1];
        geom[3][i+1] = temp;
        change = true;
      }
      else if (geom[1][i+1] == geom[1][i]) {
        if (geom[2][i+1] < geom[2][i]) {
          double temp = geom[0][i];
          geom[0][i] = geom[0][i+1];
          geom[0][i+1] = temp;
          temp = geom[1][i];
          geom[1][i] = geom[1][i+1];
          geom[1][i+1] = temp;
          temp = geom[2][i];
          geom[2][i] = geom[2][i+1];
          geom[2][i+1] = temp;
          temp = geom[3][i];
          geom[3][i] = geom[3][i+1];
          geom[3][i+1] = temp;
          change = true;
        }
        else if (geom[2][i+1] == geom[2][i])
        {
          if (geom[3][i+1] < geom[3][i])
          {
            double temp = geom[0][i];
            geom[0][i] = geom[0][i+1];
            geom[0][i+1] = temp;
            temp = geom[1][i];
            geom[1][i] = geom[1][i+1];
            geom[1][i+1] = temp;
            temp = geom[2][i];
            geom[2][i] = geom[2][i+1];
            geom[2][i+1] = temp;
            temp = geom[3][i];
            geom[3][i] = geom[3][i+1];
            geom[3][i+1] = temp;
            change = true;
          }
        }
      }
    }
  }
}
for (unsigned int i = 0; i < (geom[0].size()-1); ++i) {
  for (unsigned int j = i+1; j < geom[0].size(); ++j) {
    if (geom[0][i] == geom[0][j] && geom[1][i] == geom[1][j] && geom[2][i] == geom[2][j] && geom[3][i] == geom[3][j])
    {
      geom[0].erase(geom[0].begin() + j);
      geom[1].erase(geom[1].begin() + j);
      geom[2].erase(geom[2].begin() + j);
      geom[3].erase(geom[3].begin() + j);
      j--;
    }
  }
}

for (unsigned int i = 0; i < geom[0].size(); ++i) {
  if (geom[0][i]== 0 && geom[1][i] == 0 && geom[2][i] == 0 && geom[3][i] == 0) {
    geom[0].erase(geom[0].begin() + i);
    geom[1].erase(geom[1].begin() + i);
    geom[2].erase(geom[2].begin() + i);
    geom[3].erase(geom[3].begin() + i);
  }
}
}



void  makeparvec(std::string &Material, Eigen::VectorXd &par43 ){
double eRGB;
  VectorXd  par42Al(42);
  par42Al <<0.405204179289160,0.738862004021890,0.351631012630026,2.40065811939667,1.34694439281655,0.352260396651516,0.602137375062785,1.58082498976078,0.596442399566661,1.30981422643602,3.21443408257354,0.893016409093743,0.835332505166333,0.933176738717594,0.896076948651935,0.775053293192055,0.391719619979054,0.782601780600192,0.678572601273508,1.14716256515278,0.529386201144101,0.909044736601838,0.664018011430602,0.597206897283586,0.200371750006251,0.826325891814124,0.111228512469435,0.664039563157148,0.241537262980083,0.736315075146365,0.514591177241156,1.73804335876546,3.04687038671309,1.48989831680317,0.664965104218438,0.495035051289975,0.495402996460658,0.468878130180681,0.836548944799803,0.619285521065571,0.844685390948170,1.02295427618256;
  VectorXd  par42Cu(42);
  par42Cu<<0.405204179289160,0.738862004021890,0.351631012630026,2.40065811939667,1.34694439281655,3.37892632736175,0.602137375062785,1.58082498976078,0.710489498577995,0.737834049784765,3.21443408257354,0.893016409093743,0.835332505166333,0.933176738717594,0.896076948651935,0.775053293192055,0.509781056492307,0.782601780600192,0.762160812499734,1.10473084066580,0.529386201144101,0.909044736601838,0.664018011430602,0.597206897283586,0.200371750006251,0.826325891814124,0.0226010533470218,0.664039563157148,0.297920289861751,0.666383447163744,0.514591177241156,1.73804335876546,2.69805148576400,1.95956771207484,0.948894352912787,0.495035051289975,0.301975031994664,0.574050577702240,0.836548944799803,0.619285521065571,0.844685390948170,0.0491040633104212;


double AlCuparameter;

  if(Material.compare(0,2,"Ni")==0){
    eRGB = 1.44532834613925;
    AlCuparameter = 0.767911805073948;}

    else if(Material.compare(0,2,"Al")==0){

      eRGB = 0.547128733614891;
      AlCuparameter = 0;

    }

    else if(Material.compare(0,2,"Au")==0){
      eRGB = 0.529912885175204;
      AlCuparameter = 0.784289766313152;


    }
    else if(Material.compare(0,2,"Cu")==0){

      eRGB = 1.03669431227427;
      AlCuparameter = 1;

    }

    par43(0)=eRGB;
    par43.tail(42)<<par42Al+AlCuparameter*(par42Cu-par42Al);



  }

  void  quat2mat(Eigen::MatrixXd &m,Eigen::Vector4d &q ){

    double  e0;
    e0=q(0);
    double  e1;
    e1 = q(1);
    double  e2;
    e2 = q(2);
    double  e3;
    e3 = q(3);
    m.resize(3,3);
    m << e0*e0+e1*e1-e2*e2-e3*e3 , 2.0*(e1*e2-e0*e3) , 2.0*(e1*e3+e0*e2),
    2.0*(e1*e2+e0*e3) , e0*e0-e1*e1+e2*e2-e3*e3 , 2.0*(e2*e3-e0*e1),
    2.0*(e1*e3-e0*e2) , 2.0*(e2*e3+e0*e1) , e0*e0-e1*e1-e2*e2+e3*e3;
    m=m/(e0*e0+e1*e1+e2*e2+e3*e3);
  }

  void mat2quat(Eigen::MatrixXd &m,Eigen::Vector4d &q ){
    double t;
    double e0;
    double e1;
    double e3;
    VectorXd e(3);
    t = m(0,0)+m(1,1)+m(2,2);
    e0 = sqrt(1.0+t)/2.0;

    if(t>-0.999999999){
      e << m(1,2)-m(2,1) ,m(2,0)-m(0,2),m(0,1)-m(1,0);
      e=e/(4.0*e0);
    }else
    {
      e0=0.0;
      e3 = sqrt(-(m(0,0)+m(1,1))/2.0);
      if(std::abs(e3)>2.0e-8){
        e<< m(0,2)/(2.0*e3) , m(1,2)/(2.0*e3) , e3;
      }
      else{
        e1 = sqrt((m(0,0)+1.0)/2.0);
        if(e1!=0.0  ){
          e << e1,m(1,0)/(2*e1),0.0;
        }
        else{

          e<<0.0,1.0,0.0;
        }
      }
    }

    q<< e0,-e(0),-e(1),-e(2);
  }

  void rsw(Eigen::VectorXd &theta, double theta1, double theta2, double a,Eigen::VectorXd &en){

  double  dtheta = theta2 - theta1;


   VectorXd sins(theta.size());
    VectorXd xlogx(theta.size());
    for(int i=0; i<xlogx.size();i++){
   xlogx(i)=0.0;

    }

    for(int i=0; i<theta.size(); i++){
           theta(i) = (theta(i)-theta1)/dtheta*M_PI/2.0 ;
          }
  sins=theta.array().sin();

  for (int i=0; i<sins.size(); i++){
  xlogx(i)=sins(i)*std::log(sins(i));
  }

  en = sins-a*xlogx;
  for (int i=0;i<en.size();i++){

    if(isnan(en(i))==1){
      en(i)=0.0;

    }

}
  }


  /*void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove)
  {
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
    matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.rightCols(numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
  }*/

  void distances_to_set(Eigen::MatrixXd &P, Eigen::MatrixXd &Q ,int whichaxes, Eigen::MatrixXd &GEOM ){
    MatrixXd axes(10,10);
    MatrixXd dirs(10,10);
    MatrixXd Qp(10,10);
    Qp=Q;
std::vector<std::vector<double>>  geom;
geom.resize(4, std::vector<double>(10, 0));
    double  dismax = 0.999999;
    int switcher;
    int z=0;
    if(whichaxes==2){

      axes.resize(3,6);
      axes << 1.0 , 1.0 , 1.0,  1.0,  0.0,  0.0 ,
      1.0 ,-1.0 , 0.0 , 0.0 , 1.0 , 1.0 ,
      0.0  ,0.0 , 1.0, -1.0,  1.0, -1.0;
      axes=axes/sqrt(2);
      dirs.resize(3,6);
      dirs << 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 1.0 ,
      0.0  ,0.0  ,1.0 , 1.0,  0.0 , 0.0 ,
      1.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0  ;


    }
    else if(whichaxes==3){

      axes.resize(3,4);
      axes << 1.0 , 1.0 ,-1.0, -1.0 ,
      1.0, -1.0,  1.0 ,-1.0 ,
      1.0, -1.0, -1.0,  1.0 ;

      axes=axes/sqrt(3);
      dirs.resize(3,4);
      dirs <<1.0,  1.0,  1.0,  1.0 ,
      -1.0,  1.0,  1.0, -1.0 ,
      0.0,  0.0,  0.0,  0.0 ;

      dirs=dirs/sqrt(2) ;

    }
    else{


      axes.resize(3,3);
      axes << 1.0,  0.0,  0.0 ,
      0.0,  1.0,  0.0 ,
      0.0,  0.0,  1.0 ;
      dirs.resize(3,3);
      dirs << 0.0,  0.0,  1.0,
      1.0,  0.0,  0.0 ,
      0.0,  1.0,  0.0 ;
    }



    int naxes=axes.cols();
    double period=M_PI*naxes/6.0;




    Matrix3d rotX90;
    Matrix3d rotY90;
    Matrix3d rotZ90;
    Matrix3d rotZ90m;
    rotX90 << 1.0,  0.0,  0.0,
    0.0,  0.0, -1.0,
    0.0,  1.0,  0.0 ;

    rotY90 << 0.0,  0.0,  1.0,
    0.0,  1.0,  0.0 ,
    -1.0,  0.0,  0.0 ;

    rotZ90 <<0.0, -1.0,  0.0,
    1.0,  0.0,  0.0,
    0.0,  0.0,  1.0;

    rotZ90m<< 0.0,  1.0,  0.0 ,
    -1.0,  0.0,  0.0 ,
    0.0,  0.0 , 1.0 ;



    std::vector<Eigen::MatrixXd,Eigen::aligned_allocator<Eigen::MatrixXd> >  V(24);

    V[0]=Qp;
    V[1]  = V[0]*rotX90 ;
    V[2]  = V[1]*rotX90 ;
    V[3]  = V[2]*rotX90 ;

    for ( int j = 0; j<12;j++){
      V[j+4] = V[j]*rotY90 ;
    }

    for ( int j = 0; j<4;j++){

      V[j+16] = V[j]*rotZ90;
      V[j+20] = V[j]*rotZ90m;
    }

    VectorXd distances(24*naxes);
    VectorXd phis(24*naxes);
    VectorXd ksis(24*naxes);
    VectorXd etas(24*naxes);
    geom[0].resize(24*naxes); // distances
    geom[1].resize(24*naxes); // ksis
    geom[2].resize(24*naxes); // etas
    geom[3].resize(24*naxes); // phis
    MatrixXd R(3,3);
    MatrixXd RA(3,3);
    Vector4d q(4);
    Vector3d ax;
    Vector3d dir;
    Vector3d dir2;
    Vector3d axi;
    VectorXd n1(P.cols());
    VectorXd n2(Q.cols());
    VectorXd m1(P.cols());
    VectorXd m2(Q.cols());
    double psi;
    double eta;
    double ksi;
    double dotp;
    double dis;
    double theta;
    double phi;
    Vector4d qa;
    double  theta1;
    double theta2;
    for(int i=0; i<distances.size();i++){
      distances(i)=0.0;
      phis(i)=0.0;
      ksis(i)=0.0;
      etas(i)=0.0;
    }
    int thisindex =0;


    for(int i=0; i<naxes;i++){

      ax = axes.col(i);
      dir=dirs.col(i);
      dir2=ax.cross(dir);

      for(int j=0;j<24;j++){
        Qp=V[j];

        R=Qp.transpose()*P;

        mat2quat(R,q);

        axi=((q.tail(3)).transpose()).normalized();

        psi=2.0*std::acos(q(0));
        dotp=axi.dot(ax);


        dis=2.0*sqrt(std::abs(1.0-dotp*dotp))*std::sin(psi/2.0);

        if(dis<dismax){
          thisindex = thisindex + 1;
          theta = 2.0*std::atan(dotp*std::tan(psi/2.0));

          n1    = (P.row(0)).transpose() ;
          n2    = (Qp.row(0)).transpose() ;

          qa(0)=std::cos(theta/2.0);
          qa.tail(3)=std::sin(theta/2)*ax;

          quat2mat(RA,qa);

          m1    = n1 + RA.transpose()*n2 ;

          m1=m1.normalized();
          m2=RA*m1;


          phi   = std::real(std::acos(std::abs(m1.transpose()*ax)));

          if( std::abs(ax.transpose()*m1) >0.9999){
            theta1= -theta/2.0;
            theta2= theta/2.0;
          }
          else{
            theta1=std::atan2(dir2.transpose()*m1,dir.transpose()*m1);
            theta2=std::atan2(dir2.transpose()*m2,dir.transpose()*m2);
          }



          theta2  = theta2 - round(theta2/period)*period ;
          theta1  = theta1 - round(theta1/period)*period ;

          if (std::abs(theta2+period/2)<0.000001){
            theta2 = theta2 + period;
          }

          if (std::abs(theta1+period/2)<0.000001){
            theta1 = theta1 + period;
          }


          ksi  = std::abs(theta2 - theta1) ;
          eta   = std::abs(theta2 + theta1) ;


          distances(thisindex) = dis;

          ksis(thisindex)      = ksi;
          etas(thisindex)      = eta;

          phis(thisindex)      = phi;



          geom[0][thisindex] = 1.0e-6*std::round(dis*1.0e6);
          geom[1][thisindex] = 1.0e-6*std::round(ksi*1.0e6);
          geom[2][thisindex]= 1.0e-6*std::round(eta*1.0e6);
          geom[3][thisindex] = 1.0e-6*std::round(phi*1.0e6);


        }

    }
}



geom[0].resize(thisindex);
geom[1].resize(thisindex);
geom[2].resize(thisindex);
geom[3].resize(thisindex);

sortGeom(geom);
GEOM.resize(4,geom[0].size());
for(int i=0;i<4;i++){
  for(int j=0; j<geom[i].size();j++){
     GEOM(i,j)=geom[i][j];
  }

}


  /*  VectorXi ind;
    VectorXd sorted_vec;
    int Num=0;
    for(int i=0; i<distances.rows();i++ ){

      if(distances(i)>0.0){
        Num=Num+1;
      }

    }

    ind=VectorXi::LinSpaced(distances.size(),0,distances.size()-1);
    auto rule=[distances](int i, int j)->bool{
      return distances(i)<distances(j);
    };
    std::sort(ind.data(),ind.data()+ind.size(),rule);
    sorted_vec.resize(distances.size());
    for(int i=0;i<distances.size();i++){
      sorted_vec(i)=distances(ind(i));
    }
    distances=sorted_vec.tail(Num);

    for(int i=0;i<ksis.size();i++){
      sorted_vec(i)=ksis(ind(i));
    }
    ksis=sorted_vec.tail(Num);
    for(int i=0;i<etas.size();i++){
      sorted_vec(i)=etas(ind(i));
    }
    etas=sorted_vec.tail(Num);
    for(int i=0;i<phis.size();i++){
      sorted_vec(i)=phis(ind(i));
    }
    phis=sorted_vec.tail(Num);

    for (int i=0;i<distances.size();i++){


      break;
    }

    //geom = unique([distances',ksis',etas',phis'],'rows')';

    MatrixXd geom(4,distances.size());

  /*  geom.row(0)=distances.transpose();
    geom.row(1)=ksis.transpose();
    geom.row(2)=etas.transpose();
    geom.row(3)=phis.transpose();
*/
//std::cout.precision(6);
//std::cout << geom << '\n';
    /*MatrixXd newgeom(geom.rows(),geom.cols());
    newgeom =geom;
    int removed_count = 0;
    for (int i=0;i<geom.cols();i++){
      if (i<geom.cols()-1){
        bool flag = true;
        for (int j=0;j<geom.rows();j++){
          if (flag==true){
            if(geom(j,i) != geom(j,i+1)){
              flag = false;
            }

          }
        }

        if (flag == true){
          removeColumn(newgeom,i-removed_count);
          removed_count+=1;
        }
      }
    }


std::cout << GEOM << '\n';*/
//std::cout << GEOM << '\n';
//GEOM=geom;
}


  void twist100(  Eigen::VectorXd &ksi ,Eigen::VectorXd &   pars, Eigen::VectorXd &en){
     double a;
     a= pars(9);
     double b;
     b= pars(9)*pars(10);
     double perio;
     perio =  M_PI/2.0;
VectorXd vec(ksi.size());
VectorXd sins(vec.size());
VectorXd xlogx(sins.size());
vec=ksi;

    for(int i=0;i<vec.size();i++){
      vec(i)=std::fmod(std::abs(vec(i)),perio);
    }
for(int i=0;i<vec.size();i++){

  if(vec(i)>M_PI/2.0){
    vec(i)=perio-vec(i);
  }
}

sins=((2.0*vec).array()).sin();
xlogx=(sins.array())*(sins.array()).log();
for(int i=0;i<xlogx.size();i++){

  if(isnan(xlogx(i))==1){
    xlogx(i)=0.0;
  }
}

en=a*sins-b*xlogx;



  }


  void twist111(  Eigen::VectorXd &theta ,Eigen::VectorXd &   pars, Eigen::VectorXd &en){

     double thd= pars(36);
     double enm=pars(37);
    double en2=pars(27);
    double a1=pars(35);
    double a2=a1;
     double perio;
     VectorXd vec(theta.size());
      vec=theta;
  for(int i=0;i<vec.size();i++){
    if(vec(i)>M_PI/3.0){
      vec(i)=2.0*M_PI/3.0- vec(i);
    }
  }
  for(int i=0;i<en.size();i++){
    en(i)=0.0;

    }
VectorXd blu(en.size());
VectorXd mul1(en.size());
VectorXd mul2(en.size());
blu=(vec.array()<=thd).select(vec,0.0);
rsw(blu,0.0,thd,a1,mul1);
blu=(vec.array()>thd).select(vec,0.0);
rsw(blu,M_PI/3.0,thd,a2,mul2);
for (int i=0; i<en.size();i++){
  if(vec(i)<=thd){
    en(i)=enm*mul1(i);
  }
  if(vec(i)>thd){
    en(i)=en2+(enm-en2)*mul2(i);
  }
}

  for(int i=0;i<vec.size();i++){
      vec(i)=std::fmod(std::abs(vec(i)),perio);

    }
  }

  void twist110(  Eigen::VectorXd &th ,Eigen::VectorXd &   pars, Eigen::VectorXd &en){
       double th1;
       th1= pars(21);

    double en1,en2,en3;
     en1 = pars(22);
      en2 = pars(23);
       en3 = pars(24);


  double  a01 = 0.5;
  double  a12 = 0.5;
  double  a23 = 0.5;

  double  th2 = std::acos(1.0/3.0) ;
  double  th3 = M_PI/2.0 ;

  double  perio = M_PI ;

  VectorXd vec(th.size());
  vec=th;

  for(int i=0;i<vec.size();i++){
        vec(i)=std::fmod(std::abs(vec(i)),perio);
      }

     for (int i=0; i< vec.size();i++){

       if(vec(i)>perio/2){
         vec(i)=perio-vec(i);
       }
     }


  VectorXd mul1(th.size());
  VectorXd mul2(th.size());
  VectorXd mul3(th.size());
  VectorXd  blu(th.size());
  blu=(vec.array() <=th1).select(vec,0);



  rsw(blu,0.0,th1,a01,mul1);
  blu=(vec.array() >th1 && vec.array() <=th2).select(vec,0);
  rsw(blu,th2,th1,a12,mul2);
  blu=(vec.array() >th2).select(vec,0);
  rsw(blu,th3,th2,a23,mul3);


  for (int i=0; i<th.size(); i++){

      if(th(i)<= th1){
        en(i) = en1*mul1(i) ;
      }
      if(th(i)>th1 && th(i)<=th2){
        en(i) = en2 + (en1-en2)*mul2(i) ;
      }

      if(th(i)>th2){
        en(i) = en3 + (en2-en3)*mul3(i) ;
      }
    }


  }


void stgb100(  Eigen::VectorXd &ksi ,Eigen::VectorXd &   pars, Eigen::VectorXd &en){


  double  en2 = pars(12);
  double   en3 = pars(13);
  double   en4 = pars(14);
  double   en5 = pars(15);
  double   en6 = pars(16);

  double   th2 = pars(17);
  double   th4 = pars(18);


  double th6 = 2.0*std::acos(5.0/std::sqrt(34.0));
  double a12 = 0.5;
  double a23 = a12;
  double a34 = a12;
  double a45 = a12;
  double a56 = a12;
  double a67 = a12;

  double  en1 = 0.0 ;
  double en7 = 0.0 ;

  double th1 = 0.0 ;
  double th3 = std::acos(4.0/5.0) ;
  double th5 = std::acos(3.0/5.0) ;
  double th7 = M_PI/2.0 ;


VectorXd  mul1(ksi.size());
VectorXd  mul2(ksi.size());
VectorXd  mul3(ksi.size());
VectorXd  mul4(ksi.size());
VectorXd  mul5(ksi.size());
VectorXd  mul6(ksi.size());

for (int i=0;i<ksi.size();i++){
  en(i)=0.0;
  mul1(i)=0.0;
    mul2(i)=0.0;
      mul3(i)=0.0;
        mul4(i)=0.0;
          mul5(i)=0.0;
            mul6(i)=0.0;

}


VectorXd  blu(ksi.size());
blu=(ksi.array() <=th2).select(ksi,0);
rsw(blu,th1,th2,a12,mul1);
blu=(ksi.array() >=th2 && ksi.array()<=th3 ).select(ksi,0);

rsw(blu,th3,th2,a23,mul2);
blu=(ksi.array() >=th3 && ksi.array()<=th4 ).select(ksi,0);
rsw(blu,th3,th4,a34,mul3);
blu=(ksi.array() >=th4 && ksi.array()<=th5 ).select(ksi,0);

rsw(blu,th5,th4,a45,mul4) ;
blu=(ksi.array() >=th5 && ksi.array()<=th6 ).select(ksi,0);
rsw(blu,th6,th5,a56,mul5);
blu=(ksi.array() >=th6 && ksi.array()<=th7 ).select(ksi,0);
rsw(blu,th7,th6,a67,mul6);

  for (int i=0; i<ksi.size(); i++){

      if(ksi(i)<=th2){
        en(i) = en1 + (en2-en1)*mul1(i) ;
      }
      if(ksi(i)>=th2 && ksi(i) <= th3){

        en(i) = en3 + (en2-en3)*mul2(i);
      }
      if(ksi(i)>=th3 && ksi(i) <= th4){

        en(i) = en3 + (en4-en3)*mul3(i) ;
      }
      if(ksi(i)>=th4 && ksi(i) <= th5){
        en(i) = en5 + (en4-en5)*mul4(i);
      }
      if(ksi(i)>=th5 && ksi(i) <= th6){
        en(i) = en6 + (en5-en6)*mul5(i);
      }
      if(ksi(i)>=th6 && ksi(i) <= th7){
        en(i) = en7 + (en6-en7)*mul6(i);
      }
  }


}


void stgb110(  Eigen::VectorXd &th ,Eigen::VectorXd &   pars, Eigen::VectorXd &en){




VectorXd vec(th.size());
vec=th;


  double  en2 = pars(26);
  double   en3 = pars(27);
  double   en4 = pars(28);
  double   en5 = pars(29);
  double   en6 = pars(30);

  double   th2 = pars(31);
  double   th4 = pars(32);
  double   th6 = pars(33);



  double a12 = 0.5;
  double a23 = 0.5;
  double a34 = 0.5;
  double a45 = 0.5;
  double a56 = 0.5;
  double a67 = 0.5;

  double  en1 = 0.0 ;
  double en7 = 0.0 ;

  double th1 = 0.0 ;
  double th3 = std::acos(1.0/3.0) ;
  double th5 = std::acos(-7.0/11.0) ;
  double th7 = M_PI;

for(int i=0; i<th.size();i++){
  vec(i)=M_PI-vec(i);
}
VectorXd  mul1(th.size());
VectorXd  mul2(th.size());
VectorXd  mul3(th.size());
VectorXd  mul4(th.size());
VectorXd  mul5(th.size());
VectorXd  mul6(th.size());

for (int i=0;i<th.size();i++){
  en(i)=0.0;
  mul1(i)=0.0;
    mul2(i)=0.0;
      mul3(i)=0.0;
        mul4(i)=0.0;
          mul5(i)=0.0;
            mul6(i)=0.0;

}


VectorXd  blu(vec.size());
blu=(vec.array() <=th2).select(vec,0);
rsw(blu,th1,th2,a12,mul1);
blu=(vec.array() >=th2 && vec.array()<=th3 ).select(vec,0);

rsw(blu,th3,th2,a23,mul2);
blu=(vec.array() >=th3 && vec.array()<=th4 ).select(vec,0);
rsw(blu,th3,th4,a34,mul3);
blu=(vec.array() >=th4 && vec.array()<=th5 ).select(vec,0);

rsw(blu,th5,th4,a45,mul4) ;
blu=(vec.array() >=th5 && vec.array()<=th6 ).select(vec,0);

rsw(blu,th5,th6,a56,mul5);

blu=(vec.array() >=th6 && vec.array()<=th7 ).select(vec,0);

rsw(blu,th7,th6,a67,mul6);

  for (int i=0; i<vec.size(); i++){

      if(vec(i)<=th2){
        en(i) = en1 + (en2-en1)*mul1(i) ;

      }
      if(vec(i)>=th2 && vec(i) <= th3){

        en(i) = en3 + (en2-en3)*mul2(i);
      }
      if(vec(i)>=th3 && vec(i) <= th4){

        en(i) = en3 + (en4-en3)*mul3(i) ;
      }
      if(vec(i)>=th4 && vec(i) <= th5){

        en(i) = en5 + (en4-en5)*mul4(i);
      }
      if(vec(i)>=th5 && vec(i) <= th6){

        en(i) = en5+ (en6-en5)*mul5(i);
      }
      if(vec(i)>=th6 && vec(i) <= th7){

        en(i) = en7 + (en6-en7)*mul6(i);
      }
  }


}




void  atgb111(Eigen::VectorXd &eta, Eigen::VectorXd &ksi, Eigen::VectorXd &  pars, Eigen::VectorXd &en ){


  double    ksim  = pars(38);
  double   enmax = pars(39);
  double   enmin = pars(40);
  double   encnt = pars(41);
  double   a1    = 0.5;
  double   a2    = 0.5;
  double   etascale = pars(42);

  VectorXd blu(en.size());
  VectorXd mul1(en.size());
  VectorXd mul2(en.size());
  VectorXd mul3(en.size());
VectorXd veta(eta.size());
VectorXd vksi(ksi.size());
VectorXd chi(ksi.size());
veta=eta;
vksi=ksi;

for (int i=0;i<veta.size();i++){

  if(veta(i)>M_PI/3.0){
    veta(i)=2.0*M_PI/3-veta(i);
  }
  if(vksi(i)>M_PI/3.0){
    vksi(i)=2.0*M_PI/3-vksi(i);
  }
}
for (int i=0;i<en.size();i++){
  en(i)=0.0;
}






blu=(vksi.array()<=ksim).select(vksi,0.0);

rsw(blu,0.0,ksim,a1,mul1);

blu=(vksi.array()>ksim).select(vksi,0.0);
rsw(blu,M_PI/3.0,ksim,a2,mul2);

blu=(vksi.array()>ksim).select(veta,0.0);
rsw(blu,0.0,M_PI/(2.0*etascale),0.5,mul3);

for (int i=0;i<mul1.size();i++){

  if(isnan(mul1(i))==1){
    mul1(i)=0.0;
  }
  if(isnan(mul2(i))==1){
    mul2(i)=0.0;
  }
  if(isnan(mul3(i))==1){
    mul3(i)=0.0;
  }
}


for (int i=0; i<en.size();i++){
  if(vksi(i)<=ksim){
    en(i)=enmax*mul1(i);
  }
  if(vksi(i)>ksim){
    chi(i)=enmin+(encnt-enmin)*mul3(i);
    en(i)=chi(i)+(enmax-chi(i))*mul2(i);
  }
}





}





void  atgb100(Eigen::VectorXd &eta, Eigen::VectorXd &ksi, Eigen::VectorXd &  pars, Eigen::VectorXd &en ){


  double pwr;
  pwr=pars(11);

  double period;
  period= M_PI/2.0;
  VectorXd en1(ksi.size());
  VectorXd en2(ksi.size());
  VectorXd add(ksi.size());
  for(int i=0; i<ksi.size();i++){
    add(i)=period-ksi(i);
  }

  stgb100(ksi,pars,en1);

  stgb100(add,pars,en2);




  //std::cout << en2 << '\n';
  for (int i=0; i<en.size();i++){
    if (en1(i)>=en2(i)){

      en(i)=en1(i)-(en1(i)-en2(i))*std::pow((eta(i)/period),pwr);
    }
    else{
      en(i)=en2(i)-(en2(i)-en1(i))*std::pow(1.0-(eta(i)/period),pwr);
        }


  }
}


void  atgb110(Eigen::VectorXd &eta, Eigen::VectorXd &ksi, Eigen::VectorXd &  pars, Eigen::VectorXd &en ){


  double a;
  a=pars(25);

  double period;
  period= M_PI;
  VectorXd en1(ksi.size());
  VectorXd en2(ksi.size());
  VectorXd add(ksi.size());

  for(int i=0; i<ksi.size();i++){
    add(i)=period-ksi(i);
  }
for( int i=0; i<ksi.size();i++){

  en1(i)=0.0;
  en2(i)=0.0;
}

  stgb110(ksi,pars,en1);

  stgb110(add,pars,en2);
double zero= 0.0;


VectorXd  mul1(eta.size());
VectorXd  mul2(eta.size());
  for (int i=0;i<mul1.size();i++){
mul1(i)=0.0;
mul2(i)=0.0;
}

  VectorXd  blu(ksi.size());
  blu=(en1.array() >=en2.array() ).select(eta,0);

  rsw(blu,period,zero,a,mul1);
  blu=(en1.array() <en2.array() ).select(eta,0);

  rsw(blu,zero,period,a,mul2);
  for (int i=0;i<mul1.size();i++){

    if(isnan(mul1(i))==1){
      mul1(i)=0.0;
    }
    if(isnan(mul2(i))==1){
      mul2(i)=0.0;
    }
  }


  for (int i=0; i<en.size();i++){
    if (en1(i)>=en2(i)){

      en(i)=en2(i)+(en1(i)-en2(i))*mul1(i);

    }
    else if(en1(i)<en2(i)){
      en(i)=en1(i)+(en2(i)-en1(i))*mul2(i);
        }
      }

}


  void set111(  Eigen::MatrixXd &geom111 ,Eigen::VectorXd &  pars,Eigen::VectorXd &en){
  double a = pars(34);
  double b=a-1.0;

  VectorXd ksi(geom111.cols());
  VectorXd eta(geom111.cols());
  VectorXd phi(geom111.cols());
  VectorXd x(geom111.cols());
  VectorXd entwist(geom111.cols());
  VectorXd entilt(geom111.cols());
  VectorXd temp1(geom111.cols());
  VectorXd temp2(geom111.cols());
  ksi= geom111.row(1);
  eta= geom111.row(2);
  phi= geom111.row(3);


  twist111(ksi,pars,entwist);

  atgb111(eta,ksi,pars,entilt);

  //entilt
for(int i=0; i<en.size();i++){
  en(i)=0.0;
}
  x=phi/(M_PI/2.0);

temp1=(VectorXd::Ones(x.size())-x);

temp2=x;


  en= entwist.array()+(entilt.array()- entwist.array())*(a*x.array()-b*x.array()*x.array());
  }



  void set100(  Eigen::MatrixXd &geom100 ,Eigen::VectorXd &  pars,Eigen::VectorXd &en){
  double pwr1=pars(7);
  double pwr2=pars(8);

  VectorXd ksi(geom100.cols());
  VectorXd eta(geom100.cols());
  VectorXd phi(geom100.cols());
  VectorXd x(geom100.cols());
  VectorXd entwist(geom100.cols());
  VectorXd entilt(geom100.cols());
  VectorXd temp1(geom100.cols());
  VectorXd temp2(geom100.cols());
  ksi= geom100.row(1);
  eta= geom100.row(2);
  phi= geom100.row(3);


  twist100(ksi,pars,entwist);

  atgb100(eta,ksi,pars,entilt);

for(int i=0; i<en.size();i++){
  en(i)=0.0;
}
  x=phi/(M_PI/2.0);

temp1=(VectorXd::Ones(x.size())-x);

temp2=x;

for(int i=0;i<temp1.size();i++){
  temp1(i)=std::pow(temp1(i),pwr1);
  temp2(i)=std::pow(temp2(i),pwr2);
  en(i)= entwist(i)*temp1(i)+entilt(i)*temp2(i);
}

  }



  void set110(  Eigen::MatrixXd &geom110 ,Eigen::VectorXd &  pars,Eigen::VectorXd &en){
  double pwr1=pars(19);
  double pwr2=pars(20);

  VectorXd ksi(geom110.cols());
  VectorXd eta(geom110.cols());
  VectorXd phi(geom110.cols());
  VectorXd x(geom110.cols());
  VectorXd entwist(geom110.cols());
  VectorXd entilt(geom110.cols());
  VectorXd temp1(geom110.cols());
  VectorXd temp2(geom110.cols());
  ksi= geom110.row(1);
  eta= geom110.row(2);
  phi= geom110.row(3);



  twist110(ksi,pars,entwist);

  atgb110(eta,ksi,pars,entilt);

  //entilt
for(int i=0; i<en.size();i++){
  en(i)=0.0;
}
  x=phi/(M_PI/2.0);

temp1=(VectorXd::Ones(x.size())-x);

temp2=x;

for(int i=0;i<en.size();i++){
  temp1(i)=std::pow(temp1(i),pwr1);
  temp2(i)=std::pow(temp2(i),pwr2);
  en(i)= entwist(i)*temp1(i)+entilt(i)*temp2(i);
    }
  }







void weightedmeanenergy( Eigen::MatrixXd &geom100,Eigen::MatrixXd & geom110, Eigen::MatrixXd & geom111,Eigen::VectorXd &   pars, double &en){

double     eRGB = pars(0);
double     d0100 = pars(1);
double     d0110 = pars(2);
double     d0111 = pars(3);
double     weight100 = pars(4);
double     weight110 = pars(5);
double     weight111 = pars(6);
VectorXd  en100(geom100.cols());
VectorXd  en110(geom110.cols());
VectorXd  en111(geom111.cols());
VectorXd  d100(geom100.cols());
VectorXd  d110(geom110.cols());
VectorXd  d111(geom111.cols());
double offset = 0.00001;

set100(geom100,pars,en100);
set110(geom110,pars,en110);
set111(geom111,pars,en111);


d100=geom100.row(0);

d110=geom110.row(0);
d111=geom111.row(0);
VectorXd  s100(d100.size());
VectorXd  s110(d110.size());
VectorXd  s111(d111.size());
VectorXd  w100(d100.size());
VectorXd  w110(d110.size());
VectorXd  w111(d111.size());

s100=(d100.array()/d0100*M_PI/2.0).sin();

for (int i=0; i<d100.size();i++){
if(d100(i)>d0100){
  s100(i)=1.0;
}
}

for (int i=0; i<d100.size();i++){
if(d100(i)<d0100*offset){
  s100(i)=offset*M_PI/2.0;
}
}

w100=(1.0/(s100.array()*(1.0-0.5*s100.array().log()))-1.0)*weight100;

///////////////////
s110=(d110.array()/d0110*M_PI/2.0).sin();

for (int i=0; i<d110.size();i++){
if(d110(i)>d0110){
  s110(i)=1.0;
}
}

for (int i=0; i<d110.size();i++){
if(d110(i)<d0110*offset){
  s110(i)=offset*M_PI/2.0;
}
}

w110=(1.0/(s110.array()*(1.0-0.5*s110.array().log()))-1.0)*weight110;

///////////

s111=(d111.array()/d0111*M_PI/2.0).sin();

for (int i=0; i<d111.size();i++){
if(d111(i)>d0111){
  s111(i)=1.0;
}
}

for (int i=0; i<d111.size();i++){
if(d111(i)<d0111*offset){
  s111(i)=offset*M_PI/2.0;
}
}

w111=(1.0/(s111.array()*(1.0-0.5*s111.array().log()))-1.0)*weight111;


en=eRGB*(((en100.array()*w100.array())).sum()+(en110.array()*w110.array()).sum()+(en111.array()*w111.array()).sum()+1.0)/(w100.array().sum()+w110.array().sum()+w111.array().sum()+1.0);



}

void rotVecToZ(Eigen::Vector3d vec, Eigen::Matrix3d & Tensor){
Vector3d v1(3);
Vector3d w(3);
vec = vec.normalized();

w<<std::abs(vec(0)), std::abs(vec(1)),std::abs(vec(2));
if ((w(2) >= w(1) && w(1) >= w(0)) || (w(1) >= w(2) && w(2) >= w(0)))
       // vec(0) is the smallest component
        v1(0) = 1.0;
     else if ((w(2) >= w(0) && w(0) >= w(1)) || (w(0) >= w(2) && w(2) >= w(1)))
       // vec(1) is the smallest component
       v1(1) = 1.0;
     else
       // vec(2) is the smallest component
     v1(2) = 1.0;

     v1 -= (v1.dot(vec)) * vec;
    v1=v1.normalized();

    Vector3d v0(3);

      v0(0) = v1(1) * vec(2) - v1(2) * vec(1);
      v0(1) = v1(2) * vec(0) - v1(0) * vec(2);
    v0(2) = v1(0) * vec(1) - v1(1) * vec(0);

Tensor<<v0(0), v0(1), v0(2),
        v1(0), v1(1), v1(2),
        vec(0), vec(1), vec(2);


}

void rotVec1toVec2( Eigen::Vector3d & vec1 ,Eigen::Vector3d & vec2 ,Eigen::Matrix3d & Tensor){

  Matrix3d  rot1_to_z(3,3);
  rotVecToZ(vec1,rot1_to_z);
  Matrix3d  rot2_to_z(3,3);
  rotVecToZ(vec2,rot2_to_z);

  Tensor=rot2_to_z.transpose()*rot1_to_z;



}

void RotationMatrix(Eigen::VectorXd  angle, Eigen::MatrixXd & Tensor ){

double Phi_1= angle(0)*M_PI/180.0;
double Phi= angle(1)*M_PI/180.0;
double Phi_2= angle(2)*M_PI/180.0;

double c1= std::cos(Phi_1);
double c2= std::cos(Phi);
double c3= std::cos(Phi_2);

double s1= std::sin(Phi_1);
double s2= std::sin(Phi);
double s3= std::sin(Phi_2);



Tensor.resize(3,3);
Tensor<<  c1 * c3 - c2 * s1 * s3 ,  -c1 * s3 - c2 * c3 * s1 ,  s1 * s2 ,
           c3 * s1 + c1 * c2 * s3,  c1 * c2 * c3 - s1 * s3, -c1 * s2,
           s2 * s3, c3 * s2,   c2;


}


  void bulatov(std::vector<double> &v, double &gamma){

  double en;
  VectorXd par43(43);
  std::string Material= "Cu";
  int or1= 1;
  int or2= 2;
  int or3= 3;
  MatrixXd P(3,3);
  MatrixXd Q(3,3);
  MatrixXd geom100(3,3);
  MatrixXd geom110(3,3);
  MatrixXd geom111(3,3);
  Vector3d vec100;
  Vector3d vecGB;
  Matrix3d ROTGB;

  vec100(0)=1.0;
  vec100(1)=0.0;
  vec100(2)=0.0;
  vecGB(0)=v[0]/v.size();
  vecGB(1)=v[1]/v.size();
  vecGB(2)=0.0;


  rotVec1toVec2(vecGB,vec100,ROTGB);


geom100<<0,0,0,
0,0,0,
0,0,0;
geom110<<0,0,0,
0,0,0,
0,0,0;
geom111<<0,0,0,
0,0,0,
0,0,0;

  Vector3d angle1;
  Vector3d angle2;
  angle1<<45.0,0.0,0;
  angle2<<0.0,0.0,0.0;
 RotationMatrix(angle1, P);
 RotationMatrix(angle2, Q);

P=ROTGB*P;

Q=ROTGB*Q;

/*std::cout << P << '\n';
std::cout << Q << '\n';*/
  /*  P << 0.5774,  0.5774,  0.5774,
    0.7071,  -0.7071, 0.0,
    0.4082,  0.4082,  -0.8165 ;

    Q << 0.5774,  0.5774,  0.5774,
    -0.7071,  0.7071, 0.0,
    -0.4082,  -0.4082,  0.8165 ;*/

    distances_to_set(P, Q, or1,geom100);
  distances_to_set(P, Q, or2,geom110);
    distances_to_set(P, Q, or3,geom111);


    makeparvec(Material, par43);
    weightedmeanenergy( geom100,geom110,geom111,par43,en);

gamma=en/par43(0);

  }
