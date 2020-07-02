// ===========================================================================
// FUNCTION FOR INITIAL CONDITIONS
// ===========================================================================

template <int dim, int degree>
void customPDE<dim,degree>::setInitialCondition(const dealii::Point<dim> &p, const unsigned int index, double & scalar_IC, dealii::Vector<double> & vector_IC){
    // ---------------------------------------------------------------------
    // ENTER THE INITIAL CONDITIONS HERE
    // ---------------------------------------------------------------------
    // Enter the function describing conditions for the fields at point "p".
    // Use "if" statements to set the initial condition for each variable
    // according to its variable index
	  // The initial condition is two circles/spheres defined
	  // by a hyperbolic tangent function. The center of each circle/sphere is
	  // given by "center" and its radius is given by "rad".

      double dist;
	  scalar_IC = 0.0;


      dist = 0.0;
      for (unsigned int dir = 0; dir < dim; dir++){
          dist += (p[dir]-center1[dir])*(p[dir]-center1[dir]);
      }
      dist = std::sqrt(dist);

      // Initial condition for the concentration field
      if (index == 0){
          scalar_IC += 0.5*(1.0-std::tanh((dist-radius1)/(1.0)));
      }
      else if (index==2){
         scalar_IC += 0.5*(1.0-std::tanh((dist-radius1)/(1.0)));
     }

/////////////////////////////////////////////
      dist = 0.0;
      for (unsigned int dir = 0; dir < dim; dir++){
          dist += (p[dir]-center2[dir])*(p[dir]-center2[dir]);
      }
      dist = std::sqrt(dist);

      // Initial condition for the concentration field
      if (index == 0){
          scalar_IC += 0.5*(1.0-std::tanh((dist-radius2)/(1.0)));
      }
      else if (index==3){
         scalar_IC += 0.5*(1.0-std::tanh((dist-radius2)/(1.0)));
     }
//////////////////////////////
  /*   dist = 0.0;
     for (unsigned int dir = 0; dir < dim; dir++){
         dist += (p[dir]-center3[dir])*(p[dir]-center3[dir]);
     }
     dist = std::sqrt(dist);

     // Initial condition for the concentration field
     if (index == 0){
         scalar_IC += 0.5*(1.0-std::tanh((dist-radius3)/(1.0)));
     }
     else if (index==4){
        scalar_IC += 0.5*(1.0-std::tanh((dist-radius3)/(1.0)));
    }

///////////////////////////////////

    dist = 0.0;
    for (unsigned int dir = 0; dir < dim; dir++){
        dist += (p[dir]-center4[dir])*(p[dir]-center4[dir]);
    }
    dist = std::sqrt(dist);

    // Initial condition for the concentration field
    if (index == 0){
        scalar_IC += 0.5*(1.0-std::tanh((dist-radius4)/(1.0)));
    }
    else if (index==5){
       scalar_IC += 0.5*(1.0-std::tanh((dist-radius4)/(1.0)));
   }

//////////////////////////////////////////////////

   dist = 0.0;
   for (unsigned int dir = 0; dir < dim; dir++){
       dist += (p[dir]-center5[dir])*(p[dir]-center5[dir]);
   }
   dist = std::sqrt(dist);

   // Initial condition for the concentration field
   if (index == 0){
       scalar_IC += 0.5*(1.0-std::tanh((dist-radius5)/( 1.0)));
   }
   else if (index==6){
      scalar_IC += 0.5*(1.0-std::tanh((dist-radius5)/(1.0)));
  }
///////////////////////////////////////////////////


/*  dist = 0.0;
  for (unsigned int dir = 0; dir < dim; dir++){
      dist += (p[dir]-center6[dir])*(p[dir]-center6[dir]);
  }
  dist = std::sqrt(dist);

  // Initial condition for the concentration field
  if (index == 0){
      scalar_IC += 0.5*(1.0-std::tanh((dist-radius6)/(1.0)));
  }
  else if (index==7){
     scalar_IC += 0.5*(1.0-std::tanh((dist-radius6)/(1.0)));
 }
/////////////////////////////////////////////////////////
 dist = 0.0;
 for (unsigned int dir = 0; dir < dim; dir++){
     dist += (p[dir]-center7[dir])*(p[dir]-center7[dir]);
 }
 dist = std::sqrt(dist);

 // Initial condition for the concentration field
 if (index == 0){
     scalar_IC += 0.5*(1.0-std::tanh((dist-radius7)/(1.0)));
 }
 else if (index==8){
    scalar_IC += 0.5*(1.0-std::tanh((dist-radius7)/(1.0)));
}
/////////////////////////////////////////////
dist = 0.0;
for (unsigned int dir = 0; dir < dim; dir++){
    dist += (p[dir]-center8[dir])*(p[dir]-center8[dir]);
}
dist = std::sqrt(dist);

// Initial condition for the concentration field
if (index == 0){
    scalar_IC += 0.5*(1.0-std::tanh((dist-radius8)/(1.0)));
}
else if (index==9){
   scalar_IC += 0.5*(1.0-std::tanh((dist-radius8)/(1.0)));
}
//////////////////////////////////////////////////////
dist = 0.0;
for (unsigned int dir = 0; dir < dim; dir++){
    dist += (p[dir]-center9[dir])*(p[dir]-center9[dir]);
}
dist = std::sqrt(dist);

// Initial condition for the concentration field
if (index == 0){
    scalar_IC += 0.5*(1.0-std::tanh((dist-radius9)/(1.0)));
}
else if (index==10){
   scalar_IC += 0.5*(1.0-std::tanh((dist-radius9)/(1.0)));
}

/////////////////////////////////////////////////////////////////

dist = 0.0;
for (unsigned int dir = 0; dir < dim; dir++){
    dist += (p[dir]-center10[dir])*(p[dir]-center10[dir]);
}
dist = std::sqrt(dist);

// Initial condition for the concentration field
if (index == 0){
    scalar_IC += 0.5*(1.0-std::tanh((dist-radius10)/(1.0)));
}
else if (index==11){
   scalar_IC += 0.5*(1.0-std::tanh((dist-radius10)/(1.0)));
}

/////////////////////////////////////////////////

dist = 0.0;
for (unsigned int dir = 0; dir < dim; dir++){
    dist += (p[dir]-center11[dir])*(p[dir]-center11[dir]);
}
dist = std::sqrt(dist);

// Initial condition for the concentration field
if (index == 0){
    scalar_IC += 0.5*(1.0-std::tanh((dist-radius11)/(1.0)));
}
else if (index==12){
   scalar_IC += 0.5*(1.0-std::tanh((dist-radius11)/(1.0)));
}
//////////////////////////////////////////////////////////////
dist = 0.0;
for (unsigned int dir = 0; dir < dim; dir++){
    dist += (p[dir]-center12[dir])*(p[dir]-center12[dir]);
}
dist = std::sqrt(dist);

// Initial condition for the concentration field
if (index == 0){
    scalar_IC += 0.5*(1.0-std::tanh((dist-radius12)/(1.0)));
}
else if (index==13){
   scalar_IC += 0.5*(1.0-std::tanh((dist-radius12)/(1.0)));
}
*/

	  // --------------------------------------------------------------------------
}

// ===========================================================================
// FUNCTION FOR NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS
// ===========================================================================

template <int dim, int degree>
void customPDE<dim,degree>::setNonUniformDirichletBCs(const dealii::Point<dim> &p, const unsigned int index, const unsigned int direction, const double time, double & scalar_BC, dealii::Vector<double> & vector_BC)
{
    // --------------------------------------------------------------------------
    // ENTER THE NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS HERE
    // --------------------------------------------------------------------------
    // Enter the function describing conditions for the fields at point "p".
    // Use "if" statements to set the boundary condition for each variable
    // according to its variable index. This function can be left blank if there
    // are no non-uniform Dirichlet boundary conditions. For BCs that change in
    // time, you can access the current time through the variable "time". The
    // boundary index can be accessed via the variable "direction", which starts
    // at zero and uses the same order as the BC specification in parameters.in
    // (i.e. left = 0, right = 1, bottom = 2, top = 3, front = 4, back = 5).


    // -------------------------------------------------------------------------

}
