// \file PdShell.cpp
//
// File containing a version of the inelastic peridynamic shell.
// For details of the model, refer to:
// Behzadinasab, M., Alaydin, M., Trask, N., Bazilevs, Y. (2021). 
// A general-purpose, inelastic, rotation-free Kirchhoff-Love shell formulation for peridynamics. 
// Computer Methods in Applied Mechanics and Engineering.
// ***** If you are using this code, please cite the above reference. *****

// ************************************************************************
//
// BSD 3-Clause License
//
// Copyright (c) 2021, Masoud Behzadinasab 
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// ************************************************************************

// ***** Included Files *****
//
// The folder Geo contains a python file (create_plate_PD.py), which can be
// used to create the particle setting (initial discretization). 
// Alternative means can be used to create Geometry.dat, which includes the
// nodal x- y- z- coordinates and area in the undeformed configuration.
//
// The file user.h includes an interface where the user can modify depending on
// the problem setup. E.g., solver variables, material properties, damage model 
// parameters, body force, and initial and boundary conditions should be 
// specified there.
//

#include <iostream>
#include "math.h"
#include <vector>
#include "string.h"
#include "PdShell.h"
#include "user.h"


#undef  __FUNCT__
#define __FUNCT__ "input"
void input(PARAMETERS *par) 
{
  // This function reads the Geometry file as the input and initializes the PD system
  
  int in,j;
  POINTS *point;

  int dim = par->dim;

  char filenameCoord[256];
	sprintf(filenameCoord,"./Geo/Geometry.dat");

  if (  (fopen(filenameCoord, "rt")) == NULL)  {
    printf("Geometry.dat file NOT found!\n");
  } else {
    FILE *fL = fopen(filenameCoord, "rt");

    char temp[5];

    // read the number of PD nodes
    fscanf(fL, "%d %s %s %s %s", &par->numPoints, &temp, &temp, &temp, &temp);

    printf("numPoints %d \n",par->numPoints);

    // allocate memory for each PD node
    par->puntos = (POINTS*) malloc(sizeof(*par->puntos)*par->numPoints);

    // read the position and volume of each node in the reference/initial setting
    for (in=0; in<(par->numPoints); in++){
      point = &par->puntos[in];
      fscanf(fL, "%d %le %le %le %le",&par->puntos[in].ID, &par->puntos[in].initialCoord[0], &par->puntos[in].initialCoord[1], &par->puntos[in].initialCoord[2], &par->puntos[in].area);
      point->ID -= 1;

      for(j=0; j<dim; j++){
        point->referenceCoord[j] = point->initialCoord[j]; // undeformed configuration is used as the reference for calculating derivatives here
        point->currentCoord[0][j] = point->initialCoord[j];
        point->currentCoord[1][j] = point->initialCoord[j];
      }

      point->horizon = par->horizon; // currently, we give the same horizon to all points; later, we can input this with the Geometry file
    }

    fclose(fL);
  }
}


#undef  __FUNCT__
#define __FUNCT__ "directNeighborSearch"
void directNeighborSearch(PARAMETERS *par) 
{
  // Here we find the neighbors of each point based on the given horizon
  // then we initialize the bond-level quantities 

  int dim = par->dim;

  int n,j,l,in;
  double distance, delta;
  double xyz[3],xyzNgh[3];

  POINTS *point;
  POINTS *neighbor;

  int neighborCounter; // a counter for all neighbors 

  // first we count the number of points within the neighborhood of each point
  for(in=0;in<par->numPoints;++in){
    point = &par->puntos[in];

    for(j=0; j<dim; j++)
      xyz[j] = point->referenceCoord[j];

    neighborCounter=0;

    delta = point->horizon;

    for(n=0;n<par->numPoints;++n){
      neighbor = &par->puntos[n];

      if(in!=n){
        for(j=0; j<dim; j++)
          xyzNgh[j] = neighbor->referenceCoord[j];

        distance = 0.0;
        for(j=0; j<dim; j++)
          distance += (xyz[j]-xyzNgh[j])*(xyz[j]-xyzNgh[j]);
        distance = sqrt(distance);

        if(distance <= delta) 
          neighborCounter++;
      }
    }
    point->numNeighbors = neighborCounter;
  }


  // now we initialize bond-associated fields based on the number of neighbors
  for( in=0;in<par->numPoints;++in){
    point = &par->puntos[in];

    point->neighbors = (int*) calloc(point->numNeighbors, sizeof(int));

    point->influenceState = (double*) calloc(point->numNeighbors, sizeof(double));

    for (j=0;j<2;j++){
      point->parametricCoordinates[j] = (double*) calloc(point->numNeighbors, sizeof(double));
    }

    for (j=0;j<2;j++){
      point->gradientWeight1[j] = (double*) calloc(point->numNeighbors, sizeof(double));
    }

    for (j=0;j<3;j++){
      point->gradientWeight2[j] = (double*) calloc(point->numNeighbors, sizeof(double));
    }

    point->bondDamage[0] = (double*) calloc(point->numNeighbors, sizeof(double));
    point->bondDamage[1] = (double*) calloc(point->numNeighbors, sizeof(double));

    for(l=0;l<par->numLayers;l++){
      point->layerAssociatedBondLevelJacobianDeterminant[0][l] = (double*) calloc(point->numNeighbors, sizeof(double));
      point->layerAssociatedBondLevelJacobianDeterminant[1][l] = (double*) calloc(point->numNeighbors, sizeof(double));
      point->layerAssociatedBondLevelEquivalentPlasticStrain[0][l] = (double*) calloc(point->numNeighbors, sizeof(double));
      point->layerAssociatedBondLevelEquivalentPlasticStrain[1][l] = (double*) calloc(point->numNeighbors, sizeof(double));
      point->layerAssociatedBondLevelVonMisesStress[l] = (double*) calloc(point->numNeighbors, sizeof(double));
      point->layerAssociatedBondLevelMinPrincipalStress[l] = (double*) calloc(point->numNeighbors, sizeof(double));
      point->layerAssociatedBondLevelMaxPrincipalStress[l] = (double*) calloc(point->numNeighbors, sizeof(double));

      for (j=0;j<dim;j++){
        point->layerAssociatedBondDeformedState[l][j] = (double*) calloc(point->numNeighbors, sizeof(double));
      }

      for (j=0;j<dim*dim;j++){
        point->layerAssociatedBondLevelLeftStretchTensor[0][l][j] = (double*) calloc(point->numNeighbors, sizeof(double));
        point->layerAssociatedBondLevelLeftStretchTensor[1][l][j] = (double*) calloc(point->numNeighbors, sizeof(double));
        point->layerAssociatedBondLevelRotationTensor[0][l][j] = (double*) calloc(point->numNeighbors, sizeof(double));
        point->layerAssociatedBondLevelRotationTensor[1][l][j] = (double*) calloc(point->numNeighbors, sizeof(double));
        point->layerAssociatedBondLevelUnrotatedRateOfDeformation[l][j] = (double*) calloc(point->numNeighbors, sizeof(double));
        point->layerAssociatedBondLevelVelocityGradient[l][j] = (double*) calloc(point->numNeighbors, sizeof(double));
        point->layerAssociatedBondLevelUnrotatedCauchyStress[0][l][j] = (double*) calloc(point->numNeighbors, sizeof(double));
        point->layerAssociatedBondLevelUnrotatedCauchyStress[1][l][j] = (double*) calloc(point->numNeighbors, sizeof(double));
        point->layerAssociatedBondLevelCauchyStress[l][j] = (double*) calloc(point->numNeighbors, sizeof(double));
        point->layerAssociatedBondLevelKirchhoffStress[l][j] = (double*) calloc(point->numNeighbors, sizeof(double));
      }
    }

    // add neighbor ID
    for(j=0; j<dim; j++)
      xyz[j] = point->referenceCoord[j];

    neighborCounter=0;

    delta = point->horizon;

    for(n=0;n<par->numPoints && neighborCounter < point->numNeighbors ;++n){
      neighbor = &par->puntos[n];

      if(in!=n){
        for(j=0; j<dim; j++)
          xyzNgh[j] = neighbor->referenceCoord[j];

        distance = 0.0;
        for(j=0; j<dim; j++)
          distance += (xyz[j]-xyzNgh[j])*(xyz[j]-xyzNgh[j]);
        distance = sqrt(distance);

        if(distance <= delta){
          point->neighbors[neighborCounter] = neighbor->ID;
          neighborCounter++;
        }
      }
    }
  }
}


#undef __FUNCT__
#define __FUNCT__ "initializeFields"
void initializeFields(PARAMETERS *par)
{
  int n,j,l,in;
  POINTS *point;
  int dim = par->dim;

  double identity[9]={1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0};

  for (in=0; in<par->numPoints ; in++){
    point = &par->puntos[in];

    for(j=0; j<dim; j++){
      point->displacement[0][j] = 0.0; 
      point->displacement[1][j] = 0.0; 
      point->velocity[0][j] = 0.0;
      point->velocity[1][j] = 0.0;
      point->acceleration[j] = 0.0;
      point->forceDensity[j] = 0.0;
    }

    point->referenceShellThickness = par->shellThickness;
    point->shellThickness[0] = par->shellThickness;
    point->density2D = par->density * point->shellThickness[0];

    point->weightEvaluationFlag = true; // we need to initially calculate the gradient shape functions (weights)

    point->jacobianDeterminant[0] = 1.0;
    point->equivalentPlasticStrain[0] = 0.0;

    for (j=0;j<dim*dim; j++){ 
      point->unrotatedCauchyStress[0][j] = 0.0;
      point->leftStretchTensor[0][j] = identity[j];
      point->deformationGradient[j] = identity[j];
    }

    point->damage[0] = 0.0;
    point->damage[1] = 0.0;

    for(n=0;n<point->numNeighbors;n++){
      point->bondDamage[0][n] = 0.0; 
      point->bondDamage[1][n] = 0.0; 

      for(l=0;l<par->numLayers;l++){

        point->layerAssociatedBondLevelJacobianDeterminant[0][l][n] = 1.0; 
        point->layerAssociatedBondLevelEquivalentPlasticStrain[0][l][n] = 0.0; 

        for (j=0;j<dim*dim;j++){ 
          point->layerAssociatedBondLevelUnrotatedCauchyStress[0][l][j][n] = 0.0; 
          point->layerAssociatedBondLevelLeftStretchTensor[0][l][j][n] = identity[j];
        }
      }
    }
  }
}


#undef  __FUNCT__
#define __FUNCT__ "kernel"
double kernel(int kernelType,
              double distanceRatio)
{
  double returnValue=0.0;

  if(distanceRatio>=1.0){
    return returnValue;
  }

  if(kernelType==0){
    //Heaviside function
    returnValue = 1.0;
  }
  else if(kernelType==1){
    //Hat function
    returnValue = 1.0 - distanceRatio;
  }
  else if(kernelType==2){
    //Quadratic Spline
    if(distanceRatio > 1.0/3.0)
      returnValue = 3.0/2.0 +distanceRatio*(-3.0+1.5*distanceRatio);
    else
      returnValue = 1.0 - 3.0*distanceRatio*distanceRatio;
  }
  else if(kernelType==3){
    //Cubic Spline
    if(distanceRatio > 0.50)
      returnValue = 2 + 2*distanceRatio*(-3 + (3 - distanceRatio)*distanceRatio);
    else
      returnValue = 1 + 6*distanceRatio*distanceRatio*(-1 + distanceRatio);
  }
  else if(kernelType==4){
    //4th order Spline
    if(distanceRatio > 0.6)
      returnValue = 2.7173913043478260870 +distanceRatio*(- 10.869565217391304348 + distanceRatio*(16.304347826086956522 + distanceRatio*(-10.869565217391304348 + 2.7173913043478260870*distanceRatio)));
    else if(distanceRatio > 1.0/5.0 )
      returnValue = 0.95652173913043478261 +distanceRatio*(0.86956521739130434783 + distanceRatio*( - 13.043478260869565217 + distanceRatio*(21.739130434782608696-distanceRatio*10.869565217391304348)));
    else
      returnValue = 1.0 +distanceRatio*distanceRatio * (- 6.5217391304347826087 + 16.304347826086956522*distanceRatio*distanceRatio);
  }
  else if(kernelType==5){
    //Parabolic Decay
    if(distanceRatio > 0.50)
      returnValue = -4.0*distanceRatio*distanceRatio + 4.0*distanceRatio; 
    else
      returnValue = 1.0;
  }
  else{
    printf(" Invalid Kernel Type \n");
  }

  return returnValue;
}

/*******************************************
********** Start of PCA functions **********
*******************************************/

// Symmetric Householder reduction to tridiagonal form.
void tred2(double V[3][3], double d[3], double e[3]) {

//  This is derived from the Algol procedures tred2 by
//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
//  Fortran subroutine in EISPACK.

  for (int j = 0; j < 3; j++) {
    d[j] = V[3-1][j];
  }

  // Householder reduction to tridiagonal form.

  for (int i = 3-1; i > 0; i--) {

    // Scale to avoid under/overflow.

    double scale = 0.0;
    double h = 0.0;
    for (int k = 0; k < i; k++) {
      scale = scale + fabs(d[k]);
    }
    if (scale == 0.0) {
      e[i] = d[i-1];
      for (int j = 0; j < i; j++) {
        d[j] = V[i-1][j];
        V[i][j] = 0.0;
        V[j][i] = 0.0;
      }
    } else {

      // Generate Householder vector.

      for (int k = 0; k < i; k++) {
        d[k] /= scale;
        h += d[k] * d[k];
      }
      double f = d[i-1];
      double g = sqrt(h);
      if (f > 0) {
        g = -g;
      }
      e[i] = scale * g;
      h = h - f * g;
      d[i-1] = f - g;
      for (int j = 0; j < i; j++) {
        e[j] = 0.0;
      }

      // Apply similarity transformation to remaining columns.

      for (int j = 0; j < i; j++) {
        f = d[j];
        V[j][i] = f;
        g = e[j] + V[j][j] * f;
        for (int k = j+1; k <= i-1; k++) {
          g += V[k][j] * d[k];
          e[k] += V[k][j] * f;
        }
        e[j] = g;
      }
      f = 0.0;
      for (int j = 0; j < i; j++) {
        e[j] /= h;
        f += e[j] * d[j];
      }
      double hh = f / (h + h);
      for (int j = 0; j < i; j++) {
        e[j] -= hh * d[j];
      }
      for (int j = 0; j < i; j++) {
        f = d[j];
        g = e[j];
        for (int k = j; k <= i-1; k++) {
          V[k][j] -= (f * e[k] + g * d[k]);
        }
        d[j] = V[i-1][j];
        V[i][j] = 0.0;
      }
    }
    d[i] = h;
  }

  // Accumulate transformations.

  for (int i = 0; i < 3-1; i++) {
    V[3-1][i] = V[i][i];
    V[i][i] = 1.0;
    double h = d[i+1];
    if (h != 0.0) {
      for (int k = 0; k <= i; k++) {
        d[k] = V[k][i+1] / h;
      }
      for (int j = 0; j <= i; j++) {
        double g = 0.0;
        for (int k = 0; k <= i; k++) {
          g += V[k][i+1] * V[k][j];
        }
        for (int k = 0; k <= i; k++) {
          V[k][j] -= g * d[k];
        }
      }
    }
    for (int k = 0; k <= i; k++) {
      V[k][i+1] = 0.0;
    }
  }
  for (int j = 0; j < 3; j++) {
    d[j] = V[3-1][j];
    V[3-1][j] = 0.0;
  }
  V[3-1][3-1] = 1.0;
  e[0] = 0.0;
} 

// Symmetric tridiagonal QL algorithm.
void tql2(double V[3][3], double d[3], double e[3]) {

//  This is derived from the Algol procedures tql2, by
//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
//  Fortran subroutine in EISPACK.

  for (int i = 1; i < 3; i++) {
    e[i-1] = e[i];
  }
  e[3-1] = 0.0;

  double f = 0.0;
  double tst1 = 0.0;
  double eps = pow(2.0,-52.0);
  for (int l = 0; l < 3; l++) {

    // Find small subdiagonal element

    tst1 = fmax(tst1,fabs(d[l]) + fabs(e[l]));
    int m = l;
    while (m < 3) {
      if (fabs(e[m]) <= eps*tst1) {
        break;
      }
      m++;
    }

    // If m == l, d[l] is an eigenvalue,
    // otherwise, iterate.

    if (m > l) {
      int iter = 0;
      do {
        iter = iter + 1;  // (Could check iteration count here.)

        // Compute implicit shift

        double g = d[l];
        double p = (d[l+1] - g) / (2.0 * e[l]);
        double r = sqrt(p*p+1.0);
        if (p < 0) {
          r = -r;
        }
        d[l] = e[l] / (p + r);
        d[l+1] = e[l] * (p + r);
        double dl1 = d[l+1];
        double h = g - d[l];
        for (int i = l+2; i < 3; i++) {
          d[i] -= h;
        }
        f = f + h;

        // Implicit QL transformation.

        p = d[m];
        double c = 1.0;
        double c2 = c;
        double c3 = c;
        double el1 = e[l+1];
        double s = 0.0;
        double s2 = 0.0;
        for (int i = m-1; i >= l; i--) {
          c3 = c2;
          c2 = c;
          s2 = s;
          g = c * e[i];
          h = c * p;
          r = sqrt(p*p+e[i]*e[i]);
          e[i+1] = s * r;
          s = e[i] / r;
          c = p / r;
          p = c * d[i] - s * g;
          d[i+1] = h + s * (c * g + s * d[i]);

          // Accumulate transformation.

          for (int k = 0; k < 3; k++) {
            h = V[k][i+1];
            V[k][i+1] = s * V[k][i] + c * h;
            V[k][i] = c * V[k][i] - s * h;
          }
        }
        p = -s * s2 * c3 * el1 * e[l] / dl1;
        e[l] = s * p;
        d[l] = c * p;

        // Check for convergence.

      } while (fabs(e[l]) > eps*tst1);
    }
    d[l] = d[l] + f;
    e[l] = 0.0;
  }
  
  // Sort eigenvalues and corresponding vectors.
  for (int i = 0; i < 3-1; i++) {
    int k = i;
    double p = d[i];
    for (int j = i+1; j < 3; j++) {
      if (d[j] < p) {
        k = j;
        p = d[j];
      }
    }
    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      for (int j = 0; j < 3; j++) {
        p = V[j][i];
        V[j][i] = V[j][k];
        V[j][k] = p;
      }
    }
  }
}


#undef  __FUNCT__
#define __FUNCT__ "calculateTwoLargestEigenvectors"
void calculateTwoLargestEigenvectors
(
    const double* input3x3SymmetricMatrix,
    double* psi1,
    double* psi2
)
{
  // This function calculates the two largest eigenvectors of a 3x3 symmetric matrix

  double V[3][3], d[3], e[3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      V[i][j] = *(input3x3SymmetricMatrix+3*i+j);
    }
  }
  tred2(V, d, e);
  tql2(V, d, e);

  // grab only the two largest eigenvectors (sorted in ascending order)
  for (int i = 0; i < 3; i++) {
    *(psi1+i) = V[i][1];
    *(psi2+i) = V[i][2];
  }
}


#undef  __FUNCT__
#define __FUNCT__ "constructTwoDimensionalSpace"
void constructTwoDimensionalSpace(PARAMETERS *par) 
{
  // This function constructs 2D coordinate charts locally on each neighborhood
  // using the positions of the nodes in the point cloud
  // Principal Component Analysis (PCA) is employed to approximate
  // local tangent planes to the point cloud. 

  int dim = par->dim;

  POINTS *point;
  POINTS *neighbor;

  double *modelCoord;
  double *neighborModelCoord;
  double *psi1;
  double *psi2;
  double *eta;
  double *expNormalDir;

  double **xi;

  double undeformedBond[dim];

  double averageModelCoord[dim];

  double covariantTensor[dim*dim]={0.0};

  double temp[dim];

  double etaDotExpectedNormal;

  int in,n,j,k;

  int numNeighbors, neigh;

  for(in=0;in<par->numPoints;in++){

    point = &par->puntos[in];

    modelCoord = point->initialCoord;
    psi1 = point->parametricVector1;
    psi2 = point->parametricVector2;
    eta = point->etaVector;
    expNormalDir = point->expectedNormalVector;

    // Compute the centroid point to the set
    // Starting with the contribution of the centering node first
    for(j=0;j<dim;j++){
      *(averageModelCoord+j) = *(modelCoord+j);
    }

    numNeighbors = point->numNeighbors;
    for(n=0; n<numNeighbors; n++){

      neigh = point->neighbors[n];
      neighbor = &par->puntos[neigh];

      neighborModelCoord = neighbor->initialCoord;

      for(j=0;j<dim;j++){
        *(averageModelCoord+j) += *(neighborModelCoord+j);
      }
    }

    for(j=0;j<dim;j++){
      *(averageModelCoord+j) /= (numNeighbors+1.0);
    }


    // Compute the covariant tensor
    // Starting with the contribution of the centering node first
    for(j=0;j<dim;j++){
      *(undeformedBond+j) = *(modelCoord+j) - *(averageModelCoord+j);
    }

    for(j=0;j<dim;j++){
      for(k=0;k<dim;k++){
        *(covariantTensor+dim*j+k) = *(undeformedBond+j) * *(undeformedBond+k);
      }
    }

    //Re-iterate over the neighbor set and add each neighbor's contribution to
    //the covariant tensor
    numNeighbors = point->numNeighbors;
    for(n=0; n<numNeighbors; n++){

      neigh = point->neighbors[n];
      neighbor = &par->puntos[neigh];

      neighborModelCoord = neighbor->initialCoord;

      for(j=0;j<dim;j++){
        *(undeformedBond+j) = *(neighborModelCoord+j) - *(averageModelCoord+j);
      }

      for(j=0;j<dim;j++){
        for(k=0;k<dim;k++){
          *(covariantTensor+dim*j+k) += *(undeformedBond+j) * *(undeformedBond+k);
        }
      }
    }

    for(j=0;j<dim;j++){
      for(k=0;k<dim;k++){
        *(covariantTensor+dim*j+k) /= (numNeighbors+1.0);
      }
    }

    //Compute the two largest eigenvectors of the covariant tensor
    //Normalize and call them psi1 and psi2
    calculateTwoLargestEigenvectors(covariantTensor,
                                    psi1,
                                    psi2);
    
    // eta = psi1 x psi2
    // for 3D
    *(eta+0) = *(psi1+1) * *(psi2+2) - *(psi1+2) * *(psi2+1);
    *(eta+1) = *(psi1+2) * *(psi2+0) - *(psi1+0) * *(psi2+2);
    *(eta+2) = *(psi1+0) * *(psi2+1) - *(psi1+1) * *(psi2+0);

    // check if the calculated normal direction is in the same direction with
    // the expected normal (based on the initial geometry)
    etaDotExpectedNormal = 0.0;
    for(j=0;j<dim;j++)
      etaDotExpectedNormal += *(eta+j) * *(expNormalDir+j);

    if(etaDotExpectedNormal < 0.0){
      // need to reverse psi1 and psi2
      for(j=0;j<dim;j++){
        *(temp+j) = *(psi1+j);
        *(psi1+j) = *(psi2+j);
        *(psi2+j) = *(temp+j);
        *(eta+j) *= -1.0;
      }
    }

    // Iterate over the neighbor set and compute xi1, xi2
   
    numNeighbors = point->numNeighbors;
    xi = point->parametricCoordinates;
    for(n=0; n<numNeighbors; n++){

      neigh = point->neighbors[n];
      neighbor = &par->puntos[neigh];

      neighborModelCoord = neighbor->initialCoord;

      for(j=0;j<dim;j++){
        *(undeformedBond+j) = *(neighborModelCoord+j) - *(modelCoord+j);
      }

      // Zero out xi1, xi2
      xi[0][n] = 0.0; xi[1][n] = 0.0;
      for(j=0;j<dim;j++){
        xi[0][n] += *(undeformedBond+j) * *(psi1+j);
        xi[1][n] += *(undeformedBond+j) * *(psi2+j);
      }
    }
  }
}

/*****************************************
********** End of PCA functions **********
*****************************************/

#undef __FUNCT__
#define __FUNCT__ "computeInfluenceFunction"
void computeInfluenceFunction(PARAMETERS *par) 
{
  //This function computes the influence state using the reference positions
  //and damage state

  int dim = par->dim;
  int kernelType = par->kernelType;

  POINTS *point;
  POINTS *neighbor;
  int in,n,j;
  int numNeighbors, neigh;

  double *refCoord;
  double *neighborRefCoord;

  double delta;
  double distances[3];
  double distanceSquared;
  double ratio;

  double *omega;
  double *damage;

  for(in=0;in<par->numPoints;in++){

    point = &par->puntos[in];

    // influence state calculation is only needed initially unless damage grows for a point
    if(point->weightEvaluationFlag){

      refCoord = point->referenceCoord;
      delta = point->horizon;

      omega = point->influenceState;
      damage = point->bondDamage[1];

      numNeighbors = point->numNeighbors;
      for(n=0; n<numNeighbors; n++){
        neigh = point->neighbors[n];
        neighbor = &par->puntos[neigh];

        neighborRefCoord = neighbor->referenceCoord;

        // relative distance
        distanceSquared = 0.0;
        for(j=0; j<dim; j++){
          distances[j] = *(neighborRefCoord+j) - *(refCoord+j);
          distanceSquared += distances[j]*distances[j];
        }
        ratio = sqrt(distanceSquared)/delta;
        omega[n] = (1.0-damage[n]) * kernel(kernelType,ratio);
      }
    }
  }
}


#undef __FUNCT__
#define __FUNCT__ "computeWeightedArea"
void computeWeightedArea(PARAMETERS *par) 
{
  //This function calculates the weighted area for each PD point

  computeInfluenceFunction(par); // Compute influence state 

  POINTS *point;
  int in,n;
  int numNeighbors, neigh;

  double *weightedArea;
  double *omega;
  double neighborArea;

  for(in=0;in<par->numPoints;in++){

    point = &par->puntos[in];
    weightedArea = &(point->weightedArea);
    omega = point->influenceState;

    *weightedArea = 0.0; //initialize

    numNeighbors = point->numNeighbors;
    for(n=0; n<numNeighbors; n++, omega++){
      neigh = point->neighbors[n];
      neighborArea = par->puntos[neigh].area;
      *weightedArea += *omega * neighborArea; 
    }
  }
}


#undef __FUNCT__
#define __FUNCT__ ""
void setRotationTensorToMaterialCoordinates(PARAMETERS *par) 
{
  //This function sets the initial rotation tensor aligned with the material
  //coordinates at each point (i.e., parametric vectors)

  POINTS *point;
  int i,in;

  double *rotTensor;
  double *psi1;
  double *psi2;
  double *eta;

  for(in=0;in<par->numPoints;in++){

    point = &par->puntos[in];

    rotTensor = point->rotationTensor[0];
    psi1 = point->parametricVector1;
    psi2 = point->parametricVector2;
    eta = point->etaVector;

    for(i=0; i<3; i++){
      *(rotTensor+3*i+0) = *(psi1+i);
      *(rotTensor+3*i+1) = *(psi2+i);
      *(rotTensor+3*i+2) = *(eta+i);
    }

  }
}


#undef __FUNCT__
#define __FUNCT__ "setBondLevelRotationTensorToMaterialCoordinates"
void setBondLevelRotationTensorToMaterialCoordinates(PARAMETERS *par) 
{
  //This function sets the initial rotation tensor aligned with the material
  //coordinates at each bond (i.e., parametric vectors)

  POINTS *point;
  POINTS *neighbor;
  int in,n,i,l;
  int numNeighbors, neigh;

  double *rotTensorXX, *rotTensorXY, *rotTensorXZ;
  double *rotTensorYX, *rotTensorYY, *rotTensorYZ;
  double *rotTensorZX, *rotTensorZY, *rotTensorZZ;
  double *eta;
  double *neighborEta;

  double nn[3], eb1[3], eb2[3];
  double nnNorm;

  for(in=0;in<par->numPoints;in++){

    point = &par->puntos[in];

    eta = point->etaVector;

    numNeighbors = point->numNeighbors;
    for(n=0; n<numNeighbors; n++){
      neigh = point->neighbors[n];
      neighbor = &par->puntos[neigh];

      neighborEta = neighbor->etaVector;

      // average the normal vectors
      for(i=0; i<3; i++)
        nn[i] = ( *(eta+i) + *(neighborEta+i) )/2.0;

      nnNorm = sqrt( nn[0] * nn[0] +
                     nn[1] * nn[1] +
                     nn[2] * nn[2]);

      // normalize 
      for(i=0; i<3; i++)
        nn[i] /= nnNorm;

      // find two basis vectors in the plane normal to eta
      // take z-component = 0 for the first vector
      eb1[2] = 0.0;
      if(fabs(nn[0]) < 1.0e-14){
        eb1[1] = 0.0;
        eb1[0] = 1.0;
      }
      else{
        eb1[1] = sqrt(1.0 / (1.0 + (nn[1]*nn[1])/(nn[0]*nn[0])));
        eb1[0] = - nn[1]/nn[0] * eb1[1];
      }

      // eb2 = nn * eb1
      eb2[0] = nn[1]*eb1[2] - nn[2]*eb1[1];
      eb2[1] = nn[2]*eb1[0] - nn[0]*eb1[2];
      eb2[2] = nn[0]*eb1[1] - nn[1]*eb1[0];

      for(l=0;l<par->numLayers;l++){
        rotTensorXX = &point->layerAssociatedBondLevelRotationTensor[0][l][0][n];
        rotTensorXY = &point->layerAssociatedBondLevelRotationTensor[0][l][1][n];
        rotTensorXZ = &point->layerAssociatedBondLevelRotationTensor[0][l][2][n];
        rotTensorYX = &point->layerAssociatedBondLevelRotationTensor[0][l][3][n];
        rotTensorYY = &point->layerAssociatedBondLevelRotationTensor[0][l][4][n];
        rotTensorYZ = &point->layerAssociatedBondLevelRotationTensor[0][l][5][n];
        rotTensorZX = &point->layerAssociatedBondLevelRotationTensor[0][l][6][n];
        rotTensorZY = &point->layerAssociatedBondLevelRotationTensor[0][l][7][n];
        rotTensorZZ = &point->layerAssociatedBondLevelRotationTensor[0][l][8][n];

        *rotTensorXX = eb1[0]; *rotTensorYX = eb1[1]; *rotTensorZX = eb1[2];
        *rotTensorXY = eb2[0]; *rotTensorYY = eb2[1]; *rotTensorZY = eb2[2];
        *rotTensorXZ = nn[0]; *rotTensorYZ = nn[1]; *rotTensorZZ = nn[2];
      }
    }
  }
}


#undef __FUNCT__
#define __FUNCT__ "updateShellThickness"
void updateShellThickness(PARAMETERS *par) 
{
  //This function updates the shell thickness at each PD point

  POINTS *point;
  int in;

  double *thicknessN, *thicknessNP1;
  double *D33;

  double dt = par->timeStep;

  for(in=0;in<par->numPoints;in++){

    point = &par->puntos[in];
    thicknessN = &(point->shellThickness[0]);
    thicknessNP1 = &(point->shellThickness[1]);
    D33 = &(point->unrotatedRateOfDeformation[8]);

    // initialize
    *thicknessNP1 = *thicknessN;

    // hdot = h * D_33
    *thicknessNP1 += *thicknessN * *D33 * dt;
  }
}


#undef __FUNCT__
#define __FUNCT__ "updateGradientWeightEvaluationFlag"
void updateGradientWeightEvaluationFlag(PARAMETERS *par) 
{
  //In this function we indicate if damage grows at each point
  //The determined flag will be used to see if a recalculation of weights is needed

  POINTS *point;
  int in;

  double *damageN, *damageNP1;
  bool *weightEvaluationFlag;

  const double tol = 1.0e-12;

  for(in=0;in<par->numPoints;in++){

    point = &par->puntos[in];
    damageN = &(point->damage[0]);
    damageNP1 = &(point->damage[1]);
    weightEvaluationFlag = &(point->weightEvaluationFlag);

    if(*damageNP1 - *damageN > tol){
      // Damage grows. We need to update gradient weights
      *weightEvaluationFlag = true;
    }
  }
}


#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static double PYTHAG(double a, 
                     double b)
{
  double at = fabs(a), bt = fabs(b), ct, result;

  if (at > bt){ 
    ct = bt / at; 
    result = at * sqrt(1.0 + ct * ct); 
  } else if(bt > 0.0){
    ct = at / bt; 
    result = bt * sqrt(1.0 + ct * ct); 
  } else 
    result = 0.0;

  return(result);
}


int dsvd(double *a,
         int m,
         int n,
         double *w,
         double *v)
{
  int flag, i, its, j, jj, k, l, nm;
  double c, f, h, s, x, y, z;
  double anorm = 0.0, g = 0.0, scale = 0.0;
  double *rv1;

  if (m < n){
    fprintf(stderr, "#rows must be > #cols \n");
    return(0);
  }

  rv1 = (double *)malloc((unsigned int) n*sizeof(double));

  // Householder reduction to bidiagonal form 
  for (i = 0; i < n; i++){
    // left-hand reduction 
    l = i + 1;
    rv1[i] = scale * g;
    g = s = scale = 0.0;

    if (i < m){
      for (k = i; k < m; k++)
        scale += fabs((double)a[k*n+i]);

      if (scale){
        for (k = i; k < m; k++){
          a[k*n+i] = (double)((double)a[k*n+i]/scale);
          s += ((double)a[k*n+i] * (double)a[k*n+i]);
        }

        f = (double)a[i*n+i];
        g = -SIGN(sqrt(s), f);
        h = f * g - s;
        a[i*n+i] = (double)(f - g);

        if (i != n - 1){
          for (j = l; j < n; j++){
            for (s = 0.0, k = i; k < m; k++)
              s += ((double)a[k*n+i] * (double)a[k*n+j]);

            f = s / h;

            for (k = i; k < m; k++)
              a[k*n+j] += (double)(f * (double)a[k*n+i]);
          }
        }

        for (k = i; k < m; k++)
          a[k*n+i] = (double)((double)a[k*n+i]*scale);
      }
    }
    w[i] = (double)(scale * g);

    // right-hand reduction 
    g = s = scale = 0.0;

    if (i < m && i != n - 1){
      for (k = l; k < n; k++)
        scale += fabs((double)a[i*n+k]);

      if (scale){
        for (k = l; k < n; k++){
          a[i*n+k] = (double)((double)a[i*n+k]/scale);
          s += ((double)a[i*n+k] * (double)a[i*n+k]);
        }

        f = (double)a[i*n+l];
        g = -SIGN(sqrt(s), f);
        h = f * g - s;
        a[i*n+l] = (double)(f - g);

        for (k = l; k < n; k++)
          rv1[k] = (double)a[i*n+k] / h;

        if (i != m - 1){
          for (j = l; j < m; j++){
            for (s = 0.0, k = l; k < n; k++)
              s += ((double)a[j*n+k] * (double)a[i*n+k]);

            for (k = l; k < n; k++)
              a[j*n+k] += (double)(s * rv1[k]);
          }
        }

        for (k = l; k < n; k++)
        a[i*n+k] = (double)((double)a[i*n+k]*scale);
      }
    }

    anorm = fmax(anorm, (fabs((double)w[i]) + fabs(rv1[i])));
  }

  // accumulate the right-hand transformation 
  for (i = n - 1; i >= 0; i--){
    if (i < n - 1){
      if (g){
        for (j = l; j < n; j++)
          v[j*n+i] = (double)(((double)a[i*n+j] / (double)a[i*n+l]) / g);

        // double division to avoid underflow 
        for (j = l; j < n; j++){
          for (s = 0.0, k = l; k < n; k++)
            s += ((double)a[i*n+k] * (double)v[k*n+j]);

          for (k = l; k < n; k++)
            v[k*n+j] += (double)(s * (double)v[k*n+i]);
        }
      }
      for (j = l; j < n; j++)
        v[i*n+j] = v[j*n+i] = 0.0;
    }
    v[i*n+i] = 1.0;
    g = rv1[i];
    l = i;
  }

  // accumulate the left-hand transformation
  for (i = n - 1; i >= 0; i--){
    l = i + 1;
    g = (double)w[i];

    if (i < n - 1)
      for (j = l; j < n; j++)
        a[i*n+j] = 0.0;

    if (g){
      g = 1.0 / g;

      if (i != n - 1){
        for (j = l; j < n; j++){
          for (s = 0.0, k = l; k < m; k++)
            s += ((double)a[k*n+i] * (double)a[k*n+j]);

          f = (s / (double)a[i*n+i]) * g;

          for (k = i; k < m; k++)
            a[k*n+j] += (double)(f * (double)a[k*n+i]);
        }
      }

      for (j = i; j < m; j++)
      a[j*n+i] = (double)((double)a[j*n+i]*g);
    }
    else
    {
      for (j = i; j < m; j++)
        a[j*n+i] = 0.0;
    }
    ++a[i*n+i];
  }

  //diagonalize the bidiagonal form 
  for (k = n - 1; k >= 0; k--){                             //loop over singular values 
    for (its = 0; its < 30; its++){                         //loop over allowed iterations 
      flag = 1;

      for (l = k; l >= 0; l--){                     //test for splitting 
        nm = l - 1;

        if (fabs(rv1[l]) + anorm == anorm){
          flag = 0;
          break;
        }

        if (fabs((double)w[nm]) + anorm == anorm)
          break;
      }

      if (flag){
        s = 1.0;

        for (i = l; i <= k; i++){
          f = s * rv1[i];

          if (fabs(f) + anorm != anorm){
            g = (double)w[i];
            h = PYTHAG(f, g);
            w[i] = (double)h;
            h = 1.0 / h;
            c = g * h;
            s = (- f * h);
            for (j = 0; j < m; j++){
              y = (double)a[j*n+nm];
              z = (double)a[j*n+i];
              a[j*n+nm] = (double)(y * c + z * s);
              a[j*n+i] = (double)(z * c - y * s);
            }
          }
        }
      }

      z = (double)w[k];
      if (l == k){                  // convergence 
        if (z < 0.0){              //make singular value nonnegative 
          w[k] = (double)(-z);
          for (j = 0; j < n; j++)
            v[j*n+k] = (-v[j*n+k]);
        }
        break;
      }

      if (its >= 30) {
        free (rv1);

        fprintf(stderr, "No convergence after 30,000! iterations \n");
        return(0);
      }

      // shift from bottom 2 x 2 minor 
      x = (double)w[l];
      nm = k - 1;
      y = (double)w[nm];
      g = rv1[nm];
      h = rv1[k];
      f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
      g = PYTHAG(f, 1.0);
      f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

      // next QR transformation 
      c = s = 1.0;

      for (j = l; j <= nm; j++){
        i = j + 1;
        g = rv1[i];
        y = (double)w[i];
        h = s * g;
        g = c * g;
        z = PYTHAG(f, h);
        rv1[j] = z;
        c = f / z;
        s = h / z;
        f = x * c + g * s;
        g = g * c - x * s;
        h = y * s;
        y = y * c;
        for (jj = 0; jj < n; jj++){
          x = (double)v[jj*n+j];
          z = (double)v[jj*n+i];
          v[jj*n+j] = (double)(x * c + z * s);
          v[jj*n+i] = (double)(z * c - x * s);
        }

        z = PYTHAG(f, h);
        w[j] = (double)z;

        if (z){
          z = 1.0 / z;
          c = f * z;
          s = h * z;
        }

        f = (c * g) + (s * y);
        x = (c * y) - (s * g);

        for (jj = 0; jj < m; jj++){
          y = (double)a[jj*n+j];
          z = (double)a[jj*n+i];
          a[jj*n+j] = (double)(y * c + z * s);
          a[jj*n+i] = (double)(z * c - y * s);
        }
      }

      rv1[l] = 0.0;
      rv1[k] = f;
      w[k] = (double)x;
    }
  }

  free (rv1);

  return(1);
}


#undef  __FUNCT__
#define __FUNCT__ "invertAndCond"
void invertAndCond(double *Min,
                   double *Mout,
                   int size,
                   double thresVal)
{
  //    double conditioning;
  double *v=(double*)calloc(size*size,sizeof(double));
  double *w=(double*)calloc(size,sizeof(double));
  double *u=(double*)calloc(size*size,sizeof(double));

  int i,j,k;

  for(i=0;i<size*size;i++){
    u[i]=Min[i];
    Mout[i]=0.0;
  }
  dsvd(u, size, size, w, v);

  for(i=0 ; i < size ; i++ ){
    for(j=0 ; j < size ; j++ ){
      for( k=0 ; k < size ; k++ ){
        if(w[k] > thresVal){ // pseudo inverse approach (ignore zero eigs)
          Mout[i*size+j] += v[i*size+k]*1.0/w[k]*u[j*size+k];
        }
      }
    }
  }

  free(u);
  free(w);
  free(v);
}


#undef __FUNCT__
#define __FUNCT__ "computeGradientShapeFunctions"
void computeGradientShapeFunctions(PARAMETERS *par) 
{

  //Using RK Method to calculate the gradient shape functions (same as the PD Gradient Operator)
  //these shape functions are constructed on the parametric space (2D gradients).
  //we need first-order as well as second-order derivatives for KL shell

  updateGradientWeightEvaluationFlag(par); // Update the flag for evaluation of gradient weights based on damage growth

  computeInfluenceFunction(par); // Compute influence state (update as damage grows)
  
  POINTS *point;
  POINTS *neighbor;
  int in;

  int numPoints = par->numPoints;
  int basisOrder = par->orderOfBasis;
  int i,j,n;

  double delta;

  int numNeigh;
  int neighborIndex;
  int *neighbors;
  double *omega;

  double **xi;
  double **phi1;
  double **phi2;

  double neighborArea, temp;

  // Calculate the dimension of Q vector (parametric space is 2d)
  int Qdim = (basisOrder+1)*(basisOrder+2)/2 - 1;

  double Q[Qdim];
  double M[Qdim*Qdim];
  double Minv[Qdim*Qdim];

  // Constants for calculating gradient weights 1, 2, 11, 12, 22
  double phi1const = 1.0; int phi1ind = 0;
  double phi2const = 1.0; int phi2ind = 1;
  double phi11const = 2.0; int phi11ind = 2;
  double phi12const = 1.0; int phi12ind = 3;
  double phi22const = 2.0; int phi22ind = 4;

  int counter, thisOrder, p1, p2;

  double thresVal;

  for(in=0;in<numPoints;in++){
    point = &par->puntos[in];

    // weight calculation is only needed initially unless damage grows for a point
    if(point->weightEvaluationFlag){

      delta = point->horizon;

      numNeigh = point->numNeighbors;
      neighbors = point->neighbors;
      omega = point->influenceState;

      xi = point->parametricCoordinates;
      phi1 = point->gradientWeight1;
      phi2 = point->gradientWeight2;

      // Zero out data
      for(j=0; j<Qdim; j++)
        for(i=0; i<Qdim; i++)
          *(M+Qdim*i+j) = 0.0;

      // Calculate Moment matrix
      for(n=0;n<numNeigh;n++){

        neighborIndex = neighbors[n];
        neighbor = &par->puntos[neighborIndex];
        neighborArea = neighbor->area;

        // Calculate Q for this bond
        counter = 0;
        for(thisOrder=1; thisOrder<=basisOrder; thisOrder++){
          for(p1=thisOrder; p1>=0; p1--){ // xi1-power
            p2=thisOrder-p1; //xi2-power

            Q[counter] = 1.0;
            for(i=0; i<p1; i++)
              Q[counter] *= xi[0][n]/delta;
            for(i=0; i<p2; i++)
              Q[counter] *= xi[1][n]/delta;

            counter++;
          }
        }

        temp = omega[n] * neighborArea;
        for(j=0; j<Qdim; j++)
          for(i=0; i<Qdim; i++)
            *(M+Qdim*i+j) += temp * Q[i] * Q[j];
      }

      thresVal = 1.0e-5 * delta * delta; // this value will be used to detect a bad conditioned moment matrix 

      invertAndCond(M, Minv, Qdim, thresVal);

      // Re-iterate over the neighbor set and compute Phis
      for(n=0;n<numNeigh;n++){

        // Calculate Q for this bond
        counter = 0;
        for(thisOrder=1; thisOrder<=basisOrder; thisOrder++){
          for(p1=thisOrder; p1>=0; p1--){ // xi1-power
            p2=thisOrder-p1; //xi2-power

            Q[counter] = 1.0;
            for(i=0; i<p1; i++)
              Q[counter] *= xi[0][n]/delta;
            for(i=0; i<p2; i++)
              Q[counter] *= xi[1][n]/delta;

            counter++;
          }
        }

        // caluclate phi1
        phi1[0][n] = 0.0;
        temp = phi1const * omega[n];  
        for(j=0; j<Qdim; j++)
          phi1[0][n] += temp * *(Minv+Qdim*phi1ind+j) * Q[j]/delta;

        // caluclate phi2
        phi1[1][n] = 0.0;
        temp = phi2const * omega[n];  
        for(j=0; j<Qdim; j++)
          phi1[1][n] += temp * *(Minv+Qdim*phi2ind+j) * Q[j]/delta;

        // caluclate phi11
        phi2[0][n] = 0.0;
        temp = phi11const * omega[n];  
        for(j=0; j<Qdim; j++)
          phi2[0][n] += temp * *(Minv+Qdim*phi11ind+j) * Q[j]/delta;

        // caluclate phi12
        phi2[1][n] = 0.0;
        temp = phi12const * omega[n];  
        for(j=0; j<Qdim; j++)
          phi2[1][n] += temp * *(Minv+Qdim*phi12ind+j) * Q[j]/delta;

        // caluclate phi22
        phi2[2][n] = 0.0;
        temp = phi22const * omega[n];  
        for(j=0; j<Qdim; j++)
          phi2[2][n] += temp * *(Minv+Qdim*phi22ind+j) * Q[j]/delta;
      }

      // no need to recalculate the weights next time
      point->weightEvaluationFlag = false;
    }
  }
}


#undef __FUNCT__
#define __FUNCT__ "computeReferenceDeformationGradient"
void computeReferenceDeformationGradient(PARAMETERS *par) 
{
  // Calculate the reference deformation gradient 
  // mapping from parametric space to the undeformed physical space

  POINTS *point;
  POINTS *neighbor;

  int dim = par->dim;
  int i,j,n,in;

  int numNeigh;
  int neighborIndex;
  int *neighbors;

  double* coord;
  double* neighborCoord;
  double neighborArea;

  double* phi1;
  double* phi2;
  double* phi11;
  double* phi12;
  double* phi22;

  double* defGrad1;
  double* defGrad2;
  double* defGrad11;
  double* defGrad12;
  double* defGrad22;

  double defState[3];

  for(in=0;in<par->numPoints;in++){
    point = &par->puntos[in];

    coord = point->referenceCoord;

    defGrad1 = point->referenceParametricDeformationGradient1[0];
    defGrad2 = point->referenceParametricDeformationGradient1[1];
    defGrad11 = point->referenceParametricDeformationGradient2[0];
    defGrad12 = point->referenceParametricDeformationGradient2[1];
    defGrad22 = point->referenceParametricDeformationGradient2[2];

    // Zero out data
    for(j=0; j<dim; j++){
      *(defGrad1+j) = 0.0;
      *(defGrad2+j) = 0.0;
      *(defGrad11+j) = 0.0;
      *(defGrad12+j) = 0.0;
      *(defGrad22+j) = 0.0;
    }

    numNeigh = point->numNeighbors;
    neighbors = point->neighbors;

    phi1 = point->gradientWeight1[0];
    phi2 = point->gradientWeight1[1];
    phi11 = point->gradientWeight2[0];
    phi12 = point->gradientWeight2[1];
    phi22 = point->gradientWeight2[2];

    // Calculate the deformation gradient
    for(n=0;n<numNeigh;n++,phi1++, phi2++, phi11++, phi12++, phi22++){

      neighborIndex = neighbors[n];
      neighbor = &par->puntos[neighborIndex];
      neighborArea = neighbor->area;
      neighborCoord = neighbor->referenceCoord;

      for(i=0; i<dim; i++){
        *(defState+i) = *(neighborCoord+i) - *(coord+i);
      }

      for(i=0; i<dim; i++){
        *(defGrad1+i) += *(defState+i) * *phi1 * neighborArea;
        *(defGrad2+i) += *(defState+i) * *phi2 * neighborArea;
        *(defGrad11+i) += *(defState+i) * *phi11 * neighborArea;
        *(defGrad12+i) += *(defState+i) * *phi12 * neighborArea;
        *(defGrad22+i) += *(defState+i) * *phi22 * neighborArea;
      }
    }
  }
}



#undef __FUNCT__
#define __FUNCT__ "computeParametricGradients"
void computeParametricGradients(PARAMETERS *par) 
{
  // Calculate the parametric gradients

  POINTS *point;
  POINTS *neighbor;

  int dim = par->dim;
  int i,j,n,in;

  int numNeigh;
  int neighborIndex;
  int *neighbors;

  double* disp;
  double* neighborDisp;
  double* vel;
  double* neighborVel;
  double neighborArea;

  double* phi1;
  double* phi2;
  double* phi11;
  double* phi12;
  double* phi22;

  double* refDefGrad1;
  double* refDefGrad2;
  double* refDefGrad11;
  double* refDefGrad12;
  double* refDefGrad22;

  double* defGrad1;
  double* defGrad2;
  double* defGrad11;
  double* defGrad12;
  double* defGrad22;

  double* velGrad1;
  double* velGrad2;
  double* velGrad11;
  double* velGrad12;
  double* velGrad22;

  double dispState[3];
  double velState[3];

  for(in=0;in<par->numPoints;in++){
    point = &par->puntos[in];

    disp = point->displacement[1];
    vel = point->velocity[1];

    refDefGrad1 = point->referenceParametricDeformationGradient1[0];
    refDefGrad2 = point->referenceParametricDeformationGradient1[1];
    refDefGrad11 = point->referenceParametricDeformationGradient2[0];
    refDefGrad12 = point->referenceParametricDeformationGradient2[1];
    refDefGrad22 = point->referenceParametricDeformationGradient2[2];

    defGrad1 = point->parametricDeformationGradient1[0];
    defGrad2 = point->parametricDeformationGradient1[1];
    defGrad11 = point->parametricDeformationGradient2[0];
    defGrad12 = point->parametricDeformationGradient2[1];
    defGrad22 = point->parametricDeformationGradient2[2];

    velGrad1 = point->parametricVelocityGradient1[0];
    velGrad2 = point->parametricVelocityGradient1[1];
    velGrad11 = point->parametricVelocityGradient2[0];
    velGrad12 = point->parametricVelocityGradient2[1];
    velGrad22 = point->parametricVelocityGradient2[2];

    // initialize data
    for(j=0; j<dim; j++){
      *(defGrad1+j) = *(refDefGrad1+j);
      *(defGrad2+j) = *(refDefGrad2+j);
      *(defGrad11+j) = *(refDefGrad11+j);
      *(defGrad12+j) = *(refDefGrad12+j);
      *(defGrad22+j) = *(refDefGrad22+j);

      *(velGrad1+j) = 0.0;
      *(velGrad2+j) = 0.0;
      *(velGrad11+j) = 0.0;
      *(velGrad12+j) = 0.0;
      *(velGrad22+j) = 0.0;
    }

    numNeigh = point->numNeighbors;
    neighbors = point->neighbors;

    phi1 = point->gradientWeight1[0];
    phi2 = point->gradientWeight1[1];
    phi11 = point->gradientWeight2[0];
    phi12 = point->gradientWeight2[1];
    phi22 = point->gradientWeight2[2];

    // Calculate the deformation gradient
    for(n=0;n<numNeigh;n++,phi1++, phi2++, phi11++, phi12++, phi22++){

      neighborIndex = neighbors[n];
      neighbor = &par->puntos[neighborIndex];
      neighborArea = neighbor->area;
      neighborDisp = neighbor->displacement[1];
      neighborVel = neighbor->velocity[1];

      for(i=0; i<dim; i++){
        *(dispState+i) = *(neighborDisp+i) - *(disp+i);
        *(velState+i) = *(neighborVel+i) - *(vel+i);
      }

      for(i=0; i<dim; i++){
        *(defGrad1+i) += *(dispState+i) * *phi1 * neighborArea;
        *(defGrad2+i) += *(dispState+i) * *phi2 * neighborArea;
        *(defGrad11+i) += *(dispState+i) * *phi11 * neighborArea;
        *(defGrad12+i) += *(dispState+i) * *phi12 * neighborArea;
        *(defGrad22+i) += *(dispState+i) * *phi22 * neighborArea;

        *(velGrad1+i) += *(velState+i) * *phi1 * neighborArea;
        *(velGrad2+i) += *(velState+i) * *phi2 * neighborArea;
        *(velGrad11+i) += *(velState+i) * *phi11 * neighborArea;
        *(velGrad12+i) += *(velState+i) * *phi12 * neighborArea;
        *(velGrad22+i) += *(velState+i) * *phi22 * neighborArea;
      }
    }
  }
}


#undef __FUNCT__
#define __FUNCT__ "computeNormalVectorAndRelatedTensors"
void computeNormalVectorAndRelatedTensors(PARAMETERS *par) 
{
  POINTS *point;

  int i,j,in;

  double* defGrad1;
  double* defGrad2;
  double* defGrad11;
  double* defGrad12;
  double* defGrad22;

  double* velGrad1;
  double* velGrad2;

  double* normal;
  double* ndot;
  double* A;
  double* Ad1;
  double* Ad2;
  double* B1;
  double* B1d1;
  double* B1d2;
  double* B2;
  double* B2d1;
  double* B2d2;

  double outerProduct[3], outerProductNorm;

  for(in=0;in<par->numPoints;in++){
    point = &par->puntos[in];

    normal = point->normalVector;
    ndot = point->normalDotVector[1];
    A = point->ATensor;
    Ad1 = point->ATensorGradient[0];
    Ad2 = point->ATensorGradient[1];
    B1 = point->B1Tensor;
    B1d1 = point->B1TensorGradient[0];
    B1d2 = point->B1TensorGradient[1];
    B2 = point->B2Tensor;
    B2d1 = point->B2TensorGradient[0];
    B2d2 = point->B2TensorGradient[1];

    defGrad1 = point->parametricDeformationGradient1[0];
    defGrad2 = point->parametricDeformationGradient1[1];
    defGrad11 = point->parametricDeformationGradient2[0];
    defGrad12 = point->parametricDeformationGradient2[1];
    defGrad22 = point->parametricDeformationGradient2[2];

    velGrad1 = point->parametricVelocityGradient1[0];
    velGrad2 = point->parametricVelocityGradient1[1];

    *(outerProduct+0) = *(defGrad1+1) * *(defGrad2+2) - *(defGrad1+2) * *(defGrad2+1);
    *(outerProduct+1) = *(defGrad1+2) * *(defGrad2+0) - *(defGrad1+0) * *(defGrad2+2);
    *(outerProduct+2) = *(defGrad1+0) * *(defGrad2+1) - *(defGrad1+1) * *(defGrad2+0);
    outerProductNorm = sqrt( *(outerProduct+0) * *(outerProduct+0) +
                             *(outerProduct+1) * *(outerProduct+1) +
                             *(outerProduct+2) * *(outerProduct+2) );

    // calculate normal vector
    for(i=0; i<3; i++)
      *(normal+i) = *(outerProduct+i) / outerProductNorm;
  
    // A tensor 
    *(A+0) = 1.0 - *(normal+0) * *(normal+0); *(A+1) = - *(normal+0) * *(normal+1); *(A+2) = - *(normal+0) * *(normal+2);
    *(A+3) = - *(normal+1) * *(normal+0); *(A+4) = 1.0 - *(normal+1) * *(normal+1); *(A+5) = - *(normal+1) * *(normal+2);
    *(A+6) = - *(normal+2) * *(normal+0); *(A+7) = - *(normal+2) * *(normal+1); *(A+8) = 1.0 - *(normal+2) * *(normal+2);
    for(i=0; i<9; i++)
      *(A+i) /= outerProductNorm;

    // B1 tensor
    for(i=0; i<3; i++){
      *(B1+3*i+0) = *(A+3*i+2) * *(defGrad2+1) - *(A+3*i+1) * *(defGrad2+2);
      *(B1+3*i+1) = *(A+3*i+0) * *(defGrad2+2) - *(A+3*i+2) * *(defGrad2+0);
      *(B1+3*i+2) = *(A+3*i+1) * *(defGrad2+0) - *(A+3*i+0) * *(defGrad2+1);
    }

    // B2 tensor
    for(i=0; i<3; i++){
      *(B2+3*i+0) = *(A+3*i+1) * *(defGrad1+2) - *(A+3*i+2) * *(defGrad1+1);
      *(B2+3*i+1) = *(A+3*i+2) * *(defGrad1+0) - *(A+3*i+0) * *(defGrad1+2);
      *(B2+3*i+2) = *(A+3*i+0) * *(defGrad1+1) - *(A+3*i+1) * *(defGrad1+0);
    }

    // ndot
    for(i=0; i<3; i++){
      *(ndot+i) = 0.0;
      for(j=0; j<3; j++){
        *(ndot+i) += *(B1+3*i+j) * *(velGrad1+j);
        *(ndot+i) += *(B2+3*i+j) * *(velGrad2+j);
      }
    }

    // dAdxi1 tensor
    for(i=0; i<3; i++){
      for(j=0; j<3; j++){
        *(Ad1+3*i+j) = *(normal+j) * *(A+3*i+0) * ( *(defGrad11+1) * *(defGrad2+2) - *(defGrad11+2) * *(defGrad2+1) +
                                               *(defGrad1+1) * *(defGrad12+2) - *(defGrad1+2) * *(defGrad12+1) );
        *(Ad1+3*i+j) += *(normal+j) * *(A+3*i+1) * ( *(defGrad11+2) * *(defGrad2+0) - *(defGrad11+0) * *(defGrad2+2) +
                                                *(defGrad1+2) * *(defGrad12+0) - *(defGrad1+0) * *(defGrad12+2) );
        *(Ad1+3*i+j) += *(normal+j) * *(A+3*i+2) * ( *(defGrad11+0) * *(defGrad2+1) - *(defGrad11+1) * *(defGrad2+0) +
                                                *(defGrad1+0) * *(defGrad12+1) - *(defGrad1+1) * *(defGrad12+0) );

        *(Ad1+3*i+j) += *(normal+i) * *(A+3*j+0) * ( *(defGrad11+1) * *(defGrad2+2) - *(defGrad11+2) * *(defGrad2+1) +
                                                *(defGrad1+1) * *(defGrad12+2) - *(defGrad1+2) * *(defGrad12+1) );
        *(Ad1+3*i+j) += *(normal+i) * *(A+3*j+1) * ( *(defGrad11+2) * *(defGrad2+0) - *(defGrad11+0) * *(defGrad2+2) +
                                                *(defGrad1+2) * *(defGrad12+0) - *(defGrad1+0) * *(defGrad12+2) );
        *(Ad1+3*i+j) += *(normal+i) * *(A+3*j+2) * ( *(defGrad11+0) * *(defGrad2+1) - *(defGrad11+1) * *(defGrad2+0) +
                                                *(defGrad1+0) * *(defGrad12+1) - *(defGrad1+1) * *(defGrad12+0) );

        *(Ad1+3*i+j) += *(A+3*i+j) * *(normal+0) * ( *(defGrad11+1) * *(defGrad2+2) - *(defGrad11+2) * *(defGrad2+1) +
                                                *(defGrad1+1) * *(defGrad12+2) - *(defGrad1+2) * *(defGrad12+1) );
        *(Ad1+3*i+j) += *(A+3*i+j) * *(normal+1) * ( *(defGrad11+2) * *(defGrad2+0) - *(defGrad11+0) * *(defGrad2+2) +
                                                *(defGrad1+2) * *(defGrad12+0) - *(defGrad1+0) * *(defGrad12+2) );
        *(Ad1+3*i+j) += *(A+3*i+j) * *(normal+2) * ( *(defGrad11+0) * *(defGrad2+1) - *(defGrad11+1) * *(defGrad2+0) +
                                                *(defGrad1+0) * *(defGrad12+1) - *(defGrad1+1) * *(defGrad12+0) );
      }
    }
    for(i=0; i<9; i++)
      *(Ad1+i) /= - outerProductNorm;

    // dAdxi2 tensor
    for(i=0; i<3; i++){
      for(j=0; j<3; j++){
        *(Ad2+3*i+j) = *(normal+j) * *(A+3*i+0) * ( *(defGrad12+1) * *(defGrad2+2) - *(defGrad12+2) * *(defGrad2+1) +
                                               *(defGrad1+1) * *(defGrad22+2) - *(defGrad1+2) * *(defGrad22+1) );
        *(Ad2+3*i+j) += *(normal+j) * *(A+3*i+1) * ( *(defGrad12+2) * *(defGrad2+0) - *(defGrad12+0) * *(defGrad2+2) +
                                                *(defGrad1+2) * *(defGrad22+0) - *(defGrad1+0) * *(defGrad22+2) );
        *(Ad2+3*i+j) += *(normal+j) * *(A+3*i+2) * ( *(defGrad12+0) * *(defGrad2+1) - *(defGrad12+1) * *(defGrad2+0) +
                                                *(defGrad1+0) * *(defGrad22+1) - *(defGrad1+1) * *(defGrad22+0) );

        *(Ad2+3*i+j) += *(normal+i) * *(A+3*j+0) * ( *(defGrad12+1) * *(defGrad2+2) - *(defGrad12+2) * *(defGrad2+1) +
                                                *(defGrad1+1) * *(defGrad22+2) - *(defGrad1+2) * *(defGrad22+1) );
        *(Ad2+3*i+j) += *(normal+i) * *(A+3*j+1) * ( *(defGrad12+2) * *(defGrad2+0) - *(defGrad12+0) * *(defGrad2+2) +
                                                *(defGrad1+2) * *(defGrad22+0) - *(defGrad1+0) * *(defGrad22+2) );
        *(Ad2+3*i+j) += *(normal+i) * *(A+3*j+2) * ( *(defGrad12+0) * *(defGrad2+1) - *(defGrad12+1) * *(defGrad2+0) +
                                                *(defGrad1+0) * *(defGrad22+1) - *(defGrad1+1) * *(defGrad22+0) );

        *(Ad2+3*i+j) += *(A+3*i+j) * *(normal+0) * ( *(defGrad12+1) * *(defGrad2+2) - *(defGrad12+2) * *(defGrad2+1) +
                                                *(defGrad1+1) * *(defGrad22+2) - *(defGrad1+2) * *(defGrad22+1) );
        *(Ad2+3*i+j) += *(A+3*i+j) * *(normal+1) * ( *(defGrad12+2) * *(defGrad2+0) - *(defGrad12+0) * *(defGrad2+2) +
                                                *(defGrad1+2) * *(defGrad22+0) - *(defGrad1+0) * *(defGrad22+2) );
        *(Ad2+3*i+j) += *(A+3*i+j) * *(normal+2) * ( *(defGrad12+0) * *(defGrad2+1) - *(defGrad12+1) * *(defGrad2+0) +
                                                *(defGrad1+0) * *(defGrad22+1) - *(defGrad1+1) * *(defGrad22+0) );
      }
    }
    for(i=0; i<9; i++)
      *(Ad2+i) /= - outerProductNorm;

    // dB1dxi1 tensor
    for(i=0; i<3; i++){
      *(B1d1+3*i+0) = *(Ad1+3*i+2) * *(defGrad2+1) - *(Ad1+3*i+1) * *(defGrad2+2) + 
                      *(A+3*i+2) * *(defGrad12+1) - *(A+3*i+1) * *(defGrad12+2); 
      *(B1d1+3*i+1) = *(Ad1+3*i+0) * *(defGrad2+2) - *(Ad1+3*i+2) * *(defGrad2+0) + 
                      *(A+3*i+0) * *(defGrad12+2) - *(A+3*i+2) * *(defGrad12+0); 
      *(B1d1+3*i+2) = *(Ad1+3*i+1) * *(defGrad2+0) - *(Ad1+3*i+0) * *(defGrad2+1) + 
                      *(A+3*i+1) * *(defGrad12+0) - *(A+3*i+0) * *(defGrad12+1); 
    }

    // dB1dxi2 tensor
    for(i=0; i<3; i++){
      *(B1d2+3*i+0) = *(Ad2+3*i+2) * *(defGrad2+1) - *(Ad2+3*i+1) * *(defGrad2+2) + 
                      *(A+3*i+2) * *(defGrad22+1) - *(A+3*i+1) * *(defGrad22+2); 
      *(B1d2+3*i+1) = *(Ad2+3*i+0) * *(defGrad2+2) - *(Ad2+3*i+2) * *(defGrad2+0) + 
                      *(A+3*i+0) * *(defGrad22+2) - *(A+3*i+2) * *(defGrad22+0); 
      *(B1d2+3*i+2) = *(Ad2+3*i+1) * *(defGrad2+0) - *(Ad2+3*i+0) * *(defGrad2+1) + 
                      *(A+3*i+1) * *(defGrad22+0) - *(A+3*i+0) * *(defGrad22+1); 
    }

    // dB2dxi1 tensor
    for(i=0; i<3; i++){
      *(B2d1+3*i+0) = *(Ad1+3*i+1) * *(defGrad1+2) - *(Ad1+3*i+2) * *(defGrad1+1) + 
                      *(A+3*i+1) * *(defGrad11+2) - *(A+3*i+2) * *(defGrad11+1); 
      *(B2d1+3*i+1) = *(Ad1+3*i+2) * *(defGrad1+0) - *(Ad1+3*i+0) * *(defGrad1+2) + 
                      *(A+3*i+2) * *(defGrad11+0) - *(A+3*i+0) * *(defGrad11+2); 
      *(B2d1+3*i+2) = *(Ad1+3*i+0) * *(defGrad1+1) - *(Ad1+3*i+1) * *(defGrad1+0) + 
                      *(A+3*i+0) * *(defGrad11+1) - *(A+3*i+1) * *(defGrad11+0); 
    }

    // dB2dxi2 tensor
    for(i=0; i<3; i++){
      *(B2d2+3*i+0) = *(Ad2+3*i+1) * *(defGrad1+2) - *(Ad2+3*i+2) * *(defGrad1+1) + 
                      *(A+3*i+1) * *(defGrad12+2) - *(A+3*i+2) * *(defGrad12+1); 
      *(B2d2+3*i+1) = *(Ad2+3*i+2) * *(defGrad1+0) - *(Ad2+3*i+0) * *(defGrad1+2) + 
                      *(A+3*i+2) * *(defGrad12+0) - *(A+3*i+0) * *(defGrad12+2); 
      *(B2d2+3*i+2) = *(Ad2+3*i+0) * *(defGrad1+1) - *(Ad2+3*i+1) * *(defGrad1+0) + 
                      *(A+3*i+0) * *(defGrad12+1) - *(A+3*i+1) * *(defGrad12+0); 
    }
  }
}


void Invert3by3Matrix
(
    const double* matrix,
    double& determinant,
    double* inverse
)
{
  double minor0 =  *(matrix+4) * *(matrix+8) - *(matrix+5) * *(matrix+7);
  double minor1 =  *(matrix+3) * *(matrix+8) - *(matrix+5) * *(matrix+6);
  double minor2 =  *(matrix+3) * *(matrix+7) - *(matrix+4) * *(matrix+6);
  double minor3 =  *(matrix+1) * *(matrix+8) - *(matrix+2) * *(matrix+7);
  double minor4 =  *(matrix)   * *(matrix+8) - *(matrix+6) * *(matrix+2);
  double minor5 =  *(matrix)   * *(matrix+7) - *(matrix+1) * *(matrix+6);
  double minor6 =  *(matrix+1) * *(matrix+5) - *(matrix+2) * *(matrix+4);
  double minor7 =  *(matrix)   * *(matrix+5) - *(matrix+2) * *(matrix+3);
  double minor8 =  *(matrix)   * *(matrix+4) - *(matrix+1) * *(matrix+3);
  determinant = *(matrix) * minor0 - *(matrix+1) * minor1 + *(matrix+2) * minor2;

  *(inverse) = minor0/determinant;
  *(inverse+1) = -1.0*minor3/determinant;
  *(inverse+2) = minor6/determinant;
  *(inverse+3) = -1.0*minor1/determinant;
  *(inverse+4) = minor4/determinant;
  *(inverse+5) = -1.0*minor7/determinant;
  *(inverse+6) = minor2/determinant;
  *(inverse+7) = -1.0*minor5/determinant;
  *(inverse+8) = minor8/determinant;
}


#undef __FUNCT__
#define __FUNCT__ "computeDeformationGradientInverse"
void computeDeformationGradientInverse(PARAMETERS *par)
{
  POINTS *point;

  int i,j,l,in;
  double xi3;

  double* thickness;
  double* defGrad1;
  double* defGrad2;
  double* defGrad11;
  double* defGrad12;
  double* defGrad22;
  double* normal;
  double* B1;
  double* B2;

  double* defGradInv;

  double defGrad[9];
  double determinant;

  for(in=0;in<par->numPoints;in++){
    point = &par->puntos[in];

    thickness = &(point->shellThickness[1]);
    defGrad1 = point->parametricDeformationGradient1[0];
    defGrad2 = point->parametricDeformationGradient1[1];
    defGrad11 = point->parametricDeformationGradient2[0];
    defGrad12 = point->parametricDeformationGradient2[1];
    defGrad22 = point->parametricDeformationGradient2[2];
    normal = point->normalVector;
    B1 = point->B1Tensor;
    B2 = point->B2Tensor;

    for(l=0; l<par->numLayers; l++){

      xi3 = par->xi3list[l];
      defGradInv = point->layerAssociatedDeformationGradientInverse[l];

      for(i=0; i<3; i++){

        // F_i1
        *(defGrad+3*i+0) = *(defGrad1+i);
        for(j=0; j<3; j++){
          *(defGrad+3*i+0) += *thickness/2.0 * xi3 * ( *(B1+3*i+j) * *(defGrad11+j) );
          *(defGrad+3*i+0) += *thickness/2.0 * xi3 * ( *(B2+3*i+j) * *(defGrad12+j) );
        }

        // F_i2
        *(defGrad+3*i+1) = *(defGrad2+i);
        for(j=0; j<3; j++){
          *(defGrad+3*i+1) += *thickness/2.0 * xi3 * ( *(B1+3*i+j) * *(defGrad12+j) );
          *(defGrad+3*i+1) += *thickness/2.0 * xi3 * ( *(B2+3*i+j) * *(defGrad22+j) );
        }

        // F_i3
        *(defGrad+3*i+2) = *thickness/2.0 * *(normal+i);
      }

      Invert3by3Matrix(defGrad, determinant, defGradInv);
    }
  }
}


#undef __FUNCT__
#define __FUNCT__ "computeNodeLevelVelocityGradient"
void computeNodeLevelVelocityGradient(PARAMETERS *par)
{
  POINTS *point;

  int i,j,c,l,in;
  double xi3;

  double* thickness;
  double* velGrad1;
  double* velGrad2;
  double* velGrad11;
  double* velGrad12;
  double* velGrad22;
  double* ndot;
  double* B1;
  double* B1d1;
  double* B1d2;
  double* B2;
  double* B2d1;
  double* B2d2;

  double* defGradInv;
  double* velGrad;

  double parVelGrad[9];

  for(in=0;in<par->numPoints;in++){
    point = &par->puntos[in];

    thickness = &(point->shellThickness[1]);
    velGrad1 = point->parametricVelocityGradient1[0];
    velGrad2 = point->parametricVelocityGradient1[1];
    velGrad11 = point->parametricVelocityGradient2[0];
    velGrad12 = point->parametricVelocityGradient2[1];
    velGrad22 = point->parametricVelocityGradient2[2];
    ndot = point->normalDotVector[1];
    B1 = point->B1Tensor;
    B1d1 = point->B1TensorGradient[0];
    B1d2 = point->B1TensorGradient[1];
    B2 = point->B2Tensor;
    B2d1 = point->B2TensorGradient[0];
    B2d2 = point->B2TensorGradient[1];

    for(l=0; l<par->numLayers; l++){

      xi3 = par->xi3list[l];
      defGradInv = point->layerAssociatedDeformationGradientInverse[l];
      velGrad = point->layerAssociatedVelocityGradient[l];

      for(i=0; i<3; i++){
        // L_i1
        *(parVelGrad+3*i+0) = *(velGrad1+i);
        for(j=0; j<3; j++){
          *(parVelGrad+3*i+0) += *thickness/2.0 * xi3 * ( *(B1d1+3*i+j) * *(velGrad1+j) );
          *(parVelGrad+3*i+0) += *thickness/2.0 * xi3 * ( *(B1+3*i+j) * *(velGrad11+j) );
          *(parVelGrad+3*i+0) += *thickness/2.0 * xi3 * ( *(B2d1+3*i+j) * *(velGrad2+j) );
          *(parVelGrad+3*i+0) += *thickness/2.0 * xi3 * ( *(B2+3*i+j) * *(velGrad12+j) );
        }

        // L_i2
        *(parVelGrad+3*i+1) = *(velGrad2+i);
        for(j=0; j<3; j++){
          *(parVelGrad+3*i+1) += *thickness/2.0 * xi3 * ( *(B1d2+3*i+j) * *(velGrad1+j) );
          *(parVelGrad+3*i+1) += *thickness/2.0 * xi3 * ( *(B1+3*i+j) * *(velGrad12+j) );
          *(parVelGrad+3*i+1) += *thickness/2.0 * xi3 * ( *(B2d2+3*i+j) * *(velGrad2+j) );
          *(parVelGrad+3*i+1) += *thickness/2.0 * xi3 * ( *(B2+3*i+j) * *(velGrad22+j) );
        }

        // L_i3
        *(parVelGrad+3*i+2) = *thickness/2.0 * *(ndot+i);
      }

      // compute L_ij = L_ic Finv_cj  ; c = \xi
      for(i=0; i<3; i++){
        for(j=0; j<3; j++){
          *(velGrad+3*i+j) = 0.0;
          for(c=0; c<3; c++){
            *(velGrad+3*i+j) += *(parVelGrad+3*i+c) * *(defGradInv+3*c+j);

          }
        }
      }
    }
  }
}


#undef __FUNCT__
#define __FUNCT__ "computeBondLevelVelocityGradient"
void computeBondLevelVelocityGradient(PARAMETERS *par) 
{
  //Compute the bond-level velocity gradient for each layer

  POINTS *point;
  POINTS *neighbor;
  int i,l,in,n;
  int numNeighbors, neigh;
  double xi3;

  double* thickness;
  double* neighborThickness;
  double* coord;
  double* neighborCoord;
  double* vel;
  double* neighborVel;
  double* normal;
  double* neighborN;
  double* ndot;
  double* neighborNdot;
  double* velGrad;
  double* neighborVelGrad;
  double *bondDefX, *bondDefY, *bondDefZ;
  double *bondLevelVelGradXX, *bondLevelVelGradXY, *bondLevelVelGradXZ;
  double *bondLevelVelGradYX, *bondLevelVelGradYY, *bondLevelVelGradYZ;
  double *bondLevelVelGradZX, *bondLevelVelGradZY, *bondLevelVelGradZZ;

  double* omega;

  double deformedBondX, deformedBondY, deformedBondZ, deformedBondLengthSq;
  double velStateX, velStateY, velStateZ;
  double temp;

  std::vector<double> meanVelGradVector(9);
  double* meanVelGrad = &meanVelGradVector[0];

  for(in=0;in<par->numPoints;in++){

    point = &par->puntos[in];

    thickness = &(point->shellThickness[1]);
    coord = point->currentCoord[1];
    vel = point->velocity[1];
    normal = point->normalVector;
    ndot = point->normalDotVector[1];

    for(l=0; l<par->numLayers; l++){

      xi3 = par->xi3list[l];
      velGrad = point->layerAssociatedVelocityGradient[l];

      bondDefX = point->layerAssociatedBondDeformedState[l][0];
      bondDefY = point->layerAssociatedBondDeformedState[l][1];
      bondDefZ = point->layerAssociatedBondDeformedState[l][2];

      bondLevelVelGradXX = point->layerAssociatedBondLevelVelocityGradient[l][0];
      bondLevelVelGradXY = point->layerAssociatedBondLevelVelocityGradient[l][1];
      bondLevelVelGradXZ = point->layerAssociatedBondLevelVelocityGradient[l][2];
      bondLevelVelGradYX = point->layerAssociatedBondLevelVelocityGradient[l][3];
      bondLevelVelGradYY = point->layerAssociatedBondLevelVelocityGradient[l][4];
      bondLevelVelGradYZ = point->layerAssociatedBondLevelVelocityGradient[l][5];
      bondLevelVelGradZX = point->layerAssociatedBondLevelVelocityGradient[l][6];
      bondLevelVelGradZY = point->layerAssociatedBondLevelVelocityGradient[l][7];
      bondLevelVelGradZZ = point->layerAssociatedBondLevelVelocityGradient[l][8];

      omega = point->influenceState;

      numNeighbors = point->numNeighbors;
      for(n=0; n<numNeighbors; n++,
        omega++, bondDefX++, bondDefY++, bondDefZ++,
        bondLevelVelGradXX++, bondLevelVelGradXY++, bondLevelVelGradXZ++, 
        bondLevelVelGradYX++, bondLevelVelGradYY++, bondLevelVelGradYZ++,
        bondLevelVelGradZX++, bondLevelVelGradZY++, bondLevelVelGradZZ++){

        neigh = point->neighbors[n];
        neighbor = &par->puntos[neigh];

        // only do the calculations for the unbroken bonds
        if(*omega > 0.0){

          neighborThickness = neighbor->shellThickness;
          neighborCoord = neighbor->currentCoord[1];
          neighborVel = neighbor->velocity[1];
          neighborN = neighbor->normalVector;
          neighborNdot = neighbor->normalDotVector[1];
          neighborVelGrad = neighbor->layerAssociatedVelocityGradient[l];

          // The defState is the relative difference in current positions of the nodes at
          // each end of a bond. i.e., x_j - x_i
          *bondDefX = *(neighborCoord+0) - *(coord+0);
          *bondDefX += xi3/2.0 * ( *neighborThickness * *(neighborN+0) - *thickness * *(normal+0) );
          *bondDefY = *(neighborCoord+1) - *(coord+1);
          *bondDefY += xi3/2.0 * ( *neighborThickness * *(neighborN+1) - *thickness * *(normal+1) );
          *bondDefZ = *(neighborCoord+2) - *(coord+2);
          *bondDefZ += xi3/2.0 * ( *neighborThickness * *(neighborN+2) - *thickness * *(normal+2) );

          deformedBondX = *bondDefX; deformedBondY = *bondDefY; deformedBondZ = *bondDefZ; 
          deformedBondLengthSq = deformedBondX*deformedBondX +
                                 deformedBondY*deformedBondY +
                                 deformedBondZ*deformedBondZ;

          // The velState is the relative difference in velocities of the nodes at
          // each end of a bond. i.e., v_j - v_i
          velStateX = *(neighborVel+0) - *(vel+0);
          velStateX += xi3/2.0 * ( *neighborThickness * *(neighborNdot+0) - *thickness * *(ndot+0) );
          velStateY = *(neighborVel+1) - *(vel+1);
          velStateY += xi3/2.0 * ( *neighborThickness * *(neighborNdot+1) - *thickness * *(ndot+1) );
          velStateZ = *(neighborVel+2) - *(vel+2);
          velStateZ += xi3/2.0 * ( *neighborThickness * *(neighborNdot+2) - *thickness * *(ndot+2) );

          // average of the two points 
          for(i=0; i<9; i++)
            *(meanVelGrad+i) = 0.5 * (*(velGrad+i) + *(neighborVelGrad+i));
          
          temp = *(meanVelGrad+0) * deformedBondX + *(meanVelGrad+1) * deformedBondY + *(meanVelGrad+2) * deformedBondZ;
          *bondLevelVelGradXX = *(meanVelGrad+0) + (velStateX - temp) * deformedBondX/deformedBondLengthSq;
          *bondLevelVelGradXY = *(meanVelGrad+1) + (velStateX - temp) * deformedBondY/deformedBondLengthSq;
          *bondLevelVelGradXZ = *(meanVelGrad+2) + (velStateX - temp) * deformedBondZ/deformedBondLengthSq;

          temp = *(meanVelGrad+3) * deformedBondX + *(meanVelGrad+4) * deformedBondY + *(meanVelGrad+5) * deformedBondZ;
          *bondLevelVelGradYX = *(meanVelGrad+3) + (velStateY - temp) * deformedBondX/deformedBondLengthSq;
          *bondLevelVelGradYY = *(meanVelGrad+4) + (velStateY - temp) * deformedBondY/deformedBondLengthSq;
          *bondLevelVelGradYZ = *(meanVelGrad+5) + (velStateY - temp) * deformedBondZ/deformedBondLengthSq;

          temp = *(meanVelGrad+6) * deformedBondX + *(meanVelGrad+7) * deformedBondY + *(meanVelGrad+8) * deformedBondZ;
          *bondLevelVelGradZX = *(meanVelGrad+6) + (velStateZ - temp) * deformedBondX/deformedBondLengthSq;
          *bondLevelVelGradZY = *(meanVelGrad+7) + (velStateZ - temp) * deformedBondY/deformedBondLengthSq;
          *bondLevelVelGradZZ = *(meanVelGrad+8) + (velStateZ - temp) * deformedBondZ/deformedBondLengthSq;
        }
      }
    }
  }
}


void MatrixMultiply
(
    bool transA,
    bool transB,
    const double alpha,
    const double* A,
    const double* B,
    const int p,
    const int q,
    const int r,
    double* result
)
{
  // This function computes result = alpha * a * b
  // where alpha is a scalar and a and b are rectangular matrices
  // The result has dimension pxr and q is the other dimension of a & b
  // The arguments transA and transB denote whether or not
  // to use the transpose of a and b, respectively.

  // The default ordering is row-major:
  //

  int i, j, k;

  if(!transA && !transB){
    // result = A * B
    // A = [pxq] 
    // B = [qxr]
    
    for(i=0 ; i<p ; ++i){
      for(k=0 ; k<r ; ++k){
        *(result+i*p+k) = 0.0;
        for(j=0 ; j<q ; ++j){
          *(result+i*p+k) += *(A+i*p+j) * *(B+j*q+k);
        }
      }
    }
  }
  else if(transA && !transB){
    // result = AT * B
    // A = [qxp] 
    // B = [qxr]
    
    for(i=0 ; i<p ; ++i){
      for(k=0 ; k<r ; ++k){
        *(result+i*p+k) = 0.0;
        for(j=0 ; j<q ; ++j){
          *(result+i*p+k) += *(A+j*p+i) * *(B+j*q+k);
        }
      }
    }
  }
  else if(!transA && transB){
    // result = A * BT
    // A = [pxq] 
    // B = [rxq]
    
    for(i=0 ; i<p ; ++i){
      for(k=0 ; k<r ; ++k){
        *(result+i*p+k) = 0.0;
        for(j=0 ; j<q ; ++j){
          *(result+i*p+k) += *(A+i*p+j) * *(B+k*q+j);
        }
      }
    }
  }
  else{
    // result = AT * BT
    // A = [qxp] 
    // B = [rxq]
    
    for(i=0 ; i<p ; ++i){
      for(k=0 ; k<r ; ++k){
        *(result+i*p+k) = 0.0;
        for(j=0 ; j<q ; ++j){
          *(result+i*p+k) += *(A+j*p+i) * *(B+k*q+j);
        }
      }
    }
  }

  if(alpha != 1.0){
    for(i=0 ; i<p*r ; ++i)
      *(result+i) *= alpha;
  }
}


#undef  __FUNCT__
#define __FUNCT__ "computeUnrotatedRateOfDeformationAndRotationTensor"
void computeUnrotatedRateOfDeformationAndRotationTensor(PARAMETERS *par) 
{
  //Performs kinematic computations following Flanagan and Taylor (1987), 
  //returns the unrotated rate-of-deformation and rotation tensors.
  //This function computes the node-level values.
  //I only calculate the values on the first layer since it is only used for
  //visualization purposes.

  POINTS *point;
  int in,i;
  double dt = par->timeStep;

  double* eulerianVelGrad;
  double* leftStretchN;
  double* rotTensorN;
  double* leftStretchNP1;
  double* rotTensorNP1;
  double* unrotRateOfDef;

  double rateOfDef[9];
  double spin[9];
  double temp[9];
  double tempInv[9];
  double OmegaTensor[9];
  double QMatrix[9];
  double OmegaTensorSq[9];
  double tempA[9];
  double tempB[9];
  double rateOfStretch[9];

  double determinant;
  double omegaX, omegaY, omegaZ;
  double zX, zY, zZ;
  double wX, wY, wZ;
  double traceV, Omega, OmegaSq, scaleFactor1, scaleFactor2;

  for(in=0;in<par->numPoints;in++){
    point = &par->puntos[in];

    eulerianVelGrad = point->layerAssociatedVelocityGradient[0]; // use the values from first layer
    unrotRateOfDef = point->unrotatedRateOfDeformation;
    leftStretchN = point->leftStretchTensor[0];
    leftStretchNP1 = point->leftStretchTensor[1];
    rotTensorN = point->rotationTensor[0];
    rotTensorNP1 = point->rotationTensor[1];

    // Compute rate-of-deformation tensor, D = 1/2 * (L + Lt)
    *(rateOfDef)   = *(eulerianVelGrad);
    *(rateOfDef+1) = 0.5 * ( *(eulerianVelGrad+1) + *(eulerianVelGrad+3) );
    *(rateOfDef+2) = 0.5 * ( *(eulerianVelGrad+2) + *(eulerianVelGrad+6) );
    *(rateOfDef+3) = *(rateOfDef+1);
    *(rateOfDef+4) = *(eulerianVelGrad+4);
    *(rateOfDef+5) = 0.5 * ( *(eulerianVelGrad+5) + *(eulerianVelGrad+7) );
    *(rateOfDef+6) = *(rateOfDef+2);
    *(rateOfDef+7) = *(rateOfDef+5);
    *(rateOfDef+8) = *(eulerianVelGrad+8);

    // Compute spin tensor, W = 1/2 * (L - Lt)
    *(spin)   = 0.0;
    *(spin+1) = 0.5 * ( *(eulerianVelGrad+1) - *(eulerianVelGrad+3) );
    *(spin+2) = 0.5 * ( *(eulerianVelGrad+2) - *(eulerianVelGrad+6) );
    *(spin+3) = -1.0 * *(spin+1);
    *(spin+4) = 0.0;
    *(spin+5) = 0.5 * ( *(eulerianVelGrad+5) - *(eulerianVelGrad+7) );
    *(spin+6) = -1.0 * *(spin+2);
    *(spin+7) = -1.0 * *(spin+5);
    *(spin+8) = 0.0;
   
    //Following Flanagan & Taylor (T&F) 
    //
    //Find the vector z_i = \epsilon_{ikj} * D_{jm} * V_{mk} (T&F Eq. 13)
    //
    //where \epsilon_{ikj} is the alternator tensor.
    //
    //Components below copied from computer algebra solution to the expansion
    //above
    
    zX = - *(leftStretchN+2) *  *(rateOfDef+3) -  *(leftStretchN+5) *  *(rateOfDef+4) - 
           *(leftStretchN+8) *  *(rateOfDef+5) +  *(leftStretchN+1) *  *(rateOfDef+6) + 
           *(leftStretchN+4) *  *(rateOfDef+7) +  *(leftStretchN+7) *  *(rateOfDef+8);
    zY =   *(leftStretchN+2) *  *(rateOfDef)   +  *(leftStretchN+5) *  *(rateOfDef+1) + 
           *(leftStretchN+8) *  *(rateOfDef+2) -  *(leftStretchN)   *  *(rateOfDef+6) - 
           *(leftStretchN+3) *  *(rateOfDef+7) -  *(leftStretchN+6) *  *(rateOfDef+8);
    zZ = - *(leftStretchN+1) *  *(rateOfDef)   -  *(leftStretchN+4) *  *(rateOfDef+1) - 
           *(leftStretchN+7) *  *(rateOfDef+2) +  *(leftStretchN)   *  *(rateOfDef+3) + 
           *(leftStretchN+3) *  *(rateOfDef+4) +  *(leftStretchN+6) *  *(rateOfDef+5);

    //Find the vector w_i = -1/2 * \epsilon_{ijk} * W_{jk} (T&F Eq. 11)
    wX = 0.5 * ( *(spin+7) - *(spin+5) );
    wY = 0.5 * ( *(spin+2) - *(spin+6) );
    wZ = 0.5 * ( *(spin+3) - *(spin+1) );

    //Find trace(V)
    traceV = *(leftStretchN) + *(leftStretchN+4) + *(leftStretchN+8);

    // Compute (trace(V) * I - V) store in temp
    *(temp)   = traceV - *(leftStretchN);
    *(temp+1) = - *(leftStretchN+1);
    *(temp+2) = - *(leftStretchN+2);
    *(temp+3) = - *(leftStretchN+3);
    *(temp+4) = traceV - *(leftStretchN+4);
    *(temp+5) = - *(leftStretchN+5);
    *(temp+6) = - *(leftStretchN+6);
    *(temp+7) = - *(leftStretchN+7);
    *(temp+8) = traceV - *(leftStretchN+8);

    // Compute the inverse of the temp matrix
    Invert3by3Matrix(temp, determinant, tempInv);

    //Find omega vector, i.e. \omega = w +  (trace(V) I - V)^(-1) * z (T&F Eq. 12)
    omegaX =  wX + *(tempInv)   * zX + *(tempInv+1) * zY + *(tempInv+2) * zZ;
    omegaY =  wY + *(tempInv+3) * zX + *(tempInv+4) * zY + *(tempInv+5) * zZ;
    omegaZ =  wZ + *(tempInv+6) * zX + *(tempInv+7) * zY + *(tempInv+8) * zZ;

    //Find the tensor \Omega_{ij} = \epsilon_{ikj} * w_k (T&F Eq. 10)
    *(OmegaTensor) = 0.0;
    *(OmegaTensor+1) = -omegaZ;
    *(OmegaTensor+2) = omegaY;
    *(OmegaTensor+3) = omegaZ;
    *(OmegaTensor+4) = 0.0;
    *(OmegaTensor+5) = -omegaX;
    *(OmegaTensor+6) = -omegaY;
    *(OmegaTensor+7) = omegaX;
    *(OmegaTensor+8) = 0.0;

    //Increment R with (T&F Eq. 36 and 44) as opposed to solving (T&F 39) this
    //is desirable for accuracy in implicit solves and has no effect on
    //explicit solves (other than a slight decrease in speed).
    //
    // Compute Q with (T&F Eq. 44)
    //
    // Omega^2 = w_i * w_i (T&F Eq. 42)
    OmegaSq = omegaX*omegaX + omegaY*omegaY + omegaZ*omegaZ;
    // Omega = \sqrt{OmegaSq}
    Omega = sqrt(OmegaSq);

    // Avoid a potential divide-by-zero
    if(OmegaSq > 1.e-30){

      // Compute Q = I + sin( dt * Omega ) * OmegaTensor / Omega - (1. - cos(dt * Omega)) * omegaTensor^2 / OmegaSq
      //           = I + scaleFactor1 * OmegaTensor + scaleFactor2 * OmegaTensorSq
      scaleFactor1 = sin(dt*Omega) / Omega;
      scaleFactor2 = -(1.0 - cos(dt*Omega)) / OmegaSq;
      MatrixMultiply(false, false, 1.0, OmegaTensor, OmegaTensor, 3, 3, 3, OmegaTensorSq);
      *(QMatrix)   = 1.0 + scaleFactor1 * *(OmegaTensor)   + scaleFactor2 * *(OmegaTensorSq)   ;
      *(QMatrix+1) =       scaleFactor1 * *(OmegaTensor+1) + scaleFactor2 * *(OmegaTensorSq+1) ;
      *(QMatrix+2) =       scaleFactor1 * *(OmegaTensor+2) + scaleFactor2 * *(OmegaTensorSq+2) ;
      *(QMatrix+3) =       scaleFactor1 * *(OmegaTensor+3) + scaleFactor2 * *(OmegaTensorSq+3) ;
      *(QMatrix+4) = 1.0 + scaleFactor1 * *(OmegaTensor+4) + scaleFactor2 * *(OmegaTensorSq+4) ;
      *(QMatrix+5) =       scaleFactor1 * *(OmegaTensor+5) + scaleFactor2 * *(OmegaTensorSq+5) ;
      *(QMatrix+6) =       scaleFactor1 * *(OmegaTensor+6) + scaleFactor2 * *(OmegaTensorSq+6) ;
      *(QMatrix+7) =       scaleFactor1 * *(OmegaTensor+7) + scaleFactor2 * *(OmegaTensorSq+7) ;
      *(QMatrix+8) = 1.0 + scaleFactor1 * *(OmegaTensor+8) + scaleFactor2 * *(OmegaTensorSq+8) ;

    } else {
      *(QMatrix)   = 1.0 ; *(QMatrix+1) = 0.0 ; *(QMatrix+2) = 0.0 ;
      *(QMatrix+3) = 0.0 ; *(QMatrix+4) = 1.0 ; *(QMatrix+5) = 0.0 ;
      *(QMatrix+6) = 0.0 ; *(QMatrix+7) = 0.0 ; *(QMatrix+8) = 1.0 ;
    };

    // Compute R_STEP_NP1 = QMatrix * R_STEP_N (T&F Eq. 36)
    MatrixMultiply(false, false, 1.0, QMatrix, rotTensorN, 3, 3, 3, rotTensorNP1);

    // Compute rate of stretch, Vdot = L*V - V*Omega
    // First tempA = L*V, 
    MatrixMultiply(false, false, 1.0, eulerianVelGrad, leftStretchN, 3, 3, 3, tempA);

    // tempB = V*Omega
    MatrixMultiply(false, false, 1.0, leftStretchN, OmegaTensor, 3, 3, 3, tempB);

    //Vdot = tempA - tempB
    for(i=0 ; i<9 ; ++i)
      *(rateOfStretch+i) = *(tempA+i) - *(tempB+i);

    //V_STEP_NP1 = V_STEP_N + dt*Vdot
    for(i=0 ; i<9 ; ++i)
      *(leftStretchNP1+i) = *(leftStretchN+i) + dt * *(rateOfStretch+i);

    // Compute the unrotated rate-of-deformation, d, i.e., temp = D * R
    MatrixMultiply(false, false, 1.0, rateOfDef, rotTensorNP1, 3, 3, 3, temp);

    // d = Rt * temp
    MatrixMultiply(true, false, 1.0, rotTensorNP1, temp, 3, 3, 3, unrotRateOfDef);
  }
}


#undef  __FUNCT__
#define __FUNCT__ "computeBondLevelUnrotatedRateOfDeformationAndRotationTensor"
void computeBondLevelUnrotatedRateOfDeformationAndRotationTensor(PARAMETERS *par) 
{
  POINTS *point;
  int i,l,in,n;
  int numNeighbors;
  double dt = par->timeStep;

  double* velGradXX;
  double* velGradXY;
  double* velGradXZ;
  double* velGradYX;
  double* velGradYY;
  double* velGradYZ;
  double* velGradZX;
  double* velGradZY;
  double* velGradZZ;
  double* leftStretchXXN;
  double* leftStretchXYN;
  double* leftStretchXZN;
  double* leftStretchYXN;
  double* leftStretchYYN;
  double* leftStretchYZN;
  double* leftStretchZXN;
  double* leftStretchZYN;
  double* leftStretchZZN;
  double* rotTensorXXN;
  double* rotTensorXYN;
  double* rotTensorXZN;
  double* rotTensorYXN;
  double* rotTensorYYN;
  double* rotTensorYZN;
  double* rotTensorZXN;
  double* rotTensorZYN;
  double* rotTensorZZN;

  double* leftStretchXXNP1;
  double* leftStretchXYNP1;
  double* leftStretchXZNP1;
  double* leftStretchYXNP1;
  double* leftStretchYYNP1;
  double* leftStretchYZNP1;
  double* leftStretchZXNP1;
  double* leftStretchZYNP1;
  double* leftStretchZZNP1;
  double* rotTensorXXNP1;
  double* rotTensorXYNP1;
  double* rotTensorXZNP1;
  double* rotTensorYXNP1;
  double* rotTensorYYNP1;
  double* rotTensorYZNP1;
  double* rotTensorZXNP1;
  double* rotTensorZYNP1;
  double* rotTensorZZNP1;
  double* unrotRateOfDefXX;
  double* unrotRateOfDefXY;
  double* unrotRateOfDefXZ;
  double* unrotRateOfDefYX;
  double* unrotRateOfDefYY;
  double* unrotRateOfDefYZ;
  double* unrotRateOfDefZX;
  double* unrotRateOfDefZY;
  double* unrotRateOfDefZZ;

  double eulerianVelGrad[9];
  double leftStretchN[9];
  double rotTensorN[9];
  double leftStretchNP1[9];
  double rotTensorNP1[9];
  double unrotRateOfDef[9];

  double rateOfDef[9];
  double spin[9];
  double temp[9];
  double tempInv[9];
  double OmegaTensor[9];
  double QMatrix[9];
  double OmegaTensorSq[9];
  double tempA[9];
  double tempB[9];
  double rateOfStretch[9];

  double determinant;
  double omegaX, omegaY, omegaZ;
  double zX, zY, zZ;
  double wX, wY, wZ;
  double traceV, Omega, OmegaSq, scaleFactor1, scaleFactor2;

  for(in=0;in<par->numPoints;in++){

    point = &par->puntos[in];

    for(l=0; l<par->numLayers; l++){

      velGradXX = point->layerAssociatedBondLevelVelocityGradient[l][0];
      velGradXY = point->layerAssociatedBondLevelVelocityGradient[l][1];
      velGradXZ = point->layerAssociatedBondLevelVelocityGradient[l][2];
      velGradYX = point->layerAssociatedBondLevelVelocityGradient[l][3];
      velGradYY = point->layerAssociatedBondLevelVelocityGradient[l][4];
      velGradYZ = point->layerAssociatedBondLevelVelocityGradient[l][5];
      velGradZX = point->layerAssociatedBondLevelVelocityGradient[l][6];
      velGradZY = point->layerAssociatedBondLevelVelocityGradient[l][7];
      velGradZZ = point->layerAssociatedBondLevelVelocityGradient[l][8];

      leftStretchXXN = point->layerAssociatedBondLevelLeftStretchTensor[0][l][0];
      leftStretchXYN = point->layerAssociatedBondLevelLeftStretchTensor[0][l][1];
      leftStretchXZN = point->layerAssociatedBondLevelLeftStretchTensor[0][l][2];
      leftStretchYXN = point->layerAssociatedBondLevelLeftStretchTensor[0][l][3];
      leftStretchYYN = point->layerAssociatedBondLevelLeftStretchTensor[0][l][4];
      leftStretchYZN = point->layerAssociatedBondLevelLeftStretchTensor[0][l][5];
      leftStretchZXN = point->layerAssociatedBondLevelLeftStretchTensor[0][l][6];
      leftStretchZYN = point->layerAssociatedBondLevelLeftStretchTensor[0][l][7];
      leftStretchZZN = point->layerAssociatedBondLevelLeftStretchTensor[0][l][8];

      leftStretchXXNP1 = point->layerAssociatedBondLevelLeftStretchTensor[1][l][0];
      leftStretchXYNP1 = point->layerAssociatedBondLevelLeftStretchTensor[1][l][1];
      leftStretchXZNP1 = point->layerAssociatedBondLevelLeftStretchTensor[1][l][2];
      leftStretchYXNP1 = point->layerAssociatedBondLevelLeftStretchTensor[1][l][3];
      leftStretchYYNP1 = point->layerAssociatedBondLevelLeftStretchTensor[1][l][4];
      leftStretchYZNP1 = point->layerAssociatedBondLevelLeftStretchTensor[1][l][5];
      leftStretchZXNP1 = point->layerAssociatedBondLevelLeftStretchTensor[1][l][6];
      leftStretchZYNP1 = point->layerAssociatedBondLevelLeftStretchTensor[1][l][7];
      leftStretchZZNP1 = point->layerAssociatedBondLevelLeftStretchTensor[1][l][8];

      rotTensorXXN = point->layerAssociatedBondLevelRotationTensor[0][l][0];
      rotTensorXYN = point->layerAssociatedBondLevelRotationTensor[0][l][1];
      rotTensorXZN = point->layerAssociatedBondLevelRotationTensor[0][l][2];
      rotTensorYXN = point->layerAssociatedBondLevelRotationTensor[0][l][3];
      rotTensorYYN = point->layerAssociatedBondLevelRotationTensor[0][l][4];
      rotTensorYZN = point->layerAssociatedBondLevelRotationTensor[0][l][5];
      rotTensorZXN = point->layerAssociatedBondLevelRotationTensor[0][l][6];
      rotTensorZYN = point->layerAssociatedBondLevelRotationTensor[0][l][7];
      rotTensorZZN = point->layerAssociatedBondLevelRotationTensor[0][l][8];

      rotTensorXXNP1 = point->layerAssociatedBondLevelRotationTensor[1][l][0];
      rotTensorXYNP1 = point->layerAssociatedBondLevelRotationTensor[1][l][1];
      rotTensorXZNP1 = point->layerAssociatedBondLevelRotationTensor[1][l][2];
      rotTensorYXNP1 = point->layerAssociatedBondLevelRotationTensor[1][l][3];
      rotTensorYYNP1 = point->layerAssociatedBondLevelRotationTensor[1][l][4];
      rotTensorYZNP1 = point->layerAssociatedBondLevelRotationTensor[1][l][5];
      rotTensorZXNP1 = point->layerAssociatedBondLevelRotationTensor[1][l][6];
      rotTensorZYNP1 = point->layerAssociatedBondLevelRotationTensor[1][l][7];
      rotTensorZZNP1 = point->layerAssociatedBondLevelRotationTensor[1][l][8];

      unrotRateOfDefXX = point->layerAssociatedBondLevelUnrotatedRateOfDeformation[l][0];
      unrotRateOfDefXY = point->layerAssociatedBondLevelUnrotatedRateOfDeformation[l][1];
      unrotRateOfDefXZ = point->layerAssociatedBondLevelUnrotatedRateOfDeformation[l][2];
      unrotRateOfDefYX = point->layerAssociatedBondLevelUnrotatedRateOfDeformation[l][3];
      unrotRateOfDefYY = point->layerAssociatedBondLevelUnrotatedRateOfDeformation[l][4];
      unrotRateOfDefYZ = point->layerAssociatedBondLevelUnrotatedRateOfDeformation[l][5];
      unrotRateOfDefZX = point->layerAssociatedBondLevelUnrotatedRateOfDeformation[l][6];
      unrotRateOfDefZY = point->layerAssociatedBondLevelUnrotatedRateOfDeformation[l][7];
      unrotRateOfDefZZ = point->layerAssociatedBondLevelUnrotatedRateOfDeformation[l][8];

      numNeighbors = point->numNeighbors;
      for(n=0; n<numNeighbors; n++,
        velGradXX++, velGradXY++, velGradXZ++, 
        velGradYX++, velGradYY++, velGradYZ++, 
        velGradZX++, velGradZY++, velGradZZ++,
        leftStretchXXN++, leftStretchXYN++, leftStretchXZN++, 
        leftStretchYXN++, leftStretchYYN++, leftStretchYZN++, 
        leftStretchZXN++, leftStretchZYN++, leftStretchZZN++,
        rotTensorXXN++, rotTensorXYN++, rotTensorXZN++,
        rotTensorYXN++, rotTensorYYN++, rotTensorYZN++,
        rotTensorZXN++, rotTensorZYN++, rotTensorZZN++,
        leftStretchXXNP1++, leftStretchXYNP1++, leftStretchXZNP1++, 
        leftStretchYXNP1++, leftStretchYYNP1++, leftStretchYZNP1++, 
        leftStretchZXNP1++, leftStretchZYNP1++, leftStretchZZNP1++,
        rotTensorXXNP1++, rotTensorXYNP1++, rotTensorXZNP1++,
        rotTensorYXNP1++, rotTensorYYNP1++, rotTensorYZNP1++,
        rotTensorZXNP1++, rotTensorZYNP1++, rotTensorZZNP1++,
        unrotRateOfDefXX++, unrotRateOfDefXY++, unrotRateOfDefXZ++,
        unrotRateOfDefYX++, unrotRateOfDefYY++, unrotRateOfDefYZ++,
        unrotRateOfDefZX++, unrotRateOfDefZY++, unrotRateOfDefZZ++){

        // Store in a tensor form 
        *(eulerianVelGrad+0) = *velGradXX; *(eulerianVelGrad+1) = *velGradXY; *(eulerianVelGrad+2) = *velGradXZ;
        *(eulerianVelGrad+3) = *velGradYX; *(eulerianVelGrad+4) = *velGradYY; *(eulerianVelGrad+5) = *velGradYZ;
        *(eulerianVelGrad+6) = *velGradZX; *(eulerianVelGrad+7) = *velGradZY; *(eulerianVelGrad+8) = *velGradZZ;
        *(leftStretchN+0) = *leftStretchXXN; *(leftStretchN+1) = *leftStretchXYN; *(leftStretchN+2) = *leftStretchXZN;
        *(leftStretchN+3) = *leftStretchYXN; *(leftStretchN+4) = *leftStretchYYN; *(leftStretchN+5) = *leftStretchYZN;
        *(leftStretchN+6) = *leftStretchZXN; *(leftStretchN+7) = *leftStretchZYN; *(leftStretchN+8) = *leftStretchZZN;
        *(rotTensorN+0) = *rotTensorXXN; *(rotTensorN+1) = *rotTensorXYN; *(rotTensorN+2) = *rotTensorXZN;
        *(rotTensorN+3) = *rotTensorYXN; *(rotTensorN+4) = *rotTensorYYN; *(rotTensorN+5) = *rotTensorYZN;
        *(rotTensorN+6) = *rotTensorZXN; *(rotTensorN+7) = *rotTensorZYN; *(rotTensorN+8) = *rotTensorZZN;

        // Compute rate-of-deformation tensor, D = 1/2 * (L + Lt)
        *(rateOfDef)   = *(eulerianVelGrad);
        *(rateOfDef+1) = 0.5 * ( *(eulerianVelGrad+1) + *(eulerianVelGrad+3) );
        *(rateOfDef+2) = 0.5 * ( *(eulerianVelGrad+2) + *(eulerianVelGrad+6) );
        *(rateOfDef+3) = *(rateOfDef+1);
        *(rateOfDef+4) = *(eulerianVelGrad+4);
        *(rateOfDef+5) = 0.5 * ( *(eulerianVelGrad+5) + *(eulerianVelGrad+7) );
        *(rateOfDef+6) = *(rateOfDef+2);
        *(rateOfDef+7) = *(rateOfDef+5);
        *(rateOfDef+8) = *(eulerianVelGrad+8);

        // Compute spin tensor, W = 1/2 * (L - Lt)
        *(spin)   = 0.0;
        *(spin+1) = 0.5 * ( *(eulerianVelGrad+1) - *(eulerianVelGrad+3) );
        *(spin+2) = 0.5 * ( *(eulerianVelGrad+2) - *(eulerianVelGrad+6) );
        *(spin+3) = -1.0 * *(spin+1);
        *(spin+4) = 0.0;
        *(spin+5) = 0.5 * ( *(eulerianVelGrad+5) - *(eulerianVelGrad+7) );
        *(spin+6) = -1.0 * *(spin+2);
        *(spin+7) = -1.0 * *(spin+5);
        *(spin+8) = 0.0;
       
        //Following Flanagan & Taylor (T&F) 
        //
        //Find the vector z_i = \epsilon_{ikj} * D_{jm} * V_{mk} (T&F Eq. 13)
        //
        //where \epsilon_{ikj} is the alternator tensor.
        //
        //Components below copied from computer algebra solution to the expansion
        //above
        zX = - *(leftStretchN+2) *  *(rateOfDef+3) -  *(leftStretchN+5) *  *(rateOfDef+4) - 
               *(leftStretchN+8) *  *(rateOfDef+5) +  *(leftStretchN+1) *  *(rateOfDef+6) + 
               *(leftStretchN+4) *  *(rateOfDef+7) +  *(leftStretchN+7) *  *(rateOfDef+8);
        zY =   *(leftStretchN+2) *  *(rateOfDef)   +  *(leftStretchN+5) *  *(rateOfDef+1) + 
               *(leftStretchN+8) *  *(rateOfDef+2) -  *(leftStretchN)   *  *(rateOfDef+6) - 
               *(leftStretchN+3) *  *(rateOfDef+7) -  *(leftStretchN+6) *  *(rateOfDef+8);
        zZ = - *(leftStretchN+1) *  *(rateOfDef)   -  *(leftStretchN+4) *  *(rateOfDef+1) - 
               *(leftStretchN+7) *  *(rateOfDef+2) +  *(leftStretchN)   *  *(rateOfDef+3) + 
               *(leftStretchN+3) *  *(rateOfDef+4) +  *(leftStretchN+6) *  *(rateOfDef+5);

        //Find the vector w_i = -1/2 * \epsilon_{ijk} * W_{jk} (T&F Eq. 11)
        wX = 0.5 * ( *(spin+7) - *(spin+5) );
        wY = 0.5 * ( *(spin+2) - *(spin+6) );
        wZ = 0.5 * ( *(spin+3) - *(spin+1) );

        //Find trace(V)
        traceV = *(leftStretchN) + *(leftStretchN+4) + *(leftStretchN+8);

        // Compute (trace(V) * I - V) store in temp
        *(temp)   = traceV - *(leftStretchN);
        *(temp+1) = - *(leftStretchN+1);
        *(temp+2) = - *(leftStretchN+2);
        *(temp+3) = - *(leftStretchN+3);
        *(temp+4) = traceV - *(leftStretchN+4);
        *(temp+5) = - *(leftStretchN+5);
        *(temp+6) = - *(leftStretchN+6);
        *(temp+7) = - *(leftStretchN+7);
        *(temp+8) = traceV - *(leftStretchN+8);

        // Compute the inverse of the temp matrix
        Invert3by3Matrix(temp, determinant, tempInv);

        //Find omega vector, i.e. \omega = w +  (trace(V) I - V)^(-1) * z (T&F Eq. 12)
        omegaX =  wX + *(tempInv)   * zX + *(tempInv+1) * zY + *(tempInv+2) * zZ;
        omegaY =  wY + *(tempInv+3) * zX + *(tempInv+4) * zY + *(tempInv+5) * zZ;
        omegaZ =  wZ + *(tempInv+6) * zX + *(tempInv+7) * zY + *(tempInv+8) * zZ;

        //Find the tensor \Omega_{ij} = \epsilon_{ikj} * w_k (T&F Eq. 10)
        *(OmegaTensor) = 0.0;
        *(OmegaTensor+1) = -omegaZ;
        *(OmegaTensor+2) = omegaY;
        *(OmegaTensor+3) = omegaZ;
        *(OmegaTensor+4) = 0.0;
        *(OmegaTensor+5) = -omegaX;
        *(OmegaTensor+6) = -omegaY;
        *(OmegaTensor+7) = omegaX;
        *(OmegaTensor+8) = 0.0;

        //Increment R with (T&F Eq. 36 and 44) as opposed to solving (T&F 39) this
        //is desirable for accuracy in implicit solves and has no effect on
        //explicit solves (other than a slight decrease in speed).
        //
        // Compute Q with (T&F Eq. 44)
        //
        // Omega^2 = w_i * w_i (T&F Eq. 42)
        OmegaSq = omegaX*omegaX + omegaY*omegaY + omegaZ*omegaZ;
        // Omega = \sqrt{OmegaSq}
        Omega = sqrt(OmegaSq);

        // Avoid a potential divide-by-zero
        if(OmegaSq > 1.e-30){

          // Compute Q = I + sin( dt * Omega ) * OmegaTensor / Omega - (1. - cos(dt * Omega)) * omegaTensor^2 / OmegaSq
          //           = I + scaleFactor1 * OmegaTensor + scaleFactor2 * OmegaTensorSq
          scaleFactor1 = sin(dt*Omega) / Omega;
          scaleFactor2 = -(1.0 - cos(dt*Omega)) / OmegaSq;
          MatrixMultiply(false, false, 1.0, OmegaTensor, OmegaTensor, 3, 3, 3, OmegaTensorSq);
          *(QMatrix)   = 1.0 + scaleFactor1 * *(OmegaTensor)   + scaleFactor2 * *(OmegaTensorSq)   ;
          *(QMatrix+1) =       scaleFactor1 * *(OmegaTensor+1) + scaleFactor2 * *(OmegaTensorSq+1) ;
          *(QMatrix+2) =       scaleFactor1 * *(OmegaTensor+2) + scaleFactor2 * *(OmegaTensorSq+2) ;
          *(QMatrix+3) =       scaleFactor1 * *(OmegaTensor+3) + scaleFactor2 * *(OmegaTensorSq+3) ;
          *(QMatrix+4) = 1.0 + scaleFactor1 * *(OmegaTensor+4) + scaleFactor2 * *(OmegaTensorSq+4) ;
          *(QMatrix+5) =       scaleFactor1 * *(OmegaTensor+5) + scaleFactor2 * *(OmegaTensorSq+5) ;
          *(QMatrix+6) =       scaleFactor1 * *(OmegaTensor+6) + scaleFactor2 * *(OmegaTensorSq+6) ;
          *(QMatrix+7) =       scaleFactor1 * *(OmegaTensor+7) + scaleFactor2 * *(OmegaTensorSq+7) ;
          *(QMatrix+8) = 1.0 + scaleFactor1 * *(OmegaTensor+8) + scaleFactor2 * *(OmegaTensorSq+8) ;

        } else {
          *(QMatrix)   = 1.0 ; *(QMatrix+1) = 0.0 ; *(QMatrix+2) = 0.0 ;
          *(QMatrix+3) = 0.0 ; *(QMatrix+4) = 1.0 ; *(QMatrix+5) = 0.0 ;
          *(QMatrix+6) = 0.0 ; *(QMatrix+7) = 0.0 ; *(QMatrix+8) = 1.0 ;
        };

        // Compute R_STEP_NP1 = QMatrix * R_STEP_N (T&F Eq. 36)
        MatrixMultiply(false, false, 1.0, QMatrix, rotTensorN, 3, 3, 3, rotTensorNP1);

        // Compute rate of stretch, Vdot = L*V - V*Omega
        // First tempA = L*V, 
        MatrixMultiply(false, false, 1.0, eulerianVelGrad, leftStretchN, 3, 3, 3, tempA);

        // tempB = V*Omega
        MatrixMultiply(false, false, 1.0, leftStretchN, OmegaTensor, 3, 3, 3, tempB);

        //Vdot = tempA - tempB
        for(i=0 ; i<9 ; ++i)
          *(rateOfStretch+i) = *(tempA+i) - *(tempB+i);

        //V_STEP_NP1 = V_STEP_N + dt*Vdot
        for(i=0 ; i<9 ; ++i)
          *(leftStretchNP1+i) = *(leftStretchN+i) + dt * *(rateOfStretch+i);

        // Compute the unrotated rate-of-deformation, d, i.e., temp = D * R
        MatrixMultiply(false, false, 1.0, rateOfDef, rotTensorNP1, 3, 3, 3, temp);

        // d = Rt * temp
        MatrixMultiply(true, false, 1.0, rotTensorNP1, temp, 3, 3, 3, unrotRateOfDef);

        // Store back in element-wise format
        *leftStretchXXNP1 = *(leftStretchNP1+0); *leftStretchXYNP1 = *(leftStretchNP1+1); *leftStretchXZNP1 = *(leftStretchNP1+2);
        *leftStretchYXNP1 = *(leftStretchNP1+3); *leftStretchYYNP1 = *(leftStretchNP1+4); *leftStretchYZNP1 = *(leftStretchNP1+5);
        *leftStretchZXNP1 = *(leftStretchNP1+6); *leftStretchZYNP1 = *(leftStretchNP1+7); *leftStretchZZNP1 = *(leftStretchNP1+8);
        *rotTensorXXNP1 = *(rotTensorNP1+0); *rotTensorXYNP1 = *(rotTensorNP1+1); *rotTensorXZNP1 = *(rotTensorNP1+2);
        *rotTensorYXNP1 = *(rotTensorNP1+3); *rotTensorYYNP1 = *(rotTensorNP1+4); *rotTensorYZNP1 = *(rotTensorNP1+5);
        *rotTensorZXNP1 = *(rotTensorNP1+6); *rotTensorZYNP1 = *(rotTensorNP1+7); *rotTensorZZNP1 = *(rotTensorNP1+8);
        *unrotRateOfDefXX = *(unrotRateOfDef+0); *unrotRateOfDefXY = *(unrotRateOfDef+1); *unrotRateOfDefXZ = *(unrotRateOfDef+2);
        *unrotRateOfDefYX = *(unrotRateOfDef+3); *unrotRateOfDefYY = *(unrotRateOfDef+4); *unrotRateOfDefYZ = *(unrotRateOfDef+5);
        *unrotRateOfDefZX = *(unrotRateOfDef+6); *unrotRateOfDefZY = *(unrotRateOfDef+7); *unrotRateOfDefZZ = *(unrotRateOfDef+8);
      }
    }
  }
}


#undef  __FUNCT__
#define __FUNCT__ "radialReturnSaturationExponential"
void radialReturnSaturationExponential
(
    const double* unrotatedRateOfDeformation, 
    const double* cauchyStressN, 
    double* cauchyStressNP1, 
    double* vonMisesStress,
    const double* equivalentPlasticStrainN,
    double* equivalentPlasticStrainNP1,
    double *C33,
    const double bulkMod,
    const double shearMod,
    const double initialYieldStress,
    const double saturatedYieldStress,
    const double exponentialConstant,
    const double linearConstant,
    const double dt
)
{
  const double* rateOfDef = unrotatedRateOfDeformation;
  const double* stressN = cauchyStressN;
  double* stressNP1 = cauchyStressNP1;
  double* vmStress = vonMisesStress;
  const double* eqpsN = equivalentPlasticStrainN;
  double* eqpsNP1 = equivalentPlasticStrainNP1;

  double strainInc[9];
  double deviatoricStrainInc[9];

  double deviatoricStressNP1[9];
  double deviatoricStressMagnitudeNP1;
  double flowVector[9];
  double vmStress_trial;

  double dilatationInc;
  double sphericalStressNP1;
  double tempScalar;
  double yieldFunctionVal;
  double hardenedYieldStress;
  double yieldFunctionDerivative;
  double deltaPlasticIncrement;
  double plasticIncrement;

  double tol = 1.0e-16 * bulkMod;
  int counter;
  int maxIterations = 100;

  double hardMod;
  double beta, gamma_hat;

  //strainInc = dt * rateOfDef
  for(int i = 0; i < 9; i++){
    strainInc[i] = *(rateOfDef+i)*dt;
    deviatoricStrainInc[i] = strainInc[i];
  }

  //dilatation
  dilatationInc = strainInc[0] + strainInc[4] + strainInc[8];

  //deviatoric strain
  deviatoricStrainInc[0] -= dilatationInc/3.0;
  deviatoricStrainInc[4] -= dilatationInc/3.0;
  deviatoricStrainInc[8] -= dilatationInc/3.0;

  //Compute an elastic ``trail stress''
  for(int i = 0; i < 9; i++){
    *(stressNP1+i) = *(stressN+i) + deviatoricStrainInc[i]*2.0*shearMod;
  }
  *(stressNP1+0) += bulkMod*dilatationInc;
  *(stressNP1+4) += bulkMod*dilatationInc;
  *(stressNP1+8) += bulkMod*dilatationInc;

  sphericalStressNP1 = (*(stressNP1) + *(stressNP1+4) + *(stressNP1+8))/3.0;

  // Compute the ``trial'' von Mises stress
  for(int i = 0; i < 9; i++){
    deviatoricStressNP1[i] = *(stressNP1+i);
  }
  deviatoricStressNP1[0] -= sphericalStressNP1;
  deviatoricStressNP1[4] -= sphericalStressNP1;
  deviatoricStressNP1[8] -= sphericalStressNP1;

  // Compute \S_ij * \S_ij
  tempScalar = 0.0;
  for(int i = 0; i < 9; i++)
    tempScalar += deviatoricStressNP1[i] * deviatoricStressNP1[i];

  // Avoid divide-by-zero
  deviatoricStressMagnitudeNP1 = std::max(1.0e-20,sqrt(tempScalar));

  vmStress_trial = sqrt(3.0/2.0*tempScalar);

  //Evaluate the current yield stress
  hardenedYieldStress = initialYieldStress + (saturatedYieldStress - initialYieldStress) * ( 1.0 - exp(-exponentialConstant * *eqpsN)) + linearConstant * *eqpsN;

  //Evaluate the yield function
  yieldFunctionVal = vmStress_trial - hardenedYieldStress;

  // Elastic or plastic?
  if(yieldFunctionVal < 0.0){

    // elastic it is
    *eqpsNP1 = *eqpsN;

    *vmStress = vmStress_trial;

    // Calculate C33 from the consistent tangent
    *C33 = bulkMod + 4.0/3.0*shearMod;

  } else {
    // The step is plastic and we need to update the yield surface
    // (hardening) and return to that.
    
    counter = 0;
    plasticIncrement = 0.0;
    *eqpsNP1 = *eqpsN;

    // Flow vector's direction does not change in radial return: this is N = df/dsigma
    for (int i = 0; i < 9; i++)
      flowVector[i] = deviatoricStressNP1[i] / deviatoricStressMagnitudeNP1;

    // Newton iteration to project stress to the yield surface
    while(fabs(yieldFunctionVal) > tol){

      counter++;

      if(counter > maxIterations){
        std::cout << "radialReturnSaturationExponential: Too many iterations in the plastic radial return." << std::endl;
        exit(1);
      }

      // evaluate the hardening modulus 
      hardMod = exponentialConstant * (saturatedYieldStress - initialYieldStress) * exp(-exponentialConstant * *eqpsNP1) + linearConstant ;

      // --- (1) COMPUTE INCREMENT IN THE PLASTIC MULTIPLIER ---
      // Residual derivative
      yieldFunctionDerivative = 3.0*shearMod + hardMod;

      // Plastic increment
      deltaPlasticIncrement = (yieldFunctionVal/yieldFunctionDerivative);

      plasticIncrement += deltaPlasticIncrement;

      // --- (2) UPDATE --- 
      // Update eq. plastic strain
      *eqpsNP1 += deltaPlasticIncrement;

      // Update stress
      tempScalar = sqrt(6.0) * shearMod * deltaPlasticIncrement;
      for (int i = 0; i < 9; i++){
        *(stressNP1+i) -= tempScalar * flowVector[i];
      }

      //Evaluate the current yield stress
      hardenedYieldStress = initialYieldStress + (saturatedYieldStress - initialYieldStress) * ( 1.0 - exp(-exponentialConstant * *eqpsNP1)) + linearConstant * *eqpsNP1;

      //Evaluate the yield function
      yieldFunctionVal = vmStress_trial - 3.0 * shearMod * plasticIncrement - hardenedYieldStress;
    } 

    // Update von Mises stress
    *vmStress = vmStress_trial - 3.0 * shearMod * plasticIncrement;
     
    // 33 component of consistent tangent
    beta = hardenedYieldStress / vmStress_trial;
    gamma_hat = (1.0/(1.0+(hardMod/(3.0*shearMod)))) - (1.0-beta);
    *C33 = bulkMod + 4.0/3.0 * shearMod * beta 
          - 2.0 * shearMod * gamma_hat * flowVector[8] * flowVector[8];
  }
}


#undef  __FUNCT__
#define __FUNCT__ "updateIsotropicHardeningCauchyStress"
void updateIsotropicHardeningCauchyStress(PARAMETERS *par) 
{
  //Rotate the unrotated Cauchy stress to the current configuration using
  //the rotation tensor
  //This function computes the node-level values (only for visualization).

  POINTS *point;
  int in,i,j;

  double bulkMod = par->bulkModulus;
  double shearMod = par->shearModulus;
  double initialYieldStress = par->initialYieldStress;
  double saturatedYieldStress = par->saturatedYieldStress;
  double exponentialConstant = par->hardeningExponentialConstant;
  double linearConstant = par->hardeningLinearConstant;
  double dt = par->timeStep;

  double* rateOfDef;
  double* stressN;
  double* stressNP1;
  double* eqpsN;
  double* eqpsNP1;
  double* vmStress;
  double* minPrincipalStress;
  double* maxPrincipalStress;

  double C33;
  double tol = 1.0e-16 * bulkMod;
  int counter;
  int maxIterations = 100;

  double V[3][3], d[3], e[3];

  for(in=0;in<par->numPoints;in++){
    point = &par->puntos[in];

    rateOfDef = point->unrotatedRateOfDeformation;
    stressN = point->unrotatedCauchyStress[0];
    stressNP1 = point->unrotatedCauchyStress[1];
    eqpsN = &(point->equivalentPlasticStrain[0]);
    eqpsNP1 = &(point->equivalentPlasticStrain[1]);
    vmStress = &(point->vonMisesStress);
    minPrincipalStress = &(point->minPrincipalStress);
    maxPrincipalStress = &(point->maxPrincipalStress);

    radialReturnSaturationExponential(rateOfDef, 
                                      stressN, 
                                      stressNP1, 
                                      vmStress,
                                      eqpsN, 
                                      eqpsNP1, 
                                      &C33,
                                      bulkMod, 
                                      shearMod,
                                      initialYieldStress, 
                                      saturatedYieldStress,
                                      exponentialConstant, 
                                      linearConstant, 
                                      dt);
  

    // Enforce the zero-through-thickness stress 
    // A D33 is added to the system such that make stress33 zero
    counter = 0;
    while(fabs(*(stressNP1+8)) > tol){

      if(counter > maxIterations){
        std::cout << "updateIsotropicHardeningCauchyStress: Too many iterations to enforce zero through-thickness-stress." << std::endl;
        exit(1);
      }

      *(rateOfDef+8) -= *(stressNP1+8) / (C33*dt);

      radialReturnSaturationExponential(rateOfDef, 
                                        stressN, 
                                        stressNP1, 
                                        vmStress,
                                        eqpsN, 
                                        eqpsNP1, 
                                        &C33,
                                        bulkMod, 
                                        shearMod,
                                        initialYieldStress, 
                                        saturatedYieldStress,
                                        exponentialConstant, 
                                        linearConstant, 
                                        dt);
      counter++;
    }

    // compute eigenvalues
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        V[i][j] = *(stressNP1+3*i+j);
      }
    }
    tred2(V, d, e); // diagonalize V

    // find the max and min of eigenvalues
    *minPrincipalStress = d[0];
    *maxPrincipalStress = d[0];
    for(i=1; i<3; i++){
      if(d[i] > *maxPrincipalStress)
        *maxPrincipalStress = d[i];
      if(d[i] < *minPrincipalStress)
        *minPrincipalStress = d[i];
    }
  }
}


#undef  __FUNCT__
#define __FUNCT__ "updateBondLevelIsotropicHardeningCauchyStress"
void updateBondLevelIsotropicHardeningCauchyStress(PARAMETERS *par) 
{
  POINTS *point;
  int in,i,j,l,n;
  int numNeighbors;

  double bulkMod = par->bulkModulus;
  double shearMod = par->shearModulus;
  double initialYieldStress = par->initialYieldStress;
  double saturatedYieldStress = par->saturatedYieldStress;
  double exponentialConstant = par->hardeningExponentialConstant;
  double linearConstant = par->hardeningLinearConstant;
  double dt = par->timeStep;

  double* rateOfDefXX;
  double* rateOfDefXY;
  double* rateOfDefXZ;
  double* rateOfDefYX;
  double* rateOfDefYY;
  double* rateOfDefYZ;
  double* rateOfDefZX;
  double* rateOfDefZY;
  double* rateOfDefZZ;
  double* stressXXN;
  double* stressXYN;
  double* stressXZN;
  double* stressYXN;
  double* stressYYN;
  double* stressYZN;
  double* stressZXN;
  double* stressZYN;
  double* stressZZN;
  double* stressXXNP1;
  double* stressXYNP1;
  double* stressXZNP1;
  double* stressYXNP1;
  double* stressYYNP1;
  double* stressYZNP1;
  double* stressZXNP1;
  double* stressZYNP1;
  double* stressZZNP1;
  
  double* eqpsN;
  double* eqpsNP1;

  double* vmStress;
  double* minPrincipalStress;
  double* maxPrincipalStress;

  double C33;

  double rateOfDef[9];
  double stressN[9];
  double stressNP1[9];

  double tol = 1.0e-16 * bulkMod;
  int counter;
  int maxIterations = 100;

  double V[3][3], d[3], e[3];

  for(in=0;in<par->numPoints;in++){

    point = &par->puntos[in];

    for(l=0; l<par->numLayers; l++){

      rateOfDefXX = point->layerAssociatedBondLevelUnrotatedRateOfDeformation[l][0];
      rateOfDefXY = point->layerAssociatedBondLevelUnrotatedRateOfDeformation[l][1];
      rateOfDefXZ = point->layerAssociatedBondLevelUnrotatedRateOfDeformation[l][2];
      rateOfDefYX = point->layerAssociatedBondLevelUnrotatedRateOfDeformation[l][3];
      rateOfDefYY = point->layerAssociatedBondLevelUnrotatedRateOfDeformation[l][4];
      rateOfDefYZ = point->layerAssociatedBondLevelUnrotatedRateOfDeformation[l][5];
      rateOfDefZX = point->layerAssociatedBondLevelUnrotatedRateOfDeformation[l][6];
      rateOfDefZY = point->layerAssociatedBondLevelUnrotatedRateOfDeformation[l][7];
      rateOfDefZZ = point->layerAssociatedBondLevelUnrotatedRateOfDeformation[l][8];

      stressXXN = point->layerAssociatedBondLevelUnrotatedCauchyStress[0][l][0];
      stressXYN = point->layerAssociatedBondLevelUnrotatedCauchyStress[0][l][1];
      stressXZN = point->layerAssociatedBondLevelUnrotatedCauchyStress[0][l][2];
      stressYXN = point->layerAssociatedBondLevelUnrotatedCauchyStress[0][l][3];
      stressYYN = point->layerAssociatedBondLevelUnrotatedCauchyStress[0][l][4];
      stressYZN = point->layerAssociatedBondLevelUnrotatedCauchyStress[0][l][5];
      stressZXN = point->layerAssociatedBondLevelUnrotatedCauchyStress[0][l][6];
      stressZYN = point->layerAssociatedBondLevelUnrotatedCauchyStress[0][l][7];
      stressZZN = point->layerAssociatedBondLevelUnrotatedCauchyStress[0][l][8];

      stressXXNP1 = point->layerAssociatedBondLevelUnrotatedCauchyStress[1][l][0];
      stressXYNP1 = point->layerAssociatedBondLevelUnrotatedCauchyStress[1][l][1];
      stressXZNP1 = point->layerAssociatedBondLevelUnrotatedCauchyStress[1][l][2];
      stressYXNP1 = point->layerAssociatedBondLevelUnrotatedCauchyStress[1][l][3];
      stressYYNP1 = point->layerAssociatedBondLevelUnrotatedCauchyStress[1][l][4];
      stressYZNP1 = point->layerAssociatedBondLevelUnrotatedCauchyStress[1][l][5];
      stressZXNP1 = point->layerAssociatedBondLevelUnrotatedCauchyStress[1][l][6];
      stressZYNP1 = point->layerAssociatedBondLevelUnrotatedCauchyStress[1][l][7];
      stressZZNP1 = point->layerAssociatedBondLevelUnrotatedCauchyStress[1][l][8];

      eqpsN = point->layerAssociatedBondLevelEquivalentPlasticStrain[0][l];
      eqpsNP1 = point->layerAssociatedBondLevelEquivalentPlasticStrain[1][l];
      vmStress = point->layerAssociatedBondLevelVonMisesStress[l];
      minPrincipalStress = point->layerAssociatedBondLevelMinPrincipalStress[l];
      maxPrincipalStress = point->layerAssociatedBondLevelMaxPrincipalStress[l];

      numNeighbors = point->numNeighbors;
      for(n=0; n<numNeighbors; n++,
          rateOfDefXX++, rateOfDefXY++, rateOfDefXZ++, 
          rateOfDefYX++, rateOfDefYY++, rateOfDefYZ++, 
          rateOfDefZX++, rateOfDefZY++, rateOfDefZZ++,
          stressXXN++, stressXYN++, stressXZN++, 
          stressYXN++, stressYYN++, stressYZN++, 
          stressZXN++, stressZYN++, stressZZN++,
          stressXXNP1++, stressXYNP1++, stressXZNP1++, 
          stressYXNP1++, stressYYNP1++, stressYZNP1++, 
          stressZXNP1++, stressZYNP1++, stressZZNP1++,
          eqpsN++, eqpsNP1++, vmStress++, 
          minPrincipalStress++, maxPrincipalStress++){

        *(rateOfDef+0) = *rateOfDefXX; *(rateOfDef+1) = *rateOfDefXY; *(rateOfDef+2) = *rateOfDefXZ; 
        *(rateOfDef+3) = *rateOfDefYX; *(rateOfDef+4) = *rateOfDefYY; *(rateOfDef+5) = *rateOfDefYZ; 
        *(rateOfDef+6) = *rateOfDefZX; *(rateOfDef+7) = *rateOfDefZY; *(rateOfDef+8) = *rateOfDefZZ; 
        *(stressN+0) = *stressXXN; *(stressN+1) = *stressXYN; *(stressN+2) = *stressXZN; 
        *(stressN+3) = *stressYXN; *(stressN+4) = *stressYYN; *(stressN+5) = *stressYZN; 
        *(stressN+6) = *stressZXN; *(stressN+7) = *stressZYN; *(stressN+8) = *stressZZN; 

        radialReturnSaturationExponential(rateOfDef, 
                                          stressN, 
                                          stressNP1, 
                                          vmStress,
                                          eqpsN, 
                                          eqpsNP1, 
                                          &C33,
                                          bulkMod, 
                                          shearMod,
                                          initialYieldStress, 
                                          saturatedYieldStress,
                                          exponentialConstant, 
                                          linearConstant, 
                                          dt);
      

        // Enforce the zero-through-thickness stress 
        // A D33 is added to the system such that make stress33 zero
        counter = 0;
        while(fabs(*(stressNP1+8)) > tol){

          if(counter > maxIterations){
            std::cout << "updateBondLevelElasticIsotropicSaturationExponentialHardeningPlasticCauchyStress: Too many iterations to enforce zero through-thickness-stress." << std::endl;
            exit(1);
          }

          *(rateOfDef+8) -= *(stressNP1+8) / (C33*dt);

          radialReturnSaturationExponential(rateOfDef, 
                                            stressN, 
                                            stressNP1, 
                                            vmStress,
                                            eqpsN, 
                                            eqpsNP1, 
                                            &C33,
                                            bulkMod, 
                                            shearMod,
                                            initialYieldStress, 
                                            saturatedYieldStress,
                                            exponentialConstant, 
                                            linearConstant, 
                                            dt);
          counter++;
        }

        // compute eigenvalues
        for (i = 0; i < 3; i++) {
          for (j = 0; j < 3; j++) {
            V[i][j] = *(stressNP1+3*i+j);
          }
        }
        tred2(V, d, e); // diagonalize V

        // find the max and min of eigenvalues
        *minPrincipalStress = d[0];
        *maxPrincipalStress = d[0];
        for(i=1; i<3; i++){
          if(d[i] > *maxPrincipalStress)
            *maxPrincipalStress = d[i];
          if(d[i] < *minPrincipalStress)
            *minPrincipalStress = d[i];
        }

        *rateOfDefZZ = *(rateOfDef+8);

        *stressXXNP1 = *(stressNP1+0); *stressXYNP1 = *(stressNP1+1); *stressXZNP1 = *(stressNP1+2); 
        *stressYXNP1 = *(stressNP1+3); *stressYYNP1 = *(stressNP1+4); *stressYZNP1 = *(stressNP1+5); 
        *stressZXNP1 = *(stressNP1+6); *stressZYNP1 = *(stressNP1+7); *stressZZNP1 = *(stressNP1+8); 
      }
    }
  }
}


#undef  __FUNCT__
#define __FUNCT__ "rotateCauchyStress"
void rotateCauchyStress(PARAMETERS *par) 
{
  //Rotate the unrotated Cauchy stress to the current configuration using
  //the rotation tensor
  //This function computes the node-level values (only for visualization).

  POINTS *point;
  int in;

  double* rotTensor;
  double* unrotatedStress;
  double* rotatedStress;

  double temp[9];

  for(in=0;in<par->numPoints;in++){
    point = &par->puntos[in];

    unrotatedStress = point->unrotatedCauchyStress[1];
    rotTensor = point->rotationTensor[1];
    rotatedStress = point->cauchyStress; 

    // temp = \sigma_unrot * Rt
    MatrixMultiply(false, true, 1.0, unrotatedStress, rotTensor, 3, 3, 3, temp);

    // \sigma_rot = R * temp
    MatrixMultiply(false, false, 1.0, rotTensor, temp, 3, 3, 3, rotatedStress);
  }
}


#undef  __FUNCT__
#define __FUNCT__ "rotateBondLevelCauchyStress"
void rotateBondLevelCauchyStress(PARAMETERS *par) 
{
  POINTS *point;
  int in,l,n;
  int numNeighbors;

  double* rotTensorXX;
  double* rotTensorXY;
  double* rotTensorXZ;
  double* rotTensorYX;
  double* rotTensorYY;
  double* rotTensorYZ;
  double* rotTensorZX;
  double* rotTensorZY;
  double* rotTensorZZ;
  double* unrotatedStressXX;
  double* unrotatedStressXY;
  double* unrotatedStressXZ;
  double* unrotatedStressYX;
  double* unrotatedStressYY;
  double* unrotatedStressYZ;
  double* unrotatedStressZX;
  double* unrotatedStressZY;
  double* unrotatedStressZZ;
  double* rotatedStressXX;
  double* rotatedStressXY;
  double* rotatedStressXZ;
  double* rotatedStressYX;
  double* rotatedStressYY;
  double* rotatedStressYZ;
  double* rotatedStressZX;
  double* rotatedStressZY;
  double* rotatedStressZZ;

  double unrotatedStress[9], rotTensor[9], rotatedStress[9], temp[9];


  for(in=0;in<par->numPoints;in++){

    point = &par->puntos[in];

    for(l=0; l<par->numLayers; l++){

      rotTensorXX = point->layerAssociatedBondLevelRotationTensor[1][l][0];
      rotTensorXY = point->layerAssociatedBondLevelRotationTensor[1][l][1];
      rotTensorXZ = point->layerAssociatedBondLevelRotationTensor[1][l][2];
      rotTensorYX = point->layerAssociatedBondLevelRotationTensor[1][l][3];
      rotTensorYY = point->layerAssociatedBondLevelRotationTensor[1][l][4];
      rotTensorYZ = point->layerAssociatedBondLevelRotationTensor[1][l][5];
      rotTensorZX = point->layerAssociatedBondLevelRotationTensor[1][l][6];
      rotTensorZY = point->layerAssociatedBondLevelRotationTensor[1][l][7];
      rotTensorZZ = point->layerAssociatedBondLevelRotationTensor[1][l][8];

      unrotatedStressXX = point->layerAssociatedBondLevelUnrotatedCauchyStress[1][l][0];
      unrotatedStressXY = point->layerAssociatedBondLevelUnrotatedCauchyStress[1][l][1];
      unrotatedStressXZ = point->layerAssociatedBondLevelUnrotatedCauchyStress[1][l][2];
      unrotatedStressYX = point->layerAssociatedBondLevelUnrotatedCauchyStress[1][l][3];
      unrotatedStressYY = point->layerAssociatedBondLevelUnrotatedCauchyStress[1][l][4];
      unrotatedStressYZ = point->layerAssociatedBondLevelUnrotatedCauchyStress[1][l][5];
      unrotatedStressZX = point->layerAssociatedBondLevelUnrotatedCauchyStress[1][l][6];
      unrotatedStressZY = point->layerAssociatedBondLevelUnrotatedCauchyStress[1][l][7];
      unrotatedStressZZ = point->layerAssociatedBondLevelUnrotatedCauchyStress[1][l][8];

      rotatedStressXX = point->layerAssociatedBondLevelCauchyStress[l][0];
      rotatedStressXY = point->layerAssociatedBondLevelCauchyStress[l][1];
      rotatedStressXZ = point->layerAssociatedBondLevelCauchyStress[l][2];
      rotatedStressYX = point->layerAssociatedBondLevelCauchyStress[l][3];
      rotatedStressYY = point->layerAssociatedBondLevelCauchyStress[l][4];
      rotatedStressYZ = point->layerAssociatedBondLevelCauchyStress[l][5];
      rotatedStressZX = point->layerAssociatedBondLevelCauchyStress[l][6];
      rotatedStressZY = point->layerAssociatedBondLevelCauchyStress[l][7];
      rotatedStressZZ = point->layerAssociatedBondLevelCauchyStress[l][8];


      numNeighbors = point->numNeighbors;
      for(n=0; n<numNeighbors; n++,
          rotTensorXX++, rotTensorXY++, rotTensorXZ++, 
          rotTensorYX++, rotTensorYY++, rotTensorYZ++, 
          rotTensorZX++, rotTensorZY++, rotTensorZZ++, 
          unrotatedStressXX++, unrotatedStressXY++, unrotatedStressXZ++, 
          unrotatedStressYX++, unrotatedStressYY++, unrotatedStressYZ++, 
          unrotatedStressZX++, unrotatedStressZY++, unrotatedStressZZ++, 
          rotatedStressXX++, rotatedStressXY++, rotatedStressXZ++, 
          rotatedStressYX++, rotatedStressYY++, rotatedStressYZ++, 
          rotatedStressZX++, rotatedStressZY++, rotatedStressZZ++){

        // write in matrix form 
        rotTensor[0] = *rotTensorXX; rotTensor[1] = *rotTensorXY; rotTensor[2] = *rotTensorXZ;
        rotTensor[3] = *rotTensorYX; rotTensor[4] = *rotTensorYY; rotTensor[5] = *rotTensorYZ;
        rotTensor[6] = *rotTensorZX; rotTensor[7] = *rotTensorZY; rotTensor[8] = *rotTensorZZ;
        unrotatedStress[0] = *unrotatedStressXX; unrotatedStress[1] = *unrotatedStressXY; unrotatedStress[2] = *unrotatedStressXZ;
        unrotatedStress[3] = *unrotatedStressYX; unrotatedStress[4] = *unrotatedStressYY; unrotatedStress[5] = *unrotatedStressYZ;
        unrotatedStress[6] = *unrotatedStressZX; unrotatedStress[7] = *unrotatedStressZY; unrotatedStress[8] = *unrotatedStressZZ;

        // temp = \sigma_unrot * Rt
        MatrixMultiply(false, true, 1.0, unrotatedStress, rotTensor, 3, 3, 3, temp);
        // \sigma_rot = R * temp
        MatrixMultiply(false, false, 1.0, rotTensor, temp, 3, 3, 3, rotatedStress);

        // update bond-level data field
        *rotatedStressXX = rotatedStress[0]; *rotatedStressXY = rotatedStress[1]; *rotatedStressXZ = rotatedStress[2]; 
        *rotatedStressYX = rotatedStress[3]; *rotatedStressYY = rotatedStress[4]; *rotatedStressYZ = rotatedStress[5]; 
        *rotatedStressZX = rotatedStress[6]; *rotatedStressZY = rotatedStress[7]; *rotatedStressZZ = rotatedStress[8]; 
      }
    }
  }
}


#undef  __FUNCT__
#define __FUNCT__ "updateBondLevelJacobianDeterminant"
void updateBondLevelJacobianDeterminant(PARAMETERS *par) 
{
  POINTS *point;
  int in,l,n;
  int numNeighbors;

  double dt = par->timeStep;

  double* bondLevelJN;
  double* bondLevelJNP1;
  double* unrotRateOfDefXX;
  double* unrotRateOfDefYY;
  double* unrotRateOfDefZZ;
  double velGradTrace;

  for(in=0;in<par->numPoints;in++){

    point = &par->puntos[in];

    for(l=0; l<par->numLayers; l++){

      bondLevelJN = point->layerAssociatedBondLevelJacobianDeterminant[0][l];
      bondLevelJNP1 = point->layerAssociatedBondLevelJacobianDeterminant[1][l];
      unrotRateOfDefXX = point->layerAssociatedBondLevelUnrotatedRateOfDeformation[l][0];
      unrotRateOfDefYY = point->layerAssociatedBondLevelUnrotatedRateOfDeformation[l][4];
      unrotRateOfDefZZ = point->layerAssociatedBondLevelUnrotatedRateOfDeformation[l][8];

      numNeighbors = point->numNeighbors;
      for(n=0; n<numNeighbors; n++,
          bondLevelJN++, bondLevelJNP1++, 
          unrotRateOfDefXX++, unrotRateOfDefYY++, unrotRateOfDefZZ++){

        velGradTrace = *unrotRateOfDefXX + *unrotRateOfDefYY + *unrotRateOfDefZZ;
        *bondLevelJNP1 = *bondLevelJN * (1.0 + velGradTrace*dt);
      }
    }
  }
}


#undef  __FUNCT__
#define __FUNCT__ "computeBondLevelKirchhoffStress"
void computeBondLevelKirchhoffStress(PARAMETERS *par) 
{
  POINTS *point;
  int in,l,n;
  int numNeighbors;

  double* J;
  double* cauchyStressXX;
  double* cauchyStressXY;
  double* cauchyStressXZ;
  double* cauchyStressYX;
  double* cauchyStressYY;
  double* cauchyStressYZ;
  double* cauchyStressZX;
  double* cauchyStressZY;
  double* cauchyStressZZ;
  double* kirchhoffStressXX;
  double* kirchhoffStressXY;
  double* kirchhoffStressXZ;
  double* kirchhoffStressYX;
  double* kirchhoffStressYY;
  double* kirchhoffStressYZ;
  double* kirchhoffStressZX;
  double* kirchhoffStressZY;
  double* kirchhoffStressZZ;

  for(in=0;in<par->numPoints;in++){

    point = &par->puntos[in];

    for(l=0; l<par->numLayers; l++){

      J = point->layerAssociatedBondLevelJacobianDeterminant[1][l];
      cauchyStressXX = point->layerAssociatedBondLevelCauchyStress[l][0];
      cauchyStressXY = point->layerAssociatedBondLevelCauchyStress[l][1];
      cauchyStressXZ = point->layerAssociatedBondLevelCauchyStress[l][2];
      cauchyStressYX = point->layerAssociatedBondLevelCauchyStress[l][3];
      cauchyStressYY = point->layerAssociatedBondLevelCauchyStress[l][4];
      cauchyStressYZ = point->layerAssociatedBondLevelCauchyStress[l][5];
      cauchyStressZX = point->layerAssociatedBondLevelCauchyStress[l][6];
      cauchyStressZY = point->layerAssociatedBondLevelCauchyStress[l][7];
      cauchyStressZZ = point->layerAssociatedBondLevelCauchyStress[l][8];
      kirchhoffStressXX = point->layerAssociatedBondLevelKirchhoffStress[l][0];
      kirchhoffStressXY = point->layerAssociatedBondLevelKirchhoffStress[l][1];
      kirchhoffStressXZ = point->layerAssociatedBondLevelKirchhoffStress[l][2];
      kirchhoffStressYX = point->layerAssociatedBondLevelKirchhoffStress[l][3];
      kirchhoffStressYY = point->layerAssociatedBondLevelKirchhoffStress[l][4];
      kirchhoffStressYZ = point->layerAssociatedBondLevelKirchhoffStress[l][5];
      kirchhoffStressZX = point->layerAssociatedBondLevelKirchhoffStress[l][6];
      kirchhoffStressZY = point->layerAssociatedBondLevelKirchhoffStress[l][7];
      kirchhoffStressZZ = point->layerAssociatedBondLevelKirchhoffStress[l][8];

      numNeighbors = point->numNeighbors;
      for(n=0; n<numNeighbors; n++,
          J++,
          cauchyStressXX++, cauchyStressXY++, cauchyStressXZ++, 
          cauchyStressYX++, cauchyStressYY++, cauchyStressYZ++, 
          cauchyStressZX++, cauchyStressZY++, cauchyStressZZ++, 
          kirchhoffStressXX++, kirchhoffStressXY++, kirchhoffStressXZ++, 
          kirchhoffStressYX++, kirchhoffStressYY++, kirchhoffStressYZ++, 
          kirchhoffStressZX++, kirchhoffStressZY++, kirchhoffStressZZ++){

        *kirchhoffStressXX = *J * *cauchyStressXX;
        *kirchhoffStressXY = *J * *cauchyStressXY;
        *kirchhoffStressXZ = *J * *cauchyStressXZ;
        *kirchhoffStressYX = *J * *cauchyStressYX;
        *kirchhoffStressYY = *J * *cauchyStressYY;
        *kirchhoffStressYZ = *J * *cauchyStressYZ;
        *kirchhoffStressZX = *J * *cauchyStressZX;
        *kirchhoffStressZY = *J * *cauchyStressZY;
        *kirchhoffStressZZ = *J * *cauchyStressZZ;
      }
    }
  }
}


#undef  __FUNCT__
#define __FUNCT__ "updateNodeLevelStress"
void updateNodeLevelStress(PARAMETERS *par) 
{
  // A set of functions to be called for updating the stress tensor 
  // at the material point level (only for visualization)

  // The stress update takes place in the unrotated configuration to maintain
  // objectivity - Following Flanagan and Taylor (1987)
  
  computeUnrotatedRateOfDeformationAndRotationTensor(par); // Compute the unrotated rate of deformation tensor and update the rotation tensor


  // Evaluate the Cauchy stress using the implemented routine
  // The general idea is to compute the stress based on:
  //   1) The unrotated rate-of-deformation tensor
  //   2) The time step
  //   3) Whatever state variables are managed by the derived class
  //
  //   Input:  unrotated rate-of-deformation tensor
  //   Input:  unrotated Cauchy stress at step N
  //   Input:  internal state data (managed in the derived class)
  //   Output: unrotated Cauchy stress at step N+1
  // 
  // Currently, an isotropic hardening material law is included in the code
  updateIsotropicHardeningCauchyStress(par); // Update the unrotated stress using the unrotated rate of deformation

  rotateCauchyStress(par); // Rotate the unrotated stress to the current configuration
}


#undef  __FUNCT__
#define __FUNCT__ "updateStress"
void updateBondLevelStress(PARAMETERS *par) 
{
  // A set of functions to be called for updating the stress tensor 
  // at the bond level 

  // The stress update takes place in the unrotated configuration to maintain
  // objectivity - Following Flanagan and Taylor (1987)

  computeBondLevelUnrotatedRateOfDeformationAndRotationTensor(par); // Compute the unrotated rate of deformation tensor and update the rotation tensor
  
  // Evaluate the Cauchy stress using the implemented routine
  // The general idea is to compute the stress based on:
  //   1) The unrotated rate-of-deformation tensor
  //   2) The time step
  //   3) Whatever state variables are managed by the derived class
  //
  //   Input:  unrotated rate-of-deformation tensor
  //   Input:  unrotated Cauchy stress at step N
  //   Input:  internal state data (managed in the derived class)
  //   Output: unrotated Cauchy stress at step N+1
  // 
  // Currently, an isotropic hardening material law is included in the code
  updateBondLevelIsotropicHardeningCauchyStress(par); // Update the unrotated stress using the unrotated rate of deformation

  rotateBondLevelCauchyStress(par); // Rotate the unrotated stress to the current configuration

  updateBondLevelJacobianDeterminant(par); // Update the Jacobian determinant

  computeBondLevelKirchhoffStress(par); // Compute the Kirchhoff stress tensor
}


#undef  __FUNCT__
#define __FUNCT__ "computeAStateIntegral"
void computeAStateIntegral(PARAMETERS *par) 
{
  POINTS *point;
  POINTS *neighbor;
  int in,i,l,n;
  int numNeighbors, neigh;
  int *neighbors;

  double* neighborArea;
  double* A0;
  double* neighborA0;
  double* refThickness;
  double* neighborRefThickness;
  double* omega;
  double* bondDefX;
  double* bondDefY;
  double* bondDefZ;
  double* stressXX;
  double* stressXY;
  double* stressXZ;
  double* stressYX;
  double* stressYY;
  double* stressYZ;
  double* stressZX;
  double* stressZY;
  double* stressZZ;
  double* integral;

  double deformedBondX, deformedBondY, deformedBondZ, deformedBondLengthSq;
  double temp;

  double aState[3];

  for(in=0;in<par->numPoints;in++){

    point = &par->puntos[in];
    neighbors = point->neighbors;
    A0 = &(point->weightedArea);
    refThickness = &(point->referenceShellThickness);

    for(l=0; l<par->numLayers; l++){

      integral = point->layerAssociatedAStateIntegral[l];

      // Zero out data
      for(i=0; i<3; i++)
        *(integral+i) = 0.0;

      omega = point->influenceState;
      bondDefX = point->layerAssociatedBondDeformedState[l][0];
      bondDefY = point->layerAssociatedBondDeformedState[l][1];
      bondDefZ = point->layerAssociatedBondDeformedState[l][2];
      stressXX = point->layerAssociatedBondLevelKirchhoffStress[l][0];
      stressXY = point->layerAssociatedBondLevelKirchhoffStress[l][1];
      stressXZ = point->layerAssociatedBondLevelKirchhoffStress[l][2];
      stressYX = point->layerAssociatedBondLevelKirchhoffStress[l][3];
      stressYY = point->layerAssociatedBondLevelKirchhoffStress[l][4];
      stressYZ = point->layerAssociatedBondLevelKirchhoffStress[l][5];
      stressZX = point->layerAssociatedBondLevelKirchhoffStress[l][6];
      stressZY = point->layerAssociatedBondLevelKirchhoffStress[l][7];
      stressZZ = point->layerAssociatedBondLevelKirchhoffStress[l][8];

      numNeighbors = point->numNeighbors;
      for(n=0; n<numNeighbors; n++, omega++,
        bondDefX++, bondDefY++, bondDefZ++,
        stressXX++, stressXY++, stressXZ++,
        stressYX++, stressYY++, stressYZ++,
        stressZX++, stressZY++, stressZZ++){

        // only do the calculations for the unbroken bonds
        if(*omega > 0.0){

          neigh = neighbors[n];
          neighbor = &par->puntos[neigh];
          neighborArea = &(neighbor->area);
          neighborA0 = &(neighbor->weightedArea);
          neighborRefThickness = &(neighbor->referenceShellThickness);

          deformedBondX = *bondDefX;
          deformedBondY = *bondDefY;
          deformedBondZ = *bondDefZ;
          deformedBondLengthSq = deformedBondX*deformedBondX +
                                 deformedBondY*deformedBondY +
                                 deformedBondZ*deformedBondZ;

          // The negative sign is because it actually is the negative of the integral of aState
          temp = - *omega/2.0 * ( *refThickness / *A0 + *neighborRefThickness / *neighborA0 ) / deformedBondLengthSq;

          *(aState+0) = temp * ( *stressXX * deformedBondX + 
                                 *stressXY * deformedBondY + 
                                 *stressXZ * deformedBondZ ); 

          *(aState+1) = temp * ( *stressYX * deformedBondX + 
                                 *stressYY * deformedBondY + 
                                 *stressYZ * deformedBondZ ); 

          *(aState+2) = temp * ( *stressZX * deformedBondX + 
                                 *stressZY * deformedBondY + 
                                 *stressZZ * deformedBondZ ); 

          for(i=0; i<3; ++i)
            *(integral+i) += *(aState+i) * *neighborArea;
        }
      }
    }
  }
}


#undef  __FUNCT__
#define __FUNCT__ "computeBStateIntegral"
void computeBStateIntegral(PARAMETERS *par) 
{
  POINTS *point;
  POINTS *neighbor;
  int in,i,l,n;
  int numNeighbors, neigh;
  int *neighbors;

  double* neighborArea;
  double* A0;
  double* neighborA0;
  double* refThickness;
  double* neighborRefThickness;
  double* omega;
  double* bondDefX;
  double* bondDefY;
  double* bondDefZ;
  double* stressXX;
  double* stressXY;
  double* stressXZ;
  double* stressYX;
  double* stressYY;
  double* stressYZ;
  double* stressZX;
  double* stressZY;
  double* stressZZ;
  double* defGradInv;
  double* integral;

  double deformedBondX, deformedBondY, deformedBondZ, deformedBondLengthSq;
  double scalarTemp;

  double temp[9];
  double stress[9];
  double bState[9];
  double bInt[9];

  for(in=0;in<par->numPoints;in++){

    point = &par->puntos[in];
    neighbors = point->neighbors;
    A0 = &(point->weightedArea);
    refThickness = &(point->referenceShellThickness);

    for(l=0; l<par->numLayers; l++){

      integral = point->layerAssociatedBStateIntegral[l];
      defGradInv = point->layerAssociatedDeformationGradientInverse[l];

      // Zero out data
      for(i=0; i<9; i++)
        *(bInt+i) = 0.0;

      omega = point->influenceState;
      bondDefX = point->layerAssociatedBondDeformedState[l][0];
      bondDefY = point->layerAssociatedBondDeformedState[l][1];
      bondDefZ = point->layerAssociatedBondDeformedState[l][2];
      stressXX = point->layerAssociatedBondLevelKirchhoffStress[l][0];
      stressXY = point->layerAssociatedBondLevelKirchhoffStress[l][1];
      stressXZ = point->layerAssociatedBondLevelKirchhoffStress[l][2];
      stressYX = point->layerAssociatedBondLevelKirchhoffStress[l][3];
      stressYY = point->layerAssociatedBondLevelKirchhoffStress[l][4];
      stressYZ = point->layerAssociatedBondLevelKirchhoffStress[l][5];
      stressZX = point->layerAssociatedBondLevelKirchhoffStress[l][6];
      stressZY = point->layerAssociatedBondLevelKirchhoffStress[l][7];
      stressZZ = point->layerAssociatedBondLevelKirchhoffStress[l][8];

      numNeighbors = point->numNeighbors;
      for(n=0; n<numNeighbors; n++, omega++,
        bondDefX++, bondDefY++, bondDefZ++,
        stressXX++, stressXY++, stressXZ++,
        stressYX++, stressYY++, stressYZ++,
        stressZX++, stressZY++, stressZZ++){

        // only do the calculations for the unbroken bonds
        if(*omega > 0.0){

          neigh = neighbors[n];
          neighbor = &par->puntos[neigh];
          neighborArea = &(neighbor->area);
          neighborA0 = &(neighbor->weightedArea);
          neighborRefThickness = &(neighbor->referenceShellThickness);

          deformedBondX = *bondDefX;
          deformedBondY = *bondDefY;
          deformedBondZ = *bondDefZ;
          deformedBondLengthSq = deformedBondX*deformedBondX +
                                 deformedBondY*deformedBondY +
                                 deformedBondZ*deformedBondZ;

          scalarTemp = *omega/4.0 * (*refThickness / *A0 + *neighborRefThickness / *neighborA0);

          // write the stress in matrix form 
          stress[0] = *stressXX; stress[1] = *stressXY; stress[2] = *stressXZ; 
          stress[3] = *stressYX; stress[4] = *stressYY; stress[5] = *stressYZ; 
          stress[6] = *stressZX; stress[7] = *stressZY; stress[8] = *stressZZ; 

          // delta_jp - (y_j y_p)/|y|^2
          *(temp+0) = 1.0 - deformedBondX * deformedBondX / deformedBondLengthSq;
          *(temp+1) = - deformedBondX * deformedBondY / deformedBondLengthSq;
          *(temp+2) = - deformedBondX * deformedBondZ / deformedBondLengthSq;
          *(temp+3) = *(temp+1);
          *(temp+4) = 1.0 - deformedBondY * deformedBondY / deformedBondLengthSq;
          *(temp+5) = - deformedBondY * deformedBondZ / deformedBondLengthSq;
          *(temp+6) = *(temp+2);
          *(temp+7) = *(temp+5);
          *(temp+8) = 1.0 - deformedBondZ * deformedBondZ / deformedBondLengthSq;

          // Matrix multiply the stress and the second term to compute b
          MatrixMultiply(false, false, scalarTemp, stress, temp, 3, 3, 3, bState);

          for(i=0 ; i<9 ; ++i)
            *(bInt+i) += *(bState+i) * *neighborArea;
        }
      }

      // b_int * F_inv^T
      MatrixMultiply(false, true, 1.0, bInt, defGradInv, 3, 3, 3, integral);
    }
  }
}


#undef  __FUNCT__
#define __FUNCT__ "computeInternalForceState"
void computeInternalForceState(PARAMETERS *par) 
{
  // After computing the bond-level stress and damage variables, here we 
  // compute the internal force state in each layer and eventually integrate 
  // over the thickness to obtain the force density per unit area

  computeAStateIntegral(par); // Compute the discrete sum of a_state

  computeBStateIntegral(par); // Compute the discrete sum of b_state F_inv

  POINTS *point;
  POINTS *neighbor;
  int in,i,j,l,n;
  int numNeighbors, neigh;
  int *neighbors;

  // gauss quadrature points along the thickness
  double* xi3list = par->xi3list;
  double* w3list = par->w3list;
  double xi3, w3;

  double* A0;
  double* refThickness;
  double* thickness;
  double* B1;
  double* dB1xi1;
  double* dB1xi2;
  double* B2;
  double* dB2xi1;
  double* dB2xi2;

  double* aBar;
  double* bBar;

  double *omega;
  double* phi1;
  double* phi2;
  double* phi11;
  double* phi12;
  double* phi22;

  double* bondDefX;
  double* bondDefY;
  double* bondDefZ;
  double* stressXX;
  double* stressXY;
  double* stressXZ;
  double* stressYX;
  double* stressYY;
  double* stressYZ;
  double* stressZX;
  double* stressZY;
  double* stressZZ;

  double deformedBondX, deformedBondY, deformedBondZ, deformedBondLengthSq;

  double *forceDensity, *neighborForceDensity;
  double TX, TY, TZ, nodalArea, neighborArea;
  double alpha, temp;

  double aState[3];

  double gamma[9];

  double betaX[9];
  double betaY[9];
  double betaZ[9];

  // Zero out the force state 
  for(in=0;in<par->numPoints;in++){
    for(i=0; i<3; i++){
      par->puntos[in].internalForceDensity[i] = 0.0;
    }
  }

  // Loop over the material points and convert the Cauchy stress into pairwise peridynamic force densities
  for(in=0;in<par->numPoints;in++){

    point = &par->puntos[in];
    neighbors = point->neighbors;
    A0 = &(point->weightedArea);
    refThickness = &(point->referenceShellThickness);
    thickness = &(point->shellThickness[1]);
    B1 = point->B1Tensor;
    dB1xi1 = point->B1TensorGradient[0];
    dB1xi2 = point->B1TensorGradient[1];
    B2 = point->B2Tensor;
    dB2xi1 = point->B2TensorGradient[0];
    dB2xi2 = point->B2TensorGradient[1];
    forceDensity = point->internalForceDensity;
    nodalArea = point->area;

    for(l=0; l<par->numLayers; l++){

      xi3 = xi3list[l];
      w3 = w3list[l];

      aBar = point->layerAssociatedAStateIntegral[l];
      bBar = point->layerAssociatedBStateIntegral[l];

      omega = point->influenceState;
      phi1 = point->gradientWeight1[0];
      phi2 = point->gradientWeight1[1];
      phi11 = point->gradientWeight2[0];
      phi12 = point->gradientWeight2[1];
      phi22 = point->gradientWeight2[2];
      bondDefX = point->layerAssociatedBondDeformedState[l][0];
      bondDefY = point->layerAssociatedBondDeformedState[l][1];
      bondDefZ = point->layerAssociatedBondDeformedState[l][2];
      stressXX = point->layerAssociatedBondLevelKirchhoffStress[l][0];
      stressXY = point->layerAssociatedBondLevelKirchhoffStress[l][1];
      stressXZ = point->layerAssociatedBondLevelKirchhoffStress[l][2];
      stressYX = point->layerAssociatedBondLevelKirchhoffStress[l][3];
      stressYY = point->layerAssociatedBondLevelKirchhoffStress[l][4];
      stressYZ = point->layerAssociatedBondLevelKirchhoffStress[l][5];
      stressZX = point->layerAssociatedBondLevelKirchhoffStress[l][6];
      stressZY = point->layerAssociatedBondLevelKirchhoffStress[l][7];
      stressZZ = point->layerAssociatedBondLevelKirchhoffStress[l][8];

      // Loop over the neighbors and compute contribution to force densities
      numNeighbors = point->numNeighbors;
      for(n=0; n<numNeighbors; n++, omega++,
              phi1++, phi2++, phi11++, phi12++, phi22++,
              bondDefX++, bondDefY++, bondDefZ++,
              stressXX++, stressXY++, stressXZ++, 
              stressYX++, stressYY++, stressYZ++, 
              stressZX++, stressZY++, stressZZ++){

        // only do the calculations for the unbroken bonds
        if(*omega > 0.0){

          neigh = neighbors[n];
          neighbor = &par->puntos[neigh];
          neighborArea = neighbor->area;
          neighborForceDensity = neighbor->internalForceDensity;

          deformedBondX = *bondDefX;
          deformedBondY = *bondDefY;
          deformedBondZ = *bondDefZ;
          deformedBondLengthSq = deformedBondX*deformedBondX +
                                 deformedBondY*deformedBondY +
                                 deformedBondZ*deformedBondZ;

          // Compute a state
          alpha = *omega / *A0;
          temp = alpha/2.0 * *refThickness / deformedBondLengthSq;

          *(aState+0) = temp * ( *stressXX * deformedBondX + 
                                 *stressXY * deformedBondY + 
                                 *stressXZ * deformedBondZ ); 

          *(aState+1) = temp * ( *stressYX * deformedBondX + 
                                 *stressYY * deformedBondY + 
                                 *stressYZ * deformedBondZ ); 

          *(aState+2) = temp * ( *stressZX * deformedBondX + 
                                 *stressZY * deformedBondY + 
                                 *stressZZ * deformedBondZ ); 

          // Compute gamma
          for(i=0; i<9; i++)
            *(gamma+i) = *thickness/2.0 * xi3 * ( *(B1+i) * *phi1 + *(B2+i) * *phi2 );

          // Compute beta
          for(i=0; i<3; i++){
            *(betaX+3*i+0) = *thickness/2.0 * xi3 * ( *(B1+3*i+0) * *phi11 + *(B2+3*i+0) * *phi12 + *(dB1xi1+3*i+0) * *phi1 + *(dB2xi1+3*i+0) * *phi2 );
            *(betaY+3*i+0) = *thickness/2.0 * xi3 * ( *(B1+3*i+1) * *phi11 + *(B2+3*i+1) * *phi12 + *(dB1xi1+3*i+1) * *phi1 + *(dB2xi1+3*i+1) * *phi2 );
            *(betaZ+3*i+0) = *thickness/2.0 * xi3 * ( *(B1+3*i+2) * *phi11 + *(B2+3*i+2) * *phi12 + *(dB1xi1+3*i+2) * *phi1 + *(dB2xi1+3*i+2) * *phi2 );
          }
          *(betaX+3*0+0) += *phi1;
          *(betaY+3*1+0) += *phi1;
          *(betaZ+3*2+0) += *phi1;

          for(i=0; i<3; i++){
            *(betaX+3*i+1) = *thickness/2.0 * xi3 * ( *(B1+3*i+0) * *phi12 + *(B2+3*i+0) * *phi22 + *(dB1xi2+3*i+0) * *phi1 + *(dB2xi2+3*i+0) * *phi2 );
            *(betaY+3*i+1) = *thickness/2.0 * xi3 * ( *(B1+3*i+1) * *phi12 + *(B2+3*i+1) * *phi22 + *(dB1xi2+3*i+1) * *phi1 + *(dB2xi2+3*i+1) * *phi2 );
            *(betaZ+3*i+1) = *thickness/2.0 * xi3 * ( *(B1+3*i+2) * *phi12 + *(B2+3*i+2) * *phi22 + *(dB1xi2+3*i+2) * *phi1 + *(dB2xi2+3*i+2) * *phi2 );
          }
          *(betaX+3*0+1) += *phi2;
          *(betaY+3*1+1) += *phi2;
          *(betaZ+3*2+1) += *phi2;

          for(i=0; i<3; i++){
            *(betaX+3*i+2) = *thickness/2.0 * ( *(B1+3*i+0) * *phi1 + *(B2+3*i+0) * *phi2 );
            *(betaY+3*i+2) = *thickness/2.0 * ( *(B1+3*i+1) * *phi1 + *(B2+3*i+1) * *phi2 );
            *(betaZ+3*i+2) = *thickness/2.0 * ( *(B1+3*i+2) * *phi1 + *(B2+3*i+2) * *phi2 );
          }

          TX = *(aState+0);
          TY = *(aState+1);
          TZ = *(aState+2);

          for(j=0;j<3;j++){
            TX += *(aBar+j) * *(gamma+3*j+0);
            TY += *(aBar+j) * *(gamma+3*j+1);
            TZ += *(aBar+j) * *(gamma+3*j+2);
          }

          for(i=0;i<3;i++){
            for(j=0;j<3;j++){
              TX += *(bBar+3*i+j) * *(betaX+3*i+j);
              TY += *(bBar+3*i+j) * *(betaY+3*i+j);
              TZ += *(bBar+3*i+j) * *(betaZ+3*i+j);
            }
          }

          // multiply by the given gauss quadrature weight (thickness-level)
          TX *= w3;
          TY *= w3;
          TZ *= w3;

          *(forceDensity)   += TX * neighborArea;
          *(forceDensity+1) += TY * neighborArea;
          *(forceDensity+2) += TZ * neighborArea;
          *(neighborForceDensity)   -= TX * nodalArea;
          *(neighborForceDensity+1) -= TY * nodalArea;
          *(neighborForceDensity+2) -= TZ * nodalArea;
        }
      }
    }
  }
}


#undef  __FUNCT__
#define __FUNCT__ "computeRotationalForceState"
void computeRotationalForceState(PARAMETERS *par) 
{
  // The rotational part of the internal force is computed here
  // This part should have insignificant contributions for thin shells;
  // however, we keep it for the sake of completeness.

  double dt = par->timeStep;

  POINTS *point;
  POINTS *neighbor;
  int in,i,j,n;
  int numNeighbors, neigh;
  int *neighbors;

  double rho;
  double thickness;
  double* phi1;
  double* phi2;
  double* B1;
  double* B2;
  double* ndotN;
  double* ndotNP1;

  double nDdot[3];
  double nDdotB1[3], nDdotB2[3];
  double temp1, temp2;

  double *forceDensity, *neighborForceDensity;
  double TX, TY, TZ, nodalArea, neighborArea;

  // Zero out the force state 
  for(in=0;in<par->numPoints;in++){
    for(i=0; i<3; i++){
      par->puntos[in].rotationalForceDensity[i] = 0.0;
    }
  }

  // Loop over the material points and convert the Cauchy stress into pairwise peridynamic force densities
  for(in=0;in<par->numPoints;in++){

    point = &par->puntos[in];
    neighbors = point->neighbors;
    rho = point->density2D;
    thickness = point->shellThickness[1];
    B1 = point->B1Tensor;
    B2 = point->B2Tensor;
    ndotN = point->normalDotVector[0];
    ndotNP1 = point->normalDotVector[1];
    forceDensity = point->rotationalForceDensity;
    nodalArea = point->area;

    phi1 = point->gradientWeight1[0];
    phi2 = point->gradientWeight1[1];

    for(i=0; i<3; i++)
      nDdot[i] = ( *(ndotNP1+i) - *(ndotN+i) )/dt;

    for(j=0; j<3; j++){
      for(i=0; i<3; i++){
         nDdotB1[j] = nDdot[i] * *(B1+3*i+j);
         nDdotB2[j] = nDdot[i] * *(B2+3*i+j);
      }
    }

    // Loop over the neighbors and compute contribution to force densities
    numNeighbors = point->numNeighbors;
    for(n=0; n<numNeighbors; n++, phi1++, phi2++){

      neigh = neighbors[n];
      neighbor = &par->puntos[neigh];
      neighborArea = neighbor->area;
      neighborForceDensity = neighbor->rotationalForceDensity;

      temp1 = rho/12.0 * thickness * thickness * *phi1;
      temp2 = rho/12.0 * thickness * thickness * *phi2;

      TX = nDdotB1[0] * temp1 + nDdotB2[0] * temp2;
      TY = nDdotB1[1] * temp1 + nDdotB2[1] * temp2;
      TZ = nDdotB1[2] * temp1 + nDdotB2[2] * temp2;

      *(forceDensity)   += TX * neighborArea;
      *(forceDensity+1) += TY * neighborArea;
      *(forceDensity+2) += TZ * neighborArea;
      *(neighborForceDensity)   -= TX * nodalArea;
      *(neighborForceDensity+1) -= TY * nodalArea;
      *(neighborForceDensity+2) -= TZ * nodalArea;
    }
  }
}


#undef  __FUNCT__
#define __FUNCT__ "initializer"
void initializer(PARAMETERS *par) 
{
  // A set of functions to be called only in the initialization phase

  input(par); // Read the input geometry and initialize PD points

  directNeighborSearch(par); // Set up neighborhoods and allocate required memory

  initializeFields(par); // Initialize fields that are incrementally updated or have non-zero values

  formInitialCondition(par); // Set up the initial conditions

  constructTwoDimensionalSpace(par); // Construct parametric coordinates (2d space) from a 3d space, for the meshfree shell structure

  computeWeightedArea(par); // Calculate the weighted area 

  setRotationTensorToMaterialCoordinates(par); // Set the initial rotation tensor aligned with the material coordinates on the point level

  setBondLevelRotationTensorToMaterialCoordinates(par); // Set the initial rotation tensor aligned with the material coordinates on the bond level

  computeGradientShapeFunctions(par); // Compute the gradient weights at the initialization step

  computeReferenceDeformationGradient(par); // Compute the reference deformation gradient (mapping from parametric space to the undeformed physical space)
}


#undef  __FUNCT__
#define __FUNCT__ "computeVelocityGradient"
void computeVelocityGradient(PARAMETERS *par) 
{
  // A set of functions to be called for computing the velocity gradient tensor 

  computeGradientShapeFunctions(par); // Re-compute the gradient weights for nodes with damage growth

  updateShellThickness(par); // Update the shell thickness

  computeParametricGradients(par); // Compute the parametric gradients (with respect to the deformed configuration)

  computeNormalVectorAndRelatedTensors(par); // Compute the geometrical properties (normal vector and A, B1, B2 tensors)

  computeDeformationGradientInverse(par); // Compute the inverse of 3D deformation gradient along the thickness

  computeNodeLevelVelocityGradient(par); // Compute the 3D velocity gradient along the thickness at the point level

  computeBondLevelVelocityGradient(par); // Compute the 3D velocity gradient for each bond in all layers
}


#undef  __FUNCT__
#define __FUNCT__ "updateStress"
void updateStress(PARAMETERS *par) 
{
  // A set of functions to be called for updating the stress tensor 

  updateNodeLevelStress(par); // Update the stress at the point level (only for visualization)

  updateBondLevelStress(par); // Update the stress at the bond level
}


#undef  __FUNCT__
#define __FUNCT__ "computeForceState"
void computeForceState(PARAMETERS *par) 
{
  // The force state is computed here

  computeInternalForceState(par); // Compute the internal force state

  computeRotationalForceState(par); // Compute the rotational force state

  // Sum up the internal and rotational forces
  for(int in=0;in<par->numPoints;in++){
    for(int i=0; i<3; i++){
      par->puntos[in].forceDensity[i] = par->puntos[in].internalForceDensity[i] + par->puntos[in].rotationalForceDensity[i];
    }
  }
}


#undef  __FUNCT__
#define __FUNCT__ "computeAcceleration"
void computeAcceleration(PARAMETERS *par) 
{
  // The acceleration value is computed here by solving the E.O.M.
  
  POINTS *point;

  int dim = par->dim;
  int in,i;

  double *f, *b, *a, rho;

  for(in=0;in<par->numPoints;in++){
    point = &par->puntos[in];

    rho = point->density2D;
    f = point->forceDensity;
    b = point->bodyForceDensity;
    a = point->acceleration;

    for(i=0; i<dim; i++)
      *(a+i) = ( *(f+i) + *(b+i) ) / rho;
  }
}


#undef  __FUNCT__
#define __FUNCT__ "computeDamageByPrincipalStress"
void computeDamageByPrincipalStress(PARAMETERS *par) 
{
  // Compute damage using the principal stress failure criterion
  //
  // I use the first layer-associated bond-level variables to calculate damage
  // The same damage is applied to all layers (could be modified later)
  //
  // The bond breakage can be applied abruptly 
  // (by setting critcial and threshold damage variables to 1)
  // or can be applied with a transition (to enhance stability)
  double thresholdDamage = 0.99; 
  double criticalDamage = 1.0;

  double ultimateCompressionStrength = par->ultimateCompressionStrength;
  double ultimateTensionStrength = par->ultimateTensionStrength;

  POINTS *point;
  int in,n;
  int numNeighbors;

  double* bondDamageN;
  double* bondDamageNP1;
  double* damage;

  double* minPrincipalStress;
  double* maxPrincipalStress;

  double tempDmg, tempDmg1, tempDmg2;

  // Loop over the material points and convert the Cauchy stress into pairwise peridynamic force densities
  for(in=0;in<par->numPoints;in++){

    point = &par->puntos[in];
    damage = &(point->damage[1]);

    // Zero out fields
    *damage = 0.0;

    bondDamageN = point->bondDamage[0];
    bondDamageNP1 = point->bondDamage[1];
    minPrincipalStress = point->layerAssociatedBondLevelMinPrincipalStress[0]; // data from the first layer
    maxPrincipalStress = point->layerAssociatedBondLevelMaxPrincipalStress[0]; // data from the first layer

    // Loop over the neighbors and compute contribution to force densities
    numNeighbors = point->numNeighbors;
    for(n=0; n<numNeighbors; n++,
        bondDamageN++, bondDamageNP1++, 
        minPrincipalStress++, maxPrincipalStress++){

      if(*minPrincipalStress / ultimateCompressionStrength >= criticalDamage)
        tempDmg1  = 1.0;
      else if (*minPrincipalStress / ultimateCompressionStrength <= thresholdDamage)
        tempDmg1 = 0.0;
      else
        tempDmg1 = (*minPrincipalStress / ultimateCompressionStrength - thresholdDamage) / (criticalDamage - thresholdDamage);
      
      if(*maxPrincipalStress / ultimateTensionStrength >= criticalDamage)
        tempDmg2  = 1.0;
      else if (*maxPrincipalStress / ultimateTensionStrength <= thresholdDamage)
        tempDmg2 = 0.0;
      else
        tempDmg2 = (*maxPrincipalStress / ultimateTensionStrength - thresholdDamage) / (criticalDamage - thresholdDamage);

      tempDmg = fmax(tempDmg1, tempDmg2);
      
      *bondDamageNP1 = fmax(*bondDamageN, tempDmg); // non-decreasing damage 

      *damage += *bondDamageNP1;
    }

    *damage /= numNeighbors;
  }
}


#undef  __FUNCT__
#define __FUNCT__ "computeDamageByPlasticStrain"
void computeDamageByPlasticStrain(PARAMETERS *par) 
{
  // Compute damage using the plastic strain failure criterion
  //
  // I use the first layer-associated bond-level variables to calculate damage
  // The same damage is applied to all layers (could be modified later)
  //
  double thresholdEpsilonPlastic = par->thresholdEpsilonPlastic;
  double criticalEpsilonPlastic = par->criticalEpsilonPlastic;

  POINTS *point;
  int in,n;
  int numNeighbors;

  double* bondDamageN;
  double* bondDamageNP1;
  double* damage;

  double* eqps;

  double tempDmg;

  // Loop over the material points and convert the Cauchy stress into pairwise peridynamic force densities
  for(in=0;in<par->numPoints;in++){

    point = &par->puntos[in];
    damage = &(point->damage[1]);

    // Zero out fields
    *damage = 0.0;

    bondDamageN = point->bondDamage[0];
    bondDamageNP1 = point->bondDamage[1];
    eqps = point->layerAssociatedBondLevelEquivalentPlasticStrain[1][0]; // data from the first layer

    // Loop over the neighbors and compute contribution to force densities
    numNeighbors = point->numNeighbors;
    for(n=0; n<numNeighbors; n++,
        bondDamageN++, bondDamageNP1++, eqps++){

      if(*eqps >= criticalEpsilonPlastic)
        tempDmg  = 1.0;
      else if(*eqps <= thresholdEpsilonPlastic)
        tempDmg = 0.0;
      else
        tempDmg = (*eqps - thresholdEpsilonPlastic) / (criticalEpsilonPlastic - thresholdEpsilonPlastic);
      
      *bondDamageNP1 = fmax(*bondDamageN, tempDmg); // non-decreasing damage 

      *damage += *bondDamageNP1;
    }

    *damage /= numNeighbors;
  }
}


#undef  __FUNCT__
#define __FUNCT__ "computeDamage"
void computeDamage(PARAMETERS *par) 
{
  // Compute material damage here
  // Call the corresponding routine for different criteria

  if(!strcmp(par->failureCriterion, "NONE")){
    // no damage modeling
  }
  else if(!strcmp(par->failureCriterion, "Principal Stress")){
    // Principal Stress damage modeling (brittle fracture)
    computeDamageByPrincipalStress(par);
  }
  else if(!strcmp(par->failureCriterion, "Plastic Strain")){
    // Plastic Strain damage modeling (ductile fracture)
    computeDamageByPlasticStrain(par);
  }
  else{
    printf("Failure Criterion is not identified!\n");
    printf("Available options are 'Principal Stress', 'Principal Stress', 'NONE'.\n");
    exit(1);
  }
}


#undef  __FUNCT__
#define __FUNCT__ "updateKinematicsAtTheBeginningOfStep"
void updateKinematicsAtTheBeginningOfStep(PARAMETERS *par) 
{
  // The kinematics (velocities and positions) are evaluated at the beginning of each step
  // First step of velocity-Verlet algorithm is used to update the kinematics

  POINTS *point;

  int dim = par->dim;
  int in,i;

  double dt = par->timeStep;

  double* dispN;
  double* disp;
  double* coord;
  double* modelCoord;
  double* velN;
  double* vel;
  double* acc;

  // First update the velocities
  // V^{n+1/2} = V^{n} + (dt/2)*A^{n}
  for(in=0;in<par->numPoints;in++){
    point = &par->puntos[in];

    velN = point->velocity[0];
    vel = point->velocity[1];
    acc = point->acceleration;

    for(i=0; i<dim; i++)
      *(vel+i) = *(velN+i) + *(acc+i) * dt/2.0;
  }

  // Now prescribe the essential boundary conditions on the velocity DOFs
  prescribeEssentialBoundaryConditions(par); 

  // Update the positions now
  // U^{n+1} = U^{n} + (dt)*V^{n+1/2}
  for(in=0;in<par->numPoints;in++){
    point = &par->puntos[in];

    dispN = point->displacement[0];
    disp = point->displacement[1];
    coord = point->currentCoord[1];
    modelCoord = point->initialCoord;
    vel = point->velocity[1];

    for(i=0; i<dim; i++){
      *(disp+i) = *(dispN+i) + *(vel+i) * dt;
      *(coord+i) = *(modelCoord+i) + *(disp+i);
    }
  }
}


#undef  __FUNCT__
#define __FUNCT__ "updateKinematicsAtTheEndOfStep"
void updateKinematicsAtTheEndOfStep(PARAMETERS *par) 
{
  // The kinematics (velocities and positions) are evaluated at the end of each step
  // Second step of velocity-Verlet algorithm is used to update the kinematics
  
  POINTS *point;

  int dim = par->dim;
  int in,i;

  double dt = par->timeStep;

  double* vel;
  double* acc;

  // Velocities should be advanced from half step to the end of step
  // V^{n+1}   = V^{n+1/2} + (dt/2)*A^{n+1}
  for(in=0;in<par->numPoints;in++){
    point = &par->puntos[in];

    vel = point->velocity[1];
    acc = point->acceleration;

    for(i=0; i<dim; i++)
      *(vel+i) += *(acc+i) * dt/2.0;
  }
}


#undef  __FUNCT__
#define __FUNCT__ "updateStepNvalues"
void updateStepNvalues(PARAMETERS *par) 
{

  int in,n,j,l;
  POINTS *point;
  int dim = par->dim;

  for (in=0; in<par->numPoints ; in++){
    point = &par->puntos[in];

    point->shellThickness[0] = point->shellThickness[1]; 
    point->jacobianDeterminant[0] = point->jacobianDeterminant[1];
    point->equivalentPlasticStrain[0] = point->equivalentPlasticStrain[1];
    point->damage[0] = point->damage[1];

    for(j=0; j<dim; j++){
      point->currentCoord[0][j] = point->currentCoord[1][j]; 
      point->displacement[0][j] = point->displacement[1][j]; 
      point->velocity[0][j] = point->velocity[1][j]; 
      point->normalDotVector[0][j] = point->normalDotVector[1][j]; 
    }

    for (j=0;j<dim*dim; j++){ 
      point->unrotatedCauchyStress[0][j] = point->unrotatedCauchyStress[1][j]; 
      point->leftStretchTensor[0][j] = point->leftStretchTensor[1][j]; 
      point->rotationTensor[0][j] = point->rotationTensor[1][j]; 
    }

    for(n=0;n<point->numNeighbors;n++){
      point->bondDamage[0][n] = point->bondDamage[1][n]; 

      for(l=0;l<par->numLayers;l++){

        point->layerAssociatedBondLevelJacobianDeterminant[0][l][n] = point->layerAssociatedBondLevelJacobianDeterminant[1][l][n]; 
        point->layerAssociatedBondLevelEquivalentPlasticStrain[0][l][n] = point->layerAssociatedBondLevelEquivalentPlasticStrain[1][l][n]; 

        for (j=0;j<dim*dim; j++){ 
          point->layerAssociatedBondLevelUnrotatedCauchyStress[0][l][j][n] = point->layerAssociatedBondLevelUnrotatedCauchyStress[1][l][j][n]; 
          point->layerAssociatedBondLevelLeftStretchTensor[0][l][j][n] = point->layerAssociatedBondLevelLeftStretchTensor[1][l][j][n]; 
          point->layerAssociatedBondLevelRotationTensor[0][l][j][n] = point->layerAssociatedBondLevelRotationTensor[1][l][j][n]; 
        }
      }
    }
  }
}


#undef  __FUNCT__
#define __FUNCT__ "outputReults"
void outputReults(PARAMETERS *par, int counter) 
{
  // Write to output files 
  int numPoints = par->numPoints;

  // ##################################################
  //                  Geometry File
  // ##################################################

  FILE *fileResult;
  char filenameResults[50];

  sprintf(filenameResults,"Output/Meshless.%d.geo",counter);

  fileResult=fopen(filenameResults,"wt");
  if (fileResult == NULL){
    printf("Error opening file!\n");
    exit(1);
  }

  fprintf(fileResult,"Meshless \n");
  fprintf(fileResult,"node \n");
  fprintf(fileResult,"node id given \n");
  fprintf(fileResult,"element id given \n");
  fprintf(fileResult,"coordinates \n");

  fprintf(fileResult,"%d \n", numPoints);

  int i,j;

  for(i=0; i<numPoints;i++){
    POINTS *point = &par->puntos[i];
    fprintf(fileResult,"%d   %.4e   %.4e   %.4e \n", i+1, point->currentCoord[1][0], point->currentCoord[1][1], point->currentCoord[1][2]);
  }

  fprintf(fileResult,"part     1\n");
  fprintf(fileResult,"todo \n");
  fprintf(fileResult,"point \n");

  fprintf(fileResult,"%d \n", numPoints);

  for(j=1;j<=numPoints;j++)
    fprintf(fileResult,"%d  %d \n",j,j);

  fclose(fileResult);


  // ##################################################
  //                  Case File
  // ##################################################

  if (counter == 0){
    sprintf(filenameResults,"Output/pvFile.case");

    fileResult=fopen(filenameResults,"wt");
    if (fileResult == NULL){
      printf("Error opening file!\n");
      exit(1);
    }

    fprintf(fileResult,"#BOF: meshless.case \n");
    fprintf(fileResult,"\n");
    fprintf(fileResult,"FORMAT \n");
    fprintf(fileResult,"\n");
    fprintf(fileResult,"type: ensight \n");
    fprintf(fileResult,"\n");
    fprintf(fileResult,"GEOMETRY \n");
    fprintf(fileResult,"\n");
    fprintf(fileResult,"model: 1 Meshless.*.geo \n");
    fprintf(fileResult,"\n");
    fprintf(fileResult,"VARIABLE \n");
    fprintf(fileResult,"\n");
    fprintf(fileResult,"scalar per node: 1 Area Area.*.res \n");
    fprintf(fileResult,"vector per node: 1 Displacement Displacement.*.res \n");
    fprintf(fileResult,"vector per node: 1 Velocity Velocity.*.res \n");
    fprintf(fileResult,"vector per node: 1 Acceleration Acceleration.*.res \n");
    fprintf(fileResult,"vector per node: 1 ForceDensity ForceDensity.*.res \n");
    fprintf(fileResult,"tensor symm per node: 1 CauchyStress CauchyStress.*.res \n");
    fprintf(fileResult,"vector per node: 1 NormalVector NormalVector.*.res \n");
    fprintf(fileResult,"scalar per node: 1 PlasticStrain PlasticStrain.*.res \n");
    fprintf(fileResult,"scalar per node: 1 Damage Damage.*.res \n");
    fprintf(fileResult,"scalar per node: 1 VonMisesStress VonMisesStress.*.res \n");
    fprintf(fileResult,"scalar per node: 1 MinPrincipalStress MinPrincipalStress.*.res \n");
    fprintf(fileResult,"scalar per node: 1 MaxPrincipalStress MaxPrincipalStress.*.res \n");
    fprintf(fileResult,"\n");
    fprintf(fileResult,"TIME \n");
    fprintf(fileResult,"\n");

    fprintf(fileResult,"time set: 1 \n");
    fprintf(fileResult,"number of steps: %d \n", par->nOutputSteps);
    fprintf(fileResult,"filename start number: 0 \n");
    fprintf(fileResult,"filename increment: %d \n", 1);
    fprintf(fileResult,"time values: \n");

    for(i=0; i<par->nOutputSteps;i++){
      fprintf(fileResult,"%.5e  ", par->timeStep*par->FreqResults*i);

      if ((i+1) % 6 == 0) 
        fprintf(fileResult,"\n  ");
    }

    fclose(fileResult);
  }


  // ##################################################
  //                  Area File
  // ##################################################

  sprintf(filenameResults,"Output/Area.%d.res",counter);

  fileResult=fopen(filenameResults,"wt");
  if (fileResult == NULL){
    printf("Error opening file!\n");
    exit(1);
  }

  fprintf(fileResult,"Area \n");

  for(i=0; i<numPoints;i++){
    POINTS *point = par->puntos;
    fprintf( fileResult, "%.5e  ", point[i].area);

    if ((i+1) % 6 == 0) 
      fprintf(fileResult,"\n  ");
  }

  fclose(fileResult);


  // ##################################################
  //                  Displacement File
  // ##################################################

  sprintf(filenameResults,"Output/Displacement.%d.res",counter);

  fileResult=fopen(filenameResults,"wt");
  if (fileResult == NULL){
    printf("Error opening file!\n");
    exit(1);
  }

  fprintf(fileResult,"Displacement \n");

  for(i=0; i<numPoints;i++){
    POINTS *point = par->puntos;
    fprintf( fileResult, "%.5e   %.5e   %.5e   ", 
                 point[i].displacement[1][0],
                 point[i].displacement[1][1],
                 point[i].displacement[1][2]);

    if ((i+1) % 2 == 0) 
      fprintf(fileResult,"\n  ");
  }

  fclose(fileResult);


  // ##################################################
  //                  Velocity File
  // ##################################################

  sprintf(filenameResults,"Output/Velocity.%d.res",counter);

  fileResult=fopen(filenameResults,"wt");
  if (fileResult == NULL){
    printf("Error opening file!\n");
    exit(1);
  }

  fprintf(fileResult,"Velocity \n");

  for(i=0; i<numPoints;i++){
    POINTS *point = par->puntos;
    fprintf( fileResult, "%.5e   %.5e   %.5e   ", 
                 point[i].velocity[1][0],
                 point[i].velocity[1][1],
                 point[i].velocity[1][2]);

    if ((i+1) % 2 == 0) 
      fprintf(fileResult,"\n  ");
  }

  fclose(fileResult);


  // ##################################################
  //                  Acceleration File
  // ##################################################

  sprintf(filenameResults,"Output/Acceleration.%d.res",counter);

  fileResult=fopen(filenameResults,"wt");
  if (fileResult == NULL){
    printf("Error opening file!\n");
    exit(1);
  }

  fprintf(fileResult,"Acceleration \n");

  for(i=0; i<numPoints;i++){
    POINTS *point = par->puntos;
    fprintf( fileResult, "%.5e   %.5e   %.5e   ", 
                 point[i].acceleration[0],
                 point[i].acceleration[1],
                 point[i].acceleration[2]);

    if ((i+1) % 2 == 0) 
      fprintf(fileResult,"\n  ");
  }

  fclose(fileResult);


  // ##################################################
  //                  ForceDensity File
  // ##################################################

  sprintf(filenameResults,"Output/ForceDensity.%d.res",counter);

  fileResult=fopen(filenameResults,"wt");
  if (fileResult == NULL){
    printf("Error opening file!\n");
    exit(1);
  }

  fprintf(fileResult,"ForceDensity \n");

  for(i=0; i<numPoints;i++){
    POINTS *point = par->puntos;
    fprintf( fileResult, "%.5e   %.5e   %.5e   ", 
                 point[i].forceDensity[0],
                 point[i].forceDensity[1],
                 point[i].forceDensity[2]);

    if ((i+1) % 2 == 0) 
      fprintf(fileResult,"\n  ");
  }

  fclose(fileResult);


  // ##################################################
  //                  CauchyStress File
  // ##################################################

  sprintf(filenameResults,"Output/CauchyStress.%d.res",counter);

  fileResult=fopen(filenameResults,"wt");
  if (fileResult == NULL){
    printf("Error opening file!\n");
    exit(1);
  }

  fprintf(fileResult,"CauchyStress \n");

  for(i=0; i<numPoints;i++){
    POINTS *point = par->puntos;
    fprintf( fileResult, "%.5e  %.5e  %.5e  %.5e  %.5e  %.5e  ", point[i].cauchyStress[0], 
                                                                 point[i].cauchyStress[4],
                                                                 point[i].cauchyStress[8],
                                                                 point[i].cauchyStress[1],
                                                                 point[i].cauchyStress[2],
                                                                 point[i].cauchyStress[5]);

    fprintf(fileResult,"\n  ");
  }

  fclose(fileResult);


  // ##################################################
  //                  NormalVector File
  // ##################################################

  sprintf(filenameResults,"Output/NormalVector.%d.res",counter);

  fileResult=fopen(filenameResults,"wt");
  if (fileResult == NULL){
    printf("Error opening file!\n");
    exit(1);
  }

  fprintf(fileResult,"NormalVector \n");

  for(i=0; i<numPoints;i++){
    POINTS *point = par->puntos;
    fprintf( fileResult, "%.5e   %.5e   %.5e   ", 
                 point[i].normalVector[0],
                 point[i].normalVector[1],
                 point[i].normalVector[2]);

    if ((i+1) % 2 == 0) 
      fprintf(fileResult,"\n  ");
  }

  fclose(fileResult);


  // ##################################################
  //                  PlasticStrain File
  // ##################################################

  sprintf(filenameResults,"Output/PlasticStrain.%d.res",counter);

  fileResult=fopen(filenameResults,"wt");
  if (fileResult == NULL){
    printf("Error opening file!\n");
    exit(1);
  }

  fprintf(fileResult,"PlasticStrain \n");

  for(i=0; i<numPoints;i++){
    POINTS *point = par->puntos;
    fprintf( fileResult, "%.5e  ", point[i].equivalentPlasticStrain[1]);

    if ((i+1) % 6 == 0) 
      fprintf(fileResult,"\n  ");
  }

  fclose(fileResult);


  // ##################################################
  //                  Damage File
  // ##################################################

  sprintf(filenameResults,"Output/Damage.%d.res",counter);

  fileResult=fopen(filenameResults,"wt");
  if (fileResult == NULL){
    printf("Error opening file!\n");
    exit(1);
  }

  fprintf(fileResult,"Damage \n");

  for(i=0; i<numPoints;i++){
    POINTS *point = par->puntos;
    fprintf( fileResult, "%.5e  ", point[i].damage[1]);

    if ((i+1) % 6 == 0) 
      fprintf(fileResult,"\n  ");
  }

  fclose(fileResult);

  // ##################################################
  //                  VonMisesStress File
  // ##################################################

  sprintf(filenameResults,"Output/VonMisesStress.%d.res",counter);

  fileResult=fopen(filenameResults,"wt");
  if (fileResult == NULL){
    printf("Error opening file!\n");
    exit(1);
  }

  fprintf(fileResult,"VonMisesStress \n");

  for(i=0; i<numPoints;i++){
    POINTS *point = par->puntos;
    fprintf( fileResult, "%.5e  ", point[i].vonMisesStress);

    if ((i+1) % 6 == 0) 
      fprintf(fileResult,"\n  ");
  }

  fclose(fileResult);

  // ##################################################
  //                  MinPrincipalStress File
  // ##################################################

  sprintf(filenameResults,"Output/MinPrincipalStress.%d.res",counter);

  fileResult=fopen(filenameResults,"wt");
  if (fileResult == NULL){
    printf("Error opening file!\n");
    exit(1);
  }

  fprintf(fileResult,"MinPrincipalStress \n");

  for(i=0; i<numPoints;i++){
    POINTS *point = par->puntos;
    fprintf( fileResult, "%.5e  ", point[i].minPrincipalStress);

    if ((i+1) % 6 == 0) 
      fprintf(fileResult,"\n  ");
  }

  fclose(fileResult);

  // ##################################################
  //                  MaxPrincipalStress File
  // ##################################################

  sprintf(filenameResults,"Output/MaxPrincipalStress.%d.res",counter);

  fileResult=fopen(filenameResults,"wt");
  if (fileResult == NULL){
    printf("Error opening file!\n");
    exit(1);
  }

  fprintf(fileResult,"MaxPrincipalStress \n");

  for(i=0; i<numPoints;i++){
    POINTS *point = par->puntos;
    fprintf( fileResult, "%.5e  ", point[i].maxPrincipalStress);

    if ((i+1) % 6 == 0) 
      fprintf(fileResult,"\n  ");
  }

  fclose(fileResult);
}


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,
         char *argv[])
{

  //AppCtx user;
  //int i,j;
  PARAMETERS *par;
  par = (PARAMETERS*) malloc (sizeof(*par));
  
  // ------
  // initialize parameters
  par->dim=3; // Spatial dimension of the space
  par->kernelType = 3; // Cubic Spline (type of influence function)
  
  inputParameters(par); // read the rest of parameters from the user input file

  // Gauss points along the thickness
  if(par->numLayers == 1){
    par->xi3list[0] = 0.0; // Gaussian points along the thickness
    par->w3list[0] = 2.0; // Gaussian weights for each point
  }
  else if(par->numLayers == 2){
    par->xi3list[0] = -sqrt(1.0/3.0); par->xi3list[1] = sqrt(1.0/3.0); // Gaussian points along the thickness
    par->w3list[0] = 1.0; par->w3list[1] = 1.0; // Gaussian weights for each point
  }
  else if(par->numLayers == 3){
    par->xi3list[0] = 0.0; par->xi3list[1] = -sqrt(3.0/5.0); par->xi3list[2] = sqrt(3.0/5.0); // Gaussian points along the thickness
    par->w3list[0] = 8.0/9.0; par->w3list[1] = 5.0/9.0; par->w3list[2] = 5.0/9.0; // Gaussian weights for each point
  }
  else{
    printf("Number of discretized layers along the thickness is currently limited to 3!\n");
    exit(1);
  }
  // ------

  par->nSteps = int( floor((par->finalTime - par->initialTime)/par->timeStep) );
  par->nOutputSteps = par->nSteps/par->FreqResults + 1; // frequency for outputting results

  // we really need the bulk and shear moduli in the calculations (they can be directly input also)
  par->bulkModulus = par->youngModulus / (3.0*(1.0 - 2.0*par->poissonRatio));
  par->shearModulus = par->youngModulus / (2.0*(1.0 + par->poissonRatio));

  int displayTrigger = par->nSteps/100;
  if(displayTrigger == 0)
    displayTrigger = 1;

  // Perform functions that are called once, only in the initialization phase
  initializer(par);
  
  int outCounter = 0; // a counter for outputting
  outputReults(par,outCounter); // Output the initial settings 

  // Start the solver
  par->currentTime = par->initialTime;
  par->stepNumber = 0;
  while(par->stepNumber < par->nSteps){

    par->stepNumber++;

    if((par->stepNumber-1)%displayTrigger==0){
      printf("######################################################### \n");
      printf("%d percent complete!\n", (par->stepNumber-1)*100/par->nSteps);
      printf("Step Number: %d  Time step: %e, Time: %e \n", par->stepNumber, par->timeStep, par->currentTime);
      printf("######################################################### \n\n");
    }

    par->currentTime += par->timeStep/2.0; // Advance the time half a step
    
    computeDamage(par); // Update the damage using previous step values (for consistency purposes)

    updateKinematicsAtTheBeginningOfStep(par); // Perform the first step of velocity-Verlet integrator 

    computeVelocityGradient(par); // Compute the 3D velocity gradient 

    updateStress(par); // Update the stress tensor

    computeForceState(par); // Compute the force state

    computeAcceleration(par); // Compute the acceleration by solving the equation of motion

    par->currentTime += par->timeStep/2.0; // Advance the time another half a step

    updateKinematicsAtTheEndOfStep(par); // Perform the second step of velocity-Verlet integrator 

    updateStepNvalues(par); // Update the step N values for next time step

    if((par->stepNumber-1)%par->FreqResults==0){
      outCounter++;
      outputReults(par,outCounter); // Write to output files
    }
  }

  printf("######################################################### \n");
  printf("FINITO\n");
  printf("######################################################### \n");

  return 0;
}
