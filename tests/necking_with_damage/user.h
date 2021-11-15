// For details of this problem, look at example 4.3.1 in the following paper:
// 
// Behzadinasab et al. "A General-Purpose, Inelastic, Rotation-Free 
// Kirchhoff–Love Shell Formulation for Peridynamics", Computer Methods 
// in Applied Mechanics and Engineering, 2021.
//
// User should specify these functions with respect to a given problem

#undef __FUNCT__
#define __FUNCT__ "inputParameters"
void inputParameters(PARAMETERS *par)
{
  // Assign the input parameters here
  // User should specify with respect to the problem
  
  // neighborhood parameters
  par->orderOfBasis = 2; // Order of accuracy (polynomial reproducibility)
  par->horizon = 1.2; // PD horizon for each point

  // shell parameters
  par->shellThickness = 1.0; // Reference shell thickness
  par->numLayers = 3; // Number of discretized layers along the shell thickness (1, 2, or 3)

  // solver parameters
  par->initialTime = 0.0; //Initial time of our computation, usually 0
  par->finalTime = 8.0e-3; //Final time of our computation. Problem dependent
  par->timeStep	= 2.0e-8; //Time step, problem and mesh dependent
  par->FreqResults = 1000; // frequency for outputting results

  // material parameters
  par->density = 1.0e-9; //Undeformed material mass density per unit volume
  par->youngModulus = 189.0e3; //Young's modulus
  par->poissonRatio = 0.29; //Poisson's ratio

  // exponential saturation type hardening rule
  // Y = Y_0 + (Y_sat - Y_0)(1 - exp( - expCst * eqps)) + linCst * eqps
  par->initialYieldStress = 343.0; //Initial yield stress under which the material starts plasticizing. 
  par->saturatedYieldStress = 680.0; //Saturated yield stress in this plasticity model
  par->hardeningExponentialConstant = 16.93; //The exponential constant in this plasticity model
  par->hardeningLinearConstant = 300.0; //The linear hardening constant in this plasticity model

  //------
  // damage parameters (current damage models = "Principal Stress" "Plastic Strain"
  
  //strcpy(par->failureCriterion, "NONE"); // specify this as "NONE" if you don't want to include damage

  //strcpy(par->failureCriterion, "Principal Stress");
  //par->failureStress = 8.0e10;
  //par->ultimateCompressionStrength; //Ultimate Compression Strength for bond breakage
  //par->ultimateTensionStrength; //Ultimate Tension Strength for bond breakage

  strcpy(par->failureCriterion, "Plastic Strain");
  par->thresholdEpsilonPlastic = 0.55; //Threshold Epsilon Plastic for bond degradation
  par->criticalEpsilonPlastic = 0.8; //Critical Epsilon Plastic for bond breakage
  //------

}


#undef __FUNCT__
#define __FUNCT__ "formInitialCondition"
void formInitialCondition(PARAMETERS *par)
{
  // Assign the initial conditions here
  // User should specify with respect to the problem
  
  POINTS *point;
  int in,j;
  int dim = par->dim;

  double *modelCoord;
  double *expNormalDir;
  double *bodyForceDensity;

  for (in=0; in<par->numPoints ; in++){
    point = &par->puntos[in];

    expNormalDir = point->expectedNormalVector;
    bodyForceDensity = point->bodyForceDensity;

    // Normal vector in the initial geometry 
    // This does not be to exact
    // Just needs to be in the direction that the initial normals are pointing out
    // This will be used to make sure that PCA normals remain consistent (PCA
    // can give result in a normal N or its opposite direction -N
    expNormalDir[0] = 0.0; 
    expNormalDir[1] = 0.0; 
    expNormalDir[2] = 1.0; 

    // Body force in this problem
    // (if it changes with time, move it to boundary conditions)
    for(j=0; j<dim; j++)
      bodyForceDensity[j] = 0.0;
  }
}


#undef __FUNCT__
#define __FUNCT__ "prescribeEssentialBoundaryConditions"
void prescribeEssentialBoundaryConditions(PARAMETERS *par)
{
  // Prescribe the essential boundary conditions here
  // User should specify with respect to the problem

  POINTS *point;

  int dim = par->dim;
  int in,i;

  double t = par->currentTime;
  double t0 = 1.0e-5;
  double pi2 = 1.57079632679; // pi/2

  double* initCoord;
  double* vel;

  // With respect to the current format of the code, 
  // essential BCs should be prescribed on velocities
  for(in=0;in<par->numPoints;in++){
    point = &par->puntos[in];

    initCoord = point->initialCoord;
    vel = point->velocity[1];

    // fixed nodes in the bottom row
    if(initCoord[1]<0.5){
      vel[0] = 0.0;
      vel[1] = 0.0;
      vel[2] = 0.0;
    }

    // constant y-velocity for the top row and fixed x and z coordinates
    if(initCoord[1]>49.5){
      vel[0] = 0.0;
      vel[2] = 0.0;

      // ramp up the velocity smoothly
      if(t < t0)
        vel[1] = 1500.0 * sin(pi2*t/t0);
      else
        vel[1] = 1500.0;
    }
  }
}

