// User should specify these functions with respect to a given problem

#undef __FUNCT__
#define __FUNCT__ "inputParameters"
void inputParameters(PARAMETERS *par)
{
  // Assign the input parameters here
  // User should specify with respect to the problem
  
  par->shellThickness = 1.0; // Reference shell thickness
  par->orderOfBasis = 2; // Order of accuracy (polynomial reproducibility)
  par->horizon = 2.5; // PD horizon for each point

  // solver parameters
  par->initialTime = 0.0; //Initial time of our computation, usually 0
  par->finalTime = 5.0e-3; //Final time of our computation. Problem dependent
  par->timeStep	= 1.0e-7; //Time step, problem and mesh dependent
  par->FreqResults = 100; // frequency for outputting results

  // material parameters
  par->density = 7.8e-6; //Undeformed material mass density per unit volume
  par->youngModulus = 200.0e6; //Young's modulus
  par->poissonRatio = 0.3; //Poisson's ratio

  // exponential saturation type hardening rule
  // Y = Y_0 + (Y_sat - Y_0)(1 - exp( - expCst * eqps)) + linCst * eqps
  par->initialYieldStress = 100.0e3; //Initial yield stress under which the material starts plasticizing. 
  par->saturatedYieldStress = 200.0e3; //Saturated yield stress in this plasticity model
  par->hardeningExponentialConstant = 1.0; //The exponential constant in this plasticity model
  par->hardeningLinearConstant = 100.0e3; //The linear hardening constant in this plasticity model

  //------
  // damage parameters (current damage models = "Principal Stress" "Plastic Strain"
  
  //strcpy(par->failureCriterion, "NONE"); // specify this as "NONE" if you don't want to include damage

  //strcpy(par->failureCriterion, "Principal Stress");
  //par->failureStress = 8.0e10;
  //par->ultimateCompressionStrength; //Ultimate Compression Strength for bond breakage
  //par->ultimateTensionStrength; //Ultimate Tension Strength for bond breakage

  strcpy(par->failureCriterion, "Plastic Strain");
  par->thresholdEpsilonPlastic = 0.2; //Threshold Epsilon Plastic for bond degradation
  par->criticalEpsilonPlastic = 0.25; //Critical Epsilon Plastic for bond breakage
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

  double* initCoord;
  double* vel;

  // With respect to the current format of the code, 
  // essential BCs should be prescribed on velocities
  for(in=0;in<par->numPoints;in++){
    point = &par->puntos[in];

    initCoord = point->initialCoord;
    vel = point->velocity[1];

    // uni-axial tension in Y
    if(initCoord[1]<-25.0){
      vel[0] = 0.0;
      vel[1] = -1000.0;
      vel[2] = 0.0;
    }

    if(initCoord[1]>25.0){
      vel[0] = 0.0;
      vel[1] = 1000.0;
      vel[2] = 0.0;
    }
  }
}

