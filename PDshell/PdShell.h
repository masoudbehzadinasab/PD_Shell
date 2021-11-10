// Structure of each PD point and the fields it carries
typedef struct
{
  int  ID; //ID of the particle. 
  double initialCoord[3]; //Initial coordinates of the particle for 3D
  double referenceCoord[3]; //Reference coordinates of the particle for 3D
  double currentCoord[2][3]; //Current coordinates of the particle for 3D @ step N and N+1
  double displacement[2][3]; //Displacement of the particle between two consecutive time steps in the two directions @ step N and N+1
  double velocity[2][3]; //Velocity of the particle at the current time position @ step N and N+1
  double acceleration[3]; //Acceleration of the particle at the current time position

  int numNeighbors; //Number of neighbors for each point
  int *neighbors; //List of neighbors for each point
  double *influenceState; //Influence state 
  double *parametricCoordinates[2]; //Parametric coordinates of the neighbors
  double *gradientWeight1[2]; //First-order derivative shape functions
  double *gradientWeight2[3]; //Second-order derivative shape functions
  double *bondDamage[2]; //Bond damage associated with the neighbors @ step N and N+1

  double referenceShellThickness; //Initial shell thickness at each point
  double shellThickness[2]; //Shell thickness at each point @ step N and N+1

  double density2D; //Undeformed material mass density per unit area

  double horizon; //Support (horizon) for each point
  double area; //Area associated with each point. It is used as the integration weight in nodal integration.

  double weightedArea; //Area associated with each point. It is used as the integration weight in nodal integration.

  double parametricVelocityGradient1[2][3]; //first-order parametric velocity gradient on each particle. 3x2 for 3D. dvx/dxi1,dvx/dxi2,dvy/dxi1,dvy/dxi2,dvz/dxi1,dvz/dxi2 for indexing 0-->6
  double parametricVelocityGradient2[3][3]; //second-order parametric velocity gradient on each particle. 3x3 for 3D. dvx/dxi1dxi1,dvx/dxi1dxi2,dvx/dxi2dxi2,dvy/dxi1dxi1,dvy/dxi1dxi2,dvy/dxi2dxi2,dvz/dxi1dxi1,dvz/dxi1dxi2,dvz/dxi2dxi2 for indexing 0-->6

  double parametricDeformationGradient1[2][3]; //first-order parametric deformation gradient on each particle. 3x2 for 3D. 
  double parametricDeformationGradient2[3][3]; //second-order parametric deformation gradient on each particle. 3x3 for 3D.
  double referenceParametricDeformationGradient1[2][3]; //first-order parametric deformation gradient on each particle. 3x2 for 3D. For the undeformed configuration.
  double referenceParametricDeformationGradient2[3][3]; //second-order parametric deformation gradient on each particle. 3x3 for 3D. For the undeformed configuration.

  double parametricVector1[3]; //Psi1 vector for each point
  double parametricVector2[3]; //Psi2 vector for each point
  double etaVector[3]; //Eta vector for each point

  double expectedNormalVector[3]; //the expected direction of normal is input to make sure that PCA directions do not change sign sharply 
  //(i.e., the initial normals for neighboring points should not be in opposite directions; otherwise, it will be like a checkerboard and the method becomes unstable)
  double normalVector[3]; //normal vector on each particle
  double normalDotVector[2][3]; //normal dot vector on each particle @ step N and N+1
  double normalDDotVector[3]; //normal double dot vector on each particle

  double ATensor[9]; //A tensor on each particle
  double B1Tensor[9]; //B1 tensor on each particle
  double B2Tensor[9]; //B2 tensor on each particle
  double ATensorGradient[2][9]; //Parametric gradient of A tensor on each particle
  double B1TensorGradient[2][9]; //Parametric gradient of B1 tensor on each particle
  double B2TensorGradient[2][9]; //Parametric gradient of B2 tensor on each particle

  double deformationGradient[9]; //Deformation gradient tensor for each point
  double greenLagrangeStrain[9]; //Green Lagrange strain tensor for each point
  double unrotatedRateOfDeformation[9]; //Unrotated rate of deformation tensor for each point
  double leftStretchTensor[2][9]; //Left stretch tensor for each point @ step N and N+1
  double rotationTensor[2][9]; //Rotation tensor for each point @ step N and N+1
  double unrotatedCauchyStress[2][9]; //Unrotated Cauchy stress for each point @ step N and N+1
  double cauchyStress[9]; //Rotated Cauchy stress for each point
  double kirchhoffStress[9]; //Kirchhoff stress for each point
  double jacobianDeterminant[2]; //Jacobian determinant of volume change @ step N and N+1
  double equivalentPlasticStrain[2]; //Equivalent plastic strain. Basically a non decreasing quantity that accumulates the plastic strain in every time step @ step N and N+1
  double vonMisesStress; //Von Mises stress. This could be used as a damage criterion
  double minPrincipalStress; //Minimum principal stress. This could be used as a damage criterion (e.g. Tresca)
  double maxPrincipalStress; //Maximum principal stress. This could be used as a damage criterion (e.g. Tresca)

  double damage[2]; //Damage at previous and current stages
  bool weightEvaluationFlag; //Flag to indicate if the gradient weights needs to be (re-)computed

  double internalForceDensity[3]; //Internal force density for each point
  double rotationalForceDensity[3]; //Rotational force density for each point
  double forceDensity[3]; //Total force density for each point
  double bodyForceDensity[3]; //External body force density for each point

  // Currently 3 layers are considered for the shell
  // Layer-associated fields 
  double layerAssociatedVelocityGradient[3][9]; //3-d velocity gradient on each particle. 3x3 for 3D. dvx/dx,dvx/dy,dvx/dz,dvy/dx,dvy/dy,dvy/dz,dvz/dx,dvz/dy,dvz/dz for indexing 0-->9
  double layerAssociatedDeformationGradientInverse[3][9]; //Inverse of the 3-d deformation gradient on each particle. 3x3 for 3D. 
  double layerAssociatedAStateIntegral[3][3]; //Integral of a-state for each point in the layer
  double layerAssociatedBStateIntegral[3][9]; //Integral of b-state for each point in the layer

  // Currently 3 layers are considered for the shell
  // Layer-associated Bond-level fields 
  double *layerAssociatedBondDeformedState[3][3]; //Bond-level deformation state for each bond in each layer
  double *layerAssociatedBondLevelLeftStretchTensor[2][3][9]; //Left stretch tensor for each bond in each layer @ step N and N+1
  double *layerAssociatedBondLevelRotationTensor[2][3][9]; //Rotation tensor for each bond in each layer @ step N and N+1
  double *layerAssociatedBondLevelUnrotatedCauchyStress[2][3][9]; //Unrotated Cauchy stress for each bond in each layer @ step N and N+1
  double *layerAssociatedBondLevelCauchyStress[3][9]; //Rotated Cauchy stress for each bond in each layer
  double *layerAssociatedBondLevelKirchhoffStress[3][9]; //Kirchhoff stress for each bond in each layer
  double *layerAssociatedBondLevelUnrotatedRateOfDeformation[3][9]; //Unrotated rate of deformation tensor for each bond in each layer
  double *layerAssociatedBondLevelVelocityGradient[3][9]; //Velocity gradient for each bond in each layer
  double *layerAssociatedBondLevelJacobianDeterminant[2][3]; //Jacobian determinant of volume change for each bond in each layer @ step N and N+1
  double *layerAssociatedBondLevelEquivalentPlasticStrain[2][3]; //Equivalent plastic strain for each bond in each layer @ step N and N+1
  double *layerAssociatedBondLevelVonMisesStress[3]; //Von Mises stress for each bond in each layer
  double *layerAssociatedBondLevelMinPrincipalStress[3]; //Minimum principal stress for each bond in each layer (e.g. Tresca)
  double *layerAssociatedBondLevelMaxPrincipalStress[3]; //Maximum principal stress for each bond in each layer (e.g. Tresca)
} POINTS;


// Parameters for setting up a problem 
typedef struct
{

  int dim; //Spatial dimension of the space
  int kernelType; //Type of kernel for calculating influence function
  int orderOfBasis; //Order of basis for reproducing kernel
  double horizon; //PD horizon size

  // solver parameters
	double initialTime; //Initial time of our computation, usually 0
	double finalTime; //Final time of our computation. Problem dependent
	double currentTime; //Time currently (time in the time step we are in)
	double timeStep; //Time step, problem and mesh dependent
  int nSteps; //Number of steps for solving this problem
  int stepNumber; //The current step number

  // shell parameters
  // Currently 3 layers are considered for the shell
  //Therefore, three Gaussian point is used along the thickness
  int numLayers=3; //Discretized layers along the shell thickness
  double xi3list[3]; //Gaussian points along the thickness
  double w3list[3]; //Gaussian weights for each point
  double shellThickness; //Reference shell thickness

  // material properties
  double density; //Undeformed material mass density per unit volume
  double youngModulus; //Young's modulus
  double poissonRatio; //Poisson's ratio
  double bulkModulus; //Bulk modulus
  double shearModulus; //Shear modulus

  double initialYieldStress; //Initial yield stress under which the material starts plasticizing. 
  double saturatedYieldStress; //Saturated yield stress in this plasticity model
  double hardeningExponentialConstant; //The exponential constant in this plasticity model
  double hardeningLinearConstant; //The linear hardening constant in this plasticity model

  bool damageModeling; //Flag to indicate whether Damage Modeling is considered
  char failureCriterion[40]; //Type of failure criterion
  double ultimateCompressionStrength; //Ultimate Compression Strength for bond breakage
  double ultimateTensionStrength; //Ultimate Tension Strength for bond breakage
  double thresholdEpsilonPlastic; //Threshold Epsilon Plastic for bond degradation
  double criticalEpsilonPlastic; //Critical Epsilon Plastic for bond breakage

  int FreqResults; //Frequency with which we export results for post processing
  int nOutputSteps; //Number of times we export the outputs

  int numPoints; //total number of PD nodes
  POINTS *puntos; //structure POINT is a member of PARAMETERS

} PARAMETERS;

