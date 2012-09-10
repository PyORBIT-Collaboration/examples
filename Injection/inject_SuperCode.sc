
//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// DESCRIPTION
//
// Sarah Cousineau -June 2010.  Simulation of May 2010 benchmark data.
//
///////////////////////////////////////////////////////////////////////////

////////////////////
//MPI stuff
////////////////////

  Integer nRank=MPI_rank();
  String sRank = "_";
  MPI_String_Add_Integer(sRank,nRank);

////////////////////
// Files for output:
////////////////////

  String runName, of1;
  runName = "Run_May";
  of1 = runName + sRank + ".out";

  OFstream fio(of1, ios::out);

///////////////////////////////////////////////////////////
// Make a synchronous particle and do some initialization:
///////////////////////////////////////////////////////////

  Real TSync = 0.930;   // Kinetic Energy (GeV)
  Real mSync = 1.0;     // Mass (AMU)
  Real charge = 1.0;    // charge number

  addSyncPart(mSync,charge,TSync);

  mainHerd = addMacroHerd(10000);

  Integer TURNS = 600;
  Real INTENSITY = 7.8e13;
  Real GAMMA = 1.0 + TSync / (0.93827231 * mSync);
  Real BETA = sqrt(1.0 - 1.0 / (GAMMA * GAMMA));
  Real CVEL = 2.997924580e+08;
  Real LRING = 248.009350;
  Real TTURN = LRING / (BETA * CVEL);
  Real TWAVETOT=0.910; // 910 us of injection kicker sqrt(t) falloff


///////////////////////////////////////////
// Make  Linac X & Y Distribution functions
////////////////////////////////////////////
//Real xpaint = 21.0, xppaint = 0.0;
// Real ypaint = 25.0, yppaint = 0.0;
// Real dxpaint = 44.07 - xpaint;
// Real dypaint = 46.00 - ypaint;

//Measured spot relative to c.o
  Real xinjspot = 16.4, pxinjspot=1.0;
  Real yinjspot = 28.4, pyinjspot=-0.06;
//Kickers were flat-topped.
//This is the initial bump calculated from inj. kicker settings at turn 0
  Real xbump = 30.4, xppaint = 0.0;
  Real ybump = 20.8, yppaint = 0.0;

  betaXInj = 10.209; alphaXInj = 0.063;
//  epsXRMSInj = 0.158; MXJoho = 3.;
  epsXRMSInj = 0.00152; MYJoho = 3.;
  epsXLimInj = epsXRMSInj * 2 * (MXJoho+1);
  x0Inj = xbump + xinjspot;   xP0Inj = pxinjspot;
  addXInitializer("Gaussian-X",JohoXDist);

  betaYInj = 10.776; alphaYInj = 0.063;
//  epsYRMSInj = 0.152; MYJoho = 3.;
  epsYRMSInj = 0.00152; MYJoho = 3.;
  epsYLimInj  = epsYRMSInj * 2 * (MYJoho+1);
  y0Inj = ybump + yinjspot;   yP0Inj = pyinjspot;
  addYInitializer("Gaussian-Y",JohoYDist);

//Got the energy spread factor from Mike.  It's a measured number.
/*
  Real efac = 0.784;
  phiMaxInj = 120.;
  phiMinInj = -120.;
  EInjMean  = TSync;
  EInjSigma = 0.001*efac;
  EInjMin   = TSync - 0.0025*efac;
  EInjMax   = TSync + 0.0025*efac;
  ECMean    = 0.0;
  ECSigma   = 0.0015*efac;       // Random Centroid Jitter
  ECMin     = -0.0035*efac;      // Random Centroid Jitter
  ECMax     = 0.0035*efac;       // Random Centroid Jitter
  ECDriftI  = -0.000;
  ECDriftF  =  0.000;
  DriftTime = 1000 * TURNS * TTURN;
  ESNu      = 100.;
  ESPhase    = 0.;
  ESMax     = 0.0;
*/

MLJoho = 3;
phiLimInj = 120.;
dELimInj = 0.001;
nLongInjBunch = 100;
deltaPhiBunch = 1;
deltaPhiNotch = 0;
lTailFraction = 0;
lTailFactor = 0;

  
/*
  ECSigma   = 0.0015;       // Random Centroid Jitter
  ECSigma   = 0.000000001;  // No Centroid Jitter
  ECSigma   = 0.000388;     // Corrected Centroid Jitter
  ECMin     = -0.005;
  ECMax     = 0.005;
  ESMax     = 0.004;
*/
  addLongInitializer("Uniform -L",JohoLDist);
//addLongInitializer("Uniform -L",SNSESpreadDist);


//////////////////////////////
// Add foil node
//////////////////////////////

  nxfh = 10;
  nyfh = 11;

  Real xfoilmin = x0Inj - 100.0;
  Real yfoilmin = y0Inj - 100.0;
  Real xfoilmax = x0Inj + 100.0;
  Real yfoilmax = y0Inj + 100.0;

  useFoilScattering = 0;
  addFoil("Foil", 1, xfoilmin, xfoilmax, yfoilmin, yfoilmax,0.);
 

//////////////////////////////
// Inject particles
//////////////////////////////

  nMacrosPerTurn = 10000;
  nMaxMacroParticles = TURNS * nMacrosPerTurn;
  nReals_Macro = INTENSITY / nMaxMacroParticles;


///////////////////////////////////////////////////
// do some turns, and dump particles for later plots:
/////////////////////////////////////////////////////
  

  if( nRank == 0 ) cerr << "Start Tracking\n";
  Real et;
  timerOn();
  showNodes(cerr);

  turnToNode(mainHerd,1);
  OFstream fio12("Bm_Parts_000"+sRank, ios::out);
  dumpPartsGlobal(mainHerd, fio12);
  fio12.close();

  quit


