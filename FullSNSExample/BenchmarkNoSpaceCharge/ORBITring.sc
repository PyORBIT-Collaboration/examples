
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

  Real TSync = 1.0;   // Kinetic Energy (GeV)
  Real mSync = 1.0;     // Mass (AMU)
  Real charge = 1.0;    // charge number

  addSyncPart(mSync,charge,TSync);

  mainHerd = addMacroHerd(10000);

  Integer TURNS = 1000;
  Real INTENSITY = 7.8e13;
  Real GAMMA = 1.0 + TSync / (0.93827231 * mSync);
  Real BETA = sqrt(1.0 - 1.0 / (GAMMA * GAMMA));
  Real CVEL = 2.997924580e+08;
  Real LRING = 248.009350;
  Real TTURN = LRING / (BETA * CVEL);
  Real TWAVETOT=0.910; // 910 us of injection kicker sqrt(t) falloff


//////////////////////////////
// Make  a Ring
//////////////////////////////

  Integer  Nextr = 0;

  const Integer nstepTPD = 1;
  const Integer fringeD = 1;
  const Integer nstepTPM = 4;
  const Integer fringeM = 1;
  const Integer nstepTPQ = 4;
  const Integer fringeQ = 1;
  const Integer nstepTPB = 10;
  const Integer fringeB = 1;
  const Integer nstepTPS = 4;
  const Integer fringeS = 1;
  const Integer nstepTPK = 4;
  const Integer fringeK = 1;


buildTPlatticeNew("MAD_Lattice/RealInjection/SNSring_pyOrbitBenchmark.TP","Ring",
                 nstepTPD,
                 nstepTPM, fringeM,
                 nstepTPQ, fringeQ,
                 nstepTPB, fringeB,
                 nstepTPS,
                 nstepTPK);


  Ring::gammaTrans = 5.245869;
  Ring::nuX=6.22;
  Ring::nuY=6.17;
  
///////////////////////////////////////////
// Make  Linac X & Y Distribution functions
////////////////////////////////////////////


  betaXInj = 10.209; alphaXInj = 0.063;
  epsXRMSInj = 0.152; MYJoho = 3.;
  epsXLimInj = epsXRMSInj * 2 * (MXJoho+1); 
//x0Inj = 46.8;   xP0Inj = 0;
  x0Inj = 0;   xP0Inj = 0;
  addXInitializer("Gaussian-X",JohoXDist);

  betaYInj = 10.776; alphaYInj = 0.063;
  epsYRMSInj = 0.152; MYJoho = 3.;
  epsYLimInj  = epsYRMSInj * 2 * (MYJoho+1);
//y0Inj = 49.2;   yP0Inj = 0;
  y0Inj = 0;   yP0Inj = 0;
  addYInitializer("Gaussian-Y",JohoYDist);

//Got the energy spread factor from Mike.  It's a measured number.

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
  NullTime  = 0.0;

addLongInitializer("SNS Spread",SNSESpreadDist);

////////////////////////////////
// Add information for injection 
// painting waveform
////////////////////////////////

//These parameters correspond to the parameters used to define the 
//kicker waveforms in the Waveform3() routin of TeaPot.cc

dhwf = 0.58;  // Fractional kicker amp at the end of 1000 turns
dvwf = 0.58; 
hwf = 1-dhwf;
vwf = 1-dvwf;
tstart = 0.0;
tstop = 1.0; // Total paint time in ms.


//////////////////////////////
// Add foil node
//////////////////////////////

  nxfh = 10;
  nyfh = 11;

  Real xfoilmin = x0Inj - 8.5;
  Real xfoilmax = x0Inj + 8.5;
  Real yfoilmin = y0Inj - 8.5;  
  Real yfoilmax = y0Inj + 100.0;

  useFoilScattering = 2;
  addFoil("Foil", 1, xfoilmin, xfoilmax, yfoilmin, yfoilmax, 400.0);

///////////////////////////////////////////
// Add a Transverse Space Charge Node Set
///////////////////////////////////////////

  Real eps = 1.0e-06;
  String BPShape = "Circle";
  Real BP1 = 110., BP2 = 0., BP3 = 0., BP4 = 0.;
  Integer BPPoints = 128, BPModes = 32;
  Real Gridfact = 2.0;
  Integer nMacroSCMin = 1;

  Integer nxBins = 64, nyBins = 64;
// addFFTTransSCSet(nxBins, nyBins, eps, nMacroSCMin);


//////////////////////////////////////
// Add aperutre
//////////////////////////////////////

globalAperture = 200.0;

//////////////////////////////
// Inject particles
//////////////////////////////

  nMacrosPerTurn = 260;
  nMaxMacroParticles = TURNS * nMacrosPerTurn;
  nReals_Macro = INTENSITY / nMaxMacroParticles;


///////////////////////////////////////////////////
// do some turns, and dump particles for later plots:
/////////////////////////////////////////////////////
  

  if( nRank == 0 ) cerr << "Start Tracking\n";
  Real et;
  timerOn();
  showNodes(cerr);

doTurn(1);;
  OFstream fio12("Bm_Parts_1"+sRank, ios::out);
  dumpPartsGlobal(mainHerd, fio12);
  fio12.close();

  doTurn(99);
  OFstream fio12("Bm_Parts_100"+sRank, ios::out);
  dumpPartsGlobal(mainHerd, fio12);
  fio12.close();

  doTurn(100);
  OFstream fio12("Bm_Parts_200"+sRank, ios::out);
  dumpPartsGlobal(mainHerd, fio12);
  fio12.close();

  doTurn(100);
  OFstream fio12("Bm_Parts_300"+sRank, ios::out);
  dumpPartsGlobal(mainHerd, fio12);

  doTurn(100);
  OFstream fio12("Bm_Parts_400"+sRank, ios::out);
  dumpPartsGlobal(mainHerd, fio12);
  fio12.close();

  doTurn(100);
  OFstream fio12("Bm_Parts_500"+sRank, ios::out);
  dumpPartsGlobal(mainHerd, fio12);
  fio12.close();

  doTurn(100);
  OFstream fio12("Bm_Parts_600"+sRank, ios::out);
  dumpPartsGlobal(mainHerd, fio12);
  fio12.close();

  doTurn(100);
  OFstream fio12("Bm_Parts_700"+sRank, ios::out);
  dumpPartsGlobal(mainHerd, fio12);
  fio12.close();

  doTurn(100);
  OFstream fio12("Bm_Parts_800"+sRank, ios::out);
  dumpPartsGlobal(mainHerd, fio12);
  fio12.close();

  doTurn(100);
  OFstream fio12("Bm_Parts_900"+sRank, ios::out);
  dumpPartsGlobal(mainHerd, fio12);
  fio12.close();

  doTurn(100);
  OFstream fio12("Bm_Parts_1000"+sRank, ios::out);
  dumpPartsGlobal(mainHerd, fio12);
  fio12.close();

  quit


