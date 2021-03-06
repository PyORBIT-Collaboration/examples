
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


buildTPlatticeNew("MAD_Lattice/SNSring_pyOrbitBenchmark.TP","Ring",
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
  x0Inj = 46.8;   xP0Inj = 0;
//x0Inj = 0;   xP0Inj = 0;
  addXInitializer("Gaussian-X",JohoXDist);

  betaYInj = 10.776; alphaYInj = 0.063;
  epsYRMSInj = 0.152; MYJoho = 3.;
  epsYLimInj  = epsYRMSInj * 2 * (MYJoho+1);
  y0Inj = 49.2;   yP0Inj = 0;
//y0Inj = 0;   yP0Inj = 0;
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
  Integer nMacroSCMin = 1000;
  Real sc_path_length_min = 0.00000001;
  Integer nxBins = 64, nyBins = 64;
  addFFTTransSCSet(nxBins, nyBins, eps, nMacroSCMin);
  

//////////////////////////////////////
// Add aperture
//////////////////////////////////////

globalAperture = 200.0;

//////////////////////////////
// Inject particles
//////////////////////////////

  nMacrosPerTurn = 260;
  nMaxMacroParticles = TURNS * nMacrosPerTurn;
  nReals_Macro = INTENSITY / nMaxMacroParticles;

//////////////////////////////
// Add RF Cavities
//////////////////////////////

  RealVector harmNum(2);
  harmNum(1) = 1.; harmNum(2) = 2.;

  Integer nRFHarms1 = 1;
  RealVector volts1(nRFHarms1), harmNum1(nRFHarms1), RFPhase1(nRFHarms1);
  volts1(1) = 16.0;
  harmNum1(1) = harmNum(1);
  RFPhase1(1) = 0;

  Integer nRFHarms2 = 2;
  RealVector volts2(nRFHarms2), harmNum2(nRFHarms2), RFPhase2(nRFHarms2);
  volts2(1) = 0.; volts2(2) = -3.0;
  harmNum2(1) = harmNum(1); harmNum2(2) = harmNum(2);
  RFPhase2(1) = 0.; RFPhase2(2) = 0.;

  addRFCavity("RF 1", 4767, nRFHarms1, volts1, harmNum1, RFPhase1);
  addRFCavity("RF 2", 4787, nRFHarms2, volts2, harmNum2, RFPhase2);

//////////////////////////////
// Add a Longitudinal Impedance Node
//////////////////////////////

  ComplexVector ZImpExtr(128), ZImpRFCav(128), ZImped(128);
  ZImpExtr = Complex(0.,0.);
  ZImpRFCav = Complex(0.,0.);
  ZImped = Complex(0.,0.);

/*
  define ZImped for SNS modes 1-32.
  eg.
  values are real & imaginary components in Ohm/n
  ZImped() by mode number
*/

  ZImpExtr(1) = Complex(42., 182.);
  ZImpExtr(2) = Complex(35., 101.5);
  ZImpExtr(3) = Complex(30.3333, 74.6667);
  ZImpExtr(4) = Complex(31.5, 66.5);
  ZImpExtr(5) = Complex(32.2, 57.4);
  ZImpExtr(6) = Complex(31.5, 51.3333);
  ZImpExtr(7) = Complex(31., 49.);
  ZImpExtr(8) = Complex(31.5, 46.375);
  ZImpExtr(9) = Complex(31.8889, 43.5556);
  ZImpExtr(10) = Complex(32.9, 40.6);
  ZImpExtr(11) = Complex(32.7273, 38.1818);
  ZImpExtr(12) = Complex(32.25, 35.58333);
  ZImpExtr(13) = Complex(34.46, 32.846);
  ZImpExtr(14) = Complex(35., 30.5);
  ZImpExtr(15) = Complex(35.4667, 28.);
  ZImpExtr(16) = Complex(36.75, 25.8125);
  ZImpExtr(17) = Complex(36.647, 23.88);
  ZImpExtr(18) = Complex(36.944, 22.1667);
  ZImpExtr(19) = Complex(36.474, 20.263);
  ZImpExtr(20) = Complex(36.4, 18.55);
  ZImpExtr(21) = Complex(35.333, 17.);
  ZImpExtr(22) = Complex(35.,14.9545 );
  ZImpExtr(23) = Complex(33.478, 13.696);
  ZImpExtr(24) = Complex(32.375, 11.6667);
  ZImpExtr(25) = Complex(30.8, 10.08);
  ZImpExtr(26) = Complex(29.615, 8.077);
  ZImpExtr(27) = Complex(28.519, 6.741);
  ZImpExtr(28) = Complex(27.5, 5.);
  ZImpExtr(29) = Complex(26.552, 4.103);
  ZImpExtr(30) = Complex(25.4333, 3.26667);
  ZImpExtr(31) = Complex(24.3871, 2.7097);
  ZImpExtr(32) = Complex(23.40625, 2.1875);

  ZImpRFCav(1)  = Complex(0.000, 0.0);
  ZImpRFCav(2)  = Complex(0.750, 0.0);
  ZImpRFCav(3)  = Complex(0.333, 0.0);
  ZImpRFCav(4)  = Complex(0.250, 0.0);
  ZImpRFCav(5)  = Complex(0.200, 0.0);
  ZImpRFCav(6)  = Complex(0.167, 0.0);
  ZImpRFCav(7)  = Complex(3.214, 0.0);
  ZImpRFCav(8)  = Complex(0.188, 0.0);
  ZImpRFCav(9)  = Complex(0.167, 0.0);
  ZImpRFCav(10) = Complex(0.150, 0.0);
  ZImpRFCav(11) = Complex(1.000, 0.0);
  ZImpRFCav(12) = Complex(0.125, 0.0);
  ZImpRFCav(13) = Complex(0.115, 0.0);
  ZImpRFCav(14) = Complex(0.143, 0.0);
  ZImpRFCav(15) = Complex(0.333, 0.0);
  ZImpRFCav(16) = Complex(0.313, 0.0);
  ZImpRFCav(17) = Complex(0.294, 0.0);
  ZImpRFCav(18) = Complex(0.278, 0.0);
  ZImpRFCav(19) = Complex(0.263, 0.0);
  ZImpRFCav(20) = Complex(0.250, 0.0);
  ZImpRFCav(21) = Complex(0.714, 0.0);
  ZImpRFCav(22) = Complex(0.682, 0.0);
  ZImpRFCav(23) = Complex(0.652, 0.0);
  ZImpRFCav(24) = Complex(0.625, 0.0);
  ZImpRFCav(25) = Complex(0.600, 0.0);
  ZImpRFCav(26) = Complex(0.577, 0.0);
  ZImpRFCav(27) = Complex(0.536, 0.0);
  ZImpRFCav(28) = Complex(0.536, 0.0);
  ZImpRFCav(29) = Complex(0.517, 0.0);
  ZImpRFCav(30) = Complex(0.500, 0.0);
  ZImpRFCav(31) = Complex(0.484, 0.0);
  ZImpRFCav(32) = Complex(0.469, 0.0);

  Integer harm;
  Real realp, imagp;

  for (harm = 1; harm <= 32; harm++)
  {
    realp = ZImpExtr(harm).re / 1.75 + ZImpRFCav(harm).re;
    imagp = ZImpExtr(harm).im / 1.75 + ZImpRFCav(harm).im;
    ZImped(harm) = Complex(realp, imagp);
  }
  Integer maxModes = 32;

  nLongBins = 128;
  Real b_a = 10.0/3.0;
  Integer useAvg = 0;
  Integer nMacroLSCMin = 1000;

  Integer useSpaceCharge = 1;
  addFFTLSpaceCharge("LSC1", 100, ZImped, b_a, useAvg,
                    nMacroLSCMin, useSpaceCharge);

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


