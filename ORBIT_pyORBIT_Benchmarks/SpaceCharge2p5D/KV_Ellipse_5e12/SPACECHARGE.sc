//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// DESCRIPTION
//
// Jeff Holmes 07/28/2010 Space Charge test for Xiyin
//
// Xiyin 12/07/2010 Benchmark with py-orbit
//
///////////////////////////////////////////////////////////////////////////

runName = "KV_Ellipse_5e12_1";

OFstream fio(runName+".out", ios::out);

///////////////////////////////////////////////////////////
// Make a synchronous particle and do some initialization:
///////////////////////////////////////////////////////////

  Real TSync = 1.0;     // Kinetic Energy (GeV)
  Real mSync = 1.0;     // Mass (AMU)
  Real charge = 1.0;    // charge number

  addSyncPart(mSync,charge,TSync);

  Real INTENSITY = 5.0e+12;
  nMaxMacroParticles = 100000;
  nReals_Macro = INTENSITY / nMaxMacroParticles;
  
  mainHerd = addMacroHerd(nMaxMacroParticles);
  readParts(mainHerd, "../DISTRIBUTIONS/Bm_KV_Ellipse", nMaxMacroParticles);

  Integer TURNS = 1;

  Real GAMMA = 1.0 + TSync / (0.93827231 * mSync);
  Real BETA = sqrt(1.0 - 1.0 / (GAMMA * GAMMA));
  Real CVEL = 2.997924580e+08;
  Real LRING = 4.0;
  Real TTURN = LRING / (BETA * CVEL);

//////////////////////////////
// Make  a Ring
//////////////////////////////

  const Integer nstepTPD = 1;
  const Integer nstepTPM = 4;
  const Integer fringeM = 1;
  const Integer nstepTPQ = 4;
  const Integer fringeQ = 1;
  const Integer nstepTPB = 10;
  const Integer fringeB = 1;
  const Integer nstepTPS = 4;
  const Integer nstepTPK = 4;

  buildTPlatticeNew("../LATTICES/Q_0p125.TP",
                    "Ring",
                    nstepTPD,
                    nstepTPM, fringeM,
                    nstepTPQ, fringeQ,
                    nstepTPB, fringeB,
                    nstepTPS,
                    nstepTPK);

//////////////////////////////
// Add a Longitudinal Impedance Node
//////////////////////////////

  ComplexVector ZImped(128);
  ZImped = Complex(0.,0.);

  nLongBins = 32;
  Real b_a = 10.0;
  Integer useAvg = 1;
  Integer nMacroLSCMin = 1;
  Integer useSpaceCharge = 0;

  addFFTLSpaceCharge("LSC1", 6, ZImped, b_a, useAvg, nMacroLSCMin, useSpaceCharge);
	
///////////////////////////////////////////
// Add a Transverse Space Charge Node Set
///////////////////////////////////////////

  Real eps = 1.0e-06;
  String BPShape = "Circle";
  
  Real BP1 = 110., BP2 = 0., BP3 = 0., BP4 = 0.;
  Integer BPPoints = 128, BPModes = 32;
  Real Gridfact = 2.0;
  Integer nMacroSCMin = 1;

  Integer nxBins = 128, nyBins = 128;
  addPotentialTransSCSet(nxBins, nyBins, eps,
                         BPShape, BP1, BP2, BP3, BP4,
                         BPPoints, BPModes, Gridfact, nMacroSCMin);

//////////////////////////////
// Start Output:
//////////////////////////////

  showStart(fio);
  showErrors(fio);

/////////////////////////////////////////////////////
// do some turns, and dump particles for later plots:
/////////////////////////////////////////////////////

  Integer nTurn;

  cerr << "Start Tracking\n";

  for(nTurn = 1; nTurn <= TURNS; nTurn++)
  {
      doTurn(1);
  }


  OFstream fio98(runName+"_output.dat", ios::out);
  dumpPartsGlobal(mainHerd, fio98);
  fio98.close();
  
  showTurnInfo(fio);
  showTurnInfo(cerr);
  
  fio.close();

  quit
