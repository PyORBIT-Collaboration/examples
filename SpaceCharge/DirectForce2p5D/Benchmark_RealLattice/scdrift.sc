
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


buildTPlatticeNew("../MAD_Lattice/SNSring_pyOrbitBenchmark_noKick.TP","Ring",
                 nstepTPD,
                 nstepTPM, fringeM,
                 nstepTPQ, fringeQ,
                 nstepTPB, fringeB,
                 nstepTPS,
                 nstepTPK);


  Ring::gammaTrans = 5.245869;
  Ring::nuX=6.22;
  Ring::nuY=6.17;

//////////////////////////////
// Add Beam:
//////////////////////////////

  nMaxMacroParticles = 10000;
  nReals_Macro = 1.0e+14 / nMaxMacroParticles;
  readParts(mainHerd, "kv_orbit.dat", nMaxMacroParticles);

  nMacrosPerTurn = 0;

///////////////////////////////////////////
// Add a Transverse Space Charge Node Set
///////////////////////////////////////////

  Real eps = 1.0e-06;
  Integer nMacroSCMin = 1;
  Real length = 248.0;
  Integer nxBins = 16, nyBins = 16;
addFFTTransSCSet(nxBins, nyBins, eps, nMacroSCMin);

//////////////////////////////////////////
// Do a turn
//////////////////////////////////////////

  if( nRank == 0 ) cerr << "Start Tracking\n";
  Real et;
  timerOn();
  showNodes(cerr);

  doTurn(1);
  OFstream fio12("Bm_Parts1", ios::out);
  dumpParts(mainHerd, fio12);
  fio12.close();
  
doTurn(1);
  OFstream fio12("Bm_Parts2", ios::out);
  dumpParts(mainHerd, fio12);
  fio12.close();
  doTurn(1);

  OFstream fio12("Bm_Parts3", ios::out);
  dumpParts(mainHerd, fio12);
  fio12.close();
  
doTurn(1);
  OFstream fio12("Bm_Parts4", ios::out);
  dumpParts(mainHerd, fio12);
  fio12.close();

doTurn(1);
  OFstream fio12("Bm_Parts5", ios::out);
  dumpParts(mainHerd, fio12);
  fio12.close();
  fio.close();

doTurn(95);
OFstream fio12("Bm_Parts100", ios::out);
  dumpParts(mainHerd, fio12);
  fio12.close();
  fio.close();


  quit


