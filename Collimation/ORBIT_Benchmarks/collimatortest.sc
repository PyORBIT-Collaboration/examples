
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

  Real length = 0.5;
  addTPD("Drift", 1, 0, 0, 0, 0, 0, 0, length, "Drift", nstepTPD);

//////////////////////////////
// Add Beam:
//////////////////////////////

  nMaxMacroParticles = 10000;
  nReals_Macro = 1.0e+10 / nMaxMacroParticles;

  readParts(mainHerd, "Bm_KV_Uniform_10000", nMaxMacroParticles);

//////////////////////////////
// Add RF Cavities
//////////////////////////////


addCollimator("P1", 1, length, 4, 1, 1, 0.00, 0.00, 0., 0., 0.);


  if( nRank == 0 ) cerr << "Start Tracking\n";
  Real et;
  timerOn();
  showNodes(cerr);

  doTurn(1);
  OFstream fio12("Bm_Parts", ios::out);
  dumpParts(mainHerd, fio12);
  OFstream fio13("Lost_Parts", ios::out);
  dumpLostParts(mainHerd, fio13);

  fio12.close();
  fio13.close();
  quit


