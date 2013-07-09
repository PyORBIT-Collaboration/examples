
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
  runName = "Run_Benchmark";
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

  const Integer nstepTPD = 4;
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


buildTPlatticeNew("MAD_Lattice/RealInjection/SNSring_pyOrbitBenchmark_nokick.TP","Ring",
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
// Make a Herd
////////////////////////////////////////////

readParts(mainHerd, "Bunches/controlbunch_600_ORBIT.dat", 1000);

addFracTuneNode("FTune", 943, "FracTunes");
//////////////////////////////

activateFracTuneNodes();


addStatLatNodeSet("StatLats", "all");
addMomentNodeSet(4, "Moments", "all"); 
activateStatLatNodes();
activateMomentNodes();
///////////////////////////////////////////////////
// do some turns, and dump particles for later plots:
/////////////////////////////////////////////////////

//useSimpleTuneCalc = 1;  


if( nRank == 0 ) cerr << "Start Tracking\n";
Real et;
timerOn();
showNodes(cerr);

doTurn(2);
OFstream fio07("FracTunes", ios::out);
dumpFracTunes(fio07);
fio07.close();
deactivateFracTuneNodes();
OFstream fio12("Bm_Parts2", ios::out);

dumpParts(mainHerd, fio12);
fio12.close();


  

  quit


