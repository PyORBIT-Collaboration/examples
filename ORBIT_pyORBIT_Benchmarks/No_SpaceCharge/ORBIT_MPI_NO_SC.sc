//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// DESCRIPTION
//
// Jeff Holmes 07/28/2010 Space Charge test for Xiyin
//
///////////////////////////////////////////////////////////////////////////

	
OFstream fio("ORBIT_MPI_NO_SC.out", ios::out);

///////////////////////////////////////////////////////////
// Make a synchronous particle and do some initialization:
///////////////////////////////////////////////////////////

  Real TSync = 1.0;     // Kinetic Energy (GeV)
  Real mSync = 1.0;     // Mass (AMU)
  Real charge = 1.0;    // charge number

  addSyncPart(mSync,charge,TSync);

  Real INTENSITY = 1.0e+13;
  nMaxMacroParticles = 10000;
  nReals_Macro = INTENSITY / nMaxMacroParticles;
  
  mainHerd = addMacroHerd(nMaxMacroParticles);
  readParts(mainHerd, "orbit_mpi_bunch_input.dat", nMaxMacroParticles);

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

  buildTPlatticeNew("LATTICES/Q_0p125.TP",
                    "Ring",
                    nstepTPD,
                    nstepTPM, fringeM,
                    nstepTPQ, fringeQ,
                    nstepTPB, fringeB,
                    nstepTPS,
                    nstepTPK);

  //Ring::gammaTrans = 1.e+10;
  //Ring::nuX=0.125;
  //Ring::nuY=0.125;

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

  OFstream fio98("orbit_mpi_bunch_output.dat", ios::out);
  dumpPartsGlobal(mainHerd, fio98);
  fio98.close();
	
  fio.close();

  quit
