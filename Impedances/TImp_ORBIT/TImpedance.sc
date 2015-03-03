
//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// DESCRIPTION
//
// Test Transverse Impedance. Compare to pyORBIT.
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

  Real TSync  = 1.0;   // Kinetic Energy (GeV)
  Real mSync  = 1.0;   // Mass (AMU)
  Real charge = 1.0;   // charge number

  addSyncPart(mSync, charge, TSync);

  mainHerd = addMacroHerd(10000);

//////////////////////////////
// Make  a Ring
//////////////////////////////

  const Integer nstepTPD = 1;

  Ring::nuX    =  6.2;
  Ring::betaX  = 10.0;
  Ring::alphaX =  0.0;
  Ring::nuY    =  6.2;
  Ring::betaY  = 10.0;
  Ring::alphaY =  0.0;

  Real length = 248.0;
  addTPD("Drift", 2, Ring::betaX,   Ring::betaY,
                     Ring::alphaX, Ring::alphaY,
                     0, 0, length, "Drift", nstepTPD);

//////////////////////////////
// Add Beam:
//////////////////////////////

  nMaxMacroParticles = 10000;
  nReals_Macro = 1.0e+16 / nMaxMacroParticles;
  readParts(mainHerd, "Bm_KV_Uniform_10000", nMaxMacroParticles);

  nMacrosPerTurn = 0;

//////////////////////////////
// Add a Transverse Impedance Node
//////////////////////////////

  Integer ilist;

  nLBinsTImp = 64;
  ComplexVector ZTImpedplus(nLBinsTImp/2), ZTImpedminus(nLBinsTImp/2);
  ComplexVector ZTImpedplusX(nLBinsTImp/2), ZTImpedminusX(nLBinsTImp/2);
  ComplexVector ZTImpedplusY(nLBinsTImp/2), ZTImpedminusY(nLBinsTImp/2);

  for(ilist = 1; ilist <= nLBinsTImp/2; ilist++)
  {
    ZTImpedplusX(ilist)  = Complex(0.,0.);
    ZTImpedminusX(ilist) = Complex(0.,0.);
    ZTImpedplusY(ilist)  = Complex(0.,0.);
    ZTImpedminusY(ilist) = Complex(0.,0.);
  }

  IFstream  fio03("TImpedance.dat", ios::in);

  Integer dummy;
  Real PlusRe, PlusIm, MinusRe, MinusIm;
  for(ilist = 1; ilist <= nLBinsTImp/2; ilist++)
  {
    fio03 >> dummy
          >> PlusRe
          >> PlusIm
          >> MinusRe
          >> MinusIm;

    PlusRe  *= 100.0;
    PlusIm  *= 100.0;
    MinusRe *= 100.0;
    MinusIm *= 100.0;

//    ZTImpedplusX(ilist)  = Complex(PlusRe,  PlusIm);
//    ZTImpedminusX(ilist) = Complex(MinusRe, MinusIm);
    ZTImpedplusY(ilist)  = Complex(PlusRe,  PlusIm);
    ZTImpedminusY(ilist) = Complex(MinusRe, MinusIm);
  }
  fio03.close();

  Real b_a = 10.0/3.0;
  Integer useAvg = 0;
//  useXdimension = 1;
//  useYdimension = 0;
  useXdimension = 0;
  useYdimension = 1;

  Integer nMacroTImpMin = 1;
  addFFTTImpedance("TImp", 1, ZTImpedplusX, ZTImpedminusX,
                              ZTImpedplusY, ZTImpedminusY,
                              b_a, useAvg, nMacroTImpMin);

  if( nRank == 0 ) cerr << "Start Tracking\n";
  Real et;
  timerOn();
  showNodes(cerr);

  turnToNode(mainHerd,1);
  OFstream fio12("Bm_Parts", ios::out);
  dumpParts(mainHerd, fio12);
  OFstream fio13("Lost_Parts", ios::out);
  dumpLostParts(mainHerd, fio13);

  fio12.close();
  fio13.close();
  fio.close();
  quit


