
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

  Real length = 248.0;
  addTPD("Drift", 2, 0, 0, 0, 0, 0, 0, length, "Drift", nstepTPD);

//////////////////////////////
// Add Beam:
//////////////////////////////

  nMaxMacroParticles = 10000;
  nReals_Macro = 1.0e+10 / nMaxMacroParticles;
  readParts(mainHerd, "Bm_KV_Uniform_10000", nMaxMacroParticles);

  nMacrosPerTurn = 0;

//////////////////////////////
// Add Long SC Node
//////////////////////////////

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
   
    //Multiply by an order of magnitude to make the effect bigger for the benchmark.
    realp = 10.0 * realp;
    imagp = 10.0 * imagp;
    ZImped(harm) = Complex(realp, imagp);
  }
  Integer maxModes = 32;

  nLongBins = 128;
  Real b_a = 10.0/3.0;
  Integer useAvg = 0;
  Integer nMacroLSCMin = 1;

  Integer useSpaceCharge = 1;
  addFFTLSpaceCharge("LSC1", 1, ZImped, b_a, useAvg,
                    nMacroLSCMin, useSpaceCharge);


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


