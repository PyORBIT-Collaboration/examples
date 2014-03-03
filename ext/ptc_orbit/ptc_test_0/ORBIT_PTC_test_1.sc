///////////////////////////////////////////////////////////////////////
// ORBIT_PTC test srcript - at the end it reads the external PTC script
///////////////////////////////////////////////////////////////////////

//initialize PTC
InitPTC("ptc_data_test_0.txt");


mainHerd = addMacroHerd(1000);



nMacrosPerTurn = 0;
nMaxMacroParticles = 8;
nReals_Macro = 3.00e11/Real(nMaxMacroParticles);
readParts(mainHerd, "bunch_ini.dat", nMaxMacroParticles);

//////////////////////////////
// Start Output:
//////////////////////////////
  OFstream fio("statistics.out", ios::out);
  showStart(fio);
  fio.close();


  cerr << "Start tracking.\n";
  doTurn(1);
  cerr << "Stop tracking.\n";


  OFstream fio1("bunch_final.dat", ios::out);
  dumpPartsGlobal(mainHerd, fio1);
  fio1.close();

  ReadPTCscript("ptc_data_test_0.txt");

quit
