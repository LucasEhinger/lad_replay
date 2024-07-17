void rootlogon() {
  cout << "rootlogon.C: loading libLAD" << endl;
  // Assumes that the library can be found in (DY)LD_LIBRARY_PATH
  gSystem->Load("libLAD");
}