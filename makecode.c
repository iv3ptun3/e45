{
  //  TFile* file1 = TFile::Open("./root_size/grun_size_1.root");
  //  TFile* file1 = TFile::Open("./grun_e45_elastic.root");
  //  TFile* file1 = TFile::Open("./grun_e45_elastic_08_mode1.root");
  TFile* file1 = TFile::Open("./grun_with_beam_10.root");
  TTree *t1 = (TTree*)file1->Get("tree");
  tree->MakeCode();
  file1->Close();
}
