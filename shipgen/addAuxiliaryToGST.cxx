#include <vector>
#include <string>

#include "TTree.h"
#include "TBranch.h"

#include "Tools/Flux/GSimpleNtpFlux.h"

/*
  Copies auxiliary branches of GENIE default output to the GST ntuple. The gntpc app doesn't copy the auxiliary branches...
*/


int main(int argc, char** argv){

  if (argc != 3) {
        std::cout << "Two arguments required: path of ghep file AND path of gst file." << std::endl;
	return -1;
  }

  // 500 auxiliary variables should be enough?
  int this_auxint[500];
  double this_auxdbl[500];

  std::vector<std::string> * auxname = new std::vector<std::string>();

  TFile * f_ghep = new TFile(argv[1], "READ");
  TFile * f_gst = new TFile(argv[2], "UPDATE");

  TTree * ghep_meta = (TTree*) f_ghep->Get("meta");
  genie::flux::GSimpleNtpMeta*  meta_entry  = new genie::flux::GSimpleNtpMeta;
  ghep_meta->SetBranchAddress("meta", &meta_entry);
  ghep_meta->GetEntry(0);

  TTree * gst = (TTree*) f_gst->Get("gst");
  
  std::vector<TBranch*> aux_branches = new std::vector<TBranch*>();

  for (int i_name = 0; i_name< meta_entry->auxintname.size(); i_name++){
    aux_branches->push_back(f_gst->Branch(meta_entry->auxintname.at(i_name).c_str(), &(auxint[i_name]), (meta_entry->auxintname.at(i_name)+"/I")+.c_str()));
  }

  for (int i_name = 0; i_name< meta_entry->auxdblname.size(); i_name++){
    aux_branches->push_back(f_gst->Branch(meta_entry->auxdblname.at(i_name).c_str(), &(this_auxdbl[i_name]), (meta_entry->auxdblname.at(i_name)+"/D")+.c_str()));
  }
  
  TTree * ghep_aux = (TTree*) f_ghep->Get("aux");
  genie::flux::GSimpleNtpAux*  aux_entry  = new genie::flux::GSimpleNtpAux; 
  ghep_aux->SetBranchAddress("aux", &aux_entry);

  for (int i_entry = 0; i_entry < f_ghep->GetEntries(); f_ghep++){
    ghep_aux->GetEntry(i_entry);

    for (int i_auxint = 0; i_auxint < meta_entry->auxintname.size(); i_auxint++) this_auxint[i_auxint] = ghep_aux->auxint.at(i_auxint);
    for (int i_auxdbl = 0; i_auxdbl < meta_entry->auxdblname.size(); i_auxdbl++) this_auxdbl[i_auxdbl] = ghep_aux->auxdbl.at(i_auxdbl);

    for (int i_branch = 0; i_branch < aux_branches->size(); i_branch++) aux_branches->at(i_branch)->Fill();
    
  };

  f_gst->Write();
  f_gst->Close();
  f_ghep->Close();
  
}



