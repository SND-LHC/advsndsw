#include "TBranch.h"
#include "TFile.h"
#include "TTree.h"
#include "Tools/Flux/GSimpleNtpFlux.h"

#include <iostream>
#include <string>
#include <vector>

/*
  Copies auxiliary branches of GENIE default output to the GST ntuple. The gntpc app doesn't copy the auxiliary
  branches...
*/

int main(int argc, char** argv)
{

    if (argc != 3) {
        std::cout << "Two arguments required: path of ghep file AND path of gst file." << std::endl;
        return -1;
    }

    // 500 auxiliary variables should be enough?
    int this_auxint[500];
    double this_auxdbl[500];

    std::vector<std::string>* auxname = new std::vector<std::string>();

    TFile* f_ghep = new TFile(argv[1], "READ");
    TFile* f_gst = new TFile(argv[2], "UPDATE");

    TTree* ghep_meta = (TTree*)f_ghep->Get("meta");
    genie::flux::GSimpleNtpMeta* meta_entry = new genie::flux::GSimpleNtpMeta;
    ghep_meta->SetBranchAddress("meta", &meta_entry);
    ghep_meta->GetEntry(0);

    TTree* gst = (TTree*)f_gst->Get("gst");

    std::vector<TBranch*>* aux_branches = new std::vector<TBranch*>();

    for (int i_aux_var = 0; i_aux_var < meta_entry->auxintname.size(); i_aux_var++) {
        aux_branches->push_back(gst->Branch(meta_entry->auxintname.at(i_aux_var).c_str(),
                                            &(this_auxint[i_aux_var]),
                                            (meta_entry->auxintname.at(i_aux_var) + "/I").c_str()));
        std::cout << "Added branch " << meta_entry->auxintname.at(i_aux_var) << std::endl;
    }

    for (int i_aux_var = 0; i_aux_var < meta_entry->auxdblname.size(); i_aux_var++) {
        aux_branches->push_back(gst->Branch(meta_entry->auxdblname.at(i_aux_var).c_str(),
                                            &(this_auxdbl[i_aux_var]),
                                            (meta_entry->auxdblname.at(i_aux_var) + "/D").c_str()));
        std::cout << "Added branch " << meta_entry->auxdblname.at(i_aux_var) << std::endl;
    }

    TTree* ghep_gtree = (TTree*)f_ghep->Get("gtree");
    genie::flux::GSimpleNtpAux* aux_entry = new genie::flux::GSimpleNtpAux;
    ghep_gtree->SetBranchAddress("aux", &aux_entry);

    for (int i_entry = 0; i_entry < ghep_gtree->GetEntries(); i_entry++) {
        aux_entry->Reset();

        ghep_gtree->GetEntry(i_entry);

        for (int i_auxint = 0; i_auxint < meta_entry->auxintname.size(); i_auxint++)
            this_auxint[i_auxint] = aux_entry->auxint.at(i_auxint);
        for (int i_auxdbl = 0; i_auxdbl < meta_entry->auxdblname.size(); i_auxdbl++)
            this_auxdbl[i_auxdbl] = aux_entry->auxdbl.at(i_auxdbl);

        for (int i_branch = 0; i_branch < aux_branches->size(); i_branch++)
            aux_branches->at(i_branch)->Fill();
    };

    std::cout << "Done copying auxiliary variables. Closing files." << std::endl;

    f_gst->Write();
    f_gst->Close();
    f_ghep->Close();
}
