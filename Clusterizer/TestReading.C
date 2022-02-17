#include "TClonesArray.h"
#include "AliJHMREvent.h"

TFile *fin;
TTree *tree;
AliJHMREvent *event;
TBranch *brtrack;
TBranch *brcluster;
TClonesArray *tracks = nullptr;
TClonesArray *clusters = nullptr;

int LoadInput(TString infile);

void TestReading(TString infile="results/output.root")
{
    //gSystem->Load("AliJHMREvent_cxx.so");
    if (!LoadInput(infile)) return;
    int nev = tree->GetEntries();
    cout << "Number of events : " << nev << endl;
    for (int iev=0; iev<nev; iev++) {
    cout << "\nEvent " << iev << endl;
        tree->GetEntry(iev);
        // Go through tracks
        int ntrack = tracks->GetEntries();
        cout << "\tNumber of tracks : " << ntrack << endl;
        for (int itrack=0; itrack<ntrack; itrack++) {
            Track *track = (Track*)tracks->At(itrack);
            //cout << "\t\t" << track->GetEta() << endl;
        }

        // Go through clusters
        int ncluster = clusters->GetEntries();
        cout << "\tNumber of clusters : " << ncluster << endl;
        for (int icluster=0; icluster<ncluster; icluster++) {
            Cluster *cluster = (Cluster*)clusters->At(icluster);
            cout << "\t\t" << cluster->GetEta() << endl;
        }
    }
}

//******************************************************************************
//******************************************************************************

int LoadInput(TString infile)
{
    fin = TFile::Open(infile.Data(), "READ");
    if (!fin) {
        cout << "File " << infile << " not found!" << endl;
        return 0;
    }

    tree = (TTree*)fin->Get("Run");
    tree->SetBranchAddress("Events", &event);
    brtrack = tree->GetBranch("Events.fTracks");
    brtrack->SetAddress(&tracks);
    brcluster = tree->GetBranch("Events.fClusters");
    brcluster->SetAddress(&clusters);
    return 1;
}
