/**
 * Convert the kinematics file and clusterizer output from AliRoot framework into classes used 
 * neutral pion correlation analysis.
 * 
 * If one simulation run has several jobs this combines the files into one output file.
 * 
 * To be executed in AliRoot environment.
*/

TFile *clusterFile, *fout;
TTree *tree;

AliJHMREvent *event;

void InitCombinedData(TString outputdir);
void GetKinematics(AliRunLoader *runloader);
void GetClusters(TFile *clusterfile, int ievt);
void SaveCombinedOutput();

void CombineSimulationOutput(TString simfolder, TString clusterfolder, TString outputdir, int startFolder, int endFolder)
{

    int nprocessedfiles = 0;
    InitCombinedData(outputdir.Data());
	for (Int_t ifolder = startFolder; ifolder <= endFolder; ifolder++)
	{
		char filename[200];
		sprintf(filename, "%s/%03d/focalClusters.root", clusterfolder.Data(), ifolder);

		clusterFile = TFile::Open(filename);
		cout << "Clusters from: " << clusterFile->GetName() << endl;
		cout << "FOLDER: " << ifolder << endl;

		AliRunLoader *fRunLoader = 0;

		sprintf(filename, "%s/%03d/%s", simfolder.Data(), ifolder, "galice.root");
        cout << filename << endl;

		Long_t id = 0, size = 0, flags = 0, mt = 0;
		if (gSystem->GetPathInfo(filename, &id, &size, &flags, &mt) == 1) {
			cout << "ERROR: FOLDER: " << ifolder << endl;
			continue;
		}

		// Alice run loader
		fRunLoader = AliRunLoader::Open(filename);

		if (!fRunLoader) {
			cout << "ERROR: FOLDER: " << ifolder << endl;
			continue;
		}

		if (!fRunLoader->GetAliRun()) fRunLoader->LoadgAlice();
		if (!fRunLoader->TreeE()) fRunLoader->LoadHeader();
		if (!fRunLoader->TreeK()) fRunLoader->LoadKinematics();

		Int_t nevt = fRunLoader->GetNumberOfEvents();
		for (int ievt = 0; ievt < nevt; ievt++) {
			Int_t ie = fRunLoader->GetEvent(ievt);
			//cout << "Event: " << ievt << " folder " << ifolder << " event " << ievt << endl;

			GetKinematics(fRunLoader);
            GetClusters(clusterFile, ievt);
            SaveCombinedOutput();
		}
		clusterFile->Close();
		fRunLoader->Delete();
        nprocessedfiles++;
	}
    fout->cd();
    tree->Write("", TObject::kOverwrite);

    cout << "\nData combination done for run " << simfolder << endl;
    cout << "Output was saved to " << outputdir << endl;
    cout << "PROCESSING WAS SUCCESFUL FOR " << nprocessedfiles << "/" << endFolder << " JOBS\n" << endl;
}

//************************************************************************************************
//******************************** FUNTCIONS *****************************************************
//************************************************************************************************

void InitCombinedData(TString outputdir)
{
    gSystem->Exec(Form("mkdir %s", outputdir.Data()));
    fout = TFile::Open(Form("%s/output.root", outputdir.Data()), "RECREATE");
    tree = new TTree("Run", "MC and cluster data for FoCal");
    event = new AliJHMREvent();
    event->InitEvent();
    tree->Branch("Events", &event, 32000, 2);
}

void GetKinematics(AliRunLoader *runLoader)
{
    TTree* treeK = runLoader->TreeK();
    AliStack *stack = runLoader->Stack();

    //int npart = stack->GetNtrack();
    int npart = stack->GetNprimary(); // Save only primaries
    for (int ipart = 0; ipart < npart; ipart++) {
        auto track = stack->Particle(ipart);
        int pdg = track->GetPdgCode();
        if (TDatabasePDG::Instance()->GetParticle(pdg)==NULL) {
            continue;
        }
        bool isPrimary = stack->IsPhysicalPrimary(ipart);
        Float_t charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
        Float_t px = track->Px();
        Float_t py = track->Py();
        Float_t pz = track->Pz();
        Float_t e = track->Energy();
        Int_t motherId = track->GetFirstMother();
        Int_t motherPdg = -1;
        if (motherId!=-1) {
            auto motherTrack = stack->Particle(motherId);
            motherPdg = motherTrack->GetPdgCode();
        }
        Float_t mass = TMath::Sqrt(e*e + track->P()*track->P());
        //if (pdg==22) cout << "PDG 22, isPrimary=" << isPrimary << ", mass=" << mass << ", mass(true)=" << track->GetMass() << ", status=" << track->GetStatusCode() << ", mother=" << motherPdg << ", isSecondaryFromWeakDecay=" << stack->IsSecondaryFromWeakDecay(ipart) << ", IsSecondaryFromMaterial=" << stack->IsSecondaryFromMaterial(ipart) << endl;
        event->AddTrack(pdg, px, py, pz, e, charge, ipart, motherId, motherPdg, isPrimary);
    }
}

void GetClusters(TFile *clusterfile, int ievt)
{
    TTree *tClusters = 0;
    if (clusterfile->GetDirectory(Form("Event%i", ievt))) {
        clusterfile->GetDirectory(Form("Event%i", ievt))->GetObject("fTreeR", tClusters);
    } else {
        cout << "Cannot find event " << ievt << " in cluster file " << clusterfile->GetName() << endl;
        clusterfile->ls();
    }

    TBranch *bClusters = tClusters->GetBranch("AliFOCALCluster"); // Branch for final ECAL clusters

    TClonesArray *clustersArray = 0;
    bClusters->SetAddress(&clustersArray);
    bClusters->GetEvent(0);
    for (int iclust = 0; iclust < clustersArray->GetEntries(); iclust++)
    {
        AliFOCALCluster *cluster = (AliFOCALCluster*)clustersArray->At(iclust);
        double clustX = cluster->X();
        double clustY = cluster->Y();
        double clustZ = cluster->Z();
        double clustE = cluster->E();
        double clustEHCAL = cluster->GetHCALEnergy();
        float segE[AliFOCALCluster::N_DEPTH_FIELDS];
        float seedE[AliFOCALCluster::N_DEPTH_FIELDS];
        float width1[AliFOCALCluster::N_DEPTH_FIELDS];
        float width2[AliFOCALCluster::N_DEPTH_FIELDS];
        for (int i=0; i<AliFOCALCluster::N_DEPTH_FIELDS; i++) {
            segE[i] = cluster->GetSegmentEnergy(i);
            seedE[i] = cluster->GetSeedEnergy(i);
            width1[i] = cluster->GetWidth1(i);
            width2[i] = cluster->GetWidth2(i);
        }
        event->AddCluster(clustX, clustY, clustZ, clustE, clustEHCAL, segE, seedE, width1, width2);
    }
}

void SaveCombinedOutput()
{
    tree->Fill();
    event->Clear();
}