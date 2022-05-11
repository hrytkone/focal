//******************************************************************************
//      CONSTANTS & PARAMETERS
//******************************************************************************

const TString detectorGeometryFile      = "geometry.txt";
const TString clusterizerParametersFile = "parameters.txt";
const TString clusteringOutputFileDir   = "results";
const TString clusteringParametersTag   = "pi0";

bool saveFigures = false;

// Whether should clusterizer provide debug messages
bool clusteringDebugMode = false;

// Whether clusterizer should run in 'calibration mode'
// in which raw signal from detector is saved. Used
// for determining calibration function and its parameters
bool clusteringCalibrationRun = false;

// Default parameters for clusterizer
// LKegacy leftover; should be overridden by 'parameters.txt'
Float_t dist = 4.;
Float_t eThreshold = 100000;

// Noise settings
Float_t pixelNoiseProb = 0;
Float_t padNoiseSigma = 0;

Float_t logWeight = 4.5; // w_0 for log weights
Float_t powWeight = 2.2; // p for power weights

// Which segments to use as 'seeds'
// for creating clusters from sub-clusters
Int_t  COARSE_SEEDS[3] = {2,4,0};
Int_t  FINE_SEEDS[2]   = {1,3};
Int_t  ALL_SEEDS [5]   = {2,4,0,1,3};

// Which segments paramteres should be used for calibration
Int_t COARSE_CALIB = 0;
Int_t FINE_CALIB = 1;

// Energy weights for individual segments
Float_t  ENERGY_WEIGHTS[6] = {1,1,1,1,1,1};

//******************************************************************************
//      OBJECTS
//******************************************************************************

AliFOCALClusterizerv2 *fFOCALCluster;
AliFOCALGeometry * geometry;
AliRunLoader * fRunLoader;
AliFOCAL *fFOCAL;
AliFOCALLoader *fFOCALLoader;
TDirectory * actualDir;

// Cluster maps
ObjectMap * * clusterMap;
ObjectMap * coarseClusterMap;
ObjectMap * fineClusterMap;
ObjectMap * finalClusterMap;

Int_t nSegments;
Int_t firstFineSeed;
Int_t firstCoarseSeed;

// For saving combined output
TFile *outfile;
TTree *tree;
TClonesArray *clusterArray = 0;

// Incident particle info
Int_t pdgCode[4];
Float_t e[4];
Float_t pt[4];
Float_t phi[4];
Float_t theta[4];
Float_t Vx[4], Vy[4], Vz[4];
Float_t cVx[4], cVy[4], cVz[4];
Int_t conv_flag[4];

//******************************************************************************
//      FUNCTIONS
//******************************************************************************

int CheckFile(TString filename);
void InitCombinedData();
void SaveCombinedOutput(); // Save kinematics & clusterizer results to same file
void GetKinematics(int nevtot);
void InitClusterizer();
void CreateClusterMaps();
void LoadClusterizerHit(TString inputfile);
void LoadEvent(TString simFolder, int ifolder, int ievt);
void Clusterize(int nfolder, int ievt);
Float_t Calibration(Float_t signal, Int_t segment, Float_t * pars, bool calibrationRun, bool debug);
