#ifndef CONVRAWDATA_H_
#define CONVRAWDATA_H_

#include "FairEventHeader.h"     // for FairEventHeader
#include "FairTask.h"            // for FairTask, InitStatus
#include "MuFilterHit.h"         // for Muon Filter Hit
#include "SNDLHCEventHeader.h"   // for EventHeader
#include "Scifi.h"               // for Scifi detector
#include "sndScifiHit.h"         // for SciFi Hit

#include <TClonesArray.h>
#include <TFile.h>
#include <TMap.h>
#include <iostream>
#include <map>
#include <tuple>

using namespace std;

class ConvRawData : public FairTask
{
  public:
    /** Default constructor **/
    ConvRawData();

    /** Destructor **/
    ~ConvRawData();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);

    /** Update input raw-data file and first-to-process event **/
    void UpdateInput(int n);

  private:
    /** Start time of run **/
    void StartTimeofRun(string path);
    /** Board mapping for Scifi and MuFilter **/
    void DetMapping(string path);
    void checkBoardMapping(string path);
    void debugMapping(string board, int tofpetID, int tofpetChannel);
    /** Define calibration functions **/
    double qdc_calibration(int board_id,
                           int tofpet_id,
                           int channel,
                           int tac,
                           uint16_t v_coarse,
                           uint16_t v_fine,
                           uint16_t tf);
    double qdc_chi2(int board_id, int tofpet_id, int channel, int tac, int TDC);
    double qdc_sat(int board_id, int tofpet_id, int channel, int tac, uint16_t v_fine);
    double
        time_calibration(int board_id, int tofpet_id, int channel, int tac, int64_t t_coarse, uint16_t t_fine, int TDC);
    tuple<double, double, double, double> comb_calibration(int board_id,
                                                           int tofpet_id,
                                                           int channel,
                                                           int tac,
                                                           uint16_t v_coarse,
                                                           uint16_t v_fine,
                                                           int64_t t_coarse,
                                                           uint16_t t_fine,
                                                           double GQDC,
                                                           int TDC);
    // max gain QDC = 3.6
    map<double, pair<double, double>> calibrationReport();

    /** Define some other functions **/
    int channel_func(int tofpet_id, int tofpet_channel, int position);
    /** Read csv data files **/
    void read_csv(string path);
    /** Processing of different raw-data formats **/
    void Process0();
    void Process1();

    /** Data structures to be used in the class **/
    map<int, MuFilterHit*> digiMuFilterStore{};
    map<int, sndScifiHit*> digiSciFiStore{};
    map<vector<int>, map<string, double>> X_qdc{};
    map<vector<int>, map<string, double>> X_tdc{};
    map<string, map<string, map<string, int>>> boardMaps{};
    map<int, map<int, int>> MufiSystem{};   // <board_id_mu, <slot(tofpetID), s>
    map<int, string> slots = {{0, "A"}, {1, "A"}, {2, "B"}, {3, "B"}, {4, "C"}, {5, "C"}, {6, "D"}, {7, "D"}};
    map<int, map<int, int>> TofpetMap{};
    map<string, map<string, map<string, string>>> boardMapsMu{};
    map<string, vector<int>>
        offMap{};   // name is key, vector is first bar, number of sipm channels / bar and direction
    // For time monitoring
    map<string, double> counters = {{"N", 0},
                                    {"event", 0},
                                    {"qdc", 0},
                                    {"tdc", 0},
                                    {"chi2", 0},
                                    {"make", 0},
                                    {"storage", 0},
                                    {"createScifi", 0},
                                    {"createMufi", 0}};

    Scifi* ScifiDet;
    TFile* fOut;
    /** Input data **/
    TTree* fEventTree;
    // Board_x data
    map<string, TTree*> boards{};
    /** Input parameters **/
    int frunNumber;
    int fnStart, fnEvents;
    int fheartBeat;
    int debug, stop, makeCalibration, local, newFormat;
    string fpathCalib, fpathJSON;
    int eventNumber;
    double chi2Max, saturationLimit;
    double runStartUTC;
    /** Filling scheme per run **/
    TMap* FSmap;

    /** Output data **/
    SNDLHCEventHeader* fSNDLHCEventHeader;
    FairEventHeader* fEventHeader;
    TClonesArray* fDigiSciFi;
    TClonesArray* fDigiMuFilter;

    ConvRawData(const ConvRawData&);
    ConvRawData& operator=(const ConvRawData&);

    ClassDef(ConvRawData, 3);
};
#endif /* CONVRAWDATA_H_ */
