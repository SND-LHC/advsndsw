#include "ChargeDivision.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
#include <sstream>

ChargeDivision::ChargeDivision() {}

void ChargeDivision::ReadPulseShape(std::string PulseFileName)
{
    bool APVPeakMode = true;   // to be included in header configuration file
    if (APVPeakMode == true) {
        ifstream inputFile(PulseFileName);

        if (!inputFile.is_open()) {
            cout << "Error opening the file!" << endl;
        }
        string line;
        string res_find = "resolution:";
        std::vector<double> pulsevalues;
        float res;
        string s;

        while (getline(inputFile, line)) {
            if ((!line.empty()) && (line.substr(0, 1) != "#")) {
                stringstream ss(line);
                if (line.find(res_find) != std::string::npos) {
                    res = stof(line.substr(line.find(res_find) + res_find.length()));   // implement check for
                                                                                        // resolution
                } else {
                    string value;
                    while (getline(ss, value, ' ')) {
                        pulsevalues.push_back(
                            stod(value));   // implement check to see if the max value of the pulse shape is 1
                    }
                }
            }
        }
    }
    // get the vector of beginning to max of the pulse
    // check if pulse correction is needed
}
