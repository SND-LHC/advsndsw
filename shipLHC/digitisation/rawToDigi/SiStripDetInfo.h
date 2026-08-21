#ifndef SNDHLLHC_SISTRIP_DETINFO_H
#define SNDHLLHC_SISTRIP_DETINFO_H

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>

struct DetectorInfo {
    int detid;
    int dcuhardid;
    int fedid;
    int ccuaddress;
    int i2channel;
    int i2caddress;
    int fedchannel;
    int napvpairs;
    int layer;
    int row;
    int column;
};

inline bool CheckDetectorInfo(const DetectorInfo& detinfo) {
    // If there are only 2 apv pairs the strips per module are 512 and wider
    if (detinfo.napvpairs == 2) {
        if ((detinfo.i2caddress != 32) && (detinfo.i2caddress != 36)) {
            std::cerr << "(detinfo.i2caddress != 32) && (detinfo.i2caddress != 36) (" << detinfo.i2caddress << ")\n";
            return false;
        }
    }
    else if (detinfo.napvpairs == 3) {
        if ((detinfo.i2caddress != 32) && (detinfo.i2caddress != 34) && (detinfo.i2caddress != 36)) {
            std::cerr << "(detinfo.i2caddress != 32) && (detinfo.i2caddress != 34) && (detinfo.i2caddress != 36) (" << detinfo.i2caddress << ")\n";
            return false;
        }
    }
    else {
        std::cerr << "(detinfo.napvpairs != 2) && (detinfo.napvpairs != 3) (" << detinfo.napvpairs << ")\n";
        return false;     
    }

    if (detinfo.row > 3) {
        std::cerr << "detinfo.row > 3 (" << detinfo.row << ")\n";
        return false;
    }

    if (detinfo.column > 2) {
        std::cerr << "detinfo.column > 2 (" << detinfo.column << ")\n";
        return false;
    }
    return true;
}

inline std::vector<DetectorInfo> GetDetectorInfo(const std::string& file_name) {
    std::vector<DetectorInfo> records;
    std::ifstream file(file_name);

    if (!file.is_open()) {
        std::cerr << "Could not open file " << file_name << "\n";
        return records;
    }

    std::string line;

    // Skip header line
    std::getline(file, line);

    // Read each CSV row
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value;

        DetectorInfo data;
        try {
            std::getline(ss, value, ',');
            data.detid = std::stoi(value);

            std::getline(ss, value, ',');
            data.dcuhardid = std::stoi(value);

            std::getline(ss, value, ',');
            data.fedid = std::stoi(value);

            std::getline(ss, value, ',');
            data.ccuaddress = std::stoi(value);

            std::getline(ss, value, ',');
            data.i2channel = std::stoi(value);

            std::getline(ss, value, ',');
            data.i2caddress = std::stoi(value);

            std::getline(ss, value, ',');
            data.fedchannel = std::stoi(value);

            std::getline(ss, value, ',');
            data.napvpairs = std::stoi(value);

            std::getline(ss, value, ',');
            data.layer = std::stoi(value);

            std::getline(ss, value, ',');
            data.row = std::stoi(value);

            std::getline(ss, value, ',');
            data.column = std::stoi(value);
        }
        catch (const std::exception& e) {
            std::cerr << "CSV parse error: " << e.what() << "\n";
            continue;
        }
        
        if (CheckDetectorInfo(data)) {
            records.push_back(data);
        }
        else {
            std::cerr << "Invalid detector info found\n";
        }
    }

    return records;
}

inline int GetApvPair(const DetectorInfo& detinfo) {
    // Already checked by CheckDetectorInfo
    if (detinfo.napvpairs == 2) {
        if (detinfo.i2caddress == 32) return 0;
        else if (detinfo.i2caddress == 36) return 1;
        else {
            std::cerr << "Invalid napvpairs: " << detinfo.napvpairs << "\ti2caddress: " << detinfo.i2caddress << "\n";
            return -1;
        }
    }
    else if (detinfo.napvpairs == 3) {
        if (detinfo.i2caddress == 32) return 0;
        else if (detinfo.i2caddress == 34) return 1;
        else if (detinfo.i2caddress == 36) return 2;
        else {
            std::cerr << "Invalid napvpairs: " << detinfo.napvpairs << "\ti2caddress: " << detinfo.i2caddress << "\n";
            return -1;
        }
    }
    else {
        std::cerr << "Invalid napvpairs: " << detinfo.napvpairs << "\ti2caddress: " << detinfo.i2caddress << "\n";
        return -1;
    }
}

#endif