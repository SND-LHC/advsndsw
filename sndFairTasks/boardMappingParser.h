// This is the C++ version of shipLHC/rawData/boardMappingParser.py
#ifndef BOARDMAP_H
#define BOARDMAP_H 1

#include "nlohmann/json.hpp"

#include <iostream>
#include <map>
#include <string>
#include <tuple>

using namespace std;
using json = nlohmann::json;

tuple<map<string, map<string, map<string, int>>>, map<string, map<string, map<string, string>>>> getBoardMapping(
    json j);

#endif
