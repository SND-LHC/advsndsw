#pragma once

#include <string>

namespace sndAnalysis{
  class baseCut{
  protected :
    std::string cutName;
  public :
    virtual bool passCut() = 0;
    std::string getName() {return cutName;}
  };
};
