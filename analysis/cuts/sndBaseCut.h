#pragma once

#include <string>

namespace sndAnalysis{
  class baseCut{
  protected :
    std::string cutName;
    std::string shortName;
    std::vector<int> nbins;
    std::vector<double> range_start;
    std::vector<double> range_end;
    std::vector<double> plot_var;
  public :
    virtual bool passCut() = 0;
    std::string getName() {return cutName;}

    // For histograms
    std::string getShortName() {return shortName;}
    std::vector<int> getNbins() {return nbins;}
    std::vector<double> getRangeStart() {return range_start;}
    std::vector<double> getRangeEnd() {return range_end;}
    std::vector<double> getPlotVar() {return plot_var;}
  };
};
