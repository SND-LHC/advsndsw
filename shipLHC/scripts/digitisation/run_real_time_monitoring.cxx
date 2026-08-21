#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <memory>
#include <string>
#include <cstdint>
#include <numeric>
#include "ROOT/RDataFrame.hxx"
#include <ROOT/RVec.hxx>
#include "TCanvas.h"
#include "TH1D.h"
#include "TFile.h"
#include "TDirectory.h"
#include "SiStripIO.h"
#include "AdvHit.h"
#include "SiStripRawToDigi.h"
#include "SiStripDigiClustering.h"
#include "SiStripDetInfo.h"
#include "SiStripPosition.h"
#include "SiStripHardwareConstants.h"

template <typename T, typename Ret, typename Method>
auto ExtractRVec(Method method)
{
    return [method](const std::vector<T>& v) {
        ROOT::RVec<Ret> out;
        out.reserve(v.size());
        for (const auto& obj : v) {
            out.emplace_back((obj.*method)());
        }
        return out;
    };
}

int main(int argc, char* argv[]){
    if (argc != 6) {
        std::cerr << "5 arguments expected but " << argc - 1 << " provided\n";
        std::cerr << "Usage: " << argv[0] << " <input_root_file> <detector_info> <geometry_file> <output_histos_root_file> <n_treads>\n";
        return 1;
    }
    std::string input_root_file(argv[1]);
    std::string detector_info_path(argv[2]);
    std::string geometry_file(argv[3]);
    std::string output_root_file(argv[4]);

    ROOT::EnableImplicitMT(std::stoi(argv[5]));
    std::cout << "Requested threads: " << std::stoi(argv[5]) << '\n';
    std::cout << "ROOT pool size:    " << ROOT::GetThreadPoolSize() << '\n';

    TGeoManager::Import(geometry_file.c_str());
    gGeoManager->SetMaxThreads(ROOT::GetThreadPoolSize());

    auto df = ROOT::RDataFrame("Events", input_root_file);
    // Perform digitization + clustering
    auto df2 = df.Define("ClusteringProducts", SiStripDigiClustering(), {"FedChannelDigis"})
        .Define("FedChannelDigis_Clust", [](const SiStripClusteringProducts& p) { return p.digis; }, {"ClusteringProducts"})
        .Define("Cluster", [](const SiStripClusteringProducts& p) { return p.clusters; }, {"ClusteringProducts"});
    // Prepare columns for histos
    auto df3 = df2.Define("adc", ExtractRVec<AdvHit, uint16_t>(&AdvHit::GetSignal), {"FedChannelDigis_Clust"})
        .Define("strip", ExtractRVec<AdvHit, uint16_t>(&AdvHit::GetStrip), {"FedChannelDigis_Clust"})
        .Define("layer", ExtractRVec<AdvHit, int>(&AdvHit::GetLayer), {"FedChannelDigis_Clust"})
        .Define("row", ExtractRVec<AdvHit, int>(&AdvHit::GetRow), {"FedChannelDigis_Clust"})
        .Define("column", ExtractRVec<AdvHit, int>(&AdvHit::GetColumn), {"FedChannelDigis_Clust"})
        .Define("detector_id", ExtractRVec<AdvHit, uint32_t>(&AdvHit::GetDetectorID), {"FedChannelDigis_Clust"})
        .Define("is_vertical", ExtractRVec<AdvHit, bool>(&AdvHit::IsVertical), {"FedChannelDigis_Clust"})
        .Define("position", [](const ROOT::VecOps::RVec<uint32_t>& detids) {return ROOT::VecOps::Map(detids, GetSiStripPosition);}, {"detector_id"})
        .Define("x", [](const ROOT::VecOps::RVec<ROOT::Math::XYZPoint>& points) {return ROOT::VecOps::Map(points, [](const auto& p) { return p.X(); });}, {"position"})
        .Define("y", [](const ROOT::VecOps::RVec<ROOT::Math::XYZPoint>& points) {return ROOT::VecOps::Map(points, [](const auto& p) { return p.Y(); });}, {"position"})
        .Define("z", [](const ROOT::VecOps::RVec<ROOT::Math::XYZPoint>& points) {return ROOT::VecOps::Map(points, [](const auto& p) { return p.Z(); });}, {"position"})
        .Define("x_vertical", "x[is_vertical]")
        .Define("y_not_vertical", "y[is_vertical == false]")
        .Define("z_vertical", "z[is_vertical]")
        .Define("z_not_vertical", "z[is_vertical == false]")
        .Define("cluster_adc", ExtractRVec<SiStripCluster, uint32_t>(&SiStripCluster::GetSignal), {"Cluster"})
        .Define("cluster_size", ExtractRVec<SiStripCluster, size_t>(&SiStripCluster::GetSize), {"Cluster"})
        .Define("cluster_layer", ExtractRVec<SiStripCluster, int>(&SiStripCluster::GetLayer), {"Cluster"})
        .Define("cluster_row", ExtractRVec<SiStripCluster, int>(&SiStripCluster::GetRow), {"Cluster"})
        .Define("cluster_column", ExtractRVec<SiStripCluster, int>(&SiStripCluster::GetColumn), {"Cluster"})
        .Define("cluster_detector_id", ExtractRVec<SiStripCluster, uint32_t>(&SiStripCluster::GetDetectorID), {"Cluster"})
        .Define("cluster_is_vertical", ExtractRVec<SiStripCluster, bool>(&SiStripCluster::IsVertical), {"Cluster"})
        .Define("cluster_position", [](const ROOT::VecOps::RVec<uint32_t>& detids) {return ROOT::VecOps::Map(detids, GetSiStripPosition);}, {"cluster_detector_id"})
        .Define("cluster_x", [](const ROOT::VecOps::RVec<ROOT::Math::XYZPoint>& points) {return ROOT::VecOps::Map(points, [](const auto& p) { return p.X(); });}, {"cluster_position"})
        .Define("cluster_y", [](const ROOT::VecOps::RVec<ROOT::Math::XYZPoint>& points) {return ROOT::VecOps::Map(points, [](const auto& p) { return p.Y(); });}, {"cluster_position"})
        .Define("cluster_z", [](const ROOT::VecOps::RVec<ROOT::Math::XYZPoint>& points) {return ROOT::VecOps::Map(points, [](const auto& p) { return p.Z(); });}, {"cluster_position"});

    // Retrieve detinfo values to generate histograms dynamically
    const std::vector<DetectorInfo> detinfo = GetDetectorInfo(detector_info_path);
    const int max_row = std::max_element(detinfo.begin(), detinfo.end(), [](const DetectorInfo& x, const DetectorInfo& y) {return x.row < y.row;})->row;
    const int max_column = std::max_element(detinfo.begin(), detinfo.end(), [](const DetectorInfo& x, const DetectorInfo& y) {return x.column < y.column;})->column;

    std::map<int, std::set<std::pair<int,int>>> layer_map;
    for (const auto& d : detinfo) {
        layer_map[d.layer].insert({d.row, d.column});
    }

    std::vector<ROOT::RDF::RResultPtr<TH1D>> histos_1d_global;
    std::vector<ROOT::RDF::RResultPtr<TH2D>> histos_2d_global;

    std::map<int, std::vector<ROOT::RDF::RResultPtr<TH1D>>> histos_1d;
    std::map<int, std::vector<ROOT::RDF::RResultPtr<TH2D>>> histos_2d;
    std::vector<ROOT::RDF::RResultPtr<TProfile>> profiles;

    const int max_layer = layer_map.rbegin()->first;
    constexpr double x_min = -30.0;
    constexpr double x_max = 0.0;
    constexpr double y_min = 15.0;
    constexpr double y_max = 45.0;
    constexpr double z_min = 165.0;
    constexpr double z_max = 195.0;

    auto df4 = df3.Define("layer_vec", [max_layer]() {
            ROOT::RVec<int> l_vec(max_layer + 1);
            std::iota(l_vec.begin(), l_vec.end(), 0);
            return l_vec;
        }).Define("saturated_percentage_vec", [max_layer](const ROOT::RVec<uint16_t>& adc, const ROOT::RVec<int>& layer) {
            ROOT::RVec<float> sat_vec(max_layer + 1);
            for (int i{0}; i < max_layer + 1; ++i) {
                const auto& adc_layer = adc[layer == i];
                sat_vec[i] = adc_layer.empty() ? 0.f : 100.f * adc_layer[adc_layer > 253].size() / adc_layer.size();
            }
            return sat_vec;
        }, {"adc", "layer"})
        .Define("adc_vec", [max_layer](const ROOT::RVec<uint16_t>& adc, const ROOT::RVec<int>& layer) {
            ROOT::RVec<uint32_t> adc_vec(max_layer + 1, 0);
            for (int i{0}; i < max_layer + 1; ++i) {
                adc_vec[i] = ROOT::VecOps::Sum(adc[layer == i]);
            }
            return adc_vec;
        }, {"adc", "layer"});


    histos_1d_global.emplace_back(df3.Histo1D<ROOT::RVec<double>>({"x", "x;x [cm];Entries", 1000, x_min, x_max}, "x_vertical"));
    histos_1d_global.emplace_back(df3.Histo1D<ROOT::RVec<double>>({"y", "y;y [cm];Entries", 1000, y_min, y_max}, "y_not_vertical"));
    histos_1d_global.emplace_back(df3.Histo1D<ROOT::RVec<double>>({"z", "z;z [cm];Entries", 1000, z_min, z_max}, "z"));
    histos_2d_global.emplace_back(df4.Histo2D<ROOT::RVec<double>, ROOT::RVec<double>>({"x vs z", "x vs z;z [cm];x [cm]", 1000, z_min, z_max, 1000, x_min, x_max}, "z_vertical", "x_vertical"));
    histos_2d_global.emplace_back(df4.Histo2D<ROOT::RVec<double>, ROOT::RVec<double>>({"y vs z", "y vs z;z [cm];y [cm]", 1000, z_min, z_max, 1000, y_min, y_max}, "z_not_vertical", "y_not_vertical"));

    histos_1d_global.emplace_back(df3.Histo1D<ROOT::RVec<int>>({"layer", "layer;layer;Entries", max_layer + 1, 0, static_cast<double>(max_layer + 1)}, "layer"));
    histos_1d_global.emplace_back(df3.Histo1D<ROOT::RVec<uint32_t>>({"det_id (advsndsw)", "det_id (advsndsw);det_id;Entries", 18000, 0, 180000}, "detector_id"));
    histos_2d_global.emplace_back(df4.Histo2D<ROOT::RVec<int>, ROOT::RVec<float>>({"saturated percentage (adc > 253) vs layer", "saturated percentage (adc > 253) vs layer;layer;saturated %", max_layer + 1, 0, static_cast<double>(max_layer + 1), 101, 0, 101}, "layer_vec", "saturated_percentage_vec"));
    histos_2d_global.emplace_back(df4.Histo2D<ROOT::RVec<int>, ROOT::RVec<uint32_t>>({"adc sum vs layer", "adc sum vs layer;layer;adc sum", max_layer + 1, 0, static_cast<double>(max_layer + 1), 1000, 0, 30000}, "layer_vec", "adc_vec"));

    ROOT::RDF::TProfile1DModel saturation_model(
        "saturated percentage (adc > 253) vs layer (profiled)",
        "saturated percentage (adc > 253) vs layer;layer;saturated %",
        max_layer + 1, 0, max_layer + 1
    );
    profiles.emplace_back(df4.Profile1D<float>(saturation_model, "layer_vec", "saturated_percentage_vec"));

    ROOT::RDF::TProfile1DModel adc_model(
        "adc sum vs layer (profiled)",
        "adc sum vs layer;layer;adc sum",
        max_layer + 1, 0, max_layer + 1
    );
    profiles.emplace_back(df4.Profile1D<float>(adc_model, "layer_vec", "adc_vec"));

    auto compute_saturated_percentage = [](const ROOT::RVec<uint16_t>& adc) {return adc.empty() ? 0.f : 100.f * adc[adc > 253].size() / adc.size();};

    for (const auto& [layer, modules] : layer_map) {

        std::string h_row_vs_column_name = Form("row vs column Layer %d", layer);
        

        std::string select_row = Form("row[(layer == %d)]", layer);
        std::string select_column = Form("column[(layer == %d)]", layer);
        
        

        if (layer % 2 == 0) {
            std::string h_x_name = Form("x Layer %d", layer);
            std::string h_cluster_x_name = Form("cluster x Layer %d", layer);
            std::string select_x = Form("x[(layer == %d) && is_vertical]", layer);
            std::string select_cluster_x = Form("cluster_x[(cluster_layer == %d) && cluster_is_vertical]", layer);
            auto df_layer = df3.Define("row_layer", select_row).Define("column_layer", select_column).Define("x_vertical_layer", select_x).Define("cluster_x_vertical_layer", select_cluster_x);
            histos_1d[layer].emplace_back(df_layer.Histo1D<ROOT::RVec<double>>({h_x_name.c_str(), (h_x_name + std::string(";x [cm];Entries")).c_str(), 1000, x_min, x_max}, "x_vertical_layer"));
            histos_1d[layer].emplace_back(df_layer.Histo1D<ROOT::RVec<double>>({h_cluster_x_name.c_str(), (h_cluster_x_name + std::string(";x [cm];Entries")).c_str(), 1000, x_min, x_max}, "cluster_x_vertical_layer"));
            histos_2d[layer].emplace_back(df_layer.Histo2D<ROOT::RVec<int>, ROOT::RVec<int>>({h_row_vs_column_name.c_str(), (h_row_vs_column_name + std::string(";column;row;Entries")).c_str(), max_column + 1, 0, static_cast<double>(max_column + 1), max_row + 1, 0, static_cast<double>(max_row + 1)}, "column_layer", "row_layer"));
        }
        else {
            std::string h_y_name = Form("y Layer %d", layer);
            std::string h_cluster_y_name = Form("cluster y Layer %d", layer);
            std::string select_y = Form("y[(layer == %d) && (is_vertical == false)]", layer);
            std::string select_cluster_y = Form("cluster_y[(cluster_layer == %d) && (cluster_is_vertical == false)]", layer);
            auto df_layer = df3.Define("row_layer", select_row).Define("column_layer", select_column).Define("y_not_vertical_layer", select_y).Define("cluster_y_not_vertical_layer", select_cluster_y);
            histos_1d[layer].emplace_back(df_layer.Histo1D<ROOT::RVec<double>>({h_y_name.c_str(), (h_y_name + std::string(";y [cm];Entries")).c_str(), 1000, y_min, y_max}, "y_not_vertical_layer"));
            histos_1d[layer].emplace_back(df_layer.Histo1D<ROOT::RVec<double>>({h_cluster_y_name.c_str(), (h_cluster_y_name + std::string(";y [cm];Entries")).c_str(), 1000, y_min, y_max}, "cluster_y_not_vertical_layer"));
            histos_2d[layer].emplace_back(df_layer.Histo2D<ROOT::RVec<int>, ROOT::RVec<int>>({h_row_vs_column_name.c_str(), (h_row_vs_column_name + std::string(";column;row;Entries")).c_str(), max_column + 1, 0, static_cast<double>(max_column + 1), max_row + 1, 0, static_cast<double>(max_row + 1)}, "column_layer", "row_layer"));
        }
       
        for (const auto& [row, col] : modules) {

            std::string h_adc_name = Form("adc Layer %d Row %d Column %d", layer, row, col);
            std::string h_strip_name = Form("strip Layer %d Row %d Column %d", layer, row, col);
            std::string h_nhits_name = Form("nhits Layer %d Row %d Column %d", layer, row, col);
            std::string h_saturated_percentage_name = Form("saturated percentage (adc > 253) Layer %d Row %d Column %d", layer, row, col);
            std::string h_adc_vs_strip_name = Form("adc vs strip Layer %d Row %d Column %d", layer, row, col);
            std::string h_cluster_adc_name = Form("cluster adc Layer %d Row %d Column %d", layer, row, col);
            std::string h_cluster_size_name = Form("cluster size Layer %d Row %d Column %d", layer, row, col);
            std::string h_n_clusters_name = Form("n clusters Layer %d Row %d Column %d", layer, row, col);
            std::string h_cluster_adc_vs_size_name = Form("cluster adc vs size Layer %d Row %d Column %d", layer, row, col);

            std::string select_adc = Form("adc[(layer == %d) && (row == %d) && (column == %d)]", layer, row, col);
            std::string select_strip = Form("strip[(layer == %d) && (row == %d) && (column == %d)]", layer, row, col);
            std::string select_cluster_adc = Form("cluster_adc[(cluster_layer == %d) && (cluster_row == %d) && (cluster_column == %d)]", layer, row, col);
            std::string select_cluster_size = Form("cluster_size[(cluster_layer == %d) && (cluster_row == %d) && (cluster_column == %d)]", layer, row, col);

            auto df_module = df3.Define("adc_module", select_adc).Define("sistrip_module", select_strip).Define("nhits_module", "adc_module.size()").Define("saturated_percentage_module", compute_saturated_percentage, {"adc_module"})
                .Define("cluster_adc_module", select_cluster_adc).Define("cluster_size_module", select_cluster_size).Define("n_clusters_module", "cluster_adc_module.size()");

            histos_1d[layer].emplace_back(df_module.Histo1D<ROOT::RVec<uint16_t>>({h_adc_name.c_str(), (h_adc_name + std::string(";adc;Entries")).c_str(), 256, 0, 256}, "adc_module"));
            histos_1d[layer].emplace_back(df_module.Histo1D<ROOT::RVec<uint16_t>>({h_strip_name.c_str(), (h_strip_name + std::string(";strip;Entries")).c_str(), MAX_SISTRIPS_PER_MODULE, 0, static_cast<double>(MAX_SISTRIPS_PER_MODULE)}, "sistrip_module"));
            histos_1d[layer].emplace_back(df_module.Histo1D<std::size_t>({h_nhits_name.c_str(), (h_nhits_name + std::string(";nhits;Entries")).c_str(), MAX_SISTRIPS_PER_MODULE, 0, static_cast<double>(MAX_SISTRIPS_PER_MODULE)}, "nhits_module"));
            histos_1d[layer].emplace_back(df_module.Histo1D<ROOT::RVec<uint32_t>>({h_cluster_adc_name.c_str(), (h_cluster_adc_name + std::string(";adc;Entries")).c_str(), 3000, 0, 3000}, "cluster_adc_module"));
            histos_1d[layer].emplace_back(df_module.Histo1D<ROOT::RVec<size_t>>({h_cluster_size_name.c_str(), (h_cluster_size_name + std::string(";size;Entries")).c_str(), 100, 0, 100}, "cluster_size_module"));
            histos_1d[layer].emplace_back(df_module.Histo1D<std::size_t>({h_n_clusters_name.c_str(), (h_n_clusters_name + std::string(";n clusters;Entries")).c_str(), MAX_SISTRIPS_PER_MODULE / 8, 0, static_cast<double>(MAX_SISTRIPS_PER_MODULE / 8)}, "n_clusters_module"));
            histos_1d[layer].emplace_back(df_module.Histo1D<float>({h_saturated_percentage_name.c_str(), (h_saturated_percentage_name + std::string(";saturated %;Entries")).c_str(), 101, 0, 101}, "saturated_percentage_module"));
            histos_2d[layer].emplace_back(df_module.Histo2D<ROOT::RVec<uint16_t>, ROOT::RVec<uint16_t>>({h_adc_vs_strip_name.c_str(), (h_adc_vs_strip_name + std::string(";strip;adc;Entries")).c_str(),  MAX_SISTRIPS_PER_MODULE, 0, static_cast<double>(MAX_SISTRIPS_PER_MODULE), 256, 0, 256}, "sistrip_module", "adc_module"));
            histos_2d[layer].emplace_back(df_module.Histo2D<ROOT::RVec<size_t>, ROOT::RVec<uint32_t>>({h_cluster_adc_vs_size_name.c_str(), (h_cluster_adc_vs_size_name + std::string(";size;adc;Entries")).c_str(), 100, 0, 100, 3000, 0, 3000}, "cluster_size_module", "cluster_adc_module"));

            if (layer % 2 == 0) {
                std::string h_x_vs_strip_name = Form("x vs strip Layer %d Row %d Column %d", layer, row, col);
                auto df_debug = df_module.Define("x_module", Form("x[(layer == %d) && (row == %d) && (column == %d) && is_vertical]", layer, row, col));
                histos_2d[layer].emplace_back(df_debug.Histo2D<ROOT::RVec<uint16_t>, ROOT::RVec<double>>({h_x_vs_strip_name.c_str(), (h_x_vs_strip_name + std::string(";strip;x [cm];Entries")).c_str(), MAX_SISTRIPS_PER_MODULE, 0, static_cast<double>(MAX_SISTRIPS_PER_MODULE), 1000, x_min, x_max}, "sistrip_module", "x_module"));
            }
            else {
                std::string h_y_vs_strip_name = Form("y vs strip Layer %d Row %d Column %d", layer, row, col);
                auto df_debug = df_module.Define("y_module", Form("y[(layer == %d) && (row == %d) && (column == %d) && (is_vertical == false)]", layer, row, col));
                histos_2d[layer].emplace_back(df_debug.Histo2D<ROOT::RVec<uint16_t>, ROOT::RVec<double>>({h_y_vs_strip_name.c_str(), (h_y_vs_strip_name + std::string(";strip;y [cm];Entries")).c_str(), MAX_SISTRIPS_PER_MODULE, 0, static_cast<double>(MAX_SISTRIPS_PER_MODULE), 1000, y_min, y_max}, "sistrip_module", "y_module"));
            }
        }
    }

    // Write histos
    TFile output_file(output_root_file.c_str(), "RECREATE");

    for (auto& h : histos_1d_global) {
        h->SetFillColor(kAzure - 9);
        h->Write();
    }

    // Trigger evaluation and get raw pointers  
    auto hX = histos_1d_global[0].GetPtr();  // "x" from x_vertical  
    auto hY = histos_1d_global[1].GetPtr();  // "y" from y_not_vertical  
    
    TH2D* h2 = new TH2D("beam spot", "beam spot;x [cm];y [cm]",  
        200, hX->GetXaxis()->GetXmin(), hX->GetXaxis()->GetXmax(),  // reduced from 1000  
        200, hY->GetXaxis()->GetXmin(), hY->GetXaxis()->GetXmax()  
    );  
    
    for (int i = 1; i <= h2->GetNbinsX(); i++) {  
        for (int j = 1; j <= h2->GetNbinsY(); j++) {  
            // Interpolate since we're changing bin count from 1000 to 200  
            double x = h2->GetXaxis()->GetBinCenter(i);  
            double y = h2->GetYaxis()->GetBinCenter(j);  
            h2->SetBinContent(i, j, hX->Interpolate(x) * hY->Interpolate(y));  
        }  
    }

    h2->Write();

    for (auto& h : histos_2d_global) {
        h->Write();
    }

    for (auto& p : profiles) {
        p->Write();
    }

    for (auto& [layer, hlist] : histos_1d) {

        std::string dir_name = "layer_" + std::to_string(layer);
        TDirectory* dir = output_file.mkdir(dir_name.c_str());
        if (!dir) dir = output_file.GetDirectory(dir_name.c_str());

        for (auto& h : hlist) {
            h->SetFillColor(kAzure - 9);
            dir->WriteTObject(h.GetPtr());
        }
    }

    for (auto& [layer, hlist] : histos_2d) {

        std::string dir_name = "layer_" + std::to_string(layer);
        TDirectory* dir = output_file.GetDirectory(dir_name.c_str());
        if (!dir) dir = output_file.mkdir(dir_name.c_str());

        for (auto& h : hlist) {
            dir->WriteTObject(h.GetPtr());
        }
    }

    output_file.Close();
}
