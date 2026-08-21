#ifndef SNDHLLHC_SISTRIP_CLUSTER_H
#define SNDHLLHC_SISTRIP_CLUSTER_H

#include <cstdint>

struct Module {
    int layer, row, col;

    bool operator<(const Module& other) const {
        if (layer != other.layer) return layer < other.layer;
        if (row != other.row) return row < other.row;
        return col < other.col;
    }
};

class SiStripCluster {
public:
    SiStripCluster(uint32_t detector_id, uint32_t adc, size_t size, int layer, int row, int column, bool is_vertical) : 
        detector_id_(detector_id), adc_(adc), size_(size), layer_(layer), row_(row), column_(column), is_vertical_(is_vertical) {}

    inline uint32_t GetDetectorId() const { return detector_id_; }
    inline uint32_t GetSignal() const { return adc_; }
    inline size_t GetSize() const { return size_; }
    inline int GetLayer() const { return layer_; }
    inline int GetRow() const { return row_; }
    inline int GetColumn() const { return column_; }
    inline int IsVertical() const { return is_vertical_; }

private:
    uint32_t detector_id_; // Of the strip with most adc
    uint32_t adc_;
    size_t size_;
    int layer_;
    int row_;
    int column_;
    bool is_vertical_;
};

struct SiStripClusteringProducts {
    std::vector<SiStripDigi> digis;
    std::vector<SiStripCluster> clusters;
};

#endif