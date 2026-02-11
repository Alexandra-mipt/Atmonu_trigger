#include "Trigger.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <iomanip>
#include <string>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <array>

///////////////////////////////////////////////////////////////////////////////////////

inline uint64_t make_key(int eventID, int hitSet_id) {
    return (uint64_t(uint32_t(eventID)) << 32) | uint32_t(hitSet_id);
}

struct Features {
    float n;
    float adc;
    float adc_per_hit;
    float n_asym;
    float adc_asym;
    
    int eventID;
    int hitSet_id;
};

class PrepareClusterParameters {
public:
    PrepareClusterParameters() = default;

    void prepare_features(std::vector<Hit>& hits);

    std::vector<Features> &get_features();
private:
    std::vector<Features> vFeatures;
};

///////////////////////////////////////////////////////////////////////////////////////

// Structure for tree node
struct TreeNode {
    int feature;           // Feature index for splitting
    float threshold;      // Threshold value for splitting
    float value;          // Value in leaf node
    int left_child;        // Left child index (-1 for leaf)
    int right_child;       // Right child index (-1 for leaf)
    bool is_leaf;          // Leaf node flag
};

struct Tree {
    std::vector<TreeNode> nodes;
    
    float predict(const std::vector<float>& features) const;
};

class GradientBoostingModel {
public:
    GradientBoostingModel(const std::string& model_filename, float thr_opt);

    void filter(std::vector<Hit>& hits, const std::vector<Features> vFeatures);
private:
    std::vector<float> predict_proba_one(const std::vector<float>& features) const;
    
    std::vector<std::vector<float>> predict_proba(const std::vector<Features>& features) const;
    

    float init_lo;                // Initial log odds
    float learning_rate;          // Learning rate
    std::vector<Tree> trees;      // Trees
    float thr_opt;                // Optimal threshold for filtering
};

///////////////////////////////////////////////////////////////////////////////////////

class TriggerNonTrack {
public:
    TriggerNonTrack(
        const std::string& model_filename, float thr_opt
    );
    
    void fPrepareClusterParameters(std::vector<Hit>& hits);
    
    std::vector<Features> &GetPrepareClusterParameters();

    void fGradientBoostingModel(
        std::vector<Hit>& hits,
        const std::vector<Features>& vFeatures
    );

private:
    std::string model_filename;
    PrepareClusterParameters _PrepareClusterParameters;
    GradientBoostingModel _GradientBoostingModel;
};

///////////////////////////////////////////////////////////////////////////////////////
