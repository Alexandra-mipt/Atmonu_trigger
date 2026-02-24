#ifndef GUARD_DDT_ATMNU_TRIGGER_HH
#define GUARD_DDT_ATMNU_TRIGGER_HH



#include "DDTBaseDataProducts/BaseProducts.h"
#include "DDTBaseDataProducts/DAQHit.h"
#include "DDTBaseDataProducts/HitList.h"
#include "DDTBaseDataProducts/TriggerDecision.h"

#include <art/Framework/Principal/Handle.h>
#include <fhiclcpp/ParameterSet.h>

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

#include "DDTBaseDataProducts/TriggerDecision.h"
#include "Tracktrigger.h"

///////////////////////////////////////////////////////////////////////////////////////

struct Features {
    float n;
    float adc;
    float adc_per_hit;
    float n_asym;
    float adc_asym;
    
    int hitSet_id;
};

class PrepareClusterParameters {
public:
    PrepareClusterParameters() = default;

    std::vector<Features>& PrepareFeatures(std::vector<Hit>& hits);

private:
    std::vector<Features> _features;
};

///////////////////////////////////////////////////////////////////////////////////////

// Структура для узла дерева
struct TreeNode {
    int feature;           // Индекс признака для разделения
    float threshold;       // Порог для разделения
    float value;           // Значение в листе
    int left_child;        // Индекс левого потомка (-1 для листа)
    int right_child;       // Индекс правого потомка (-1 для листа)
    bool is_leaf;          // Флаг листа
};

// Структура для дерева
struct Tree {
    std::vector<TreeNode> nodes;
    
    float predict(const std::vector<float>& features) const;
};

// Структура для хранения модели градиентного бустинга
class GradientBoostingModel {
public:
    GradientBoostingModel(const std::string& model_filename, float thr_opt);

    void operator()(std::vector<Hit>& hits, std::vector<Features>& features);
private:
    // Предсказание вероятностей для одного образца
    std::vector<float> predict_proba_one(const std::vector<float>& features) const;
    
    // Предсказание вероятностей для матрицы
    std::vector<std::vector<float>> predict_proba(const std::vector<Features>& features) const;
    
    // Предсказание классов
    // std::vector<int> predict(const std::vector<std::vector<float>>& X) const;

    float init_lo;                // Начальный логарифм шансов
    float learning_rate;          // Скорость обучения
    std::vector<Tree> trees;      // Деревья
    float thr_opt;                // Оптимальный порог для фильтрации
};


///////////////////////////////////////////////////////////////////////////////////////

class TriggerNonTrack {
public:
    TriggerNonTrack(fhicl::ParameterSet const& p);
    
    bool run_algorithm(std::vector<Hit>& hits, const art::Event& event);

    std::vector<novaddt::TriggerDecision> TriggerDecisions() const;

    size_t sum=0;
private:
    std::string model_filename;
    PrepareClusterParameters _PrepareClusterParameters;
    GradientBoostingModel    _GradientBoostingModel;
    std::vector<novaddt::TriggerDecision> _triggerDecisions;
};

///////////////////////////////////////////////////////////////////////////////////////

#endif 
