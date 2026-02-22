#include "Nontracktrigger.h"

///////////////////////////////////////////////////////////////////////////////////////

TriggerNonTrack::TriggerNonTrack(fhicl::ParameterSet const& p):
    _PrepareClusterParameters(),
    _GradientBoostingModel(
        p.get<std::string>("model_filename"),
        p.get<float>("thr_opt")
    )
{ }


bool TriggerNonTrack::run_algorithm(std::vector<Hit>& hits, const art::Event& event)
{
    std::cout << "\n=== Event ID: " << event.id().event() << " ===\n" << std::endl;

    std::cout << "sum: " << hits.size()  << "\n";
    sum+=hits.size();
    auto& features = _PrepareClusterParameters.PrepareFeatures(hits);

    std::cout << "HITS: " << hits.size()   << "\n";    

    _GradientBoostingModel(hits, features);
    
    std::cout << "HITS_after_gb: " << hits.size()   << "\n";

    return !hits.empty();
}

///////////////////////////////////////////////////////////////////////////////////////

std::vector<Features>& PrepareClusterParameters::PrepareFeatures(std::vector<Hit>& hits)
{
    _features.clear();

    struct Accumulator {
        int n = 0;
        float adc = 0.0;
        int n_x = 0;
        int n_y = 0;
        float adc_x = 0.0;
        float adc_y = 0.0;
    };

    std::unordered_map<uint64_t, Accumulator> acc;
    acc.reserve(hits.size());

    for (const auto& h : hits) {
        auto& a = acc[h.hitSet_id];

        a.n++;
        a.adc += h.adc;

        if (h.view == 1) {
            a.n_x++;
            a.adc_x += h.adc;
        } else {
            a.n_y++;
            a.adc_y += h.adc;
        }
    }

    _features.reserve(acc.size());

    for (const auto& kv : acc) {
        uint64_t key = kv.first;
        const auto& a = kv.second;

        Features f{};
        
        f.hitSet_id = static_cast<int>(key & 0xFFFFFFFF);
        
        f.n = static_cast<float>(a.n);
        f.adc = a.adc;

        f.adc_per_hit = (a.n > 0) ? a.adc / a.n : 0.0;

        int n_sum = a.n_x + a.n_y;
        f.n_asym = (n_sum > 0)
            ? float(a.n_x - a.n_y) / n_sum
            : 0.0;

        float adc_sum = a.adc_x + a.adc_y;
        f.adc_asym = (adc_sum > 0.0)
            ? (a.adc_x - a.adc_y) / adc_sum
            : 0.0;

        _features.push_back(f);
    }

    return _features;
}

///////////////////////////////////////////////////////////////////////////////////////

GradientBoostingModel::GradientBoostingModel(const std::string& model_filename, float thr_opt):
    thr_opt(thr_opt)
{
    std::ifstream file(model_filename);
    
    std::string line;
    int n_trees = 0;
    
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string key;
        
        if (std::getline(iss, key, ':')) {
            if (key == "init_lo") {
                iss >> init_lo;
            } 
            else if (key == "learning_rate") {
                iss >> learning_rate;
            } 
            else if (key == "n_trees") {
                iss >> n_trees;
                trees.resize(n_trees);
            }
            else if (key.find("tree_") == 0) {
                // Формат: tree_X_nodes:N
                size_t pos = key.find("_nodes");
                if (pos != std::string::npos) {
                    // Это заголовок дерева - размерность уже установлена
                }
            }
            else if (key == "node") {
                // Формат: node:tree_idx:node_idx:feature:threshold:value:left:right:is_leaf
                std::string token;
                std::vector<std::string> tokens;
                
                while (std::getline(iss, token, ':')) {
                    tokens.push_back(token);
                }
                
                if (tokens.size() >= 8) {
                    int tree_idx = std::stoi(tokens[0]);
                    int node_idx = std::stoi(tokens[1]);
                    
                    if (tree_idx >= 0 && tree_idx < trees.size()) {
                        TreeNode node;
                        node.feature = std::stoi(tokens[2]);
                        node.threshold = std::stof(tokens[3]);
                        node.value = std::stof(tokens[4]);
                        node.left_child = std::stoi(tokens[5]);
                        node.right_child = std::stoi(tokens[6]);
                        node.is_leaf = (tokens[7] == "1" || tokens[7] == "true");
                        
                        if (node_idx >= trees[tree_idx].nodes.size()) {
                            trees[tree_idx].nodes.resize(node_idx + 1);
                        }
                        trees[tree_idx].nodes[node_idx] = node;
                    }
                }
            }
        }
    }
    
    file.close();
}

float Tree::predict(const std::vector<float>& x_test) const {
    if (nodes.empty()) return 0.0;
    
    int current_node = 0;
    
    while (true) {
        const TreeNode& node = nodes[current_node];
        
        if (node.is_leaf) {
            return node.value;
        }
        
        // may be optimized, because not all vector x_test is needed within Tree::predict()
        if (x_test[node.feature] <= node.threshold) {
            current_node = node.left_child;
        } else {
            current_node = node.right_child;
        }
    }
}

// Предсказание вероятностей для одного образца
std::vector<float> GradientBoostingModel::predict_proba_one(const std::vector<float>& features) const {
    float total = init_lo;
    
    // Суммируем предсказания всех деревьев
    for (const auto& tree : trees) {
        total += learning_rate * tree.predict(features);
    }
    
    // Преобразуем в вероятности через сигмоид
    float p1 = 1.0f / (1.0f + std::exp(-total));
    
    return {1.0f - p1, p1};
}

// Предсказание вероятностей для матрицы
std::vector<std::vector<float>> GradientBoostingModel::predict_proba(const std::vector<Features>& features) const 
{
    std::vector<std::vector<float>> result;
    result.reserve(features.size());
    
    for (const auto& f : features) {
        result.push_back(predict_proba_one({ f.n, f.adc, f.adc_per_hit, f.n_asym, f.adc_asym }));
    }
    
    return result;
}

void GradientBoostingModel::operator()(std::vector<Hit>& hits, std::vector<Features>& features)
{
    std::vector<Hit> out;

    std::vector<std::vector<float>> proba = predict_proba(features);
    
    for (size_t i = 0; i < features.size(); ++i) {
        if (proba[i][1] > thr_opt) { // 0.63875492f
            out = std::move(hits);
            //Hit h;
            //h.hitSet_id = features[i].hitSet_id;
            //out.push_back(h);
        }
    }

    hits = std::move(out);
}

///////////////////////////////////////////////////////////////////////////////////////
