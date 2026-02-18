#include "Tracktrigger.h"
#include "art/Framework/Principal/Event.h"
//#include <cmath>
//#include <algorithm>


double DistanceToTrack(
double x1, double x2,
double y1, double y2,
double x, double y)
{
double Rx = x2 - x1;
double Ry = y2 - y1;

double rx = x - x1;
double ry = y - y1;

double Rabs = std::sqrt(Rx*Rx + Ry*Ry);
if (Rabs == 0) return 1e9;

double Lx = (Rx*rx + Ry*ry) / Rabs;
double Ly2 = rx*rx + ry*ry - Lx*Lx;
double Ly = (Ly2 >= 0) ? std::sqrt(Ly2) : 100.0;

double distance_to_end = std::max({Lx - Rabs, -Lx, 0.0});

return std::sqrt(Ly*Ly + distance_to_end*distance_to_end);
}

// Trigger

Trigger::Trigger(fhicl::ParameterSet const & p):
_ParametersCuts(
    p.get<double>("gaps_min"), 
    p.get<double>("gaps_max"),
    p.get<double>("length_min"),
    p.get<double>("length_max"),
    p.get<double>("max_gap"), 
    p.get<double>("epsilon"),
    p.get<double>("total_hits_min"),
    p.get<double>("total_hits_max")),
_Score(        
    p.get<double>("length_weight"),
    p.get<double>("hits_weight"),
    p.get<double>("gaps_weight"),
    p.get<double>("score_threshold")),
_AngleCut(
    p.get<double>("xz_cut"),
    p.get<double>("yz_cut"),
    p.get<double>("xy_cut")),
_SeparateHitsWithTracksRestActivity(
    p.get<double>("epsilon")),
_RestActivityCut(p.get<double>("min_hits")
)
{
}

bool Trigger::run_algorithm(
std::vector<Hit>& hits,
std::vector<Track>& tracks,
const art::Event& event)
{
/*    std::cout << "\n=== Event ID: " << event.id().event() << " ===\n" << std::endl;
    std::cout << "========== ВХОДНЫЕ ДАННЫЕ ==========" << std::endl;
    std::cout << "  Всего треков: " << tracks.size() << std::endl;
    std::cout << "  Всего хитов: " << hits.size() << std::endl;
    
    // Вывод первых 5 треков
    if (!tracks.empty()) {
        std::cout << "\n  Первые 5 треков:" << std::endl;
        std::cout << "  idx\tfView\tStartX\tStartY\tStartZ\tEndX\tEndY\tEndZ\tsliceID" << std::endl;
        for (size_t i = 0; i < std::min(size_t(5), tracks.size()); ++i) {
            const auto& t = tracks[i];
            std::cout << "  " << i << "\t"
                      << t.fView << "\t"
                      << std::fixed << std::setprecision(1) << t.StartX << "\t"
                      << t.StartY << "\t"
                      << t.StartZ << "\t"
                      << t.EndX << "\t"
                      << t.EndY << "\t"
                      << t.EndZ << "\t"
                      << t.sliceID << std::endl;
        }
    }
    
    // Вывод первых 10 хитов
    if (!hits.empty()) {
        std::cout << "\n  Первые 10 хитов:" << std::endl;
        std::cout << "  idx\thitSet_id\tplane\tcell\tadc\ttdc\tview\tused" << std::endl;
        for (size_t i = 0; i < std::min(size_t(100), hits.size()); ++i) {
            const auto& h = hits[i];
            std::cout << "  " << i << "\t"
                      << h.hitSet_id << "\t\t"
                      << h.plane << "\t"
                      << h.cell << "\t"
                      << h.adc << "\t"
                      << h.tdc << "\t"
                      << h.view << "\t"
                      << h.used << std::endl;
        }
    }
    std::cout << "====================================\n" << std::endl;
    // ===== КОНЕЦ ВЫВОДА ВХОДНЫХ ДАННЫХ =====
 */    
 //   std::cout << "STEP 1: Parameter Cuts" << std::endl;
//    std::cout << "  Input tracks: " << tracks.size() << std::endl;
//    std::cout << "  Input hits: " << hits.size() << std::endl;
    
    auto& paramTracks = _ParametersCuts(hits, tracks);
    
//    std::cout << "  Output tracks after parameter cuts: " << paramTracks.size() << std::endl;
 //   if (!paramTracks.empty()) {
 //   std::cout << "  Все треки:" << std::endl;
 //   for (size_t i = 0; i < paramTracks.size(); ++i) {
 //       const auto& t = paramTracks[i];
 //       std::cout << "    Трек " << i << ": sliceID=" << t.sliceID 
 //                 << ", length=" << t.track_length
  //                << ", total_hits=" << t.total_hits
 //                 << ", gaps=" << t.total_gaps << std::endl;
 //   }
//}
//    std::cout << std::endl;     
//    std::cout << "STEP 2: Score" << std::endl;
    const auto& scoredTracks = _Score(paramTracks);
    
//    std::cout << "  Output tracks after score cut: " << scoredTracks.size() << std::endl;
//    if (!scoredTracks.empty()) {
//        const auto& t = scoredTracks[0];
//        std::cout << "  First track: sliceID=" << t.sliceID 
//                  << ", score=" << t.score << std::endl;
//    }
//    std::cout << std::endl;

//    std::cout << "STEP 3: Angle Cut" << std::endl;
    const auto& angleTracks = _AngleCut(scoredTracks);
    
//    std::cout << "  Output tracks after angle cut: " << angleTracks.size() << std::endl;
//    if (!angleTracks.empty()) {
//        const auto& t = angleTracks[0];
//        std::cout << "  First track: sliceID=" << t.sliceID << std::endl;
//    }
//    std::cout << std::endl;

//    std::cout << "STEP 4: Choose Track Hits" << std::endl;
    _ChooseTrackHits.Filter(hits, angleTracks);
    
//    std::cout << "  Selected hits: " << _ChooseTrackHits.get_hits().size() << std::endl;
//    std::cout << "  Selected tracks: " << _ChooseTrackHits.get_tracks().size() << std::endl;
//    std::cout << std::endl;

//    std::cout << "STEP 5: Separate Hits With Tracks Rest Activity" << std::endl;
    _SeparateHitsWithTracksRestActivity.Filter(
        _ChooseTrackHits.get_hits(),
        _ChooseTrackHits.get_tracks());
    
//    std::cout << "  Rest hits: " << _SeparateHitsWithTracksRestActivity.get_rest_hits().size() << std::endl;
//    std::cout << std::endl;

//    std::cout << "STEP 6: Rest Activity Cut" << std::endl;
    _RestActivityCut.Filter(_SeparateHitsWithTracksRestActivity.get_rest_hits());
    
//    std::cout << "  Final filtered hits: " << _RestActivityCut.get_rest_hits().size() << std::endl;
//    std::cout << std::endl;

    return !_RestActivityCut.get_rest_hits().empty();
}

// ParameterCuts

ParameterCuts::ParameterCuts(double gaps_min,       double gaps_max,                                                         
                         double length_min,     double length_max,
                         double max_gap,        double epsilon,                                                          
                         double total_hits_min, double total_hits_max):
gaps_min(gaps_min),
gaps_max(gaps_max),
length_min(length_min),
length_max(length_max),
max_gap(max_gap),
epsilon(epsilon),
total_hits_min(total_hits_min),
total_hits_max(total_hits_max)
{
}

std::vector<TrackParams>&
ParameterCuts::operator()(
const std::vector<Hit>& hits,
const std::vector<Track>& tracks)
{
table.clear();

for (const auto& tr : tracks)
{
    TrackParams tp =
        CalculateTrackParams(hits, tr);

    if (FilterTrack(tp))
        table.push_back(tp);
}

return table;
}

TrackParams ParameterCuts::CalculateTrackParams(
const std::vector<Hit>& hits,
const Track& track)
{
TrackParams p;
p.track = track;
p.sliceID = track.sliceID;

std::vector<Hit> corridor;
for (const auto& h : hits)
{

    if (h.hitSet_id != track.sliceID) continue;

    double y1 = (h.view == 1) ? track.StartX : track.StartY;
    double y2 = (h.view == 1) ? track.EndX : track.EndY;
    double d = DistanceToTrack(
        track.StartZ, track.EndZ,
        y1, y2,
        h.plane, h.cell);

    if (d <= epsilon)
        corridor.push_back(h);
}

std::sort(corridor.begin(), corridor.end(),
    [](const Hit& a, const Hit& b){ return a.plane < b.plane; });

p.total_hits = corridor.size();

int total_gaps = 0;
int max_gap = 0;

for (size_t i = 1; i < corridor.size(); ++i)
{
    int gap = corridor[i].plane - corridor[i-1].plane - 1;
    if (gap > 0)
    {
        total_gaps += gap;
        if (gap > max_gap) max_gap = gap;
    }
}

p.total_gaps = total_gaps;
p.max_gap_size = max_gap;

double dx = track.EndX - track.StartX;
double dy = track.EndY - track.StartY;
double dz = track.EndZ - track.StartZ;

p.track_length = std::sqrt(dx*dx + dy*dy + dz*dz);
p.gaps_per_length = total_gaps / (std::abs(dz) + 1.0);
//    std::cout << "  Все треки:" << std::endl;
//        std::cout  << ": sliceID=" << p.sliceID
//                  << ", length=" << p.track_length
//                  << ", total_hits=" << p.total_hits
//                  << ", gaps=" << p.total_gaps <<"Track!"<<p.track.StartX<<"\n";
//std::cout << std::endl;


return p;
}


bool ParameterCuts::FilterTrack(
const TrackParams& t)
{
if (t.gaps_per_length < gaps_min)
    return false;

if (t.gaps_per_length > gaps_max)
    return false;

if (t.track_length < length_min)
    return false;

if (t.track_length > length_max)
    return false;

if (t.max_gap_size > max_gap)
    return false;

if (t.total_hits < total_hits_min)
    return false;

if (t.total_hits > total_hits_max)
    return false;

return true;
}

// Score

Score::Score(double length_weight, double hits_weight,
         double gaps_weight,   double score_threshold):
length_weight(length_weight),
hits_weight(hits_weight),
gaps_weight(gaps_weight),
score_threshold(score_threshold)

{
}
//////
const std::vector<TrackParams>&
Score::operator()(std::vector<TrackParams>& in)
{
    table.clear();
    
    // Шаг 1: вычисляем скоры для всех треков и сохраняем во временный вектор
    
//    std::vector<double> lengths, hits, gaps;

//    for (const auto& t : in) {
     //   lengths.push_back(t.track_length);
     //   hits.push_back(static_cast<double>(t.total_hits));
    //    gaps.push_back(t.gaps_per_length);
   // }

   // std::sort(lengths.begin(), lengths.end());
   // std::sort(hits.begin(), hits.end());
   // std::sort(gaps.begin(), gaps.end());

    double m_len = 70.285133;/*Quantile(lengths, 0.5);*/
    double q1_len = 54.196402;/*Quantile(lengths, 0.25);*/
    double q3_len = 96.341832;/*Quantile(lengths, 0.75);*/

    double m_hit = 47;/*Quantile(hits, 0.5);*/
    double q1_hit = 35;/*Quantile(hits, 0.25);*/
    double q3_hit = 63;/*Quantile(hits, 0.75);*/

    double m_gap = 0.2;/*Quantile(gaps, 0.5);*/
    double q1_gap = 0.058823529; /*Quantile(gaps, 0.25);*/
    double q3_gap = 0.33333333; /*Quantile(gaps, 0.75);*/

    // Вычисляем score для каждого трека
    for (auto& t : in) {
        double norm_len = RobustNormalize(t.track_length, m_len, q1_len, q3_len);
        double norm_hit = RobustNormalize(t.total_hits,  m_hit, q1_hit, q3_hit);
        double norm_gap = -RobustNormalize(t.gaps_per_length, m_gap, q1_gap, q3_gap);


        t.score = length_weight * norm_len +
                  hits_weight   * norm_hit +
                  gaps_weight   * norm_gap;
 //       std::cout<<t.score<<"\n";
    
    }

    // Шаг 2: для каждого slice_id выбираем трек с максимальным score
    std::map<int, TrackParams> best_track_per_slice;  // ключ = slice_id
    
    for (const auto& t : in) {
        auto it = best_track_per_slice.find(t.sliceID);
        if (it == best_track_per_slice.end()) {
            // Первый трек в этом slice_id
            best_track_per_slice[t.sliceID] = t;
        } else if (t.score > it->second.score) {
            // Нашли трек с большим score в этом slice_id
            it->second = t;
        }
    }
    
    // Шаг 3: добавляем в table только те треки, у которых score > threshold
    for (const auto& pair : best_track_per_slice) {
        const auto& t = pair.second;
        if (t.score > score_threshold) {  // используем >, не >=
            table.push_back(t);
        }
    }
    
//    for (const auto& t: table){
 //   std::cout<<t.score<<"\n";
//}



    return table;
}

/////

/*
double Score::Quantile(const std::vector<double>& sorted, double q)
{
    if (sorted.empty()) return 0.0;

    double pos = (sorted.size() - 1) * q;
    size_t idx = static_cast<size_t>(std::floor(pos));
    double frac = pos - idx;

    if (idx + 1 < sorted.size())
        return sorted[idx] + frac * (sorted[idx + 1] - sorted[idx]);
    else
        return sorted[idx];
}
*/

double Score::RobustNormalize(
    double value, double median, double q1, double q3
)
{
    double iqr = q3 - q1;
    if (iqr == 0.0) return 0.0;
    return (value - median) / iqr;
}





// AngleCut

AngleCut::AngleCut(double xz_cut, double yz_cut, double xy_cut):
    xz_cut(xz_cut),
    yz_cut(yz_cut),
    xy_cut(xy_cut) 
{
}

const std::vector<TrackParams>&
AngleCut::operator()(const std::vector<TrackParams>& in)
{
table.clear();

for (const auto& t : in)
{
    double dx = t.track.EndX - t.track.StartX;
    double dy = t.track.EndY - t.track.StartY;
    double dz = t.track.EndZ - t.track.StartZ;

    double angle_xz = std::atan2(dz, dx) * 180.0/ M_PI;
    double angle_yz = std::atan2(dz, dy) * 180.0/ M_PI;
    double angle_xy = std::atan2(dx, dy) * 180.0/ M_PI;

 //   std::cout<<"Angle xz"<<angle_xz<<"\n";
 //   std::cout<<"Angle yz"<<angle_yz<<"\n";
 //   std::cout<<"Angle xy"<<angle_xy<<"\n";
   
   
 //  std::cout<<"xz_cut"<<xz_cut<<"yz_cut"<<yz_cut<<"xy_cut"<<xy_cut<<"\n";


    if ((angle_xz <= xz_cut) ||
     (angle_xz >= (180.0 - xz_cut))) continue;

    if ((angle_yz <=  yz_cut) ||
       (angle_yz >= (180.0 - yz_cut))) continue;
    
    if (std::abs(angle_xy) >= (180.0 - xy_cut)) continue;

    table.push_back(t);
}

return table;
}

// ChooseTrackHits


void ChooseTrackHits::Filter(
const std::vector<Hit>& hits,
const std::vector<TrackParams>& tracks)
{
selected_hits.clear();
selected_tracks.clear();

    // Только sliceID (hitSet_id) из треков
    std::unordered_set<int> track_slice_ids;
    track_slice_ids.reserve(tracks.size());

    for (const auto& t : tracks) {
        selected_tracks.push_back(t.track);
        track_slice_ids.insert(t.sliceID);
    }

    // фильтрация хитов
    selected_hits.reserve(hits.size());

    for (const auto& h : hits) {
        if (track_slice_ids.count(h.hitSet_id)) {
            selected_hits.push_back(h);
        }
    }

 
}

/*
void ChooseTrackHits::Filter(
const std::vector<Hit>& hits,
const std::vector<TrackParams>& tracks)
{
selected_hits.clear();
selected_tracks.clear();

for (const auto& t : tracks)
{
    selected_tracks.push_back(t.track);
    for (const auto& h : hits){
        if (h.hitSet_id == t.sliceID)
            selected_hits.push_back(h);
    }
}
}
*/

std::vector<Hit>&
ChooseTrackHits::get_hits()
{
return selected_hits;
}

const std::vector<Track>&
ChooseTrackHits::get_tracks() const
{
return selected_tracks;
}

// SeparateHitsWithTracksRestActivity

SeparateHitsWithTracksRestActivity::SeparateHitsWithTracksRestActivity(double epsilon):
epsilon(epsilon)
{
}
void SeparateHitsWithTracksRestActivity::Filter(
    std::vector<Hit>& hits,
    const std::vector<Track>& tracks
)
    
{
    rest_hits.clear();
    for (const auto& track : tracks) {
        auto hit_indices =
            FilterHitsAroundTrack(hits, track);

        for (int idx : hit_indices) {
            hits[idx].used = true;
        }
    }

    for (const auto& h : hits) {
        if (!h.used)
            rest_hits.push_back(h);
    }   
}

std::vector<int> SeparateHitsWithTracksRestActivity::FilterHitsAroundTrack(
    const std::vector<Hit>& hits,
    const Track& track
)
{
    std::vector<int> indices;

    for (size_t i = 0; i < hits.size(); ++i) {
        if (hits[i].used) continue;
        if (hits[i].hitSet_id != track.sliceID) continue;

        double d = HitDistanceToTrack(hits[i], track);
        if (d <= epsilon)
            indices.push_back(i);
    }

    std::sort(indices.begin(), indices.end(),
              [&](int a, int b){
                  return hits[a].plane < hits[b].plane;
              });

    return indices;
}


double SeparateHitsWithTracksRestActivity::HitDistanceToTrack(const Hit& h, const Track& t)
{
    double y1 = (h.view == 1) ? t.StartX : t.StartY;
    double y2 = (h.view == 1) ? t.EndX   : t.EndY;

    return DistanceToTrack(
        t.StartZ, t.EndZ,
        y1, y2,
        h.plane, h.cell
    );
}
std::vector<Hit>& SeparateHitsWithTracksRestActivity::get_rest_hits()
{
    return rest_hits;
}


// RestActivityCut

RestActivityCut::RestActivityCut(double min_hits):
min_hits(min_hits)
{
}
void RestActivityCut::Filter(const std::vector<Hit>& hits)
{
    filtered_hits.clear();
    std::map<int, int> group_counts;
    std::map<int, int> hits_counts;
    
    for (const auto& h : hits) group_counts[h.hitSet_id]++;

    for (const auto& h : hits) {
        if (group_counts[h.hitSet_id] > min_hits) {filtered_hits.push_back(h);
     hits_counts[h.hitSet_id] = group_counts[h.hitSet_id];
             }
                 }

                     int tot_groups = group_counts.size();

                         int passed = hits_counts.size();

                             std::cout << "Hits:\n";
                                 std::cout << "  All: " << tot_groups << '\n';
                                     std::cout << "  Events with >" << min_hits << " hits: "
                                                   << passed << " ("
                                                                 << std::fixed << std::setprecision(3)
                                                                               << 100.0 * passed / tot_groups << "%)\n";
    
    }



const RestHits&
RestActivityCut::get_rest_hits() const
{
return filtered_hits;
}

