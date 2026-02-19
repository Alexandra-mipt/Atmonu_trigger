#include "Tracktrigger.h"

///////////////////////////////////////////////////////////////////////////////////////

Trigger::Trigger(fhicl::ParameterSet const& p):
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
{ }

bool Trigger::run_algorithm(
    std::vector<Hit>& hits,std::vector<Track>& tracks, const art::Event& event
)
{
    std::cout << "\n=== Event ID: " << event.id().event() << " ===\n" << std::endl;
    
    auto& paramTracks = _ParametersCuts(hits, tracks);

    auto& scoreTracks = _Score(paramTracks);
        
    auto& angleTracks = _AngleCut(scoreTracks);
        
    _ChooseTrackHits(hits, angleTracks);
        
    auto& filteredHits   = _ChooseTrackHits.GetHits();
    auto& filteredTracks = _ChooseTrackHits.GetTracks();
    _SeparateHitsWithTracksRestActivity(filteredHits, filteredTracks);
    
    auto& restHits = _SeparateHitsWithTracksRestActivity.GetHits();
    _RestActivityCut(restHits);

    return !_RestActivityCut.GetHits().empty();
}

///////////////////////////////////////////////////////////////////////////////////////

double DistanceToTrack(double x1, double x2, double y1, double y2, double x, double y)
{
    double Rx = x2 - x1;
    double Ry = y2 - y1;
    double r_x = x - x1;
    double r_y = y - y1;

    double Rabs = std::sqrt(Rx*Rx + Ry*Ry);
    if (Rabs == 0) return 1e9;

    double Lx = (Rx*r_x + Ry*r_y) / Rabs;
    double Ly2 = r_x*r_x + r_y*r_y - Lx*Lx;
    double Ly = (Ly2 >= 0) ? std::sqrt(Ly2) : 100.0;

    double dist_end = std::max({Lx - Rabs, -Lx, 0.0});
    return std::sqrt(Ly*Ly + dist_end*dist_end);
}

///////////////////////////////////////////////////////////////////////////////////////

ParameterCuts::ParameterCuts(
    double gaps_min,       double gaps_max,                                                         
    double length_min,     double length_max,
    double max_gap,        double epsilon,                                                          
    double total_hits_min, double total_hits_max
):
    gaps_min(gaps_min),             gaps_max(gaps_max),
    length_min(length_min),         length_max(length_max),
    max_gap(max_gap),               epsilon(epsilon),
    total_hits_min(total_hits_min), total_hits_max(total_hits_max)
{ }

std::vector<TrackParams>& ParameterCuts::operator()(const std::vector<Hit>& hits, const std::vector<Track>& tracks)
{
    _paramTracks.clear();

    for (auto& tr : tracks) {
        TrackParams params = CalculateTrackParams(hits, tr);
        
        if (FilterTrack(params)) _paramTracks.push_back(params);
    }

    return _paramTracks;
}

TrackParams ParameterCuts::CalculateTrackParams(const std::vector<Hit>& hits, const Track& track)
{
    TrackParams p;
    p.sliceID = track.sliceID;
    p.track   = track;

    std::vector<Hit> corridor;
    for (const auto& h : hits) {
        if (h.hitSet_id != track.sliceID) continue;

        double y1 = (h.view == 1) ? track.StartX : track.StartY;
        double y2 = (h.view == 1) ? track.EndX : track.EndY;

        double d = DistanceToTrack(track.StartZ, track.EndZ, y1, y2, h.plane, h.cell);

        if (d <= epsilon) corridor.push_back(h);
    }

    std::sort(corridor.begin(), corridor.end(),
        [](const Hit& a, const Hit& b){ return a.plane < b.plane; });

    p.total_hits = corridor.size();

    int total_gaps = 0;
    int max_gap    = 0;
    for (size_t i = 1; i < corridor.size(); ++i) {
        int gap = corridor[i].plane - corridor[i-1].plane - 1;
        if (gap > 0) {
            total_gaps += gap;
            if (gap > max_gap) max_gap = gap;
        }
    }

    p.total_gaps = total_gaps;
    p.max_gap_size = max_gap;

    double dx = track.StartX - track.EndX;
    double dy = track.StartY - track.EndY;
    double dz = track.StartZ - track.EndZ;
    p.track_length = std::sqrt(dx*dx + dy*dy + dz*dz);

    p.gaps_per_length = total_gaps / (std::abs(dz) + 1.0);

    return p;
}

bool ParameterCuts::FilterTrack(const TrackParams& t)
{
    if (t.gaps_per_length < gaps_min      || t.gaps_per_length > gaps_max    ) return false;
    if (t.track_length    < length_min       || t.track_length    > length_max     ) return false;
    if (t.max_gap_size    > max_gap                                          ) return false;
    if (t.total_hits      < total_hits_min  || t.total_hits      > total_hits_max) return false;
    return true;
}

///////////////////////////////////////////////////////////////////////////////////////

Score::Score(
    double length_weight, double hits_weight,
    double gaps_weight,   double score_threshold
):
    length_weight(length_weight), hits_weight(hits_weight),
    gaps_weight(gaps_weight),     score_threshold(score_threshold)
{ }

std::vector<TrackParams>& Score::operator()(std::vector<TrackParams>& tracks)
{
    _scoreTracks.clear();

    ScoreTracks(tracks);
    
    SelectBestTracks(tracks);

    return _scoreTracks;
}

void Score::ScoreTracks(std::vector<TrackParams>& tracks)
{
/*
    std::vector<double> lengths, hits, gaps;

    for (const auto& t : tracks) {
        lengths.push_back(t.track_length);
        hits.push_back(static_cast<double>(t.total_hits));
        gaps.push_back(t.gaps_per_length);
    }

    std::sort(lengths.begin(), lengths.end());
    std::sort(hits.begin(), hits.end());
    std::sort(gaps.begin(), gaps.end());

    double m_len  = Quantile(lengths, 0.5);
    double q1_len = Quantile(lengths, 0.25);
    double q3_len = Quantile(lengths, 0.75);

    double m_hit  = Quantile(hits, 0.5);
    double q1_hit = Quantile(hits, 0.25);
    double q3_hit = Quantile(hits, 0.75);

    double m_gap  = Quantile(gaps, 0.5);
    double q1_gap = Quantile(gaps, 0.25);
    double q3_gap = Quantile(gaps, 0.75);
*/

    double m_len  = 70.285133;
    double q1_len = 54.196402;
    double q3_len = 96.341832;

    double m_hit  = 47;
    double q1_hit = 35;
    double q3_hit = 63;

    double m_gap  = 0.2;
    double q1_gap = 0.058823529;
    double q3_gap = 0.33333333;

    std::vector<TrackParams> scored;

    for (auto& t : tracks) {
        double norm_len =  RobustNormalize(t.track_length, m_len, q1_len, q3_len);
        double norm_hit =  RobustNormalize(t.total_hits,  m_hit, q1_hit, q3_hit);
        double norm_gap = -RobustNormalize(t.gaps_per_length, m_gap, q1_gap, q3_gap);

        t.score = length_weight  * norm_len + hits_weight * norm_hit + gaps_weight * norm_gap;
    }
}

void Score::SelectBestTracks(std::vector<TrackParams>& tracks)
{
    std::map<int, TrackParams> best;

    for (const auto& t : tracks) {
        if (t.score <= score_threshold) continue;

        auto key = t.sliceID;

        if (!best.count(key) || t.score > best[key].score)
            best[key] = t;
    }

    _scoreTracks.reserve(best.size());
    for (auto& pair : best) {
        _scoreTracks.push_back(pair.second);
    }
}

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

double Score::RobustNormalize(double value, double median, double q1, double q3)
{
    double iqr = q3 - q1;
    if (iqr == 0.0) return 0.0;
    return (value - median) / iqr;
}

///////////////////////////////////////////////////////////////////////////////////////

AngleCut::AngleCut(double xz_cut, double yz_cut, double xy_cut):
    xz_cut(xz_cut), yz_cut(yz_cut), xy_cut(xy_cut) 
{ }

std::vector<TrackParams>& AngleCut::operator()(std::vector<TrackParams>& tracks)
{
    _angleTracks.clear();

    double dx, dy, dz;
    double angle_xz, angle_yz, angle_yx;

    for (const auto& t : tracks) {
        dx = t.track.EndX - t.track.StartX;
        dy = t.track.EndY - t.track.StartY;
        dz = t.track.EndZ - t.track.StartZ;

        angle_xz = std::atan2(dz, dx) * 180.0 / M_PI;
        angle_yz = std::atan2(dz, dy) * 180.0 / M_PI;
        angle_yx = std::atan2(dx, dy) * 180.0 / M_PI;

        if ((angle_xz <= xz_cut) || (angle_xz >= (180.0 - xz_cut))) continue;
        if ((angle_yz <= yz_cut) || (angle_yz >= (180.0 - yz_cut))) continue;
        if (std::abs(angle_yx) >= (180.0 - xy_cut))                 continue;

        _angleTracks.push_back(t);
    }

    return _angleTracks;
}

///////////////////////////////////////////////////////////////////////////////////////

void ChooseTrackHits::operator()(
    const std::vector<Hit>& hits, const std::vector<TrackParams>& tracks
)
{
    _selected_hits.clear();
    _selected_tracks.clear();

    std::unordered_set<int> track_slice_ids;
    track_slice_ids.reserve(tracks.size());

    for (const auto& t : tracks) {
        _selected_tracks.push_back(t.track);
        track_slice_ids.insert(t.sliceID);
    }

    _selected_hits.reserve(hits.size());
    for (const auto& h : hits) if (track_slice_ids.count(h.hitSet_id)) _selected_hits.push_back(h);
}

std::vector<Hit>& ChooseTrackHits::GetHits()
{
    return _selected_hits;
}

std::vector<Track>& ChooseTrackHits::GetTracks()
{
    return _selected_tracks;
}

///////////////////////////////////////////////////////////////////////////////////////

SeparateHitsWithTracksRestActivity::SeparateHitsWithTracksRestActivity(double epsilon):
    epsilon(epsilon)
{ }

void SeparateHitsWithTracksRestActivity::operator()(
    std::vector<Hit>& hits, std::vector<Track>& tracks
)
{
    _restHits.clear();

    for (const auto& track : tracks) {
        auto hit_indices = FilterHitsAroundTrack(hits, track);
        
        for (int idx : hit_indices) hits[idx].used = true;
    }
    
    for (const auto& h : hits) if (!h.used) _restHits.push_back(h);
}

std::vector<int> SeparateHitsWithTracksRestActivity::FilterHitsAroundTrack(
    const std::vector<Hit>& hits, const Track& track
)
{
    std::vector<int> indices;

    for (size_t i = 0; i < hits.size(); ++i) {
        if (hits[i].used)                       continue;
        if (hits[i].hitSet_id != track.sliceID) continue;

        double y1 = (hits[i].view == 1) ? track.StartX : track.StartY;
        double y2 = (hits[i].view == 1) ? track.EndX   : track.EndY;

        double d = DistanceToTrack(track.StartZ, track.EndZ, y1, y2, hits[i].plane, hits[i].cell);
        if (d <= epsilon) indices.push_back(i);
    }

    std::sort(indices.begin(), indices.end(),
              [&](int a, int b){ return hits[a].plane < hits[b].plane;});
    
    return indices;
}

std::vector<Hit>& SeparateHitsWithTracksRestActivity::GetHits()
{
    return _restHits;
}

///////////////////////////////////////////////////////////////////////////////////////

RestActivityCut::RestActivityCut(double min_hits):
    min_hits(min_hits)
{ }

void RestActivityCut::operator()(std::vector<Hit>& hits)
{
    _restHits.clear();
    
    std::map<int, int> group_counts;
    std::map<int, int> hits_counts;

    for (const auto& h : hits) group_counts[h.hitSet_id]++;

    for (const auto& h : hits) {
        if (group_counts[h.hitSet_id] > min_hits) {
            _restHits.push_back(h);
            hits_counts[h.hitSet_id] = group_counts[h.hitSet_id];
        }
    }

    int tot_groups = group_counts.size();
    int passed     = hits_counts.size();

    std::cout << "Hits:\n";
    std::cout << "  All: " << tot_groups << '\n';
    std::cout << "  Events with >" << min_hits << " hits: "
              << passed << " ("
              << std::fixed << std::setprecision(1)
              << 100.0 * passed / tot_groups << "%)\n";
}

std::vector<Hit>& RestActivityCut::GetHits()
{
    return _restHits;
}

///////////////////////////////////////////////////////////////////////////////////////
