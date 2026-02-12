#include "Tracktrigger.h"
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
if (Rabs == 0.0) return 1e9;

double Lx = (Rx*rx + Ry*ry) / Rabs;
double Ly2 = rx*rx + ry*ry - Lx*Lx;
double Ly = (Ly2 > 0.0) ? std::sqrt(Ly2) : 0.0;

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
std::vector<Track>& tracks)
{
const auto& paramTracks =
    _ParametersCuts(hits, tracks);

const auto& scoredTracks =
    _Score(paramTracks);

const auto& angleTracks =
    _AngleCut(scoredTracks);

_ChooseTrackHits.Filter(hits, angleTracks);

_SeparateHitsWithTracksRestActivity.Filter(
    _ChooseTrackHits.get_hits(),
    _ChooseTrackHits.get_tracks());

_RestActivityCut.Filter(
_SeparateHitsWithTracksRestActivity.get_rest_hits());

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

const std::vector<TrackParams>&
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

const std::vector<TrackParams>&
Score::operator()(const std::vector<TrackParams>& in)
{
table.clear();

for (auto t : in)
{
    double score =
        length_weight * t.track_length +
        hits_weight   * t.total_hits  -
        gaps_weight   * t.total_gaps;

    if (score >= score_threshold)
    {
        t.score = score;
        table.push_back(t);
    }
}

return table;
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

    double angle_xz = std::atan2(dx, dz);
    double angle_yz = std::atan2(dy, dz);
    double angle_xy = std::atan2(dy, dx);

    if (std::abs(angle_xz) > xz_cut)
        continue;
    if (std::abs(angle_yz) > yz_cut)
        continue;
    if (std::abs(angle_xy) > xy_cut)
        continue;

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

for (const auto& t : tracks)
{
    selected_tracks.push_back(t.track);
    for (const auto& h : hits)
        selected_hits.push_back(h); 
}
}


const std::vector<Hit>&
ChooseTrackHits::get_hits() const
{
return selected_hits;
}

const std::vector<Track>&
ChooseTrackHits::get_tracks() const
{
return selected_tracks;
}

// SeparateHitsWithTracksRestActivity

SeparateHitsWithTracksRestActivity::
SeparateHitsWithTracksRestActivity(double epsilon):
epsilon(epsilon)
{
}

void SeparateHitsWithTracksRestActivity::Filter(
const std::vector<Hit>& hits,
const std::vector<Track>&)
{
rest_hits.clear();

for (const auto& h : hits)
    rest_hits.push_back(h);
}

const RestHits&
SeparateHitsWithTracksRestActivity::
get_rest_hits() const
{
return rest_hits;
}

// RestActivityCut

RestActivityCut::RestActivityCut(double min_hits):
min_hits(min_hits)
{
}

void RestActivityCut::Filter(
const RestHits& hits)
{
filtered_hits.clear();

if (hits.size() >= min_hits)
{
    filtered_hits = hits;
}
}

const RestHits&
RestActivityCut::get_rest_hits() const
{
return filtered_hits;
}

