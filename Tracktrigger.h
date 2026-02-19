#include "fhiclcpp/ParameterSet.h"
#include "DDTBaseDataProducts/DAQHit.h"
#include "DDTBaseDataProducts/Track3D.h"
#include "art/Framework/Principal/Event.h"
#include <vector>
#include <cmath>
#include <string>
#include <map>
#include <unordered_set>
#include <iomanip>
#include <iostream>
#include <algorithm>

///////////////////////////////////////////////////////////////////////////////////////

double DistanceToTrack(double x1, double x2, double y1, double y2, double x, double y);

struct Hit {
    int hitSet_id;
    int plane;
    int cell;
    int adc;
    long long tdc;
    int view;
    bool used;

    Hit(const novaddt::DAQHit& h, int slice_id);
};

struct Track {
    int fView;
    double StartX, StartY, StartZ;
    double EndX, EndY, EndZ;
    int sliceID;
    int hitSet_id;
   
    Track() = default;
    Track(const novaddt::Track3D& t, int slice_id);
};

struct TrackParams {
    int sliceID;

    double gaps_per_length{};
    int total_hits{};
    int total_gaps{};
    double track_length{};
    int max_gap_size{};

    double score{};
    Track track;
};

///////////////////////////////////////////////////////////////////////////////////////

class ParameterCuts {
public:
    ParameterCuts(
        double gaps_min,       double gaps_max,                                                         
        double length_min,     double length_max,
        double max_gap,        double epsilon,                                                          
        double total_hits_min, double total_hits_max
    );

    std::vector<TrackParams>& operator()(
        const std::vector<Hit>& hits, const std::vector<Track>& tracks
    );

private:
    TrackParams CalculateTrackParams(
        const std::vector<Hit>& hits, const Track& track
    );

    bool FilterTrack(const TrackParams& t);

    double gaps_min;                                                         
    double gaps_max;                                                         
    double length_min;                                                       
    double length_max;
    double max_gap;                                                       
    double epsilon;                                                          
    double total_hits_min;                                                   
    double total_hits_max;
    
    std::vector<TrackParams> _paramTracks;
};

///////////////////////////////////////////////////////////////////////////////////////

class Score {
public:
    Score(
        double length_weight, double hits_weight,
        double gaps_weight,   double score_threshold
    );

    std::vector<TrackParams>& operator()(std::vector<TrackParams>& tracks);

private:
    void SelectBestTracks(std::vector<TrackParams>& tracks);

    void ScoreTracks(std::vector<TrackParams>& tracks);

    double Quantile(const std::vector<double>& sorted, double q);
    
    double RobustNormalize(double value, double median, double q1, double q3);

    double length_weight;
    double hits_weight;
    double gaps_weight;
    double score_threshold;
    std::vector<TrackParams> _scoreTracks;
};

///////////////////////////////////////////////////////////////////////////////////////

class AngleCut {
public:
    AngleCut(double xz_cut, double yz_cut, double xy_cut);

    std::vector<TrackParams>& operator() (std::vector<TrackParams>& tracks);

private:
    double xz_cut;
    double yz_cut;
    double xy_cut;
    std::vector<TrackParams> _angleTracks;
};

///////////////////////////////////////////////////////////////////////////////////////

class ChooseTrackHits {
public:
    void operator()(
        const std::vector<Hit>& hits, const std::vector<TrackParams>& tracks
    );

    std::vector<Hit>& GetHits();

    std::vector<Track>& GetTracks();

private:
    std::vector<Hit>   _selected_hits;
    std::vector<Track> _selected_tracks;
};

///////////////////////////////////////////////////////////////////////////////////////

class SeparateHitsWithTracksRestActivity {
public:
    SeparateHitsWithTracksRestActivity(double epsilon);

    void operator()(
        std::vector<Hit>& hits, std::vector<Track>& tracks
    );

    std::vector<Hit>& GetHits();

private:
    std::vector<int> FilterHitsAroundTrack(
        const std::vector<Hit>& hits, const Track& track
    );

    double epsilon;
    std::vector<Hit> _restHits;
};

///////////////////////////////////////////////////////////////////////////////////////

class RestActivityCut {
public:
    RestActivityCut(double min_hits);

    void operator()(std::vector<Hit>& hits);

    std::vector<Hit>& GetHits();

private:
    double min_hits;
    std::vector<Hit> _restHits;
};

///////////////////////////////////////////////////////////////////////////////////////

class Trigger {
public:
    Trigger(const fhicl::ParameterSet& p);

    bool run_algorithm(
        std::vector<Hit>& hits,std::vector<Track>& tracks, const art::Event& event
    );

private:
    ParameterCuts                      _ParametersCuts;
    Score                              _Score;
    AngleCut                           _AngleCut;
    ChooseTrackHits                    _ChooseTrackHits;
    SeparateHitsWithTracksRestActivity _SeparateHitsWithTracksRestActivity;
    RestActivityCut                    _RestActivityCut;
};

///////////////////////////////////////////////////////////////////////////////////////
