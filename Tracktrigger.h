#include "fhiclcpp/ParameterSet.h"
#include "DDTBaseDataProducts/DAQHit.h"
#include "DDTBaseDataProducts/Track3D.h"
#include <vector>
#include <cmath>
#include <string>
#include <map>

double DistanceToTrack(
    double x1, double x2,
    double y1, double y2,
    double x, double y);

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

// Track parameters

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

// ParameterCuts

class ParameterCuts {
public:
    ParameterCuts(double gaps_min, 
                  double gaps_max,
                  double length_min,
                  double length_max,
                  double max_gap,
                  double epsilon,
                  double total_hits_min,
                  double total_hits_max);

    const std::vector<TrackParams>& operator()(
        const std::vector<Hit>& hits,
        const std::vector<Track>& tracks);

private:
    TrackParams CalculateTrackParams(
        const std::vector<Hit>& hits,
        const Track& track);

    bool FilterTrack(const TrackParams& t);

    double gaps_min;                                                         
    double gaps_max;                                                         
    double length_min;                                                       
    double length_max;
    double max_gap;                                                       
    double epsilon;                                                          
    double total_hits_min;                                                   
    double total_hits_max;
    
    std::vector<TrackParams> table;
};

// Score

class Score {
public:
    Score(
        double length_weight,
        double hits_weight,
        double gaps_weight,
        double score_threshold);

    const std::vector<TrackParams>& operator()(
        const std::vector<TrackParams>& in);

private:
    double length_weight;
    double hits_weight;
    double gaps_weight;
    double score_threshold;
    std::vector<TrackParams> table;

    };

// AngleCut

class AngleCut {
public:
    AngleCut(double xz_cut, double yz_cut, double xy_cut);

    const std::vector<TrackParams>& operator()(
        const std::vector<TrackParams>& in);

private:
    double xz_cut;
    double yz_cut;
    double xy_cut;
    std::vector<TrackParams> table;
};

// ChooseTrackHits

class ChooseTrackHits {
public:
    void Filter(
        const std::vector<Hit>& hits,
        const std::vector<TrackParams>& tracks);

    const std::vector<Hit>& get_hits() const;
    const std::vector<Track>& get_tracks() const;

private:
    std::vector<Hit> selected_hits;
    std::vector<Track> selected_tracks;
};

// Rest activity

using RestHits = std::vector<Hit>;

class SeparateHitsWithTracksRestActivity {
public:
    SeparateHitsWithTracksRestActivity(double epsilon);

    void Filter(
        const std::vector<Hit>& hits,
        const std::vector<Track>& tracks);

    const RestHits& get_rest_hits() const;

private:
    double epsilon;
    RestHits rest_hits;
};

// RestActivityCut

class RestActivityCut {
public:
    RestActivityCut(double min_hits);

    void Filter(const RestHits& hits);

    const RestHits& get_rest_hits() const;

private:
    double min_hits;
    RestHits filtered_hits;
};

// Trigger

class Trigger {
public:
   // Trigger(ParameterSet const& p);
    Trigger(const fhicl::ParameterSet& p);

    bool run_algorithm(
         std::vector<Hit>& hits,
         std::vector<Track>& tracks);

private:
    ParameterCuts _ParametersCuts;
    Score _Score;
    AngleCut _AngleCut;
    ChooseTrackHits _ChooseTrackHits;
    SeparateHitsWithTracksRestActivity _SeparateHitsWithTracksRestActivity;
    RestActivityCut _RestActivityCut;
};

