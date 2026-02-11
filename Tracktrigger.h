#include <vector>
#include <cmath>
#include <string>
#include <map>

double DistanceToTrack(
    double x1, double x2,
    double y1, double y2,
    double x, double y);

/*struct Hit {
    int hitSet_id;
    int plane;
    int cell;
    int adc;
    long long tdc;
    int view;
    bool used;
};

struct Track {
    int fView;
    double StartX, StartY, StartZ;
    double EndX, EndY, EndZ;
    int sliceID;
    int hitSet_id;
};
*/

// ParameterSet

class ParameterSet {
public:
    ParameterSet(std::vector<std::pair<std::string,double>>& input);
    double get(const std::string& key) const;

private:
    std::map<std::string,double> params;
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
    ParameterCuts(ParameterSet const& p);

    const std::vector<TrackParams>& operator()(
        const std::vector<Hit>& hits,
        const std::vector<Track>& tracks);

private:
    TrackParams CalculateTrackParams(
        const std::vector<Hit>& hits,
        const Track& track);

    bool FilterTrack(const TrackParams& t);

    ParameterSet const& pset;
    std::vector<TrackParams> table;
};

// Score

class Score {
public:
    Score(ParameterSet const& p);

    const std::vector<TrackParams>& operator()(
        const std::vector<TrackParams>& in);

private:
    ParameterSet const& pset;
    std::vector<TrackParams> table;
};

// AngleCut

class AngleCut {
public:
    AngleCut(ParameterSet const& p);

    const std::vector<TrackParams>& operator()(
        const std::vector<TrackParams>& in);

private:
    ParameterSet const& pset;
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
    SeparateHitsWithTracksRestActivity(ParameterSet const& p);

    void Filter(
        const std::vector<Hit>& hits,
        const std::vector<Track>& tracks);

    const RestHits& get_rest_hits() const;

private:
    ParameterSet const& pset;
    RestHits rest_hits;
};

// RestActivityCut

class RestActivityCut {
public:
    RestActivityCut(ParameterSet const& p);

    void Filter(const RestHits& hits);

    const RestHits& get_rest_hits() const;

private:
    ParameterSet const& pset;
    RestHits filtered_hits;
};

// Trigger

class Trigger {
public:
    Trigger(ParameterSet const& p);

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

