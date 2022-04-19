//
// Created by Cloud on 2020/8/21.
//

#ifndef POLYFIT_CLION_VERDICTDB_H
#define POLYFIT_CLION_VERDICTDB_H

#include <vector>
#include <cstdlib>
#include <cmath>
#include <map>
#include <chrono>
#include "utils.h"

/**
 * implement VerdictDB's variational subsampling (for simple query)
 * @param dataset: the original dataset
 * @param ns: the total size of records to sampled, by default, by default it's square root of original dataset's size
 * @param b: the number of subsamples, by default it's 100
 */
template <typename Tk>
void VariationalSampling(const std::vector<Tk> dataset, int ns, int b, std::map<int, std::vector<Tk>> &subsamples){

    srand((unsigned)time(NULL));
    int size_data = dataset.size();
    double sample_threshold = double(ns) / double(size_data);

    typedef std::map<int, std::vector<Tk>> MapType;
    subsamples.clear();

    for(int i = 0; i < size_data; i++){
        double sample_prob = ((double) rand() / (RAND_MAX));
        if(sample_prob <= sample_threshold){
            int sid = rand() % b;
            // insert if not exist
            typename MapType::iterator lb = subsamples.lower_bound(sid);
            if(lb != subsamples.end() && !(subsamples.key_comp()(sid, lb->first))){
                // key already exists
                lb->second.push_back(dataset[i]);
            } else {
                // the key does not exist, add it to the map, so it can avoid another lookup
                std::vector<Tk> v{dataset[i]};
                subsamples.insert(lb, typename MapType::value_type(sid, v));    // Use lb as a hint to insert,
            }
        }
    }
}

template <typename Tk>
void VariationalSampling2D(const std::vector<Tk> keys1, const std::vector<Tk> keys2, int ns, int b, std::map<int, std::vector<pair<Tk, Tk>>> &subsamples){

    srand((unsigned)time(NULL));
    int size_data = keys1.size();
    double sample_threshold = double(ns) / double(size_data);

    typedef std::map<int, std::vector<pair<Tk, Tk>>> MapType;
    subsamples.clear();

    for(int i = 0; i < size_data; i++){
        double sample_prob = ((double) rand() / (RAND_MAX));
        if(sample_prob <= sample_threshold){
            int sid = rand() % b;
            // insert if not exist
            typename MapType::iterator lb = subsamples.lower_bound(sid);
            if(lb != subsamples.end() && !(subsamples.key_comp()(sid, lb->first))){
                // key already exists
                lb->second.push_back(pair<Tk, Tk>(keys1[i], keys2[i]));
            } else {
                // the key does not exist, add it to the map, so it can avoid another lookup
                std::vector<pair<Tk, Tk>> v{pair<Tk, Tk>(keys1[i], keys2[i])};
                subsamples.insert(lb, typename MapType::value_type(sid, v));    // Use lb as a hint to insert,
            }
        }
    }
}

/**
 * performing count aggregate on single subsample
 */
template <typename Tk, typename Tv>
void CountAggSingle(const std::vector<Tk> &subsample, double L, double U, Tv &result){
    result = 0;
    for(int i = 0 ; i < subsample.size(); i++){
        if(subsample[i] >= L and subsample[i] <= U){
            result += 1;
        }
    }
}

template <typename Tk, typename Tv>
void CountAggSingle2D(const std::vector<pair<Tk, Tk>> &subsample, double L1, double L2, double U1, double U2, Tv &result){
    result = 0;
    for(int i = 0 ; i < subsample.size(); i++){
        if(subsample[i].first >= L1 and subsample[i].first <= U1 and subsample[i].second >= L2 and subsample[i].second <= U2){
            result += 1;
        }
    }
}

// sort ascending according to the first value
bool sortpair(const pair<double,int> &a, const pair<double,int> &b){
    return (a.first < b.first);
}

/**
 * perform count aggregate on all provided queries, this one do not have refinement!
 * @tparam Tk
 * @tparam Tv
 * @param queryset_low
 * @param queryset_up
 * @param results
 * @param key_v
 * @param RealResultPath
 * @param ns: total subsample size
 * @param b: number of subsamples
 * @param confidence : integer, confidence %
 * @return
 */
 template<typename Tk, typename Tv>
QueryResult VerdictDBCountAggregate(const std::vector<Tk> &queryset_low, const std::vector<Tk> &queryset_up, std::vector<Tv> &results,
                                    const std::vector<Tk> &key_v, std::string RealResultPath, int ns = -1, int b = -1, double confidence = 0.9){

    if(ns == -1) ns = int(sqrt(key_v.size())); // use the default sampling size mentioned by VerdictDB
    if(b == -1) b = min(100, ns / 50); // at least 50 elements per group

    double sample_threshold = double(ns) / double(key_v.size());
    double amplified_factor = 1.0 / sample_threshold;

    // perform variational sampling first
    std::map<int, std::vector<Tk>> subsamples;
    VariationalSampling(key_v, ns, b, subsamples);
    results.clear();

    auto t0 = std::chrono::steady_clock::now();

    vector<pair<double, int>> subsample_agg_size(subsamples.size()); // estimated count and subsample size

    // aggregate on subsamples
    for(int i = 0; i < queryset_low.size(); i++){
        typename std::map<int, std::vector<Tk>>::iterator it;
        Tv count = 0;
        Tv total_count = 0;
        double subsample_agg = 0;
        subsample_agg_size.clear();

        // first round, calculate original sample aggregate and each subsample aggregate
        for(it = subsamples.begin(); it != subsamples.end(); it++){
            CountAggSingle(it->second, queryset_low[i], queryset_up[i], count);
            subsample_agg = count * (key_v.size() / it->second.size()); // gj
            subsample_agg_size.push_back(pair<double, int>(subsample_agg, it->second.size()));
            total_count += (count);
        }

        // second round, calculate confidence interval
        total_count *= amplified_factor; // g0

        // sort pair according to the modified subsample aggregate value, ascending
        sort(subsample_agg_size.begin(), subsample_agg_size.end(), sortpair);

        double lower_side = (1.0 - confidence) / 2;
        double upper_side = confidence + (1.0 - confidence) / 2;
        int lower_index = lower_side * subsample_agg_size.size();
        int upper_index = upper_side * subsample_agg_size.size();
        double to_minus_lower = (total_count - subsample_agg_size[lower_index].first) * sqrt(subsample_agg_size[lower_index].second / key_v.size());
        double to_minus_upper = (total_count - subsample_agg_size[upper_index].first) * sqrt(subsample_agg_size[upper_index].second / key_v.size());
        double lower_threshold = total_count - to_minus_lower;
        double upper_threshold = total_count - to_minus_upper;

        // calculate average inbetween
        int average = 0;
        int valid_count = 0;
        for(int j = 0; j < subsample_agg_size.size(); j++){
            if(subsample_agg_size[j].first >= lower_threshold || subsample_agg_size[j].first <= upper_threshold){
                valid_count += 1;
                average += subsample_agg_size[j].first;
            }
        }
        average /= valid_count;
        results.push_back(average);
    }

    auto t1 = std::chrono::steady_clock::now();
    auto average_time = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() / queryset_low.size();
    auto total_time = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();

    double MEabs, MErel;
    MeasureAccuracy(results, RealResultPath, MEabs, MErel);

    QueryResult query_result;
    query_result.average_query_time = average_time;
    query_result.total_query_time = total_time;
    query_result.measured_absolute_error = MEabs;
    query_result.measured_relative_error = MErel;

    return query_result;
}

template<typename Tk, typename Tv>
QueryResult VerdictDBCountAggregate2D(const std::vector<Tk> &keys1, const std::vector<Tk> &keys2, const std::vector<Tk> &queryset_low1,
                                      const std::vector<Tk> &queryset_low2, const std::vector<Tk> &queryset_up1, const std::vector<Tk> &queryset_up2,
                                      std::vector<Tv> &results, std::string RealResultPath, int ns = -1, int b = -1, double confidence = 0.9) {

    if(ns == -1) ns = int(sqrt(keys1.size())); // use the default sampling size mentioned by VerdictDB
    if(b == -1) b = min(100, ns / 50); // at least 50 elements per group
    double sample_threshold = double(ns) / double(keys1.size());
    double amplified_factor = 1.0 / sample_threshold;

    // perform variational sampling first
    std::map<int, std::vector<pair<Tk, Tk>>> subsamples;
    VariationalSampling2D(keys1, keys2, ns, b, subsamples);
    results.clear();

    auto t0 = std::chrono::steady_clock::now();

    vector<pair<double, int>> subsample_agg_size(subsamples.size()); // estimated count and subsample size

    // aggregate on subsamples
    for(int i = 0; i < queryset_low1.size(); i++){
        typename std::map<int, std::vector<pair<Tk, Tk>>>::iterator it;
        Tv count = 0;
        Tv total_count = 0;
        double subsample_agg = 0;
        subsample_agg_size.clear();

        // first round, calculate original sample aggregate and each subsample aggregate
        for(it = subsamples.begin(); it != subsamples.end(); it++){
            CountAggSingle2D(it->second, queryset_low1[i], queryset_low2[i], queryset_up1[i], queryset_up2[i], count);
            subsample_agg = count * (keys1.size() / it->second.size()); // gj
            subsample_agg_size.push_back(pair<double, int>(subsample_agg, it->second.size()));
            total_count += (count);
        }

        // second round, calculate confidence interval
        total_count *= amplified_factor; // g0

        // sort pair according to the modified subsample aggregate value, ascending
        sort(subsample_agg_size.begin(), subsample_agg_size.end(), sortpair);

        double lower_side = (1.0 - confidence) / 2;
        double upper_side = confidence + (1.0 - confidence) / 2;
        int lower_index = lower_side * subsample_agg_size.size();
        int upper_index = upper_side * subsample_agg_size.size();
        double to_minus_lower = (total_count - subsample_agg_size[lower_index].first) * sqrt(subsample_agg_size[lower_index].second / keys1.size());
        double to_minus_upper = (total_count - subsample_agg_size[upper_index].first) * sqrt(subsample_agg_size[upper_index].second / keys1.size());
        double lower_threshold = total_count - to_minus_lower;
        double upper_threshold = total_count - to_minus_upper;

        // calculate average inbetween
        int average = 0;
        int valid_count = 0;
        for(int j = 0; j < subsample_agg_size.size(); j++){
            if(subsample_agg_size[j].first >= lower_threshold || subsample_agg_size[j].first <= upper_threshold){
                valid_count += 1;
                average += subsample_agg_size[j].first;
            }
        }
        average /= valid_count;
        results.push_back(average);
    }

    auto t1 = std::chrono::steady_clock::now();
    auto average_time = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() / queryset_low1.size();
    auto total_time = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();

    double MEabs, MErel;
    MeasureAccuracy(results, RealResultPath, MEabs, MErel);

    QueryResult query_result;
    query_result.average_query_time = average_time;
    query_result.total_query_time = total_time;
    query_result.measured_absolute_error = MEabs;
    query_result.measured_relative_error = MErel;

    return query_result;

}
#endif //POLYFIT_CLION_VERDICTDB_H
