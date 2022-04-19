//
// Created by Cloud on 2020/8/17.
//

#ifndef POLYFIT_CLION_PGM_H
#define POLYFIT_CLION_PGM_H

#include <vector>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <chrono>
#include <cmath>
#include "pgm_index.hpp"
#include "utils.h"

using namespace std;

// seems the official implementation cannot handle repeatted data
template<typename Tk, typename Tv>
class PGMWrapper{
public:
    // notice the keys should be already sorted
    PGMWrapper(vector<Tk> &keys, size_t absolute_error, double relative_error_threshold){
        this->absolute_error = absolute_error;
        this->t_rel = relative_error_threshold;
        pgm_index = PGMIndex<Tk, 100>(keys, absolute_error);
        cout << "construct PGM finished, total key size: " << keys.size() << endl;
    }

    QueryResult RangeCountQuery(vector<Tk> &queryset_low, vector<Tk> &queryset_up, vector<Tv> &results, vector<Tk> &keys, bool DoRefinement,
            string RealResultPath, bool statistics = true){

        auto t0 = chrono::steady_clock::now();

        results.clear();
        int result, result_low, result_up;
        double max_err_rel = 0; // the estimated maximum possible relative error
        int count_refinement = 0;
        for (int i = 0; i < queryset_low.size(); i++) {
            // calculate the lower key position
            result_low = pgm_index.ApproxSearch(queryset_low[i]);
            result_up = pgm_index.ApproxSearch(queryset_up[i]);
            result = result_up - result_low;
            if (DoRefinement && !RelativeErrorCheck_1D_Count(result, absolute_error, t_rel)) {

                count_refinement++;
                // do refinement
                auto range_l = pgm_index.search(queryset_low[i]);
                auto lo_l = keys.begin() + range_l.lo;
                auto hi_l = keys.begin() + range_l.hi;
                result_low = *std::lower_bound(lo_l, hi_l, queryset_low[i]);

                auto range_u = pgm_index.search(queryset_low[i]);
                auto lo_u = keys.begin() + range_u.lo;
                auto hi_u = keys.begin() + range_u.hi;
                result_low = *std::lower_bound(lo_u, hi_u, queryset_up[i]);

                result = result_up - result_low;
            }
            results.push_back(result);
        }
        auto t1 = chrono::steady_clock::now();
        auto average_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count() / queryset_low.size();
        auto total_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();

        double MEabs, MErel;
        cout << "result size: " << results.size() << endl;
        if(statistics)
            MeasureAccuracy(results, RealResultPath, MEabs, MErel);

        QueryResult query_result;
        query_result.average_query_time = average_time;
        query_result.total_query_time = total_time;
        query_result.measured_absolute_error = MEabs;
        query_result.measured_relative_error = MErel;
        query_result.refinement_count = count_refinement;
        query_result.hit_count = queryset_low.size() - count_refinement;
        query_result.model_amount = pgm_index.segments_count();
        query_result.hit_count = queryset_low.size() - count_refinement;
        query_result.total_bytes = pgm_index.size_in_bytes();
        return query_result;
    }

    size_t absolute_error;
    double t_rel; // relative error threshold
    PGMIndex<Tk, 100> pgm_index;
};


void TryPGM(){

    // dataset
    std::vector<double> K, V;
    LoadTweetDataset(K, V);

    // Construct the PGM-index
    const int epsilon = 128; // space-time trade-off parameter
    // PGMIndex<double, 100> index(data, 120); // this is OK
    PGMIndex<double, 100> index(K, 100); // this is OK

    // Query the PGM-index
    auto t0 = std::chrono::steady_clock::now();

    auto q = 42;
    auto pos = index.ApproxSearch(q);

    auto t1 = std::chrono::steady_clock::now();
    std::cout << "processing time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "ns" << std::endl;

    std::cout << "pos: " << pos << std::endl;
    std::cout << "leaf level size: " << index.segments_count() << std::endl;
}

void TryPGM2(){

    // dataset
    std::vector<double> K, V, QL, QU;
    LoadOSM_Dataset_Queryset_1D_SUBSET(K, V, QL, QU, 2, false); // true for largeInt version

    // Construct the PGM-index
    const int epsilon = 128; // space-time trade-off parameter
    // PGMIndex<double, 100> index(data, 120); // this is OK
    PGMIndex<double, 100> index(K, 100); // this is OK

    // Query the PGM-index
    auto t0 = std::chrono::steady_clock::now();

    auto q = 42;
    auto pos = index.ApproxSearch(q);

    auto t1 = std::chrono::steady_clock::now();
    std::cout << "processing time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() << "ns" << std::endl;

    std::cout << "pos: " << pos << std::endl;
    std::cout << "leaf level size: " << index.segments_count() << std::endl;
}

void TryPGM3(){
    std::vector<int> K, V;

    arma::mat dataset;
//    bool loaded = mlpack::data::Load("/mnt/d/Polyfit_Dataset/TWEET/test.csv", dataset); // deduplicated, works
//    bool loaded = mlpack::data::Load("/mnt/d/Polyfit_Dataset/TWEET/1M_all.csv", dataset); // original, do not work
    bool loaded = mlpack::data::Load("/mnt/d/Polyfit_Dataset/TWEET/1M_keys_toLargeInt.csv", dataset); // transformed to large int, works

    arma::rowvec trainingset = dataset.row(0);
    arma::rowvec responses = dataset.row(dataset.n_rows - 1);
    RowvecToVector(trainingset, K);
//    for (int i = 0; i < 3; i++) {
//        K.push_back(trainingset[i]);
//    }

    std::cout << "K size: " << K.size() << " " << K[0] << " " << K[1]  << " " << K[2] << endl;

    PGMIndex<int, 100> index(K,50);
    std::cout << "leaf level size: " << index.segments_count() << std::endl;

    auto q = 4195000;
    auto range = index.search(q);
    auto lo = K.begin() + range.lo;
    auto hi = K.begin() + range.hi;
    std::cout << *std::lower_bound(lo, hi, q) << endl;
    std::cout << "approxmate pos: " << range.pos << endl;
    std::cout << "lower pos: " << range.lo << endl;
    std::cout << "upper pos: " << range.hi << endl;
}

#endif //POLYFIT_CLION_PGM_H
