//
// Created by Cloud on 2020/8/28.
//

#ifndef POLYFIT_CLION_RADIXSPLINE_H
#define POLYFIT_CLION_RADIXSPLINE_H

#include <vector>
#include "rs/builder.h" // for sorted data
#include "rs/radix_spline.h"
#include "utils.h"

using namespace std;
/**
 * Notice RadixSpline only support uint_32 and uint_64, the keys should be converted!!!
 * @tparam Tk : type of key
 */
class RSWrapper{
public:
    RSWrapper(double t_abs, double t_rel):t_abs(t_abs),t_rel(t_rel){}

    /**
     * build the RadixSpline
     * @param keys
     * @param Tabs : absolute error threshold
     */
    void build(vector<uint64_t> &keys){
        uint64_t min = keys.front();
        uint64_t max = keys.back();
        rs::Builder<uint64_t> rsb(min, max, 18, t_abs);
        for (const auto& key : keys) rsb.AddKey(key);
        rs = rsb.Finalize();
    }

    QueryResult RSCountPrediction(const vector<uint64_t> &queryset_low, const vector<uint64_t> &queryset_up, vector<uint64_t> &results,
                                  const vector<uint64_t> &key_v, bool refinement, string RealResultPath, bool statistics = true){

        results.clear();
        int count_refinement = 0;
        auto t0 = chrono::steady_clock::now();

        for(int i = 0; i < queryset_low.size(); i++){
            auto pred_l = rs.GetEstimatedPosition(queryset_low[i]);
            auto pred_u = rs.GetEstimatedPosition(queryset_up[i]);
            auto result = pred_u - pred_l;
            if(refinement && !RelativeErrorCheck_1D_Count(result, t_abs, t_rel)){
                // refinement, find the exact value for lower key and upper key
                count_refinement++;
                int shift_l_l = max(int(pred_l - t_abs), 0), shift_l_u = min(int(pred_l + t_abs), (int)key_v.size());
                int shift_u_l = max(int(pred_u - t_abs), 0), shift_u_u = min(int(pred_u + t_abs), (int)key_v.size());
                auto l = prev(upper_bound(next(key_v.begin(), shift_l_l), next(key_v.begin(), shift_l_u), queryset_low[i]));
                auto u = prev(upper_bound(next(key_v.begin(), shift_u_l), next(key_v.begin(), shift_u_u), queryset_low[i]));
                result = std::distance(l, u); // since it's count
            }
            results.push_back(result);
        }

        auto t1 = chrono::steady_clock::now();

        QueryResult query_result;
        // measure query performance
        auto average_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count() / queryset_low.size();
        auto total_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
        double MEabs, MErel;

        if(statistics)
            MeasureAccuracy(results, RealResultPath, MEabs, MErel);

        query_result.average_query_time = average_time;
        query_result.total_query_time = total_time;
        query_result.measured_absolute_error = MEabs;
        query_result.measured_relative_error = MErel;
        query_result.hit_count = queryset_low.size() - count_refinement;
        query_result.model_amount = rs.GetSegmentAmount();
        query_result.total_bytes = rs.GetSize();

        return query_result;
    }

    void TestRS(){
        //rs::SearchBound bound = rs.GetSearchBound(8128);
        auto pred = rs.GetEstimatedPosition((42+90)*100000);
        cout << "pred: " << pred << endl;
    }

    double t_abs; // absolute error threshold
    double t_rel; // relative error threshold
    rs::RadixSpline<uint64_t> rs;
};

#endif //POLYFIT_CLION_RADIXSPLINE_H
