//
// Created by Cloud on 2020/9/4.
//

#ifndef POLYFIT_CLION_ICTAI1997_H
#define POLYFIT_CLION_ICTAI1997_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <stx/btree.h>
#include "utils.h"

using namespace std;

template <typename Tk, typename Tv>
class ICTAI1997 {
public:
    ICTAI1997(const vector<Tk> &key_v, const vector<Tv> &value_v){
        BottomUp(key_v, value_v);
    }

    // y = (Sxy / Sx^2)(x - x_avg) + y_avg
    // Sxy = ( sum(xy) - n * x_avg * y_avg ) / n-1
    // Sxx = ( sum(xi^2) - n(x_avg)^2 ) / n-1    square of standard deviation
    struct Segment{
        double slope;
        double intercept;

        Segment(double slope, double intercept): slope(slope), intercept(intercept) {}

        inline Tv operator()(const Tk &k) const {
            auto pred = slope * k + intercept;
            return pred;
        }

        friend inline bool operator<(const Segment &s, const Tk &k) { return s.key < k; }
        friend inline bool operator<(const Tk &k, const Segment &s) { return k < s.key; }
    };

    void BottomUp(const vector<Tk> &key_v, const vector<Tv> &value_v){

        vector<Segment> segs;
        vector<double> nor_errs;
        vector<pair<int, int>> intervals;

        // create initial finest segments, using n/3 segments
        int start_index = 0, end_index = 0;
        while(end_index < key_v.size()){
            end_index = start_index + 6 <= key_v.size()? start_index + 3 : key_v.size(); // should be 3
            Segment s = ClassicRegression(key_v, value_v, start_index, end_index);
            segs.push_back(s);
            double normalized_error = NormalizedError(key_v, value_v, start_index, end_index, s);
            nor_errs.push_back(normalized_error);
            intervals.push_back(pair<int, int>(start_index, end_index));
            start_index = end_index;
        }

        // calculate the merget cost
        vector<double> Bk_over_iterations;

        // find the lowest merge cost
        double min_bk;
        int min_index = 0, iteration_count = 0;
        double sum_error_square = -1, sum_error = -1, min_merged_nor_err;

        while(segs.size() > 1){

            CalculateMergedBk(key_v, value_v, intervals, nor_errs, min_bk, min_index, min_merged_nor_err, sum_error_square, sum_error);

            //int min_index = std::min_element(merge_costs.begin(), merge_costs.end()) - merge_costs.begin();
            if(Bk_over_iterations.size() > 0 && min_bk > Bk_over_iterations.back()){
                segments = segs;
                segment_keys.clear();
                for(int i = 0; i < intervals.size(); i++){
                    segment_keys.push_back(key_v[intervals[i].first]);
                }
                seg_idx.clear();
                seg_idx.bulk_load(segment_keys.begin(), segment_keys.end());
                break;
            }

            Bk_over_iterations.push_back(min_bk);

            // merge the segments
            Segment merged_s = ClassicRegression(key_v, value_v, intervals[min_index].first, intervals[min_index+1].second);

            segs[min_index] = merged_s;
            segs.erase(segs.begin() + min_index + 1);

            sum_error_square = sum_error_square - nor_errs[min_index] * nor_errs[min_index] - nor_errs[min_index+1] * nor_errs[min_index+1] + min_merged_nor_err * min_merged_nor_err;
            sum_error = sum_error - nor_errs[min_index] - nor_errs[min_index+1] + min_merged_nor_err;

            nor_errs[min_index] = min_merged_nor_err;
            nor_errs.erase(nor_errs.begin() + min_index + 1);

            intervals[min_index] = pair<int, int>(intervals[min_index].first, intervals[min_index+1].second);
            intervals.erase(intervals.begin() + min_index + 1);

            if(iteration_count % 100 == 0)
                cout << "iteration: " << iteration_count << " segment amount: " << segs.size() << endl;
            iteration_count++;
        }
    }

    QueryResult CountPrediction(vector<Tk> &queryset_low, vector<Tk> &queryset_up, vector<Tv> &results, const vector<Tk> &key_v, string RealResultPath){
        results.clear();
        int count_refinement = 0;
        auto t0 = chrono::steady_clock::now();

        for(int i = 0; i < queryset_low.size(); i++){
            auto pred_l = segments[prev(seg_idx.upper_bound(queryset_low[i]))->second](queryset_low[i]); // pred pos L
            auto pred_u = segments[prev(seg_idx.upper_bound(queryset_up[i]))->second](queryset_up[i]); // pred pos U
            Tv result = pred_u - pred_l;
            results.push_back(result);
        }

        auto t1 = chrono::steady_clock::now();

        QueryResult query_result;

        // measure query performance
        auto average_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count() / queryset_low.size();
        auto total_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
        double MEabs, MErel;
        MeasureAccuracy(results, RealResultPath, MEabs, MErel);

        query_result.average_query_time = average_time;
        query_result.total_query_time = total_time;
        query_result.measured_absolute_error = MEabs;
        query_result.measured_relative_error = MErel;
        query_result.hit_count = queryset_low.size() - count_refinement;
        query_result.model_amount = this->segments.size();
        query_result.tree_paras = this->seg_idx.CountParametersNewPrimary(false);
        query_result.total_paras = this->segments.size() * 3 + query_result.tree_paras;

        return query_result;
    }

    // y = (Sxy / Sx^2)(x - x_avg) + y_avg
    // Sxy = ( sum(xy) - n * x_avg * y_avg ) / n-1
    // Sxx = ( sum(xi^2) - n(x_avg)^2 ) / n-1    square of standard deviation
    Segment ClassicRegression(const vector<Tk> &key_v, const vector<Tv> &value_v, int start_index, int end_index){
        double error = 0;

        double n = end_index - start_index; // [ , )
        double Sxy, Sxx;
        double sum_xy = 0, sum_xi2 = 0, sum_x = 0, sum_y = 0, x_avg, y_avg;
        for(int i = start_index; i < end_index; i++){
            sum_xy += key_v[i] * value_v[i];
            sum_xi2 += key_v[i] * key_v[i];
            sum_x += key_v[i];
            sum_y += value_v[i];
        }
        x_avg = sum_x / n;
        y_avg = sum_y / n;
        Sxy = (sum_xy - n * x_avg * y_avg) / (n-1);
        Sxx = (sum_xi2 - n * x_avg * x_avg) / (n-1);

        double slope, intercept;
        slope = Sxy / Sxx;
        intercept = y_avg - slope * x_avg;
        Segment s(slope, intercept);
        return s;
    }

    double NormalizedError(const vector<Tk> &key_v, const vector<Tv> &value_v, int start_index, int end_index, Segment s){
        double err, dis;
        for(int i = start_index; i < end_index; i++){
            dis = s(key_v[i]) - value_v[i];
            err += dis * dis;
        }
        err /= (end_index - start_index);
        return err;
    }

    double BalanceError(vector<double> &nor_errs){
        double sum_err = 0, avg_err;
        for(int i = 0; i < nor_errs.size(); i++){
            sum_err += nor_errs[i];
        }
        avg_err = sum_err / nor_errs.size();
        double variance = 0;
        for(int i = 0; i < nor_errs.size(); i++){
            variance += (nor_errs[i] - avg_err) * (nor_errs[i] - avg_err);
        }
        variance /= (nor_errs.size() - 1);
        double std_var_err = sqrt(variance);
        return std_var_err;
    }

    double inline IncrementalBalanceError(const vector<double> &nor_errs, const int &merged_idx1, const int &merged_idx2,
                                          const double &merged_err, const double &err_sum_square, const double &err_sum){
        double new_err_std = 0; // instead of standard deviation
        new_err_std = err_sum_square - nor_errs[merged_idx1] * nor_errs[merged_idx1] - nor_errs[merged_idx2] * nor_errs[merged_idx2] + merged_err * merged_err;
        double new_avg = (err_sum - nor_errs[merged_idx1] - nor_errs[merged_idx2] + merged_err) / (nor_errs.size() - 1);
        new_err_std -= (nor_errs.size() - 1) * new_avg * new_avg;
        new_err_std /= (nor_errs.size() - 2);
        new_err_std = sqrt(new_err_std);
        return new_err_std;
    }

    void CalculateMergedBk(const vector<Tk> &key_v, const vector<Tv> &value_v, vector<pair<int, int>> &intervals, vector<double> &nor_errs,
                           double &min_bk, int &min_index, double &min_merged_err, double &sum_error_square, double &sum_error){
        if(sum_error_square == -1 || sum_error == -1){
            for(int i = 0; i < nor_errs.size(); i++){
                sum_error_square += nor_errs[i] * nor_errs[i];
                sum_error += nor_errs[i];
            }
        }

        bool init_tag = true;
        for(int i = 0; i < nor_errs.size() - 1; i++){
            Segment merged_s = ClassicRegression(key_v, value_v, intervals[i].first, intervals[i+1].second);
            double merged_nor_err = NormalizedError(key_v, value_v, intervals[i].first, intervals[i+1].second, merged_s);
            //vector<double> merged_nor_errs(nor_errs);
            //merged_nor_errs[i] = merged_nor_err;
            //merged_nor_errs.erase(merged_nor_errs.begin() + i + 1);
            double merged_bk = IncrementalBalanceError(nor_errs, i, i+1, merged_nor_err, sum_error_square, sum_error);
            //double merged_bk2 = BalanceError(merged_nor_errs);
            if(init_tag){
                min_bk = merged_bk;
                min_index = i;
                min_merged_err = merged_nor_err;
                init_tag = false;
            } else if (merged_bk < min_bk){
                min_bk = merged_bk;
                min_index = i;
                min_merged_err = merged_nor_err;
            }
        }
    }

    vector<Segment> segments;
    vector<Tk> segment_keys;
    stx::btree<Tk, int> seg_idx; // used to index the segments
    Tv t_abs; // the absolute error threshold, default 100
};

#endif //POLYFIT_CLION_ICTAI1997_H
