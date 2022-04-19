//
// Created by Cloud on 2020/8/26.
//

#ifndef POLYFIT_CLION_POLYFIT_H
#define POLYFIT_CLION_POLYFIT_H

#include <vector>
#include <string>
#include "ilcplex/ilocplex.h"
#include "pgm_index.hpp" // use it for the non-leaf layer (replace STX B+tree)
#include <stx/btree.h>
#include "utils.h"
#include "SimpleRTree.h"

using namespace std;

template<class DATATYPE, class ELEMTYPE, int NUMDIMS, class ELEMTYPEREAL, int TMAXNODES, int TMINNODES, class AGGTYPE, typename Tk, typename Tv, int degree>
class RTree; // forward declare to solve the circular include, alternatively, define Segment in another file such as utils
// https://stackoverflow.com/questions/1281641/circular-c-header-includes
// http://www.cplusplus.com/articles/Gw6AC542/

/**
 * Polyfit 1D
 * @tparam Tk : type of key
 * @tparam Tv : type of value
 * @tparam degree : the highest degree
 */
template<typename Tk, typename Tv, int degree = 1>
class Polyfit{
public:

    Polyfit(Tv t_abs = 100, double t_rel = 0.01): t_abs(t_abs), t_rel(t_rel) {}
    ~Polyfit(){}

    Polyfit(const vector<Tk> &key_v, const vector<Tv> &value_v, Tv t_abs = 100, double t_rel = 0.01): t_abs(t_abs), t_rel(t_rel) {
        this->build(key_v, value_v);
        this->build_non_leaf();
        //this->build_non_leaf_with_pgm();
    }

    /**
     * build the bottm layer, i.e., optimal polynomial functions
     * @param key_v
     * @param value_v
     */
    void build(const vector<Tk> &key_v, const vector<Tv> &value_v){
        int idx_l = 0; // segment start position, inclusive
        int idx_u = 1; // segment end position, exclusive

        vector<double> segment_paras;
        double loss;

        while(idx_u < key_v.size()){
            // find the initial position for upper index within twice absolute error
            idx_u = idx_l+1;
            while(value_v[idx_u]-value_v[idx_l] <= 2*t_abs && idx_u < key_v.size()) idx_u++;

            // what if idx_u is now end of the pos
            SolveMaxlossLP(key_v, value_v, idx_l, idx_u, segment_paras, loss);

            // exponential search for the upper index position
            int segment_size = idx_u - idx_l;
            while(loss <= t_abs && idx_u < key_v.size()){
                segment_size *= 2;
                idx_u = min(idx_l + segment_size, (int)key_v.size());
                SolveMaxlossLP(key_v, value_v, idx_l, idx_u, segment_paras, loss);
            }

            // binary search for the exact upper index position
            int L = idx_l + segment_size / 2;
            int U = idx_u;
            int mid;
            while(L < U - 1){ // as U is exclusive
                mid = (L + U) / 2;
                SolveMaxlossLP(key_v, value_v, idx_l, mid, segment_paras, loss);
                if(loss <= t_abs) L = mid;
                else U = mid - 1;
            }
            idx_u = U;
            SolveMaxlossLP(key_v, value_v, idx_l, idx_u, segment_paras, loss);

            // save this segment
            Segment segment(key_v[idx_l], segment_paras);
            segments.push_back(segment);
            intervals.push_back(pair<int, int>(idx_l, idx_u)); // in intervals, the last element is exclusive

            idx_l = idx_u;
            //cout << "current index_upper: " << idx_u << endl;
        }
    }

    /**
     * build the non-leaf layer, assume the bottom layer segments have been built
     */
    void build_non_leaf(){
        vector<pair<Tk, int>> keys(segments.size());
        for(int i = 0; i < segments.size(); i++) keys[i] = pair<Tk, int>(segments[i].key, i);
        seg_idx.bulk_load(keys.begin(), keys.end());
    }

    void build_non_leaf_with_pgm(){
        vector<Tk> keys(segments.size());
        for(int i = 0; i < segments.size(); i++) keys[i] = segments[i].key;
        PGMIndex<Tk, 4> index(keys, 4);
        auto upper_size = index.size_in_bytes();
        cout << " = = = upper layer size of polyfit using pgm: " << upper_size << endl;
        cout << "= = = total size of polyfit using pgm: " << keys.size()*(degree+2)*4 << endl;
    }

    void build_aggregate_for_max(const std::vector<Tk> &keys, const std::vector<Tv> &values){
        //artree.Reset();
        double seg_max = 0;
        for (int i = 0; i < intervals.size(); i++){
            seg_max = values[intervals[i].first];
            for(int j = intervals[i].first; j < intervals[i].second; j++){
                if(values[j] > seg_max)
                    seg_max = values[j];
            }
            seg_maxs.push_back(seg_max);

            Tk min[1]{keys[intervals[i].first]};
            Tk max[1]{keys[intervals[i].second]-0.001}; // so that it will not include the end

            artree.Insert(min, max, i);
        }
        std::cout << "finsih inserting segments to ARTree for Polyfit" << std::endl;
        artree.GenerateMaxAggregate(artree.m_root, seg_maxs);
        std::cout << "finsih generating Max aggreate in ARTree for Polyfit" << std::endl;
    }

    void testRealMaxFromLeft(const Tk &k, int id, const vector<Tk> &keys, const vector<Tv> &values){
        Tk max_agg = values[intervals[id].first];
        for(int i = intervals[id].first; i < intervals[id].second; i++){
            if(keys[i] > k)
                break;
            if(values[i] > max_agg)
                max_agg = values[i];
        }
        cout << "right part true max value: " << max_agg << endl;
        //return max_agg;
    }

    void testRealMaxFromRight(const Tk &k, int id, const vector<Tk> &keys, const vector<Tv> &values){
        Tk max_agg = values[intervals[id].second-1];
        for(int i = intervals[id].second-1; i >= intervals[id].first; i--){
            if(keys[i] < k)
                break;
            if(values[i] > max_agg)
                max_agg = values[i];
        }
        cout << "left part true max value: " << max_agg << endl;
        //return max_agg;
    }

    Tv MaxRefinement(const Tk &kl, const Tk &ku, const vector<Tk> &keys, const vector<Tv> &values, int segid){
        Tv max_agg = numeric_limits<Tv>::min();
        for(int i = intervals[segid].first; i < intervals[segid].second; i++){
            if(keys[i] > ku)
                break;
            else if(keys[i] >= kl && values[i] > max_agg)
                max_agg = values[i];
        }
        return max_agg;
    }

    QueryResult MaxPrediction(vector<Tk> &queryset_low, vector<Tk> &queryset_up, vector<Tv> &results, const vector<Tk> &key_v,
                              const vector<Tv> &values_v, bool refinement, string RealResultPath, bool statistics = true){
        results.clear();
        double result;
        int count_refinement = 0;
        auto t0 = chrono::steady_clock::now();
        int left_model_id, right_model_id;
        typename SimpleRtreeType::Rect rect;
        for(int i = 0; i < queryset_low.size(); i++){
            //cout << "= = = query " << i << " = = =" << endl;
            rect.m_min[0] = queryset_low[i];
            rect.m_max[0] = queryset_up[i];
            result = artree.ProgressiveMaxForPolyfit(artree.m_root, &rect, seg_maxs, segments, t_abs, left_model_id, right_model_id);
            if(refinement && !RelativeErrorCheck_1D_Max(result, t_abs, t_rel)){
                // refinement
                count_refinement++;
            }
            results.push_back(result);
        }

        auto t1 = chrono::steady_clock::now();

        QueryResult query_result;
        if(!statistics) return query_result;

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
        //query_result.tree_paras = this->seg_idx.CountParametersNewPrimary(false);
        //query_result.total_paras = this->segments.size() * (2 + degree) + query_result.tree_paras;

        return query_result;
    }

//    QueryResult MaxPrediction(vector<Tk> &queryset_low, vector<Tk> &queryset_up, vector<Tv> &results, const vector<Tk> &key_v,
//                              const vector<Tv> &values_v, bool refinement, string RealResultPath, bool statistics = true){
//        results.clear();
//        double result;
//        int count_refinement = 0;
//        auto t0 = chrono::steady_clock::now();
//        int left_model_id, right_model_id;
//        bool refined_label = false;
//        typename SimpleRtreeType::Rect rect;
//        for(int i = 0; i < queryset_low.size(); i++){
//            //cout << "= = = query " << i << " = = =" << endl;
//            rect.m_min[0] = queryset_low[i];
//            rect.m_max[0] = queryset_up[i];
//            result = artree.ProgressiveMaxForPolyfitWithRelErrGuarantee(artree.m_root, &rect, key_v, values_v, seg_maxs, segments, t_abs, t_rel,
//                                                                        left_model_id, right_model_id, intervals, refinement, refined_label);
//            if(refined_label){
//                // refinement
//                count_refinement++;
//
//                refined_label = false;
//            }
//            results.push_back(result);
//        }
//
//        auto t1 = chrono::steady_clock::now();
//
//        QueryResult query_result;
//        if(!statistics) return query_result;
//
//        // measure query performance
//        auto average_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count() / queryset_low.size();
//        auto total_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
//        double MEabs, MErel;
//        MeasureAccuracy(results, RealResultPath, MEabs, MErel);
//
//        query_result.average_query_time = average_time;
//        query_result.total_query_time = total_time;
//        query_result.measured_absolute_error = MEabs;
//        query_result.measured_relative_error = MErel;
//        query_result.hit_count = queryset_low.size() - count_refinement;
//        query_result.model_amount = this->segments.size();
//        //query_result.tree_paras = this->seg_idx.CountParametersNewPrimary(false);
//        //query_result.total_paras = this->segments.size() * (2 + degree) + query_result.tree_paras;
//
//        return query_result;
//    }

    /**
     * CountPrediction, assume the segments have been built
     * @param queryset_low
     * @param queryset_up
     * @param results
     * @param refinement
     */
    QueryResult CountPrediction(vector<Tk> &queryset_low, vector<Tk> &queryset_up, vector<Tv> &results, const vector<Tk> &key_v, bool refinement,
                                string RealResultPath, bool statistics = true){
        results.clear();
        int count_refinement = 0;
        auto t0 = chrono::steady_clock::now();

        for(int i = 0; i < queryset_low.size(); i++){
            auto pred_l = segments[prev(seg_idx.upper_bound(queryset_low[i]))->second](queryset_low[i]); // pred pos L
            auto pred_u = segments[prev(seg_idx.upper_bound(queryset_up[i]))->second](queryset_up[i]); // pred pos U
            Tv result = pred_u - pred_l;
            if(refinement && !RelativeErrorCheck_1D_Count(result, t_abs, t_rel)){
                // refinement, find the exact value for lower key and upper key
                count_refinement++;
                int shift_l_l = max((int)(pred_l - t_abs), 0), shift_l_u = min((int)(pred_l + t_abs), (int)key_v.size());
                int shift_u_l = max((int)(pred_u - t_abs), 0), shift_u_u = min((int)(pred_u + t_abs), (int)key_v.size());
                auto l = prev(upper_bound(next(key_v.begin(), shift_l_l), next(key_v.begin(), shift_l_u), queryset_low[i]));
                auto u = prev(upper_bound(next(key_v.begin(), shift_u_l), next(key_v.begin(), shift_u_u), queryset_low[i]));
                result = std::distance(l, u); // since it's count
            }
            results.push_back(result);
        }

        auto t1 = chrono::steady_clock::now();

        // measure query performance
        auto average_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count() / queryset_low.size();
        auto total_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
        double MEabs, MErel;

        if(statistics)
            MeasureAccuracy(results, RealResultPath, MEabs, MErel);

        QueryResult query_result;
        query_result.average_query_time = average_time;
        query_result.total_query_time = total_time;
        query_result.measured_absolute_error = MEabs;
        query_result.measured_relative_error = MErel;
        query_result.hit_count = queryset_low.size() - count_refinement;
        query_result.refinement_count = count_refinement;
        query_result.model_amount = this->segments.size();
        query_result.tree_paras = this->seg_idx.CountParametersNewPrimary(false);
        query_result.total_paras = this->segments.size() * (2 + degree) + query_result.tree_paras;
        query_result.total_bytes = query_result.total_paras * 4;

        return query_result;
    }

    void test_query(){
        auto it = seg_idx.upper_bound(9000000);
        cout << "query key for 90: " << std::prev(it)->first << " " << std::prev(it)->second << endl; // last one, 797
        it = seg_idx.upper_bound(-9000000);
        cout << "query key for -90: " << std::prev(it)->first << " " << std::prev(it)->second << endl; // first one, 0

        int key = 4200000;
        auto segment_id = prev(seg_idx.upper_bound(key))->second;
        auto pred = segments[segment_id](key);
        cout << "segment id for 42: " << segment_id << endl;
        cout << "segment key: " << segments[segment_id].key << endl;
        cout << "segment paras: " << segments[segment_id].paras[0] << " " << segments[segment_id].paras[1] << endl;
        cout << "pred: " << pred << endl;
    }

    /**
     * the linear programming solver to find the lowest absoluter error of an segment
     * @param key_v
     * @param value_v
     * @param lower_index : segment start position in key_v, inclusive
     * @param upper_index : segment end position in key_v, exclusive
     * @param paras
     * @param loss
     */
    void SolveMaxlossLP(const vector<Tk> &key_v, const vector<Tv> &value_v, int lower_index, int upper_index, vector<double> &paras, double &loss) {
        IloEnv env;
        IloModel model(env);
        IloCplex cplex(model);
        IloObjective obj(env);
        IloNumVarArray vars(env);
        IloRangeArray ranges(env);

        env.setOut(env.getNullStream());
        cplex.setOut(env.getNullStream());

        //cplex.setParam(IloCplex::NumericalEmphasis, CPX_ON);
        //cplex.setParam(IloCplex::Param::Advance, 0); // turnoff advanced start
        //cplex.setParam(IloCplex::Param::Preprocessing::Presolve, false); // turnoff presolve

        //cplex.setParam(IloCplex::RootAlg, IloCplex::Primal); // using simplex
        //cplex.setParam(IloCplex::RootAlg, IloCplex::Dual); // using dual simplex
        //cplex.setParam(IloCplex::RootAlg, IloCplex::Barrier); // set optimizer used interior point method

        //cplex.setParam(IloCplex::Param::Parallel, CPX_PARALLEL_DETERMINISTIC); // multi thread
        //cplex.setParam(IloCplex::Param::MIP::Strategy::KappaStats, 1); // https://perso.ensta-paris.fr/~diam/ro/online/cplex/cplex1271/CPLEX/Parameters/topics/MIPKappaStats.html

        // set variable type, IloNumVarArray starts from 0.
        for (int i = 0; i <= degree; i++) {
            vars.add(IloNumVar(env, -INFINITY, INFINITY,ILOFLOAT));  // the weight, i.e., a2 for x^2, start from the lowest, i.e., a0, a1,...,an
        }
        IloNumVar target(env, 0.0, INFINITY, ILOFLOAT); // the max loss value

        // declare objective
        obj.setExpr(target);
        obj.setSense(IloObjective::Minimize);
        model.add(obj);

        vector<Tk> key_term;
        Tk key; // x
        Tk base_key = key_v[lower_index]; // x0 for x - x0
        Tk key_diff; // x - x0
        Tk current_key_term; // (x - x0)^n

        // add constraint for each record
        for (int i = lower_index; i < upper_index; i++) {
            IloExpr model_term(env);
            key = key_v[i];
            key_diff = key - base_key;
            current_key_term = 1.0;
            key_term.clear();

            // 0, (x-x0), (x-x0)^2, ... , (x-x0)^n
            for (int j = 0; j <= degree; j++) {
                key_term.push_back(current_key_term);
                current_key_term *= key_diff;
            }

            // a0*0 + a1*(x-x0) + a2*(x-x0)^2 + ... + an*(x-x0)^n
            for (int j = 0; j <= degree; j++) {
                model_term += vars[j] * key_term[j];
            }

            model.add(model_term - value_v[i] <= target);
            model.add(model_term - value_v[i] >= -target);
        }

        IloNum starttime_ = cplex.getTime();
        cplex.solve();
        IloNum endtime_ = cplex.getTime();

        //cplex.exportModel("path../model.sav");
        //cplex.exportModel("path../model.lp");

        //IloNum violation = cplex.getQuality(IloCplex::MaxPrimalInfeas);
        //cout << "violation: " << violation << endl;
        //IloNum Kappa = cplex.getQuality(IloCplex::Kappa);
        //cout << "kappa: " << Kappa << endl;

        loss = cplex.getObjValue();

        paras.clear();
        try{
            for (int i = 0; i <= degree; i++) {
                paras.push_back(cplex.getValue(vars[i]));
            }
        } catch (IloAlgorithm::NotExtractedException e){
            cout << "catch NotExtractedException from index" << lower_index  << " to " << upper_index << endl;
            paras.clear();
            paras.push_back(value_v[lower_index]);
            for (int i = 1; i <= degree; i++)
                paras.push_back(0);
        }

        env.end();
    }

    // each leaf layer model, i.e., the polynomial function
    struct Segment{
        Tk key; // the first key of this segment
        std::array<double, degree+1> paras; // a0, a1, ... , an

        Segment(Tk key, vector<double> &pas): key(key){
            for(int i = 0; i < degree+1; i++){
                paras[i] = pas[i];
            }
        }

        inline Tv operator()(const Tk &k) const {
            auto key_diff = k - key;
            Tk current_key_term = 1.0;
            Tv pred = paras[0]; // intercept
            for(int i = 1; i <= degree; i++){
                current_key_term *= key_diff;
                pred += paras[i] * current_key_term;
            }
            return pred;
        }

        friend inline bool operator<(const Segment &s, const Tk &k) { return s.key < k; }
        friend inline bool operator<(const Tk &k, const Segment &s) { return k < s.key; }

        inline Tv max_from_left(const Tk &k) const {
            // find the 0 derivative point for deg 1 and 2
            if(degree == 1){
                if(paras[1] > 0)
                    return (*this)(k);
                else
                    return paras[0];
            } else if (degree == 2){
                Tk center = -(2 * paras[2]) / paras[1];
                if (paras[2] > 0){
                    if(center < key)
                        return (*this)(k);
                    else if (center < k)
                        return max((Tv)paras[0], (*this)(k));
                    else
                        return paras[0];
                } else {
                    if(center < key)
                        return paras[0];
                    else if (center < k)
                        return (*this)(center);
                    else
                        return (*this)(k);
                }
            } else {
                cout << "currently no supported degree for max" << endl;
            }
        }

        inline Tv max_from_right(const Tk &k, const Tk &right_border) const {
            if(degree == 1){
                if(paras[1] > 0)
                    return (*this)(right_border);
                else
                    return (*this)(k);
            } else if (degree == 2){
                Tk center = -(2 * paras[2]) / paras[1];
                if (paras[2] > 0){
                    if(center < k)
                        return (*this)(right_border);
                    else if (center < right_border)
                        return max((*this)(k), (*this)(right_border));
                    else
                        return (*this)(k);
                } else {
                    if(center < k)
                        return (*this)(k);
                    else if (center < right_border)
                        return (*this)(center);
                    else
                        return (*this)(right_border);
                }
            } else {
                cout << "currently no supported degree for max" << endl;
            }
        }

        // note: kl and ku should within this segment's range
        inline Tv region_max(const Tk &kl, const Tk &ku) const{
            if(degree == 1){
                if(paras[1] > 0)
                    return (*this)(ku);
                else
                    return (*this)(kl);
            } else if (degree == 2){
                Tk center = -(2 * paras[2]) / paras[1];
                if (paras[2] > 0){
                    if(center < kl)
                        return (*this)(ku);
                    else if (center < ku)
                        return max((*this)(kl), (*this)(ku));
                    else
                        return (*this)(kl);
                } else {
                    if(center < kl)
                        return (*this)(kl);
                    else if (center < ku)
                        return (*this)(center);
                    else
                        return (*this)(ku);
                }
            } else if (degree == 3){
                Tv border_val1 = (*this)(kl), border_val2 = (*this)(ku);
                Tv max_to_return = max(border_val1, border_val2);
                double a = 3 * paras[3], b = 2 * paras[2], c = paras[1]; // term when derivative = 0
                Tk deriv_zero_pos_1 = (-b - sqrt(b*b - 4*a*c))/(2*a); // -b - sqrt(b^2-4ac)/2a
                Tk deriv_zero_pos_2 = (-b + sqrt(b*b - 4*a*c))/(2*a); // -b + sqrt(b^2-4ac)/2a
                if (deriv_zero_pos_1 >= kl && deriv_zero_pos_1 <= ku){
                    max_to_return = max(max_to_return, (*this)(deriv_zero_pos_1));
                }
                if (deriv_zero_pos_2 >= kl && deriv_zero_pos_2 <= ku){
                    max_to_return = max(max_to_return, (*this)(deriv_zero_pos_2));
                }
                return max_to_return;
            } else {
                cout << "currently no supported degree for max" << endl;
            }
        }
    };

    vector<Segment> segments;
    stx::btree<Tk, int> seg_idx; // used to index the segments

    vector<pair<int, int>> intervals; // each segment's corresponding interval [ , ) exclusive

    typedef RTree<int, Tk, 1, double, 64, 32, Tv, Tk, Tv, degree> SimpleRtreeType;
    SimpleRtreeType artree; // based on simple Rtree, collect aggregate info
    vector<Tv> seg_maxs; // aggregate max for each segments

    Tv t_abs; // absolute error threshold
    double t_rel; // relative error threshold
};
#endif //POLYFIT_CLION_POLYFIT_H
