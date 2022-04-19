//
// Created by Cloud on 2020/9/3.
//

#ifndef POLYFIT_CLION_MULTIRESOLUTIONTREE_H
#define POLYFIT_CLION_MULTIRESOLUTIONTREE_H

#include <vector>
#include "SimpleRTree.h"
#include "utils.h"

template<typename Tk, typename Tv, int dims = 2>
class MRTree{
public:
    MRTree(Tv t_abs = 100, double t_rel = 0.01): t_abs(t_abs), t_rel(t_rel) {
        cout << "MRTree constructor: enter" << endl;
    }
    ~MRTree(){
        // mrtree.Reset(); // the destructor of SimpleRTree already call Reset()
        cout << "MRTree destructor: free mrtree" << endl;
    }

    void build1DCount(const std::vector<Tk> &keys){
        for (int i = 0; i < keys.size(); i++) {
            Tk min[1]{keys[i]};
            Tk max[1]{keys[i]};
            mrtree.Insert(min, max, i);
        }
        std::cout << "finsih inserting data to MRTree" << std::endl;
        mrtree.GenerateCountAggregate(mrtree.m_root);
        std::cout << "finsih generating Count aggreate in MRTree" << std::endl;
        mrtree.PropagateRect(mrtree.m_root);
        std::cout << "finsih generating Node Rect in MRTree" << std::endl;
    }

    QueryResult Query1DCount(const std::vector<Tk> &queryset_low, const std::vector<Tk> &queryset_up, std::vector<Tv> &results, std::string RealResultPath, bool isproblem2 = false){
        results.clear();
        double result;

        auto t0 = chrono::steady_clock::now();

        typename SimpleRtreeType::Rect rect;
        for(int i = 0; i < queryset_low.size(); i++){
            rect.m_min[0] = queryset_low[i];
            rect.m_max[0] = queryset_up[i];
            result = mrtree.ProgressiveCount(mrtree.m_root, &rect, t_abs, t_rel, isproblem2);
            results.push_back(result);
        }

        auto t1 = chrono::steady_clock::now();

        // measure query performance
        auto average_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count() / queryset_low.size();
        auto total_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
        double MEabs, MErel;
        MeasureAccuracy(results, RealResultPath, MEabs, MErel);

        QueryResult query_result;
        query_result.average_query_time = average_time;
        query_result.total_query_time = total_time;
        query_result.measured_absolute_error = MEabs;
        query_result.measured_relative_error = MErel;
        query_result.tree_paras = mrtree.CountBasicParameters();

        return query_result;
    }

    // dims should be 2
    void build2DCount(const std::vector<Tk> &keys1, const std::vector<Tk> &keys2){
        for (int i = 0; i < keys1.size(); i++) {
            //Rect rect(keys1[i], keys2[i], keys1[i], keys2[i]);
            //mrtree.Insert(rect.min, rect.max, i);
            Tk min[2]{keys1[i], keys2[i]};
            Tk max[2]{keys1[i], keys2[i]};
            mrtree.Insert(min, max, i);
            if (i % 1000000 == 0)
                cout << i << endl;
        }
        std::cout << "finsih inserting data to MRTree" << std::endl;
        mrtree.GenerateCountAggregate(mrtree.m_root);
        std::cout << "finsih generating Count aggreate in MRTree" << std::endl;
        mrtree.PropagateRect(mrtree.m_root);
        std::cout << "finsih generating Node Rect in MRTree" << std::endl;
    }

    QueryResult Query2DCount(const std::vector<Tk> &d1_low, const std::vector<Tk> &d2_low, const std::vector<Tk> &d1_up, const std::vector<Tk> &d2_up,
                             std::vector<Tv> &results, std::string RealResultPath, bool isproblem2 = false){
        results.clear();
        double result;

        auto t0 = chrono::steady_clock::now();

        typename SimpleRtreeType::Rect rect;
        for(int i = 0; i < d1_low.size(); i++){
            rect.m_min[0] = d1_low[i];
            rect.m_min[1] = d2_low[i];
            rect.m_max[0] = d1_up[i];
            rect.m_max[1] = d2_up[i];
            result = mrtree.ProgressiveCount(mrtree.m_root, &rect, t_abs, t_rel, isproblem2);
            results.push_back(result);
        }

        auto t1 = chrono::steady_clock::now();

        // measure query performance
        auto average_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count() / d1_low.size();
        auto total_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
        double MEabs, MErel;
        MeasureAccuracy(results, RealResultPath, MEabs, MErel);

        QueryResult query_result;
        query_result.average_query_time = average_time;
        query_result.total_query_time = total_time;
        query_result.measured_absolute_error = MEabs;
        query_result.measured_relative_error = MErel;

        return query_result;
    }

    // dims should be 1
    void build1DMax(const std::vector<Tk> &keys, const std::vector<Tv> &values){
        mrtree.Reset();
        for (int i = 0; i < keys.size(); i++) {
            Tk min[1]{keys[i]};
            Tk max[1]{keys[i]};
            mrtree.Insert(min, max, i);
        }
        std::cout << "finsih inserting data to MRTree" << std::endl;
        mrtree.GenerateMaxAggregate(mrtree.m_root, values);
        std::cout << "finsih generating Max aggreate in MRTree" << std::endl;
    }

    QueryResult Query1DMax(const std::vector<Tk> &queryset_low, const std::vector<Tk> &queryset_up, std::vector<Tv> &results,
                           const std::vector<Tv> &values, std::string RealResultPath, bool isproblem2 = false){

        results.clear();
        double result;

        auto t0 = chrono::steady_clock::now();

        typename SimpleRtreeType::Rect rect;
        for(int i = 0; i < queryset_low.size(); i++){
            rect.m_min[0] = queryset_low[i];
            rect.m_max[0] = queryset_up[i];
            result = mrtree.ProgressiveMax(mrtree.m_root, &rect, values, t_abs, t_rel, isproblem2);
            results.push_back(result);
        }

        auto t1 = chrono::steady_clock::now();

        // measure query performance
        auto average_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count() / queryset_low.size();
        auto total_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
        double MEabs, MErel;
        MeasureAccuracy(results, RealResultPath, MEabs, MErel);

        QueryResult query_result;
        query_result.average_query_time = average_time;
        query_result.total_query_time = total_time;
        query_result.measured_absolute_error = MEabs;
        query_result.measured_relative_error = MErel;

        return query_result;
    }

    typedef RTree<int, Tk, dims, double, 16, 8, Tv> SimpleRtreeType;
    //typedef RTree<int, Tk, dims, double, 256, 128, Tv> SimpleRtreeType; // used for max
    SimpleRtreeType mrtree; // based on simple Rtree
    Tv t_abs; // absolute error threshold
    double t_rel; // relative error threshold
};


#endif //POLYFIT_CLION_MULTIRESOLUTIONTREE_H
