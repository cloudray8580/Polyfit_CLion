//
// Created by Cloud on 2020/10/3.
//

#ifndef POLYFIT_CLION_ARTREE_H
#define POLYFIT_CLION_ARTREE_H

#include <vector>
#include "SimpleRTree.h"
#include "utils.h"

using namespace std;

template<typename Tk, typename Tv, int dims = 1>
class aRtree{
public:
    aRtree(Tv t_abs = 100, double t_rel = 0.01): t_abs(t_abs), t_rel(t_rel) {
        cout << "aRtree constructor: enter" << endl;
    }
    ~aRtree(){
        // mrtree.Reset(); // the destructor of SimpleRTree already call Reset()
        cout << "aRtree destructor: free mrtree" << endl;
    }

    void build1DCount(const std::vector<Tk> &keys){
        for (int i = 0; i < keys.size(); i++) {
            Tk min[1]{keys[i]};
            Tk max[1]{keys[i]};
            artree.Insert(min, max, i);
        }
        std::cout << "finsih inserting data to aRtree" << std::endl;
        artree.GenerateCountAggregate(artree.m_root);
        std::cout << "finsih generating Count aggreate in aRtree" << std::endl;
    }

    void buildTest(){
        for (int i = 0; i < 1000; i++) {
            Tk min[1]{i};
            Tk max[1]{i};
            artree.Insert(min, max, i);
        }
        std::cout << "basic paras of aRtree" << artree.CountBasicParameters() << std::endl;
    }

    // no difference for Problem 1 or 2, as it always find the exact
    QueryResult Query1DCount(const std::vector<Tk> &queryset_low, const std::vector<Tk> &queryset_up, std::vector<Tv> &results, std::string RealResultPath){
        results.clear();
        double result;

        auto t0 = chrono::steady_clock::now();

        typename SimpleRtreeType::Rect rect;
        for(int i = 0; i < queryset_low.size(); i++){
            int leafcount = 0;
            rect.m_min[0] = queryset_low[i];
            rect.m_max[0] = queryset_up[i];
            result = artree.Aggregate(artree.m_root, &rect, leafcount);
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
        query_result.tree_paras = artree.CountBasicParameters();

        return query_result;
    }

    void build2DCount(const std::vector<Tk> &keys1, const std::vector<Tk> &keys2){
        for (int i = 0; i < keys1.size(); i++) {
            //Rect rect(keys1[i], keys2[i], keys1[i], keys2[i]);
            //mrtree.Insert(rect.min, rect.max, i);
            Tk min[2]{keys1[i], keys2[i]};
            Tk max[2]{keys1[i], keys2[i]};
            artree.Insert(min, max, i);
            if (i % 1000000 == 0)
                cout << i << endl;
        }
        std::cout << "finsih inserting data to aRtree" << std::endl;
        artree.GenerateCountAggregate(artree.m_root);
        std::cout << "finsih generating Count aggreate in aRtree" << std::endl;
    }

    QueryResult Query2DCount(const std::vector<Tk> &d1_low, const std::vector<Tk> &d2_low, const std::vector<Tk> &d1_up, const std::vector<Tk> &d2_up,
                             std::vector<Tv> &results, std::string RealResultPath){
        results.clear();
        double result;

        auto t0 = chrono::steady_clock::now();

        typename SimpleRtreeType::Rect rect;
        for(int i = 0; i < d1_low.size(); i++){
            int leafcount = 0;
            rect.m_min[0] = d1_low[i];
            rect.m_min[1] = d2_low[i];
            rect.m_max[0] = d1_up[i];
            rect.m_max[1] = d2_up[i];
            result = artree.Aggregate(artree.m_root, &rect, leafcount);
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


    typedef RTree<int, Tk, dims, double, 16, 8, Tv> SimpleRtreeType;
    SimpleRtreeType artree; // based on simple Rtree
    Tv t_abs; // absolute error threshold
    double t_rel; // relative error threshold
};


#endif //POLYFIT_CLION_ARTREE_H
