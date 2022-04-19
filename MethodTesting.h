//
// Created by Cloud on 2020/8/17.
//

#ifndef POLYFIT_CLION_METHODTESTING_H
#define POLYFIT_CLION_METHODTESTING_H

#include <string>
#include <armadillo>
#include "utils.h"
#include "S2Sampling.h"
#include "FITingTree.h"
#include "RMI.h"
#include "PGM.h"
#include "RadixSpline.h"
#include "DBest.h"
#include "VerdictDB.h"
#include "ICTAI1997.h"
#include "MultiResolutionTree.h"
#include "aRtree.h"
#include "Polyfit_1D.h"
#include "Polyfit_2D.h"
#include "Polyfit.h"
#include "Polyfit2D.h"


using namespace std;

/**
 * Run methods for count aggregate on 1 dimensional data
 * @tparam Tk : type of keys
 * @tparam Tv : type of values
 * @tparam degree : for Polyfit, the maximum degreee used in polynomial function
 * @param keys : the keys
 * @param values : the corresponding values of keys
 * @param queryset_L : queries, lower border
 * @param queryset_U : queries, upper border
 * @param methods : methods from enum defined in utils.h
 * @param real_result_path : the real result file path
 * @param result_save_path : save the result in file path
 * @param Tabs : absolute error threshold
 * @param Trel : relative error threshold
 * @param DoRefinement : if specified relative error threshold, refine those fail the relative error check
 * @param repeat_times : each methods run this many times to retrieve average result
 */
template<typename Tk, typename Tv, int degree>
void RunMethods_1D_Count(vector<Tk> &keys, vector<Tv> &values, vector<Tk> queryset_L, vector<Tk> queryset_U,
                         set<METHODS> methods, string real_result_path, string result_save_path,
                         double Tabs = 100, double Trel = 0.01, bool DoRefinement = true, int repeat_times = 1){

    vector<Tv> predicted_results;
    QueryResult qr;
    std::ofstream header;
    header.open(result_save_path, std::ios::app);
    header << endl;
    header << "Tabs " << Tabs << "," << "Trel " << Trel << endl;
    header << "Name" << "," << "average_query_time" << "," << "total_query_time" << ","
               << "measured_absolute_error" << "," << "measured_relative_error" << ","
               << "hit_count" << "," << "refinement_count" << "," << "total_query_count" << ","
               << "tree_paras" << "," << "leaf_model_amount" << "," << "total_paras" << "," << "total_bytes" << std::endl;
    header.close();

    if(methods.find(METHODS::RMI) != methods.end()) {
        auto construction_start = chrono::steady_clock::now();
        vector<int> arch;
        arch.push_back(1);
        arch.push_back(10);
        arch.push_back(100);
        arch.push_back(1000);
        arma::rowvec dataset;
        arma::rowvec responses;
        VectorToRowvec(dataset, keys);
        VectorToRowvec(responses, values);
        StageModel stage_model(dataset, responses, arch, Tabs, Trel);
        stage_model.DumpParameters();
        auto construction_stop = chrono::steady_clock::now();
        auto construction_time = chrono::duration_cast<chrono::nanoseconds>(construction_stop - construction_start).count();
        cout << "RMI construction time (ns): " << construction_time << endl;

        qr = {}; // reset struct
        auto avg_time = 0;
        for(int i  = 0 ; i < repeat_times; i++){
            predicted_results.clear();
            qr = stage_model.CountPrediction(queryset_L, queryset_U, predicted_results, keys, DoRefinement, real_result_path);
            avg_time += qr.average_query_time;
        }
        avg_time /= repeat_times;
        qr.average_query_time = avg_time;
        WriteQueryResult2File(qr, result_save_path, "RMI");
    }

    if(methods.find(METHODS::FITINGTREE) != methods.end()){
        auto construction_start = chrono::steady_clock::now();
        ATree<Tk, Tv> atree(Tabs, Trel);
        atree.TrainAtree(keys, values);
        auto construction_stop = chrono::steady_clock::now();
        auto construction_time = chrono::duration_cast<chrono::nanoseconds>(construction_stop - construction_start).count();
        cout << "FITingTree construction time (ns): " << construction_time << endl;

        qr = {}; // reset struct
        auto avg_time = 0;
        for(int i  = 0 ; i < repeat_times; i++){
            predicted_results.clear();
            qr = atree.CountPrediction2(queryset_L, queryset_U, predicted_results, keys, DoRefinement, real_result_path);
            avg_time += qr.average_query_time;
        }
        avg_time /= repeat_times;
        qr.average_query_time = avg_time;
        WriteQueryResult2File(qr, result_save_path, "FITingTree");
    }

    if(methods.find(METHODS::PGM) != methods.end()){
        auto construction_start = chrono::steady_clock::now();
        std::vector<int> int_keys(keys.begin(), keys.end());
        std::vector<int> int_ql(queryset_L.begin(), queryset_L.end());
        std::vector<int> int_qu(queryset_U.begin(), queryset_U.end());
        PGMWrapper<Tk, Tv> PGM(int_keys, Tabs, Trel);
        auto construction_stop = chrono::steady_clock::now();
        auto construction_time = chrono::duration_cast<chrono::nanoseconds>(construction_stop - construction_start).count();
        cout << "PGM construction time (ns): " << construction_time << endl;

        qr = {}; // reset struct
        auto avg_time = 0;
        for(int i  = 0 ; i < repeat_times; i++){
            predicted_results.clear();
            qr = PGM.RangeCountQuery(int_ql, int_qu, predicted_results, int_keys, DoRefinement, real_result_path);
            avg_time += qr.average_query_time;
        }
        avg_time /= repeat_times;
        qr.average_query_time = avg_time;
        WriteQueryResult2File(qr, result_save_path, "PGM");
    }

    if(methods.find(METHODS::RADIXSPLINE) != methods.end()){
        vector<uint64_t> unkeys, unvalues, unqueryset_L, unqueryset_U, unpredicted_results;
        if (keys[0] < 0){
            uint64_t base = abs(keys[0]) + 1; // the keys are sorted
            Convert2Unsign(keys, unkeys, base);
            Convert2Unsign(values, unvalues, base);
            Convert2Unsign(queryset_L, unqueryset_L, base);
            Convert2Unsign(queryset_U, unqueryset_U, base);
        } else {
            Convert2Unsign(keys, unkeys, 0);
            Convert2Unsign(values, unvalues, 0);
            Convert2Unsign(queryset_L, unqueryset_L, 0);
            Convert2Unsign(queryset_U, unqueryset_U, 0);
        }
        RSWrapper rsw(Tabs, Trel);
        rsw.build(unkeys);
        qr = {}; // reset struct
        auto avg_time = 0;
        for(int i  = 0 ; i < repeat_times; i++){
            unpredicted_results.clear();
            qr = rsw.RSCountPrediction(unqueryset_L, unqueryset_U, unpredicted_results, unkeys, DoRefinement, real_result_path);
            avg_time += qr.average_query_time;
        }
        avg_time /= repeat_times;
        qr.average_query_time = avg_time;
        WriteQueryResult2File(qr, result_save_path, "RadixSpline");
    }

    if(methods.find(METHODS::POLYFIT) != methods.end()){
        auto construction_start = chrono::steady_clock::now();
        Polyfit<Tk, Tv, degree> polyfit(keys, values, Tabs, Trel);
        auto construction_stop = chrono::steady_clock::now();
        auto construction_time = chrono::duration_cast<chrono::nanoseconds>(construction_stop - construction_start).count();
        cout << "Polyfit construction time (ns): " << construction_time << endl;

        qr = {}; // reset struct
        auto avg_time = 0;
        for(int i  = 0 ; i < repeat_times; i++){
            predicted_results.clear();
            qr = polyfit.CountPrediction(queryset_L, queryset_U, predicted_results, keys, DoRefinement, real_result_path);
            avg_time += qr.average_query_time;
        }
        avg_time /= repeat_times;
        qr.average_query_time = avg_time;
        WriteQueryResult2File(qr, result_save_path, "Polyfit_"+std::to_string(degree));
    }

    if(methods.find(METHODS::DBEST) != methods.end()){
        auto construction_start = chrono::steady_clock::now();
        //const char* model_uri = "/mnt/c/Users/Cloud/iCloudDrive/Polyfit_Python/model_500_3.bin"; // 500 trees, each 3 level
        //const char* model_uri = "/mnt/c/Users/Cloud/iCloudDrive/Polyfit_Python/model.bin"; // 1000 trees, each 3 level
        //const char* model_uri = "/mnt/c/Users/Cloud/iCloudDrive/Polyfit_Python/model_2000_4.bin"; // 2000 trees, each 4 level
        const char* model_uri = "/mnt/c/Users/Cloud/iCloudDrive/Polyfit_Python/model_5000_5.bin"; // 5000 trees, each 5 level
        float *ql = new float[queryset_L.size()];
        float *qu = new float[queryset_U.size()];
        for(int i = 0; i < queryset_L.size(); i++){
            ql[i] = queryset_L[i];
            qu[i] = queryset_U[i];
        }
        auto construction_stop = chrono::steady_clock::now();
        auto construction_time = chrono::duration_cast<chrono::nanoseconds>(construction_stop - construction_start).count();
        cout << "DBest construction time (ns): " << construction_time << endl;

        qr = {}; // reset struct
        auto avg_time = 0;
        for(int i  = 0 ; i < repeat_times; i++){
            predicted_results.clear();
            qr = DBestCountPrediction(model_uri, ql, qu, (int)queryset_L.size(), real_result_path);
            avg_time += qr.average_query_time;
        }
        avg_time /= repeat_times;
        qr.average_query_time = avg_time;
        WriteQueryResult2File(qr, result_save_path, "DBest");
    }

    if(methods.find(METHODS::VERDICTDB) != methods.end()){
        qr = {}; // reset struct
        auto avg_time = 0;
        for(int i  = 0 ; i < repeat_times; i++){
            predicted_results.clear();
            qr = VerdictDBCountAggregate(queryset_L, queryset_U, predicted_results, keys, real_result_path);
            avg_time += qr.average_query_time;
        }
        avg_time /= repeat_times;
        qr.average_query_time = avg_time;
        WriteQueryResult2File(qr, result_save_path, "VerdictDB");
    }

    if(methods.find(METHODS::S2) != methods.end()){
        qr = {}; // reset struct
        auto avg_time = 0;
        for(int i  = 0 ; i < repeat_times; i++){
            predicted_results.clear();
            qr = TestS2Sampling1D<Tk, Tv>(keys, queryset_L, queryset_U, 0.9, Trel, Tabs, predicted_results, real_result_path);
            avg_time += qr.average_query_time;
        }
        avg_time /= repeat_times;
        qr.average_query_time = avg_time;
        WriteQueryResult2File(qr, result_save_path, "S2");
    }

    if(methods.find(METHODS::ENTROPYHIST) != methods.end()){
        // tested in legacy code
    }

    if(methods.find(METHODS::STX) != methods.end()){
        // tested in legacy code
    }

    if(methods.find(METHODS::PLATO) != methods.end()){
        auto construction_start = chrono::steady_clock::now();
        ICTAI1997<Tk, Tv> time_series(keys, values);
        auto construction_stop = chrono::steady_clock::now();
        auto construction_time = chrono::duration_cast<chrono::nanoseconds>(construction_stop - construction_start).count();
        cout << "PLATO construction time (ns): " << construction_time << endl;

        qr = {}; // reset struct
        auto avg_time = 0;
        for(int i  = 0 ; i < repeat_times; i++){
            predicted_results.clear();
            qr = time_series.CountPrediction(queryset_L, queryset_U, predicted_results, keys, real_result_path);
            avg_time += qr.average_query_time;
        }
        avg_time /= repeat_times;
        qr.average_query_time = avg_time;
        WriteQueryResult2File(qr, result_save_path, "Plato");
    }

    if(methods.find(METHODS::MRTREE) != methods.end()){
        auto construction_start = chrono::steady_clock::now();
        MRTree<Tk, Tv, 1> mrt(Tabs, Trel);
        cout << "before function build " << endl;
        mrt.build1DCount(keys);
        cout << "after function build " << endl;
        auto construction_stop = chrono::steady_clock::now();
        auto construction_time = chrono::duration_cast<chrono::nanoseconds>(construction_stop - construction_start).count();
        cout << "MRTree construction time (ns): " << construction_time << endl;

        qr = {}; // reset struct
        auto avg_time = 0;
        for(int i  = 0 ; i < repeat_times; i++){
            predicted_results.clear();
            qr = mrt.Query1DCount(queryset_L, queryset_U, predicted_results, real_result_path, DoRefinement); // use DoRefinement as a mark for prob1 or 2
            avg_time += qr.average_query_time;
        }
        avg_time /= repeat_times;
        qr.average_query_time = avg_time;
        WriteQueryResult2File(qr, result_save_path, "MRTree");
    }

    if(methods.find(METHODS::ARTREE) != methods.end()){
        auto construction_start = chrono::steady_clock::now();
        aRtree<Tk, Tv, 1> art(Tabs, Trel);
        cout << "before function build " << endl;
        art.build1DCount(keys);
        cout << "after function build " << endl;
        auto construction_stop = chrono::steady_clock::now();
        auto construction_time = chrono::duration_cast<chrono::nanoseconds>(construction_stop - construction_start).count();
        cout << "ARTree construction time (ns): " << construction_time << endl;

        qr = {}; // reset struct
        auto avg_time = 0;
        for(int i  = 0 ; i < repeat_times; i++){
            predicted_results.clear();
            qr = art.Query1DCount(queryset_L, queryset_U, predicted_results, real_result_path);
            avg_time += qr.average_query_time;
        }
        avg_time /= repeat_times;
        qr.average_query_time = avg_time;
        WriteQueryResult2File(qr, result_save_path, "aRtree");
    }
}


template<typename Tk, typename Tv, int degree>
void RunMethods_1D_Max(vector<Tk> &keys, vector<Tv> &values, vector<Tk> queryset_L, vector<Tk> queryset_U,
                         set<METHODS> methods, string real_result_path, string result_save_path,
                         double Tabs = 100, double Trel = 0.01, bool DoRefinement = true, int repeat_times = 1){
    vector<Tv> predicted_results;
    QueryResult qr;
    std::ofstream header;
    header.open(result_save_path, std::ios::app);
    header << endl;
    header << "Tabs " << Tabs << "," << "Trel " << Trel << endl;
    header << "Name" << "," << "average_query_time" << "," << "total_query_time" << ","
           << "measured_absolute_error" << "," << "measured_relative_error" << ","
           << "hit_count" << "," << "refinement_count" << "," << "total_query_count" << ","
           << "tree_paras" << "," << "leaf_model_amount" << "," << "total_paras" << "," << "total_bytes" << std::endl;
    header.close();

    if(methods.find(METHODS::POLYFIT) != methods.end()){
        Polyfit<Tk, Tv, degree> polyfit(keys, values, Tabs, Trel);
        polyfit.build_aggregate_for_max(keys, values);
        qr = {}; // reset struct
        auto avg_time = 0;
        for(int i  = 0 ; i < repeat_times; i++){
            predicted_results.clear();
            qr = polyfit.MaxPrediction(queryset_L, queryset_U, predicted_results, keys, values, DoRefinement, real_result_path);
            avg_time += qr.average_query_time;
        }
        avg_time /= repeat_times;
        qr.average_query_time = avg_time;
        WriteQueryResult2File(qr, result_save_path, "Polyfit_max_"+std::to_string(degree));
    }

    if(methods.find(METHODS::MRTREE) != methods.end()){
        cout << "found function MRTREE " << endl;
        MRTree<Tk, Tv, 1> mrt(Tabs, Trel);
        cout << "before function build " << endl;
        mrt.build1DMax(keys, values);
        cout << "after function build " << endl;
        qr = {}; // reset struct
        auto avg_time = 0;
        for(int i  = 0 ; i < repeat_times; i++){
            predicted_results.clear();
            qr = mrt.Query1DMax(queryset_L, queryset_U, predicted_results, values, real_result_path, DoRefinement);
            avg_time += qr.average_query_time;
        }
        avg_time /= repeat_times;
        qr.average_query_time = avg_time;
        WriteQueryResult2File(qr, result_save_path, "MRTree_max");
    }
}

template<typename Tk, typename Tv, int degree>
void RunMethods_2D_Count(vector<Tk> &keys1, vector<Tk> &keys2, vector<Tk> &accu_keys1, vector<Tk> &accu_keys2, vector<Tv> &accu_values,
                         vector<Tk> queryset_L1, vector<Tk> queryset_L2, vector<Tk> queryset_U1, vector<Tk> queryset_U2,
                         set<METHODS> methods, string real_result_path, string result_save_path, tuple<Tk, Tk, Tk, Tk> domain,
                         double Tabs = 100, double Trel = 0.01, bool DoRefinement = true, int repeat_times = 1){

    vector<Tv> predicted_results;
    QueryResult qr;
    std::ofstream header;
    header.open(result_save_path, std::ios::app);
    header << endl;
    header << "Tabs " << Tabs << "," << "Trel " << Trel << endl;
    header << "Name" << "," << "average_query_time" << "," << "total_query_time" << ","
           << "measured_absolute_error" << "," << "measured_relative_error" << ","
           << "hit_count" << "," << "refinement_count" << "," << "total_query_count" << ","
           << "tree_paras" << "," << "leaf_model_amount" << "," << "total_paras" << "," << "total_bytes" << std::endl;
    header.close();

    if(methods.find(METHODS::POLYFIT) != methods.end()){
        Polyfit2D<Tk, Tv, degree>  polyfit2d(Tabs, Trel, domain);
        polyfit2d.build(accu_keys1, accu_keys2, accu_values, 3);
        polyfit2d.build_non_leaf();
        qr = {}; // reset struct
        auto avg_time = 0;
        for(int i  = 0 ; i < repeat_times; i++){
            predicted_results.clear();
            qr = polyfit2d.CountPrediction2D(queryset_L1, queryset_L2, queryset_U1, queryset_U2, predicted_results, DoRefinement, real_result_path);
            avg_time += qr.average_query_time;
        }
        avg_time /= repeat_times;
        qr.average_query_time = avg_time;
        WriteQueryResult2File(qr, result_save_path, "Polyfit_2D_"+std::to_string(degree));
    }

    if(methods.find(METHODS::DBEST) != methods.end()){
        //const char* model_uri = "/mnt/c/Users/Cloud/iCloudDrive/Polyfit_Python/model_2d_5K_3.bin"; // 5K trees, depth 3
        const char* model_uri = "/mnt/c/Users/Cloud/iCloudDrive/Polyfit_Python/model_2d_2W_3.bin"; // 2W trees, depth 3
        //const char* model_uri = "/mnt/c/Users/Cloud/iCloudDrive/Polyfit_Python/model_2d.bin"; // 1W trees, depth 4
        //const char* model_uri = "/mnt/c/Users/Cloud/iCloudDrive/Polyfit_Python/model_2d_5W_4.bin"; // 5W trees, depth 4
        qr = {}; // reset struct
        auto avg_time = 0;
        for(int i  = 0 ; i < repeat_times; i++){
            predicted_results.clear();
            qr = DBestCountPrediction2D(model_uri, queryset_L1, queryset_L2, queryset_U1, queryset_U2, (int)queryset_L1.size(), real_result_path);
            avg_time += qr.average_query_time;
        }
        avg_time /= repeat_times;
        qr.average_query_time = avg_time;
        WriteQueryResult2File(qr, result_save_path, "DBest_2D");
    }

    if(methods.find(METHODS::VERDICTDB) != methods.end()){ // should use the original dataset!!!
        qr = {}; // reset struct
        auto avg_time = 0;
        for(int i  = 0 ; i < repeat_times; i++){
            predicted_results.clear();
            qr = VerdictDBCountAggregate2D(keys1, keys2, queryset_L1, queryset_L2, queryset_U1, queryset_U2, predicted_results, real_result_path);
            avg_time += qr.average_query_time;
        }
        avg_time /= repeat_times;
        qr.average_query_time = avg_time;
        WriteQueryResult2File(qr, result_save_path, "VerdictDB_2D");
    }

    if(methods.find(METHODS::S2) != methods.end()){
        qr = {}; // reset struct
        double p = 0.9;
        auto avg_time = 0;
        for(int i  = 0 ; i < repeat_times; i++){
            predicted_results.clear();
            qr = TestS2Sampling2D<Tk, Tv>(keys1, keys2, queryset_L1, queryset_L2, queryset_U1, queryset_U2, 0.9, Trel, Tabs, predicted_results, real_result_path);
            avg_time += qr.average_query_time;
        }
        avg_time /= repeat_times;
        qr.average_query_time = avg_time;
        WriteQueryResult2File(qr, result_save_path, "S2_2D");
    }

    if(methods.find(METHODS::MRTREE) != methods.end()){
        auto construction_start = chrono::steady_clock::now();
        MRTree<Tk, Tv, 2> mrt(Tabs, Trel);
        mrt.build2DCount(keys1, keys2);
        auto construction_stop = chrono::steady_clock::now();
        auto construction_time = chrono::duration_cast<chrono::nanoseconds>(construction_stop - construction_start).count();
        cout << "MRTree 2d construction time (ns): " << construction_time << endl;

        qr = {}; // reset struct
        auto avg_time = 0;
        for(int i  = 0 ; i < repeat_times; i++){
            predicted_results.clear();
            qr = mrt.Query2DCount(queryset_L1, queryset_L2, queryset_U1, queryset_U2, predicted_results, real_result_path, DoRefinement);
            avg_time += qr.average_query_time;
        }
        avg_time /= repeat_times;
        qr.average_query_time = avg_time;
        WriteQueryResult2File(qr, result_save_path, "MRTree_2D");
    }

    if(methods.find(METHODS::ARTREE) != methods.end()){
        auto construction_start = chrono::steady_clock::now();
        aRtree<Tk, Tv, 2> art(Tabs, Trel);
        art.build2DCount(keys1, keys2);
        auto construction_stop = chrono::steady_clock::now();
        auto construction_time = chrono::duration_cast<chrono::nanoseconds>(construction_stop - construction_start).count();
        cout << "ARTree 2d construction time (ns): " << construction_time << endl;

        qr = {}; // reset struct
        auto avg_time = 0;
        for(int i  = 0 ; i < repeat_times; i++){
            predicted_results.clear();
            qr = art.Query2DCount(queryset_L1, queryset_L2, queryset_U1, queryset_U2, predicted_results, real_result_path);
            avg_time += qr.average_query_time;
        }
        avg_time /= repeat_times;
        qr.average_query_time = avg_time;
        WriteQueryResult2File(qr, result_save_path, "aRtree_2D");
    }
}


// = = = = = = For methods that can not be configured with Tabs and Trel = = = = = =

/**
 * Since DBest cannot provide any error guarantee, it's tested by measured error and response time
 * @param model_paths
 * @param real_result_path
 * @param result_save_path
 * @return
 */
template<typename Tk, typename Tv>
QueryResult RunDBest(const vector<Tk> &queryset_L, const vector<Tk> &queryset_U, const vector<string> &model_paths,
                     string real_result_path, string result_save_path, int repeat_times = 1){

    QueryResult qr;
    vector<Tv> predicted_results(queryset_L.size());
    for (string mp : model_paths){
        const char* model_uri = mp.c_str();

        float *ql = new float[queryset_L.size()];
        float *qu = new float[queryset_U.size()];
        for(int i = 0; i < queryset_L.size(); i++){
            ql[i] = queryset_L[i];
            qu[i] = queryset_U[i];
        }
        qr = {}; // reset struct
        auto avg_time = 0;
        for(int i  = 0 ; i < repeat_times; i++){
            predicted_results.clear();
            qr = DBestCountPrediction(model_uri, ql, qu, (int)queryset_L.size(), real_result_path);
            avg_time += qr.average_query_time;
        }
        avg_time /= repeat_times;
        qr.average_query_time = avg_time;
        WriteQueryResult2File(qr, result_save_path, "DBest");
    }
}
template<typename Tk, typename Tv>
QueryResult RunVerdictDB(const vector<Tk> &queryset_L, const vector<Tk> &queryset_U, const vector<Tk> &keys, int ns, int b, double confidence,
                         string real_result_path, string result_save_path, int repeat_times = 1){
    QueryResult qr = {}; // reset struct
    vector<Tv> predicted_results;
    auto avg_time = 0;
    for(int i  = 0 ; i < repeat_times; i++){
        predicted_results.clear();
        qr = VerdictDBCountAggregate<Tk, Tv>(queryset_L, queryset_U, predicted_results, keys, real_result_path, ns, b, confidence);
        avg_time += qr.average_query_time;
    }
    avg_time /= repeat_times;
    qr.average_query_time = avg_time;
    //cout << qr.average_query_time << "," << qr.total_query_time << "," << qr.measured_absolute_error << "," << qr.measured_relative_error << endl;
    WriteQueryResult2File(qr, result_save_path, "VerdictDB");
}

// = = = = = = Legacy Contents = = = = = =


template<typename Tk, typename Tv, int degree>
QueryResult TestPolyfit(vector<Tk> &keys, vector<Tv> &values, vector<Tk> queryset_L, vector<Tk> queryset_U, string RealResultPath,
                        double Tabs = 100, double Trel = 0.01, bool DoRefinement = true, int repeat_times = 1) {

    vector<Tv> predicted_results;
    Polyfit<Tk, Tv, degree> polyfit(keys, values, Tabs, Trel);

    QueryResult qr;
    auto avg_time = 0;
    for(int i  = 0 ; i < repeat_times; i++){
        qr = polyfit.CountPrediction(queryset_L, queryset_U, predicted_results, keys, DoRefinement, RealResultPath);
        avg_time += qr.average_query_time;
    }
    avg_time /= repeat_times;
    qr.average_query_time = avg_time;
    return qr;
}


QueryResult TestPolyfit_MAX(vector<double> &keys, vector<double> &values, vector<double> queryset_L, vector<double> queryset_U, double Trel = 0.01, double Tabs = 100, int highest_term = 1, bool DoRefinement = true) {

    double result;
    vector<double> predicted_results;
    ReverseMaxlossOptimal RMLO(Tabs, Trel, highest_term);
    RMLO.SegmentOnTrainMaxLossModel(keys, values);
    RMLO.BuildNonLeafLayerWithBtree();
    RMLO.PrepareMaxAggregateTree(keys, values);
    RMLO.PrepareExactAggregateMaxTree(keys, values);
    QueryResult query_result = RMLO.MaxPrediction(queryset_L, queryset_U, predicted_results, DoRefinement);
    //RMLO.ExportDatasetRangeAndModels();
    return query_result;
}

QueryResult TestS2Sampling2D(vector<double> &keys1, vector<double> &keys2, vector<double> &queryset_L1, vector<double> &queryset_L2, vector<double> &queryset_U1, vector<double> &queryset_U2, double p = 0.9, double Trel = 0.01, double Tabs = 100) {

    double result;
    vector<int> predicted_results;

    auto t0 = chrono::steady_clock::now();

    for (int i = 0; i < queryset_L1.size(); i++) {
        result = SequentialSampling2D(keys1, keys2, queryset_L1[i], queryset_L2[i], queryset_U1[i], queryset_U2[i], p, Trel, Tabs);
        predicted_results.push_back(int(result));
    }

    auto t1 = chrono::steady_clock::now();

    auto average_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count() / queryset_L1.size();
    auto total_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();

    double MEabs, MErel;

    // check correctness
    MeasureAccuracy(predicted_results, "C:/Users/Cloud/iCloudDrive/LearnedAggregate/VLDB_Final_Experiments/RealQueryResults/OSM_2D.csv", MEabs, MErel);

    //ErrorGuaranteedSampling();

    QueryResult query_result;
    query_result.average_query_time = average_time;
    query_result.total_query_time = total_time;
    query_result.measured_absolute_error = MEabs;
    query_result.measured_relative_error = MErel;

    return query_result;
}

// call when the Abs error is different
QueryResult TestPolyfit_COUNT2D(vector<double> &key1, vector<double> &key2, vector<double> &queryset_L1, vector<double> &queryset_L2, vector<double> &queryset_U1, vector<double> &queryset_U2, double Trel = 0.01, double Tabs = 100, bool DoRefinement = true) {

    double result;
    vector<double> predicted_results;

    Maxloss2D_QuadDivide model2d(Tabs, Trel, -180.0, 180.0, -90.0, 90.0);
    model2d.GenerateKeysAndAccuFromFile("C:/Users/Cloud/iCloudDrive/LearnedAggregate/Sampled2D_100M_1000_1000_ADJUSTED.csv");
    model2d.TrainModel();

    int AbsErr = int(Tabs);
    string AbsErrStr = to_string(AbsErr);
    string filename = "C:/Users/Cloud/Desktop/LearnedAggregateData/2D_LP_models_100M_1000_1000_" + AbsErrStr + ".csv";

    // try to save models to file
    model2d.WriteTrainedModelsToFile(filename);

    // try to read models from file
    //model2d.ReadTrainedModelsFromFile(filename);
    //model2d.LoadRtree();

    model2d.PrepareExactAggregateRtree(key1, key2);
    QueryResult query_result = model2d.CountPrediction2(queryset_L1, queryset_L2, queryset_U1, queryset_U2, predicted_results, DoRefinement);

    return query_result;
}

#endif //POLYFIT_CLION_METHODTESTING_H
