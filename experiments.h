//
// Created by Cloud on 2020/8/18.
//

#ifndef POLYFIT_CLION_EXPERIMENTS_H
#define POLYFIT_CLION_EXPERIMENTS_H

#include "MethodTesting.h"
#include "utils.h"

// = = = = = Absolute Error = = = = =

// without refinement
void experiment_abserr_time_1D_count(){
    vector<int> keys, queryset_L, queryset_U;
    vector<int> values;
    LoadTweetDataset(keys, values);
    LoadTweetQuerySet(queryset_L, queryset_U);
    //LoadOSM_Dataset_Queryset_1D_SUBSET(keys, values, queryset_L, queryset_U, 3, false);

    set<METHODS> methods{ARTREE, MRTREE, RMI};
    //set<METHODS> methods{RMI};
    //set<METHODS> methods{PLATO};
    string real_result_path = "/mnt/c/Users/Cloud/iCloudDrive/LearnedAggregate/VLDB_Final_Experiments/RealQueryResults/TWEET_1D.csv";
    string result_save_path = "/mnt/d/Polyfit_RunResults/EDBT2021/abserr_time_1D_count.csv";

    vector<double> Eabs_collection = { 25, 50, 100, 200, 500, 1000 };
    //vector<double> Eabs_collection = { 100 };
    double Erel = 0.01;
    //bool DoRefinement = false;
    bool DoRefinement = true; // only used for measure dataset size. i.e., experiment 4
    int repeat_times = 10;

    for(auto Eabs : Eabs_collection) {
        RunMethods_1D_Count<int, int, 2>(keys, values, queryset_L, queryset_U, methods, real_result_path, result_save_path, Eabs, Erel, DoRefinement, repeat_times);
    }
}

// without refinement
void experiment_abserr_time_2D_count(){
    vector<double> keys1, keys2;
    LoadOSMDataset(keys1, keys2);
    //LoadOSMDataset2DSubset(keys1, keys2, 3);
    vector<double> accu_keys1, accu_keys2, accu_values;
    LoadOSMAccumulationSurface(accu_keys1, accu_keys2, accu_values);
    vector<double> ql1, ql2, qu1, qu2, results;
    LoadOSMQuerySet(ql1, ql2, qu1, qu2);

    set<METHODS> methods{ARTREE, MRTREE};
    //set<METHODS> methods{POLYFIT};
    string real_result_path = "/mnt/c/Users/Cloud/iCloudDrive/LearnedAggregate/VLDB_Final_Experiments/RealQueryResults/OSM_2D.csv";
    //string real_result_path = "/mnt/d/Polyfit_Dataset/OSM/100M_Long_Lat_queryset_withnoise_result.csv";
    string result_save_path = "/mnt/d/Polyfit_RunResults/EDBT2021/abserr_time_2D_count.csv";

    tuple<double, double, double, double> domain(-180.0, -90.0, 180.0, 90.0);
    vector<double> Eabs_collection = { 50 };
    //vector<double> Eabs_collection = { 25, 50, 125, 250 };
    double Erel = 0.01;
    bool DoRefinement = false;
    int repeat_times = 10;

    for(auto Eabs : Eabs_collection) {
        RunMethods_2D_Count<double, double, 4>(keys1, keys2, accu_keys1, accu_keys2, accu_values, ql1, ql2, qu1, qu2, methods, real_result_path, result_save_path, domain, Eabs, Erel, DoRefinement, repeat_times);
    }
}

void experiment_abserr_time_1D_max(){
    vector<long> keys, values, queryset_L, queryset_U;
    LoadHKIDataset(keys, values);
    LoadHKIQuerySet(queryset_L, queryset_U);

    //set<METHODS> methods{PGM, RADIXSPLINE, POLYFIT, VERDICTDB, DBEST};
    set<METHODS> methods{POLYFIT};
    string real_result_path = "/mnt/c/Users/Cloud/iCloudDrive/LearnedAggregate/VLDB_Final_Experiments/RealQueryResults/HKI_MAX.csv";
    string result_save_path = "/mnt/d/Polyfit_RunResults/EDBT2021/abserr_time_1D_max.csv";

    //vector<double> Eabs_collection = { 50, 100, 200, 500, 1000 };
    vector<double> Eabs_collection = { 100 };
    double Erel = 0.01;
    bool DoRefinement = false;
    int repeat_times = 10;

    for(auto Eabs : Eabs_collection) {
        RunMethods_1D_Max<long, long, 1>(keys, values, queryset_L, queryset_U, methods, real_result_path, result_save_path, Eabs, Erel, DoRefinement, repeat_times);
    }
}

// = = = = = Relative Error = = = = =

void experiment_relerr_time_1D_count(){
    vector<int> keys, queryset_L, queryset_U;
    vector<int> values;
    LoadTweetDataset(keys, values);
    LoadTweetQuerySet(queryset_L, queryset_U);

    set<METHODS> methods{PGM, FITINGTREE, POLYFIT};
    //set<METHODS> methods{MRTREE};
    string real_result_path = "/mnt/c/Users/Cloud/iCloudDrive/LearnedAggregate/VLDB_Final_Experiments/RealQueryResults/TWEET_1D.csv";
    string result_save_path = "/mnt/d/Polyfit_RunResults/EDBT2021/relerr_time_1D_count.csv";

    double Eabs = 50;
    vector<double> Erel_collection = {0.005, 0.01, 0.05, 0.10, 0.20};
    bool DoRefinement = true;
    int repeat_times = 10;

    for(auto Erel : Erel_collection) {
        RunMethods_1D_Count<int, int, 2>(keys, values, queryset_L, queryset_U, methods, real_result_path, result_save_path, Eabs, Erel, DoRefinement, repeat_times);
    }
}

void experiment_relerr_time_2D_count(){
    vector<double> keys1, keys2;
    LoadOSMDataset(keys1, keys2);
    vector<double> accu_keys1, accu_keys2, accu_values;
    LoadOSMAccumulationSurface(accu_keys1, accu_keys2, accu_values);
    vector<double> ql1, ql2, qu1, qu2, results;
    LoadOSMQuerySet(ql1, ql2, qu1, qu2);

    set<METHODS> methods{MRTREE, ARTREE};
    //set<METHODS> methods{MRTREE};
    string real_result_path = "/mnt/c/Users/Cloud/iCloudDrive/LearnedAggregate/VLDB_Final_Experiments/RealQueryResults/OSM_2D.csv";
    string result_save_path = "/mnt/d/Polyfit_RunResults/EDBT2021/relerr_time_2D_count.csv";

    tuple<double, double, double, double> domain(-180.0, -90.0, 180.0, 90.0);
    double Eabs = 100;
    vector<double> Erel_collection = {0.005, 0.01, 0.05, 0.10, 0.20};
    bool DoRefinement = true;
    int repeat_times = 10;

    for(auto Erel : Erel_collection) {
        RunMethods_2D_Count<double, double, 2>(keys1, keys2, accu_keys1, accu_keys2, accu_values, ql1, ql2, qu1, qu2, methods, real_result_path, result_save_path, domain, Eabs, Erel, DoRefinement, repeat_times);
    }

    // the below code is optimized for polyfit only
//    Polyfit2D<double, double, 3> polyfit2d(Eabs, 0.01, domain);
//    polyfit2d.build(accu_keys1, accu_keys2, accu_values, 3);
//    polyfit2d.build_non_leaf();
//
//    QueryResult qr;
//    vector<double> predicted_results;
//    for(auto Erel : Erel_collection) {
//        polyfit2d.t_rel = Erel;
//        qr = {}; // reset struct
//        auto avg_time = 0;
//        for(int i  = 0 ; i < repeat_times; i++){
//            predicted_results.clear();
//            qr = polyfit2d.CountPrediction2D(ql1, ql2, qu1, qu2, predicted_results, DoRefinement, real_result_path);
//            avg_time += qr.average_query_time;
//        }
//        avg_time /= repeat_times;
//        qr.average_query_time = avg_time;
//        WriteQueryResult2File(qr, result_save_path, "Polyfit_2d_"+std::to_string(2)+"_erel_"+std::to_string(Erel));
//    }
}

void experiment_relerr_time_1D_max() {
    vector<long> keys, values, queryset_L, queryset_U;
    LoadHKIDataset(keys, values);
    LoadHKIQuerySet(queryset_L, queryset_U);

    //set<METHODS> methods{PGM, RADIXSPLINE, POLYFIT, VERDICTDB, DBEST};
    set<METHODS> methods{MRTREE};
    string real_result_path = "/mnt/c/Users/Cloud/iCloudDrive/LearnedAggregate/VLDB_Final_Experiments/RealQueryResults/HKI_MAX.csv";
    string result_save_path = "/mnt/d/Polyfit_RunResults/EDBT2021/relerr_time_1D_max.csv";

    double Eabs = 100;
    vector<double> Erel_collection = {0.005, 0.01, 0.05, 0.10, 0.20};
    bool DoRefinement = true;
    int repeat_times = 10;

    for(auto Erel : Erel_collection) {
        RunMethods_1D_Max<long, long, 2>(keys, values, queryset_L, queryset_U, methods, real_result_path, result_save_path, Eabs, Erel, DoRefinement, repeat_times);
    }

    // the below code is optimized for polyfit only
//    Polyfit<long, long, 1> polyfit(keys, values, Eabs, 0.01);
//    polyfit.build_aggregate_for_max(keys, values);
//    QueryResult qr;
//    vector<long> predicted_results;
//    for(auto Erel : Erel_collection) {
//        polyfit.t_rel = Erel;
//        qr = {}; // reset struct
//        auto avg_time = 0;
//        for(int i  = 0 ; i < repeat_times; i++){
//            predicted_results.clear();
//            qr = polyfit.MaxPrediction(queryset_L, queryset_U, predicted_results, keys, values, DoRefinement, real_result_path);
//            avg_time += qr.average_query_time;
//        }
//        avg_time /= repeat_times;
//        qr.average_query_time = avg_time;
//        WriteQueryResult2File(qr, result_save_path, "Polyfit_max_"+std::to_string(1)+"_erel_"+std::to_string(Erel));
//    }
}

// = = = = = Degree = = = = =

void experiment_vary_degree() {

    vector<int> keys, values;
    LoadTweetDataset(keys, values);

    string result_save_path = "/mnt/d/Polyfit_RunResults/EDBT2021/vary_degree.csv";

    std::ofstream header;
    header.open(result_save_path, std::ios::app);
    header << endl;

    // measure construction time and query time for deg 1, 2, and 3
    auto t0 = chrono::steady_clock::now();
    Polyfit<int, int, 1> polyfit_deg1(keys, values, 100, 0.01);
    auto t1 = chrono::steady_clock::now();
    auto construction_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
    header << "construction_time" << "," << construction_time << "," << " for degree" << "," << 1 << endl;

    t0 = chrono::steady_clock::now();
    Polyfit<int, int, 2> polyfit_deg2(keys, values, 100, 0.01);
    t1 = chrono::steady_clock::now();
    construction_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
    header << "construction_time" << "," << construction_time << "," << " for degree" << "," << 2 << endl;

    t0 = chrono::steady_clock::now();
    Polyfit<int, int, 3> polyfit_deg3(keys, values, 100, 0.01);
    t1 = chrono::steady_clock::now();
    construction_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
    header << "construction_time" << "," << construction_time << "," << " for degree" << "," << 3 << endl;

    t0 = chrono::steady_clock::now();
    Polyfit<int, int, 4> polyfit_deg4(keys, values, 100, 0.01);
    t1 = chrono::steady_clock::now();
    construction_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
    header << "construction_time" << "," << construction_time << "," << " for degree" << "," << 4 << endl;
}

void experiment_vary_degree_2D() {

    vector<double> keys1, keys2;
    //LoadOSMDataset(keys1, keys2);
    vector<double> accu_keys1, accu_keys2, accu_values;
    //LoadOSMAccumulationSurface(accu_keys1, accu_keys2, accu_values);
    vector<double> ql1, ql2, qu1, qu2, results;
    LoadOSMQuerySet(ql1, ql2, qu1, qu2);

    LoadOSMAccumulationSurface2DSubset(accu_keys1, accu_keys2, accu_values, 4);
    LoadOSMDataset2DSubset(keys1, keys2, 4);

    string real_result_path = "/mnt/c/Users/Cloud/iCloudDrive/LearnedAggregate/VLDB_Final_Experiments/RealQueryResults/OSM_2D.csv";
    string result_save_path = "/mnt/d/Polyfit_RunResults/EDBT2021/vary_degree.csv";

    std::ofstream header;
    header.open(result_save_path, std::ios::app);
    header << endl;

    tuple<double, double, double, double> domain(-180.0, -90.0, 180.0, 90.0);

    // measure construction time and query time for deg 1, 2, and 3
    auto t0 = chrono::steady_clock::now();
    Polyfit2D<double, double, 1> polyfit2d_deg1(1000, 0.01, domain);
    polyfit2d_deg1.build(accu_keys1, accu_keys2, accu_values, 3);
    polyfit2d_deg1.build_non_leaf();
    auto t1 = chrono::steady_clock::now();
    auto construction_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
    header << "polyfit2d construction_time" << "," << construction_time << "," << " for degree" << "," << 1 << "," << "sgements" << "," << polyfit2d_deg1.models.size() << endl;

    t0 = chrono::steady_clock::now();
    Polyfit2D<double, double, 2> polyfit2d_deg2(1000, 0.01, domain);
    polyfit2d_deg2.build(accu_keys1, accu_keys2, accu_values, 4);
    polyfit2d_deg2.build_non_leaf();
    t1 = chrono::steady_clock::now();
    construction_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
    header << "polyfit2d construction_time" << "," << construction_time << "," << " for degree" << "," << 2 << "," << "sgements" << "," << polyfit2d_deg2.models.size() << endl;

    t0 = chrono::steady_clock::now();
    Polyfit2D<double, double, 3> polyfit2d_deg3(1000, 0.01, domain);
    polyfit2d_deg3.build(accu_keys1, accu_keys2, accu_values, 5);
    polyfit2d_deg3.build_non_leaf();
    t1 = chrono::steady_clock::now();
    construction_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
    header << "polyfit2d construction_time" << "," << construction_time << "," << " for degree" << "," << 3 << "," << "sgements" << "," << polyfit2d_deg3.models.size() << endl;

    t0 = chrono::steady_clock::now();
    Polyfit2D<double, double, 4> polyfit2d_deg4(1000, 0.01, domain);
    polyfit2d_deg4.build(accu_keys1, accu_keys2, accu_values, 6);
    polyfit2d_deg4.build_non_leaf();
    t1 = chrono::steady_clock::now();
    construction_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
    header << "polyfit2d construction_time" << "," << construction_time << "," << " for degree" << "," << 4 << "," << "sgements" << "," << polyfit2d_deg4.models.size() << endl;
}

// = = = = = Dataset Size = = = = =

void experiment_vary_datasize(){
    vector<double> keys, queryset_L, queryset_U;
    vector<double> values, predicted_results;

    constexpr int degree = 2;
    string real_result_path = "/mnt/c/Users/Cloud/iCloudDrive/LearnedAggregate/VLDB_Final_Experiments/RealQueryResults/TWEET_1D.csv";
    string result_save_path = "/mnt/d/Polyfit_RunResults/EDBT2021/vary_datasize.csv";

    std::ofstream header;
    header.open(result_save_path, std::ios::app);
    header << endl;
    header << "Name" << "," << "average_query_time" << "," << "total_query_time" << ","
           << "measured_absolute_error" << "," << "measured_relative_error" << ","
           << "hit_count" << "," << "refinement_count" << "," << "total_query_count" << ","
           << "tree_paras" << "," << "leaf_model_amount" << "," << "total_paras" << "," << "total_bytes" << std::endl;
    header.close();

    QueryResult qr = {};
    auto avg_time = 0;
    for(int i = 4; i <= 4; i++){ // should test from 1 to 4
        LoadOSM_Dataset_Queryset_1D_SUBSET(keys, values, queryset_L, queryset_U, i, false); // true for largeInt version

//        auto t0 = chrono::steady_clock::now();
//
//        Polyfit<int, int, degree> polyfit(keys, values, 100, 0.01);
//
//        auto t1 = chrono::steady_clock::now();
//        auto construction_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
//
//        cout << "construction time for polyfit degree "+std::to_string(degree)+"option " << i << " is: " << construction_time << endl;
//
//        header.open(result_save_path, std::ios::app);
//        header << "construction_time" << "," << construction_time << "," << "dataset size opt" << "," << i << endl;
//        header.close();
//
//        qr = {}; // reset struct
//        avg_time = 0;
//        for(int i  = 0 ; i < 5; i++){
//            predicted_results.clear();
//            qr = polyfit.CountPrediction(queryset_L, queryset_U, predicted_results, keys, false, real_result_path, false);
//            avg_time += qr.average_query_time;
//        }
//        avg_time /= 5;
//        qr.average_query_time = avg_time;
//        WriteQueryResult2File(qr, result_save_path, "Polyfit_"+std::to_string(degree)+"_dataset_"+std::to_string(i));

//        PGMWrapper<double, double> PGM(keys, 100, 0.01);
//        qr = {};
//        avg_time = 0;
//        for(int i  = 0 ; i < 5; i++){
//            predicted_results.clear();
//            qr = PGM.RangeCountQuery(queryset_L, queryset_U, predicted_results, keys, false, real_result_path, true);
//            avg_time += qr.average_query_time;
//        }
//        avg_time /= 5;
//        qr.average_query_time = avg_time;
//        WriteQueryResult2File(qr, result_save_path, "PGM_dataset_"+std::to_string(i));


        ATree<double, double> atree(100, 0.01);
        atree.TrainAtree(keys, values);
        qr = {}; // reset struct
        avg_time = 0;
        for(int i  = 0 ; i < 5; i++){
            predicted_results.clear();
            qr = atree.CountPrediction2(queryset_L, queryset_U, predicted_results, keys, false, real_result_path);
            avg_time += qr.average_query_time;
        }
        avg_time /= 5;
        qr.average_query_time = avg_time;
        WriteQueryResult2File(qr, result_save_path, "FITingTree");


//        vector<uint64_t> unkeys, unvalues, unqueryset_L, unqueryset_U, unpredicted_results;
//        if (keys[0] < 0){
//            uint64_t base = abs(keys[0]) + 1; // the keys are sorted
//            Convert2Unsign(keys, unkeys, base);
//            Convert2Unsign(values, unvalues, base);
//            Convert2Unsign(queryset_L, unqueryset_L, base);
//            Convert2Unsign(queryset_U, unqueryset_U, base);
//        } else {
//            Convert2Unsign(keys, unkeys, 0);
//            Convert2Unsign(values, unvalues, 0);
//            Convert2Unsign(queryset_L, unqueryset_L, 0);
//            Convert2Unsign(queryset_U, unqueryset_U, 0);
//        }
//        RSWrapper rsw(100, 0.01);
//        rsw.build(unkeys);
//        qr = {}; // reset struct
//        avg_time = 0;
//        for(int i  = 0 ; i < 5; i++){
//            unpredicted_results.clear();
//            qr = rsw.RSCountPrediction(unqueryset_L, unqueryset_U, unpredicted_results, unkeys, false, real_result_path, false);
//            avg_time += qr.average_query_time;
//        }
//        avg_time /= 5;
//        qr.average_query_time = avg_time;
//        WriteQueryResult2File(qr, result_save_path, "RadixSpline_dataset_"+std::to_string(i));
    }
}

void experiment_vary_datasize_2D() {

    vector<double> keys1, keys2;
    vector<double> accu_keys1, accu_keys2, accu_values;
    vector<double> ql1, ql2, qu1, qu2, results;
    LoadOSMQuerySet(ql1, ql2, qu1, qu2);

    string real_result_path = "/mnt/c/Users/Cloud/iCloudDrive/LearnedAggregate/VLDB_Final_Experiments/RealQueryResults/TWEET_1D.csv"; // in case of optimization
    string result_save_path = "/mnt/d/Polyfit_RunResults/EDBT2021/vary_datasize_2d.csv";

    std::ofstream header;
    header.open(result_save_path, std::ios::app);
    header << endl;

    tuple<double, double, double, double> domain(-180.0, -90.0, 180.0, 90.0);

    set<METHODS> methods{ ARTREE, MRTREE };
    double Eabs = 100;
    double Erel = 0.01;
    bool DoRefinement = false;
    int repeat_times = 5;

    QueryResult qr;
    vector<double> predicted_results;
    for(int i = 1 ; i <= 4; i ++){
        LoadOSMAccumulationSurface2DSubset(accu_keys1, accu_keys2, accu_values, i);
        LoadOSMDataset2DSubset(keys1, keys2, i);
        //RunMethods_2D_Count<double, double, 3>(keys1, keys2, accu_keys1, accu_keys2, accu_values, ql1, ql2, qu1, qu2, methods, real_result_path, result_save_path, domain, Eabs, Erel, DoRefinement, repeat_times);

        auto t0 = chrono::steady_clock::now();
        Polyfit2D<double, double, 3> polyfit2d_deg3(100, 0.01, domain); // previously it's 1000
        polyfit2d_deg3.build(accu_keys1, accu_keys2, accu_values, i+3);
        polyfit2d_deg3.build_non_leaf();
        auto t1 = chrono::steady_clock::now();
        auto construction_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();

        auto avg_time = 0;
        for(int j  = 0 ; j < repeat_times; j++){
            predicted_results.clear();
            qr = polyfit2d_deg3.CountPrediction2D(ql1, ql2, qu1, qu2, predicted_results, DoRefinement, real_result_path);
            avg_time += qr.average_query_time;
        }
        avg_time /= repeat_times;
        qr.average_query_time = avg_time;

        header << "polyfit2d construction_time" << "," << construction_time << "," << " for degree" << "," << 3 << "," << "sgements" <<
        "," << polyfit2d_deg3.models.size() << " size opt" << "," << i << " avg query time" << "," << avg_time << endl;
    }
}

// = = = = = Selectivity = = = = =

void experiment_vary_selectivity(){
    //build once, multi queryset
    vector<int> keys, queryset_L, queryset_U;
    vector<int> values;
    LoadTweetDataset(keys, values);
    LoadTweetQuerySet(queryset_L, queryset_U);

    //set<METHODS> methods{PGM, RADIXSPLINE, POLYFIT, VERDICTDB, DBEST};
    set<METHODS> methods{ ARTREE, MRTREE, PGM, FITINGTREE, RMI, POLYFIT };
    string real_result_path = "/mnt/c/Users/Cloud/iCloudDrive/LearnedAggregate/VLDB_Final_Experiments/RealQueryResults/TWEET_1D.csv"; // in case of optimization
    string result_save_path = "/mnt/d/Polyfit_RunResults/EDBT2021/selectivity_1D_count.csv"; //

    double Eabs = 100;
    double Erel = 0.01;
    bool DoRefinement = false;
    int repeat_times = 10;

    for(int i = 1; i <= 4; i++) {
        // load dataset for different selectivity, 0.01%, 0.1%, 1%, 10%
        LoadTweetQuerySetWithSelectivity(queryset_L, queryset_U, i);
        RunMethods_1D_Count<int, int, 2>(keys, values, queryset_L, queryset_U, methods, real_result_path, result_save_path, Eabs, Erel, DoRefinement, repeat_times);
    }
}

void experiment_vary_selectivity_2D(){

    vector<double> keys1, keys2;
    LoadOSMDataset(keys1, keys2);
    vector<double> accu_keys1, accu_keys2, accu_values;

    vector<double> ql1, ql2, qu1, qu2, results;

    //set<METHODS> methods{PGM, RADIXSPLINE, POLYFIT, VERDICTDB, DBEST};
    set<METHODS> methods{ ARTREE };

    string real_result_path = "/mnt/c/Users/Cloud/iCloudDrive/LearnedAggregate/VLDB_Final_Experiments/RealQueryResults/OSM_2D.csv";
    string result_save_path = "/mnt/d/Polyfit_RunResults/EDBT2021/selectivity_2D_count.csv";

    tuple<double, double, double, double> domain(-180.0, -90.0, 180.0, 90.0);
    double Eabs = 100;
    double Erel = 0.01;
    bool DoRefinement = false;
    int repeat_times = 10;

    for(int i = 1; i <= 4; i++) {
        LoadOSMQuerySetWithSelectivity(ql1, ql2, qu1, qu2, i);
        RunMethods_2D_Count<double, double, 4>(keys1, keys2, accu_keys1, accu_keys2, accu_values, ql1, ql2, qu1, qu2, methods, real_result_path, result_save_path, domain, Eabs, Erel, DoRefinement, repeat_times);
    }
}

// = = = = = Others = = = = =

void experiment_VerdictDB_1D(){
    // test different ns and b

    vector<int> keys, values, queryset_L, queryset_U;
    LoadTweetDataset(keys, values);
    LoadTweetQuerySet(queryset_L, queryset_U);

    string real_result_path = "/mnt/c/Users/Cloud/iCloudDrive/LearnedAggregate/VLDB_Final_Experiments/RealQueryResults/TWEET_1D.csv";
    string result_save_path = "/mnt/d/Polyfit_RunResults/EDBT2021/VerdictDB_1D.csv";

//    std::ofstream header;
//    header.open(result_save_path, std::ios::app);
//    header << endl;

    int b = 100;
    for (int ns = 10000; ns < 100000; ns += 10000){
        RunVerdictDB<int, int>(queryset_L, queryset_U, keys, ns, b, 0.9, real_result_path, result_save_path, 1);
//        for (int b = 50; b < 200; b += 10){
//            RunVerdictDB<int, int>(queryset_L, queryset_U, keys, ns, b, 0.9, real_result_path, result_save_path, 1);
//        }
    }
}


void experiment_VerdictDB_2D(){
    // test different ns and b
    vector<double> keys1, keys2;
    LoadOSMDataset(keys1, keys2);
    vector<double> ql1, ql2, qu1, qu2, results;
    LoadOSMQuerySet(ql1, ql2, qu1, qu2);

    string real_result_path = "/mnt/c/Users/Cloud/iCloudDrive/LearnedAggregate/VLDB_Final_Experiments/RealQueryResults/OSM_2D.csv";
    string result_save_path = "/mnt/d/Polyfit_RunResults/EDBT2021/VerdictDB_2D.csv";

    QueryResult qr;
    int b = 10;
    for (int ns = 10000; ns < 10000000; ns += 10000){
        qr = {}; // reset struct
        auto avg_time = 0;
        for(int i  = 0 ; i < 3; i++){
            results.clear();
            qr = VerdictDBCountAggregate2D(keys1, keys2, ql1, ql2, qu1, qu2, results, real_result_path, ns, b);
            avg_time += qr.average_query_time;
        }
        avg_time /= 3;
        qr.average_query_time = avg_time;
        WriteQueryResult2File(qr, result_save_path, "VerdictDB_2D");
    }
}

void experiment_MRTree2D(){
    vector<double> keys1, keys2;
    LoadOSMDataset(keys1, keys2);
    vector<double> ql1, ql2, qu1, qu2, results;
    LoadOSMQuerySet(ql1, ql2, qu1, qu2);

    string real_result_path = "/mnt/c/Users/Cloud/iCloudDrive/LearnedAggregate/VLDB_Final_Experiments/RealQueryResults/OSM_2D.csv";
    string result_save_path = "/mnt/d/Polyfit_RunResults/EDBT2021/abserr_time_2D_count.csv";

    tuple<double, double, double, double> domain(-180.0, -90.0, 180.0, 90.0);
    vector<double> Eabs_collection = { 50, 100, 200, 500, 1000 };
    double Erel = 0.01;
    bool DoRefinement = false;
    int repeat_times = 3;

    MRTree<double, double, 2> mrt(100);
    mrt.build2DCount(keys1, keys2);

    QueryResult qr = {}; // reset struct

    for (auto Eabs : Eabs_collection){
        mrt.t_abs = Eabs;
        auto avg_time = 0;
        for(int i  = 0 ; i < repeat_times; i++){
            results.clear();
            qr = mrt.Query2DCount(ql1, ql2, qu1, qu2, results, real_result_path);
            avg_time += qr.average_query_time;
        }
        avg_time /= repeat_times;
        qr.average_query_time = avg_time;
        WriteQueryResult2File(qr, result_save_path, "MRTree_2D_"+std::to_string(Eabs));
    }
}

#endif //POLYFIT_CLION_EXPERIMENTS_H
