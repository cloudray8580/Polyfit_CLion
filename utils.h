//
// Created by Cloud on 2020/8/15.
//

#ifndef POLYFIT_CLION_UTILS_H
#define POLYFIT_CLION_UTILS_H

#include <vector>
#include <string>
#include <armadillo>
#include <mlpack/core.hpp>
#include <string>

// = = = = = Return Struct = = = = =

struct Rect
{
    Rect() {}
    Rect(double a_minX, double a_minY, double a_maxX, double a_maxY)
    {
        min[0] = a_minX;
        min[1] = a_minY;
        max[0] = a_maxX;
        max[1] = a_maxY;
    }
    double min[2];
    double max[2];
};

// MRTREE = Multi-Resulution Tree
enum METHODS {RMI, FITINGTREE, PGM, RADIXSPLINE, POLYFIT, DBEST, VERDICTDB, S2, ENTROPYHIST, STX, ARTREE, PLATO, MRTREE};

struct QueryResult {
    //std::chrono::duration<double> average_query_time;
    //std::chrono::duration<double> total_query_time;
    unsigned long long average_query_time;
    unsigned long long total_query_time;
    double measured_absolute_error;
    double measured_relative_error;
    int hit_count;
    int refinement_count;
    int total_query_count;
    int tree_paras; // btree or Rtree
    int model_amount;
    int total_paras; // the total minimum paras
    int total_bytes; // total size in bytes
};

/**
 * save in csv format
 * @param qr
 * @param result_save_path
 */
void WriteQueryResult2File(QueryResult qr, std::string result_save_path, std::string tag="NAME"){
    std::ofstream run_result;
    run_result.open(result_save_path, std::ios::app);

    run_result << tag << "," << qr.average_query_time << "," << qr.total_query_time << ","
    << qr.measured_absolute_error << "," << qr.measured_relative_error << ","
    << qr.hit_count << "," << qr.refinement_count << "," << qr.total_query_count << ","
    << qr.tree_paras << "," << qr.model_amount << "," << qr.total_paras << "," << qr.total_bytes << std::endl;

    run_result.close();
}

// = = = = = Conversion = = = = =

template <typename T>
void VectorToRowvec(arma::mat& rv, const std::vector<T> &v) {
    rv.clear();
    rv.set_size(v.size());
    int count = 0;
    rv.imbue([&]() { return v[count++]; });
}

//template <typename T>
//void RowvecToVector(const arma::mat& matrix, std::vector<T> &vec) {
//    vec.clear();
//    for (int i = 0; i < matrix.n_cols; i++) {
//        vec.push_back(matrix[i]);
//    }
//}

template <typename T>
void RowvecToVector(const arma::rowvec& matrix, std::vector<T> &vec) {
    vec.clear();
    for (int i = 0; i < matrix.n_cols; i++) {
        vec.push_back(matrix[i]);
    }
}

// = = = = = Dataset and Queryset = = = = =

template <typename Tk, typename Tv>
void LoadTweetDataset(std::vector<Tk> &keys, std::vector<Tv> &values) {
    arma::mat dataset;
    //bool loaded = mlpack::data::Load("/mnt/c/Users/Cloud/Desktop/LearnIndex/data/SortedSingleDimPOIs2.csv", dataset); // original
    bool loaded = mlpack::data::Load("/mnt/d/Polyfit_Dataset/TWEET/1M_keys_toLargeInt.csv", dataset); // transformed to large int (*10^5)
    arma::rowvec trainingset = dataset.row(0);
    arma::rowvec responses = dataset.row(dataset.n_rows - 1);
    RowvecToVector(trainingset, keys);
    RowvecToVector(responses, values);
}

template <typename Tk>
void LoadTweetQuerySet(std::vector<Tk> &Querykey_L, std::vector<Tk> &Querykey_U) {
    arma::mat queryset;
    //bool loaded2 = mlpack::data::Load("/mnt/c/Users/Cloud/Desktop/LearnIndex/data/SortedSingleDimQuery2.csv", queryset); // original
    bool loaded = mlpack::data::Load("/mnt/d/Polyfit_Dataset/TWEET/1M_queryset_toLargeInt.csv", queryset); // transformed to large int (*10^5)
    arma::rowvec query_x_low = queryset.row(0);
    arma::rowvec query_x_up = queryset.row(1);
    std::vector<double> queryset_x_up_v, queryset_x_low_v;
    RowvecToVector(query_x_low, Querykey_L);
    RowvecToVector(query_x_up, Querykey_U);
}

/**
 * dedicated designed for RadixSpline
 * @tparam Tk
 * @param vec : the original vector
 * @param unvec : the unsigned vector to 'return'
 * @param base : add each original element with base
 */
template <typename Tk>
void Convert2Unsign(const std::vector<Tk> &vec, std::vector<uint64_t> &unvec, uint64_t base){
    unvec.clear();
    unvec.reserve(vec.size());
    for (int i = 0 ; i < vec.size(); i++){
        unvec.push_back((uint64_t)(vec[i] + base));
    }
}

// original OSM, for 2D
template <typename Tk>
void LoadOSMDataset(std::vector<Tk> &keys1, std::vector<Tk> &keys2) {
    arma::mat dataset;
    mlpack::data::Load("/mnt/d/Polyfit_Dataset/OSM/100M_Long_Lat_all.csv", dataset); // first key longitude, second key latitude
    arma::rowvec key1_row = dataset.row(0);
    arma::rowvec key2_row = dataset.row(1);
    RowvecToVector(key1_row, keys1);
    RowvecToVector(key2_row, keys2);
}

// the processed accumulation surface
template <typename Tk, typename Tv>
void LoadOSMAccumulationSurface(std::vector<Tk> &keys1, std::vector<Tk> &keys2, std::vector<Tv> &values){
    arma::mat dataset;
    mlpack::data::Load("/mnt/d/Polyfit_Dataset/OSM/100M_Long_Lat_accu.csv", dataset); // lat lon
    arma::rowvec key1_row = dataset.row(0);
    arma::rowvec key2_row = dataset.row(1);
    arma::rowvec vals_row = dataset.row(2);
    RowvecToVector(key1_row, keys1);
    RowvecToVector(key2_row, keys2);
    RowvecToVector(vals_row, values);
}

template <typename Tk>
void LoadOSMQuerySet(std::vector<Tk> &d1_low, std::vector<Tk> &d2_low, std::vector<Tk> &d1_up, std::vector<Tk> &d2_up) {
    arma::mat queryset;
    mlpack::data::Load("/mnt/d/Polyfit_Dataset/OSM/100M_Long_Lat_queryset.csv", queryset); // d1_low, d2_low, d1_up, d2_up
    //mlpack::data::Load("/mnt/d/Polyfit_Dataset/OSM/100M_Long_Lat_queryset_withnoise.csv", queryset);
    RowvecToVector(queryset.row(0), d1_low);
    RowvecToVector(queryset.row(1), d2_low);
    RowvecToVector(queryset.row(2), d1_up);
    RowvecToVector(queryset.row(3), d2_up);
}

// HKI for Max
template <typename Tk, typename Tv>
void LoadHKIDataset(std::vector<Tk> &keys, std::vector<Tv> &values) {
    arma::mat dataset;
    bool loaded = mlpack::data::Load("/mnt/d/Polyfit_Dataset/HKI/1M_all.csv", dataset);
    arma::rowvec trainingset = dataset.row(0);
    arma::rowvec responses = dataset.row(2);
    RowvecToVector(trainingset, keys);
    RowvecToVector(responses, values);
}

template <typename Tk>
void LoadHKIQuerySet(std::vector<Tk> &Querykey_L, std::vector<Tk> &Querykey_U) {
    arma::mat queryset;
    bool loaded2 = mlpack::data::Load("/mnt/d/Polyfit_Dataset/HKI/1M_queryset.csv", queryset);
    arma::rowvec query_x_low = queryset.row(0);
    arma::rowvec query_x_up = queryset.row(1);
    RowvecToVector(query_x_up, Querykey_U);
    RowvecToVector(query_x_low, Querykey_L);
}

// For varying dataset size, first key latitude
template <typename Tk, typename Tv>
void LoadOSM_Dataset_Queryset_1D_SUBSET(std::vector<Tk> &keys, std::vector<Tv> &values, std::vector<Tk> &q_low, std::vector<Tk> &q_up, int option=1, int toint = false) {

    keys.clear();
    values.clear();
    q_low.clear();
    q_up.clear();

    arma::mat dataset;
    arma::mat queryset;

    if(toint){
        switch (option){
            case 1:
                //mlpack::data::Load("/mnt/d/Polyfit_Dataset/MAP/1M_dataset_toLargeInt.csv", dataset);
                mlpack::data::Load("/mnt/d/Polyfit_Dataset/MAP/1M_dataset_toLargeInt_deduplicated.csv", dataset);
                mlpack::data::Load("/mnt/d/Polyfit_Dataset/MAP/1M_queryset_toLargeInt.csv", queryset);
                break;
            case 2:
                //mlpack::data::Load("/mnt/d/Polyfit_Dataset/MAP/10M_dataset_toLargeInt.csv", dataset);
                mlpack::data::Load("/mnt/d/Polyfit_Dataset/MAP/10M_dataset_toLargeInt_deduplicated.csv", dataset);
                mlpack::data::Load("/mnt/d/Polyfit_Dataset/MAP/10M_queryset_toLargeInt.csv", queryset);
                break;
            case 3:
                //mlpack::data::Load("/mnt/d/Polyfit_Dataset/MAP/30M_dataset_toLargeInt.csv", dataset);
                mlpack::data::Load("/mnt/d/Polyfit_Dataset/MAP/30M_dataset_toLargeInt_deduplicated.csv", dataset);
                mlpack::data::Load("/mnt/d/Polyfit_Dataset/MAP/30M_queryset_toLargeInt.csv", queryset);
                break;
            case 4:
                //mlpack::data::Load("/mnt/d/Polyfit_Dataset/MAP/100M_dataset_toLargeInt.csv", dataset);
                mlpack::data::Load("/mnt/d/Polyfit_Dataset/MAP/100M_dataset_toLargeInt_deduplicated.csv", dataset);
                mlpack::data::Load("/mnt/d/Polyfit_Dataset/MAP/100M_queryset_toLargeInt.csv", queryset);
                break;
            default:
                mlpack::data::Load("/mnt/d/Polyfit_Dataset/MAP/1M_dataset_toLargeInt.csv", dataset);
                mlpack::data::Load("/mnt/d/Polyfit_Dataset/MAP/1M_queryset_toLargeInt.csv", queryset);
                break;
        }
    } else {
        switch (option){
            case 1:
                //mlpack::data::Load("/mnt/c/Users/Cloud/Desktop/LearnedAggregateData/MapData_1M_Sorted_Value.csv", dataset);
                mlpack::data::Load("/mnt/d/Polyfit_Dataset/MAP/1M_dataset_deduplicated.csv", dataset);
                mlpack::data::Load("/mnt/c/Users/Cloud/Desktop/LearnedAggregateData/MapData_1M_Query_1D.csv", queryset);
                break;
            case 2:
                //mlpack::data::Load("/mnt/c/Users/Cloud/Desktop/LearnedAggregateData/MapData_10M_Sorted_Value.csv", dataset);
                mlpack::data::Load("/mnt/d/Polyfit_Dataset/MAP/10M_dataset_deduplicated.csv", dataset);
                mlpack::data::Load("/mnt/c/Users/Cloud/Desktop/LearnedAggregateData/MapData_10M_Query_1D.csv", queryset);
                break;
            case 3:
                //mlpack::data::Load("/mnt/c/Users/Cloud/Desktop/LearnedAggregateData/MapData_30M_Sorted_Value.csv", dataset);
                mlpack::data::Load("/mnt/d/Polyfit_Dataset/MAP/30M_dataset_deduplicated.csv", dataset);
                mlpack::data::Load("/mnt/c/Users/Cloud/Desktop/LearnedAggregateData/MapData_30M_Query_1D.csv", queryset);
                break;
            case 4:
                //mlpack::data::Load("/mnt/c/Users/Cloud/Desktop/LearnedAggregateData/MapData_100M_Sorted_Value.csv", dataset);
                mlpack::data::Load("/mnt/d/Polyfit_Dataset/MAP/100M_dataset_deduplicated.csv", dataset);
                mlpack::data::Load("/mnt/c/Users/Cloud/Desktop/LearnedAggregateData/MapData_100M_Query_1D.csv", queryset);
                break;
            default:
                mlpack::data::Load("/mnt/c/Users/Cloud/Desktop/LearnedAggregateData/MapData_1M_Sorted_Value.csv", dataset);
                mlpack::data::Load("/mnt/c/Users/Cloud/Desktop/LearnedAggregateData/MapData_1M_Query_1D.csv", queryset);
                break;
        }
    }

    arma::rowvec key_row = dataset.row(0);
    arma::rowvec value_row = dataset.row(dataset.n_rows - 1);
    RowvecToVector(key_row, keys);
    RowvecToVector(value_row, values);
    RowvecToVector(queryset.row(0), q_low);
    RowvecToVector(queryset.row(1), q_up);
}

template <typename Tk, typename Tv>
void LoadOSMAccumulationSurface2DSubset(std::vector<Tk> &keys1, std::vector<Tk> &keys2, std::vector<Tv> &values, int option = 1){
    arma::mat dataset;

    switch (option){
        case 1:
            mlpack::data::Load("/mnt/d/Polyfit_Dataset/OSM/100M_Long_Lat_accu.csv", dataset); // lat lon
            break;
        case 2:
            mlpack::data::Load("/mnt/d/Polyfit_Dataset/OSM/100M_Long_Lat_accu_9M.csv", dataset); // lat lon
            break;
        case 3:
            mlpack::data::Load("/mnt/d/Polyfit_Dataset/OSM/100M_Long_Lat_accu_25M.csv", dataset); // lat lon
            break;
        case 4:
            mlpack::data::Load("/mnt/d/Polyfit_Dataset/OSM/100M_Long_Lat_accu_100M.csv", dataset); // lat lon
            break;
    }

    arma::rowvec key1_row = dataset.row(0);
    arma::rowvec key2_row = dataset.row(1);
    arma::rowvec vals_row = dataset.row(2);
    RowvecToVector(key1_row, keys1);
    RowvecToVector(key2_row, keys2);
    RowvecToVector(vals_row, values);
}

// 1: 1M, 2: 10M, 3: 30M, 4: 100M
template <typename Tk>
void LoadOSMDataset2DSubset(std::vector<Tk> &keys1, std::vector<Tk> &keys2, int option) {
    arma::mat dataset;
    mlpack::data::Load("/mnt/d/Polyfit_Dataset/OSM/100M_Long_Lat_subset_"+std::to_string(option)+".csv", dataset); // first key longitude, second key latitude
    arma::rowvec key1_row = dataset.row(0);
    arma::rowvec key2_row = dataset.row(1);
    RowvecToVector(key1_row, keys1);
    RowvecToVector(key2_row, keys2);
}


// for accumulative sum
template <typename Tk, typename Tv>
void LoadHKIAccumulativeSumDataset(std::vector<Tk> &keys, std::vector<Tv> &values) {
    arma::mat dataset;
    bool loaded = mlpack::data::Load("/mnt/c/Users/Cloud/iCloudDrive/SampledFinancialAccu.csv", dataset);
    arma::rowvec trainingset = dataset.row(0);
    arma::rowvec responses = dataset.row(1);
    RowvecToVector(trainingset, keys);
    RowvecToVector(responses, values);
}

// 1: 0.01%, 2: 0.1%, 3: 1%, 4: 10%
template <typename Tk>
void LoadTweetQuerySetWithSelectivity(std::vector<Tk> &Querykey_L, std::vector<Tk> &Querykey_U, int option) {
    arma::mat queryset;
    bool loaded = mlpack::data::Load("/mnt/d/Polyfit_Dataset/TWEET/1M_queryset_toLargeInt_selectivity_"+std::to_string(option)+".csv", queryset); // large int (*10^5)
    arma::rowvec query_x_low = queryset.row(0);
    arma::rowvec query_x_up = queryset.row(1);
    std::vector<double> queryset_x_up_v, queryset_x_low_v;
    RowvecToVector(query_x_low, Querykey_L);
    RowvecToVector(query_x_up, Querykey_U);
}

// 1: 0.01%, 2: 0.1%, 3: 1%, 4: 10%
template <typename Tk>
void LoadOSMQuerySetWithSelectivity(std::vector<Tk> &d1_low, std::vector<Tk> &d2_low, std::vector<Tk> &d1_up, std::vector<Tk> &d2_up, int option) {
    arma::mat queryset;
    mlpack::data::Load("/mnt/d/Polyfit_Dataset/OSM/100M_Long_Lat_queryset_selectivity_"+std::to_string(option)+".csv", queryset); // d1_low, d2_low, d1_up, d2_up
    //mlpack::data::Load("/mnt/d/Polyfit_Dataset/OSM/100M_Long_Lat_queryset_withnoise.csv", queryset);
    RowvecToVector(queryset.row(0), d1_low);
    RowvecToVector(queryset.row(1), d2_low);
    RowvecToVector(queryset.row(2), d1_up);
    RowvecToVector(queryset.row(3), d2_up);
}

// = = = = = Query Measurement = = = = =

/**
 * this is only for CountPrediction
 * @param result
 * @return true if no refinement need
 */
template <typename Tv>
bool inline RelativeErrorCheck_1D_Count(Tv result, double t_abs, double t_rel){
    double max_err_rel = (2 * t_abs) / (result - 2 * t_abs);
    if (max_err_rel > t_rel || max_err_rel < 0)
        return false;
    return true;
}

template <typename Tv>
bool inline RelativeErrorCheck_2D_Count(Tv result, double t_abs, double t_rel){
    double max_err_rel = 4 * t_abs / (result - 4 * t_abs);
    if (max_err_rel > t_rel || max_err_rel < 0)
        return false;
    return true;
}

template <typename Tv>
bool inline RelativeErrorCheck_1D_Max(Tv result, double t_abs, double t_rel){
    double max_err_rel = t_abs / (result - t_abs);
    if (max_err_rel > t_rel || max_err_rel < 0)
        return false;
    return true;
}

template <typename Tv>
void MeasureAccuracy(std::vector<Tv> &predicted_results, std::string filepath, double &MEabs, double &MErel) {

    arma::mat results;
//    results.load(filepath);
    bool loaded = mlpack::data::Load(filepath, results);
    arma::rowvec real_results_row = results.row(0);
    std::vector<double> real_results;
    RowvecToVector(real_results_row, real_results);

    double absolute_error = 0;
    double relative_error;
    double accu_relative = 0;
    double accu_absolute = 0;
    double est_rel_err = 0;
    int total_size = predicted_results.size();
    for (int i = 0; i < total_size; i++) {
        if (real_results[i] == 0) {
            total_size--;
            continue;
        }
        absolute_error = abs(predicted_results[i] - real_results[i]);
        relative_error = abs(double(predicted_results[i] - real_results[i])) / abs(real_results[i]);
        accu_absolute += absolute_error;
        accu_relative += relative_error;
    }

    MEabs = accu_absolute / total_size;
    MErel = accu_relative / total_size;
    //cout << "measured average relative error: " << accu / total_size << endl;
    //cout << "measured average absolute error: " << accu_absolute / total_size << endl;
}

#endif //POLYFIT_CLION_UTILS_H
