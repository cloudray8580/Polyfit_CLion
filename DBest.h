//
// Created by Cloud on 2020/8/20.
//

#ifndef POLYFIT_CLION_DBEST_H
#define POLYFIT_CLION_DBEST_H


#include <cmath>
#include <xgboost/c_api.h>
#include "utils.h"
#define safe_xgboost(call) {                                                                    \
int err = (call);                                                                               \
if (err != 0) {                                                                                 \
  fprintf(stderr, "%s:%d: error in %s: %s\n", __FILE__, __LINE__, #call, XGBGetLastError());    \
  exit(1);                                                                                      \
}                                                                                               \
}

QueryResult DBestCountPrediction(const char* model_uri, float* query_low, float* query_up, int querysize, string RealResultPath){

    // load the model
    BoosterHandle booster;
    XGBoosterCreate(NULL, 0, &booster);
    safe_xgboost(XGBoosterLoadModel(booster, model_uri));

    // convert query L U into the required format
    DMatrixHandle d_l, d_u;
    safe_xgboost(XGDMatrixCreateFromMat(query_low, querysize, 1, 0, &d_l));
    safe_xgboost(XGDMatrixCreateFromMat(query_up, querysize, 1, 0, &d_u));

    auto t0 = chrono::steady_clock::now();
    // prediction
    bst_ulong out_len_l = 0, out_len_u = 0;
    const float* out_result_l = NULL;
    const float* out_result_u = NULL;
    safe_xgboost(XGBoosterPredict(booster, d_l, 0, 0, 0, &out_len_l, &out_result_l));
    std::vector<float> pred_l(out_result_l, out_result_l + out_len_l);

    safe_xgboost(XGBoosterPredict(booster, d_u, 0, 0, 0, &out_len_u, &out_result_u));
    std::vector<float> pred_u(out_result_u, out_result_u + out_len_u);

    // calculate the count
    std::vector<float> results(querysize);
    for (int i = 0; i < querysize; ++i) {
        results[i] = pred_u[i] - pred_l[i];
    }

    auto t1 = chrono::steady_clock::now();

    QueryResult query_result;

    // measure query performance
    auto average_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count() / querysize;
    auto total_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
    double MEabs, MErel;
    MeasureAccuracy(results, RealResultPath, MEabs, MErel);

    query_result.average_query_time = average_time;
    query_result.total_query_time = total_time;
    query_result.measured_absolute_error = MEabs;
    query_result.measured_relative_error = MErel;

    return query_result;
}

template<typename Tk>
QueryResult DBestCountPrediction2D(const char* model_uri, const vector<Tk> query_low_1, const vector<Tk> query_low_2, const vector<Tk> query_up_1, const vector<Tk> query_up_2,
                                   int querysize, string RealResultPath){
    // load the model
    BoosterHandle booster;
    XGBoosterCreate(NULL, 0, &booster);
    safe_xgboost(XGBoosterLoadModel(booster, model_uri));

    // convert float** to float* produce incorrect logic, but you can conver float[][] to float* correctly
//    float **q1 = new float*[querysize];
//    float **q2 = new float*[querysize];
//    float **q3 = new float*[querysize];
//    float **q4 = new float*[querysize];
//    for(int i = 0; i < querysize; i++){
//        q1[i] = new float[2]{query_low_1[i], query_low_2[i]};
//        q2[i] = new float[2]{query_up_1[i], query_up_2[i]};
//        q3[i] = new float[2]{query_low_1[i], query_up_2[i]};
//        q4[i] = new float[2]{query_up_1[i], query_low_2[i]};
//    }

    float *q1 = new float[querysize*2];
    float *q2 = new float[querysize*2];
    float *q3 = new float[querysize*2];
    float *q4 = new float[querysize*2];
    // first row, then columns
    for(int i = 0; i < querysize; i++){
        q1[2*i] = query_low_1[i];
        q1[2*i+1] = query_low_2[i];
        q2[2*i] = query_up_1[i];
        q2[2*i+1] = query_up_2[i];
        q3[2*i] = query_low_1[i];
        q3[2*i+1] = query_up_2[i];
        q4[2*i] = query_up_1[i];
        q4[2*i+1] = query_low_2[i];
    }

    DMatrixHandle d_1, d_2, d_3, d_4;
    safe_xgboost(XGDMatrixCreateFromMat(q1, querysize, 2, 0, &d_1));
    safe_xgboost(XGDMatrixCreateFromMat(q2, querysize, 2, 0, &d_2));
    safe_xgboost(XGDMatrixCreateFromMat(q3, querysize, 2, 0, &d_3));
    safe_xgboost(XGDMatrixCreateFromMat(q4, querysize, 2, 0, &d_4));

    auto t0 = chrono::steady_clock::now();
    // prediction
    bst_ulong out_len_1 = 0, out_len_2 = 0, out_len_3 = 0, out_len_4 = 0;
    const float* out_result_1 = NULL;
    const float* out_result_2 = NULL;
    const float* out_result_3 = NULL;
    const float* out_result_4 = NULL;
    safe_xgboost(XGBoosterPredict(booster, d_1, 0, 0, 0, &out_len_1, &out_result_1));
    std::vector<float> pred_1(out_result_1, out_result_1 + out_len_1);

    safe_xgboost(XGBoosterPredict(booster, d_2, 0, 0, 0, &out_len_2, &out_result_2));
    std::vector<float> pred_2(out_result_2, out_result_2 + out_len_2);

    safe_xgboost(XGBoosterPredict(booster, d_3, 0, 0, 0, &out_len_3, &out_result_3));
    std::vector<float> pred_3(out_result_3, out_result_3 + out_len_3);

    safe_xgboost(XGBoosterPredict(booster, d_4, 0, 0, 0, &out_len_4, &out_result_4));
    std::vector<float> pred_4(out_result_4, out_result_4 + out_len_4);

    // calculate the count
    std::vector<float> results(querysize);
    for (int i = 0; i < querysize; ++i) {
        results[i] = pred_1[i] + pred_2[i] - pred_3[i] - pred_4[i];
    }

    auto t1 = chrono::steady_clock::now();

    QueryResult query_result;

    // measure query performance
    auto average_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count() / querysize;
    auto total_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
    double MEabs, MErel;
    MeasureAccuracy(results, RealResultPath, MEabs, MErel);

    query_result.average_query_time = average_time;
    query_result.total_query_time = total_time;
    query_result.measured_absolute_error = MEabs;
    query_result.measured_relative_error = MErel;

    return query_result;
}

// ref: https://github.com/dmlc/xgboost/blob/master/demo/c-api/c-api-demo.c
void TestXGBoost(){
    int silent = 1; // print during loading, 0: not print

    // load the data
    DMatrixHandle keys;
    safe_xgboost(XGDMatrixCreateFromFile("/mnt/d/Polyfit_Dataset/TWEET/1M_keys.csv?format=csv", silent, &keys));

    // load the model
    BoosterHandle booster;
    XGBoosterCreate(NULL, 0, &booster);
    safe_xgboost(XGBoosterLoadModel(booster, "/mnt/c/Users/Cloud/iCloudDrive/Polyfit_Python/model.bin"));

    // predict
    bst_ulong out_len = 0;
    const float* out_result = NULL;
    int n_print = 10;

    auto t0 = chrono::steady_clock::now();
    safe_xgboost(XGBoosterPredict(booster, keys, 0, 0, 0, &out_len, &out_result));
    auto t1 = chrono::steady_clock::now();
    printf("y_pred: ");
    for (int i = 0; i < n_print; ++i) {
        printf("%1.4f ", out_result[i]);
    }

    printf("\n");
    printf("processing time: %d", t1-t0);
}

#endif //POLYFIT_CLION_DBEST_H
