//
// Created by Cloud on 2020/8/17.
//

#ifndef POLYFIT_CLION_S2SAMPLING_H
#define POLYFIT_CLION_S2SAMPLING_H

#include <iostream>
#include <vector>
#include <cmath>
#include <random>

using namespace std;

double RationalApproximation(double t) {
    double c0 = 2.515517;
    double c1 = 0.802853;
    double c2 = 0.010328;
    double d1 = 1.432788;
    double d2 = 0.189269;
    double d3 = 0.001308;

    return t - ((c2 * t + c1)*t + c0) / (((d3 * t + d2)*t + d1)*t + 1.0);
}

// for normal distribution
double InverseCDF(double p) {
    if (p < 0.5)
    {
        // F^-1(p) = - G^-1(p)
        return -RationalApproximation(sqrt(-2.0*log(p)));
    }
    else
    {
        // F^-1(p) = G^-1(1-p)
        return RationalApproximation(sqrt(-2.0*log(1 - p)));
    }
}

// If using this one, two sample is 0 will terminate this algorithm
template<typename Tk>
double SequentialSampling(vector<Tk> keys, double lower, double upper, double p = 0.9, double Trel=0.01, double Tabs=100) {

    double Vn = 0; // used to substitue sigma^2
    double Xn_sum = 0;
    double Xn_average = 0; // used to substitue Mu
    double sampled_value = 0;
    double n = 1;
    double square_sum = 0;
    double Zp;
    srand((unsigned)time(NULL));

    double d = Tabs / keys.size(); // the absolute error of each partition

    // calculate Zp
    Zp = InverseCDF((1 + p) / 2);

    // perform the first sample
    double index = (rand() % (keys.size() - 1));
    if (keys[index] >= lower && keys[index] <= upper) {
        sampled_value = 1;
    }
    else {
        sampled_value = 0;
    }
    Xn_sum += sampled_value;
    n = 1;

    double pre = sampled_value;
    int max_index = 0;

    // used for large random number generation
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<unsigned long long> dis(0, keys.size());

    // start from n>=2
    while (true) {

        // perform sampling
        // random take 1 key
        //double index = (rand() % (keys.size()-1)); // by default, rand generate values from 0 to rand_max, rand max is equal or larger than 32767, in my pc is 32767
        int index = dis(gen) % (keys.size() - 1);
        /*if (index > max_index)
            max_index = index;*/

        if (keys[index] >= lower && keys[index] <= upper) {
            sampled_value = 1;
        }
        else {
            sampled_value = 0;
        }

        n++;
        Xn_sum += sampled_value;
        Xn_average = Xn_sum/n;

        // calculate Vn
        if (n == 2) {
            square_sum += (pre - Xn_average)*(pre - Xn_average);
        }
        square_sum += (sampled_value - Xn_average)*(sampled_value - Xn_average);
        Vn = square_sum / (n - 1);

        // check stop condition:

        double left_part = Xn_sum;
        if (n*d >= Xn_sum) {
            left_part = n * d;
        }

        double diff = Trel * left_part - Zp * sqrt(n*Vn);
        //cout << "diff: " << diff << endl;
        if(diff >= 0 && Vn > 0){
            break;
        }
        /*if (int(n) % 10000 == 0) {
            cout <<"current n: " << n  << "  diff: " << diff <<  "  Vn: " << Vn  << " index: " << index  << "  sum: " << Xn_sum << " max index: " << max_index << endl;
        }*/
    }
    //cout << "number of samplling: " << n << endl;
    double estimated_result = keys.size() * Xn_sum / n;
    return estimated_result;
}

// the 2D version
// keys1 has the same length as keys2, which is the same record's two keys
template<typename Tk>
double SequentialSampling2D(const vector<Tk> &keys1, const vector<Tk> &keys2, const Tk lower1, const Tk lower2, const Tk upper1, const Tk upper2,
                            double p = 0.9, double Trel = 0.01, double Tabs = 100) {

    double Vn = 0; // used to substitue sigma^2
    double Xn_sum = 0;
    double Xn_average = 0; // used to substitue Mu
    double sampled_value = 0;
    double n = 1;
    double square_sum = 0;
    double Zp;
    srand((unsigned)time(NULL));

    double d = Tabs / keys1.size(); // the absolute error of each partition

    // calculate Zp
    Zp = InverseCDF((1 + p) / 2);

    // perform the first sample
    double index = (rand() % (keys1.size() - 1));
    if (keys1[index] >= lower1 && keys1[index] <= upper1 && keys2[index] >= lower2 && keys2[index] <= upper2) {
        sampled_value = 1;
    }
    else {
        sampled_value = 0;
    }
    Xn_sum += sampled_value;
    n = 1;

    double pre = sampled_value;
    int max_index = 0;

    // used for large random number generation
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<unsigned long long> dis(0, keys1.size());

    // start from n>=2
    while (true) {

        // perform sampling
        // random take 1 key
        //double index = (rand() % (keys.size()-1)); // by default, rand generate values from 0 to rand_max, rand max is equal or larger than 32767, in my pc is 32767
        int index = dis(gen) % (keys1.size() - 1);
        /*if (index > max_index)
            max_index = index;*/

        if (keys1[index] >= lower1 && keys1[index] <= upper1 && keys2[index] >= lower2 && keys2[index] <= upper2) {
            sampled_value = 1;
        }
        else {
            sampled_value = 0;
        }

        n++;
        Xn_sum += sampled_value;
        Xn_average = Xn_sum / n;

        // calculate Vn
        if (n == 2) {
            square_sum += (pre - Xn_average)*(pre - Xn_average);
        }
        square_sum += (sampled_value - Xn_average)*(sampled_value - Xn_average);
        Vn = square_sum / (n - 1);

        // check stop condition:

        double left_part = Xn_sum;
        if (n*d >= Xn_sum) {
            left_part = n * d;
        }

        double diff = Trel * left_part - Zp * sqrt(n*Vn);
        //cout << "diff: " << diff << endl;
        if (diff >= 0 && Vn > 0) {
            break;
        }
        /*if (int(n) % 10000 == 0) {
            cout <<"current n: " << n  << "  diff: " << diff <<  "  Vn: " << Vn  << " index: " << index  << "  sum: " << Xn_sum << " max index: " << max_index << endl;
        }*/
    }
    //cout << "number of samplling: " << n << endl;
    double estimated_result = keys1.size() * Xn_sum / n;
    return estimated_result;
}


// The same, abandoned
template <typename Tk>
double S2Sampling(vector<Tk> keys, double lower, double upper, double p = 0.9, double Trel = 0.01) {
    int n = 1;
    double s = 0;
    double sample_value = 0;
    double w = 0;

    // perform the first sample
    double index = (rand() % (keys.size() - 1));
    if (keys[index] >= lower && keys[index] <= upper) {
        sample_value = 1;
    }
    else {
        sample_value = 0;
    }
    s = sample_value;

    double Zp = InverseCDF((1 + p) / 2);

    while (!(w >= 0) && !(Trel*s >= Zp * sqrt(n*w / (n - 1)))) {
        double index = (rand() % (keys.size() - 1));
        if (keys[index] >= lower && keys[index] <= upper) {
            sample_value = 1;
        }
        else {
            sample_value = 0;
        }

        w += (s - n * sample_value)*(s - n * sample_value) / (n*(n - 1));
        s += sample_value;
        n += 1;
    }
    cout << "number of samplling: " << n << endl;
    return keys.size() * s / n;
}

template<typename Tk, typename Tv>
QueryResult TestS2Sampling2D(const vector<Tk> &keys1, const vector<Tk> &keys2, const vector<Tk> &queryset_L1, const vector<Tk> &queryset_L2,
                             const vector<Tk> &queryset_U1, const vector<Tk> &queryset_U2, double p, double Trel,
                             double Tabs, vector<Tv> predicted_results, string RealResultPath){
    double result;
    predicted_results.clear();

    auto t0 = chrono::steady_clock::now();

    for (int i = 0; i < queryset_L1.size(); i++) {
        result = SequentialSampling2D(keys1, keys2, queryset_L1[i], queryset_L2[i], queryset_U1[i], queryset_U2[i], p, Trel, Tabs);
        predicted_results.push_back(int(result));
    }

    auto t1 = chrono::steady_clock::now();

    auto average_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count() / queryset_L1.size();
    auto total_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();

    double MEabs, MErel;
    MeasureAccuracy(predicted_results, RealResultPath, MEabs, MErel);

    QueryResult query_result;
    query_result.average_query_time = average_time;
    query_result.total_query_time = total_time;
    query_result.measured_absolute_error = MEabs;
    query_result.measured_relative_error = MErel;

    return query_result;

}

template<typename Tk, typename Tv>
QueryResult TestS2Sampling1D(const vector<Tk> &keys, const vector<Tk> &queryset_L, const vector<Tk> &queryset_U, double p, double Trel, double Tabs,
                             vector<Tv> predicted_results, string RealResultPath) {

    double result;
    predicted_results.clear();

    auto t0 = chrono::steady_clock::now();

    for (int i = 0; i < queryset_L.size(); i++) {
        result = SequentialSampling(keys, queryset_L[i], queryset_U[i], 0.9, Trel, Tabs);
        predicted_results.push_back(int(result));
    }

    auto t1 = chrono::steady_clock::now();

    auto average_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count() / queryset_L.size();
    auto total_time = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();

    double MEabs, MErel;

    // check correctness
    MeasureAccuracy(predicted_results, RealResultPath, MEabs, MErel);

    //ErrorGuaranteedSampling();

    QueryResult query_result;
    query_result.average_query_time = average_time;
    query_result.total_query_time = total_time;
    query_result.measured_absolute_error = MEabs;
    query_result.measured_relative_error = MErel;

    return query_result;
}

//void ErrorGuaranteedSampling() {
//    mat dataset;
//    bool loaded = mlpack::data::Load("C:/Users/Cloud/Desktop/LearnIndex/data/SortedSingleDimPOIs2.csv", dataset);
//    arma::rowvec trainingset = dataset.row(0);
//    arma::rowvec responses = dataset.row(dataset.n_rows - 1);
//    vector<double> keys;
//    RowvecToVector(trainingset, keys);
//
//    mat queryset;
//    bool loaded2 = mlpack::data::Load("C:/Users/Cloud/Desktop/LearnIndex/data/SortedSingleDimQuery2.csv", queryset);
//    arma::rowvec query_x_low = queryset.row(0);
//    arma::rowvec query_x_up = queryset.row(1);
//    vector<double> queryset_x_up_v, queryset_x_low_v;
//    RowvecToVector(query_x_up, queryset_x_up_v);
//    RowvecToVector(query_x_low, queryset_x_low_v);
//    vector<int> predicted_results, real_results;
//    vector<double> key_v;
//    RowvecToVector(trainingset, key_v);
//
//    double result;
//    auto t00 = chrono::steady_clock::now();
//    for (int i = 0; i < queryset_x_low_v.size(); i++) {
//        result = SequentialSampling(keys, queryset_x_low_v[i], queryset_x_up_v[i]);
//        //result = S2Sampling(keys, queryset_x_low_v[i], queryset_x_up_v[i]);
//        predicted_results.push_back(int(result));
//    }
//    auto t11 = chrono::steady_clock::now();
//    cout << "total query time in chrono: " << chrono::duration_cast<chrono::nanoseconds>(t11 - t00).count() << " in ns    " << chrono::duration_cast<chrono::nanoseconds>(t11 - t00).count() / (1000 * 1000 * 1000) << "in s" << endl;
//    cout << "average query time in chrono: " << chrono::duration_cast<chrono::nanoseconds>(t11 - t00).count() / queryset_x_up_v.size() << " in ns    " << chrono::duration_cast<chrono::nanoseconds>(t11 - t00).count() / queryset_x_up_v.size() / (1000 * 1000 * 1000) << "in s" << endl;
//
//    // check correctness
//    CalculateRealCountWithScan1D(queryset, real_results);
//    MeasureAccuracy(predicted_results, real_results);
//
//    system("pause");
//}

#endif //POLYFIT_CLION_S2SAMPLING_H
