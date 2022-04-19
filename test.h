//
// Created by Cloud on 2020/8/27.
//

#ifndef POLYFIT_CLION_TEST_H
#define POLYFIT_CLION_TEST_H

#include <vector>
#include "Polyfit.h"

using namespace std;

void BuildSingleModel(){
    vector<double> keys, values;
    LoadHKIAccumulativeSumDataset(keys, values);

    Polyfit<double, double, 3> polyfit(20000, 0.01);

    vector<double> paras;
    double loss;

    polyfit.SolveMaxlossLP(keys, values, 0, keys.size(), paras, loss);

    cout.precision(11);
    cout << "max loss: " << loss << endl;
    for (int i = 0; i < paras.size(); i++) { // start from a0
        cout << paras[i] << " ";
    }
}

//template <int deg, typename T>
//class Test{
//public:
//    Test(vector<T> &aa){
//        this->a = aa;
//        cout << a[0] << endl;
//    }
//
//    vector<T> a;
//};
//
//template <int deg, typename T>
//void tt(T key){
//    cout << key << endl;
//}

#endif //POLYFIT_CLION_TEST_H
