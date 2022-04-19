//
// Created by Cloud on 2020/8/28.
//

#ifndef POLYFIT_CLION_POLYFIT2D_H
#define POLYFIT_CLION_POLYFIT2D_H

#include <vector>
#include <array>
#include "ilcplex/ilocplex.h"
#include "utils.h"
#include "SimpleRTree.h"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<double, 2, bg::cs::cartesian> point;
typedef bg::model::box<point> box;
typedef std::pair<box, unsigned> value;

using namespace std;

/**
 * Polyfit 1D
 * @tparam Tk : type of key
 * @tparam Tv : type of value
 * @tparam degree : the highest degree
 */
template<typename Tk, typename Tv, int degree = 2>
class Polyfit2D {
public:

    /**
     *
     * @param t_abs
     * @param t_rel
     * @param boundary
     */
    Polyfit2D(Tv t_abs, double t_rel, tuple<Tk, Tk, Tk, Tk> domain): t_abs(t_abs), t_rel(t_rel), global_domain(domain) {}

    inline bool inside(const std::tuple<Tk, Tk, Tk, Tk> &domain, const Tk &k1, const Tk &k2){
        if(k1 >= get<0>(domain) and k1 < get<2>(domain) and k2 >= get<1>(domain) and k2 < get<3>(domain))
            return true;
        return false;
    }

    int called_times = 0;

    void build_subregion(const vector<Tk> &key1_v, const vector<Tk> &key2_v, const vector<Tv> &value_v,
                         std::tuple<Tk, Tk, Tk, Tk> domain, vector<int> indexes, int start_level = 0, int level = 0){
        double loss = t_abs + 1;
        vector<double> paras;
        if(level >= start_level){
            SolveMaxlossLP2D(key1_v, key2_v, value_v, std::get<0>(domain), std::get<1>(domain), indexes, paras, loss);
//            cout << "call times: " << called_times++ << " processing level " << level << "  domain: " << get<0>(domain) << "," << get<1>(domain) << ","
//                 << get<2>(domain) << "," << get<3>(domain) << "  indexes size: " << indexes.size() << "  loss: " << loss <<
//                 "  paras: " << paras[0] << "," << paras[1] << "," << paras[2] << "," << paras[3] << "," << paras[4] << "..." << endl;
        }

        if(loss <= t_abs){
            Model model(get<0>(domain), get<1>(domain), get<2>(domain), get<3>(domain), paras);
            models.push_back(model);
        } else {
            Tk d1_middle = (get<0>(domain) + get<2>(domain)) / 2;
            Tk d2_middle = (get<1>(domain) + get<3>(domain)) / 2;
            tuple<Tk, Tk, Tk, Tk> domain1(get<0>(domain), get<1>(domain), d1_middle, d2_middle); // lower left
            tuple<Tk, Tk, Tk, Tk> domain2(get<0>(domain), d2_middle, d1_middle, get<3>(domain)); // upper left
            tuple<Tk, Tk, Tk, Tk> domain3(d1_middle, get<1>(domain), get<2>(domain), d2_middle); // lower right
            tuple<Tk, Tk, Tk, Tk> domain4(d1_middle, d2_middle, get<2>(domain), get<3>(domain)); // upper right
            vector<int> indexes1, indexes2, indexes3, indexes4;
            // distributes points to each index
            for (int i : indexes){
                if (inside(domain1, key1_v[i], key2_v[i]))
                    indexes1.push_back(i);
                else if (inside(domain2, key1_v[i], key2_v[i]))
                    indexes2.push_back(i);
                else if (inside(domain3, key1_v[i], key2_v[i]))
                    indexes3.push_back(i);
                else if (inside(domain4, key1_v[i], key2_v[i]))
                    indexes4.push_back(i);
            }
            if(indexes1.size() > 0) build_subregion(key1_v, key2_v, value_v, domain1, indexes1, start_level,level+1);
            if(indexes2.size() > 0) build_subregion(key1_v, key2_v, value_v, domain2, indexes2, start_level,level+1);
            if(indexes3.size() > 0) build_subregion(key1_v, key2_v, value_v, domain3, indexes3, start_level,level+1);
            if(indexes4.size() > 0) build_subregion(key1_v, key2_v, value_v, domain4, indexes4, start_level,level+1);
        }
    }

    void build(const vector<Tk> &key1_v, const vector<Tk> &key2_v, const vector<Tv> &value_v, int start_level = 0){
        vector<int> indexes(value_v.size());
        std::iota(indexes.begin(), indexes.end(), 0);
        build_subregion(key1_v, key2_v, value_v, global_domain, indexes, start_level);
        cout << "total models: " << models.size() << " with degree " << degree << endl;
    }

    /**
     * build the non-leaf layer, assume the bottom layer segments have been built
     */
    void build_non_leaf(){
        // boost rtree
        rtree.clear();
        for ( unsigned i = 0 ; i < models.size() ; i++ ) {
            // create a box
            box b(point(models[i].key_L1, models[i].key_L2), point(models[i].key_U1, models[i].key_U2)); // lower point and upper point
            // insert new value
            rtree.insert(std::make_pair(b, i));
        }
    }

    /**
     * check whether the KNN return the correct result and compare which methods is faster
     * Result: KNN is faster
     */
    void TestKNNAndIntersectPerformance(){
        box query_box(point(0, 0), point(0, 0));
        std::vector<value> result_s;
        auto t0 = chrono::steady_clock::now();
        rtree.query(bgi::intersects(query_box), std::back_inserter(result_s));
        auto t1 = chrono::steady_clock::now();
        cout << "query time using intersection: " <<  chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count() << "   box id: " << result_s[0].second << endl;

        std::vector<value> result_n;
        auto t2 = chrono::steady_clock::now();
        rtree.query(bgi::nearest(point(0, 0), 1), std::back_inserter(result_n));
        auto t3 = chrono::steady_clock::now();
        cout << "query time using KNN: " <<  chrono::duration_cast<chrono::nanoseconds>(t3 - t2).count() << "   box id: " << result_n[0].second << endl;

    }

    void SolveMaxlossLP2D(const vector<Tk> &key1_v, const vector<Tk> &key2_v, const vector<Tv> &value_v, Tk base_key1, Tk base_key2,
                          const vector<int> &indexes, vector<double> &paras, double &loss) {
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

        cplex.setParam(IloCplex::RootAlg, IloCplex::Primal); // using simplex
        //cplex.setParam(IloCplex::RootAlg, IloCplex::Dual); // using dual simplex
        //cplex.setParam(IloCplex::RootAlg, IloCplex::Barrier); // set optimizer used interior point method

        //cplex.setParam(IloCplex::ObjLLim, 100);
        cplex.setParam(IloCplex::Param::TimeLimit, 20); // set the maximum time to solve a segment to 20 seconds

        // set variable type, IloNumVarArray starts from 0.
        for (int i = 0; i < (degree+1)*(degree+1); i++) {
            vars.add(IloNumVar(env, -INFINITY, INFINITY, ILOFLOAT));  // the weights, i.e., a_0_0, a_0_1, .. a_0_n, a_1_0, .. a_1_n, ... , a_n_0, ... a_n_n, n is deg
        }
        IloNumVar target(env, 0.0, INFINITY, ILOFLOAT); // the max loss value

        // declare objective
        obj.setExpr(target);
        obj.setSense(IloObjective::Minimize);
        model.add(obj);

        vector<Tk> key_term_1, key_term_2;
        Tk key_1, key_2; // x
        Tk key_diff_1, key_diff_2; // x - x0
        Tk current_key_term_1, current_key_term_2; // (x - x0)^n

        // add constraint for each record
        for (auto i : indexes) {
            key_1 = key1_v[i];
            key_2 = key2_v[i];
            key_diff_1 = key_1 - base_key1;
            key_diff_2 = key_2 - base_key2;
            current_key_term_1 = 1.0;
            current_key_term_2 = 1.0;
            key_term_1.clear();
            key_term_2.clear();
            key_term_1.push_back(current_key_term_1);
            key_term_2.push_back(current_key_term_2);
            // 0, (x-x0), (x-x0)^2, ... , (x-x0)^n
            for (int j = 1; j <= degree; j++) {
                current_key_term_1 *= key_diff_1;
                current_key_term_2 *= key_diff_2;
                key_term_1.push_back(current_key_term_1);
                key_term_2.push_back(current_key_term_2);
            }

            // a0*0 + a1*(x-x0) + a2*(x-x0)^2 + ... + an*(x-x0)^n
            IloExpr model_term(env);
            for(int m = 0 ; m < degree + 1; m++){
                for(int n = 0; n < degree + 1; n++){
                    model_term += vars[m*(degree+1)+n] * key_term_1[m] * key_term_2[n];
                }
            }

            model.add(model_term - value_v[i] <= target);
            model.add(model_term - value_v[i] >= -target);
        }

        IloNum starttime_ = cplex.getTime();
        cplex.solve();
        IloNum endtime_ = cplex.getTime();

        //cplex.exportModel("path../model2d.sav");
        //cplex.exportModel("path../model2d.lp");

        loss = cplex.getObjValue();

        paras.clear();
        for (int i = 0; i < (degree+1)*(degree+1); i++) {
            paras.push_back(cplex.getValue(vars[i]));
        }

        env.end();
    }

    QueryResult CountPrediction2D(vector<Tk> &d1_low, vector<Tk> &d2_low, vector<Tk> &d1_up, vector<Tk> &d2_up, vector<Tv> &results,
                                  bool refinement, string RealResultPath){
        results.clear();
        Tv pred1, pred2, pred3, pred4, result;
        std::vector<value> result_n1, result_n2, result_n3, result_n4;
        int count_refinement = 0;
        auto t0 = chrono::steady_clock::now();
        for(int i = 0; i < d1_low.size(); i ++){
            rtree.query(bgi::nearest(point(d1_low[i], d2_low[i]), 1), std::back_inserter(result_n1));
            rtree.query(bgi::nearest(point(d1_up[i], d2_up[i]), 1), std::back_inserter(result_n2));
            rtree.query(bgi::nearest(point(d1_low[i], d2_up[i]), 1), std::back_inserter(result_n3));
            rtree.query(bgi::nearest(point(d1_up[i], d2_low[i]), 1), std::back_inserter(result_n4));
            pred1 = models[result_n1.back().second](d1_low[i], d2_low[i]);
            pred2 = models[result_n2.back().second](d1_up[i], d2_up[i]);
            pred3 = models[result_n3.back().second](d1_low[i], d2_up[i]);
            pred4 = models[result_n4.back().second](d1_up[i], d2_low[i]);
            result = pred1 + pred2 - pred3 - pred4;
            if(refinement && !RelativeErrorCheck_1D_Count(result, t_abs, t_rel)){
                count_refinement++;
                // if this is too slow, use aggregate Rtree instead
                box query_box(point(d1_low[i], d2_low[i]), point(d1_up[i], d2_up[i]));
                auto r = boost::make_iterator_range(bgi::qbegin(data_rtree, bgi::intersects(query_box)), {});
                result = boost::distance(r); // https://stackoverflow.com/questions/29261563/boost-r-tree-counting-elements-satisfying-a-query
            }
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
        query_result.hit_count = d1_low.size() - count_refinement;
        query_result.model_amount = this->models.size();

        return query_result;
    }

    // each leaf layer model, i.e., the polynomial function
    struct Model{
        Tk key_L1, key_L2, key_U1, key_U2; // the bounding rectangle of this model region
        array<array<double,degree+1>, degree+1> paras;

        Model(Tk key_L1, Tk key_L2, Tk key_U1, Tk key_U2, array<array<double,degree+1>, degree+1> &pas):
            key_L1(key_L1), key_L2(key_L2), key_U1(key_U1), key_U2(key_U2){
            for(int i = 0; i <= degree; i++){
                for(int j = 0; j <= degree; j++){
                    paras[i][j] = pas[i][j];
                }
            }
        }

        // @param paras : the paras should have size (deg+1)*(deg+1)
        Model(Tk key_L1, Tk key_L2, Tk key_U1, Tk key_U2, vector<double> &pas):
            key_L1(key_L1), key_L2(key_L2), key_U1(key_U1), key_U2(key_U2){
            for(int i = 0; i <= degree; i++){
                for(int j = 0; j <= degree; j++){
                    paras[i][j] = pas[i*(degree+1)+j];
                }
            }
        }

        inline Tv operator()(const Tk &k1, const Tk &k2) const {
            auto key_diff_1 = k1 - key_L1;
            auto key_diff_2 = k2 - key_L2;
            array<Tk, degree+1> key_term_1, key_term_2;
            //vector<Tk> key_term_1(degree+1), key_term_2(degree+1);
            Tk current_key_term_1 = 1.0, current_key_term_2 = 1.0;
            key_term_1[0] = 1.0;
            key_term_2[0] = 1.0;
            for(int i = 1; i <= degree; i++){
                current_key_term_1 *= key_diff_1;
                current_key_term_2 *= key_diff_2;
                key_term_1[i] = current_key_term_1;
                key_term_2[i] = current_key_term_2;
            }
            // calculate prediction
            Tv pred = 0;
            for(int i = 0; i <= degree; i++){
                for(int j = 0; j <= degree; j++){
                    pred += paras[i][j] * key_term_1[i] * key_term_2[j];
                }
            }
            return pred;
        }
    };

    std::tuple<Tk, Tk, Tk, Tk> global_domain; // indicate the domain, L1, L2, U1, U2

    vector<Model> models;
    bgi::rtree<value, bgi::quadratic<32>> rtree; // Boost Rtree, used to index bottom layer models, quadratic<max elements, min elements>
    bgi::rtree<value, bgi::quadratic<32>> data_rtree; // full rtee that store every records
    //RTree<int, double, 2, float> rtree; // simple Rtree

    Tv t_abs; // absolute error threshold
    double t_rel; // relative error threshold
};

#endif //POLYFIT_CLION_POLYFIT2D_H
