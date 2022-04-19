//
// Created by Cloud on 2020/9/3.
//

#ifndef POLYFIT_CLION_SWAB_H
#define POLYFIT_CLION_SWAB_H

#endif //POLYFIT_CLION_SWAB_H

#include <vector>
#include <deque>
using namespace std;

/**
 * the online segmentation algorithm from time-series database
 * Sliding Window and Bottom-up
 * @tparam Tk
 * @tparam Tv
 */
template <typename Tk, typename Tv>
class SWAB{
public:
    SWAB(Tv t_abs, int buffer_seg_num): t_abs(t_abs),  buffer_seg_num(buffer_seg_num) {}

    struct Segment{
        Tk key; // the first key of this segment
        double slope;
        double intercept;

        Segment(Tk key, double slope, double intercept): key(key), slope(slope), intercept(intercept) {}

        inline Tv operator()(const Tk &k) const {
            auto key_diff = k - key;
            auto pred = slope * key_diff + intercept;
            return pred;
        }

        friend inline bool operator<(const Segment &s, const Tk &k) { return s.key < k; }
        friend inline bool operator<(const Tk &k, const Segment &s) { return k < s.key; }
    };

    void build(const vector<Tk> &key_v, const vector<Tv> &value_v){
        int buffer_size = buffer_seg_num * t_abs * 2;
        int lower_bound = buffer_size / 2;
        int upper_bound = buffer_size * 2;
        deque<Tk> buffer(key_v.begin(), key_v.begin()+upper_bound); // read in the initial data points

        segments.clear();
        int start_index = upper_bound, end_index = upper_bound;
        while(end_index < key_v.size()){
            Segment s = BottomUp(buffer);
            segments.push_back(s);

            BestLine(key_v, start_index, end_index);
            // append start_index to end_index to buffer
            for(int i = start_index; i < end_index; i++){
                buffer.push_back(key_v[i]);
            }
            start_index = end_index;
        }

    }

    Segment BottomUp(deque<Tk> &buffer){

        // create initial finest segments, using n/2 segments
        vector<Segment> temp_segments;

        // calculate the merget cost
        vector<double> merge_cost(temp_segments.size());

        // find the lowest merge cost
        while(*min_element(temp_segments) < t_abs){
            int minElementIndex = min_element(temp_segments.begin(), temp_segments.end()) - temp_segments.begin();
            // merge the ith and i+1 segments


            temp_segments.erase(temp_segments.begin() + minElementIndex+1);

            // update the merge cost


        }
    }

    void BestLine(const vector<Tk> &key_v, int start_index, int &end_index){
        double error = 0;
        end_index = start_index;
        while(error < t_abs){
            error = MeasureApproxError(key_v, start_index, end_index);
        }
    }

    void MeasureApproxError(const vector<Tk> &key_v, int start_index, int end_index){

    }



    vector<Segment> segments;

    Tv t_abs; // the absolute error threshold, default 100
    int buffer_seg_num; // maximum number of segments can stored in a buffer, default 5
};