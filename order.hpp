#pragma once

#include <iostream>
#include <limits.h>
#include <parallel/algorithm>
#include <parallel/numeric>
#include <algorithm>
#include <random>

#include "graph.hpp"
#include "vec2d.hpp"
#include "global.hpp"
#include <boost/timer/timer.hpp>
using namespace boost::timer;

class Orderer {
   // typedef std::vector<std::vector<unsigned>> Vector2D;
    typedef std::pair<unsigned, unsigned> degree_id_pair;
    unsigned num_vertex;
    unsigned num_edges;
    unsigned num_partitions;
    unsigned num_levels;
    unsigned average_degree;
    std::vector<unsigned> levels;

    Graph* graph;
public:
    std::vector<unsigned> new_id;

    Orderer(Graph* g) {
        this->graph = g;
        num_vertex = graph->num_vertex;
        num_edges = graph->num_edges;
        num_partitions = params::num_partitions;
        new_id = std::vector<unsigned>(num_vertex, 0);
        average_degree = num_edges/num_vertex;

        // level: aver_degree/2, aver_degree, aver_degree*2
        num_levels = (unsigned) log2((float) num_partitions) + 2;

        for(int i = 0; i < num_levels; i++)
            levels.push_back(average_degree * pow(2, i - 1));

        levels.back() = UINT_MAX;
    }

    void fastRandom() {        
        std::vector<unsigned> index(num_vertex);
        #pragma omp parallel for
        for(unsigned i = 0; i < num_vertex; i++)
            index[i] = i; 

        auto rng = std::default_random_engine {};
        std::shuffle(std::begin(index), std::end(index), rng);

        #pragma omp parallel for
        for(unsigned i = 0; i < num_vertex; i++)
            new_id[index[i]] = i;  

    }

    void fastSort() {
        std::vector<unsigned>& out_degree = graph->out_degree;
        std::vector<degree_id_pair> degree_vs_id(num_vertex);

        #pragma omp parallel for
        for(unsigned i = 0; i < num_vertex; i++)
            degree_vs_id[i] = std::make_pair(out_degree[i], i);

         __gnu_parallel::sort(degree_vs_id.begin(), degree_vs_id.end(), std::greater<degree_id_pair>());
        
        std::cout << degree_vs_id.front().first << " " << degree_vs_id.front().second << '\n';

        #pragma omp parallel for
        for(unsigned i = 0; i < num_vertex; i++)
            new_id[degree_vs_id[i].second] = i;
    }

    // hubcluster
    void fastDBG(unsigned num_levels) {

        levels.clear();
        for(int i = 0; i < num_levels; i++)
            levels.push_back(average_degree * pow(2, i - 1));
        levels.back() = UINT_MAX;

        const auto& out_degree = graph->out_degree;

        std::vector<unsigned>segment[params::num_threads][num_levels];
        
        #pragma omp parallel for schedule(static) num_threads(params::num_threads)
        for(unsigned i = 0; i < num_vertex; i++) {
            for(unsigned j = 0; j < num_levels; j++) {
                if(out_degree[i] <= levels[j]) {
                    segment[omp_get_thread_num()][j].push_back(i);
                    break;
                }
            }
        }
        unsigned tmp = 0;
        unsigned seg_offset[params::num_threads][num_levels];
        for(int j = num_levels - 1; j >= 0; j--)
            for(unsigned t = 0; t < params::num_threads; t++) {
                seg_offset[t][j] = tmp;
                tmp += segment[t][j].size();
            }
        
        #pragma omp parallel for schedule(static) num_threads(params::num_threads)
        for(unsigned t = 0; t < params::num_threads; t++)
            for(int j = num_levels - 1; j >= 0; j--) {
                unsigned offset = seg_offset[t][j];
                const std::vector<unsigned>& curr_seg = segment[t][j];
                for(auto id: curr_seg)
                    new_id[id] = offset++;
            }
    }

    
    void fastFBC() {
        new_id.clear();
        new_id = std::move(std::vector<unsigned>(num_vertex, UINT_MAX));

        Vector2d<degree_id_pair> local_degree_vs_id(params::num_threads);
        unsigned slice = num_vertex / params::num_threads;
        std::vector<unsigned> start(params::num_threads); 
        std::vector<unsigned> end(params::num_threads);
        std::vector<unsigned> hub_count(params::num_threads);
        std::vector<unsigned> non_hub_count(params::num_threads);
        std::vector<unsigned> index(params::num_threads);
        std::vector<unsigned>& out_degree = graph->out_degree;
        unsigned sum_hub_count = 0;

        for(int t = 0; t < params::num_threads; t++) {
            start[t] = t * slice;
            end[t] = (t + 1) * slice;
            hub_count[t] = 0;
        }
        end[params::num_threads - 1] = num_vertex;

        #pragma omp parallel for schedule(static) num_threads(params::num_threads)
        for(unsigned t = 0; t < params::num_threads; t++) 
            for(unsigned i = start[t]; i < end[t]; i++)
                if(out_degree[i] > average_degree) 
                    local_degree_vs_id[t].push_back(std::make_pair(out_degree[i], i));

        std::vector<unsigned> hub_offset(params::num_threads + 1, 0);
        for(int t = 0; t < params::num_threads; t++) {
            hub_count[t] = local_degree_vs_id[t].size();
            sum_hub_count += hub_count[t];
            non_hub_count[t] = end[t] - start[t] - hub_count[t];
            hub_offset[t + 1] =  hub_offset[t] + hub_count[t];
        }        
        index[0] = sum_hub_count;
        for(int t = 1; t < params::num_threads; t++) 
            index[t] = index[t - 1] + non_hub_count[t - 1];

        std::vector<degree_id_pair> degree_vs_id(sum_hub_count);

        unsigned tmp = 0;
        #pragma omp parallel for schedule(static) num_threads(params::num_threads)
        for(int i = 0; i < params::num_threads; i++) {
            for(unsigned j = 0; j < hub_count[i]; j++)
                degree_vs_id[hub_offset[i]++] = local_degree_vs_id[i][j];
            std::vector<degree_id_pair>().swap(local_degree_vs_id[i]);
        }
        __gnu_parallel::sort(degree_vs_id.begin(), degree_vs_id.end(),
            std::greater<degree_id_pair>());

        #pragma omp parallel for
        for(unsigned i = 0; i < sum_hub_count; i++)
            new_id[degree_vs_id[i].second] = i;
        
        std::vector<degree_id_pair>().swap(degree_vs_id);

    #pragma omp parallel for schedule(static) num_threads(params::num_threads)
        for(int t = 0; t < params::num_threads; t++)
            for(unsigned i = start[t]; i < end[t]; i++)
                if(new_id[i] == UINT_MAX) 
                    new_id[i] = index[t]++;

    }
   
    void fastHC() {
        levels.clear();
        num_levels = 2;
        levels.push_back(average_degree);
        levels.push_back(UINT_MAX);

        const auto& out_degree = graph->out_degree;

        std::vector<unsigned>segment[params::num_threads][num_levels];
        
        #pragma omp parallel for schedule(static) num_threads(params::num_threads)
        for(unsigned i = 0; i < num_vertex; i++) {
            for(unsigned j = 0; j < num_levels; j++) {
                if(out_degree[i] <= levels[j]) {
                    segment[omp_get_thread_num()][j].push_back(i);
                    break;
                }
            }
        }
        unsigned tmp = 0;
        unsigned seg_offset[params::num_threads][num_levels];
        for(int j = num_levels - 1; j >= 0; j--)
            for(unsigned t = 0; t < params::num_threads; t++) {
                seg_offset[t][j] = tmp;
                tmp += segment[t][j].size();
            }
        
        #pragma omp parallel for schedule(static) num_threads(params::num_threads)
        for(unsigned t = 0; t < params::num_threads; t++)
            for(int j = num_levels - 1; j >= 0; j--) {
                unsigned offset = seg_offset[t][j];
                const std::vector<unsigned>& curr_seg = segment[t][j];
                for(auto id: curr_seg)
                    new_id[id] = offset++;
            }
   }

   void Corder() {
     
        unsigned max_threads = omp_get_max_threads();
        std::vector<unsigned> segment_large;
        segment_large.reserve(num_vertex);
        std::vector<unsigned> segment_small;
        segment_small.reserve(num_vertex/2);

        for(unsigned i = 0; i < num_vertex; i++)
            if(graph->out_degree[i] > 1 * average_degree)
                segment_large.push_back(i);
            else
                segment_small.push_back(i);

        unsigned num_large_per_seg = ceil((float) segment_large.size() / num_partitions);
        params::overflow_ceil = num_large_per_seg;

        unsigned num_small_per_seg = params::partition_size - num_large_per_seg;

        std::cout << "partition size: " << params::partition_size  << " num of large: " << num_large_per_seg << " num of small: " << num_small_per_seg << '\n';
        unsigned last_cls = num_partitions - 1;

        #pragma omp parallel for schedule(static) num_threads(max_threads)
        for(unsigned i = 0; i < last_cls; i++) {
            unsigned index = i * params::partition_size;
            for(unsigned j = 0; j < num_large_per_seg; j++) {
                new_id[segment_large[i * num_large_per_seg + j]] = index++;
            }
            for(unsigned j = 0; j < num_small_per_seg; j++)
                new_id[segment_small[i * num_small_per_seg + j]] = index++;
        }

        auto last_large = num_large_per_seg * last_cls;
        auto last_small = num_small_per_seg * last_cls;
        unsigned index = last_cls * params::partition_size;

        for(unsigned i = last_large; i < segment_large.size(); i++) {
            new_id[segment_large[i]] = index++;
        }
        for(unsigned i = last_small; i < segment_small.size(); i++) {
            new_id[segment_small[i]] = index++;
        }
   }
    
   void fastCorder() {
        unsigned max_threads = omp_get_max_threads();

        Vector2d<unsigned> large_segment(max_threads);
        Vector2d<unsigned> small_segment(max_threads);

        const auto average_degree = num_edges/num_vertex;
        
        // classifying hot/cold vertices
        // the static scheduler ensures the relative order 
        // (static, 1024) disrupts the relative order but achieves very good performance
        #pragma omp parallel for schedule(static) num_threads(max_threads) 
        for(unsigned i = 0; i < num_vertex; i++)
            if(graph->out_degree[i] > average_degree) 
                large_segment[omp_get_thread_num()].push_back(i);
            else
                small_segment[omp_get_thread_num()].push_back(i);

        std::vector<unsigned> large_offset(max_threads + 1, 0);
        std::vector<unsigned> small_offset(max_threads + 1, 0);

        large_offset[1] = large_segment[0].size();
        small_offset[1] = small_segment[0].size(); 
        for(unsigned i = 0; i < max_threads ; i++) {
            large_offset[i+1] = large_offset[i] + large_segment[i].size();
            small_offset[i+1] = small_offset[i] + small_segment[i].size();
        }

       unsigned total_large = large_offset[max_threads];
       unsigned total_small = small_offset[max_threads]; 
     
        unsigned num_large_per_seg = ceil((float) total_large  / num_partitions);
        unsigned num_small_per_seg = params::partition_size - num_large_per_seg;

        unsigned last_cls = num_partitions - 1;

      //  constructing partitions based on the classified hot/cold vertices
        #pragma omp parallel for schedule(static) num_threads(max_threads)
        for(unsigned i = 0; i < num_partitions; i++) {
            unsigned index = i * params::partition_size;
            unsigned num_large =  (i != num_partitions - 1) ? (i + 1) * num_large_per_seg: total_large;
            unsigned large_start_t = 0;
            unsigned large_end_t = 0;
            unsigned large_start_v = 0;
            unsigned large_end_v = 0;
            unsigned large_per_seg = (i != num_partitions - 1) ? num_large_per_seg: total_large - i * num_large_per_seg;

            unsigned num_small =  (i != num_partitions - 1) ? (i + 1) * num_small_per_seg: total_small;
            unsigned small_start_t = 0;
            unsigned small_end_t = 0;
            unsigned small_start_v = 0;
            unsigned small_end_v = 0;
            unsigned small_per_seg = (i != num_partitions - 1) ? num_small_per_seg: total_small - i * num_small_per_seg;
            //HOT find the starting segment and starting vertex
            for(unsigned t = 0; t < max_threads; t++) {
                if(large_offset[t+1] > num_large - large_per_seg) {
                    large_start_t = t;
                    large_start_v = num_large - large_per_seg - large_offset[t];
                    break;
                }
            }
            //HOT find the ending segment and ending vertex
            for(unsigned t = large_start_t; t < max_threads; t++) {
                if(large_offset[t+1] >= num_large) {
                    large_end_t = t;
                    large_end_v =  num_large - large_offset[t] - 1;
                    break;
                }
            }

            //COLD find the starting segment and starting vertex
            for(unsigned t = 0; t < max_threads; t++) {
                if(small_offset[t+1] > num_small - small_per_seg) {
                    small_start_t = t;
                    small_start_v = num_small - small_per_seg - small_offset[t];
                    break;
                }
            }
            //COLD find the ending segment and ending vertex
           for(unsigned t = small_start_t; t < max_threads; t++) {
                if(small_offset[t+1] >= num_small) {
                    small_end_t = t;
                    small_end_v =  num_small - small_offset[t] - 1;
                    break;
                }
            }
     // HOT move the vertices form hot segment(s) to a partition
            if(large_start_t == large_end_t) {
                for(unsigned j = large_start_v; j <= large_end_v; j++) {
                    new_id[large_segment[large_start_t][j]] = index++;
                }
            } else {
                for(unsigned t = large_start_t; t < large_end_t; t++) {
                    if(t!=large_start_t)
                        large_start_v = 0;
                    for(unsigned j = large_start_v; j < large_segment[t].size(); j++) {
                        new_id[large_segment[t][j]] = index++;
                    }
                }
                for(unsigned j = 0; j <= large_end_v; j++) {
                    new_id[large_segment[large_end_t][j]] = index++;
                }
            }
    // COLD move the vertices form cold segment(s) to a partition
            if(small_start_t == small_end_t) {
                for(unsigned j = small_start_v; j <= small_end_v; j++) {
                    new_id[small_segment[small_start_t][j]] = index++;
                }
            } else {
                for(unsigned t = small_start_t; t < small_end_t; t++) {
                    if(t!=small_start_t)
                        small_start_v = 0;
                    for(unsigned j = small_start_v; j < small_segment[t].size(); j++) {
                        new_id[small_segment[t][j]] = index++;
                    }
                }
                for(unsigned j = 0; j <= small_end_v; j++) {
                    new_id[small_segment[small_end_t][j]] = index++;
                }
            }
        }
   }

    void getNewGraph() {
        cpu_timer timer;
        float time = 0.0;
        unsigned max_threads = omp_get_max_threads();
      //  Graph new_graph(num_vertex, num_edges);

        std::vector<unsigned> new_degree(num_vertex, 0);
        timer.start();
        //Assign the outdegree to new id
        #pragma omp parallel for schedule(static) num_threads(max_threads)
        for(unsigned i = 0; i < num_vertex; i++) 
            new_degree[new_id[i]] = graph->out_degree[i];
        float tm = timer.elapsed().wall/(1e9); 

        // Build new row_index array
        std::vector<unsigned> new_row(num_vertex + 1, 0);
        __gnu_parallel::partial_sum(new_degree.begin(), new_degree.end(), new_row.begin() + 1);
        std::vector<unsigned> new_col(num_edges, 0);
        tm = timer.elapsed().wall/(1e9); 

        #ifdef WEIGHTED
         std::vector<unsigned> new_wei(num_edges, 0);
        #endif
        //Build new col_index array
        #pragma omp parallel for schedule(dynamic) num_threads(max_threads)
        for(unsigned i = 0; i < num_vertex; i++) {
            unsigned count = 0;
            for(unsigned j = graph->row_index[i]; j < graph->row_index[i + 1]; j++) {
                new_col[new_row[new_id[i]] + count] = new_id[graph->col_index[j]];
                #ifdef WEIGHTED
                new_wei[new_row[new_id[i]] + count] = graph->edge_weight[j];
                #endif
                count++;
            }
         }
        tm = timer.elapsed().wall/(1e9); 

        this->graph->out_degree.swap(new_degree);
        this->graph->row_index.swap(new_row);
        this->graph->col_index.swap(new_col);

        #ifdef WEIGHTED
        new_graph.edge_weight.swap(new_wei);
        #endif
    }

    std::vector<unsigned> getLevels() {
        return levels;
    }

    void reorder(Algo algo){
        switch(algo) {
            case Algo::original: 
                std::cout << "original order is maintained" << '\n';
                break;
            case Algo::randm:
                std::cout << "reordering method: random" << '\n';
                fastRandom();
                break;
            case Algo::sort:
                std::cout << "reordering method: sort" << '\n';
                fastSort();
                break;
            case Algo::fbc:
                std::cout << "reordering method: fbc" << '\n';
                fastFBC();
                break;
            case Algo::dbg:
                std::cout << "reordering method: dbg" << '\n';
                fastDBG(8);
                break;
            case Algo::hc:
                std::cout << "reordering method: hc" << '\n';
                fastHC();
                break;
            case Algo::corder:
                std::cout << "reordering method: corder" << '\n';
                Corder();
                break;
            case Algo::fastCorder:
                std::cout << "reordering method: parallized corder" << '\n';
                fastCorder();
                break;
            default:
                std::cout << "choose a correct algorithm!" << '\n';
        }
    }
};