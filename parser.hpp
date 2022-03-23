#pragma once


#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "graph.hpp"
#include <algorithm>


bool parseGraph(std::string filename, Graph& graph) {
    std::ifstream csr_file;
    csr_file.open(filename, std::ios::binary);
    if(!csr_file.is_open()) {
        std::cout << "cannot open csr file!" << std::endl;
        return false;
    }

    csr_file.read(reinterpret_cast<char*>(&graph.num_vertex), sizeof(unsigned));
    csr_file.read(reinterpret_cast<char*>(&graph.num_edges), sizeof(unsigned));

    std::vector<unsigned> local_row(graph.num_vertex);
    std::vector<unsigned> local_col(graph.num_edges);

    csr_file.read(reinterpret_cast<char*>(local_row.data()), graph.num_vertex * sizeof(unsigned));
    csr_file.read(reinterpret_cast<char*>(local_col.data()), graph.num_edges * sizeof(unsigned));

    local_row.push_back(graph.num_edges);
    graph.row_index = std::move(local_row);
    graph.col_index = std::move(local_col);

#ifdef WEIGHTED
    std::vector<unsigned> local_wei(graph.num_edges);
    csr_file.read(reinterpret_cast<char*>(local_wei.data()), graph.num_edges * sizeof(unsigned));
    graph.edge_weight = std::move(local_wei);
#endif
    csr_file.close();
    return true;
};

bool writeEdgelist(std::string filename, Graph& graph) {
    std::ofstream output(filename);
    if(!output.is_open()) {
        std::cout << "cannot open txt file!" << std::endl;
        return false;
    }
    for(unsigned i = 0; i < graph.num_vertex; i++) {
        for(unsigned j = graph.row_index[i]; j < graph.row_index[i+1]; j++) {
            output << i << " " << graph.col_index[j] << '\n';
        }
    }
        output.close();
}




bool writeGraph(std::string filename, Graph& graph) {
    std::ofstream  output_file(filename, std::ios::out | std::ios::binary);
    if(!output_file.is_open()) {
        std::cout << "cannot open output file!" << '\n';
        return false;
    }
    output_file.write(reinterpret_cast<char*>(&graph.num_vertex), sizeof(unsigned));
    output_file.write(reinterpret_cast<char*>(&graph.num_edges), sizeof(unsigned));

    output_file.write(reinterpret_cast<char*>(graph.row_index.data()), graph.num_vertex * sizeof(unsigned));
    output_file.write(reinterpret_cast<char*>(graph.col_index.data()), graph.num_edges * sizeof(unsigned));

#ifdef WEIGHTED
    output_file.write(reinterpret_cast<char*>(graph.edge_weight.data()), graph.num_edges * sizeof(unsigned));
#endif

    output_file.close();

    std::cout << filename << " is writen" << std::endl;
    return true;
};


void calculate_skew_locality(Graph& graph) {
        unsigned num_clusters = params::num_partitions;
        std::vector<unsigned> sum_degree(num_clusters - 1, 0);
        unsigned tmp = 0;
        #pragma omp parallel for schedule(static) num_threads(params::num_threads)
        for(unsigned i = 0; i < num_clusters - 1; i++) {
            for(unsigned j = 0; j < params::partition_size; j++) {
                sum_degree[i] += graph.out_degree[i * params::partition_size + j];
            }
        }
        std::sort(sum_degree.begin(), sum_degree.end());

        std::vector<unsigned> count;
        std::vector<float> locality_skew;
        std::vector<std::pair<unsigned, unsigned>> hotness;

        count.push_back(ceil(params::num_partitions * 0.01));
        count.push_back(ceil(params::num_partitions * 0.1));
        count.push_back(ceil(params::num_partitions * 0.2));
        count.push_back(ceil(params::num_partitions * 0.3));
        count.push_back(ceil(params::num_partitions * 0.4));
        count.push_back(ceil(params::num_partitions * 0.5));
        for(auto c : count) {
            float first = 0, second = 0;
            for(auto it = sum_degree.begin(); it != sum_degree.begin() + c; it++)
                first += *it;
            for(auto it = sum_degree.end() - c; it != sum_degree.end(); it++)
                second += *it;

            locality_skew.push_back(second/first);
            hotness.emplace_back(first, second);
        }

        std::cout << " 1%, 10%, 20%, 30%, 40%, 50% Locality Skew: " << " ";
        for(auto lk: locality_skew)
            std::cout << lk << " ";
        std::cout << std::endl;
        std::cout << " 1%, 10%, 20%, 30%, 40%, 50% Hotness Pairs: " << " ";
        for(auto h: hotness)
            std::cout << h.first << " " << h.second << " / ";
        std::cout << std::endl;
    }

void isOverflow(Graph& graph) {
    std::vector<unsigned> segment_large;
    segment_large.reserve(graph.num_vertex);
    std::vector<unsigned> segment_small;
    segment_small.reserve(graph.num_vertex/2);
    unsigned average_degree = graph.num_edges / graph.num_vertex;
    for(unsigned i = 0; i < graph.num_vertex; i++)
    if(graph.out_degree[i] > 1 * average_degree)
        segment_large.push_back(i);
    else
        segment_small.push_back(i);

    unsigned num_large_per_seg = ceil((float) segment_large.size() / params::num_partitions);
    params::overflow_ceil = num_large_per_seg;
    
    std::vector<unsigned> count(params::num_partitions - 1, 0);
    for(unsigned i = 0; i < params::num_partitions - 1; i++) {
        unsigned offset = i * params::partition_size;
        for(unsigned j = 0; j < params::partition_size; j++) {
            if(graph.out_degree[offset + j] > average_degree)
                count[i]++;
        }
    }

    for(unsigned i = 0; i < params::num_partitions - 1; i++) {
        int extra = count[i] - params::overflow_ceil;
        if(extra > 0) {
            std::cout << "subgraph " << i << " is overflowed by " << extra << " vertices." << std::endl;
        }
    }
}