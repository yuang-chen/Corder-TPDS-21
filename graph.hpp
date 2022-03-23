#pragma once


#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <vector>
#include <cmath>
#ifndef SORT_HEADER
#include "sort.hpp"
#include "global.hpp"
#endif



class Graph{
public:
    unsigned num_vertex;
    unsigned num_edges;
    std::vector<unsigned> row_index;
    std::vector<unsigned> col_index;
    std::vector<unsigned> out_degree;
    std::vector<float> attr;
#ifdef WEIGHTED
    std::vector<unsigned> edge_weight;
#endif
    Graph(unsigned num_vertex = 0, unsigned num_edges = 0):
            num_vertex(num_vertex), num_edges(num_edges){
                row_index.reserve(num_vertex + 1);
                col_index.reserve(num_edges);
                attr.reserve(num_vertex);
            }
    ~Graph(){}
    
    void computeOutDegree() {
        out_degree = std::vector<unsigned>(num_vertex, 0);
        uint64_t total = 0;
        float avr = num_edges / num_vertex;
        
        #pragma omp parallel for schedule(static, 1024 * 256) num_threads(params::num_threads)
        for(unsigned i = 0; i < num_vertex; i++) {
            out_degree[i] = row_index[i + 1] - row_index[i];
        }

    }
   
    void printGraph(bool all = false) {
        std::cout << "num vertex: " << num_vertex << std::endl;
        std::cout << "num edges: " << num_edges <<std::endl;
        #ifdef WEIGHTED
        std::cout << "weighted graph" <<std::endl;
        #endif
        if(all) {
            for(auto it = row_index.begin(); it != row_index.end(); ++it) 
                std::cout << *it << " ";
            std::cout << std::endl;
            for(auto it = col_index.begin(); it != col_index.end(); ++it) 
                std::cout << *it << " ";
            std::cout << std::endl;
            for(auto it = attr.begin(); it != attr.end(); ++it) 
                std::cout << *it << " ";
            std::cout << std::endl;
        }
    }

    void initAttribute()
    {
        attr = std::vector<float>(num_vertex, 1.0);
        #pragma omp parallel for schedule(static) num_threads(params::num_threads)
        for (unsigned int i=0; i<num_vertex; i++)
        {
            attr[i] = 1;
            if (out_degree[i] > 0)
                attr[i] = (1.0)/out_degree[i];  
        }
        return;
    }
};

