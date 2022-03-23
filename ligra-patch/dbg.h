/*
 * @author Priyank Faldu <Priyank.Faldu@ed.ac.uk> <http://faldupriyank.com>
 *
 * Copyright 2019 The University of Edinburgh
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <vector>
#include <parallel/algorithm>
#include "ligra.h"
#include "pvector.h"
#include "math.h"

typedef std::pair<uintT, uintE> degree_nodeid_t;

template <typename T>
using Vector2d = typename std::vector<std::vector<T>>;

extern int rand_gran;
extern string map_file;

enum ReorderingAlgo {
    ORIGINAL = 0,
    Corder = 1,
    HubCluster = 2,
    HubSort = 3,
    DBG = 4,
    Sort = 5,
    Random = 6, 
    HubSortDBG = 8,
    HubClusterDBG = 7,
    MAP = 10, //TODO: MAP file format
    MAP_ORIGINAL = MAP,
    MAP_Sort = 12,
    MAP_HubSort = 13,
    MAP_HubCluster = 14,
    MAP_DBG = 15,
    MAP_HubSortDBG = 16,
    MAP_HubClusterDBG = 17
    //Random = 18
};

const string ReorderingAlgoStr(ReorderingAlgo type) {
    switch(type) {
        case Corder:
            return "Corder";
        case HubSort:
            return "HubSort";
        case DBG:
            return "DBG";
        case HubClusterDBG:
            return "HubClusterDBG";
        case HubSortDBG:
            return "HubSortDBG";
        case HubCluster:
            return "HubCluster";
        case Random:
            return "Random-" + to_string(rand_gran);
        case ORIGINAL:
            return "Original";
        case Sort:
            return "Sort";
        case MAP:
            return "MAP";
        case MAP_HubSort:
            return "MAP_HubSort";
        case MAP_HubCluster:
            return "MAP_HubCluster";
        case MAP_DBG:
            return "MAP_DBG";
        case MAP_Sort:
            return "MAP_Sort";
        case MAP_HubSortDBG:
            return "MAP_HubSortDBG";
        case MAP_HubClusterDBG:
            return "MAP_HubClusterDBG";
        default:
            std::cout << "Unknown PreprocessTypeStr type: " << type << std::endl;
            abort();
    }
    assert(0);
    return "";
}

template <class vertex>
int verify_mapping(const graph<vertex>& GA, const pvector<uintE>& new_ids, long numVertices) {
    uintE* hist = newA(uintE, numVertices);
    {parallel_for(long i = 0 ; i < numVertices ; i++ ) {
            hist[i] = new_ids[i];
    }}

	__gnu_parallel::sort(&hist[0], &hist[numVertices]); 

    uintE count = 0;
    {parallel_for(long i = 0 ; i < numVertices ; i++ ) {
        if ( hist[i] != i) {
            __sync_fetch_and_add(&count, 1);
        }
    }}
    if ( count != 0 ){
        std::cout << "Num of vertices did not match: " << count << std::endl;
        std::cout << "Mapping is invalid.!" << std::endl;
        abort();
    } else {
        std::cout << "Mapping is valid.!" << std::endl;
    }
    free(hist);
}

template <class vertex>
void generateRandomMapping(const graph<vertex>& GA, bool isSym, bool useOutdeg, pvector<uintE>& new_ids, bool isPageRank, bool isDenseWrite) {
    Timer t; 
    t.Start();

    auto num_vertices = GA.n;
    auto num_edges = GA.m;

    uintE granularity = rand_gran;
    uintE slice = (num_vertices - granularity + 1) / granularity;
    uintE artificial_num_vertices = slice * granularity;
    assert(artificial_num_vertices <= num_vertices);
    pvector<uintE> slice_index;
    slice_index.resize(slice);

    {parallel_for(long i=0;i<slice;i++) {
        slice_index[i] = i;
    }}

    std::random_shuffle(slice_index.begin(), slice_index.end());

    {parallel_for(long i = 0 ; i < slice ; i++) {
        long new_index = slice_index[i] * granularity;
        for(long j = 0 ; j < granularity ; j++ ) {
            long v = (i * granularity) + j;
            if ( v < artificial_num_vertices) {
                new_ids[v] = new_index + j;
            }
        }
    }}
    for (long i = artificial_num_vertices ; i < num_vertices ; i++ ) {
        new_ids[i] = i;
    }
    slice_index.clear();

    t.Stop();
    t.PrintTime("Random Map Time", t.Seconds());
}

template <class vertex>
void generateHubSortDBGMapping(const graph<vertex>& GA, bool isSym, bool useOutdeg, pvector<uintE>& new_ids, bool isPageRank, bool isDenseWrite) {
    Timer t; 
    t.Start();

    auto numVertices = GA.n;
    auto numEdges    = GA.m;
    vertex *origG    = GA.V;
    uintT avgDegree  = numEdges / numVertices;
    uintT hubCount {0};

    const int num_threads = omp_get_max_threads();
    pvector<degree_nodeid_t> local_degree_id_pairs[num_threads];
    uintE slice = numVertices / num_threads;
    uintE start[num_threads];
    uintE end[num_threads];
    uintE hub_count[num_threads];
    uintE non_hub_count[num_threads];
    uintE new_index[num_threads];
    for ( int t = 0 ; t < num_threads ; t++ ) {
        start[t] = t * slice;
        end[t] = (t+1) * slice;
        hub_count[t] = 0;
    }
    end[num_threads-1] = numVertices;

#pragma omp parallel for schedule(static) num_threads(num_threads)
    for ( uintE t = 0 ; t < num_threads ; t++ ) {
        for (uintE v = start[t]; v < end[t]; ++v) {
            vertex vtx = origG[v];
            if (useOutdeg) {
                if (vtx.getOutDegree() > avgDegree) {
                    local_degree_id_pairs[t].push_back(std::make_pair(vtx.getOutDegree(), v));
                }
            } else {
                if (vtx.getInDegree() > avgDegree) {
                    local_degree_id_pairs[t].push_back(std::make_pair(vtx.getInDegree(), v));
                }
            }
        }
    }
    for ( int t = 0 ; t < num_threads ; t++ ) {
        hub_count[t] = local_degree_id_pairs[t].size();
        hubCount += hub_count[t];
        non_hub_count[t] = end[t] - start[t] - hub_count[t];
    }
    new_index[0] = hubCount;
    for ( int t = 1 ; t < num_threads ; t++ ) {
        new_index[t] = new_index[t-1] + non_hub_count[t-1];
    }
    pvector<degree_nodeid_t> degree_id_pairs(hubCount);

    long k = 0;
    for ( int i = 0 ; i < num_threads ; i++ ) {
        for ( long j = 0 ; j < local_degree_id_pairs[i].size() ; j++ ) {
            degree_id_pairs[k++] = local_degree_id_pairs[i][j];
        }
        local_degree_id_pairs[i].clear();
    }
    assert(degree_id_pairs.size() == hubCount);
    assert(k == hubCount);

    __gnu_parallel::sort(degree_id_pairs.begin(), degree_id_pairs.end(),
            std::greater<degree_nodeid_t>());

#pragma omp parallel for
    for (uintE n = 0; n < hubCount; ++n) {
        new_ids[degree_id_pairs[n].second] = n;
    }
    pvector<degree_nodeid_t>().swap(degree_id_pairs);

#pragma omp parallel for schedule(static) num_threads(num_threads)
    for ( uintE t = 0 ; t < num_threads ; t++ ) {
        for (uintE v = start[t]; v < end[t]; ++v) {
            if ( new_ids[v] == UINT_E_MAX ) {
                new_ids[v] = new_index[t]++;
            }
        }
    }

    t.Stop();
    t.PrintTime("HubSortDBG Map Time", t.Seconds());
}

template <class vertex>
void generateHubClusterDBGMapping(const graph<vertex>& GA, bool isSym, bool useOutdeg, pvector<uintE>& new_ids, bool isPageRank, bool isDenseWrite) {
    Timer t;
    t.Start();

    auto num_vertices = GA.n;
    auto num_edges = GA.m;
    vertex *origG    = GA.V;

    uint32_t avg_vertex = num_edges / num_vertices;

    const int num_buckets = 2;
    avg_vertex = avg_vertex;
    uint32_t bucket_threshold[] = {avg_vertex, static_cast<uint32_t>(-1)};
    
    vector<uint32_t> bucket_vertices[num_buckets];
    const int num_threads = omp_get_max_threads();
    vector<uint32_t> local_buckets[num_threads][num_buckets];

    if ( useOutdeg ) {
    // This loop relies on a static scheduling
        #pragma omp parallel for schedule(static)
        for ( uint64_t i = 0 ; i < num_vertices ; i++ ) {
            for ( unsigned int j = 0 ; j < num_buckets ; j++ ) {
                const uintE& count = origG[i].getOutDegree();
                if ( count <= bucket_threshold[j] ) {
                    local_buckets[omp_get_thread_num()][j].push_back(i);
                    break;
                }
            }
        }
    } else {
        #pragma omp parallel for schedule(static)
        for ( uint64_t i = 0 ; i < num_vertices ; i++ ) {
            for ( unsigned int j = 0 ; j < num_buckets ; j++ ) {
                const uintE& count = origG[i].getInDegree();
                if ( count <= bucket_threshold[j] ) {
                    local_buckets[omp_get_thread_num()][j].push_back(i);
                    break;
                }
            }
        }
    }

    int temp_k = 0;
    uint32_t start_k[num_threads][num_buckets];
    for ( int32_t j = num_buckets-1 ; j >= 0 ; j-- ) {
        for ( int t = 0 ; t < num_threads ; t++ ) {
            start_k[t][j] = temp_k;
            temp_k += local_buckets[t][j].size();
        }
    }

    #pragma omp parallel for schedule(static)
    for ( int t = 0 ; t < num_threads; t++ ) {
        for ( int32_t j = num_buckets-1 ; j >= 0 ; j-- ) {
            const vector<uint32_t>& current_bucket = local_buckets[t][j];
            int k = start_k[t][j];
            const size_t& size = current_bucket.size();
            for ( uint32_t i = 0 ; i < size ; i++ ) {
                new_ids[ current_bucket[i] ] = k++;
            }
        }
    }

    for ( uint64_t i = 0 ; i < num_threads ; i++ ) {
        for ( unsigned int j = 0 ; j < num_buckets ; j++ ) {
            local_buckets[i][j].clear();
        }
    }

    t.Stop();
    t.PrintTime("HubClusterDBG Map Time", t.Seconds());
}

template <class vertex>
void generateDBGMapping(const graph<vertex>& GA, bool isSym, bool useOutdeg, pvector<uintE>& new_ids, bool isPageRank, bool isDenseWrite) {
    Timer t; 
    t.Start();

    auto num_vertices = GA.n;
    auto num_edges = GA.m;
    vertex *origG    = GA.V;


    uint32_t avg_vertex = num_edges / num_vertices;
    const uint32_t& av = avg_vertex;
    
    uint32_t bucket_threshold[] = {av/2, av, av*2, av*4, av*8, av*16, av*32, av*64, av*128, av*256, av*512, static_cast<uint32_t>(-1)};
    int num_buckets = 8;
    if ( num_buckets > 11 ) {
        // if you really want to increase the bucket count, add more thresholds to the bucket_threshold above.
            std::cout << "Unsupported bucket size: " << num_buckets << std::endl;
            assert(0);
    }
    bucket_threshold[num_buckets-1] = static_cast<uint32_t>(-1);

    vector<uint32_t> bucket_vertices[num_buckets];
    const int num_threads = omp_get_max_threads();
    vector<uint32_t> local_buckets[num_threads][num_buckets];

    if ( useOutdeg ) {
    // This loop relies on a static scheduling
        #pragma omp parallel for schedule(static)
        for ( uint64_t i = 0 ; i < num_vertices ; i++ ) {
            for ( unsigned int j = 0 ; j < num_buckets ; j++ ) {
                const uintE& count = origG[i].getOutDegree();
                if ( count <= bucket_threshold[j] ) {
                    local_buckets[omp_get_thread_num()][j].push_back(i);
                    break;
                }
            }
        }
    } else {
        #pragma omp parallel for schedule(static)
        for ( uint64_t i = 0 ; i < num_vertices ; i++ ) {
            for ( unsigned int j = 0 ; j < num_buckets ; j++ ) {
                const uintE& count = origG[i].getInDegree();
                if ( count <= bucket_threshold[j] ) {
                    local_buckets[omp_get_thread_num()][j].push_back(i);
                    break;
                }
            }
        }
    }

    int temp_k = 0;
    uint32_t start_k[num_threads][num_buckets];
    
    for ( int32_t j = num_buckets-1 ; j >= 0 ; j-- ) {
        for ( int t = 0 ; t < num_threads ; t++ ) {
            start_k[t][j] = temp_k;
            temp_k += local_buckets[t][j].size();
        }
    }

    #pragma omp parallel for schedule(static)
    for ( int t = 0 ; t < num_threads; t++ ) {
        for ( int32_t j = num_buckets-1 ; j >= 0 ; j-- ) {
            const vector<uint32_t>& current_bucket = local_buckets[t][j];
            int k = start_k[t][j];
            const size_t& size = current_bucket.size();
            for ( uint32_t i = 0 ; i < size ; i++ ) {
                new_ids[ current_bucket[i] ] = k++;
            }
        }
    }

    for ( uint64_t i = 0 ; i < num_threads ; i++ ) {
        for ( unsigned int j = 0 ; j < num_buckets ; j++ ) {
            local_buckets[i][j].clear();
        }
    }

    t.Stop();
    t.PrintTime("DBG Map Time", t.Seconds());
}

template <class vertex>
void generateSortMapping(const graph<vertex>& GA, bool isSym, bool useOutdeg, pvector<uintE>& new_ids, bool isPageRank, bool isDenseWrite) {
    Timer t; 
    t.Start();

    auto numVertices = GA.n;
    auto numEdges    = GA.m;
    vertex *origG    = GA.V;
    pvector<degree_nodeid_t> degree_id_pairs(numVertices);

    if (useOutdeg) {
#pragma omp parallel for
        for (uintE v = 0; v < numVertices; ++v) {
            degree_id_pairs[v] = std::make_pair(origG[v].getOutDegree(), v);
        }
    } else {
#pragma omp parallel for
        for (uintE v = 0; v < numVertices; ++v) {
            degree_id_pairs[v] = std::make_pair(origG[v].getInDegree(), v);
        }
    }

    __gnu_parallel::sort(degree_id_pairs.begin(), degree_id_pairs.end(),
            std::greater<degree_nodeid_t>());

#pragma omp parallel for
    for (uintE n = 0; n < numVertices; ++n) {
        new_ids[degree_id_pairs[n].second] = n;
    }

    pvector<degree_nodeid_t>().swap(degree_id_pairs);

    t.Stop();
    t.PrintTime("Sort Map Time", t.Seconds());
}

template <class vertex>
void generateCorder(const graph<vertex>& GA, bool isSym, bool useOutdeg, pvector<uintE>& new_id, bool isPageRank, bool isDenseWrite) {
     Timer t; 
    t.Start();

    auto num_vertex = GA.n;
    auto num_edges = GA.m;
    vertex *origG    = GA.V;


    uint32_t average_degree = num_edges / num_vertex;
    
    const int max_threads = omp_get_max_threads();

    Vector2d<unsigned> large_segment(max_threads);
    Vector2d<unsigned> small_segment(max_threads);

    #pragma omp parallel for schedule(static, 1024) num_threads(max_threads) 
    for(unsigned i = 0; i < num_vertex; i++)
        if(origG[i].getOutDegree() > average_degree) 
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

    unsigned cluster_size = 1024 * 1024 / sizeof(float);
    unsigned num_clusters = (num_vertex-1)/ cluster_size + 1;
    unsigned num_large_per_seg = ceil((float) total_large  / num_clusters);
    unsigned num_small_per_seg = cluster_size - num_large_per_seg;

    // Parallelize constructing partitions based on the classified hot/cold vertices
    #pragma omp parallel for schedule(static) num_threads(max_threads)
    for(unsigned i = 0; i < num_clusters; i++) {
        unsigned index = i * cluster_size;
        unsigned num_large =  (i != num_clusters - 1) ? (i + 1) * num_large_per_seg: total_large;
        unsigned large_start_t = 0;
        unsigned large_end_t = 0;
        unsigned large_start_v = 0;
        unsigned large_end_v = 0;
        unsigned large_per_seg = (i != num_clusters - 1) ? num_large_per_seg: total_large - i * num_large_per_seg;

        unsigned num_small =  (i != num_clusters - 1) ? (i + 1) * num_small_per_seg: total_small;
        unsigned small_start_t = 0;
        unsigned small_end_t = 0;
        unsigned small_start_v = 0;
        unsigned small_end_v = 0;
        unsigned small_per_seg = (i != num_clusters - 1) ? num_small_per_seg: total_small - i * num_small_per_seg;
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
    t.Stop();
    t.PrintTime("Corder Map Time", t.Seconds());
}


//Supported format is text file
//line1: V (number of vertices)
//line2: E (number of edges)
//line3+i: a pair of IDs i and i' separated by space. (where i is from 0 to V-1)
//
//Example map file with 3 vertices and 4 edges may look like below.
//3
//4
//0 2
//1 0
//2 1
template <class vertex>
void LoadMappingFromFile(const graph<vertex>& GA, bool isSym, bool useOutdeg, pvector<uintE>& new_ids, bool isPageRank, bool isDenseWrite) {
    Timer t; 
    t.Start();

    auto num_vertex = GA.n;
    auto num_edges = GA.m;

    ifstream ifs(map_file.c_str(), std::ifstream::in);
    if (!ifs.good()) {
        cout << "File " << map_file << " does not exist!" << endl;
        exit(-1);
    }
    unsigned long int num_vertex_1, num_edges_1;
    ifs >> num_vertex_1;
    ifs >> num_edges_1;
    cout << " num_vertex: " << num_vertex_1 << " num_edges: "  << num_edges_1 << endl;
    cout << " num_vertex: " << num_vertex << " num_edges: "  << num_edges << endl;
    if ( num_vertex != num_vertex_1 ) {
        cout << "Mismatch: " << num_vertex << " " << num_vertex_1  << endl;
        exit (-1);
    }
    if ( num_vertex != new_ids.size() ) {
        cout << "Mismatch: " << num_vertex << " " << new_ids.size() << endl;
        exit (-1);
    }
    if ( num_edges != num_edges_1 ) {
        cout << "Warning! Potential mismatch: " << num_edges << " " << num_edges_1 << endl;
    }
    char c;
    unsigned long int st, v;
    bool tab = true;
    if ( tab ) {
        for ( unsigned int i = 0 ; i < num_vertex ; i++ ) {
            ifs >> st >> v;
            new_ids[st] = v;
        }
    } else {
        for ( unsigned int i = 0 ; i < num_vertex ; i++ ) {
            ifs >> c >> st >> c >> v >> c;
            new_ids[st] = v;
        }
    }
    ifs.close();

    t.Stop();
    t.PrintTime("Load Map Time", t.Seconds());
}

template <class vertex>
void generateMapping(const graph<vertex>& GA, ReorderingAlgo reordering_algo, bool isSym, bool useOutdeg, pvector<uintE>& new_ids, bool isPageRank, bool isDenseWrite) {
    switch(reordering_algo) {
        case Corder:
            generateCorder(GA, isSym, useOutdeg, new_ids, isPageRank, isDenseWrite);
            break;
        case HubSort:
            generateHubSortMapping(GA, isSym, useOutdeg, new_ids, isPageRank, isDenseWrite);
            break;
        case Sort:
            generateSortMapping(GA, isSym, useOutdeg, new_ids, isPageRank, isDenseWrite);
            break;
        case DBG:
            generateDBGMapping(GA, isSym, useOutdeg, new_ids, isPageRank, isDenseWrite);
            break;
        case HubSortDBG:
            generateHubSortDBGMapping(GA, isSym, useOutdeg, new_ids, isPageRank, isDenseWrite);
            break;
        case HubClusterDBG:
            generateHubClusterDBGMapping(GA, isSym, useOutdeg, new_ids, isPageRank, isDenseWrite);
            break;
        case HubCluster:
            generateHubClusterMapping(GA, isSym, useOutdeg, new_ids, isPageRank, isDenseWrite);
            break;
        case Random:
            generateRandomMapping(GA, isSym, useOutdeg, new_ids, isPageRank, isDenseWrite);
            break;
        case MAP:
            LoadMappingFromFile(GA, isSym, useOutdeg, new_ids, isPageRank, isDenseWrite);
            break;
        case ORIGINAL:
            std::cout << "Should not be here!" << std::endl;
            abort();
            return;
            break;
        default:
            std::cout << "Unknown generateMapping type: " << reordering_algo << std::endl;
            abort();
    }
#ifdef _DEBUG
    verify_mapping(GA, new_ids, GA.n);
    exit(-1);
#endif
}
