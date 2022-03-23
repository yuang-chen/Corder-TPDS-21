/**
 *
 * This code implements work optimized propagation blocking with
 * transposed bin graph to reduce cache misses in scatter
 */
#include <pthread.h>
#include <time.h>
#include "graph.hpp"
#include "parser.hpp"
#include <boost/timer/timer.hpp>
#include <boost/program_options.hpp>
#include "vec2d.hpp"
#include "global.hpp"
#include "order.hpp"

using namespace std;
using namespace boost::timer;
using namespace boost::program_options;

//////////////////////////////////////////
//main function
//////////////////////////////////////////
int main(int argc, char** argv)
{

    options_description desc{"Options"};
    desc.add_options() 
        ("data,d", value<string>()->default_value(""), "input data path")
        ("output,o", value<string>()->default_value(""), "output file path")
        ("size,s", value<int>()->default_value(1024), "partition size")
        ("algorithm,a", value<int>()->default_value(0), "reordering algorithm")
        ("thread,t", value<int>()->default_value(20), "threads");


    variables_map vm;
    try
    {
        store(parse_command_line(argc, argv, desc), vm);
        notify(vm);
    } catch (error& e) {
        cerr << "ERROR: " << e.what() << '\n' << '\n' << desc << '\n';
        return 1;
    }

    string data_file = vm["data"].as<string>();
    string write_file = vm["output"].as<string>();
    int size = vm["size"].as<int>();
    int reorder = vm["algorithm"].as<int>();
    params::num_threads = vm["thread"].as<int>();

    if (data_file.empty()) {
        cout << desc << '\n';
        exit(1);
    } 

    Algo algo = static_cast<Algo>(reorder);
    params::partition_size = size * 1024 / sizeof(float);
     
    // graph object
    Graph graph;

    /**************************************************
     Compute the preprocessing time
     *************************************************/
    cpu_timer timer;
    cpu_times times;

    if(parseGraph(data_file, graph)) {
        times = timer.elapsed();
        cout << times.wall/(1e9) << "s: parse done for " << data_file << '\n';
        graph.printGraph();
    }
    
    //////////////////////////////////////////
    // read csr file
    //////////////////////////////////////////
  //  writeEdgelist("edgelist.txt", graph);
    params::num_partitions = (graph.num_vertex-1)/params::partition_size + 1;

    cout << "number of partitions used for reordering: " << params::num_partitions
         << " with partition size " <<  params::partition_size * sizeof(float) / 1024 << " KB" << '\n';
    //////////////////////////////////////////
    // output Degree array
    //////////////////////////////////////////
    graph.computeOutDegree();
    times = timer.elapsed();
    cout << times.wall/(1e9) << "s: outdegree is computed "  << '\n';

    ///////////////////////////////
    // re-order the graph 
    //////////////////////////////
    Orderer orderer(&graph);

   if(algo != Algo::original) {
        cpu_timer reorder_timer;
        float order_time = 0.0;

        orderer.reorder(algo);
        cout << "reordering time is: "  << reorder_timer.elapsed().wall/(1e9) << '\n';

        orderer.getNewGraph();
        cout << times.wall/(1e9) << "s: a new graph is constructed, total time is: "  << reorder_timer.elapsed().wall/(1e9) << '\n';

        if(!write_file.empty())
           writeGraph(write_file, graph);
   }

}
