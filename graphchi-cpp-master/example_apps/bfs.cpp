
/**
 * @file
 * @author  Aapo Kyrola <akyrola@cs.cmu.edu>
 * @version 1.0
 *
 * @section LICENSE
 *
 * Copyright [2012] [Aapo Kyrola, Guy Blelloch, Carlos Guestrin / Carnegie Mellon University]
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 
 *
 * @section DESCRIPTION
 *
 * Template for GraphChi applications. To create a new application, duplicate
 * this template.
 */

#include <string>

#include "graphchi_basic_includes.hpp"
#include "util/labelanalysis.hpp"

using namespace graphchi;

/**
  * Type definitions. Remember to create suitable graph shards using the
  * Sharder-program. 
  */
typedef unsigned VertexDataType;
typedef unsigned EdgeDataType;

unsigned single_source = 0;
bool converged = false;
unsigned maxlevel = 10000000;
bool scheduler = true;
bool *BFSflag = new bool[100];
/**
  * GraphChi programs need to subclass GraphChiProgram<vertex-type, edge-type> 
  * class. The main logic is usually in the update function.
  */
struct bfs : public GraphChiProgram<VertexDataType, EdgeDataType>
{

    /**
     *  Vertex update function.
     */
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext)
    {
        //std:: cout<<vertex.id()<<" "<<vertex.get_data()<<"\n";
        if (gcontext.iteration == 0)
        {
            /* On first iteration, initialize vertex (and its edges). This is usually required, because
               on each run, GraphChi will modify the data files. To start from scratch, it is easiest
               do initialize the program in code. Alternatively, you can keep a copy of initial data files. */
            // vertex.set_data(init_value);
            //	if(vertex.num_inedges() == 0)
            if (vertex.id() == single_source)
            {
                vertex.set_data(single_source);
                for (int id = 0; id < vertex.num_outedges(); id++){
                    if (!BFSflag[vertex.outedge(id)->vertexid])
                    {
                        BFSflag[vertex.outedge(id)->vertexid] = true;
                        gcontext.scheduler->add_task(vertex.outedge(id)->vertexid);
                        vertex.outedge(id)->set_data(vertex.get_data());
                        //std::cout << vertex.outedge(id)->vertex_id() << "\n";
                    }
                }
            }
            else
            {
                vertex.set_data(maxlevel);
                for (int id = 0; id < vertex.num_outedges(); id++)
                {
                    vertex.outedge(id)->set_data(vertex.get_data());
                    // if(scheduler)
                    // {
                    // 	gcontext.scheduler->add_task(vertex.outedge(id)->vertex_id());
                    // }
                }
            }
            std:: cout<<vertex.id()<<" "<<vertex.get_data()<<"\n";
        }
        else
        {
            /* Do computation */

            /* Loop over in-edges (example) */
            for (int i = 0; i < vertex.num_inedges(); i++)
            {
                // Do something
                //    value += vertex.inedge(i).get_data();
                unsigned tmpval = vertex.inedge(i)->get_data() + 1;
                if (tmpval < vertex.get_data())
                {
                    vertex.set_data(tmpval);
                }
            }
            std:: cout<<vertex.id()<<" "<<vertex.get_data()<<"\n";
            /* Loop over out-edges (example) */
            for (int i = 0; i < vertex.num_outedges(); i++)
            {
                // Do something
                // vertex.outedge(i).set_data(x)
                vertex.outedge(i)->set_data(vertex.get_data());
                if (!BFSflag[vertex.outedge(i)->vertexid])
                {
                    BFSflag[vertex.outedge(i)->vertexid] = true;
                    gcontext.scheduler->add_task(vertex.outedge(i)->vertexid);
                    //std::cout << vertex.outedge(i)->vertex_id() << "\n";
                }
            }
            /* Loop over all edges (ignore direction) */
            //for(int i=0; i < vertex.num_edges(); i++) {
            // vertex.edge(i).get_data()
            //}

            // v.set_data(new_value);
        }
    }

    /**
     * Called before an iteration starts.
     */
    void before_iteration(int iteration, graphchi_context &gcontext)
    {
        converged = iteration > 0;
    }

    /**
     * Called after an iteration has finished.
     */
    void after_iteration(int iteration, graphchi_context &gcontext)
    {
        // if (converged)
        // {
        //     logstream(LOG_INFO) << "bfs program has converged" << std::endl;
        //     gcontext.set_last_iteration(iteration);
        // }
    }

    /**
     * Called before an execution interval is started.
     */
    void before_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext)
    {
    }

    /**
     * Called after an execution interval has finished.
     */
    void after_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext)
    {
    }
};

int main(int argc, const char **argv)
{
    for (int i = 0; i < 100; i++){
        BFSflag[i] = false;
    }
    /* GraphChi initialization will read the command line 
       arguments and the configuration file. */
    graphchi_init(argc, argv);

    /* Metrics object for keeping track of performance counters
       and other information. Currently required. */
    metrics m("my-application-name");

    /* Basic arguments for application */
    std::string filename = get_option_string("file");   // Base filename
    int niters = get_option_int("niters", 1000);        // Number of iterations
    bool scheduler = get_option_int("scheduler", true); // Whether to use selective scheduling

    /* Detect the number of shards or preprocess an input to create them */
    int nshards = convert_if_notexists<EdgeDataType>(filename,
                                                     get_option_string("nshards", "auto"));

    /* Run */
    bfs program;
    graphchi_engine<VertexDataType, EdgeDataType> engine(filename, nshards, scheduler, m);
    engine.run(program, niters);
    /*analyze result*/
    m.start_time("label-analysis");
    analyze_labels<unsigned>(filename);
    /* Report execution metrics */
    metrics_report(m);
    return 0;
}
