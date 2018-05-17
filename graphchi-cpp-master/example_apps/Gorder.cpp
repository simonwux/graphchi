
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
#define __builtin_expect(a, b) b // for gcc
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>
#include <ctime>
#include "graphchi_basic_includes.hpp"
#include "Util.h"
#include "Util.cpp"
#include "UnitHeap.h"
#include "UnitHeap.cpp"
using namespace graphchi;

/**
  * Type definitions. Remember to create suitable graph shards using the
  * Sharder-program. 
  */

struct my_vertex_type
{
    int outdegree; //num_outedges()
    int indegree;  //num_inedges()
};
struct my_edge_type
{
    int first;
    int second;
};
typedef my_vertex_type VertexDataType;
typedef vid_t EdgeDataType;

int vsize;
const int maxvsize = 1000000;
VertexDataType *graph;
std::vector<int> degreevertex;
bool *BFSflag;
std::queue<int> que;
std::vector<int> order;
int retorder[maxvsize];
int reverseorder[maxvsize];
std::vector<int> tmp;
int now;
bool nonewtaskflag;
bool rcmorderflag;
bool phase1startflag;
bool phase1flag;
bool phase1subflag;
bool phase2flag;
UnitHeap unitheap(0);
std::vector<bool> popvexist;
std::vector<int> order2;
std::vector<int> retorder2;
int u1 = 0;
int tmpindex, tmpweight;
std::vector<int> zero;
int hugevertex;
clock_t time1, time2, time3, time4;
clock_t sum1 = 0, sum2 = 0, sum3 = 0;
int window = 5;
int counttt = 0;
int popv;
int popv1;
bool whileflag;
bool popvflag;
bool binaryflag;
bool popvphase1flag;
bool popvstartflag;
bool afterpopvflag;
bool phase3flag;
bool phase4flag;
bool phase4subflag;
bool phase4startflag;
bool phase4orderflag;
int phase41;
int vvv;
bool outputflagstart;
bool outputflag;
std::stack<int> st2;
ofstream Savefile;
std::string outputfilename = "";
int outputi = 0;
std::vector<int> tmporder2;
std::vector<int> tmporder3;
int rcmtmporder1[maxvsize];
int rcmtmporder2[maxvsize];
clock_t start, end, exec_start, exec_end;
double exec_sum;
int k1;
bool sumtimeflag;
bool compare(const int &a, const int &b)
{
    if (graph[a].outdegree + graph[a].indegree < graph[b].outdegree + graph[b].indegree)
        return true;
    else
        return false;
}

bool compare2(const int &a, const int &b)
{
    if (retorder2[retorder[a]] < retorder2[retorder[b]])
        return true;
    else
        return false;
}

bool compare3(const int &a, const int &b)
{
    if (retorder[a] < retorder[b])
        return true;
    else
        return false;
}
/**
  * GraphChi programs need to subclass GraphChiProgram<vertex-type, edge-type> 
  * class. The main logic is usually in the update function.
  */
struct MyGraphChiProgram : public GraphChiProgram<VertexDataType, EdgeDataType>
{

    /**
     *  Vertex update function.
     */
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &v, graphchi_context &ginfo)
    {
        if (ginfo.iteration != 0 && !outputflagstart)
            exec_start = clock();

        if (ginfo.iteration == 0)
        {
            /* On first iteration, initialize vertex (and its edges). This is usually required, because
               on each run, GraphChi will modify the data files. To start from scratch, it is easiest
               do initialize the program in code. Alternatively, you can keep a copy of initial data files. */
            graph[v.id()].indegree = v.num_inedges();
            graph[v.id()].outdegree = v.num_outedges();
            if (v.id() > vsize)
                vsize = v.id();
        }
        else
        {
            if (!rcmorderflag)
            {
                nonewtaskflag = true;
                tmp.clear();
                for (int it = 0, limit = v.num_outedges(); it < limit; it++)
                {
                    tmp.push_back(v.outedge(it)->vertex_id());
                }
                std::sort(tmp.begin(), tmp.end(), compare);
                if (tmp.size() != v.num_outedges())
                    std::cout << "tmp.size()!=graph[now].outdegree"
                              << "\n";
                for (int i = 0; i < tmp.size(); i++)
                {
                    if (BFSflag[tmp[i]] == false)
                    {
                        que.push(tmp[i]);
                        BFSflag[tmp[i]] = true;
                        order.push_back(tmp[i]);

                        //nonewtaskflag = false;
                        //ginfo.scheduler->add_task(tmp[i]);
                    }
                }
                if (!que.empty())
                {
                    nonewtaskflag = false;
                    int now = que.front();
                    que.pop();
                    ginfo.scheduler->add_task(now);
                }
            }
            else if (!phase1flag)
            {
                if (phase1subflag)
                {
                    for (int j = 0; j < v.num_outedges(); j++)
                    {
                        rcmtmporder2[j] = v.outedge(j)->vertex_id();
                    }
                    std::sort(rcmtmporder2, rcmtmporder2 + v.num_outedges(), compare3);

                    for (int j = 0; j < v.num_outedges(); j++)
                    {
                        //int w = v.outedge(j)->vertex_id();
                        int w = rcmtmporder2[j];
                        if (unitheap.update[retorder[w]] == 0)
                        {
                            unitheap.IncrementKey(retorder[w]);
                        }
                        else
                        {
#ifndef Release
                            if (unitheap.update[retorder[w]] == INT_MAX)
                                unitheap.update[retorder[w]] = INT_MAX / 2;
#endif
                            unitheap.update[retorder[w]]++;
                        }
                    }
                    phase1subflag = false;
                    ginfo.scheduler->add_task(tmpindex);
                }
                else
                {
                    if (!phase1startflag)
                    {
                        phase1startflag = true;
                        for (int i = 0; i < v.num_inedges(); i++)
                        {
                            rcmtmporder1[i] = v.inedge(i)->vertex_id();
                        }
                        std::sort(rcmtmporder1, rcmtmporder1 + v.num_inedges(), compare3);
                    }

                    phase1flag = true;
                    for (int i = u1; i < v.num_inedges(); i++)
                    {
                        phase1flag = false;
                        //int u = v.inedge(i)->vertex_id();
                        int u = rcmtmporder1[i];
                        if (i == v.num_inedges() - 1)
                            phase1flag = true;
                        if (graph[u].outdegree <= hugevertex)
                        {
                            if (unitheap.update[retorder[u]] == 0)
                            {
                                unitheap.IncrementKey(retorder[u]);
                            }
                            else
                            {
#ifndef Release
                                if (unitheap.update[retorder[u]] == INT_MAX)
                                    unitheap.update[retorder[u]] = INT_MAX / 2;
#endif
                                unitheap.update[retorder[u]]++;
                            }

                            if (i == v.num_inedges() - 1)
                                phase1flag = true;
                            if (graph[u].outdegree > 1)
                            {
                                phase1flag = false;
                                u1 = u + 1;
                                phase1subflag = true;
                                ginfo.scheduler->add_task(u);
                                break;
                            }
                        }
                    }
                }
            }
            else if (!phase2flag)
            {
                for (int i = 0; i < v.num_outedges(); i++)
                {
                    rcmtmporder1[i] = v.outedge(i)->vertex_id();
                }
                std::sort(rcmtmporder1, rcmtmporder1 + v.num_outedges(), compare3);

                for (int i = 0; i < v.num_outedges(); i++)
                {
                    //int w = v.outedge(i)->vertex_id();
                    int w = rcmtmporder1[i];
                    if (unitheap.update[retorder[w]] == 0)
                    {
                        unitheap.IncrementKey(retorder[w]);
                    }
                    else
                    {
#ifndef Release
                        if (unitheap.update[retorder[w]] == INT_MAX)
                            unitheap.update[retorder[w]] = INT_MAX / 2;
#endif
                        unitheap.update[retorder[w]]++;
                    }
                }
                phase2flag = true;
            }
            else if (popvflag && !binaryflag)
            {
                if (graph[v.id()].outdegree <= hugevertex && !popvphase1flag)
                {
                    popvphase1flag = true;

                    for (int i = 0; i < v.num_outedges(); i++)
                    {
                        rcmtmporder1[i] = v.outedge(i)->vertex_id();
                    }
                    std::sort(rcmtmporder1, rcmtmporder1 + v.num_outedges(), compare3);

                    for (int i = 0; i < v.num_outedges(); i++)
                    {
                        if (i < 0 || i >= v.num_outedges())
                        {
                            std::cout << "popvnotbinassert" << i << "\n";
                        }
                        //int w = v.outedge(i)->vertex_id();
                        int w = rcmtmporder1[i];
                        unitheap.update[retorder[w]]--;
#ifndef Release
                        if (unitheap.update[retorder[w]] == 0)
                            unitheap.update[retorder[w]] = INT_MAX / 2;
#endif
                    }
                }
                popvflag = false;

                if (!popvstartflag)
                {
                    popvstartflag = true;
                    for (int i = 0; i < v.num_inedges(); i++)
                    {
                        rcmtmporder1[i] = v.inedge(i)->vertex_id();
                    }
                    std::sort(rcmtmporder1, rcmtmporder1 + v.num_inedges(), compare3);
                }

                for (int i = popv1; i < v.num_inedges(); i++)
                {
                    popvflag = true;
                    int u = v.inedge(i)->vertex_id();//something wrong
                    //int u = rcmtmporder1[i];
                    if (i == v.num_inedges() - 1)
                        popvflag = false;
                    if (graph[u].outdegree <= hugevertex)
                    {
                        unitheap.update[retorder[u]]--;
#ifndef Release
                        if (unitheap.update[retorder[u]] == 0)
                            unitheap.update[retorder[u]] = INT_MAX / 2;
#endif
                        if (graph[u].outdegree > 1)
                        {
                            popvflag = true;
                            popv1 = i + 1;
                            binaryflag = true;
                            ginfo.scheduler->add_task(u);
                            break;
                        }
                    }
                }
            }
            else if (popvflag && binaryflag)
            {
                bool findflag = false;
                int low = 0;
                int high = v.num_outedges() - 1;
                while (low <= high)
                {
                    int mid = (low + high) / 2;
                    if (v.outedge(mid)->vertex_id() > reverseorder[vvv])
                        high = mid - 1;
                    else if (v.outedge(mid)->vertex_id() < reverseorder[vvv])
                        low = mid + 1;
                    else
                    {
                        findflag = true;
                        break;
                    };
                }
                if (!findflag)
                {
                    for (int j = 0; j < v.num_outedges(); j++)
                    {
                        rcmtmporder2[j] = v.outedge(j)->vertex_id();
                    }
                    std::sort(rcmtmporder2, rcmtmporder2 + v.num_outedges(), compare3);
                    for (int j = 0; j < v.num_outedges(); j++)
                    {

                        if (j < 0 || j >= v.num_outedges())
                        {
                            std::cout << "popvbinassert" << j << "\n";
                        }
                        //int w = v.outedge(j)->vertex_id();
                        int w = rcmtmporder2[j];
                        unitheap.update[retorder[w]]--;
#ifndef Release
                        if (unitheap.update[retorder[w]] == 0)
                            unitheap.update[retorder[w]] = INT_MAX / 2;
#endif
                    }
                }
                else
                {
                    popvexist[retorder[v.id()]] = true;
                }
                ginfo.scheduler->add_task(popv);
                binaryflag = false;
            }
            else if (!phase3flag)
            {
                for (int i = 0; i < v.num_outedges(); i++)
                {
                    rcmtmporder1[i] = v.outedge(i)->vertex_id();
                }
                std::sort(rcmtmporder1, rcmtmporder1 + v.num_outedges(), compare3);

                for (int i = 0; i < v.num_outedges(); i++)
                {
                    if (i < 0 || i >= v.num_outedges())
                    {
                        std::cout << "phase3assert" << i << "\n";
                    }
                    //int w = v.outedge(i)->vertex_id();
                    int w = rcmtmporder1[i];
                    if (__builtin_expect(unitheap.update[retorder[w]] == 0, 0))
                    {
                        unitheap.IncrementKey(retorder[w]);
                    }
                    else
                    {
#ifndef Release
                        if (unitheap.update[retorder[w]] == INT_MAX)
                            unitheap.update[retorder[w]] = INT_MAX / 2;
#endif
                        unitheap.update[retorder[w]]++;
                    }
                }
                phase3flag = true;
            }
            else if (!phase4flag)
            {
                if (phase4subflag)
                {
                    for (int j = 0; j < v.num_outedges(); j++)
                    {
                        rcmtmporder2[j] = v.outedge(j)->vertex_id();
                    }
                    std::sort(rcmtmporder2, rcmtmporder2 + v.num_outedges(), compare3);

                    for (int j = 0; j < v.num_outedges(); j++)
                    {

                        if (j < 0 || j >= v.num_outedges())
                        {
                            std::cout << "phase4assert" << j << "\n";
                        }
                        //int w = v.outedge(j)->vertex_id();
                        int w = rcmtmporder2[j];
                        if (__builtin_expect(unitheap.update[retorder[w]] == 0, 0))
                        {
                            unitheap.IncrementKey(retorder[w]);
                        }
                        else
                        {
#ifndef Release
                            if (unitheap.update[retorder[w]] == INT_MAX)
                                unitheap.update[retorder[w]] = INT_MAX / 2;
#endif
                            unitheap.update[retorder[w]]++;
                        }
                    }
                    phase4subflag = false;
                    ginfo.scheduler->add_task(reverseorder[vvv]);
                }
                else
                {
                    if (!phase4orderflag)
                    {
                        phase4orderflag = true;
                        for (int i = 0; i < v.num_inedges(); i++)
                        {
                            rcmtmporder1[i] = v.inedge(i)->vertex_id();
                        }
                        std::sort(rcmtmporder1, rcmtmporder1 + v.num_inedges(), compare3);
                    }

                    phase4flag = true;
                    for (int i = phase41; i < v.num_inedges(); i++)
                    {
                        phase4flag = false;
                        //int u = v.inedge(i)->vertex_id();
                        int u = rcmtmporder1[i];
                        if (i == v.num_inedges() - 1)
                            phase4flag = true;
                        if (graph[u].outdegree <= hugevertex)
                        {
                            if (__builtin_expect(unitheap.update[retorder[u]] == 0, 0))
                            {
                                unitheap.IncrementKey(retorder[u]);
                            }
                            else
                            {
#ifndef Release
                                if (unitheap.update[retorder[u]] == INT_MAX)
                                    unitheap.update[retorder[u]] = INT_MAX / 2;
#endif
                                unitheap.update[retorder[u]]++;
                            }

                            if (popvexist[retorder[u]] == false)
                            {
                                if (graph[u].outdegree > 1)
                                {
                                    phase41 = i + 1;
                                    phase4subflag = true;
                                    phase4flag = false;
                                    ginfo.scheduler->add_task(u);
                                    break;
                                }
                            }
                            else
                            {
                                popvexist[retorder[u]] = false;
                            }
                        }
                    }
                }
            }
            else if (!outputflag)
            {
                tmporder3.clear();
                tmporder3.resize(v.num_outedges());
                for (int i = 0; i < v.num_outedges(); i++)
                    tmporder3[i] = v.outedge(i)->vertex_id();
                std::sort(tmporder3.begin(), tmporder3.end(), compare2);
                Savefile.open(outputfilename.c_str(), ios::app);
                for (int j = 0; j < v.num_outedges(); j++)
                {
                    int u = retorder2[retorder[v.id()]];
                    int v2 = retorder2[retorder[tmporder3[j]]];
                    //int v2 = retorder2[retorder[v.outedge(j)->vertex_id()]];
                    Savefile << u << '\t' << v2 << "\n";
                }
                Savefile.close();
                if (outputi < vsize - 1)
                {
                    outputi++;
                    ginfo.scheduler->add_task(tmporder2[outputi]);
                }
            }
        }

        if (ginfo.iteration != 0 && !outputflagstart)
        {
            exec_end = clock();
            exec_sum = exec_sum + (double)(exec_end - exec_start) / CLOCKS_PER_SEC;
        }
    }

    /**
     * Called before an iteration starts.
     */
    void before_iteration(int iteration, graphchi_context &ginfo)
    {
        if (iteration == 0)
        {
            /* initialize  each vertex with its own lable */
            graph = new VertexDataType[maxvsize];
            for (int i = 0; i < maxvsize; i++)
            {
                rcmtmporder1[i] = 0;
                rcmtmporder2[i] = 0;
                graph[i] = {
                    outdegree : 0,
                    indegree : 0
                };
            }
            vsize = 0;
            nonewtaskflag = true;
            rcmorderflag = false;
            phase1startflag = false;
            phase1flag = false;
            phase1subflag = false;
            phase2flag = false;
            whileflag = false;
            popvflag = false;
            binaryflag = false;
            popvphase1flag = false;
            popvstartflag = false;
            afterpopvflag = false;
            phase3flag = false;
            phase4flag = false;
            phase4subflag = false;
            phase4startflag = false;
            phase4orderflag = false;
            phase41 = 0;
            outputflagstart = false;
            outputflag = false;
            exec_sum = 0;
            k1 = 1;
            sumtimeflag = false;
        }
        //nonewtaskflag = true;
    }

    /**
     * Called after an iteration has finished.
     */
    void after_iteration(int iteration, graphchi_context &gcontext)
    {
        if (!outputflagstart)
            exec_start = clock();

        if (iteration == 0)
        {
            vsize = vsize + 1;
            std::cout << "vsize: " << vsize << "\n";

            BFSflag = new bool[vsize];
            memset(BFSflag, 0, sizeof(bool) * vsize);
            degreevertex.resize(vsize);
            for (int i = 0; i < vsize; i++)
            {
                degreevertex[i] = i;
            }
            std::sort(degreevertex.begin(), degreevertex.end(), compare);

            int i = degreevertex[0];
            BFSflag[i] = true;
            order.push_back(i);
            gcontext.scheduler->add_task(i);
        }
        else
        {
            if (!rcmorderflag)
            {
                if (nonewtaskflag)
                {
                    bool flag = true;
                    for (int k = k1; k < vsize; k++)
                    {
                        int i = degreevertex[k];
                        if (BFSflag[i] == false)
                        {
                            k1 = k + 1;
                            flag = false;
                            nonewtaskflag = false;
                            BFSflag[i] = true;
                            order.push_back(i);
                            gcontext.scheduler->add_task(i);
                            break;
                        }
                    }
                    if (flag)
                    {
                        if (order.size() != vsize)
                        {
                            std::cout << "order.size()!=vsize"
                                      << "\n";
                        }

                        for (int i = 0; i < order.size(); i++)
                        {
                            retorder[order[i]] = order.size() - 1 - i;
                        }
                        for (int i = 0; i < vsize; i++)
                        {
                            reverseorder[retorder[i]] = i;
                        }
                        rcmorderflag = true;

                        // Savefile.open("wikircmordergc.txt");
                        // for (int i = 0;i<vsize;i++){
                        //     Savefile<<i<<"\t"<<retorder[i]<<"\n";
                        // }
                        // Savefile.close();

                        std::cout << "readGraph is complete." << endl;
                        end = clock();
                        std::cout << "Time Cost: " << (double)(end - start) / CLOCKS_PER_SEC << "\n";

                        //greedyorderstart
                        start = clock();
                        
                        unitheap.reshape(vsize);
                        popvexist.reserve(vsize);
                        for (int i = 0; i < popvexist.size(); i++)
                            popvexist[i] = false;
                        zero.reserve(10000);
                        order2.reserve(vsize);
                        hugevertex = sqrt((double)vsize);
                        for (int j = 0; j < vsize; j++)
                        {
                            int i = reverseorder[j];
                            unitheap.LinkedList[j].key = graph[i].indegree;
                            unitheap.update[j] = -graph[i].indegree;
                        }
                        unitheap.ReConstruct();
                        
                        tmpweight = -1;
                        for (int j = 0; j < vsize; j++)
                        {
                            int i = reverseorder[j];
                            if (graph[i].indegree > tmpweight)
                            {
                                tmpweight = graph[i].indegree;
                                tmpindex = j;
                            }
                            else if (graph[i].indegree + graph[i].outdegree == 0)
                            {
                                unitheap.update[j] = INT_MAX / 2;
                                zero.push_back(j);
                                unitheap.DeleteElement(j);
                            }
                        }

                        order2.push_back(tmpindex);
                        unitheap.update[tmpindex] = INT_MAX / 2;
                        unitheap.DeleteElement(tmpindex);
                        gcontext.scheduler->add_task(reverseorder[tmpindex]);
                        std::cout<<"heapsize"<<unitheap.heapsize<<"\n";
                    }
                }
            }
            else if (phase1flag && !phase2flag)
            {
                if (graph[reverseorder[tmpindex]].outdegree <= hugevertex)
                {
                    gcontext.scheduler->add_task(reverseorder[tmpindex]);
                }
                else
                    phase2flag = true;

                if ((counttt < vsize - 1 - zero.size()))
                {
                    whileflag = false;
                }
                else
                {
                    whileflag = true;
                }
            }
            if (phase1flag && phase2flag && !whileflag)
            {
                if (!popvflag && afterpopvflag && phase3flag && phase4startflag && phase4flag)
                {
#ifndef Release
                    time4 = clock();
                    sum1 += time2 - time1;
                    sum2 += time3 - time2;
                    sum3 += time4 - time3;
#endif

                    if ((counttt >= vsize - 1 - zero.size()))
                        whileflag = true;
                    if (!whileflag)
                    {
                        afterpopvflag = false;
                        popvstartflag = false;
                        popvflag = false;
                        popvphase1flag = false;
                        popv1 = 0;
                        binaryflag = false;
                        phase3flag = false;
                        phase4flag = false;
                        phase4subflag = false;
                        phase4startflag = false;
                        phase4orderflag = false;
                        phase41 = 0;
                    }
                    
                }
                if (!afterpopvflag)
                {
#ifndef Release
                    if (counttt % 1000000 == 0)
                    {
                        std::cout << counttt << "\n";
                        std::cout << "sum1: " << sum1 << "\n";
                        std::cout << "sum2: " << sum2 << "\n";
                        std::cout << "sum3: " << sum3 << "\n";
                        std::cout << endl;
                        sum1 = sum2 = sum3 = 0;
                    }

                    time1 = clock();
#endif
                    //if (iteration > 150000) std:: cout<<"extractmaxstart\n";
                    vvv = unitheap.ExtractMax();
                    counttt++;

                    //if (iteration > 150000) std:: cout<<"extractmaxends\n";

#ifndef Release
                    time2 = clock();
#endif

                    order2.push_back(vvv);
                    unitheap.update[vvv] = INT_MAX / 2;

                    if (counttt - window >= 0)
                        popv = order2[counttt - window];
                    else
                        popv = -1;

                    if (popv >= 0)
                    {
                        popv1 = 0;
                        popvflag = true;
                        popvstartflag = false;
                        popvphase1flag = false;
                        binaryflag = false;

                        gcontext.scheduler->add_task(reverseorder[popv]);
                    }

                    afterpopvflag = true;

                }
                if (!popvflag && afterpopvflag && !phase3flag)
                {

#ifndef Release
                    time3 = clock();
#endif
                    if (graph[reverseorder[vvv]].outdegree <= hugevertex)
                    {
                        gcontext.scheduler->add_task(reverseorder[vvv]);
                    }
                    else
                        phase3flag = true;
                }
                if (!popvflag && afterpopvflag && phase3flag && !phase4startflag)
                {
                    phase4startflag = true;
                    phase4orderflag = false;
                    phase41 = 0;
                    
                    gcontext.scheduler->add_task(reverseorder[vvv]);
                }
            }

            if (phase2flag && whileflag && !outputflagstart)
            {
                if (!outputflagstart)
                {
                    exec_end = clock();
                    exec_sum = exec_sum + (double)(exec_end - exec_start) / CLOCKS_PER_SEC;
                }

                std::cout << "whilefinished\n";
                order2.insert(order2.end() - 1, zero.begin(), zero.end());

#ifndef Release
                tmporder2 = order2;
                std::sort(tmporder2.begin(), tmporder2.end());
                for (int i = 0; i < tmporder2.size() - 1; i++)
                {
                    if (tmporder2[i] == tmporder2[i + 1])
                    {
                        std::cout << "same elements: " << tmporder2[i] << endl;
                    }
                }
                for (int i = 0; i < tmporder2.size(); i++)
                {
                    if (tmporder2[i] != i)
                    {
                        std::cout << "tmporder" << tmporder2[i] << '\t' << i << endl;
                    }
                }
                vector<int>().swap(tmporder2);
#endif
                retorder2.clear();
                retorder2.resize(vsize);
                for (int i = 0; i < vsize; i++)
                {
                    retorder2[order2[i]] = i;
                }
                end = clock();
                std::cout << "ReOrdered Time Cost: " << (double)(end - start) / CLOCKS_PER_SEC << "\n";
                std::cout << "Begin Output the Reordered Graph"
                          << "\n";
                
                outputflagstart = true;
                

                tmporder2.clear();
                tmporder2.resize(vsize);
                for (int i = 0; i < vsize; i++)
                    tmporder2[i] = i;
                std::sort(tmporder2.begin(), tmporder2.end(), compare2);
                Savefile.open(outputfilename.c_str());
                Savefile.close();
                outputi = 0;
                for (int i = 0; i < 1; i++)
                {
                    gcontext.scheduler->add_task(tmporder2[i]);
                }
                
            }
        }

        if (!outputflagstart)
        {
            exec_end = clock();
            exec_sum = exec_sum + (double)(exec_end - exec_start) / CLOCKS_PER_SEC;
        }
        if (outputflagstart && !sumtimeflag)
        {
            sumtimeflag = true;
            std::cout << "exec_sum: " << exec_sum << '\n';
        }
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

    start = clock();
    /* GraphChi initialization will read the command line 
       arguments and the configuration file. */
    graphchi_init(argc, argv);

    /* Metrics object for keeping track of performance counters
       and other information. Currently required. */
    metrics m("Gorder");

    /* Basic arguments for application */
    std::string filename = get_option_string("file"); // Base filename
    //std::cout<<filename<<"\n";
    for (int i = 0; i < filename.size(); i++)
    {
        if (filename[i] == '.')
            break;
        else
            outputfilename = outputfilename + filename[i];
    }
    outputfilename = outputfilename + "_gorder.txt";
    std::cout << "outputfile: " << outputfilename << "\n";
    int niters = get_option_int("niters", 100000000); // Number of iterations
    bool scheduler = true;                            // Whether to use selective scheduling

    /* Detect the number of shards or preprocess an input to create them */
    int nshards = convert_if_notexists<EdgeDataType>(filename,
                                                     get_option_string("nshards", "auto"));
    /* Run */
    MyGraphChiProgram program;
    graphchi_engine<VertexDataType, EdgeDataType> engine(filename, nshards, scheduler, m);
    engine.run(program, niters);

    /* Report execution metrics */
    metrics_report(m);
    return 0;
}
