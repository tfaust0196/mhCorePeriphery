//
//  main.cpp
//  Code for "Inference of Hierarchical Core–Periphery Structure in Temporal Networks"
//
//  Created by Theodore Faust
//

#include <iostream>
#include <fstream>
#include <random>
#include <stdexcept>
#include <vector>
#include <cmath>
#include "MathMatrix.h"
#include "MathVector.h"
#include <algorithm>
#include <complex>
#include <numbers>
#include <chrono>
#include <cmath>
using namespace std::chrono;
using namespace std::complex_literals;
#ifdef _MSC_VER
#include "iso646.h"      // So "and" is equivalenced to &&
typedef unsigned int uint;   // Define uint to be unsigned int
#endif
long double logprobratio(MathMatrix A, std::vector<MathMatrix> G, int node, int group, int layer);
long double logbinom(double n, double k);
std::vector<std::complex<double>> coeffs(int n);
std::vector<std::vector<double>> corrections(int n);
long double logNewCorrection(MathVector prevGroups, MathVector nextGroups, int numGroups, std::vector<std::vector<double>> corrections);
long double lognewcorrectionNovel(MathMatrix G, int node, int layer, int numGroups, std::vector<std::vector<double>> corrections);
bool stepNoGpSizeBiasNew(std::vector<MathMatrix> A, std::vector<MathMatrix> G, int node, int group, int layer, int numGroups, std::vector<std::vector<double>> corrections);
long double logprobrationew(MathMatrix A, std::vector<MathMatrix> G, std::vector<MathMatrix> Gtemp, std::vector<std::vector<double>> corrections);
double logfactorial2(double n);
int main(int argc, const char * argv[]) // reads a network and runs the MCMC method on it

{
    // used to generate uniform random numbers in [0,1]
    std::random_device rd3;
    std::mt19937 gen3(rd3());
    std::uniform_real_distribution<> dis3(0.0, 1.0);
    int numGroups = 4; //number of groups
    size_t numNodes = 76; //number of nodes
    std::vector<std::vector<double>> corrs = corrections(numNodes);
    size_t sizeA = numNodes;
    size_t layersA = 20; // number of layers
    size_t numLayers = layersA;
    int runs = 1; //deprecated (performed multiple runs via multiple simultaneous instantiations of the code)
    std::string dataset = "luke"; //start of filename
    double mnMoveProb = 1e-3; // probability p of performing a multi-node move
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(0,numLayers-1);
    // used to generate random group
    std::random_device rd4;
    std::mt19937 gen4(rd4());
    std::uniform_int_distribution<> distrib4(0,numNodes-1);
    long numMCSteps = 1e6; // number of Monte Carlo steps for Markov chain Monte Carlo process
    
    
    std::vector<MathMatrix> A(layersA,MathMatrix(sizeA,sizeA));
    for (size_t l = 0; l < layersA; ++l)
    {
        std::string myText;
        std::string delimiter = " ";
        std::ifstream myfile3(dataset + "layer" + std::to_string(l+1) + ".txt");
        for (size_t i = 0; i < numNodes; ++i)
        {
            getline(myfile3,myText);
            for (size_t j = 0; j < numNodes; ++j)
            {
                A[l](i,j) = std::stod(myText.substr(0, myText.find(delimiter)));
                myText.erase(0, myText.find(delimiter) + delimiter.length());
            }
        }
        myfile3.close();
    }
    std::vector<MathMatrix> groupAssignsFinal(numLayers,MathMatrix(numNodes,runs));
    for (size_t run = 0; run < runs; ++run)
    {
        std::vector<MathMatrix> G(numGroups,MathMatrix(numNodes,numLayers)); // list of matrices of indicator variables: G[r](i,l) = g^r_{(i,l)}
        G[0] = MathMatrix(numNodes,numLayers,1.0); // every node-layer is always an element of group 0
        for (size_t g = 1; g < numGroups; ++g)
        {
            for (size_t k = 0; k < numNodes; ++k)
            {
                for (size_t l = 0; l < numLayers; ++l)
                {
                    if (dis3(gen3) > 0.5)
                    {
                        G[g](k,l) = 1.0;
                    }
                }
            }
        }
        std::vector<MathMatrix> intGroupAssignML(numLayers,MathMatrix(numNodes,numMCSteps/1e4));
        for (long i = 0; i < numMCSteps; ++i)
        {
            if (i % 10000 == 0)
            {
                std::cout << i << "\n";
                if (numGroups > 1.5)
                {
                    std::cout << "G[1]: " << G[1];
                }
                if (numGroups > 2.5)
                {
                    std::cout << "G[2]: " << G[2];
                }
            }
            std::random_device rd2;
            std::mt19937 gen2(rd2());
            double randomNum2 = dis3(gen3);
            std::uniform_int_distribution<> distrib2(1,std::max(numGroups,2)-1);
            int layer = distrib(gen); // generate random layer
            int group = distrib2(gen2); // generate random group
            if (dis3(gen3) < mnMoveProb)
            {
                // perform a multi-node move
                std::random_device rdev;
                std::mt19937 rggen(rdev());
                std::uniform_int_distribution<> rgdis(0,1);
                MathVector g1(numGroups,1.0);
                MathVector g2(numGroups,1.0);
                for (size_t z = 1; z < numGroups; ++z)
                {
                    g1(z) = (double)rgdis(rggen);
                    g2(z) = (double)rgdis(rggen);
                }
                std::vector<MathMatrix> Gtemp = G;
                size_t lay = layer;
                for (size_t y = 0; y < numNodes; ++y)
                {
                    double temp1 = 0;
                    double temp2 = 0;
                    for (size_t x = 0; x < numGroups; ++x)
                    {
                        temp1 += fabs(G[x](y,lay) - g1(x));
                        temp2 += fabs(G[x](y,lay) - g2(x));
                    }
                    if (temp1 < 0.5)
                    {
                        for (size_t x = 0; x < numGroups; ++x)
                        {
                            Gtemp[x](y,lay) = g2(x);
                        }
                    }
                    else if (temp2 < 0.5)
                    {
                        for (size_t x = 0; x < numGroups; ++x)
                        {
                            Gtemp[x](y,lay) = g1(x);
                        }
                    }
                }
                long double logAcceptanceProb = logprobrationew(A[0], G, Gtemp, corrs);
                if (dis3(gen3) < std::exp(logAcceptanceProb))
                {
                    G = Gtemp;
                }
            }
            else
            {
                if (layer == 0)
                {
                    double randomNum2 = dis3(gen3);
                    if (randomNum2 < 1.0/(2.0*numGroups*(numNodes+1)))
                    {
                        double randomNum = dis3(gen3);
                        MathVector zeros = MathVector(numNodes,0.0);
                        double corrEmptyGroup = std::exp((numLayers-1)*logNewCorrection(zeros, zeros,  numGroups, corrs));
                        if (randomNum < corrEmptyGroup)
                        {
                            numGroups++;
                            G.insert(G.begin() + group, MathMatrix(numNodes,numLayers,0.0));
                        }
                    }
                    else
                    {
                        if (numGroups > 1.5)
                        {
                            double randomNum = dis3(gen3); //generate uniform random number in [0,1]
                            if (randomNum > 0.5) //remove node from group with prob 0.5
                            {
                                std::vector<int> nodesInGroup; // vector to hold indices of node-layers for the selected layer which are in the selected group
                                
                                // loop over nodes in selected layer and check which are in the selected group
                                for (size_t n = 0; n < numNodes; ++n)
                                {
                                    if (G[group](n,layer) > 0.9)
                                    {
                                        nodesInGroup.push_back(n);
                                    }
                                }
                                if (nodesInGroup.size() > 0) // can only rmove a node if there is a node in the group
                                {
                                    // choose node uniformly at random
                                    std::random_device rd5;
                                    std::mt19937 gen5(rd5());
                                    std::uniform_int_distribution<> distrib5(0,nodesInGroup.size()-1);
                                    int randIndex = distrib5(gen5);
                                    int node = nodesInGroup[randIndex];
                                    // determine whether or not to accept move
                                    if (stepNoGpSizeBiasNew(A,G,node,group,layer,numGroups,corrs))
                                    {
                                        G[group](node,layer) = -1.0*(G[group](node,layer)-1.0); // if accepted, "flip" indicator variable of selected node-layer for given group
                                    }
                                }
                                else
                                {
                                    double randomNum = dis3(gen3);
                                    MathVector zeros = MathVector(numNodes,0.0);
                                    double probEmptyGroup = std::exp(-1.0*(numLayers-1)*logNewCorrection(zeros, zeros,  numGroups, corrs));
                                    if ((randomNum < probEmptyGroup) && (numGroups > 1))
                                    {
                                        numGroups--;
                                        G.erase(G.begin() + group);
                                    }
                                }
                            }
                            else //adding node (essentially the same as the above)
                            {
                                std::vector<int> nodesNotInGroup;
                                for (size_t n = 0; n < numNodes; ++n)
                                {
                                    if (G[group](n,layer) < 0.1)
                                    {
                                        nodesNotInGroup.push_back(n);
                                    }
                                }
                                if (nodesNotInGroup.size() > 0)
                                {
                                    std::random_device rd5;
                                    std::mt19937 gen5(rd5());
                                    std::uniform_int_distribution<> distrib5(0,nodesNotInGroup.size()-1);
                                    int randIndex = distrib5(gen5);
                                    int node = nodesNotInGroup[randIndex];
                                    if (stepNoGpSizeBiasNew(A,G,node,group,layer,numGroups,corrs))
                                    {
                                        G[group](node,layer) = -1.0*(G[group](node,layer)-1.0);
                                    }
                                }
                            }
                        }
                    }
                }
                else // if the random layer is beyond the first
                {
                    if (numGroups > 1.5)
                    {
                        int node = distrib4(gen4); // choose node uniformly at random (from all possible nodes)
                        if (stepNoGpSizeBiasNew(A,G,node,group,layer,numGroups,corrs))
                        {
                            G[group](node,layer) = -1.0*(G[group](node,layer)-1.0);
                        }
                    }
                }
                if (i % 10000 == 0)
                {
                    for (size_t n = 0; n < numNodes; ++n)
                    {
                        for (size_t l = 0; l < numLayers; ++l)
                        {
                            int groupNum = 0;
                            for (size_t g = 1; g < numGroups; ++g)
                            {
                                if (G[g](n,l) > 0.9)
                                {
                                    groupNum += std::pow(2,g-1);
                                }
                            }
                            intGroupAssignML[l](n,i/10000) = groupNum;
                        }
                    }
                }
                else
                {
                    for (size_t n = 0; n < numNodes; ++n)
                    {
                        for (size_t l = 0; l < numLayers; ++l)
                        {
                            int groupNum = 0;
                            for (size_t g = 1; g < numGroups; ++g)
                            {
                                if (G[g](n,l) > 0.9)
                                {
                                    groupNum += std::pow(2,g-1);
                                }
                            }
                        }
                    }
                }
            }
            for (size_t l = 0; l < numLayers; ++l)
            {
                std::ofstream myfile3(dataset + "layer" + std::to_string(l) + "groupAssigns.txt", std::ios::trunc);
                myfile3 << intGroupAssignML[l];
            }
        }
    }
    return 0;
}

long double logprobratio(MathMatrix A, std::vector<MathMatrix> G, int node, int group, int layer)
/*
 Computes a term necessary in the computation of P(A|G',k)/P(A|G,k) where G' contains the same group assignments as G with the exception of G{group}[node,layer]
 A: adjacency matrix for given layer
 G: group assignments (for all node-layers)
 node: given node (of the node-layer whose group assignment will change)
 group: given group (for which the group assignment changes)
 layer: given layer (of the node-layer whose group assignment will change)
 */
{
    size_t numGroups = G.size();
    size_t numNodes = G[0].getRowSize();
    MathVector oldt = MathVector(numGroups);
    MathVector newt = MathVector(numGroups);
    MathVector oldm = MathVector(numGroups);
    MathVector newm = MathVector(numGroups);
    for (size_t i = 0; i < numNodes; ++i)
    {
        for (size_t j = 0; j < i; ++j)
        {
            bool stoppingCond = false;
            size_t g = numGroups - 1;
            if ((i != node) && (j != node))
            {
                while (not(stoppingCond))
                {
                    if ((G[g](i,layer) > 0.9) && (G[g](j,layer) > 0.9))
                    {
                        oldt(g) = oldt(g) + 1.0;
                        newt(g) = newt(g) + 1.0;
                        if (A(i,j) > 0.9)
                        {
                            oldm(g) = oldm(g) + 1.0;
                            newm(g) = newm(g) + 1.0;
                        }
                        stoppingCond = true;
                    }
                    else
                    {
                        g--;
                    }
                }
            }
            else
            {
                bool stoppingCondOld = false;
                bool stoppingCondNew = false;
                size_t g = numGroups - 1;
                if (i == node)
                {
                    while (not(stoppingCondOld) || not(stoppingCondNew))
                    {
                        if (g != group)
                        {
                            if ((G[g](i,layer) > 0.9) && (G[g](j,layer) > 0.9))
                            {
                                if (not(stoppingCondOld))
                                {
                                    oldt(g) = oldt(g) + 1.0;
                                }
                                if (not(stoppingCondNew))
                                {
                                    newt(g) = newt(g) + 1.0;
                                }
                                if (A(i,j) > 0.9)
                                {
                                    if (not(stoppingCondOld))
                                    {
                                        oldm(g) = oldm(g) + 1.0;
                                    }
                                    if (not(stoppingCondNew))
                                    {
                                        newm(g) = newm(g) + 1.0;
                                    }
                                }
                                stoppingCondOld = true;
                                stoppingCondNew = true;
                            }
                            else
                            {
                                g--;
                            }
                        }
                        else
                        {
                            if ((G[g](i,layer) > 0.9) && (G[g](j,layer) > 0.9))
                            {
                                if (not(stoppingCondOld))
                                {
                                    oldt(g) = oldt(g) + 1.0;
                                }
                                if (A(i,j) > 0.9)
                                {
                                    if (not(stoppingCondOld))
                                    {
                                        oldm(g) = oldm(g) + 1.0;
                                    }
                                }
                                stoppingCondOld = true;
                            }
                            else if ((G[g](i,layer) < 0.1) && (G[g](j,layer) > 0.9))
                            {
                                if (not(stoppingCondNew))
                                {
                                    newt(g) = newt(g) + 1.0;
                                }
                                if (A(i,j) > 0.9)
                                {
                                    if (not(stoppingCondNew))
                                    {
                                        newm(g) = newm(g) + 1.0;
                                    }
                                }
                                stoppingCondNew = true;
                            }
                            g--;
                        }
                    }
                }
                if (j == node)
                {
                    while (not(stoppingCondOld) || not(stoppingCondNew))
                    {
                        if (g != group)
                        {
                            if ((G[g](i,layer) > 0.9) && (G[g](j,layer) > 0.9))
                            {
                                if (not(stoppingCondOld))
                                {
                                    oldt(g) = oldt(g) + 1.0;
                                }
                                if (not(stoppingCondNew))
                                {
                                    newt(g) = newt(g) + 1.0;
                                }
                                if (A(i,j) > 0.9)
                                {
                                    if (not(stoppingCondOld))
                                    {
                                        oldm(g) = oldm(g) + 1.0;
                                    }
                                    if (not(stoppingCondNew))
                                    {
                                        newm(g) = newm(g) + 1.0;
                                    }
                                }
                                stoppingCondOld = true;
                                stoppingCondNew = true;
                            }
                            else
                            {
                                g--;
                            }
                        }
                        else
                        {
                            if ((G[g](i,layer) > 0.9) && (G[g](j,layer) > 0.9))
                            {
                                if (not(stoppingCondOld))
                                {
                                    oldt(g) = oldt(g) + 1.0;
                                }
                                if (A(i,j) > 0.9)
                                {
                                    if (not(stoppingCondOld))
                                    {
                                        oldm(g) = oldm(g) + 1.0;
                                    }
                                }
                                stoppingCondOld = true;
                            }
                            else if ((G[g](j,layer) < 0.1) && (G[g](i,layer) > 0.9))
                            {
                                if (not(stoppingCondNew))
                                {
                                    newt(g) = newt(g) + 1.0;
                                }
                                if (A(i,j) > 0.9)
                                {
                                    if (not(stoppingCondNew))
                                    {
                                        newm(g) = newm(g) + 1.0;
                                    }
                                }
                                stoppingCondNew = true;
                            }
                            g--;
                        }
                    }
                }
            }
        }
    }
    long double out = 0.0;
    for (size_t g = 0; g < numGroups; ++g)
    {
        out += (logfactorial2(newm(g)) + logfactorial2(newt(g)-newm(g))) - logfactorial2(newt(g)+1);
        out -= (logfactorial2(oldm(g)) + logfactorial2(oldt(g)-oldm(g))) - logfactorial2(oldt(g)+1);
    }
    return out;
}

long double logbinom(double n, double k)
//uses expression for log(n!) to compute an approximation to log(n choose k)
{
    return logfactorial2(n) - logfactorial2(k) - logfactorial2(n-k);
}
long double logNewCorrection(MathVector prevGroups, MathVector nextGroups, int numGroups, std::vector<std::vector<double>> corrections)
/*
 Computes log P(g^r_l|g^r_{l-1}) for the novel approach by using (4)
 prevGroups: vector containing the community assignments in layer l-1 (i.e., g^r_{l-1})
 nextGroups: vector containing the community assignments in layer l (i.e., g^r_l)
 numGroups: given number of communities k
 corrections: the values of J(k_1,k_2) for all 0 <= k_1 < k_2
 */
{
    long double out = 0.0;
    size_t numNodes = prevGroups.size();
    for (size_t i = 0; i < numGroups; ++i)
    {
        long double temp = 0.0;
        int groupSize = 0;
        int numSame = 0;
        MathVector groupSizes = MathVector(numGroups-1);
        for (size_t j = 0; j < numNodes; ++j)
        {
            if (fabs(prevGroups(j) - i) < 0.1)
            {
                groupSize++;
                if (fabs(nextGroups(j) - i) < 0.1)
                {
                    numSame++;
                }
                else if (nextGroups(j) < i)
                {
                    groupSizes((int)nextGroups(j)) = groupSizes((int)nextGroups(j)) + 1.0;
                }
                else
                {
                    groupSizes((int)nextGroups(j)-1) = groupSizes((int)nextGroups(j)-1) + 1.0;
                }
            }
        }
        out += corrections[groupSize][groupSize-numSame] - logbinom(groupSize-numSame+numGroups-2, numGroups-2) - logfactorial2(groupSize);
        for (size_t k = 0; k < numGroups-1; ++k)
        {
            out += logfactorial2(groupSizes(k));
        }
        out += logfactorial2(numSame);
    }
    return out;
}
bool stepNoGpSizeBiasNew(std::vector<MathMatrix> A, std::vector<MathMatrix> G, int node, int group, int layer, int numGroups, std::vector<std::vector<double>> corrections)
/*
 Determines whether to accept or reject a standard move which changes the group-r assignment of a given node in a given layer
 A: adjacency matrix for given layer
 G: group assignments (for all node-layers)
 node: given node (of the node-layer whose group assignment will change)
 group: given group (for which the group assignment changes)
 layer: given layer (of the node-layer whose group assignment will change)
 numGroups: given number of communities k
 corrections: the values of J(k_1,k_2) for all 0 <= k_1 < k_2
 */
{
    long double logpr = logprobratio(A[layer], G, node, group, layer);
    long double logcr = lognewcorrectionNovel(G[group], node, layer, 2, corrections);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    long double logstepProb = logpr + logcr;
    long double stepProb = std::exp(logstepProb);
    double randomNum = dis(gen);
    return (randomNum < stepProb);
}
//
long double lognewcorrectionNovel(MathMatrix G, int node, int layer, int numGroups, std::vector<std::vector<double>> corrections)
/*
 Computes a term necessary in the computation of P(A|G',k)/P(A|G,k) where G' contains the same group assignments as G with the exception of G{group}[node,layer]
 G: group assignments (for all node-layers)
 node: given node (of the node-layer whose group assignment will change)
 layer: given layer (of the node-layer whose group assignment will change)
 numGroups: given number of communities k
 corrections: the values of J(k_1,k_2) for all 0 <= k_1 < k_2
 */
{
    size_t numNodes = G.getRowSize();
    size_t numLayers = G.getColSize();
    long double logcorrectionnumprev = 0;
    long double logcorrectiondenomprev = 0;
    long double logcorrectionnumnext = 0;
    long double logcorrectiondenomnext = 0;
    MathVector prevLayer = MathVector(numNodes);
    MathVector currLayerOld = MathVector(numNodes);
    MathVector currLayerNew = MathVector(numNodes);
    MathVector nextLayer = MathVector(numNodes);
    if (layer > 0)
    {
        for (size_t i = 0; i < numNodes; ++i)
        {
            prevLayer(i) = G(i,layer-1);
        }
    }
    if (layer < numLayers-1)
    {
        for (size_t i = 0; i < numNodes; ++i)
        {
            nextLayer(i) = G(i,layer+1);
        }
    }
    for (size_t i = 0; i < numNodes; ++i)
    {
        if (i == node)
        {
            currLayerOld(i) = G(i,layer);
            currLayerNew(i) = -1.0*(G(i,layer)-1.0);
        }
        else
        {
            currLayerOld(i) = G(i,layer);
            currLayerNew(i) = G(i,layer);
        }
    }
    if (layer > 0)
    {
        logcorrectionnumprev = logNewCorrection(prevLayer, currLayerNew, numGroups, corrections);
        logcorrectiondenomprev = logNewCorrection(prevLayer, currLayerOld, numGroups, corrections);
    }
    if (layer < numLayers-1)
    {
        logcorrectionnumnext = logNewCorrection(currLayerNew, nextLayer, numGroups, corrections);
        logcorrectiondenomnext = logNewCorrection(currLayerOld, nextLayer, numGroups, corrections);
    }
    long double logcorrection = 0.0;
    if ((layer > 0) && (layer < numLayers-1))
    {
        logcorrection = (logcorrectionnumprev-logcorrectiondenomprev) + (logcorrectionnumnext-logcorrectiondenomnext);
    }
    else if (layer <= 0)
    {
        logcorrection = (logcorrectionnumnext-logcorrectiondenomnext);
    }
    else
    {
        logcorrection = (logcorrectionnumprev-logcorrectiondenomprev);
    }
    return logcorrection;
}

std::vector<std::complex<double>> coeffs(int n)
/*
 Returns values used in the computation of J(k_1,k_2).
 */
{
    std::vector<std::complex<double>> out(n+1,{0,0});
    std::complex<double> last = {1,0};
    for (size_t k = 0; k <= n-1; ++k)
    {
        std::complex<double> sum = {0,0};
        for (size_t r = 1; r <= n; ++r)
        {
            std::complex<double> temp = {1,0};
            for (size_t m = 1; m <= n; ++m)
            {
                if (m != r)
                {
                    temp /= std::polar(1.0, 2.0 * M_PI * r / (n+1)) - std::polar(1.0, 2.0 * M_PI * m / (n+1));
                }
            }
            sum += 1.0 * temp * std::polar(1.0,2.0 * M_PI * r * k / (n+1)) * (std::log(1.0 - std::polar(1.0,2.0 * M_PI * r  / (n+1))) - std::log(-1.0 * std::polar(1.0, 2.0 * M_PI * r  / (n+1))));
        }
        out[k] = sum;
        last -= out[k];
    }
    out[n] = last;
    return out;
}

std::vector<std::vector<double>> corrections(int n)
/*
 Returns the values of J(k_1,k_2) for all 0 <= k_1 < k_2; in particular, the value of J(k_1,k_2) for a given k_1, k_2 is out[k_2][k_1]
 */
{
    std::vector<double> outTemp(n+1,0);
    std::vector<std::vector<double>> out(n+1,outTemp);
    for (size_t k = 1; k <= n; ++k)
    {
        bool stoppingCond = false;
        std::vector<std::complex<double>> coeffsN = coeffs(k);
        for (size_t a = 0; a < k + 1; ++a)
        {
            if (!stoppingCond)
            {
                out[k][a] = std::log(coeffsN[a].real());
            }
            else
            {
                out[k][a] = out[k][a-1] - std::log(1.0);
            }
            stoppingCond = (out[k][a] < -16.0);
        }
    }
    return out;
}

long double logprobrationew(MathMatrix A, std::vector<MathMatrix> G, std::vector<MathMatrix> Gtemp, std::vector<std::vector<double>> corrections)
/*
 Computes a term necessary in the computation of P(A|Gtemp,k)/P(A|G,k)
 A: adjacency matrix for given layer
 G: old group assignments (for all node-layers)
 Gtemp: new group assignments (for all node-layers)
 corrections: the values of J(k_1,k_2) for all 0 <= k_1 < k_2
 */
{
    long double out = 0.0;
    size_t numGroups = G.size();
    size_t numNodes = G[0].getRowSize();
    size_t numLayers = G[0].getColSize();
    for (size_t layer = 0; layer < numLayers; ++layer)
    {
        MathVector oldt = MathVector(numGroups);
        MathVector newt = MathVector(numGroups);
        MathVector oldm = MathVector(numGroups);
        MathVector newm = MathVector(numGroups);
        for (size_t i = 0; i < numNodes; ++i)
        {
            for (size_t j = 0; j < i; ++j)
            {
                bool stoppingCond = false;
                size_t g = numGroups - 1;
                while (not(stoppingCond))
                {
                    if ((G[g](i,layer) > 0.9) && (G[g](j,layer) > 0.9))
                    {
                        oldt(g) = oldt(g) + 1.0;
                        if (A(i,j) > 0.9)
                        {
                            oldm(g) = oldm(g) + 1.0;
                        }
                        stoppingCond = true;
                    }
                    else
                    {
                        g--;
                    }
                }
            }
        }
        for (size_t i = 0; i < numNodes; ++i)
        {
            for (size_t j = 0; j < i; ++j)
            {
                bool stoppingCond = false;
                size_t g = numGroups - 1;
                while (not(stoppingCond))
                {
                    if ((Gtemp[g](i,layer) > 0.9) && (Gtemp[g](j,layer) > 0.9))
                    {
                        newt(g) = newt(g) + 1.0;
                        if (A(i,j) > 0.9)
                        {
                            newm(g) = newm(g) + 1.0;
                        }
                        stoppingCond = true;
                    }
                    else
                    {
                        g--;
                    }
                }
            }
        }
        for (size_t g = 0; g < numGroups; ++g)
        {
            out += (logfactorial2(newm(g)) + logfactorial2(newt(g)-newm(g))) - logfactorial2(newt(g)+1);
            out -= (logfactorial2(oldm(g)) + logfactorial2(oldt(g)-oldm(g))) - logfactorial2(oldt(g)+1);
        }
    }
    for (size_t g = 1; g < numGroups; ++g)
    {
        for (size_t layer = 1; layer < numLayers; ++layer)
        {
            long double logcorrectionnumprev = 0;
            long double logcorrectiondenomprev = 0;
            long double logcorrectionnumnext = 0;
            long double logcorrectiondenomnext = 0;
            MathVector prevLayerOld = MathVector(numNodes);
            MathVector prevLayerNew = MathVector(numNodes);
            MathVector currLayerOld = MathVector(numNodes);
            MathVector currLayerNew = MathVector(numNodes);
            for (size_t i = 0; i < numNodes; ++i)
            {
                prevLayerOld(i) = G[g](i,layer-1);
                prevLayerNew(i) = Gtemp[g](i,layer-1);
            }
            for (size_t i = 0; i < numNodes; ++i)
            {
                
                currLayerOld(i) = G[g](i,layer);
                currLayerNew(i) = Gtemp[g](i,layer);
            }
            logcorrectionnumprev = logNewCorrection(prevLayerNew, currLayerNew, numGroups, corrections);
            logcorrectiondenomprev = logNewCorrection(prevLayerOld, currLayerOld, numGroups, corrections);
            out += (logcorrectionnumprev-logcorrectiondenomprev);
        }
    }
    return out;
}

double logfactorial2(double n)
//Computes log(n!) via log-gamma function
{
    long double out = 0.0;
    if (n > 1)
    {
        out = lgamma(1.0 + (long double)n);
    }
    return out;
}
