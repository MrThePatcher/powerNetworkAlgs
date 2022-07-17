import networkx as nx
import matplotlib.pyplot as plt
import random
import string
import timeit

import numpy
import numpy as np
from math import sqrt, pow
from scipy.stats import poisson

class Node:
    def __init__(self,name,x,y):
        self.name=name
        self.x=x
        self.y=y
    def getX(self):
        return self.x
    def getY(self):
        return self.y
    def getName(self):
        return self.name

def createArray(c, r):
    array = [[0] * c for i in range(r)]
    return array

# eigenvalues distribution of adjazenz matrix
#
# RT nested networks
# lcoal neigborhood distance threshhold
# k element of geometric distribution
#
# select potential neigborhood nodes:
# 	d_0= distance threshhold
# 	i,j Node numbers
# 	|j-i| < d_0
# 	N_d_0= potential close nodes
# select k links
#
# then rewire
def generateRTNested(numNodes, dMin, dMax):
    # generate array with all nodes
    G=nx.Graph()
    nodesArray=[]
    d_0=4
    for i in range(numNodes):
        #nodesArrayUnnoetig?
        nodesArray.append(Node(i,0,0))
        G.add_node(i)
    # generate adjacency array
    adj=createArray(numNodes,numNodes)
    for i_n in range(len(nodesArray)):
        k=numpy.random.geometric(0.63)
        locNbh=localNeighborhoodSelection(k,d_0,i_n,nodesArray)
    ### WIRING
        for i in range(len(locNbh)):
            adj[i_n][locNbh[i]]=1
            adj[locNbh[i]][i_n] = 1
            if i_n<locNbh[i]:
                G.add_edge(i_n,locNbh[i])
            else:
                G.add_edge(locNbh[i],i_n)
    components=[]
    pos= nx.shell_layout(G)
    nx.draw(G, pos=pos, with_labels=True, node_size=200)
    # nx.draw(G, with_labels=True, node_size=200)
    plt.show()
    print(nx.number_connected_components(G))
    for compnodes in nx.connected_components(G):
        print(compnodes)
    ### REWIRING
    markovSample=markovChain(0.4,0.8,G.number_of_nodes())
    rewireEdges=rewiringv3(markovSample[0],G)
    G.add_edges_from(rewireEdges[0])
    G.remove_edges_from(rewireEdges[1])
    nx.draw(G, pos=pos,with_labels=True, node_size=200)
    #nx.draw(G, with_labels=True, node_size=200)
    plt.show()
    ###LATTICE CONNECTIONS
    latticeConns(G)
    #nx.draw(G, with_labels=True, node_size=200)
    nx.draw(G, pos=pos, with_labels=True, node_size=200)
    plt.show()
    # generate incidence array
    sumDegree = 0
    for (node, degree) in G.degree:
        sumDegree += degree
    avgDegree = sumDegree / len(G.nodes)
    print("Average degree: " + str(avgDegree))
    nx.draw(G, with_labels=True, node_size=200)
    plt.show()
    return None

# select k edges at random from local neighborhood of n_i
#returns List of NodeNumbers to connect with
def localNeighborhoodSelection(k,d,n_i,nArr):
    minNbh=0
    if n_i-d >= 0:
        minNbh=n_i-d
    maxNbh=len(nArr)-1
    if n_i+d <len(nArr):
        maxNbh=n_i+d
    locNbh=[]
    for i in nArr[minNbh:maxNbh]:
        if not i.getName()==n_i:
            locNbh.append(i.getName())

    try:
        edgeSample=random.sample(locNbh,k)
    except ValueError:
        edgeSample=random.sample(locNbh,len(locNbh))
    # selectedEdges=[]
    # # Check for existence maybe use try w
    # if edgeSample:
    #     for i in range(len(edgeSample)):
    #         if locNbh[i]:
    #             selectedEdges.append(locNbh[i])
    return edgeSample




# beta = 1/averageClusersize
# p1= prob( nodes having rewire links)
# alpha = beta*p1/1-p1
# prw= rewiriring probability
##INPUT Alpha Beta and number of Nodes
##OUTPUT: LIST [ MarkovChainList, NodeNumbersOfRewires, NodeNumberOfNotRewires]
def markovChain(alpha,beta,numNodes):
    states = [0,1]
    transitionName=[['00','01'],['10','11']]
    transitionMatrix=[[1-alpha,alpha],[beta,1-beta]]
    #starting State is 0
    currState= 0
    sampleList=[]
    nodeRewireList=[]
    nodeNotRewireList=[]
    prob = 1
    for i in range(numNodes):
        if currState == 0:
            nextState = np.random.choice(transitionName[0],replace=True,p=transitionMatrix[0])
            if nextState=='00':
                prob=prob*transitionMatrix[0][0]
                sampleList.append(0)
                nodeNotRewireList.append(i)
            elif nextState=='01':
                prob = prob * transitionMatrix[0][1]
                currState=1
                sampleList.append(1)
                nodeRewireList.append(i)
            else:
                print("Hier sollt ich net sein")
        elif currState == 1:
            nextState = np.random.choice(transitionName[1],replace=True,p=transitionMatrix[1])
            if nextState=='10':
                prob=prob*transitionMatrix[1][0]
                currState=0
                sampleList.append(0)
                nodeNotRewireList.append(i)
            elif nextState=='11':
                prob = prob * transitionMatrix[1][1]
                sampleList.append(1)
                nodeRewireList.append(i)
            else:
                print("Hier sollt ich auch net sein")
    #returns the list of indices of nodes that have state one
    #return all lists
    return [nodeRewireList,nodeNotRewireList,sampleList]

#gets list
#rewires each 1 node to another 1 node
def rewiring(nodeRewireList):

    rewiringEdges=[]
    for i in nodeRewireList:
        numRewireEdges=random.uniform(1,len(nodeRewireList))
        numRewireEdges=1
        rngElem=random.choice(nodeRewireList)
        if i<rngElem:
            rewiringEdges.append((i,rngElem))
        elif rngElem<i:
            rewiringEdges.append((rngElem, i))
        else:
            pass
    return rewiringEdges
#rewires each 1 node to another 1 node
def rewiringv2(nodeRewireList,notRewireList,G):

    rewiringEdges=[]
    for i in nodeRewireList:
        v_edges=G.edges(i)
        for (u,v) in v_edges:
            choice= random.choices([0,1],[0.43,0.57],1)
            if choice==1:
                G.remove_edge((u,v))
                rngNode=random.choice(nodeRewireList)
                G.add_edge(u,rngNode)
            else:
                pass
    return None

    #takes markovchain and returns which edges need to be added and which edges need to be removed
##INPUT NodeNumbersOfRewires, Graph to Rewire
##OUTPUT: LIST [ EdgesAddedThroughRewirin, EdgesRemovedThroughRewiring]
def rewiringv3(nodeRewireList, G):
    #create list of 1clusters
    oneClusters=[]
    i = 0
    while i < len(nodeRewireList):
        number=G.number_of_nodes() - i
        clusterlist = []
        j = 0
        while j < number:
            if j==0:
                clusterlist.append(nodeRewireList[i+j])
                j=j+1
            elif j+i<len(nodeRewireList):
                index=i+j
                #if element now is one bigger than element before
                if (nodeRewireList[index]-nodeRewireList[index-1])==1:
                    clusterlist.append(nodeRewireList[index])
                    j=j+1
                else:
                    i=i+j
                    break
            else:
                i=i+1
                break
        oneClusters.append(clusterlist)
    #if first element in first cluster and last element in last cluster are 0 and 12 thann connect this cluster
    if oneClusters[0][0]==0 :
        lastCluster=len(oneClusters)-1
        lastElem=len(oneClusters[lastCluster])-1
        if oneClusters[lastCluster][lastElem]==(G.number_of_nodes() - 1):
            for elem in oneClusters.pop():
                oneClusters[0].append(elem)
    rewiringEdges = []
    rewireRemEdges = []
    for i in oneClusters:
        v_edges = G.edges(i)
        for (u, v) in v_edges:
            choice = random.choices([0, 1], weights=[0.43, 0.57],k=1)
            if choice[0] == 1:
                rewireRemEdges.append((u,v))
                otherClusters=[x for x in oneClusters if x != i]
                rngCluster = random.choice(otherClusters)
                rngNode= random.choice(rngCluster)
                rewiringEdges.append((u,rngNode))
            else:
                pass
    return [rewiringEdges,rewireRemEdges]
#These Lines make no Sens here right now
        # numRewireEdges=random.uniform(1,len(nodeRewireList))
        # numRewireEdges=1
        # rngElem=random.choice(nodeRewireList)
        # if i<rngElem:
        #     rewiringEdges.append((i,rngElem))
        # elif rngElem<i:
        #     rewiringEdges.append((rngElem, i))
        # else:
        #     pass


#testComponentsBeforehand
def latticeConns(G):
    if nx.number_connected_components(G)== 1:
        print("Nothing to LatticeConnect")
        return None
    else:
    #number of lattice connections around <k>
        numComponents = nx.number_connected_components(G)
        components=[]
        for compnodes in nx.connected_components(G):
            components.append(compnodes)
            print(compnodes)
        numLatConns = int(random.gauss(2.67, 1))

        for i in range(numComponents):

            j=1
            while j+i < numComponents:
                allEdges = []
                for k in components[i]:
                    for l in components[i+j]:
                        #print("another")
                        if k<l:
                            allEdges.append((k,l))
                        elif l<k:
                            allEdges.append((l,k))
                        else:
                            print("k = l ???")
                #all edges is filled with all edges between component i and component i+1
                sampleLatEdges=random.sample(allEdges,numLatConns)
                for (v,u) in sampleLatEdges:
                    G.add_edge(v,u)
                j=j+1
        return None

def testLatticeNoArg():
    G = nx.Graph()
    numNodes=13
    nodesArray = []
    d_0 = 4
    for i in range(numNodes):
        # nodesArrayUnnoetig?
        nodesArray.append(Node(i, 0, 0))
        G.add_node(i)
    G.add_edge(0,1)
    G.add_edge(1, 3)
    G.add_edge(1, 2)
    g.add_edge(2,4)

    G.add_edge(5,6)
    G.add_edge(6,7)
    G.add_edge(7, 9)

    G.add_edge(10, 11)
    G.add_edge(10,12)
    G.add_edge(8,10)
    G.add_edge(11,12)
    print("edges in 1 adn 2")
    print(G.edges([1,2]))
    print("edges on 1")
    print(G.edges([1]))
    pos = nx.shell_layout(G)
    nx.draw(G, pos=pos,with_labels=True, node_size=200)
    # nx.draw(G, with_labels=True, node_size=200)
    plt.show()
    latticeConns(G)
    # nx.draw(G, with_labels=True, node_size=200)
    nx.draw(G, pos=pos, with_labels=True, node_size=200)
    plt.show()
    return None

def testLattice(G):
    pos = nx.shell_layout(G)
    nx.draw(G, pos=pos, with_labels=True, node_size=200)
    # nx.draw(G, with_labels=True, node_size=200)
    plt.show()
    latticeConns(G)
    # nx.draw(G, with_labels=True, node_size=200)
    nx.draw(G, pos=pos, with_labels=True, node_size=200)
    plt.show()
    return None

def testRewire():
    G = nx.Graph()
    numNodes=13
    nodesArray = []
    d_0 = 4
    for i in range(numNodes):
        # nodesArrayUnnoetig?
        nodesArray.append(Node(i, 0, 0))
        G.add_node(i)
    print("Number of Nodes = "+ str(G.number_of_nodes()))
    G.add_edge(0, 1)
    G.add_edge(1, 3)
    G.add_edge(1, 2)
    G.add_edge(2, 4)

    G.add_edge(5, 6)
    G.add_edge(6, 7)
    G.add_edge(7, 9)

    G.add_edge(10, 11)
    G.add_edge(10, 12)
    G.add_edge(8, 10)
    G.add_edge(11, 12)
    pos = nx.shell_layout(G)
    nx.draw(G, pos=pos,with_labels=True, node_size=200)
    # nx.draw(G, with_labels=True, node_size=200)
    plt.show()
    markovSample = markovChain(0.4, 0.8, G.number_of_nodes())
    rewireEdges = rewiringv3(markovSample[0],G)
    G.add_edges_from(rewireEdges[0])
    G.remove_edges_from(rewireEdges[1])
    # nx.draw(G, with_labels=True, node_size=200)
    nx.draw(G, pos=pos, with_labels=True, node_size=200)
    plt.show()
    testLattice(G)
    return None
#testRewire()
generateRTNested(20,0,0)
