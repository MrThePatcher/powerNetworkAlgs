import networkx as nx
import matplotlib.pyplot as plt
import random
import string
import timeit

import numpy
import numpy as np
from math import sqrt, pow
from scipy.stats import poisson


def fillArray(fields):
    for n in range(len(fields)):
        fields[n][0] = n
        fields[n][1] = random.uniform(0, 1)
        fields[n][2] = random.uniform(0, 1)


# Creates an Array with independent fields, with r rows and c columns
def createArray(c, r):
    array = [[0] * c for i in range(r)]
    return array


def calculateDistance(nodes):
    # distanceArray = createArray(3, 0)
    distanceArray = []
    for i in range(len(nodes) - 1):
        nodedistanceList = createArray(3, 0)
        nodeV_x = nodes[i][1]
        nodeV_y = nodes[i][2]
        for j in range(len(nodes)):
            # adding i to j should remove the double counting from A-B and B-a C-D D-C (wrong)
            # changing it to while smaller than length of nodes is necessary
            nodeW_x = nodes[j][1]
            nodeW_y = nodes[j][2]
            # if i != j:
            # d(v,w) = square-root of (x_v - x_w)^2 + (y_v - y_w)^2
            # distance = sqrt(pow(nodeV_x - nodeW_x, 2) + pow(nodeV_y - nodeW_y, 2))
            # toAppend = [nodes[i][0], nodes[j][0], distance]
            # distanceArray.append(toAppend)
            if i < j:
                # d(v,w) = square-root of (x_v - x_w)^2 + (y_v - y_w)^2
                distance = sqrt(pow(nodeV_x - nodeW_x, 2) + pow(nodeV_y - nodeW_y, 2))
                toAppend = [nodes[i][0], nodes[j][0], distance]
                nodedistanceList.append(toAppend)
        distanceArray.append(nodedistanceList)
    return distanceArray


def filterDistances(d_min, d_max, distanceList):
    outputlist = []
    for nodelist in distanceList:
        appendableNodeList = createArray(3, 0)
        for i in range(len(nodelist)):
            dist = nodelist[i][2]
            if dist < d_min or dist > d_max:
                continue
            else:
                appendableNodeList.append(nodelist[i])
        outputlist.append(appendableNodeList)
    return outputlist


# returns new graph with edges?
# problem of edge selection: doestn make sence have to rethink
def edging(k, distancelist, graph):
    for nodelist in distancelist:
        # draw new poisson distribution sample
        r = random.randint(1, len(nodelist))
        data_binom = poisson.rvs(mu=2.67, size=r)
        # creates list like [2 3 2 2 0]
        for elem in data_binom:
            node_v = nodelist[elem][0]
            node_w = nodelist[elem][1]
            graph.add_edge(nodelist[elem][0], nodelist[elem][1])
    return


def paperedging(distancelist, graph):
    for nodelist in distancelist:
        if nodelist:
            k = np.random.exponential(2.67)
            kInt = int(k)
            uniformList = np.random.uniform(0, len(nodelist), kInt)
            for i in uniformList:
                indexi = int(i)
                graph.add_edge(nodelist[indexi][0], nodelist[indexi][1], length=nodelist[indexi][2])
    return


def generateOneGraphPaper(numNodes, dMin, dMax):
    G = nx.Graph()
    Nodes = createArray(3, numNodes)
    fillArray(Nodes)
    # print(Nodes)
    # Add Nodes To Graph
    for [name, x, y] in Nodes:
        G.add_node(name, pos=(x, y))

    # This creates Plot
    # nx.draw(G, with_labels=True, node_size=200)
    # nx.draw(G, pos=G.nodes('pos'), with_labels=True, node_size=50)
    # This makes created Plot visible
    # plt.savefig('nodes.png')
    # plt.show()
    distances = calculateDistance(Nodes)
    # created Array: Array[i] is list of distances from i to all nodes in range
    distances = filterDistances(dMin, dMax, distances)
    paperedging(distances, G)
    return G


# List of Nodes
# Create Array with first field is name and second field is X and third field is Y second number defines number of Nodes

def fG4Me(numNodes, dmin, dmax):
    initGraph = generateOneGraphPaper(numNodes, dmin, dmax)
    if nx.is_connected(initGraph):
        print("instant Graph!")
    else:
        generations = 0
        while not nx.is_connected(initGraph):
            generations += 1
            if generations % 100 == 0:
                print("Try number " + str(generations))
            initGraph = generateOneGraphPaper(numNodes, dmin, dmax)
        print("only took " + str(generations) + " tries !")
    return initGraph


def wrapper(func, *args, **kwargs):
    def wrapped():
        return func(*args, **kwargs)

    return wrapped


def gridLineImpedances(graph, z_0, eps):
    graph = nx.Graph()
    # Z_0 =r_jx line impedance
    # IEEE30
    impCalcNumber = random.uniform(-eps, eps)
    impedances = {}
    for (u, v, d) in graph.edges.data():
        key = (u, v)
        imp = z_0 * d + impCalcNumber;
        impedances[key]=imp
    return impedances
def perUnitImpedances(impedances={}):
    #S_B = 100MVA
    powBaseS_B=100
    #V_B = XkV
    voltageRatiov_b=1
    puImpedances={}
    for key in impedances:
        puImpedances[key]=(impedances[key]*powBaseS_B)/pow(voltageRatiov_b,2)
    return puImpedances
def connectivityMatrix(graph):

    admMatA=numpy.zeros((len(graph.edges,len(graph.nodes))))
    index1=0 # is the index which line we are looking at
    for (u,v,d) in graph.edges.data():
        admMatA[index1][u]=1
        admMatA[index1][v] =-1
    return admMatA
def networkAdmittance(graph,epsYMax):
    graph = nx.Graph()
    impedances=gridLineImpedances(graph,1,1)
    puImpedances=perUnitImpedances(impedances)
    lineAdmittance=[]
    for imp in impedances.values():
        lineAdmittance.append(1/imp)
    diagAdm=numpy.diag(lineAdmittance)
    epsY=[]
    for n in graph.nodes:
        epsY.append(random.uniform(0,epsYMax))
    diagEps=numpy.diag(epsY)
    A=connectivityMatrix(graph)
    Y_bus=A.T()@diagAdm@A+diagEps
    return Y_bus

def assignLoadP(pMin,pMax):
    return random.uniform(pMin,pMax)
def assignLoadQ(qMin,qMax):
    return random.uniform(qMin,qMax)
def assignInertiaM(mMin,mMax):
    return random.uniform(mMin, mMax)
def assignImpedanceX(xMin,xMax):
    return random.uniform(xMin, xMax)
def assignVoltagePotentialE(eMin,eMax):
    return random.uniform(eMin, eMax)

PaperGraph=fG4Me(30, 0, 0.3)

nx.draw(PaperGraph, pos=PaperGraph.nodes('pos'), with_labels=True, node_size=200)
#  plt.savefig('nodes.png')
plt.show()
sumDegree = 0
for (node, degree) in PaperGraph.degree:
    sumDegree += degree
avgDegree = sumDegree / len(PaperGraph.nodes)
print("Average degree: "+str(avgDegree))
print(PaperGraph.adj)



# wrapped = wrapper(fG4Me,40,0,0.25)
# for i in range(15):
#    print(timeit.timeit(wrapped,number=1))

# mu = lambda should be 2.67 from paper, size = k debateble
