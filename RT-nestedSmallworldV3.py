import networkx as nx
import matplotlib.pyplot as plt
import random
import string
import timeit

import numpy
import numpy as np
from math import sqrt, pow
from scipy.spatial import distance
from scipy.stats import poisson


def createArray(c, r):
    array = [[0] * c for i in range(r)]
    return array

def assignLatEdgeAdm(admList,latEdges,graph):
    if(latEdges):
        latEdges.reverse()
        for (u, v) in latEdges:
            elem=admList.pop()
            graph.add_edge(u, v, adm=elem)
    return graph

def assignRewEdgeAdm(admList,rewEdges,graph):
    if(rewEdges):
        rewEdges.reverse()
        for (u, v) in rewEdges:
            elem=admList.pop()
            graph.add_edge(u,v,adm=elem)
    return graph
def getNormalEdges(graph):
    normEdges=[]
    for (u,v) in graph.edges():
        normEdges.append((u,v))
    return normEdges


def generateRTNested(numNodes, dMin,pGeom):
    # generate array with all nodes
    G=nx.Graph()
    d_0=dMin
    for i in range(numNodes):
        G.add_node(i)
    # generate adjacency array
    adj=createArray(numNodes,numNodes)

    for i_n in range(G.number_of_nodes()):

        k=numpy.random.geometric(pGeom)
        # print('K is now '+str(k))
        locNbh=localNeighborhoodSelection(k,d_0,i_n,G)
    ### WIRING
        for i in range(len(locNbh)):
            adj[i_n][locNbh[i]]=1
            adj[locNbh[i]][i_n] = 1
            if i_n<locNbh[i]:
                G.add_edge(i_n,locNbh[i])
            else:
                G.add_edge(locNbh[i],i_n)
    ### REWIRING
    markovSample=markovChain(0.6,0.8,G.number_of_nodes())
    rewireEdges=rewiringv5(markovSample[0],G)
    G.remove_edges_from(rewireEdges[1])
    normalEdges=getNormalEdges(G)
    G.add_edges_from(rewireEdges[0])
    ###LATTICE CONNECTIONS
    latEdges=latticeConns(G)
    pos = nx.shell_layout(G)
    grAdm = impedanceCalcv2(G,pos)
    ###Assigning right admittancse
    G=assignLatEdgeAdm(grAdm,latEdges,G)
    G=assignRewEdgeAdm(grAdm,rewireEdges[0],G)
    G = assignLatEdgeAdm(grAdm, normalEdges, G)
    grAdm=[]
    for (u,v,d) in G.edges.data():
        grAdm.append(d['adm'])

    ybusPaper = networkAdmittance(G, 1, grAdm)
    # nx.draw(G, with_labels=True, node_size=200)
    # plt.show()
    return [G,ybusPaper,pos]


def calcAvgDegree(graph):
    sumDegree = 0
    for (node, degree) in graph.degree:
        sumDegree += degree
    avgDegree = sumDegree / len(graph.nodes)
    return avgDegree
# select k edges at random from local neighborhood of n_i
#returns List of NodeNumbers to connect with
def localNeighborhoodSelection(k,d,n_i,Graph):
    nArr=list(Graph.nodes)
    minNbh=0
    if n_i-d >= 0:
        minNbh=n_i-d
    maxNbh=Graph.number_of_nodes()-1
    if n_i+d <Graph.number_of_nodes():
        maxNbh=n_i+d
    locNbh=[]
    #determines the nodes that are in the local neighborhood
    for i in nArr[minNbh:maxNbh]:
        if not i==n_i:
            locNbh.append(i)
    #selects a random sample of k nodes, that are wired to node n_i
    try:
        edgeSample=random.sample(locNbh,k)
    except ValueError:
        edgeSample=random.sample(locNbh,len(locNbh))
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
def rewiringv5(nodeRewireList, G):
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
                ##we approached the last element
                i=i+j
                break
        oneClusters.append(clusterlist)
    #if first element in first cluster and last element in last cluster are 0 and 12 thann connect this cluster
    if oneClusters[0][0]==0 :
        lastCluster=len(oneClusters)-1
        lastElem=len(oneClusters[lastCluster])-1
        if oneClusters[lastCluster][lastElem]==(G.number_of_nodes() - 1):
            for elem in oneClusters.pop():
                oneClusters[0].append(elem)
    rewiringEdges = {}
    rewireRemEdges = {}
    ##choose the rewiring edges
    for i in oneClusters:
        v_edges = G.edges(i)
        for (u, v) in v_edges:
            choice = random.choices([0, 1], weights=[0.43, 0.57],k=1)
            if choice[0] == 1:
                    if (v,u) in rewireRemEdges.values() or (u,v) in rewireRemEdges.values():
                        pass
                    else:
                        ## die Kante kann entfernt werden
                        otherClusters = [x for x in oneClusters if x != i]
                        rngCluster = random.choice(otherClusters)
                        rngNode = random.choice(rngCluster)
                        if(u,rngNode) in rewiringEdges.values() or (rngNode,u) in rewiringEdges.values():
                            j=0
                            while j<10:
                                rngCluster = random.choice(otherClusters)
                                rngNode = random.choice(rngCluster)
                                if not((u,rngNode) in rewiringEdges.values() or (rngNode,u) in rewiringEdges.values()):
                                    break
                                j=j+1
                            if j<10:
                                rewiringEdges[u, rngNode] = (u, rngNode)
                                rewireRemEdges[(v, u)] = (v, u)
                        else:
                            rewiringEdges[u,rngNode]=(u,rngNode)
                            rewireRemEdges[(v, u)] = (v, u)
            else:
                pass
    return [list(rewiringEdges),list(rewireRemEdges)]


#testComponentsBeforehand
def latticeConns(G):
    if nx.number_connected_components(G)== 1:
        print("Nothing to LatticeConnect")
        return None
    else:
        latEdges=[]
    #number of lattice connections around <k>
        numComponents = nx.number_connected_components(G)
        components=[]
        for compnodes in nx.connected_components(G):
            components.append(compnodes)
        numLatConns=0
        while numLatConns==0:
            numLatConns = int(random.gauss(2.67, 1))
        for i in range(numComponents):
            j=1
            while i+j < numComponents:
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
                try:
                    sampleLatEdges=random.sample(allEdges,numLatConns)
                except ValueError:
                    sampleLatEdges=random.sample(allEdges,len(allEdges))
                for (u,v) in sampleLatEdges:
                    latEdges.append((u,v))
                    G.add_edge(u,v)
                j=j+1
        return latEdges
def impedanceCalcv2(G,pos):
    #a=2.1468 b=0.10191 Thet = 1/b
    sample=numpy.random.gamma(2.1468,9.812579,G.number_of_edges()+2)
    y_list=[]
    # mean and deviation
    #x=np.random.normal(0.1999756097560976,1)
    #r=np.random.normal(0.07569756097560976,1)
    r = 1
    x = 0.1999756097560976
    # S_B = 100MVA
    powBaseS_B = 100
    # V_B = XkV
    voltageRatiov_b = 1
    sample.sort()
    its=0
    edglst=list(G.edges)
    for z in sample:
        try:
            edge=edglst[its]
        except IndexError:
            edge=random.choice(edglst)
        posu=pos[edge[0]]
        posv=pos[edge[1]]
        lgth=distance.euclidean(posu,posv)
        img=(x/r)*z
        imp=z+img*1j
        impCalcNumber = random.uniform(-1, 1)
        imlgth=imp+lgth+impCalcNumber
        puimp= (imlgth*powBaseS_B)/pow(voltageRatiov_b,2)
        y=1/puimp
        y_list.append(y)
        its=its+1
    # y_list.sort()
    return y_list

def connectivityMatrix(graph):
    admMatA=numpy.zeros((len(graph.edges),len(graph.nodes)))
    admMatA = numpy.zeros((graph.number_of_edges(),graph.number_of_nodes()))
    index1=0 # is the index which line we are looking at
    for (u,v,d) in graph.edges.data():
        admMatA[index1][u]=1
        admMatA[index1][v] =-1
        index1=index1+1
    return admMatA
def networkAdmittance(graph,epsYMax,y_list):
    #removed code to calc impedances see this method in randomMatrices
    #y_list=[]
    # for (u, v, d) in graph.edges.data():
    #     y_list.append(d)
    diagAdm=numpy.diag(y_list)
    epsY=[]
    for n in graph.nodes:
        epsY.append(random.uniform(0,epsYMax))
    diagEps=numpy.diag(epsY)
    A=connectivityMatrix(graph)
    Y_bus=A.T@diagAdm@A+diagEps
    return Y_bus

def oneRunWithPrints():
    ###RTNested
    retlist=generateRTNested(30,2,0.37)
    rtgraph=retlist[0]
    ybusPaper=retlist[1]

    ###Eigenvalues Calculation
    eigenvalues,eigenvectors=numpy.linalg.eig(ybusPaper)
    invEigVec=np.linalg.inv(eigenvectors)
    diagEV=invEigVec.dot(ybusPaper).dot(eigenvectors)
    diagValues=numpy.diag(diagEV)
    diagV2=diagValues*(numpy.identity(rtgraph.number_of_nodes()))
    ### Topological Chars
    avgcc=nx.average_clustering(rtgraph)
    dpcc=nx.degree_pearson_correlation_coefficient(rtgraph)
    avgl=nx.average_shortest_path_length(rtgraph)
    cctext='<Cc>:'+str(avgcc)
    dpcctext='pearson Coefficient:'+str(dpcc)
    avgltext='<L>'+str(avgl)
    avgDegree=calcAvgDegree(rtgraph)
    avgtext="Average degree: " + str(avgDegree)
    print(cctext)
    print(dpcctext)
    print(avgltext)
    ##printPlots
    fig, (ax1, ax2,ax3) = plt.subplots(1, 3, figsize=(15, 5), tight_layout=True)
    x=[ele.real for ele in eigenvalues]
    y=[ele.imag for ele in eigenvalues]
    ax1.hist(y,linewidth=9,edgecolor="white")
    nx.draw(rtgraph, pos=retlist[2], with_labels=True, node_size=200,ax=ax2)
    ax2.text(0.3,1,avgtext,bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 3},fontsize='small')
    ax2.text(0.3,0.90,cctext,bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 3},fontsize='small')
    ax2.text(0.3,0.80,dpcctext,bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 3},fontsize='small')
    ax2.text(0.3,0.70,avgltext,bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 3},fontsize='small')
    ax3.scatter(x,y)
    plt.show()
    return None
def runGraphforTopological(numNodes,tries):
    sumAvgDegree=0
    sumAvgcc=0
    sumpcc=0
    sumpathLength=0
    sumEdges=0
    for i in range(tries):
        foo=generateRTNested(numNodes, 2, 0.37)
        randTopgraph = foo[0]
        avgDegree = calcAvgDegree(randTopgraph)
        avgcc = nx.average_clustering(randTopgraph)
        dpcc = nx.degree_pearson_correlation_coefficient(randTopgraph)
        avgl = nx.average_shortest_path_length(randTopgraph)
        sumAvgDegree = sumAvgDegree+avgDegree
        sumAvgcc = sumAvgcc+avgcc
        sumpcc = sumpcc+dpcc
        sumpathLength = sumpathLength+avgl
        sumEdges=sumEdges+randTopgraph.number_of_edges()
    avgDegree = sumAvgDegree/tries
    avgcc = sumAvgcc/tries
    dpcc = sumpcc/tries
    avgl = sumpathLength/tries
    avgEdges=sumEdges/tries
    print("avgavgDELTA: "+str(avgDegree))
    print("avgavgPATHL: " + str(avgl))
    print("avgavgPEARSON: " + str(dpcc))
    print("avgavgCLUSTER: " + str(avgcc))
    print("avgEdges: " + str(avgEdges))
    return None
oneRunWithPrints()


