#### Script adapted from https://github.com/ratschlab/scim

import numpy as np
import scipy as sp
import pandas as pd
import networkx as nx
import itertools
import time
from scipy.spatial import cKDTree


def get_knn_onesided(source, target, knn_k, knn_n_jobs, st=True):
    """Get kNN connections between source and target
    source: pandas dataframe which rows will be placed at the source of the graph (larger dataset)
    target: pandas dataframe which rows will be placed at the target of the graph (smaller dataset)
    knn_k: k neighbors to be found
    knn_n_jobs: how many processors to use for the job
    st: source to target (True/False)
    """
    knn = cKDTree(target).query(x=source, k=knn_k, n_jobs=knn_n_jobs)
    if st:
        knn_source_idx = ['source_'+str(x) for x in np.array(range(source.shape[0])).repeat(knn_k)]
        knn_target_idx = ['target_'+str(x) for x in knn[1].reshape((-1,1)).flatten()]
    else:
        knn_target_idx = ['target_'+str(x) for x in np.array(range(source.shape[0])).repeat(knn_k)]
        knn_source_idx = ['source_'+str(x) for x in knn[1].reshape((-1,1)).flatten()]
    knn_dist = knn[0].reshape((-1,1)).flatten()

    return np.array(knn_source_idx), np.array(knn_target_idx), np.array(knn_dist)

def get_knn_union(source, target, knn_k, knn_n_jobs):
    """Get kNN connections between source and target (union)
    source: pandas dataframe which rows will be placed at the source of the graph (larger dataset)
    target: pandas dataframe which rows will be placed at the target of the graph (smaller dataset)
    knn_k: k neighbors to be found
    knn_n_jobs: how many processors to use for the job

    output: nodes at the source of the graph, nodes at the target of the graph, Euclidean distance between the nodes
    """
    knn_source_idx, knn_target_idx, knn_dist = get_knn_onesided(source, target, knn_k, knn_n_jobs, st=True)
    knn2_source_idx, knn2_target_idx, knn2_dist = get_knn_onesided(target, source, knn_k, knn_n_jobs, st=False)

    knn_source_idx = np.append(knn_source_idx, knn2_source_idx)
    knn_target_idx = np.append(knn_target_idx, knn2_target_idx)
    knn_dist = np.append(knn_dist, knn2_dist)
    # remove duplicated connections
    knn_df = pd.DataFrame({'source':knn_source_idx, 'target':knn_target_idx, 'dist':knn_dist})
    knn_df = knn_df.drop_duplicates()
    knn_source_idx = knn_df['source'].to_list()
    knn_target_idx = knn_df['target'].to_list()
    knn_dist = np.array(knn_df['dist'].to_list())

    return knn_source_idx, knn_target_idx, knn_dist

def get_null_cost(cost, null_cost_percentile=95):
    """Compute the cost of matching to the null node based on all other costs
    cost: vector of costs (weights on all source->target edges)
    null_cost_percentile: percentile of the cost that should correspond to null node matching
    """
    null_cost = int(np.ceil(np.percentile(cost, null_cost_percentile)))

    return null_cost

def extend_graph_null(G, source_idx, null_cost):
    """Extend the graph by adding a null node on target side (linked to source nodes and sink)
    G: directed graph object
    source_idx: nodes at the source of the graph
    null_cost: null match penalty

    output: directed graph object
    """
    null_capacity = len(source_idx)
    G.add_node('target_null')
    source_null_edges = list(itertools.product(source_idx, ['target_null'], [{'capacity': 1, 'weight':null_cost}]))
    G.add_edges_from(source_null_edges)
    null_sink_edges = list(itertools.product(['target_null'], ['sink'], [{'capacity': null_capacity, 'weight':0}]))
    G.add_edges_from(null_sink_edges)

    return G

def get_target_capacities(source_idx, target_idx, capacity_method='uniform', seed=456):
    """Compute a vector of capacities from target nodes to sink
    source_idx: nodes at the source of the graph
    target_idx: nodes at the target of the graph
    capacity_method: how to set capacities on the target to sink edged {uniform, inf, top, 1to1}
    seed: random seed

    output: vector of capacities (of length=len(target_idx))
    """
    if(capacity_method=='uniform'):
        capacity = int(np.floor(len(source_idx)/len(target_idx)))
        capacities = [capacity]*len(target_idx)
        # randomly distribute remaining cells
        n_remaining = len(source_idx) - np.sum(capacities)
        np.random.seed(seed)
        to_add_idx = np.random.choice(range(len(capacities)),n_remaining, replace=False)
        for i in to_add_idx:
            capacities[i] = capacities[i] + 1
    elif(capacity_method=='inf'):
        capacity = np.inf
        capacities = [capacity]*len(target_idx)
    elif(capacity_method=='top'):
        capacity = len(source_idx) + 1000
        capacities = [capacity]*len(target_idx)
    elif(capacity_method=='1to1'):
        capacity = 1
        capacities = [capacity]*len(target_idx)
    else:
        raise NotImplementedError

    return capacities

def build_graph_base(source_idx, target_idx, capacity_method='uniform', seed=456):
    """Build a graph base
    source_idx: nodes at the source of the graph
    target_idx: nodes at the target of the graph
    method: how to set capacities on the target to sink edged {uniform, inf, top, 1to1}
    seed: random seed

    output: directed graph object
    """
    G = nx.DiGraph()
    # add initial nodes
    G.add_node('root')
    G.add_node('sink')
    G.add_nodes_from(source_idx)
    G.add_nodes_from(target_idx)
    # add edges
    source_root_edges = list(itertools.product(['root'], source_idx, [{'capacity': 1, 'weight':0}]))
    G.add_edges_from(source_root_edges)

    capacities = get_target_capacities(source_idx, target_idx, capacity_method=capacity_method, seed=seed)
    target_sink_edges = [(target_idx[i], 'sink', {'capacity':capacities[i], 'weight':0}) for i in range(len(target_idx))]
    G.add_edges_from(target_sink_edges)

    return G


def build_graph(source_idx, target_idx, knn_source_idx, knn_target_idx, knn_dist,
                capacity_method='uniform', add_null=True, null_cost_percentile=95, seed=456):
    """Build a graph based on knn indices
    source_idx: nodes at the source of the graph
    target_idx: nodes at the target of the graph
    knn_source_idx: nodes at the source of the graph with knn connections
    knn_target_idx: nodes at the target of the graph with knn connections
    knn_dist: distance (cost) on the knn connections
    capacity_method: how to set capacities on the target to sink edged {uniform, inf, top, 1to1}
    add_null: whether to add the null node ot the graph
    null_cost_percentile: which percentile of the overall costs should correspond to the null match penalty
    seed: random seed

    output: directed graph object
    """
    G = build_graph_base(source_idx, target_idx, capacity_method=capacity_method, seed=seed)
    # add inter-technology connections
    source_target_edges = [(knn_source_idx[i], knn_target_idx[i],
                           {'capacity':1, 'weight':knn_dist[i]}) for i in range(len(knn_dist))]
    G.add_edges_from(source_target_edges)
    if(add_null):
        null_cost = get_null_cost(knn_dist, null_cost_percentile)
        G = extend_graph_null(G, source_idx, null_cost)

    return G

def convert_to_int(cost, factor):
    """Convert float values into integers by multiplying by factor and cropping the floating points
    cost: vector of float values
    factor: factor to multiply the values before cutting the decimals
    """
    cost = factor*cost
    cost = [int(x) for x in cost]

    return cost

def get_cost_knn_graph(source, target, factor=100, cost_type='distance', knn_k=10, knn_n_jobs=100,
                       capacity_method='uniform', add_null=True,
                       null_cost_percentile=95, seed=456):
    """Build an extended graph based on knn graph
    source: pandas dataframe which rows will be placed at the source of the graph (larger dataset)
    target: pandas dataframe which rows will be placed at the target of the graph (smaller dataset)
    factor: factor to multiply the cost values before cutting the decimals
    cost_type: {percentile} if costs should be converted into percentiles, {distance}: cost (Euclidean distance)
    knn_k: k neighbors to be found
    knn_n_jobs: how many processors to use for the job
    capacity_method: how to set capacities on the target to sink edged {uniform, inf, top, 1to1}
    add_null: whether to add the null node ot the graph
    null_cost_percentile: which percentile of the overall costs should correspond to the null match penalty
    seed: random seed

    output: directed graph object
    """
    knn_source_idx, knn_target_idx, knn_dist = get_knn_union(source, target, knn_k, knn_n_jobs)

    if(cost_type=='percentile'):
        len_p = int(len(knn_dist))
        p = np.linspace(min(knn_dist),max(knn_dist),100)
        knn_dist = np.digitize(knn_dist, bins=p)

    knn_dist = convert_to_int(knn_dist, factor)
    print('Max dist: ', np.max(knn_dist))

    source_idx = ['source_'+str(x) for x in range(source.shape[0])]
    target_idx = ['target_'+str(x) for x in range(target.shape[0])]

    G = build_graph(source_idx, target_idx, knn_source_idx, knn_target_idx, knn_dist,
                    capacity_method=capacity_method, add_null=add_null,
                    null_cost_percentile=null_cost_percentile, seed=seed)

    print('Number of nodes: ', len(G.nodes))
    print('Number of edges: ', len(G.edges))

    return G

def extract_matches_flow(flow_dict, keep_only='source'):
    """Extract the matched pairs
    flow_dict: ooutput from the max_flow_min_cost algorithm
    keep_only: discard all matches where keep_only doesn't appear
               (eg if keep_only=='source', then only inter-technology matches are reported)
    """
    matches_source = []
    matches_target = []
    for f in flow_dict.keys():
        if(keep_only is not None):
            if(keep_only not in f):
                continue
        matches_f = [x for x in flow_dict[f].keys() if flow_dict[f][x]>0]
        matches_source.extend([f]*len(matches_f))
        matches_target.extend(matches_f)
    matches = pd.DataFrame({'source': matches_source,
                            'target': matches_target})
    return matches


def mcmf(G):
    """ Use the Min-Cost Max-Flow algorithm to find the best matches between cells across technologies
    G: directed graph with 'root' and 'sink' nodes and {capacity, weights} attributes on edges

    output: {row indices, column indices}, corresponding to source and target datasets, respectively
    """
    t_start = time.process_time()

    flow_dict = nx.max_flow_min_cost(G,'root', 'sink', capacity='capacity', weight='weight')

    t_stop = time.process_time()
    t = (t_stop-t_start)
    print('MCMF took [s]: '+str(t))

    matches = extract_matches_flow(flow_dict, keep_only='source')
    source_idx = [x.split('_')[-1] for x in matches['source']]
    source_idx = [np.nan if x=='null' else int(x) for x in source_idx]
    target_idx = [x.split('_')[-1] for x in matches['target']]
    target_idx = [np.nan if x=='null' else int(x) for x in target_idx]

    return np.array(source_idx), np.array(target_idx)
