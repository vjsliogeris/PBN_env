"""PBN.py
The environment that runs PBNs.
"""
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import copy
import time
from .Node import Node
from .utils import *

class PBN():
    def __init__(self, PBN_data):
        """Construct a PBN from given PBN data.

        Attributes:
            state (bool]: curent state of the PBN.
            PBN_size (int): Number of nodes in the PBN.
            nodes [Node]: Node objects representing each of the variables.

        Args:
            PBN_data (list): data representing the PBN.

        returns:
            PBN
        """

        self.state = None
        self.PBN_size = len(PBN_data)
        self.nodes = np.empty((self.PBN_size), dtype=object)
        self.PBN = None
        self.STG = None

        for i in range(self.PBN_size):
            mask, function = PBN_data[i]
            self.nodes[i] = Node(mask, function, i)

    def unclip(self):
        """Un-clip each of the nodes
        """
        for node in self.nodes:
            node.unclip()

    def reset(self, state = None):
        """Set the state of the PBN to a particular one.

        args:
            state [bool]: The state to be set to. If left empty, defaults to a random state.
        """
        if type(state) == type(None):
            self.state = np.random.rand(self.PBN_size) > 0.5
        else:
            if state.shape[0] != self.PBN_size:
                raise Exception('The length of the state given ({0}) is different from the PBN size ({1}).'.format(state.shape[0], self.PBN_size))
            self.state = state.astype(bool)


    def name_nodes(self, names):
        for i in range(self.PBN_size):
            self.nodes[i].name = names[i]

    def print_PBN(self):
        """Construct a networkx graph representing the connetcivities of the PBN.

        returns: networkx di-graph.
        """
        if type(self.PBN) == type(None):
            G = nx.DiGraph()
            for i in range(self.PBN_size):
                G.add_node(self.nodes[i].name)
            for i in range(self.PBN_size):
                #For each target node
                node = self.nodes[i]
                inps = []
                mask = node.mask
                weights = node.input_weights
                if type(weights) == type(None):
                    for inp in range(len(mask)):
                        if mask[inp]:
                            inps += [inp]
                    for inp in inps:
                        G.add_edge(self.nodes[inp].name,self.nodes[i].name)
                else:
                    for inp in range(len(mask)):
                        if mask[inp]:
                            inps += [inp]
                    for inp in inps:
                        G.add_edge(self.nodes[inp].name,self.nodes[i].name, weight = weights[inp])
            self.PBN = G
        return self.PBN

    def plot_PBN(self, path):
        #self.generate_weights()
        G = self.print_PBN()
        plt.figure(1, figsize = (20,20))
        pos = nx.spring_layout(G)
        weights = [G[u][v]['weight'] for u,v in G.edges()]
        nx.draw(G, pos, node_size = 900, with_labels = True, width = weights)
#        nx.draw(G, pos, node_size = 900, with_labels = True)
        plt.savefig(path, dpi=160)
        plt.clf()

    def flip(self, index):
        """Flip the value of a gene at index.

        args:
            index (int): gene index to flip.
        """
        self.state[index] = not self.state[index]

    def get_funcs(self):
        """Print the functions of the PBN to inspect visually.
        """
        for i in range(self.PBN_size):
            print(self.nodes[i].function)

    def clip(self, gene_i, value):
        """Clip a gene to a particular value.

        The clipping of a gene makes the gene ignore inputs, staying at a constant value.
        Used to represent mutations.

        args:
            gene_i (int): gene index
            value (bool): value to clip the gene to
        """

        if type(value) != type(True):
            raise Exception('Keep it as a bool!.')

        i = 0
        for node in self.nodes:
            if gene_i == i:
                node.clip(value)
                self.state[i] = value
            i += 1

    def step(self):
        """Perform a step of natural evolution.
        """
        old_state = self.state

        if type(old_state) == type(None):
            raise Exception('No state initialised. (Run reset())')

        new_state = np.empty((self.PBN_size), dtype=bool)

        for i in range(self.PBN_size):
            new_state[i] = self.nodes[i].step(old_state)

        self.state = new_state

    def gen_STG(self):
        """Generate the State Transition Graph (STG) of the PBN.

        Go through each possible state.
        Compute the probabilities of going to next states.

        returns:
            networkx DiGraph.
        """
        if type(self.STG) == type(None):
            N_states = 2**(self.PBN_size)
            G = nx.DiGraph()
            start = time.time()
            for state_index in range(N_states):
                state = booleanize(state_index, self.PBN_size)
                G.add_node(str(state.astype(int)))
                next_states = self._compute_next_states(state)
                G.add_weighted_edges_from(next_states)
                end = time.time()
                est = N_states*(end-start)/(state_index+1)
                print("\rComputing STG: {0}%. Est duration: {1}s, OR {2} mins, OR {3} hrs".format(state_index*100 / N_states, est, est/60, est/3600), end="")
            self.STG = G
        return self.STG

    def plot_STG(self, path):
        STG = self.gen_STG()
        plt.figure(1, figsize = (20,20))
        pos = nx.spring_layout(STG)
        weights = {(u,v):STG[u][v]['weight'] for u,v in STG.edges()}
        #nx.draw(G, pos, node_size = 900, with_labels = True, width = weights)
        nx.draw(STG, pos, node_size = 900, with_labels = True)
        nx.draw_networkx_edge_labels(STG,pos,edge_labels = weights)
        plt.savefig(path, dpi=160)
        plt.clf()


    def compute_attr(self):
        if type(self.STG) == type(None):
            self.STG = self.gen_STG()
        STG = self.STG
        generator = nx.algorithms.components.attracting_components(STG)
        attractors = []
        for att_set in generator:
            att_list = []
            for state_str in att_set:
                att_list += [[int(s) for s in state_str[1:-1].split(" ")]]
            attractors += [att_list]
        return attractors


    def generate_weights(self):
        """Compute weights
        """
        for node_i in range(self.PBN_size):
            node = self.nodes[node_i]
            function = node.function
            mask = node.mask
            node.compute_input_weights()

    def apply_weights(self, weights):
        """Apply the weights provided (save them)
        """
        for i in range(self.PBN_size):
            node_weights = np.around(weights[i,:], 4)
            relevant_node = self.nodes[i]
            relevant_node.input_weights = node_weights
            self.nodes[i] = relevant_node

    def _compute_next_states(self, state):
        """Compute the probabilities of going to all next possible states from current state.

        Go through each gene. Compute the probability of each gene being True after current state.
        Convert those probabilities to next possible states.

        args:
            state [bool]: State to calculate probabilities from.

        returns:
            list of triplets. (Current State, next-possible-state, probability.)

        """
        probabilities = np.zeros((2, self.PBN_size), dtype=float)
        N_states = 2**(self.PBN_size)

        output = []
        for i in range(self.PBN_size):
            prob_true = self.nodes[i].get_probs(state)
            probs = np.array([1-prob_true, prob_true])
            probabilities[:,i] = probs
        protostates = self._probs_to_states(probabilities)
        for prostate, proprob in protostates:
            output += [(str(state.astype(int)), str(prostate.astype(int)), proprob)]
        return output

    def _probs_to_states(self, probs):
        """Compute the next possible states to go to, and their probabilities, given a set of probabilities of being true for each gene.

        Set the next states as a list of 0.5 with probability 1.
        A gene can not be at value 0.5, so it is used to signify an uncomputed value.

        Go through each gene.
        If probability is 1 or 0
            set the value of all next states at that index to the particular value. Leave probability unaffected.
        Else,
            Make two copies of all next states - one for each state of the gene. Compute probabilities accordingly.

        args:
            probs [float]: Probabilities of each gene being true in the next state.

        returns:
           List of tuples. Each tuple is a possible next state with according probability. 
        """
        _, n_genes = probs.shape

        protostate = np.ones(n_genes, dtype=float) * 0.5
        protoprob = 1
        
        prostate = [(protostate, protoprob)]
        for gene_i in range(n_genes):
            p = probs[:,gene_i]
            if p[0] == 1 or p[0] == 0:
                #Deterministic. Mainly for optimisation.
                for pro in prostate:
                    protostate, protoprob = pro #Go through each next-state already computed, unpack them
                    protostate[gene_i] = p[1] #Set the value of the gene to the corresponding value.
            else:
                prostate_copy = []
                for pro in prostate:
                    pro_1, prob_1 = copy.deepcopy(pro)
                    pro_2, prob_2 = copy.deepcopy(pro)

                    pro_1[gene_i] = 0 #Set value to 0
                    pro_2[gene_i] = 1 #Set value to 1
                    prob_1 *= p[0] #Set probability to that value being 0
                    prob_2 *= p[1] #^^^
                    #Put them back in.
                    protostate_1 = (pro_1, prob_1)
                    protostate_2 = (pro_2, prob_2)
                    prostate_copy += [protostate_1]
                    prostate_copy += [protostate_2]
                prostate = prostate_copy
        return prostate


