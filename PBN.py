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
import pickle

class PBN():
    def __init__(self, PBN_data = None):
        """Construct a PBN from given PBN data.

        Args:
            PBN_data (list): data representing the PBN.

        returns:
            PBN
        """
        self.PBN_size = len(PBN_data)
        self.nodes = np.empty((self.PBN_size), dtype=object)
        self.PBN = None
        self.STG = None

        for i in range(self.PBN_size):
            _, function = PBN_data[i]
            self.nodes[i] = Node(function, i)

        for i in range(self.PBN_size):
            mask, _ = PBN_data[i]
            input_nodes = self.nodes[mask]
            self.nodes[i].input_nodes = input_nodes

    def reset(self, state = None):
        """Set the state of the PBN to a particular one.

        args:
            state [bool]: The state to be set to. If left empty, defaults to a random state.
        """
        if type(state) == type(None):
            for node in self.nodes:
                node.value = np.random.rand() > 0.5
        else:
            if state.shape[0] != self.PBN_size:
                raise Exception('The length of the state given ({0}) is different from the PBN size ({1}).'.format(state.shape[0], self.PBN_size))
            for i in range(self.PBN_size):
                self.nodes[i].value = state[i].astype(bool)


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
                node = self.nodes[i] #Current node object
                input_nodes = self.nodes[i].input_nodes
                inps = [] #List of names of input nodes
                inps_i = [] #List of indexes of input nodes

                weights = node.input_weights
                for inp in input_nodes:
                    inps += [inp.name]
                    inps_i += [inp.i]

                if type(weights) == type(None):
                    for inp in input_nodes:
                        G.add_edge(inp.name,node.name)
                else:
                    for inp in input_nodes:
                        G.add_edge(inp.name,node.name, weight = weights[inp.i])
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

    def step(self):
        """Perform a step of natural evolution.
        """

        for node in self.nodes:
            node.compute_next_value()

        for node in self.nodes:
            node.apply_next_value()

    def get_state(self):
        """Get a state from the values of all the nodes
        """
        state = np.empty(self.PBN_size, dtype=bool)
        for i in range(self.PBN_size):
            state[i] = self.nodes[i].value
        return state
        

    def gen_STG(self):
        """Generate the State Transition Graph (STG) of the PBN.

        Go through each possible state.
        Compute the probabilities of going to next states.

        returns:
            networkx DiGraph.
        """
        if type(self.STG) == type(None):
            N_states = 2**(self.PBN_size)
            print("Total number of states: {0}".format(N_states))
            G = nx.DiGraph()
            start = time.time()
            for state_index in range(N_states):
                state = booleanize(state_index, self.PBN_size)
                G.add_node(str(state.astype(int)))
                next_states = self._compute_next_states(state)
                G.add_weighted_edges_from(next_states)
                end = time.time()
                est = N_states*(end-start)/(state_index+1)
                print("\rComputing STG: At index {4} {0}%. Est duration: {1}s, OR {2} mins, OR {3} hrs".format(state_index*100 / N_states, est, est/60, est/3600, state_index), end="")
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

    def get_node_by_name(self, nodename):
        """Get the appropriate node object given the name of the node.
        """
        for node in self.nodes:
            if node.name == nodename:
                return node
        raise Exception(f'Node with name \'{nodename}\' not found.')

    def get_output_nodes(self, node_obj):
        """Get all the nodes that take node_obj as input.
        May deal poorly if there are duplicate names.
        NOTE: UNTESTED
        """
        output_nodes = []
        for node in self.nodes:
            inputs = [n.name for n in self.nodes[node.mask]]
            if node_obj.name in inputs:
                output_nodes += [node]
        return output_nodes

    def get_subPBN(self, nodes):
        """Get the subgraph of the PBN which includes the given nodes.
        """
        sub_PBN = PBN()
        PBN_nodes = []
        for node in nodes:
            node_object = self.get_node_by_name(node)
            PBN_nodes += [node_object]

        sub_PBN.nodes = PBN_nodes
        sub_PBN.resolved = False
        sub_PBN.PBN_size = len(PBN_nodes)
        return sub_PBN

    def all_parents_resolved(self):
        """Go through all parents, check if they are resolved
        """
        all_resolved = True
        for parent in self.parent_SG:
            if not parent.resolved:
                all_resolved = False
        return all_resolved


    def compute_attractors(self, expand = False):
        """Compute attractors without explicitly computing STG.
        """
        template = '*'*self.PBN_size #Generate rule template.

        all_symbolic_rules = []
        for node in self.nodes:
            symbolic_rules = node.generate_symbolic_rules(list(template))
            all_symbolic_rules+= [symbolic_rules]

        print("Original ruleset")
        print(all_symbolic_rules)
        flat_ruleset = flatten_ruleset(all_symbolic_rules)
        print("Flat rules")
        print(flat_ruleset)
        input()
        reached_states = template
#        reached_states = '0111' #Testing
#        reached_states = '100*' #Testing
#        reached_states = '010*' #Testing
#        reached_states = '0*1*' #Testing
#        reached_states = '**1*' #Testing
#        reached_states = '*1**' #Testing
#        input()
        reached_states = apply_ruleset(reached_states, flat_ruleset)
        print("Reached {0}".format(reached_states))
        loopy_states = []
        for state in reached_states:
            loopy_state = find_successors(state, [template], flat_ruleset)
            if not loopy_state in loopy_states:
                loopy_states += [loopy_state]
        print()
        print("Loopy states: {0}".format(loopy_states))
        input()
        attractors = []
        for loopy_state in loopy_states:
            print("Loopy state deepening: {0}".format(loopy_state))
            attractor_states = copy.deepcopy([loopy_state])
            attractors_old = []
            while not attractor_states == attractors_old:
                next_states = []
                for state in attractor_states:
                    next_state = apply_ruleset(state, flat_ruleset)
                    for n_s in next_state:
                        if not n_s in next_states and not check_subset(n_s, loopy_state):
                            next_states += [n_s]
                next_states += [loopy_state]
                next_states.sort()
                print(next_states)
                attractors_old = copy.deepcopy(attractor_states)
                attractor_states = next_states
#                input()
            if not attractor_states in attractors:
                attractors += [attractor_states]
        print("Loopy state + successors: {0}".format(attractors))
        attractors_compact = copy.deepcopy(attractors)
        attractors = []
        for attractor_compact in attractors_compact:
            print("Loopy state + succ: {0}".format(attractor_compact))
            attractor = []
            for state_compact in attractor_compact:
                expanded = expand_notation(state_compact)
                for s in expanded:
                    if not s in attractor:
                        attractor += [s]
            attractor.sort()
            if not attractor in attractors:
                attractors += [attractor]
                print(attractor)
        attractors.sort()

        print(attractors)
        pickle.dump(self, open("last_PBN.pkl","wb"))
#        raise Exception('')
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


