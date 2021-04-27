"""Node.py
Represents nodes in a PBN.
"""
import random
from .utils import *
import numpy as np

class Node():
    def __init__(self, function, i, name = None):
        """represents node in a PBN.

        args:
            mask [Node]: List of node objects that are inputs of this node.
            function [float]: matrix representation of function
            name (String): Name of the gene
        """
        self.input_nodes = None
        self.function = function
        self.i = i

        if type(name) == type(None):
            self.name = "G{0}".format(i)
        else:
            self.name = name

        self.state = None
        self.input_weights = None
        self.value = None

    def compute_next_value(self):
        """Return own next-state given the particular state according to own function and states of input genes.
        Wowza this is data-type madness
        """
        input_state = []

        for i in range(len(self.input_nodes)):
            input_node = self.input_nodes[i]
            input_state += [int(input_node.value)]
        input_state = tuple(input_state)

        prob_true = self.function.item(input_state)
        u = random.uniform(0,1) #Sample

        self.potential_value = u < prob_true

    def apply_next_value(self):

        if self.potential_value == type(None):
            raise Exception('Finishing transaction without computing next value.')
        self.value = self.potential_value
        self.potential_value = None



    def compute_input_weights(self):
        """BROKEN.BROKEN.BROKEN.BROKEN.BROKEN.BROKEN.BROKEN.BROKEN.BROKEN.BROKEN.BROKEN.BROKEN.BROKEN.
        """
        weights = np.zeros(np.sum(self.mask), dtype=float)
        N_states = 2**(np.sum(self.mask))
        for input_considered in range(np.sum(self.mask)):
            PT = 0
            PF = 0
            on_states = []
            off_states = []
            for state_index in range(N_states):
                state = booleanize(state_index, np.sum(self.mask))
                if state[input_considered]:
                    on_states += [state]
                else:
                    off_states += [state]
            for on_state in on_states:
                state_index = integerize(on_state)
                PT += self.function[state_index]
            for off_state in off_states:
                state_index = integerize(off_state)
                PF += self.function[state_index]
            PT /= np.sum(self.mask)
            PF /= np.sum(self.mask)
            weight = abs(PT-PF)
            weights[input_considered] = weight
        
        input_weights = np.zeros(self.mask.shape[0], dtype=float)

        i = 0
        for input_considered in range(self.mask.shape[0]):
            if self.mask[input_considered]:
                input_weights[input_considered] = weights[i]
                i += 1
        self.input_weights = input_weights

    def __str__(self):
        return self.name


