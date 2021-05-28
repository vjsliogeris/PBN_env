"""Node.py
Represents nodes in a PBN.
"""
import random
import numpy as np
import copy

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

    def get_probs(self, state):
        input_indices = []
        for i in range(len(self.input_nodes)):
            input_indices += [self.input_nodes[i].i]
        input_state = state[input_indices].astype(int)
        input_state = tuple(input_state)
        prob_true = self.function.item(input_state)
        return prob_true

    def apply_next_value(self):

        if self.potential_value == type(None):
            raise Exception('Finishing transaction without computing next value.')
        self.value = self.potential_value
        self.potential_value = None


    def __str__(self):
        return "{0}: {1}".format(self.name, self.value)
