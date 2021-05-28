"""Node.py
Represents nodes in a PBN.
"""
import random
from .utils import *
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


    def generate_symbolic_rules(self, template):
        """Generate symbolic rules for this nodes functions.
        """
        #Find indexes where the output is deterministic (0 or 1, guaranteed output.)
        F = copy.deepcopy(self.function)
        np.mod(F, 1, out=F)
        determined_outputs = (F == 0)
        index_det = np.where(determined_outputs)

        n_inputs = len(index_det)
        n_determined_combinations = len(index_det[0])
        condition_matrix = np.empty((n_inputs, n_determined_combinations), dtype=bool)
        for input_index in range(n_inputs):
            #Iterate over inputs
            for combination_index in range(n_determined_combinations):#Iterate over values where output is determined
                input_value_at_combination = index_det[input_index][combination_index]
                condition_matrix[input_index, combination_index] = input_value_at_combination
        actual_inputs_determined = []
        for i in range(n_determined_combinations):
            actual_inputs_determined += [condition_matrix[:,i]]

        rules = []
        apparent_input_indexes = []
        for child in self.input_nodes:
            apparent_input_indexes += [child.i]
        #Compute actual rules
        temp_template_inp = copy.deepcopy(template)
        temp_template_out = copy.deepcopy(template)
        for actual_input in actual_inputs_determined:
            for input_i in range(n_inputs):
                input_value = int(actual_input[input_i])
                input_position_apparent = apparent_input_indexes[input_i]
                temp_template_inp[input_position_apparent] = str(input_value)
            output_value = int(self.function.item(tuple(actual_input.astype(int))))
            apparent_index_output = self.i
            temp_template_out[apparent_index_output] = str(output_value)
            r = [("".join(temp_template_inp),"".join(temp_template_out))]
            rules += r
            #rules["".join(temp_template_inp)] = ["".join(temp_template_out)]
        return rules

    def __str__(self):
        return "{0}: {1}".format(self.name, self.value)
