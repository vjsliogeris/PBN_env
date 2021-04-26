"""Node.py
Represents nodes in a PBN.
"""
import random
from .utils import *

class Node():
    def __init__(self, mask, function, i, name = None):
        """represents node in a PBN.

        attribute:
            clipped (bool): Signified if gene is clipped.
                A clipped gene stays constant.

        args:
            mask [bool]: input mask. Signifies inputs
            function [float]: vector representation of function
            name (String): Name of the gene. Not actually used for now.
        """
        self.mask = mask
        self.function = function
        if type(name) == type(None):
            self.name = "G{0}".format(i)
        else:
            self.name = name
        self.clipped = False
        self.clip_value = None
        self.input_weights = None

    def step(self, state):
        """Return own next-state given the particular state according to own function and states of input genes.

        args:
            state [int]: Current state of entire PBN

        returns:
            bool

        """
        if self.clipped:
            return self.clip_value

        inputs = state[self.mask] #Get state of input genes
        index = integerize(inputs) #Make state into index
        prob_true = self.function[index] #Find probability of being True given the state
        u = random.uniform(0,1) #Sample
        return prob_true > u

    def get_probs(self, state):
        """Return the probability of being True given state.

        args:
            state [int]: State

        returns:
            float

        """
        if self.clipped:
            #Very non-elegant.
            if self.clip_value:
                return 1
            else:
                return 0

        inputs = state[self.mask]
        index = integerize(inputs)
        prob_true = self.function[index]
        return prob_true

    def compute_input_weights(self):
        """Compute the importance of each input after the fact.
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

    def clip(self, clip_value):
        """Clip current gene.
        args:
            clip_value (bool): value to be clipped at
        """
        if not self.clipped:
            self.clipped = True
            self.clip_value = clip_value
        else:
            print('Already Clipped!')

    def unclip(self):
        """Unclip current gene
        """
        if self.clipped:
            self.clipped = False
            self.clip_value = None


