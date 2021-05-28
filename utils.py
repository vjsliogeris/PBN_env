"""utils.py
"""
import copy
import numpy as np
import itertools


def integerize(X):
    """Turn a list of booleans in to an integer representation.
    """
    n_in = X.shape[0]
    s = 0
    for i in range(n_in):
        s += int(X[n_in-i-1]) * 2**i
    return s

def booleanize(X, l):
    """Turn an int to a boolean string of length l
    """
    output = np.zeros((l), dtype=bool)
    for i in range(l):
        h = 2**(l-i-1)
        if X >= h:
            X -= h
            output[i] = 1
    return output

def apply_until_exhaustion(ruleset_linked, function):
    """Repeat the manipulation until the ruleset doesn't change.
    This means all the possible manipulations have been applied.
    """
    ruleset = copy.deepcopy(ruleset_linked)
    can_generate_new = True
    ruleset_old = copy.deepcopy(ruleset)
    i = 0
    while can_generate_new:
        ruleset = function(copy.deepcopy(ruleset_old))
        can_generate_new = not (ruleset == ruleset_old)
        ruleset_old = ruleset
        i += 1

    return ruleset


def _union(list_1, list_2):
    out = copy.deepcopy(list_1)
    for state in list_2:
        if not state in out:
            out += [state]
    return out

def _list_intersection(list_1, list_2):
    conjunctions = []
    for state_1 in list_1:
        for state_2 in list_2:
            conjunctions += [_intersection(state_1, state_2)]
    tidy = _tidy_states(conjunctions)
    return tidy

def _union_rules(ruleset):
    fitting_rules = copy.deepcopy(ruleset)
    used_rules = []
    NU_rules = []
    for inp_1, out_1 in fitting_rules:
        for inp_2, out_2 in fitting_rules:
            if not (inp_1 == inp_2 and out_1 == out_2):
                NU_exists, NU = _nice_union(inp_1, inp_2)
                if NU_exists:
                    out = _union(out_1, out_2)
                    out.sort()
                    if not (inp_1, out_1) in used_rules:
                        used_rules += [(inp_1, out_1)]
                    if not (inp_2, out_2) in used_rules:
                        used_rules += [(inp_2, out_2)]
                    if not (NU, out) in NU_rules:
                        NU_rules += [(NU, out)]
                    #input()

    for used_rule in used_rules:
        fitting_rules.remove(used_rule)
    fitting_rules += NU_rules
    return fitting_rules

def _specify_rules(ruleset):
    fitting_rules = copy.deepcopy(ruleset)
    specific_rules = []
    used_rules = []
    for inp_1, out_1 in fitting_rules:
        working_output = copy.deepcopy(out_1)
        working_input = copy.deepcopy(inp_1)
        for inp_2, out_2 in fitting_rules:
            if not (inp_1 == inp_2 and out_1 == out_2):
                if _check_subset(inp_1, inp_2):#if rule_1 applies to state, rule 2 applies to state
                    working_input = _intersection(working_input, inp_2)
                    working_output = _list_intersection(out_2, working_output)
                    if not (inp_1, out_1) in used_rules:
                        used_rules += [(inp_1, out_1)]
                    if not (inp_2, out_2) in used_rules:
                        used_rules += [(inp_2, out_2)]

        working_output.sort()
        if not (working_input, working_output) in specific_rules:
            specific_rules += [(working_input, working_output)]
        #print("Finishing rule:")
        #_print_rule(inp_1, out_1)
        #print(working_output)
        #input()
    fitting_rules = specific_rules
    return fitting_rules

def shake_em_up(ruleset):
    fitting_rules = copy.deepcopy(ruleset)
    fitting_rules = _specify_rules(fitting_rules)
    print("Intersectioned rules")
    print(fitting_rules)
    fitting_rules = _union_rules(fitting_rules)
    print("Unioned_rules")
    print(fitting_rules)
    return fitting_rules

def apply_ruleset(state, ruleset, debug = False):
    """Apply the ruleset to the states

    Returns the resulting states after applying the ruleset to the input states.
    Recursive
    WIP
    """
    outputs = []
    print("Applying ruleset to state {0}".format(state))
    fitting_rules = []
    for inp, out in ruleset:
        intersection = _intersection(inp, state)
        if not type(intersection) == type(None):
#                _print_rule(inp, out)
            fitting_rules += [(intersection, out)]
    fitting_rules.sort(key = lambda x: _count_asterisks(x[0]), reverse = True)
    print("Fitting rules")
    print(fitting_rules)
    fitting_rules = apply_until_exhaustion(fitting_rules, lambda x: shake_em_up(x))
#    fitting_rules = _union_rules(fitting_rules)
#    print("Unioned_rules")
#    print(fitting_rules)
#
#    fitting_rules = _specify_rules(fitting_rules)

#    print("Intersectioned rules")
#    print(fitting_rules)
    outputs = []
    for inp, out in fitting_rules:
        if not out in outputs:
            outputs += out
    unaccounted_states = copy.deepcopy(state)
    for inp, out in fitting_rules:
        unaccounted_states = _subtract_states(unaccounted_states, inp)
    if len(unaccounted_states) > 0:
        print("Got unaccounted-for states")
        print(unaccounted_states)
        print(ruleset)
        unaccounted_state_outputs = apply_ruleset(unaccounted_states[0], ruleset)
        print("Outputs are {0}".format(unaccounted_state_outputs))
        raise Exception('Got unaccounted-for states')
    outputs = _tidy_states(outputs)
    outputs.sort(key = lambda x: _count_asterisks(x))
    print("Reached {0}".format(outputs))
    if debug:
        input()
    return outputs


def find_successors(state, parents, ruleset, debug = False):
#    print("_find_successors")
#    print(f"on {state}")
#    print(f" parents until now: {parents}")
    next_states = apply_ruleset(state, ruleset, debug = debug)
#    print(f" s_t+1{next_states}")
    goes_to_parents = True
    for n_s in next_states:
        if not n_s in parents:
            goes_to_parents = False
#    print(goes_to_parents)
#    input()
    if goes_to_parents:
        return state
    else:
        loopy_states = []
        for n_s in next_states:
            if not n_s in parents:
                loopy_state = find_successors(n_s, parents + [state], ruleset)
                if not type(loopy_state) == type(None):
                    return loopy_state
        raise Exception('Did not find loopy states??')

def _tidy_states(states):
    output = []
    for state_1 in states:
        redundant = False
        for state_2 in states:
            if not state_1 == state_2:
                if _check_subset(state_1, state_2):
                    redundant = True
        if not redundant:
            if not state_1 in output:
                output += [state_1]
    return output

def _subtract_states(A, B):
    """ Return A/B

    Find where there are non-asterisks.
    Get all combinations of the non-asterisks possible.

    """

    B_relevant = _intersection(A,B)

    n_elements = len(A)

    important_positions = []

    for i in range(n_elements):
        s_A = A[i]
        s_B = B[i]
        if s_A == '*' and not s_B == '*':
            important_positions += [(i, s_B)]
    combinations = list(itertools.product(['0', '1'], repeat = len(important_positions)))
    combinations = [''.join(x) for x in combinations]
    present_combination = [''.join(x[1]) for x in important_positions]
    present_combination = ''.join(present_combination)
    combinations.remove(present_combination)
    output = []
    for combination in combinations:
        c_expression = list(copy.deepcopy(A))
        for position_index in range(len(important_positions)):
            particular_position = important_positions[position_index]
            c_expression[particular_position[0]] = combination[position_index]
        c_expression = ''.join(c_expression)
        output += [c_expression]
    return output

def _intersection(temp_1, temp_2):
    """Return the intersection of two expressions
    """
    if type(temp_1) == type(None) or type(temp_2) == type(None):
        #Intersection with an empty set is an empty set.
        return None
    n_symbols = len(temp_1)
    merged_temp = ['*']*n_symbols
    for i in range(n_symbols): #Compare element-wise.
        s_1 = temp_1[i]
        s_2 = temp_2[i]
        if not s_1 == '*' and not s_2 == '*' and not s_1 == s_2:#If elements aren't asterisks but they mismatch
            #Means that there can't be an intersection because the node at that position must be at different values.
            return None
        if not s_1 == "*":
            merged_temp[i] = s_1
        if not s_2 == "*":
            merged_temp[i] = s_2
    merged_temp = ''.join(merged_temp)
    return merged_temp

def check_subset(A,B):
    return _check_subset(A,B)

def _check_subset(A, B):
    """Check if A is a subset of B (All states in A belong in B)
    """
    if type(A) == type(None): #The empty set is the subset of any set.
        return True
    if type(B) == type(None): #A non-empty set is not a subset of an empty set
        return False

    is_subset = True
    n_elements = len(A)
    for i in range(n_elements):
        a = A[i]
        b = B[i]
        if not b == '*' and not a == b:
            is_subset = False
    return is_subset

def flatten_ruleset(ruleset):
    """Remove the sub-listings in the ruleset.
    Sort it.
    """
    output = []
    for rule in ruleset:
        output += rule

    output.sort(key = lambda x: _count_asterisks(x[0]), reverse = True)
    return output

def _count_asterisks(temp):
    """Return the number of asterisks within the expression
    """
    n_elements = len(temp)
    n_asterisks = 0
    for i in range(n_elements):
        if temp[i] == '*':
            n_asterisks += 1
    return n_asterisks

def get_nice_specific(ruleset):
    """If two rules have the same input, merge outputs, since both of the rules must apply.
    """
    print("Generating specifics")
    ruleset = copy.deepcopy(ruleset)
    mergees = copy.deepcopy(ruleset[-1])
    n_output_ruleset = len(ruleset) #Number of output nodes considered by the ruleset
    for i in range(n_output_ruleset-1): #The final index of that ruleset is for mergees
        for j in range(i+1, n_output_ruleset):
            ruleset_i = ruleset[i]
            ruleset_j = ruleset[j]
            for inp_1, out_1 in ruleset_i:
                for inp_2, out_2 in ruleset_j:
                    if inp_1 == inp_2 and not out_1 == out_2:
                        out_merged = _intersection(out_1, out_2)
#                        _print_rule(inp_1, out_merged)
                        if not (inp_1, out_merged) in mergees:
                            mergees += [(inp_1, out_merged)]
    ruleset[-1] = mergees
    return ruleset

def remove_redundancies(ruleset):
    """Remove redundant rules:
    If inp_1 == inp_2, but out_1 is subset of out_2, then replace it with inp_1 --> out_1
    If inp_1 is subset of inp_2 and out_1 is subset of out_2, then replace it with inp_2 --> out_1 (more general to more specific)
    """
    print("Removing rules that are already described with other rules")
    ruleset = copy.deepcopy(ruleset)
    for inp_1, out_1 in ruleset:
        for inp_2, out_2 in ruleset:
            if not (inp_1, out_1) == (inp_2, out_2):
                if _check_subset(inp_1, inp_2): #Inp 1 is subset of inp 2
                    if _check_subset(out_1, out_2):
                        if (inp_1, out_1) in ruleset:
                            ruleset.remove((inp_1, out_1))
                        if (inp_2, out_2) in ruleset:
                            ruleset.remove((inp_2, out_2))
                        if not (inp_2, out_1) in ruleset:
                            ruleset += [(inp_2, out_1)]
                    if _check_subset(out_2, out_1):
                        if (inp_1, out_1) in ruleset:
                            ruleset.remove((inp_1, out_1))
                        if (inp_2, out_2) in ruleset:
                            ruleset.remove((inp_2, out_2))
                        if not (inp_2, out_2) in ruleset:
                            ruleset += [(inp_2, out_2)]
                if _check_subset(inp_2, inp_1): #inp 2 is subset of inp 1
                    if _check_subset(out_1, out_2):
                        if (inp_1, out_1) in ruleset:
                            ruleset.remove((inp_1, out_1))
                        if (inp_2, out_2) in ruleset:
                            ruleset.remove((inp_2, out_2))
                        if not (inp_1, out_1) in ruleset:
                            ruleset += [(inp_1, out_1)]
                    if _check_subset(out_2, out_1):
                        if (inp_1, out_1) in ruleset:
                            ruleset.remove((inp_1, out_1))
                        if (inp_2, out_2) in ruleset:
                            ruleset.remove((inp_2, out_2))
                        if not (inp_1, out_2) in ruleset:
                            ruleset += [(inp_1, out_2)]
    ruleset.sort(key = lambda x: _count_asterisks(x[0]), reverse = True)
    return ruleset

def get_nice_general_rules(ruleset):
    """Given a starting ruleset, extrapolate any ,,nice`` general rules.

    Nice rules are rules that only have one expression in input and output.
    
    For two rules to be candidates for being nice, they have to have the same output.

    This means that the rulesets for the same output can only be considered.
    """
    print("Generalizing rules")
    ruleset = copy.deepcopy(ruleset)
    n_output_ruleset = len(ruleset) #Number of output nodes considered by the ruleset
    for i in range(n_output_ruleset):
        considered_ruleset = ruleset[i]
        n_rules = len(considered_ruleset)
        for j in range(n_rules):
            for k in range(j+1, n_rules):
                inp_1, out_1 = considered_ruleset[j]
                inp_2, out_2 = considered_ruleset[k]
                if out_1 == out_2:#Outputs match
                    is_nice, union = _nice_union(inp_1, inp_2)
                    if is_nice and not (union, out_1) in considered_ruleset:
#                        _print_rule(union, out_1)
                        considered_ruleset += [(union, out_1)]
    return ruleset

def _nice_union(temp_1, temp_2):
    """Check if a ,,nice'' union between the two expressions is possible, and perform it if so.

    A ,,nice`` expression is such that it only requires one expression to show, i.e. the new 
    expression may be doable by replacing a character with an asterisk, thus making it more
    general, yet applicable to both.

    
    """
    if type(temp_1) == type(None) and type(temp_2) == type(None):#If both are empty sets
        return None #Return the empty set.
    #If one of the sets is empty, then the union will be the other set.
    if type(temp_1) == type(None):
        return temp_2
    if type(temp_2) == type(None):
        return temp_1
    #Both sets here are non-empty then
    n_elements = len(temp_1)
    
    differing_elements = []
    for i in range(n_elements):
        s_1 = temp_1[i]
        s_2 = temp_2[i]
        if not s_1 == s_2:
            differing_elements += [(i, s_1, s_2)]
    if len(differing_elements) == 1:
        #Only one differing elements there.
        position, e_1, e_2 = differing_elements[0]
        if not (e_1 == '*' or e_2 == '*'):#Both are not asterisks
            output = list(copy.deepcopy(temp_2))
            output[position] = '*'
            output = ''.join(output)
            return True, output
    return False, None


    raise Exception('')

def _print_rule(inp, out):
    """I've been using this so much I'm writing a seperate function for it.
    """
    print("{0} --> {1}".format(inp, out))

def _check_if_universe(states):
    """Check if union of set of states composes all possible states
    """
#    print("_check_if_universe({0})".format(states))
    n_elements = len(states[0])
    starting_state = ['*' * n_elements]
    for f_state in states:
        s_state_old = starting_state
        starting_state = _subtract_states(starting_state, [f_state])
#        print("{0} \ {1} = {2}".format(s_state_old, [f_state], starting_state))
#        input()
        if len(starting_state) == 0:
            return True
    print("Example reachable states: {0}".format(starting_state))
    return False




def _check_complement(f_1, f_2):
    n_elements = len(f_1)
    diff_indices = []
    for i in range(n_elements):
        if not f_1[i] == f_2[i]:
            diff_indices += [i]
    if not len(diff_indices) == 1:
        return None, False
    diff_index = diff_indices[0]
    if not f_1[diff_index] == "*" and not f_2[diff_index] == "*":
        return diff_index, True
    return None, False

def expand_notation(attractor):
    #input()
    #print("Calling _expand_notation on {0}".format(attractor))
    first_asterisk = None
    n_elements = len(attractor)
    for i in range(n_elements):
        if attractor[i] == '*':
            first_asterisk = i
            break
    if type(first_asterisk) == type(None):
        #print("Reached leaf. Returning {0}".format([attractor]))
        return [attractor]
    else:
        #print("{0}: Going deeper.".format(attractor))
        att_0 = copy.deepcopy(list(attractor))
        att_1 = copy.deepcopy(list(attractor))
        att_0[first_asterisk] = '0'
        att_1[first_asterisk] = '1'
        att_0 = ''.join(att_0)
        att_1 = ''.join(att_1)
        att_0 = expand_notation(att_0)
        att_1 = expand_notation(att_1)
        #print("{1}: This became a leaf. Returning {0}".format(att_0+ att_1, attractor))
        return att_0+ att_1


def expand_att(att):
    output = []
#    print("att to expand: {0}".format(att))
    for state in att:
        n_asterisks = _count_asterisks(state)
        if n_asterisks == 0:
            output += [state]
        else:
#            print(state)
            all_combos = _unasterisk(state)
#            print(all_combos)
            output += all_combos
            #print(all_combos)
#            print(output)
#    print("expanded att: {0}".format(output))
    return output

