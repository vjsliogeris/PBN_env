"""utils.py
"""
import copy
import numpy as np


def _check_subset(A, B):
    """Check if A is a subset of B (All states in A belong in B)
    """
    is_subset = True
    n_elements = len(A)
    for i in range(n_elements):
        a = A[i]
        b = B[i]
        if not b == '*' and not a == b:
            is_subset = False
    return is_subset

def connect_symbolic_ruleset(ruleset):
    """Connect rules
    These rules, unlike the previous ones, mean that the input will eventually, yet invariably reach the target.
    """
    #print("connect_symbolic_ruleset")
    #print(ruleset)
    for inp_1, out_1 in ruleset:
        for inp_2, out_2 in ruleset:
#            if not (inp_1 == inp_2 and out_1 == out_2): #Check if not comparing to selfi
            is_subset = _check_subset(out_1, inp_2)
            if is_subset:
                '''
                print("{0} -> {1}".format(inp_1, out_1))
                print("{0} -> {1}".format(inp_2, out_2))
                print("{0} is subset of {1}".format(out_1, inp_2))
                print("{0} -> {1}".format(inp_1, out_1))
                print("+")
                print("{0} -> {1}".format(inp_2, out_2))
                print("=")
                print("{0} -> {1}".format(inp_1, out_2))
                print("-----------------")
                #'''
                if not (inp_1, out_2) in ruleset:
                    ruleset += [(inp_1, out_2)]
#                    new_rules[inp_1] = [out_2]
                if not (out_1, out_2) in ruleset:
                    ruleset += [(out_1, out_2)]
    #print(ruleset)
#    ruleset = _update_ruleset(ruleset, new_rules)
    #raise Exception('')
    return ruleset

def apply_until_exhaustion(ruleset_linked, function):
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

def concretise_symbolic_ruleset(ruleset):
    #print("concretise_symbolic_ruleset")
#    print(ruleset)
    ruleset_new = copy.deepcopy(ruleset)
    for inp_1, out_1 in ruleset:
        n_ast_out = _count_asterisks(out_1) 
        n_ast_inp = _count_asterisks(inp_1)
        compatibility = _check_subset(out_1, inp_1)
        if n_ast_out < n_ast_inp and compatibility: #Larger set goes into smaller
#            print("{0} -> {1}".format(inp_1, out_1))
            if not (out_1, out_1) in ruleset_new:
                ruleset_new += [(out_1, out_1)]
    return ruleset_new

def remove_redundant_rules(ruleset):
    #print("remove_reduntant_rules")
#    print(ruleset)
    new_ruleset = []
    for inp_1, out_1 in ruleset:
        mega_out = copy.deepcopy(out_1)
        for inp_2, out_2 in ruleset:
            if not (inp_1 == inp_2 and out_1 == out_2): #Check if not comparing to self
                if inp_1 == inp_2 and not out_1 == out_2:
#                    print("Same input, output mismatch")
#                    print("{0} -> {1}".format(inp_1, out_1))
#                    print("{0} -> {1}".format(inp_2, out_2))
#                    print(ruleset)
                    if _check_compatibility(out_1, out_2):
                        mega_out = _merge(mega_out, out_2)
                    else:
                        raise Exception('Same input, incompatible outputs')
#                            print(ruleset)
#        print("{0} -> {1}".format(inp_1, mega_out))
        if not (inp_1, mega_out) in new_ruleset:
            new_ruleset += [(inp_1, mega_out)]
#    print(new_ruleset)
#    raise Exception('')
    return new_ruleset

def simplify_symbolic_rules(ruleset):
    """Simplify pairs of rules
    If two rules are complement, and they have the same output, they can be merged into one (generalized)
    """
    #print("simplify_symbolic_ruleset")
    #If an input is not in the ruleset, we can say that the input leads to all other states.
#    print("vvv")
#    print(ruleset)
    for inp_1, out_1 in ruleset:
        for inp_2, out_2 in ruleset:
            if not (inp_1 == inp_2 and out_1 == out_2): #Check if not comparing to self
                diff_index, complements = _check_complement(inp_1, inp_2)
                if complements and (out_1 == out_2):
                    #print("Simplifiable")
                    #print("{0} -> {1}".format(inp_1, out_1))
                    #print("{0} -> {1}".format(inp_2, out_2))
                    inp_merged = list(inp_1)
                    inp_merged[diff_index] = '*'
                    inp_merged = ''.join(inp_merged)
                    if not (inp_merged, out_1) in ruleset:
                        ruleset += [(inp_merged, out_1)]
#    print(ruleset)
#    print(new_rules)
    #ruleset = _update_ruleset(ruleset, new_rules)
    #ruleset += new_rules
#    print(ruleset)
    #input()
    return ruleset

def merge_symbolic_rules(ruleset):
    """Merge rules that we can merge.
    If inputs and outputs are compatible, then the new rule would be the intersection of both rules.
    """
#    print("merge_symbolic_rules")
#    print(ruleset)
    n_rules = len(ruleset)
    n_oper = 0
    ruleset_new = []
    print("Ruleset size: {0}".format(n_rules))
    print("Total computations to do: {0}".format(n_rules**2))
    for inp_1, out_1 in ruleset:
        for inp_2, out_2 in ruleset:
            if not (inp_1 == inp_2 and out_1 == out_2): #Check if not comparing to self
                n_rules = len(ruleset)
                print("{0}".format(n_oper/(n_rules**2)), end="\r")
                compatible = _check_compatibility(inp_1, inp_2)
                if compatible:
                    inp_merged = _merge(inp_1, inp_2)
                    out_merged = _merge(out_1, out_2)
                    '''
                    print()
                    print("{0} -> {1}".format(inp_1, out_1))
                    print("+")
                    print("{0} -> {1}".format(inp_2, out_2))
                    print("=")
                    print("{0} -> {1}".format(inp_merged, out_merged))
                    print("^^^^^")
                    input()
                    '''
                    if not (inp_merged, out_merged) in ruleset:
                        ruleset += [(inp_merged, out_merged)]
            n_oper += 1
    #print("Total computaitons done: {0}".format(n_oper))
    #ruleset = _update_ruleset(ruleset, full_ruleset)
    #ruleset += full_ruleset
    #print()
#    print(ruleset)
#    print(full_ruleset)
    #input()
    return ruleset

def _update_ruleset(ruleset, additions):
#    print("_update_ruleset")
    #print(ruleset)
    #print(additions)
    for inp_a, out_a in additions:
        already_in = False
        if not inp_a in ruleset.keys():
            #This is a new input rule
            ruleset[inp_a] = [out_a[0]]
        else:
            #This rule already exist. Get the one that is more accurate
            out_r = ruleset[inp_a][0]
            if not out_r == out_a:
                #Outputs are different
                n_ast_a = _count_asterisks(out_a)
                n_ast_r = _count_asterisks(out_r)
                if n_ast_r < n_ast_a:
                    ruleset[inp_a] = [out_r[0]]
                elif n_ast_r > n_ast_a:
                    ruleset[inp_a] = [out_a[0]]
                else:
                    #print("Original: {0} -> {1}".format(inp_a, out_r))
                    #print("Addition: {0} -> {1}".format(inp_a, out_a))
                    #print(n_ast_r)
                    #print(n_ast_a)
                    out_merged = _merge(out_r, out_a)
                    #print(out_merged)
                    ruleset[inp_a] = [out_merged]
#    print(ruleset)
#    input()

    return ruleset

def _count_asterisks(notation):
    n_elements = len(notation)
    n_asterisks = 0
    for i in range(n_elements):
        if notation[i] == '*':
            n_asterisks += 1
    return n_asterisks

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

def _check_related(temp_1, temp_2):
    n_symbols = len(temp_1)
    related = True
    for i in range(n_symbols):
        e_1 = temp_1[i]
        e_2 = temp_2[i]
        if not e_1 == '*' and not e_1 == e_2:
            related = False
    return related

def get_att_from_ruleset(ruleset, ruleset_timeless):
    #print("get_att_from_ruleset")
#    print(ruleset)
#    print(ruleset_timeless)
    states_in_att = []
    for inp,out in ruleset_timeless:
        if inp == out:
            states_in_att +=[inp]
    states_in_att.sort(key = lambda x: _count_asterisks(x))
    print("States deemed to be in attractors: {0}".format(states_in_att))
    #Pair-wise compare and merge.
    combos = []
    for a_1 in states_in_att:
        for a_2 in states_in_att:
            compatible = _check_compatibility(a_1, a_2)
            if compatible:
                product = _merge(a_1, a_2)
                if not product in combos:
                    combos += [product]
    print("combos: {0}".format(combos))
    attractors = []
    while len(combos) > 0:
        s = combos[0]
#        print(s)
        condensed_att = _condense_state([s], ruleset)
        condensed_att = list(dict.fromkeys(condensed_att))
        print("Got the entire attractor as {0}".format(condensed_att))
        for c_att in condensed_att:
            if c_att in combos:
#                print("removing {0} from states to consider".format(c_att))
                combos.remove(c_att)
#        print("Remaining states to consider: {0}".format(states_in_att))
#        input()
        attractors += [condensed_att]
    print("attractors: {0}".format(attractors))
    n_elements = len(ruleset[0][0])
    real_attractors = []
    for a in attractors:
        valid = True
        for s in a:
            if s == '*'*n_elements:
                valid =False
        if valid:
            real_attractors += [a]
    print("real_attractors: {0}".format(real_attractors))
    attractors = copy.deepcopy(real_attractors)
    real_attractors = []
    for a_1 in attractors: #Find if one of the states here is already in another state
        is_final = True
        for a_2 in attractors:
            if not a_1 == a_2:
                for s_1 in a_1:
                    for s_2 in a_2:
                        if _check_subset(s_2, s_1) and not s_1 == s_2:
                            print("{0} is subset of {1}".format(s_2, s_1))
#                            print("or")
#                            print("{0} is subset of {1}".format(s_1, s_2))
                            is_final = False
        if is_final:
            real_attractors += [a_1]
    print("realest attractors: {0}".format(real_attractors))
    if len(real_attractors) == 0:
        real_attractors = ['*'*n_elements]
    return real_attractors

def _condense_state(states, ruleset):
    """Go through ruleset to get all states.
    Since the states above must be in attractors, the states following them are also in attratcors.
    """
#    print("_condense_state")
#    print(states)
    #print(ruleset)
    next_states = []
    for state in states: #Go through each state here
        relevant_rules = []
##        print("Considering state {0}".format(state))
        rule_found = False
        for inp, out in ruleset: #Go thgourh all the rules
            if inp == state: #If there exists a rule for this state
                rule_found = True
#                print("Rule found for {0}".format(state))
#                print("{0} -> {1}".format(inp, out))
                if not out in states: #If the output of the rule is not considered
                    #print("Output deemed not to be in states")
                    s = states + [out]
                    states = _condense_state(s, ruleset)
#                    print(states)
                    next_states += [states]
                else: #The output of the rule is already considered.
             #       print("Output is already in states")
                    states = states
                    next_states += [states]
#        print(next_states)
        if not rule_found:
            print("Rule for {0} not found in ruleset".format(state))
            #Need second best?
            outstate = []
            for inp, out in ruleset:
                if _check_subset(state, inp): #If inp applies to state
                    outstate += [out]
                    print("{0} -> {1} fits".format(inp, out))
            print(outstate)
            resulting_outstate = outstate.pop()
            while len(outstate) > 0:
                resulting_outstate = _merge(outstate.pop(), resulting_outstate)
            #print(resulting_outstate)
            outstate = resulting_outstate
            print("Going with {0}".format(outstate))
            n_elements = len(states[0])
#            s = states + ['*'*n_elements]
            s = states + [outstate]
            #print('*'*n_elements)
            #input()
            return s
            raise Exception('Rule for state {0} not found in ruleset'.format(state))

    return states

def _expand_notation(attractor):
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
        att_0 = _expand_notation(att_0)
        att_1 = _expand_notation(att_1)
        #print("{1}: This became a leaf. Returning {0}".format(att_0+ att_1, attractor))
        return att_0+ att_1

def unmap(attractor, mapping):
    reverse_mapping = _reverse_map(mapping)
    n_elements = len(attractor)
    output = ''
    for i in range(n_elements):
        output += attractor[mapping[i]]
    return output

def _reverse_map(mapping):
    output = {}
    for key, value in mapping.items():
        output[value] = key
    return output

def _merge(temp_1, temp_2):
    """Merge two compatible sequences
    """
#    print(temp_1)
#    print(temp_2)
    n_symbols = len(temp_1)
    merged_temp = ['*']*n_symbols
    for i in range(n_symbols):
        s_1 = temp_1[i]
        s_2 = temp_2[i]
        if not s_1 == "*":
            merged_temp[i] = s_1
        if not s_2 == "*":
            merged_temp[i] = s_2
    merged_temp = ''.join(merged_temp)
#    print(merged_temp)
#    input()
    return merged_temp

def _check_compatibility(temp_1, temp_2):
    """Check if two input sequences could be merged.
    """
    compatible = True
    n_symbols = len(temp_1)
    for i in range(n_symbols):
        s_1 = temp_1[i]
        s_2 = temp_2[i]
        if not (s_1 == '*' or s_2 == '*'):
            if not s_1 == s_2:
                compatible = False
    return compatible

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

def _unasterisk(state):
    #print("working {0}".format(state))
    n_asterisks = _count_asterisks(state)
    output = []
    if n_asterisks:
        first_asterisk_index = _get_first_asterisk(state)
        state_0 = copy.deepcopy(list(state))
        state_1 = copy.deepcopy(list(state))
        state_0[first_asterisk_index] = '0'
        state_1[first_asterisk_index] = '1'
        state_0 = ''.join(state_0)
        state_1 = ''.join(state_1)
        state_0 = _unasterisk(state_0)
        state_1 = _unasterisk(state_1)
        output = state_0 + state_1
        #print("returning {0}".format(output))
        return output
    else:
        #print("returning {0}".format([state]))
        return [state]

def _get_first_asterisk(state):
    n_elements = len(state)
    for i in range(n_elements):
        if state[i] == '*':
            break
    return i
