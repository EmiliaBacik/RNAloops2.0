#!/usr/bin/python
import os
BP_O = "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ"
BP_C = ")]}>abcdefghijklmnopqrstuvwxyz"

def find_matching_parenthesis(text, opening_position):
    closing_position = opening_position
    counter = 1
    while counter > 0:
        closing_position += 1
        char = text[closing_position]
        if char == "(":
            counter += 1
        if char == ")":
            counter -= 1
    return closing_position

def find_matching_char(char, text, initial_position):
    print(char)
    openList = ["[","{","<","A","B","C","D","E"]
    closeList = ["]","}",">","a","b","c","d","e"]
    if char in openList:
        openChar = char
        closeChar = closeList[openList.index(openChar)]
        char_position = initial_position
        counter = 1
        while counter > 0:
            char_position += 1
            if char_position < 0 or char_position >= len(text):
                raise ValueError(
                    f"No matching bracket for '{char}' at position {initial_position}"
                )
            char = text[char_position]
            if char == openChar:
                counter += 1
            if char == closeChar:
                counter -= 1
    elif char in closeList:
        closeChar = char
        openChar = openList[closeList.index(closeChar)]
        char_position = initial_position
        counter = 1
        while counter > 0:
            char_position -= 1
            if char_position < 0 or char_position >= len(text):
                raise ValueError(
                    f"C No matching bracket for '{char}' at position {initial_position} cuz we're in {char_position} and we see {text[char_position]}"
                )
            char = text[char_position]
            if char == closeChar:
                counter += 1
            if char == openChar:
                counter -= 1
    else:
        char_position = initial_position

    return char_position

def fill_secondary(conectors, text):
    connectors_ss = []
    chars = ['.','(',')']
    for connector in conectors:
        if connector[1]-1 == connector[0]:
            #return connectors_ss, text
            connectors_ss.append('')
        else:
            for nucl_no in range(connector[0]+1,connector[1]):
                if text[nucl_no-1] not in chars:
                    paired_bracket_loc = find_matching_char(text[nucl_no-1], text, nucl_no-1)
                    if all(paired_bracket_loc not in range(pair[0],pair[1]) for pair in conectors):
                        text = text[:nucl_no-1] + '.' + text[nucl_no:]
                        text = text[:paired_bracket_loc] + '.' + text[paired_bracket_loc+1:]
            connectors_ss.append(text[connector[0]-1:connector[1]])
    return connectors_ss, text


def identify_loop_outermost_bps(structure_2d, bps_dict, outer_bps, mode, loop_outermost_bps):
    complexity = 0
    out_list = []
    for i, bp in enumerate(outer_bps):
        pairs = []
        pairs1 = []
        complexity = 0
        bpo = bpc = bp[0]
        wrong = False
        continous_fragments = []
        while(bpc + 1 != bp[1]):
            for bpc in range(bpo,len(structure_2d)):
                oidx = BP_O.find(structure_2d[bpc])
                cidx = BP_C.find(structure_2d[bpc])
                if (mode == 0 and (oidx == 0 or cidx == 0)) or (mode != 0  and ((oidx >= 0 and oidx <= bp[2]) or (cidx >= 0 and cidx <= bp[2]))):
                    break

            if bpc in bps_dict:
                fragment = (bpo, bpc)
                if fragment not in continous_fragments:
                    continous_fragments.append(fragment)
                    pairs1.append(sorted([fragment[0], bps_dict[fragment[0]-1]+1]))
                else:
                    wrong = True
                    break
                bpo = bps_dict[bpc] + 1
                pairs.append(sorted([bpc + 1, bps_dict[bpc] + 1]))

            else:
                wrong = True
                break
            complexity = complexity + 1


        if complexity >= 3 and not wrong:
            loop_outermost_bps.append((bp, complexity))
            out_list.append([complexity, pairs])
    return out_list, loop_outermost_bps

def check(list1, val):
    return(all(x < val for x in list1))

def identify_base_pairs(structure_2d, mode, bps, bps_dict):
    stacks = []
    for i in range(len(BP_O)):
        stacks.append([])
    for i, v in enumerate(structure_2d):
        oidx = BP_O.find(v)
        cidx = BP_C.find(v)
        if (mode == 0 and oidx == 0) or (mode != 0 and oidx >= 0):
            stacks[oidx].append(i)
        elif (mode == 0 and cidx == 0) or (mode != 0 and cidx >= 0):
            j = stacks[cidx].pop()
            bps.append((j+1,i+1,cidx))
            bps_dict[j] = i
            bps_dict[i] = j
    bps.sort(key = lambda x: x[0])

def takeFirst(elem):
    return elem[0]

def find_junction(db_sequence, common_stem_length_calc = True):
    mode = 1
    bps = []
    fragments = []
    bps_dict = {}
    bps_dict.clear()
    list_of_junctions =[]
    output = []
    identify_base_pairs(db_sequence, mode, bps, bps_dict)
    outer_bps = []
    identify_outer_bps(bps, outer_bps)
    loop_outermost_bps = []
    list_of_junctions, outermost_bps = identify_loop_outermost_bps(db_sequence, bps_dict, outer_bps, mode, loop_outermost_bps)
    order_list = []
    for bp in outermost_bps:
        ((i, j, order), complexity) = bp
        fragments.append(create_fragments(complexity, (i, j, order), mode, db_sequence, bps_dict))
        order_list.append(order)
    if list_of_junctions:
        for junction in list_of_junctions:
            junction[1].sort(key=takeFirst)
            output += [[*junction, common_stem_length(db_sequence, junction[1])]]
    for f in fragments:
        f.sort(key=takeFirst)
    #for f in fragments:
        #f = find_minimal(f)
    return output, fragments, order_list

def find_minimal(fragments):
    min_val = fragments[0][0]
    k = 0
    for x, el in enumerate(fragments):
        if el[0] < min_val:
            min_val = el[0]
            k = x

    for i in range(0, k):
        fragments.append(fragments.pop(0))
    return fragments

def create_fragments(complexity, outermost_bp, mode, structure_2d, bps_dict):
    bpo = outermost_bp[0]
    summary = []
    #print(f"Creating fragments for outermost_bp: {outermost_bp}, complexity={complexity}")
    while(complexity > 0):
        for bpc in range(bpo,len(structure_2d)):
            oidx = BP_O.find(structure_2d[bpc])
            cidx = BP_C.find(structure_2d[bpc])
            if (mode == 0 and (oidx == 0 or cidx == 0)) or (mode != 0 and ((oidx >= 0 and oidx <= outermost_bp[2]) or (cidx >= 0 and cidx <= outermost_bp[2]))):
                summary.append([bpo, bpc+1])
                #print(f"  Fragment added: start={bpo}, end={bpc+1}, oidx={oidx}, cidx={cidx}")
                complexity = complexity - 1
                break
        bpo = bps_dict[bpc]+1
    return summary

def find_junction_mode_single(db_sequence, junction_start, common_stem_length_calc = True):
    stems_identified = 0
    list_of_pairs = []
    if db_sequence[junction_start] == "(":
        last = find_matching_parenthesis(db_sequence, junction_start)
        list_of_pairs.append([junction_start+1, last+1])
        junction_start += 1
        stems_identified += 1

        while junction_start < last:
            if db_sequence[junction_start] == "(":
                end = find_matching_parenthesis(db_sequence, junction_start)
                list_of_pairs.append([junction_start+1, end+1])

                junction_start = end + 1
                stems_identified += 1

            else:
                junction_start += 1
    else:
        junction_start += 1
    if stems_identified > 2:
        if common_stem_length_calc:
            return True, stems_identified, list_of_pairs, common_stem_length(db_sequence, list_of_pairs)
        else:
            return True, stems_identified, list_of_pairs
    else:
        return False, stems_identified, list_of_pairs

def identify_outer_bps(bps, outer_bps):
    for i, v in enumerate(bps):
        if i == len(bps)-1:
            outer_bps.append(v)
        if i < len(bps)-1:
            (i1,j1,l1) = v
            (i2,j2,l2) = bps[i+1]
            if not ((i2 == i1 + 1) and (j2 == j1 - 1) and l1 == l2):
                outer_bps.append(v)

def common_stem_length(db_sequence, list_of_pairs):
    length_table = []
    pair_length_table = []
    for i,pair in enumerate(list_of_pairs):
        char_start = db_sequence[pair[0]-1]
        char_end = db_sequence[pair[1]-1]
        start = pair[0]-1
        open_count = 0
        end = pair[1]-1
        end_count = 0
        if i == 0:
            while db_sequence[start] == char_start and start >= 0:
                start -=1
                open_count += 1
        else:
            while db_sequence[start] == char_start:
                start +=1
                open_count += 1
        length_table.append(open_count)
        if i == 0:
            while db_sequence[end] == char_end and end <= len(db_sequence)-1:
                end +=1
                end_count += 1
        else:
            while db_sequence[end] == char_end:
                end -=1
                end_count += 1
        length_table.append(end_count)
        pair_length_table.append([open_count, end_count])
    pair_common_length_table = [min(pair) for pair in pair_length_table]
    return pair_common_length_table
