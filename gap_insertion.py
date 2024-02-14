import numpy as np
import parasail 

def old_cons_reconstruction(compopos) : 
    fin = ""
    for index in range(len(compopos)) : 
        list_compo = [compopos[index][0],compopos[index][1],compopos[index][2],compopos[index][3],compopos[index][4],compopos[index][5]]
        fin += deter_cons(list_compo)
    return fin

def deter_cons(list_compo) : 
    COMPO = ["A", "C", "G", "T", "-", "N"]
    nucl = 0
    hprop = 0
    for prop in range(len(list_compo)) : 
        if list_compo[prop] > hprop : 
            hprop = list_compo[prop]
            nucl = prop
    return COMPO[nucl]

def align_oldcons_cons(cons, compopos) : 
    old_cons = old_cons_reconstruction(compopos)
    par = parasail.nw_trace_scan(old_cons, cons, 10, 1, parasail.dnafull) 
    traceback = par.get_traceback().query 
    while len(compopos) != len(traceback) : 
        chang = -1
        for i in range(len(compopos)) : 
            if deter_cons(compopos[i]) != traceback[i] : 
                change = i
        if change != -1 :
            compopos.insert(change, [0, 0, 0, 0, 2, 0])
    return compopos 





test = [[2, 0, 0, 0, 0, 0],
        [2, 0, 0, 0, 0, 0],
        [2, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 2, 0],
        [2, 0, 0, 0, 0, 0],
        [2, 0, 0, 0, 0, 0],
        [2, 0, 0, 0, 0, 0],
        [2, 0, 0, 0, 0, 0],
        [0, 2, 0, 0, 0, 0]]



cons = "AAA--AAAA--"
print(old_cons_reconstruction(test))
print(old_cons_reconstruction(align_oldcons_cons(cons, test)))