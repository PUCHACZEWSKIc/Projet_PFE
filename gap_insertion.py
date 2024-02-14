import parasail 

def old_cons_reconstruction(compopos) : 
    #Cette fonction permet de recréer en str la séquence consensus fournie par compopos (qui est une liste de liste qui contient la composition en acide aminé à chaque position)
    fin = ""
    for index in range(len(compopos)) : #Pour chaque index dans la séquence, on récupère la quantité de chaque nucléotide, et on détermine à l'aide de dter_cons le nucléotide consensus. 
        list_compo = [compopos[index][0],compopos[index][1],compopos[index][2],compopos[index][3],compopos[index][4],compopos[index][5]] 
        fin += deter_cons(list_compo) #On rajoute le nucléotide consensus à la séquence en str, et on la renvoie après la boucle. 
    return fin

def deter_cons(list_compo) : 
    #Cette fonction permet de trouver quel est le nucléotide majoritaire à partir d'une liste de leur proportions. 
    COMPO = ["A", "C", "G", "T", "-", "N"]
    nucl = 0
    hprop = 0
    for prop in range(len(list_compo)) : 
        if list_compo[prop] > hprop : #Pour chaque nucléotide, on en regarde le nombre, le nucléotide le plus présent deviens le nucléotide consensus
            hprop = list_compo[prop]
            nucl = prop
    return COMPO[nucl]

def align_oldcons_cons(cons, compopos) : 
    #Cette fonction permet d'ajouter les gaps de la séquence consensus dans les listes de composition de séquences
    old_cons = old_cons_reconstruction(compopos) 
    par = parasail.nw_trace_scan(old_cons, cons, 10, 1, parasail.dnafull) #On utilise parasail pour aligner les séquences consensus finale et intermédiaires, pour savoir où insérer nos gaps.
    traceback = par.get_traceback().query 
    while len(compopos) != len(traceback) : #On stocke les positions de changement et à chacun on rajoute un gap. On s'arrête que la composition et la consensus finale font la même taille.
        change = -1
        for i in range(len(compopos)) : 
            if deter_cons(compopos[i]) !=  "-" and traceback[i] == "-" : 
                change = i
        if change != -1 :
            compopos.insert(change-1, [0, 0, 0, 0, 2, 0])
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