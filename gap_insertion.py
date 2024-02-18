import numpy as np

class Gap_Insertion :
    def __init__(self) : 
        pass

    def deter_cons(self, list_compo) : 
        #Cette fonction permet de trouver quel est le nucléotide majoritaire à partir d'une liste de leur proportions. 
        COMPO = ["A", "C", "G", "T", "-", "N"]
        nucl = 0
        hprop = 0
        for prop in range(len(list_compo)) : 
            if list_compo[prop] > hprop : #Pour chaque nucléotide, on en regarde le nombre, le nucléotide le plus présent deviens le nucléotide consensus
                hprop = list_compo[prop]
                nucl = prop
        return COMPO[nucl]

    def align_oldcons_cons(self, cons, compopos) : 
        #Cette fonction permet d'ajouter les gaps de la séquence consensus dans les listes de composition de séquences
        compopos = np.array(compopos)
        while len(compopos[0]) != len(cons) : #On stocke les positions de changement et à chacun on rajoute un gap. On s'arrête que la composition et la consensus finale font la même taille.
            change = -1
            for i in range(len(compopos[0])) : 
                if self.deter_cons(compopos[:,i]) !=  "-" and cons[i] == "-" : 
                    change = i
            if change != -1 :
                compopos = np.insert(compopos, change, [0, 0, 0, 0, 2, 0], axis = 1)
        return compopos #La fonction renvoie la composition mise à jour.
    
    def fusion(self, list_compopos) : 
        list_pos = []
        list_pos.append(list_compopos[0])
        for i in range(len(list_compopos[1:len(list_compopos)])+1) : 
            for a in range(i) : 
                list_pos[a-1] += list_compopos[i][a-1]
        return list_pos
