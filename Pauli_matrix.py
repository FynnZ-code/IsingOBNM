"""
    Programm to handle big sums of Pauli Matrices
"""


"""___________________ - Imports - _____________________________

"""



import numpy as np
import copy as copy
import time



"""___________________ - Global Variables - ________________________
    
"""

count_time = 0
sort_time = 0
var_time = 0



"""____________________ - Programm - _______________________________
    Classes
"""

class Pauli:
    """
        Pauli | represents a Pauli Matrix
            @site - Site it acts on (0-L)
            @kind - Kind of Pauli Matrix (0-3)
    """


    def __init__(self, site, kind):

        self.pauli = np.array([site, kind])

        return


    def get_string(self):
        """
            returns the pauli matrix in form of a string
        """
        
        return '\u03C3' + '(' + str(self.pauli[0]) + ',' + str(self.pauli[1]) + ') '



class Term:
    """
        Represents a Term in a Sum
            @matrices - List of Pauli matrices in the Term
            @coeficient - coeficient of the Term
    """

    def __init__(self, matrices, coeficient = 1):
        """
            Term initializer
                - coeficient has default value 1 
        """
        self.term = np.array([[coeficient],matrices], dtype=object)
        self.sorted = False

        return
    

    def append(self, liste):
        """
            appends the given list of pauli matrices onto the instance
        """
        self.sorted = False
        return Term(np.append(self.term[1], liste))



    def get_coeficient_shift(self, pos_aim, pos_mat):       # | NOT IN USE
        """
            calculates the sign of the coeficient when shifting the position of pauli matrix at 
            pos_mat to pos_aim
        """
       
        distance = abs(pos_aim - pos_mat)
        counter = 0

        if pos_aim > pos_mat:
            
            for iii in range(distance): # counts the number of paulis with the same site in between
                if self.matrices[pos_mat + iii + 1].site == self.matrices[pos_mat].site:
                    counter += 1

        else:
            
            for iii in range(distance): # counts the number of paulis with the same site in between
                if self.matrices[pos_aim + iii].site == self.matrices[pos_mat].site:
                    counter += 1 
        
        return  (-1)**counter


    
    @classmethod 
    def sort(cls, term):
        Flag = True
        
        matrices = term.term[1]
        coeficient = term.term[0][0]

        while Flag:
            Flag = False
            i = 0

            
            while i < (len(matrices) - 1):
                site1, kind1, site2, kind2 = matrices[i].pauli[0], matrices[i].pauli[1], matrices[i+1].pauli[0], matrices[i+1].pauli[1]
                
                
                
                
                if site1 == site2 and kind1 == kind2:
                    
                    matrices = Term.fold(i, i+1 , matrices)
                    

                    Flag = True

                elif site1 > site2 or (site1 == site2 and kind1 > kind2):
                    
                    matrices, coeficient = Term.swap(i, i+1, matrices, coeficient)
                    
                    
                    Flag = True
                    i+= 1

                else:
                    i += 1
                
                
            
            
                

        sorted = True
        return Term(matrices, coeficient=coeficient)

    @classmethod
    def fold(cls, i, j, matrices):
        
        return np.concatenate((matrices[:i],matrices[j+1:]))

    @classmethod
    def swap(cls, i, j, matrices, coeficient):
        matrices[i], matrices[j] = matrices[j], matrices[i]

        if matrices[i].pauli[0] == matrices[j].pauli[0]:
            coeficient *= -1

        return matrices, coeficient


    def print(self, printe = True):
        """
            Creates a string that represents the Term
            If printe is True the function prints the sting itself
        """

        string1 = ''
        string1 = string1 + str(self.coeficient)

        for matrix in self.matrices:
            string1 = string1 + matrix.get_string()

        if printe:
            print(string1)

        return string1


    def compare(self, term):            # | NOT IN USE
        """
            compares the given term with self
            if they are the same return True
        """

        if not len(self.matrices) == len(term.matrices):
            return False

        for i in range(len(self.matrices)):
            if not (self.matrices[i].kind == term.matrices[i].kind and self.matrices[i].site == term.matrices[i].site):
                return False

        return True

    def to_string(self):
        string = ''
        for matrix in self.term[1]:

            string = string + str(matrix.pauli[1]) +  ',' + str(matrix.pauli[0]) + ';'


        return string


    @classmethod
    def str_to_term(cls, str, coef):
        """
            Converts the given string into a term object and returns it
        """
        
        if str == '':
            return Term([], float(coef))
        
        matrices = str.split(';')
        matrices.pop()
        new_term = Term([], float(coef))

        if str == '':
            return new_term

        for matrix in matrices:
            new_term.term[1] = np.append(new_term.term[1], Pauli(int(matrix.split(',')[1]), int(matrix.split(',')[0])))

        return new_term



class Sum:
    """
        Represents a Sum of Terms
            @terms - List of terms in the sum
    """


    def __init__(self, terms):

        self.terms = terms

        return


    def set_coeficients(self, coeficient, alpha = 1):
        """
            Sets the coeficient of all terms in the sum to
                coeficient * alpha
            default  alpha = 1
        """

        for i in range(len(self.terms)):

            self.terms[i].term[0][0] = self.terms[i].term[0][0] * coeficient * alpha

        return


    def get_coeficients(self):
        """
            returns the coeficients of all terms in the sum as array
        """

        coeficients = np.array([])

        for term in self.terms:
            coeficients = np.append(coeficients, term.term[0][0])
        
        return coeficients
    


    def Hash_count(self, L, sort = True):
        
        count_st = time.time()

        hash_map = {}
        i = 0


        self.terms = [Term.sort(term) for term in self.terms if abs(self.terms[i].term[0][0]) != 0 and sort and not term.sorted]

    
        
        while i<len(self.terms):

            key = self.terms[i].to_string()
            value = self.terms[i].term[0][0]

            

            if key in hash_map:
                hash_map[key] += value

            else:
                hash_map[key] = value

            i += 1
        
        

        counted_terms = [Term.str_to_term(key, value) for key, value in hash_map.items() if value != 0]
        


        self.terms = counted_terms

        count_et = time.time()
        global count_time
        count_time += count_et - count_st

        return


    def sort(self, L):
        """
            calls the sort method for all terms in the sum
        """

        for i in range(len(self.terms)):
            self.terms[i].sort(L)

        return



    def print(self, printer = True):
        """
            Generates the sum as a string 
            if printer True the function prints out the string itself
        """

        string_sum = ''

        for term in self.terms:
            coeficient = term.term[0][0]
            matrices = term.term[1]
            string_term = str(coeficient)

            for matrix in matrices:
                string_term = string_term + matrix.get_string()

            if len(string_sum) == 0:
                string_sum = string_sum + string_term

            else:
                string_sum = string_sum + ' + ' + string_term

        if printer:
            print(string_sum)

        return string_sum



class Product:
    """
        Represents a Product of sums
            @sums - List of sums in the Product
    """



    def __init__(self, sums):

        self.sums = sums

        return


    def execute(self):
        """
            Executes the Product
            returns the resulting sum
        """
        
        sum = Sum([])
        count = 0

        for term1 in self.sums[0].terms:
            for term2 in self.sums[1].terms:

                sum.terms.append(term1.append(term2.term[1]))
                sum.terms[count].term[0][0] = term1.term[0][0] * term2.term[0][0]
                count += 1

        return sum



class System:
    """

    """

    def __init__(self, L, q, sys):
        

        self.L = L
        self.q = q

        self.sys = sys

        return

    

    def Hamiltonian(self):
        """
            Generates the Hamiltonian term of a 1-d Ising chain with a magnetic field in x direction of the length L
        """


        Hamiltonian = Sum(np.array([]))

        for l in range(self.L):
            
            if l+1 == self.L:
                Hamiltonian.terms= np.append(Hamiltonian.terms, np.array([Term([Pauli(l,3),Pauli(0,3)]), Term([Pauli(l, 1)]), Term([Pauli(l, 3)], coeficient= 0.5)], dtype=object))

            else:
                Hamiltonian.terms= np.append(Hamiltonian.terms,np.array([Term([Pauli(l,3),Pauli(l+1,3)]), Term([Pauli(l, 1)]), Term([Pauli(l, 3)], coeficient= 0.5)], dtype=object))


        if not(self.sys == None):
            
            E_max = max(self.sys.Hamiltonian.Ew)
            E_min = min(self.sys.Hamiltonian.Ew)

            a = (E_max-E_min)/2
            b = (E_max+E_min)/2
            
            b_dummy = copy.deepcopy(Hamiltonian)
            b_dummy.set_coeficients
            Hamiltonian.terms = np.append(Hamiltonian.terms, np.array([Term([], coeficient = -b)]))
            Hamiltonian.set_coeficients(1/a)

        return Hamiltonian



    def Observable(self):
        """
            returns the Pauli Operator of the Observable
        """

        Observable = Sum([Term([Pauli(self.q,3)])])

        return Observable



    def POP(self, a, power = 2):
        """
            Generates all terms of the projected Operator POP
        """
        H_filter = self.Hamiltonian()

        for i in range(power - 1):
            H_filter = Product([H_filter,self.Hamiltonian()]).execute()

        #Hsquared = Product([self.Hamiltonian(), self.Hamiltonian()])
        Observable = self.Observable()

        H_filter.Hash_count(self.L)
        
        
        POP_symbol = 'POP'
        i = 0

        POP = Sum([Term([])])
        
        while i < 3:
            
            if POP_symbol[i] == 'P':

                P = Product([POP, H_filter]).execute()
                P.set_coeficients(-1, a)
                P.Hash_count(self.L)

                POP.terms = np.append(POP.terms, P.terms)


            elif POP_symbol[i] == 'O':

                POP = Product([POP, Observable]).execute()
            
            POP.Hash_count(self.L)
            i += 1

        
        return POP

    

    def POP_squared(self, a, power = 2):
        """ recycle second mom 
            Generates all terms of the squared projected Operator (POP)^2
        """



        H_filter = self.Hamiltonian()

        for i in range(power - 1):
            H_filter = Product([H_filter,self.Hamiltonian()]).execute()

        #Hsquared = Product([self.Hamiltonian(), self.Hamiltonian()])
        Observable = self.Observable()

        H_filter.Hash_count(self.L)
        

        POP_symbol = 'POP'
        i = 0
        flag = False

        POP_squared = Sum([Term(np.array([]))])
        
        while True:


            if i == 3:

                if flag:
                    return POP_squared

                else:

                    i=0
                    flag = True


            if POP_symbol[i] == 'P':

                P = Product([POP_squared, H_filter]).execute()
                P.set_coeficients(-1, a)
                P.Hash_count(self.L)

                POP_squared.terms = np.concatenate((POP_squared.terms , P.terms))

            
            elif POP_symbol[i] == 'O':

                POP_squared = Product([POP_squared, Observable]).execute()
            

            POP_squared.Hash_count(self.L)

            i += 1

        
    def count_terms(self, O):
        """
            Calls counting Method on O
        """
        return O.Hash_count(self.L)


    def Moment(self, O):
        """
            calculates the second moment through the metric
            of POP in the Operator Basis
        """

        coeficients = O.get_coeficients()
        sum = 0

        for coeficient in coeficients:
            sum += coeficient**2

        return sum




"""____________________ - Programm - _______________________________
    Functions
"""


def analytical_solution(L, alpha):


    return alpha *(48*L**2-80*L+48)
