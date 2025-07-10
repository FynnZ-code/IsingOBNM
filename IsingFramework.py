import numpy as np
from numpy import array, linalg as LA, matrix, real
import scipy as sp
import matplotlib.pyplot as plt
import random as random
import copy as copy
import math


class Hamiltonian:
    """
    Klasse Hamiltonian.
        @Parameter      Hält Koeffizienten der verschiedenen Terme des Hamiltonians

        @__innit__      berechnet bei übergabe der Spin Basis und der Kettenlänge die Matrixstrucktur des Hamiltonians
                        Außerdem werden die Energieeigenwerte und Zustände berechnet und als Attribute der Klasse abgespeichert

        @zeromatrix     Erzeugt Matrix mit Nullen gefüllt mit Dimension des Hamiltonians

        @in_energy_basis Führt eine Basis transformatione der übergebenen Matrix in die Energie ecmigenbasis durch
    """


    J = 1
    hx = 1
    hz = 0.5

    lift_defects = False
    h_2 = 0.165
    h_5 = -0.24


    def __init__(self, basis, L, diagonalize = True, psi1 = None, psi2 = None):
        """
        Konstruktor:
            - Festlegen der Parameter Kettenlänge, Dimension und Hamiltonian matrix
            - Iteriert über die Matrixelemente des Hamiltonians und darin über die gesamte Länge der Kette
            - Berechnet dann das Matrixelement durch summieren über die Ergebnisse aus den Methoden die die Terme berechnen
            - Berechnen der Eigenvektoren und werte des Hamiltonians und speichert diese als Attribute ab
        """


        self.d = len(basis)
        self.L = L
        self.H = self.zeromatrix()




        for i in range(self.d):
            for ii in range(self.d):


                for l in range(self.L):

                    self.H[i][ii] += self.J * self.zz(basis[i], basis[ii],l)
                    self.H[i][ii] += self.hx * self.xmagnet(basis[i], basis[ii],l)
                    self.H[i][ii] += self.hz * self.zmagnet(basis[i], basis[ii],l)

                if self.lift_defects:
                    self.H[i][ii] += self.h_2 * self.defect_tran(basis[i], basis[ii])
                    self.H[i][ii] += self.h_5 * self.defect_par(basis[i], basis[ii])

                if not(psi1 is None or psi2 is None):
                    self.H[i][ii] = psi1[i]*psi2[ii]*self.H[i][ii]
        
        if diagonalize:
            w, v = sp.linalg.eigh(self.H)
        else:
            w, v = None, None


        self.Ev = v
        self.Ew = w

        

        return



    def zz(self, state1, state2, l):
        """
        ZZ - Term der Ising Chain
            - Berechnet nächste Nachbar wechselwirkung des Z-Spins auf einer Ising chain
            - Perioische Randbedingung
                -> Fallunterscheidung zwischen Bulk und Rand
        """
        error = 0

        for site in range(self.L):
            if state1[site] != state2[site]:
                error += 1



        if error == 0:
            if l+1 == self.L:  # fallunterschiedung periodische randbedingung

                matrixelement = (2*state2[l]-1)*(2*state2[0]-1)


            else:
                matrixelement = (2*state2[l]-1)*(2*state2[l+1]-1)


            return matrixelement



        return 0.0


    def xmagnet(self, state1, state2,l):
        """
        Berechnet den Term der Wechselwrikung der Spins mit einem Magnetfeld in X-Richtung
            - Perioische Randbedingung
                -> Fallunterscheidung zwischen Bulk und Rand
        """

        sum = 0
        error = 0
        flag = False



        for i in range(self.L):
            if not state1[i] == state2[i]:

                if i == l:
                    flag = True

                else:
                    error += 1



        if error == 0 and flag:
            sum += 1


        return sum



    def zmagnet(self, state1, state2, l):
        """
        Berechnet den Term der Wechselwrikung der Spins mit einem Magnetfeld in Z-Richtung
            - Perioische Randbedingung
                -> Fallunterscheidung zwischen Bulk und Rand
        """

        sum = 0
        error = 0

        for site in range(self.L):
            if state1[site] != state2[site]:
                error += 1


        if error == 0:
            sum += 2*state2[l] - 1



        return sum


    def defect_tran(self, state1, state2):

        error = 0


        for site in range(self.L):
            if state1[site] != state2[site]:
                error += 1

        if error == 0:
            matrixelement = (2*state2[1]-1)


            return matrixelement            

        return 0.0


    def defect_par(self, state1, state2):

        error = 0


        for site in range(self.L):
            if state1[site] != state2[site]:
                error += 1

        if error == 0:
            matrixelement = (2*state2[4]-1)


            return matrixelement            

        return 0.0



    def zeromatrix(self):
        """
            Setzt ein Array auf, welches vollständig mit Nullen gefüllt ist und die Dimensionen des Hamiltonians hat
        """

        M = []

        for i in range(self.d):
            M.append([])

            for ii in range(self.d):
                M[i].append(0)

        return M



    def in_energy_basis(self, O):
        """
            Transformiert die Übergebene Matrix O in die Energie Eigenbasis
        """
        
        trafo1 = self.Ev
        trafo2 = np.linalg.inv(self.Ev)

        O_transformed = np.dot(np.dot(trafo2, O),trafo1)
       

        
        return O_transformed


    def scale_hamilton(self):
        
        E_max = max(self.Ew)
        E_min = min(self.Ew)


        a = (E_max-E_min)/2
        b = (E_max+E_min)/2

        self.H = np.dot(1/a,np.array(self.H) - np.identity(self.d)*b)

        return

    def rev_scale_hamilton(self, E):

        E_max = max(self.Ew)
        E_min = min(self.Ew)


        a = (E_max-E_min)/2
        b = (E_max+E_min)/2
        

        rev_E = E*a 

        return rev_E

    def scale_Eigenvalues(self):


        E_max = max(self.Ew)
        E_min = min(self.Ew)


        a = (E_max-E_min)/2
        b = (E_max+E_min)/2
        
        for w in range(self.d):
            self.Ew[w] = (self.Ew[w] - b)/a

        return   
    
    def scaled_Eigenvalues(self):


        E_max = max(self.Ew)
        E_min = min(self.Ew)


        a = (E_max-E_min)/2
        b = (E_max+E_min)/2
        #b = 0
        E_scaled = []
        for w in range(self.d):
            E_scaled.append((self.Ew[w] - b)/a)

        return  E_scaled



class Observable:



    def __init__(self, Obs, basis, H, q , p, energy_basis = True):
        self.L = H.L
        self.Obs = Obs


        if Obs == 'A':
            self.Matrix = self.AObserv(q, basis, H)

        elif Obs == 'B':
            self.Matrix = self.BObserv(basis, H)

        elif Obs == 'C':
            self.Matrix = self.CObserv(p, basis, H)

        elif Obs == 'T':
            self.Matrix = self.TObserv(basis, H)


        if energy_basis:
            self.Matrix = H.in_energy_basis(self.Matrix)

        self.Matrix = np.array(self.Matrix)
        return


    def AObserv(self, q, basis, H):
        """
            Observable  A
        """

        A = H.zeromatrix()



        for i in range(len(basis)):
            for ii in range(len(basis)):



                value = 0
                for l in range(self.L):
                    sum = 0

                    sum += H.J * H.zz(basis[i], basis[ii],l)
                    sum += H.hx * H.xmagnet(basis[i], basis[ii],l)
                    sum += H.hz * H.zmagnet(basis[i], basis[ii],l)
                    value += sum*(1/np.sqrt(self.L))*np.cos((2*np.pi/self.L)*q*l)


                A[i][ii]= value


        return A


    def BObserv(self, basis ,H ):
        """
            Observable B
        """
        B = H.zeromatrix()


        for i in range(len(basis)):
            for ii in range(len(basis)):
                sum = 0

                for l in range(self.L): # summe in observable
                    error = 0
                    flag = False


                    for ll in range(self.L): # checken der gleichheit der Zustände ausserhalb der Position l
                        if not basis[i][ll] == basis[ii][ll]:

                            if ll == l:
                                flag = True

                            else:
                                error +=1



                    if error == 0 and flag:
                        sum += 1


                B[i][ii] = (1/(np.sqrt(self.L)) *sum)


        return B


    def CObserv(self, p, basis, H):

        C = H.zeromatrix()

        for state in range(len(basis)):
            C[state][state] += (2*basis[state][p]-1)

        return C

    
    def TObserv(self, basis, H):
        
        T = H.zeromatrix()

        for state1 in range(len(basis)):
            for state2 in range(len(basis)):
                error = 0   

                for l in range(self.L):
                    if l == (self.L-1):
                        if not(basis[state1][l] == basis[state2][0]):
                            error += 1
                    else:
                        if not(basis[state1][l] == basis[state2][l+1]):
                            error += 1
                
                if (error == 0):
                    T[state1][state2] = 1
                else:
                    T[state1][state2] = 0
                
        return T
                 
                
                

        return



    def second_moment(self, H, E0, delta):
        sum = 0
        d = H.d

        E = H.Ew

        start = 0
        end = 0

        for i in range(len(E)):
            if E[i] >= (E0-delta/2):
                start = i
                break
            
        for i in range(len(E)):
            end = i
            if E[i] >= (E0 + delta/2):
                end = i-1
                break
        
        E = E[start:end]
        O = self.Matrix[start:end,start:end]
        d = len(E)
       
        
        sum = np.trace(O @ O)

        return  sum/d


    def second_moment_a(self, H, E0, delta):

        sum = 0
        d = H.d

        for n in range(d):

            En = H.Ew[n]
            for m in range(d):
                Em = H.Ew[m]

                if (En<(E0+delta) and En>(E0-delta)) and (Em<(E0+delta) and Em>(E0-delta)):
                    sum += self.Matrix[n][m] * self.Matrix[m][n]


        return sum/d

    def fourth_moment(self, H, E0, delta):
        sum = 0
        d = H.d

        E = H.Ew

        start = 0
        end = 0

        for i in range(len(E)):
            if E[i] >= (E0-delta/2):
                start = i
                break
            
        for i in range(len(E)):
            if E[i] >= (E0 + delta/2)or i == len(E)-1:
                end = i
                break
        
        
        E = E[start:end]
        d = len(E)
        O = self.Matrix[start:end,start:end]

        sum = np.trace(O @ O @ O @ O)

        return  sum/d


    def second_moment_approx(self, H, E0, alpha):
        sum = 0
        d = H.d

        for n in range(d):

            En = H.Ew[n]
            for m in range(d):
                Em = H.Ew[m]
                P = 1 - alpha * En**2
                P_squared = 1 - 2*alpha*Em**2 + alpha**2 * Em**4
                sum += P * self.Matrix[n][m] * P_squared * self.Matrix[m][n] * P

            


        return sum/d

    def fourth_moment_approx(self, H, E0, alpha):
        sum = 0
        d = H.d

        for n in range(d):

            En = H.Ew[n]
            for m in range(d):
                Em = H.Ew[m]
                for o in range(d):
                    Eo = H.Ew[o]
                    for p in range(d):
                        Ep = H.Ew[p]

                        P = 1 - alpha * En**2
                        P_squared_I = 1 - 2*alpha*Em**2 + alpha**2 * Em**4
                        P_squared_II = 1 - 2*alpha*Eo**2 + alpha**2 * Eo**4
                        P_squared_III = 1 - 2*alpha*Ep**2 + alpha**2 * Ep**4
                        sum += P * self.Matrix[n][m] * P_squared_I * self.Matrix[m][o]*P_squared_II*self.Matrix[o][p]*P_squared_III*self.Matrix[p][n] * P

    

        return sum/d



    def second_moment_approx_matrix(self, H, alpha, power = 2):

        H_filter = H.H
        
        for i in range(power - 1):
            H_filter = H_filter @ H.H

        P = np.identity(H.d) - np.dot(alpha , H_filter)
        
        POP = (P @ np.array(self.Matrix)) @ P

        POP_squared = POP @ POP


        return POP_squared
    
    def fourth_moment_approx_matrix(self, H, alpha, power = 2):

        
        POP_squared = self.second_moment_approx_matrix(H, alpha, power= power)

        POP_fourth = POP_squared @ POP_squared
        

        return POP_fourth
    
    
    
    def second_moment_approx_matrix_gaussian(self, H, sigma):

        H_filter = H.Ew#np.diag(H.in_energy_basis(H.H)) old mit rescaled

        
        P = []
        for h in H_filter:
            sum = 0
            for n in range(1):
                sum+= 1/sp.special.factorial(n) *(- 1/2*(h/sigma)**2)**n
            
            P.append(sum)
        P1 = np.diag(np.exp(- (1/2)* (H_filter/sigma)**2))

    

        

        POP = (P1 @ np.array(self.Matrix)) @ P1

        eff_d = np.trace(P1)

        POP_squared = POP @ POP


        return POP_squared , eff_d
    

    def fourth_moment_approx_matrix_gaussian(self, H, sigma):

        
        POP_squared, dummy = self.second_moment_approx_matrix_gaussian(H, sigma)

        POP_fourth = POP_squared @ POP_squared
        

        return POP_fourth
    
    

    def second_moment_approx_fit_gauss(self, H, sigma, depth):

        H_filter = H.Ew
        P = []

        steps = np.arange(min(H_filter),max(H_filter), 0.01 )
        gaus = np.exp(- (1/2)* (steps/sigma)**2)
        poly = np.polyfit(steps, gaus, deg = depth)



        P = np.diag(np.polyval(poly,H_filter))

        POP = (P @ np.array(self.Matrix)) @ P

        POP_squared = POP @ POP

        eff_delta =2* 0.674490 * sigma
        eff_d = np.trace(P)
        print("eff_d: " + str(eff_d))

        return POP_squared, eff_delta, eff_d
    
    
    def fourth_moment_approx_fit_gauss(self, H, sigma, depth):

        
        POP_squared, dummy, dummy_1 = self.second_moment_approx_fit_gauss(H, sigma, depth)

        POP_fourth = POP_squared @ POP_squared
        

        return POP_fourth
    

    def second_moment_approx_fit_step(self, H, d, k, depth):

        H_filter = H.Ew
        P = []

        # Poly fit to window function
        steps = np.arange(min(H_filter),max(H_filter), 0.01 )
        step = step_func(steps, d, k)
        poly = np.polyfit(steps, step, deg = depth)

        # Calculate standard deviation
        P = np.polyval(poly,steps)
        eff_delta = np.sqrt(np.var(P))
        eff_delta = abs(steps[np.argmin(abs(eff_delta-P))])*2
            
        # rescaling to peak 1 value
        P = np.polyval(poly,H_filter)
        
        one_rescaling = max(P)
        if one_rescaling <1: 
            P = P* 1/one_rescaling


        P = np.diag(P)
        # calculate POP^2
        POP = (P @ np.array(self.Matrix)) @ P
        POP_squared = POP @ POP
        eff_d = np.trace(P)

        return POP_squared, eff_delta, eff_d
    

    def fourth_moment_approx_fit_step(self, H, d, k, depth):

        
        POP_squared, dummy, dummy_1= self.second_moment_approx_fit_step(H, d, k, depth)

        POP_fourth = POP_squared @ POP_squared
        

        return POP_fourth
    
    
    


    


    





    def fourth_moment_a(self, H, E0, delta):

        sum = 0
        d = H.d

        for n in range(d):
            En = H.Ew[n]

            for m in range(d):
                Em = H.Ew[m]

                for o in range(d):
                    Eo = H.Ew[o]

                    for p in range(d):

                        Ep = H.Ew[p]

                        if (En<(E0+delta) and En>(E0-delta)) and (Em<(E0+delta) and Em>(E0-delta)) and (Eo<(E0+delta) and Eo>(E0-delta)) and (Ep<(E0+delta) and Ep>(E0-delta)):
                            sum += self.Matrix[n][m] * self.Matrix[m][o] * self.Matrix[o][p] * self.Matrix[p][n]


        return sum/d


    def average(self, H, E0, delta):
        avg = 0

        sum = 0
        for n in range(H.d):
            if H.Ew[n] > E0-delta and H.Ew[n] < E0 + delta:
                sum += 1
                avg += self.Matrix[n][n]

        avg = avg/sum

        return avg/sum


    def center(self, H, E0, delta):

        O_centered = H.zeromatrix()
        O_avg = self.average(H, E0, delta)

        for n in range(H.d):
            for m in range(H.d):
                O_centered[n][m] = self.Matrix[n][m] - O_avg

        self.Matrix = O_centered

        return



class System:

    def __init__(self, L, Obs, q = 0, p = 0, energy_basis = True, symmetry_subspace = False, symmetry_value = 5):

        self.L = L
        self.d = 2**L

        self.basis = self.base()
        
        self.Hamiltonian = Hamiltonian(self.basis, self.L)

        if symmetry_subspace:
            self.symmetry_subspace(symmetry_value)
            self.d = len(self.basis)


        print(f"Energie Eigenbasis: f{self.Hamiltonian.Ew}")
        self.Hamiltonian.scale_hamilton()


        if not symmetry_subspace:
            self.Observable = Observable(Obs, self.basis, self.Hamiltonian, q, p, energy_basis)
        
        return
    
    


    def symmetry_subspace(self, symmetry_value):
        T_obs = Observable('T', self.basis, self.Hamiltonian, q = 0, p = 0, energy_basis = False)

        w, v = sp.linalg.eig(T_obs.Matrix)
        w = w.real
        v = v.real
        print(w)
        indices = []
        subbasis = []

        for i in range(len(w)):
            if w[i] < symmetry_value+1 and w[i] > symmetry_value-1:
                indices.append(i)
        
        for i in range(len(indices)):
            subbasis.append(v[i])

        self.Hamiltonian.d = len(subbasis)
        new_H = self.Hamiltonian.zeromatrix()

        for i in range(len(subbasis)):
            for ii in range(len(subbasis)):
                state1 = subbasis[i]
                state2 = subbasis[ii]
                
                parameter_matrix = np.outer(state1, state2)
    
                new_H[i][ii] = np.multiply(parameter_matrix, self.Hamiltonian.H).sum().sum()

        
        self.Hamiltonian.H = new_H
        w, v = np.linalg.eigh(self.Hamiltonian.H)
        self.Hamiltonian.Ew = w
        self.Hamiltonian.Ev = v
        

        self.basis = subbasis
        return 


    def base(self):
        """
        Erzeugt eine Spin Basis zu einer 1-D Spin Chain der Länge L
        """

        basis = []


        for i in range(self.d): # Iteriert über alle dimensionen

            basis.append([])

            composition = bin(i)[2:] # Alle binären zahlen von 0-d gibt jede mögliche kombination von basis zuständen


            if len(composition) != (self.L): # Füllt binäre zahlen zu der länge L mit Nullen auf
                composition = composition.zfill(self.L)


            for ii in range(self.L): # Setzt einen Basiszustand aus einteilchen Zuständen gemäß komposition zusammen
                basis[i].append(int(composition[ii]))


        return basis



    def autocorrelationf(self, E0, delta, t):
        """
            Berechnung der Autokorrelationsfunktion der Observablen O in der Energei eigenbasis
            E0 und delta geben das Energie Fenster an aus dem die Energie Eigenzustände gewählt werden
            t ist der Zeitpunkt
        """

        sum = 0
        d = self.Hamiltonian.d


        for n in range(d):
            En = self.Hamiltonian.Ew[n]

            for m in range(d):

                Em = self.Hamiltonian.Ew[m]
                w = En-Em

                if abs(w) < 0.0001:
                    w = 0

                if (En<(E0+delta) and En>(E0-delta)) and (Em<(E0+delta) and Em>(E0-delta)):
                    sum += np.cos((w)*t) * self.Observable.Matrix[n][m] * self.Observable.Matrix[m][n]


        return sum/d



    def auto_print(self, E0, delta, t):

        array = []
        autoNull = self.autocorrelationf(E0, delta, 0)

        for time in t:
            array.append(self.autocorrelationf(E0, delta, time)/autoNull)

        sum = 0
        for element in array:
            sum += element
        time_avg = sum/len(array)

        
        t_th = self.thermalization_time(array, t)

        plt.figure()
        plt.plot(t, array, color = '#00aedb',linewidth = 0.8)
        plt.vlines(t_th, -0.1, 1.1,colors=['#8ec127'], linestyles='dashed', label= r'$\tau_{th}$' )
        plt.xlabel('t')
        plt.ylabel('C(t)/C(0)')
        plt.legend()
        plt.show()
        print('thermalization time: ' + str(t_th))
        return

    def thermalization_time(self, array,  time_ref):
        time = np.arange(0, len(time_ref))
        thermalization_time = [t for t in time if(array[t] - array[-1]) / (array[0] - array[-1]) < 0.01 ]

        return thermalization_time[0]


    def lambda_print(self, E0, delta):

        lambdas = []

        for d in delta:
            lambdas.append(self.lambda_calc(E0,  d))
        return lambdas



    def lambda_calc(self, E0, a):
        lambdat = 0
        
        
        sec = self.Observable.second_moment(self.Hamiltonian, E0, a)
        fou = self.Observable.fourth_moment(self.Hamiltonian, E0, a)
        
        if sec == 0 or fou == 0:
            lambdat = 0
        else:
            lambdat = sec**2/fou

        return lambdat
    


    def lambda_calc_step_fit(self, E0, delta, k = 200, depth = 60):
        lambdat = 0
        
        
        sec, eff_delta, eff_d = self.Observable.second_moment_approx_fit_step(self.Hamiltonian, delta, k, depth)
        fou = self.Observable.fourth_moment_approx_fit_step(self.Hamiltonian, delta, k, depth)

        d = eff_d

        
        secm = np.trace(sec)/d
        foum = np.trace(fou)/d


        if secm == 0 or foum == 0:
            lambdat = 0
        else:
            lambdat = secm**2/foum

        
        return lambdat, eff_delta



    def lambda_calc_alt(self, a, power = 2):
        lambdat = 0
        
       
        sec = self.Observable.second_moment_approx_matrix(self.Hamiltonian, a, power = power)
        fou = self.Observable.fourth_moment_approx_matrix(self.Hamiltonian, a, power = power)

        secm = np.trace(sec)/self.d
        foum = np.trace(fou)/self.d


        if secm == 0 or foum == 0:
            lambdat = 0
        else:
            lambdat = secm**2/foum

        eff_delta = 0
        return lambdat, eff_delta
    


    def lambda_calc_alt_gaussian(self, delta):
        lambdat = 0

        g = (delta)/(2 * 0.674490)
        
       
        sec, eff_d = self.Observable.second_moment_approx_matrix_gaussian(self.Hamiltonian, g)
        fou = self.Observable.fourth_moment_approx_matrix_gaussian(self.Hamiltonian, g) 
        

        eff_delta =2 * g * 0.674490

        secm = np.trace(sec)/eff_d
        foum = np.trace(fou)/eff_d


        if secm == 0 or foum == 0:
            lambdat = 0
        else:
            lambdat = secm**2/foum

        
        
        return lambdat, eff_delta


    def lambda_calc_fit_gauss(self, delta, depth):
        lambdat = 0


        g = (delta)/(2* 0.674490)#* a)

       
        sec, eff_delta, eff_d_1 = self.Observable.second_moment_approx_fit_gauss(self.Hamiltonian, g, depth)
        fou = self.Observable.fourth_moment_approx_fit_gauss(self.Hamiltonian, g, depth) 
        
        
        E_0=0
        eff_d = 0
        
        for e in self.Hamiltonian.Ew:
            if (e > E_0- eff_delta) and (e < E_0 + eff_delta):
                eff_d += 1
        


        secm = np.trace(sec)/eff_d_1 
        foum = np.trace(fou)/eff_d_1 


        if secm == 0 or foum == 0:
            lambdat = 0
        else:
            lambdat = secm**2/foum

        
        
        return lambdat, eff_delta
    

    

    def get_effective_delta(self, a , b , c , d, power):

        H_filter = np.diag(self.Hamiltonian.in_energy_basis(self.Hamiltonian.H))

        E = self.Hamiltonian.Ew

        values = [[], []]

        for i in range(len(E)):

            values[0].append(E[i])
            values[1].append(1-(a*H_filter[i]**2 +b*H_filter[i]**4 +c*H_filter[i]**6 + d*H_filter[i]**8))
        
        sum = 0
        for i in range(len(E)):
            sum += values[1][i] * values[0][i]
        avg = sum / self.d
        sum = 0
        for i in range(len(E)):

            sum += values[1][i] * (values[0][i] - avg)**2
        var = sum/self.d



        return np.sqrt(var), avg


    def draw_step_fit(self, delta, k, depth):
        
        E = self.Hamiltonian.Ew

        steps = np.arange(min(E), max(E), 0.01)

        step = step_func(steps, delta, k)

        poly = np.polyfit(steps, step, deg = depth)

        P = np.polyval(poly,steps)
        
        one_rescaling = max(P)
        if one_rescaling <1: 
            P = P* 1/one_rescaling

        eff_delta = np.sqrt(np.var(P))
        eff_delta = abs(steps[np.argmin(abs(eff_delta-P))]) * 2


        return step_func(steps, delta, k), P, steps, eff_delta
    

    def draw_gauss_fit(self, delta, k, depth):

        E_max = max(self.Hamiltonian.Ew)
        E_min = min(self.Hamiltonian.Ew)


        a = (E_max-E_min)/2
        b = (E_max+E_min)/2

        sigma = (delta)/(2*0.674490)

        E = self.Hamiltonian.Ew

        steps = np.arange(min(E), max(E), 0.01)
        step = np.exp(- (1/2)* (steps/sigma)**2)
        poly = np.polyfit(steps, step, deg = depth)

        P = np.polyval(poly,steps)
        eff_delta = delta
        return np.exp(- (1/2)* (steps/sigma)**2), P, steps, eff_delta
    
    def draw_gauss(self, delta, k):

        E_max = max(self.Hamiltonian.Ew)
        E_min = min(self.Hamiltonian.Ew)


        a = (E_max-E_min)/2
        b = (E_max+E_min)/2

        sigma = (delta)/(2* 0.674490 )

        E = self.Hamiltonian.Ew

        steps = np.arange(min(E), max(E), 0.01)

        return step_func(steps, delta, k), np.exp(- (1/2)* (steps/sigma)**2), steps 


    def plot_pol_filter(self, sigma, power = 2):
        
        H_in_energy_basis = self.Hamiltonian.in_energy_basis(self.Hamiltonian.H)
        dia = np.diag(H_in_energy_basis)

        test = []
        for i in range(len(dia)):
            test.append(self.Hamiltonian.rev_scale_hamilton(dia[i]))

        c2 = [['#00aedb','#00aedb'], ['#8ec127','#8ec127'], ['#d41243','#d41243']]
        c1 = ['#00aedb', '#8ec127', '#d41243', '#f47835']
        counter = 0
        for g in sigma:
            result1 = 0
            for n in range(250):
                result1 += e_pol(-1/2*(dia/g)**2, n)

            result = np.exp(-(1/2)*(dia/g)**2)
            eff_delta = 2*0.674490 * self.Hamiltonian.rev_scale_hamilton(g)


            avg = (max(self.Hamiltonian.Ew)+min(self.Hamiltonian.Ew))/2
            plt.plot(self.Hamiltonian.Ew, result1, c = c1[counter], label = r'$\Delta_{eff}$=' + str(round(eff_delta, 2)))
            plt.vlines([avg + eff_delta/2,avg -eff_delta/2], [0, 0], [1, 1], colors=c2[counter], linestyles='dashed')
            counter+= 1
        plt.xlim([-20, +20])
        plt.ylim([-0.1, 1.1])
        plt.xlabel(r'$E_n$')
        plt.ylabel(r'$\mathcal{P}(E_n)$')
        plt.legend()
        plt.show()

        return
    
    def plot_par_filter(self, alpha, power = 2):
        
        H_in_energy_basis = self.Hamiltonian.in_energy_basis(self.Hamiltonian.H)
        dia = np.diag(H_in_energy_basis)

       

        c2 = [['#00aedb','#00aedb'], ['#8ec127','#8ec127'], ['#d41243','#d41243']]
        c1 = ['#00aedb', '#8ec127', '#d41243', '#f47835', '#a200ff']
        counter = 0
        for a in alpha:
            result = 1 - a * dia**power
            plt.plot(self.Hamiltonian.Ew, result, c = c1[counter], label = r'$\alpha$=' + str(a))
            counter+= 1
        plt.xlabel(r'$E_n$')
        plt.ylabel(r'$\mathcal{P}(E_n)$')
        plt.legend()
        plt.show()

        return
    
    def plot_par_issues(self, alpha, power = 2):
        
        H_in_energy_basis = self.Hamiltonian.in_energy_basis(self.Hamiltonian.H)
        dia = np.diag(H_in_energy_basis)

        power = [2,  4, 6 , 8, 10 ]
        a = 1
        c2 = [['#00aedb','#00aedb'], ['#8ec127','#8ec127'], ['#d41243','#d41243']]
        c1 = ['#00aedb', '#8ec127', '#d41243', '#f47835', '#a200ff']
        counter = 0
        for p in power:
            result = 1 - a * dia**p
            plt.plot(self.Hamiltonian.Ew, result, c = c1[counter], label = str(p))
            counter+= 1
        plt.xlabel(r'$E_n$')
        plt.ylabel(r'$\mathcal{P}(E_n)$')
        plt.legend()
        plt.show()

        return
    
    def gamma_plt(self, w_s, w_d, E0, E_d):

        E= self.Hamiltonian.Ew 
        O = self.Observable.Matrix
        O_mn_avg = self.omega_avg(O, E0, E_d, w_s, w_d)
        
        if(self.omega_avg(abs(O), E0, E_d, w_s, w_d))**2 == 0:
            Gamma = 0
        else:
            Gamma = self.omega_avg(abs(O)**2, E0, E_d, w_s, w_d)/(self.omega_avg(abs(O), E0, E_d, w_s, w_d))**2

        return Gamma, O_mn_avg
    
    def omega_avg(self, O, E_a, E_d, w_s, w_d):
        E =self.Hamiltonian.Ew
        N = 0
        sum = 0
        for n in range(len(O)):
            for m in range(len(O)):
                if (abs(E[n] - E[m]) > w_s - w_d/2 and abs(E[n] - E[m]) < w_s + w_d/2 ) and  ( (E[n] + E[m]) > E_a-E_d and (E[n] + E[m]) < E_a+E_d ):
                    
                    sum += O[n][m] 
                    N += 1
        if N == 0 :
            return 0
        
        return sum/N
    
    def off_diagonal_hisogram(self, w, w_d, bins = 20):
        O = self.Observable.Matrix
        E = self.Hamiltonian.Ew
        
        

        for i in w:
            values = []
            for n in range(len(O)):
                for m in range(len(O[n])):
                    if (abs(E[n] - E[m]) > i - w_d/2 and  abs(E[n] - E[m]) < i + w_d/2 ):
                        values.append(O[n][m])
            
            plt.hist(values, bins = bins, histtype='step', label = r"$\omega=$ " + str(round(i,2)))
            plt.xlabel(r'$O_{nm}$')
            plt.ylabel(r"P$(O_{nm})$")
            plt.legend()
        
        plt.show()

        return

def e_pol(x, depth):
    sum = 0
    for n in range(depth):
        sum += 1/math.factorial(n) * (x)**n
    return sum

def step_func(E, d, k, E_0 = 0):

    return 1/(1 + np.exp(-2*k * (-E+(E_0 + d/2))) + np.exp(-2*k * (E-(E_0 - d/2))) + np.exp(-2*k*d))

def step_func_approx(E, d, k, depth = 10, E_0 = 0):

    return 1/(1 + e_pol(-2*k * (-E+(E_0 + d/2)), depth) + e_pol(-2*k * (E-(E_0 - d/2)), depth) + e_pol(-2*k*d, depth))



