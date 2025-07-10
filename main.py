import Pauli_matrix as pm
import IsingFramework as JW
import numpy as np
import time
import matplotlib.pyplot as plt
import gc
import copy as copy
import argparse





# --------------------- Scripts --------------------- #

def Script_classical_method(power = 2):
    
    # Basic parameters
    p = np.arange(6, 12, 2)
    offset = 11
    L = np.arange(offset, 12, 1)
    q = 0 # Site of the Observable Pauli Matrix
    alpha = np.arange(0, 1.05, 0.05)
    

    
    
    lambdaT2_p6 = []
    lambdaT2_p8 = []
    lambdaT2_p10 = []

    c = 0


    for l in L:
        lambdaT2_p6.append([])
        lambdaT2_p8.append([])
        lambdaT2_p10.append([])



        sys = JW.System(l, 'C', q = 0, p = 0, energy_basis = False)
        



        for a in alpha:
            lambda_p6, dummy = sys.lambda_calc_alt(a, power = p[0])
            lambda_p8, dummy = sys.lambda_calc_alt(a, power = p[1])
            lambda_p10, dummy = sys.lambda_calc_alt(a, power = p[2])

            lambdaT2_p6[c].append(lambda_p6)
            lambdaT2_p8[c].append(lambda_p8)
            lambdaT2_p10[c].append(lambda_p10)
            
        

        del sys
        gc.collect()

        print('L=' + str(c + offset) + ' finished')
        c += 1


    for l in range(len(L)):


        plt.plot(alpha, lambdaT2_p6[l], color = '#00aedb', linewidth = 0.8, marker = 'v', label = str(p[0]),)
        plt.plot(alpha, lambdaT2_p8[l], color = '#8ec127', linewidth = 0.8, marker = 'o', label = str(p[1]))
        plt.plot(alpha, lambdaT2_p10[l], color = '#d41243', linewidth = 0.8, marker = 's', label = str(p[2]))




    # Reference of lambda = 1/2
    ref = []
    for a in alpha:
        ref.append(1/2)

    plt.plot(alpha, ref,color = '#a200ff', linewidth = 0.8)
    plt.xlabel(r'$\\alpha$')
    plt.ylabel(r'$\Lambda^T$')
    plt.legend()
    plt.show() 

    return




def Script_classical_method_gaussian(L_t, E0, q, point_1, point_2): # L = 11 : a = 16.44156513786688
    
    # Basic parameters
    p = np.arange(6, 12, 2)
    offset = L_t
    L = np.arange(offset, offset + 2, 1)
    sigma = np.arange(0.02, 2.2, 0.02)
    delta_lambda = np.arange(0.2,30.2,0.2)
    
    lambdas = []
    lambdaT2_p10 = []

    eff_delta_p10 = []

    c = 0


    for l in L:
        lambdas.append([])

        lambdaT2_p10.append([])
        eff_delta_p10.append([])



        sys_not_in_e = JW.System(l, 'C', q = 0, p = 0, energy_basis = True)
        sys_in_e = JW.System(l, 'C', q = 0, p = 0, energy_basis = True)
        


        # Alternative Method
        for d in delta_lambda:
            lambda_p10, eff_delta = sys_not_in_e.lambda_calc_alt_gaussian(d)
            
            lambdaT2_p10[c].append(lambda_p10)
            eff_delta_p10[c].append(eff_delta)
            window_1, gauss_1, steps_1 = sys_not_in_e.draw_gauss(point_1, 50)
            window_2, gauss_2, steps_2 = sys_not_in_e.draw_gauss(point_2, 50)


        # Exact result
        for d in delta_lambda:
            lambdas[c].append(sys_in_e.lambda_calc(E0,  d))
        
            
        

        del sys_not_in_e
        del sys_in_e
        gc.collect()

        print('L=' + str(c + offset) + ' finished')
        c += 1
    


    colors_exact = ['#00aedb' , '#a200ff']
    colors_gauss = ['#d41243' , '#f47835']
    for l in range(len(L)):
        plt.plot(delta_lambda, lambdas[l], color = colors_exact[l], linewidth = 0.8, marker = 'v', label = "L = " + str(offset +l) + ' | Hard-cut', ms = 4, markerfacecolor = 'None')
        plt.plot(eff_delta_p10[l], lambdaT2_p10[l],color = colors_gauss[l],  linewidth = 0.8, marker = 's', label = "L = " + str(offset + l) + ' | Gaussian filter', ms =4 , markerfacecolor = 'None')




    # Reference of lambda = 1/2
    ref = []
    for d in delta_lambda:
        ref.append(1/2)

    plt.plot(delta_lambda, ref,color = '#8ec127', label = r"$\Lambda^T = 1/2$", linewidth = 0.8)
    plt.xlabel(r'$\Delta_{(eff)}$')
    plt.ylabel(r'$\Lambda^T$')
    plt.legend()
    plt.show() 


    fig = plt.figure(figsize = (11, 7))

    ax0 = fig.add_subplot(221)
    ax0.plot(steps_1, window_1, color = '#00aedb', label = 'Hard-cut window')
    ax0.plot(steps_1, gauss_1, color = '#a200ff', label = 'Gaussian filter')
    ax0.vlines([point_1/2, -point_1/2], [-0.5, -0.5], [1.5, 1.5], colors=['#a200ff', '#a200ff'], label =r" $ \Delta_{eff} = $" + str(point_1) ,linestyles='dashed')
    ax0.set(xlabel = r'$E_n$', ylabel = r"$\mathcal{P}(E_n)$", xlim= [min(steps_1), max(steps_1)], ylim = [-0.5, 1.5])
    ax0.legend()

    ax1 = fig.add_subplot(222)
    ax1.plot(steps_2, window_2, color = '#00aedb', label = 'Hard-cut window')
    ax1.plot(steps_2, gauss_2, color = '#8ec127', label = 'Gaussian filter')
    ax1.vlines([point_2/2, -point_2/2], [-0.5, -0.5], [1.5, 1.5], colors=['#8ec127', '#8ec127'], label =r" $ \Delta_{eff} = $" + str(point_2) ,linestyles='dashed')
    ax1.set(xlabel = r'$E_n$', ylabel = r"$\mathcal{P}(E_n)$", xlim= [min(steps_2), max(steps_2)], ylim = [-0.5, 1.5])
    ax1.legend()

    ax2 = fig.add_subplot(212)
    ax2.plot(delta_lambda, lambdas[0], color = '#00aedb', linewidth = 0.8, marker = 'v',  label = 'Exact') # markerfacecolor='none',
    ax2.plot(eff_delta_p10[0], lambdaT2_p10[0],color = '#d41243',  linewidth = 0.8, marker = 's', label = 'Gaussian filter')
    ax2.plot(delta_lambda, ref,color = '#a200ff', linewidth = 0.8)
    ax2.set(xlabel = r'$\Delta_{(eff)}$', ylabel = r"$\Lambda^T$",  xlim= [-1, +30], ylim = [0.4, 1.2])
    ax2.vlines([point_1, point_2], [0.4, 0.4], [1.2, 1.2], colors=['#a200ff', '#8ec127'], linestyles='dashed') 
    ax2.legend()
    plt.show()

    return







def Script_classical_method_fit(L_t, f_shape = "gauss", q=0, k=50, depth=50, delta_1=1, delta_2=2): # L = 11 : a = 16.44156513786688
    
    # Basic parameters
    offset = L_t
    L = np.arange(offset, offset + 1, 1)
    q = 0 # Site of the Observable Pauli Matrix
    sigma = np.arange(0.2, 30.2, 0.2)
    delta_lambda = np.arange(0.01, 30.2, 0.1)
    k = 50
    depth = 50 
    delta_1 = 1
    delta_2 = 3.3
    

    
    lambdas = []
    lambdaT2_p10 = []
    eff_delta_p10 = []


    c = 0
    E0 = 0


    for l in L:
        lambdas.append([])

        lambdaT2_p10.append([])
        eff_delta_p10.append([])



        sys_not_in_e = JW.System(l, 'C', q = 0, p = 0, energy_basis = True)
        sys_in_e = JW.System(l, 'C', q = 0, p = 0, energy_basis = True)
        


        # Alternative Method
        if f_shape == "gauss":
            for g in sigma:
                lambda_p10, eff_delta = sys_not_in_e.lambda_calc_fit_gauss(g, depth)
                
                lambdaT2_p10[c].append(lambda_p10)
                eff_delta_p10[c].append(eff_delta)

            
            exact_1, poly_1, Ew_1, eff_d_1  = sys_not_in_e.draw_gauss_fit(delta_1, k, depth)
            exact_2, poly_2, Ew_2, eff_d_2 = sys_not_in_e.draw_gauss_fit(delta_2, k, depth)


        elif f_shape == "window":
            for d in delta_lambda:
                lambda_p10, eff_delta = sys_not_in_e.lambda_calc_step_fit(0, d, k = k, depth = depth)
                
                lambdaT2_p10[c].append(lambda_p10)
                eff_delta_p10[c].append(eff_delta)
            
            exact_1, poly_1, Ew_1, eff_d_1 = sys_not_in_e.draw_step_fit(delta_1, k, depth)
            exact_2, poly_2, Ew_2, eff_d_2 = sys_not_in_e.draw_step_fit(delta_2, k, depth)


        # Exact result
        for d in delta_lambda:
            lambdas[c].append(sys_in_e.lambda_calc(E0,  d))
        
            
        

        del sys_not_in_e
        del sys_in_e
        gc.collect()

        print('L=' + str(c + offset) + ' finished')
        c += 1


    colors_exact = ['#00aedb' , '#a200ff']
    colors_gauss = ['#d41243' , '#d41243']
    for l in range(len(L)):
        plt.plot(delta_lambda, lambdas[l], color = colors_exact[l], linewidth = 0.8, marker = 'v', label = 'Hard-cut energy window')
        plt.plot(eff_delta_p10[l], lambdaT2_p10[l], color = colors_gauss[l],  linewidth = 0.8, marker = 's', label = 'Gaussian filter')




    # Reference of lambda = 1/2
    ref = []
    for d in delta_lambda:
        ref.append(1/2)

    plt.plot(delta_lambda, ref,color = '#a200ff', linewidth = 0.8)
    plt.xlim([-1, +30])
    plt.ylim([0.4, 1.2])
    plt.xlabel(r'$\Delta_{(eff)}$')
    plt.ylabel(r'$\Lambda^T$')
    plt.legend()
    plt.show() 

    fig = plt.figure(figsize = (13, 9))

    eff_d_1 = round(eff_d_1, 2)
    eff_d_2 = round(eff_d_2, 2)

    ax0 = fig.add_subplot(221)
    ax0.plot(Ew_1, exact_1, color = '#00aedb', label = 'Exact ' + str(f_shape) + r' filter| $\Delta = $' + str(delta_1))
    ax0.plot(Ew_1, poly_1, color = '#a200ff', label = 'Polynomial fit')
    ax0.vlines([eff_d_1/2, -eff_d_1/2], [-0.5, -0.5], [1.5, 1.5], colors=['#a200ff', '#a200ff'], label =r" $ \Delta_{eff} = $" + str(eff_d_1) ,linestyles='dashed')
    ax0.set(xlabel = '$E_n$', ylabel = r"$\mathcal{P}(E_n)$", xlim= [min(Ew_1), max(Ew_1)], ylim = [-0.5, 1.5])
    ax0.legend()

    ax1 = fig.add_subplot(222)
    ax1.plot(Ew_2, exact_2, color = '#00aedb', label = 'Exact ' + str(f_shape) + r' filter| $\Delta = $' + str(delta_2))
    ax1.plot(Ew_2, poly_2, color = '#8ec127', label = 'Polynomial fit')
    ax1.vlines([eff_d_2/2, -eff_d_2/2], [-0.5, -0.5], [1.5, 1.5], colors=['#8ec127', '#8ec127'], label =r" $ \Delta_{eff} = $" + str(eff_d_2) ,linestyles='dashed')
    ax1.set(xlabel = r'$E_n$', ylabel = r"$\mathcal{P}(E_n)$", xlim= [min(Ew_1), max(Ew_1)], ylim = [-0.5, 1.5])
    ax1.legend()

    ax2 = fig.add_subplot(212)
    ax2.plot(delta_lambda, lambdas[0], color = '#00aedb', linewidth = 0.8, markersize = 5, marker = 'v', label = 'Hard-cut energy window')
    ax2.plot(eff_delta_p10[0], lambdaT2_p10[0],color = '#d41243',  linewidth = 0.8, markersize = 5, marker = 's', label = f_shape + ' filter fit')
    ax2.plot(delta_lambda, ref,color = '#a200ff', linewidth = 0.8)
    ax2.set(xlabel = r'$\Delta_{(eff)}$', ylabel = r"$\Lambda^T$" ,  xlim= [-1, +30], ylim = [0.4, 1.2])
    ax2.vlines([eff_d_1, eff_d_2], [0.4, 0.4], [1.2, 1.2], colors=['#a200ff', '#8ec127'], linestyles='dashed') 
    ax2.legend()

    plt.show() 


    return


def Script_with_eff_delta(power = 2):
    
    # Basic parameters
    p = np.arange(6, 12, 2)
    offset = 8
    L = np.arange(offset, 9, 1)
    q = 0 # Site of the Observable Pauli Matrix
    alpha = np.arange(0, 1, 0.008)
    delta_lambda = np.arange(0.2,25.2,0.2)
    
    
    lambdas = []
    lambdaT2_p10 = []

    eff_delta_p10 = []

    c = 0
    E0 = 0


    for l in L:
        lambdas.append([])

        lambdaT2_p10.append([])
        eff_delta_p10.append([])



        sys_not_in_e = JW.System(l, 'C', q = 0, p = 0, energy_basis = False)
        sys_in_e = JW.System(l, 'C', q = 0, p = 0, energy_basis = True)
        


        # Alternative Method
        for a in alpha:
            lambda_p10, eff_delta = sys_not_in_e.lambda_calc_alt(a, power = 2)

            lambdaT2_p10[c].append(lambda_p10)
            eff_delta_p10[c].append(eff_delta)


        # Exact result
        for d in delta_lambda:
            lambdas[c].append(sys_in_e.lambda_calc(E0,  d))
        
            
        

        del sys_not_in_e
        del sys_in_e
        gc.collect()

        print('L=' + str(c + offset) + ' finished')
        c += 1


    for l in range(len(L)):
        plt.plot(delta_lambda, lambdas[l], color = '#00aedb', linewidth = 0.8, marker = 'v', label = str(p[0]))
        plt.plot(eff_delta_p10[l], lambdaT2_p10[l], color = '#d41243', linewidth = 0.8, marker = 's', label = str(p[2]))




    # Reference of lambda = 1/2
    ref = []
    for a in alpha:
        ref.append(1/2)

    plt.plot(alpha, ref,color = '#a200ff', linewidth = 0.8)
    plt.xlabel(r'$\\alpha$')
    plt.ylabel(r'$\Lambda^T$')
    plt.legend()
    plt.show() 


    return



def Script(L_t, q, power = 2):
    """
        This Function Creates instances of classical Method and Pauli Matrix System Objects with the respected size.
        With that the Relation of Moment Indicator is calculated for the Pauli Matrix Method and the classical Method.
        A Plot for the Two Methods is then printed
    """
        
    # Basic parameters
    offset = L_t
    L = np.arange(offset, offset + 3, 1)
    alpha = np.arange(0, 1.1, 0.1)


    JW_execute = True
    
    JW_time = 0
    sec_time = 0
    fou_time = 0
    

    st = time.time()



    lambdaT1 = []
    lambdaT2 = []

    c = 0


    for l in L:
        lambdaT1.append([])
        lambdaT2.append([])


        if JW_execute:

            stJw = time.time()

            sys = JW.System(l, 'C', q = 0, p = 0, energy_basis = False)
            
            etJw = time.time()
            JW_time += etJw-stJw

        else:
            sys = None


        System1 = pm.System(l, q, sys)


        for a in alpha:
            

            st_sec = time.time()
            
            # Second Moment
            POP = System1.POP(a, power= power)
            second_moment = System1.Moment(POP)

            et_sec = time.time()
            sec_time += et_sec - st_sec


            st_fou = time.time()

            # Fourth Moment
            POP_squared = System1.POP_squared(a, power = power)
            fourth_moment = System1.Moment(POP_squared)

            et_fou = time.time()
            fou_time += et_fou - st_fou


            lambdaT1[c].append(second_moment**2 / fourth_moment)
            

            if JW_execute:

                stJw = time.time()
                #print(sys.Hamiltonian.Ew)
                lambda_dummy, dummy = sys.lambda_calc_alt(a, power = power)
                lambdaT2[c].append(lambda_dummy)

                
                etJw = time.time()
                JW_time += etJw-stJw
        

        del System1
        del POP
        del POP_squared
        del sys
        gc.collect()

        print('L=' + str(c + offset) + ' finished')
        c += 1
        


    et = time.time()
    elebaration_time = et-st

    print('runtime overall: ' + str(elebaration_time) + ' s')
    print('count_time : ' + str(pm.count_time) + ' s')
    print('sort_time : ' + str(pm.sort_time) + ' s')
    print('JW_time : ' + str(JW_time) + ' s')
    print('sec_time : ' + str(sec_time) + ' s')
    print('fou_time : ' + str(fou_time) + ' s')
    print('var_time : ' + str(pm.var_time) + ' s')


    # Reference of lambda = 1/2
    ref = []
    for a in alpha:
        ref.append(1/2)



    # Plot
    fig, axs = plt.subplots(2, sharex=True, sharey=True)


    axs[0].plot(alpha, lambdaT1[0], linewidth = 0.8, marker = 'v', c = '#00aedb', label = 'L=' + str(0 + offset))
    axs[1].plot(alpha, lambdaT2[0], linewidth = 0.8, marker = 'v', c = '#00aedb')

    axs[0].plot(alpha, lambdaT1[1], linewidth = 0.8, marker = 'o', c = '#d41243', label = 'L=' + str(1 + offset))
    axs[1].plot(alpha, lambdaT2[1], linewidth = 0.8, marker = 'o', c = '#d41243')

    axs[0].plot(alpha, lambdaT1[2], linewidth = 0.8, marker = 's', c = '#a200ff', label = 'L=' + str(2 + offset))
    axs[1].plot(alpha, lambdaT2[2], linewidth = 0.8, marker = 's', c = '#a200ff')


    axs[0].plot(alpha, ref, linewidth = 0.8, c = '#8ec127')
    axs[1].plot(alpha, ref, linewidth = 0.8, c = '#8ec127')

    axs[1].set(xlabel = '\u03B1')


    for figure in axs:
        figure.set(ylabel = r'$\Lambda^T$')

    fig.legend(loc="upper right")
    plt.show() 


    return


def classical_method(L_0, E0, delta):
    """
        Creates a classical Method System Object and calculates Quantitys with it 
    """

    t = np.arange(0, 50, 0.5)
    delta = 2
    delta_lambda = np.arange(0.2,30.2,0.2)
    lambdas = []
    L = [L_0, L_0+1, L_0 + 2]

    for l in L:
        sys = JW.System(l, 'C', q = 0, p = 0, energy_basis = True)
        sys.auto_print(E0, delta, t)
        lambdas.append(sys.lambda_print(E0, delta_lambda))

    ref_1 = []
    for a in delta_lambda:
        ref_1.append(1/2)
    c = ['#00aedb', '#d41243', '#a200ff' ]
    m = ['v', 'x', 's']
    plt.figure()
    for l in range(len(L)):
        plt.plot(delta_lambda, lambdas[l], linewidth = 0.8, label = "L= " + str(L_0 + l), c = c[l], marker = m[l])
    plt.plot(delta_lambda, ref_1, linewidth = 0.8, label = '1/2', c = '#8ec127')
    plt.legend(loc="upper left")
    plt.xlabel('\u0394')
    plt.ylabel(r'$\Lambda^T$')
    plt.show()
    return



def classical_method_energy_gap_histogram(L):
    """
        Calculates an Histogram of the Energy Gaps of the Hamiltonian of the System
    """

    sys = JW.System(L, 'C', q = 0, p = 0, energy_basis = False, symmetry_subspace=True, symmetry_value= int(L/2))

    
    Eigenvalues = sys.Hamiltonian.Ew
    
    gaps = []
    for i in range(len(Eigenvalues) - 1):
        gaps.append( abs(Eigenvalues[i] - Eigenvalues[i+1]))

    
    plt.hist(gaps, bins=80, range = (0, 0.1))
    plt.show()
    return



def Hamiltonian_powers(loops, L):
    """
        Calculates Powers of the Hamiltonian in the auli Matrix Method
    """
    System1 = pm.System(L, 0, sys = None)

    H = System1.Hamiltonian()
    power = H
    exp = []
    a = np.arange(1, loops, 1)

    i = 1
    diff =  1
    recent = 0
    while i<loops :

        print(i)

        power = pm.Product([power, H])
        power = power.execute()
        power.Hash_count(L)

        exp.append(System1.Moment(power))
        diff = abs(recent -len(power.terms))

        recent = len(power.terms)
        print("Power: " + str(i) + " Coeffizienten: ")
        print("differenz: " + str(diff))
        
        i += 1
    return


def eth_classic_diag(L_t, E0, delta):

    
    delta_sq = []
    L_arr = []
    Ed = delta

    fig, ax1 = plt.subplots()
    


    L = L_t
    sys = JW.System(L, 'C', q = 0, p = 0, energy_basis = True)
    diag = np.diag(sys.Observable.Matrix)
    x = np.multiply(1/L,sys.Hamiltonian.Ew)
    ax1.scatter(x, diag,c = '#00aedb', marker = 'v', label = 'L=' + str(L))

    O = sys.Observable.Matrix
    E = sys.Hamiltonian.Ew

    sum1 = 0
    N1 = 0
    for n in range(len(E)):
        if E[n] >= E0-Ed and E[n] <= E0+Ed:
            sum1 += O[n][n]**2
            N1 +=1
    sum1 = sum1/N1

    sum2 = 0
    N2 = 0
    for n in range(len(E)):
        if E[n] >= E0-Ed and E[n] <= E0+Ed:
            sum2 += O[n][n]
            N2 +=1
    sum2 = (sum2/N2)**2

    delta_sq.append(sum1 - sum2)
    L_arr.append(L)



    L = L_t + 1
    sys = JW.System(L, 'C', q = 0, p = 0, energy_basis = True)
    diag = np.diag(sys.Observable.Matrix)
    x = np.multiply(1/L,sys.Hamiltonian.Ew)
    ax1.scatter(x, diag,c = '#8ec127', marker = 'o', label = 'L=' + str(L))

    O = sys.Observable.Matrix
    E = sys.Hamiltonian.Ew

    sum1 = 0
    N1 = 0
    for n in range(len(E)):
        if E[n] >= E0-Ed and E[n] <= E0+Ed:
            sum1 += O[n][n]**2
            N1 +=1
    sum1 = sum1/N1

    sum2 = 0
    N2 = 0
    for n in range(len(E)):
        if E[n] >= E0-Ed and E[n] <= E0+Ed:
            sum2 += O[n][n]
            N2 +=1
    sum2 = (sum2/N2)**2

    delta_sq.append(sum1 - sum2)
    L_arr.append(L)




    L = L_t + 2
    sys = JW.System(L, 'C', q = 0, p = 0, energy_basis = True)
    diag = np.diag(sys.Observable.Matrix)
    x = np.multiply(1/L,sys.Hamiltonian.Ew)
    ax1.scatter(x, diag,c = '#d41243', marker = 's', label = 'L=' + str(L))

    O = sys.Observable.Matrix
    E = sys.Hamiltonian.Ew

    sum1 = 0
    N1 = 0
    for n in range(len(E)):
        if E[n] >= E0-Ed and E[n] <= E0+Ed:
            sum1 += O[n][n]**2
            N1 +=1
    sum1 = sum1/N1

    sum2 = 0
    N2 = 0
    for n in range(len(E)):
        if E[n] >= E0-Ed and E[n] <= E0+Ed:
            sum2 += O[n][n]
            N2 +=1
    sum2 = (sum2/N2)**2

    delta_sq.append(sum1 - sum2)
    L_arr.append(L)

    

    ax1.set_xlabel(r'$E_m/L$')
    ax1.set_ylabel(r'$E_{mm}$')
    ax1.legend(loc = 'lower right')

    left, bottom, width, height = [0.3, 0.7, 0.25, 0.17]
    ax2 = fig.add_axes([left, bottom, width, height])
    ax2.plot( L_arr, delta_sq, marker = 's', c ='#a200ff')
    ax2.set_xlabel('L')
    ax2.set_ylabel(r'$\sigma_d^2$')
    ax2.set_yscale(value = 'log')

    plt.show()
    return


def chaos_indic(L, E0, delta):

    sys = JW.System(L, 'C', q = 0, p = 0, energy_basis = False, symmetry_subspace=True, symmetry_value=np.real(get_symmetry_values(L)[int(5/2)]))

    E = sys.Hamiltonian.Ew


    s = []

    for n in range(len(E)- 1):
        s.append(abs(E[n + 1] - E[n]))

    rt = []
    n = 1
    while n < len(s):
        if not s[n-1] == 0:
            rt.append(s[n]/s[n-1])
        else:
            rt.append(0)
        n += 1

    r = []
    for i in range(len(rt)):
        if not(rt[i] == 0): 
            r.append(min([rt[i],1/rt[i]]))
        else:
            r.append(0)

    r_avg = np.average(r)

    print(r_avg)

    return


def Gamma(E0, delta, L_t):

    w_s = np.arange(0.2, 5, 0.2)
    w_d = 0.05
    goal = np.pi / 2
    L = np.arange(L_t, L_t + 3, 1)
    

    sys = JW.System(L[0], 'C', q = 0, p = 0, energy_basis = True)
    
    gamma_array= []
    ref = []
    for w in w_s:
        gamma, omn = sys.gamma_plt(w, w_d, E0, delta)
        gamma_array.append(gamma)
        ref.append(goal)
    
    plt.plot(w_s, gamma_array, c = '#00aedb', marker = 'v', label = 'L=' + str(L[0]))
    

    sys = JW.System(L[1], 'C', q = 0, p = 0, energy_basis = True)
    
    gamma_array= []
    ref = []
    for w in w_s:
        gamma, omn = sys.gamma_plt(w, w_d, E0, delta)
        gamma_array.append(gamma)
        ref.append(goal)
    
    plt.plot(w_s, gamma_array, c = '#8ec127', marker = 'o', label = 'L=' + str(L[1]))

    sys = JW.System(L[2], 'C', q = 0, p = 0, energy_basis = True)
    
    gamma_array= []
    ref = []
    for w in w_s:
        gamma, omn = sys.gamma_plt(w, w_d, E0, delta)
        gamma_array.append(gamma)
        ref.append(goal)
    
    plt.plot(w_s, gamma_array, c = '#d41243', marker = 's', label = 'L=' + str(L[2]))
    plt.plot(w_s, ref)

    plt.ylabel(r"$\Gamma(\\tilde{\omega})$")
    plt.xlabel(r"$\\tilde{\omega}$")
    plt.legend()
    plt.xlim([0, 5])
    plt.ylim([1, 2])
    
    plt.show()


    return







# --------------------- Help Functions --------------------- #
def get_symmetry_values(L):

    symmetry_values = []

    for l in range(L):
        symmetry_values.append(np.exp(-1j*2*np.pi/L * l))
    return symmetry_values


def print_pol_filter(L):

    alpha = np.arange(0.1, 1.23, 0.5)
    sys = JW.System(L, 'C', q = 0, p = 0, energy_basis = False)

    sys.plot_pol_filter(alpha, power = 10)
    return

def print_par_filter():

    alpha = np.arange(0, 1.1, 0.25)
    sys = JW.System(8, 'C', q = 0, p = 0, energy_basis = False)

    sys.plot_par_filter(alpha, power = 2)
    return

def print_par_issues():

    alpha = np.arange(0, 1.1, 0.25)
    sys = JW.System(8, 'C', q = 0, p = 0, energy_basis = False)

    sys.plot_par_issues(alpha, power = 2)
    return

def print_off_diag_hist():
    w_d = 0.05
    w = np.arange(0.2, 2, 0.6)

    sys = JW.System(11, 'C', q = 0, p = 0, energy_basis = True)

    sys.off_diagonal_hisogram(w, w_d, bins = 500)

    return


def polynomial_fit():
    x = np.linspace(-15, 15, 1000)
    sigma = 4
    gaus = np.exp(- (1/2)* (x/sigma)**2)

    poly = np.polyfit(x, gaus, deg = 10)

    fig, ax = plt.subplots()
    ax.plot(gaus, label = "gaus")
    ax.plot(np.polyval(poly, x), label = "fit")
    ax.legend()
    plt.show()
    return


def step_fit():
    x = np.linspace(-15, 15, 1000)
    sigma = 4
    step_f = JW.step_func(x)

    poly = np.polyfit(x, step_f, deg = 40)

    fig, ax = plt.subplots()
    ax.plot(step_f, label = "step_f")
    ax.plot(np.polyval(poly, x), label = "fit")
    ax.legend()
    plt.show()
    return















# --------------------- Testing --------------------- #

def Heap_sort_test(L):
    """
        Tests for the Heap count Method

        Calculates a Sum Object and then counts it with the Heap count method and the old method
        Then prints out both so they can be compared
    """

    System1 = pm.System(L, 0, None)

    Hsquared = System1.POP_squared(1)
    #pm.Product([System1.Hamiltonian(), System1.Hamiltonian()]).execute()
    
    print('before')
    #Hsquared.print()
    H_counter = copy.deepcopy(Hsquared)
    
    Hsquared.Hash_count(System1.L)
    H_counter.Heap_count(System1.L)
    print('counter')
    print(System1.Moment(Hsquared))
    #Hsquared.print()
    print('H_counter')
    print(System1.Moment(H_counter))
    #H_counter.print()

    return




def translation_test(L):
    sys = JW.System(L, 'C', q = 0, p = 0, energy_basis = False, symmetry_subspace=True, symmetry_value=np.real(get_symmetry_values(5)[int(5/2)]))

    return


def performance_testing(I):
    L = 10


    System1 = pm.System(L, 0, None)

    Hamiltonian = System1.Hamiltonian()

    for i in range(I):
        Hamiltonian = pm.Product([Hamiltonian, System1.Hamiltonian()]).execute()
    
    print('dimension of sum')
    print(len(Hamiltonian.terms))

    st = time.time()
    Hamiltonian.Hash_count(L, sort = True)
    et = time.time()
        
    print('runtime: ')
    print(et-st)
    print('dimension of sum')
    print(len(Hamiltonian.terms))
    print('count_time: ' + str(pm.count_time))
    print('sort_time: ' + str(pm.sort_time))
    # second_moment = System1.Moment(POP)


    return













# --------------------- Entrypoint --------------------- #
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Simulationswerkzeug für Spin-Ketten-Modelle",
        usage="%(prog)s -M METHOD [options]",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""
Verfügbare Methoden für -M / --Method:

MethodBasic                  Einfache Simulation mit Magnetfeld an Position q
ClassicalMethod              Klassische Approximation im Energiefenster um E0
ClassicalGaussian            Klassische Methode mit Gauß-Filter, gesteuert über Delta1/2
ClassicalFit                 Klassische Methode mit anpassbarer Filterform (f), Tiefe (D)
Gamma                        Simulation zur Gamma-Berechnung bei gegebener Energie
ChaosIndicator               Berechnung eines Chaos-Indikators für die Kette
ClassicalETHDiag             Klassische ETH-Diagonalisierung im Energiefenster
ClassicalEnergyGapHistogram  Histogramm der Energieabstände für gegebene Länge

Beispiel:
python main.py -M ClassicalGaussian -L 6 --E0 0 --Delta 1
    """
    )

    parser.add_argument("-M", "--Method", help="Methode.", required=True)
    parser.add_argument("-L", "--Length", help="Länge der Spin Kette.", type=int, required=False, default=3)

    parser.add_argument("-E", "--E0", help="Zentrum des Energiefensters", type=float, required=False, default=0)
    parser.add_argument("-q", "--MagnetPos", help="Spin Position an der das Magnetfeld wirkt.", type=int, required=False, default=0)
    parser.add_argument("-f", "--Function", help="Filter Funktion.", required=False, default="gauss")
    parser.add_argument("-k", "--K", help="Spin Position an der das Magnetfeld wirkt.", type=float, required=False, default=0)
    parser.add_argument("-D", "--Depth", help="Approximationstiefe.", type=int, required=False, default=50)
    parser.add_argument( "--Delta", help="Delta", type=float, required=False, default=0.5)
    parser.add_argument( "--Delta1", help="Inset Filter 1 (Effekte Filterbreite)", type=float, required=False, default=1)
    parser.add_argument( "--Delta2", help="Inset Filter 2 (Effekte Filterbreite)", type=float, required=False, default=2)
    

    try:
        args = parser.parse_args()
    except SystemExit as e:
        if e.code != 0:
            parser.print_help()
            exit(1)
        raise
    
    print(f"Ising Modell Kettenlänge: {args.Length}")
    print(f"Methode: {args.Method}")

    match args.Method :
        case "MethodBasic":
            Script(
                L_t=args.Length,
                q = args.MagnetPos
            )
        case "ClassicalMethod":
            classical_method(
                L_0=args.Length,
                E0=args.E0,
                delta=args.Delta
            )
        case "ClassicalGaussian":    
            Script_classical_method_gaussian(
                L_t=args.Length,
                q=args.MagnetPos,
                point_1=args.Delta1,
                point_2=args.Delta2,
                E0=args.E0
            )
        case "ClassicalFit":
            Script_classical_method_fit(
                L_t=args.Length, 
                f_shape=args.Function, 
                q=args.MagnetPos,
                k=args.K,
                depth=args.Depth,
                delta_1=args.Delta1,
                delta_2=args.Delta2
            )
        case "Gamma":
            Gamma(
                L_t=args.Length,
                E0=args.E0,
                delta=args.Delta,
            )
        case "ChaosIndicator":
            chaos_indic(
                L=args.Length,
                E0=args.E0,
                delta=args.Delta
            )
        case "ClassicalETHDiag":
            eth_classic_diag(
                L_t=args.Length,
                E0=args.E0,
                delta=args.Delta
            )
        case "ClassicalEnergyGapHistogram":
            classical_method_energy_gap_histogram(
                L=args.Length
            )
        case _:
            print("Unbekannte Methode")

    