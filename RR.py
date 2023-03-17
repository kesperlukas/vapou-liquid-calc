import numpy as np
import csv

def f(phi, K, z):
    # RR function to determine phi = V/F, depending on feed composition z_i and K-factors, with N components
    return sum(z[:-1]*(1-K[:-1])/(1+phi*(K[:-1]-1)))

def der_f(phi, K, z):
    # return the derivate of the function
    return sum(z[:-1]*(1-K[:-1])**2/(1+phi*(K[:-1]-1))**2)

def y_i(phi,K,z,i):
    return z[i]*K[i]/(1+phi*(K[i]-1))

def x_i(phi,K,z,i):
    return z[i]/(1+phi*(K[i]-1))

def check(K, z, T):
    # Initial check, since K_ammonia is below 1 and the rest is above:
    print('----------------------------------------------------------')
    print('Temperature ', T, ' °Celsius')
    print('Inital check: phi = ', f(1, K, z))

def display(res):
    # V_mole, V_i_mole, y_mole, V_mass, V_i_mass, y_mass, L_mole, L_i_mole, x_mole, L_mass, L_i_mass, x_mass
    print('Vapour mole flow:\n', np.asmatrix(res[0]))
    print('Vapour single mole flows:\n', np.asmatrix(res[1]))
    print('Vapour molar composition:\n', np.asmatrix(res[2]))
    print('Vapour mass flow:\n', np.asmatrix(res[3]))
    print('Vapour single mass flows:\n', np.asmatrix(res[4]))
    print('Vapour mass composition:\n', np.asmatrix(res[5]))
    print('Liquid mole flow:\n', np.asmatrix(res[6]))
    print('Liquid single mole flows:\n', np.asmatrix(res[7]))
    print('Liquid molar composition:\n', np.asmatrix(res[8]))
    print('Liquid mass flow:\n', np.asmatrix(res[9]))
    print('Liquid single mass flows:\n', np.asmatrix(res[10]))
    print('Liquid mass composition:\n', np.asmatrix(res[11]))

def save(res, T):
    data = [['Vapour mole flow:', res[0]],
    ['Vapour single mole flows:', res[1]],
    ['Vapour molar composition:', res[2]],
    ['Vapour mass flow:', res[3]],
    ['Vapour single mass flows:', res[4]],
    ['Vapour mass composition:', res[5]],
    ['Liquid mole flow:', res[6]],
    ['Liquid single mole flows:', res[7]],
    ['Liquid molar composition:', res[8]],
    ['Liquid mass flow:', res[9]],
    ['Liquid single mass flows:', res[10]],
    ['Liquid mass composition:', res[11]]]

    filename = "data_"+str(T)+".csv"
    # Open the file in write mode and create a csv writer object
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        # Write the data to the CSV file row by row
        for row in data:
            writer.writerow(row)

def solve(K, z, F_mole, M):
    # Find phi using the newton algorithm:
    phi = 0.5
    phi_new = 1
    threshold = 0.5
    steps = 0
    while threshold > 0.00001:
        phi_new = phi - f(phi, K, z)/der_f(phi, K, z)
        threshold = abs ((phi_new-phi)/phi)
        phi = phi_new
        steps += 1
    print('Final value for phi: ', phi_new, ' after ', steps, ' steps')

    # Solve equations and output all values:
    # Calculate streams in vapour phase
    V_mole = phi_new * F_mole
    # Use the RR definitions of y (gas composition) to define the stream:
    y_mole = []
    for i in [0,1,2,3]:
        y_mole.append(y_i(phi_new, K, z, i))
    y_mole.append(1-sum(y_mole))
    # Calculate single streams
    V_i_mole = []
    for i in [0,1,2,3,4]:
        V_i_mole.append(V_mole * y_mole[i])
    V_i_mass = []
    for i in [0,1,2,3,4]:
        V_i_mass.append(M[i]*V_i_mole[i]*1000)
    V_mass = sum(V_i_mass)
    y_mass = []
    for i in [0,1,2,3]:
        y_mass.append(V_i_mass[i]/V_mass)
    y_mass.append(1 - sum(y_mass))
    # Calculate streams in liquid phase
    L_mole = F_mole - V_mole
    # Use the RR definitions of x (liquid composition) to define the stream:
    x_mole = []
    for i in [0,1,2,3]:
        x_mole.append(x_i(phi_new, K, z, i))
    x_mole.append(1-sum(x_mole))
    # Calculate single streams
    L_i_mole = []
    for i in [0,1,2,3,4]:
        L_i_mole.append(L_mole * x_mole[i])
    L_i_mass = []
    for i in [0,1,2,3,4]:
        L_i_mass.append(M[i]*L_i_mole[i]*1000)
    L_mass = sum(L_i_mass)
    x_mass = []
    for i in [0,1,2,3]:
        x_mass.append(L_i_mass[i]/L_mass)
    x_mass.append(1 - sum(x_mass))
    return V_mole, V_i_mole, y_mole, V_mass, V_i_mass, y_mass, L_mole, L_i_mole, x_mole, L_mass, L_i_mass, x_mass,

if __name__ == '__main__':
    # Feed conditions:
    F = 1759.9 # kmol/h
    # Ammonia, Hydrogen, Nitrogen, Argon, Oxygen
    M = np.array([0.017031, 0.00100784, 0.0280134, 0.039948, 0.015999])
    z = np.array([0.160408, 0.630451, 0.209096, 0.000019559, 0.0000245245])
    # Initial check, since K_ammonia is below 1 and the rest is above:

    # For 60 degrees:
    T = 60
    K = np.array([ 0.300048, 54.7147, 38.3366, 5.84422, 9.71476])
    check(K, z, T)
    res = solve(K, z, F, M)
    display(res)
    save(res, T)

    # For 50 degrees:
    T = 50
    K = np.array([ 0.246429, 70.7681, 50.9996, 6.34951, 11.081])
    check(K, z, T)
    res = solve(K, z, F, M)
    display(res)
    save(res, T)

    # For 40 degrees:
    T = 40
    K = np.array([ 0.199187, 88.2431, 65.3795, 6.7159, 12.2772])
    check(K, z, T)
    res = solve(K, z, F, M)
    display(res)
    save(res, T)

    # For 30 degrees:
    T = 30
    K = np.array([ 0.157345, 108.005, 82.2143, 6.9679, 13.3482])
    check(K, z, T)
    res = solve(K, z, F, M)
    display(res)
    save(res, T)

    # For 20 degrees:
    T = 20
    K = np.array([ 0.116572, 130.943, 101.222, 6.99178, 14.108])
    check(K, z, T)
    res = solve(K, z, F, M)
    display(res)
    save(res, T)

    # For 10 degrees:
    T = 10
    K = np.array([ 0.0846036, 158.605, 125.077, 6.97037, 14.8498])
    check(K, z, T)
    res = solve(K, z, F, M)
    display(res)
    save(res, T)

    # For 0 degrees:
    T = 0
    K = np.array([ 0.0600113, 192.388, 155.517, 6.90837, 15.584])
    check(K, z, T)
    res = solve(K, z, F, M)
    display(res)
    save(res, T)

    # For -10 degrees:
    T = -10
    K = np.array([0.0414827, 234.197, 195.039, 6.8087, 16.3192])
    check(K, z, T)
    res = solve(K, z, F, M)
    display(res)
    save(res, T)

    # For -20 degrees:
    T = -20
    K = np.array([0.0278461, 286.669, 247.312, 6.67291, 17.0626])
    check(K, z, T)
    res = solve(K, z, F, M)
    display(res)
    save(res, T)

    # For -30 degrees:
    T = -30
    K = np.array([0.0180764, 353.519, 317.85, 6.50167, 17.8204])
    check(K, z, T)
    res = solve(K, z, F, M)
    display(res)
    save(res, T)

    # For -40 degrees:
    T = -40
    K = np.array([0.0112921, 440.094, 415.163, 6.29513, 18.5988])
    check(K, z, T)
    res = solve(K, z, F, M)
    display(res)
    save(res, T)

    # For -50 degrees:
    T = -50
    K = np.array([0.00674916, 554.253, 552.781, 6.05333, 19.4041])
    check(K, z, T)
    res = solve(K, z, F, M)
    display(res)
    save(res, T)

    # For -60 degrees:
    T = -60
    K = np.array([0.00383349, 707.853, 752.927, 5.77653, 20.2443])
    check(K, z, T)
    res = solve(K, z, F, M)
    display(res)
    save(res, T)

    # For- 70 degrees:
    T = -70
    K = np.array([ 0.00205289, 919.268, 1053.5, 5.4656, 21.1303])
    check(K, z, T)
    res = solve(K, z, F, M)
    display(res)
    save(res, T)

    # For -80 degrees:
    T = -80
    K = np.array([ 0.00102684, 1217.88, 1521.96, 5.1223, 22.0774])
    check(K, z, T)
    res = solve(K, z, F, M)
    display(res)
    save(res, T)

    # For -90 degrees:
    T = -90
    K = np.array([ 0.00047446, 1652.34, 2284.35, 4.74963, 23.109])
    check(K, z, T)
    res = solve(K, z, F, M)
    display(res)
    save(res, T)

    # For -100 degrees:
    T = -100
    K = np.array([ 0.000199872, 2306.52, 3590.15, 4.35201, 24.2609])
    check(K, z, T)
    res = solve(K, z, F, M)
    display(res)
    save(res, T)