import matplotlib.pyplot as plt
import numpy as np

def plotsf(t, fn, xb, Q2, E, label, color):
    xs = fn(-t, xb, Q2, E)
    plt.plot(t, xs, color, label=label)

def plotw(c, fn, w, Q2, E, label, color):
    xs = fn(c, w, Q2, E)
    plt.plot(c, xs, color, label=label)

def plot_structure_functions():
    fig = plt.figure(figsize=(12, 10), dpi=100)
    plt.title("Structure Function vs -t", fontsize='16')
    plt.xlabel("-t", fontsize='13')
    plt.ylabel("xs, nb", fontsize='13')
    plt.legend(['sigma_T'], loc='upper center')
    plt.grid()

    Ebeam = 10.6
    Ebeam = 24.0
    q2 = 1.14
    xb = 0.131
    t = np.linspace(0., 2., 1000)
    plt.plot([0., 2.], [0., 0.], 'black')

    plotsf(t, xsigma_T, xb, q2, Ebeam, "$\\sigma_T$", 'red')
    plotsf(t, xsigma_TT, xb, q2, Ebeam, "$\\sigma_{TT}$", 'blue')
    plotsf(t, xsigma_LT, xb, q2, Ebeam, "$\\sigma_{LT}$", 'green')

    plt.legend(loc='upper right', prop={'size': 24})
    plt.show(block=False)
    plt.savefig('ds_dxB.png')

    fig = plt.figure(figsize=(12, 10), dpi=100)
    plt.title("Structure Function vs -t", fontsize='16')
    plt.xlabel("-t", fontsize='13')
    plt.ylabel("xs, nb", fontsize='13')
    plt.legend(['sigma_T'], loc='upper center')
    plt.grid()

    cost = 0.9
    w = 2.2
    cst = np.linspace(-1., 1., 1000)
    plt.xlabel("CosT", fontsize='13')
    plt.plot([0., 2.], [0., 0.], 'black')

    plotw(cst, sigma_T, w, q2, Ebeam, "$\\sigma_T$", 'red')
    plotw(cst, sigma_TT, w, q2, Ebeam, "$\\sigma_{TT}$", 'blue')
    plotw(cst, sigma_LT, w, q2, Ebeam, "$\\sigma_{LT}$", 'green')

    plt.legend(loc='upper center', prop={'size': 24})
    plt.show(block=False)
    plt.savefig('ds_dW.png')
