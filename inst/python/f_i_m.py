import scipy
import numpy
import math
from scipy import special
from scipy.special import hyp2f1
from scipy.stats import nbinom



def get_f_i_m(M, I, f_n, f_o, P_i):
    M, I = int(M), int(I)
    F_i, f_n, f_o, P_i = numpy.zeros((M, I + 1)), numpy.array(f_n), numpy.array(f_o), numpy.array(P_i)
    for m in range(1, M + 1):
        n = numpy.arange(0, m + 1)
        thetaN = n / m
        thetaO = 1 - thetaN
        fn, fo = f_n[:(m+1)], numpy.flip(f_o[:(m+1)])
        fnxo = fn * fo
        SigmaNxO = numpy.sum(fnxo)
        fnxo = fnxo / SigmaNxO
        f_i = [numpy.sum(fnxo * (thetaN * P_i[0, i] + thetaO * P_i[1, i])) for i in range(I + 1)]
        F_i[m - 1] = f_i
    return numpy.transpose(F_i)
