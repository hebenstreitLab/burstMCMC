import scipy
import numpy
import math
from scipy import special
from scipy.special import hyp2f1
from scipy.stats import nbinom

def calculateHypergeom(a, b, tau, M):
    if all([scipy.special.hyp2f1(-n,-a,1-a-n,(1+b)/(math.exp(tau)+b)) > 0 for n in range(int(M + 1))]):
        f_n = [math.exp((n * (math.log(b) - math.log(1 + b)) + a * (math.log(1 + b * math.exp(-tau)) - math.log(1 + b))) + math.lgamma(a + n) - math.lgamma(a) - math.lgamma(n + 1) + math.log(scipy.special.hyp2f1(-n,-a,1-a-n,(1+b)/(math.exp(tau)+b)))) for n in range(int(M + 1))]
    else:
        f_n = [0]
    return f_n
