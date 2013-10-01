#   This file is part of ssa.py.
#
#    ssa.py is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    ssa.py is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with ssa.py.  If not, see <http://www.gnu.org/licenses/>.
#-------------------------------------------------------------------------------
# Name:         ssa_root.py
# Purpose:      Perform fundamental Singular Spectrum Analysis methods in Python
# Context:      These functions were originally developed in Matlab by
#                 Eric Breitenberger circa 1995. With the exeption of functions
#                 that were not ideal for Python, the credits for each
#                 translated function are provided in each function's
#                 doc-string.
# Pythonista:   KSedmera
# License:      GPLv3 (see COPYING txt file)
# Version:      2013/09
#-------------------------------------------------------------------------------
import os, sys, time, math
try:
    from numpy import *
    import scipy
    import scipy.linalg as linalg
    import scipy.signal as sig
    from scipy.stats import chi2
    from scipy.optimize import fmin_powell
    import pylab as pl
except ImportError as exc:
    sys.stderr.write("Error: {}. Closing in 5 sec...\n".format(exc))
    print "Note: These tools require Python 2.6.5 - 2.7.5,"
    print "      AND several free science-related Python libraries. See which one your missing above."
    time.sleep(5);  sys.exit()

def eofsym(E, tol=10**4*spacing(1)):
    """
    EOFSYM - check symmetry of EOFs.
    Syntax: s=eofsym(E); s=eofsym(E,tol); s=eofsym(E,'1 or 0');

    Given a matrix of EOFs E, eofsym returns a vector
    s containing the symmetry of the EOF.

    If the i-th EOF is symmetric,  s(i)=1,
    if anti-symmetric,             s(i)=0.
    if neither sym. or anti-sym.,  s(i)=-1.

    s(i)=-1 is only possible if a non-Toeplitz (BK type)
    covariance matrix was used, or if the tolerance 'tol'
    is not set high enough. 'tol' is set to tol=10^4*eps
    by default, or it can be specified as the second argument.

    If 'tol' is set to the string '1 or 0' or ('0 or 1') the
    output will be forced to give only ones and zeros. EOFSYM
    will decide whether the EOFs are symmetric or anti-symmetric
    based on which assumption gives the lowest rms error.
    """

    (M,K)=shape(E)
    s=zeros(K)
    sign = lambda x: math.copysign(1, x)

    if not isinstance(tol,str):
        for k in range(K):    # Pick off left and right halves of EOF:
            if not M%2:
                L=E[:M//2-1][:,k]; R=E[-1:(M//2):-1][:,k]
            else:
                L=E[:(M-1)//2-1,k]; R=E[-1:((M+1)//2):-1][:,k]

            if max(abs(L-R)) < tol: s[k] = 1
            elif max(abs(L+R)) < tol:   s[k] = 0
            else:   s[k] = -1
            #disp(['Warning: EOF ' num2str(k) ' is neither symmetric or antisymmetric.'])

    elif tol == '1 or 0' or tol == '0 or 1':
        Esym = sum(sqrt((E[::-1,:]-E)**2))
        Easym = sum(sqrt((-E[::-1,:]-E)**2))
        s=(sign(Easym-Esym)+1)/2.
        ties=nonzero(s==1/2.)   # Resolve any ties by calling these symmetric
        s[ties]=ones((1,len(ties)))
    else:   raise ValueError('EOFSYM: Improper specification of tolerance.')

    return s

def eoffreq(E, n=500):
    """
    EOFFREQ - find the dominant frequencies of EOFs.
    Syntax: [mrft,f]=eoffreq(E); [mrft,f]=eoffreq(E,n);

    Given the eigenvector matrix E, eoffreq computes a
    normalized reduced Fourier transform (RFT) for each
    EOF. The maximum values of the RFTs are returned in
    'mrft', and the frequencies at which the maxima
    occur are returned in 'f'.

    The second (optional) argument 'n' gives the number
    of frequencies which will be checked - the default
    is 500. The frequencies are determined by splitting
    the Nyquist interval into n sections.
    """

    s=eofsym(E) # eigenvector symmetries

    (M,K)=shape(E)

    # center the eigenfunctions:
    E=E-ones((M,1))*mean(E)

    f=vstack(linspace(0,.5,n))
    F=zeros((n,K))
    j=r_[1:M+1]
    j2=j-(M+1)/2.

    Cc=f*j2
    Cs=sin(2.*pi*Cc)
    Cc=cos(2.*pi*Cc)
    r2c=sum((Cc**2).conj().T, axis=0)/float(M)
    r2s=1.-r2c
    r2s[0]=spacing(1) # to avoid divide by zero.
    if r2s[n-1] == 0:   r2s[n-1]=spacing(1)
    for k in range(K):
        if s[k] == 1: F[:,k]=abs(dot(Cc,E[:,k]))**2/r2c.conj().T
        elif s[k] == 0:   F[:,k]=abs(dot(Cs,E[:,k]))**2/r2s.conj().T
        elif s[k] == -1:  F[:,k]=abs(dot(Cc,E[:,k]))**2+abs(dot(Cs,E[:,k]))**2
        else:   raise ValueError('Elements of s must be 1, 0, or -1.')

    mrft = F.max(0)
    #print mrft.shape
    f = [F[:][:,i].tolist().index(j) for i,j in enumerate(mrft.tolist())]
    mrft=mrft/float(M)
    f=hstack([0.5*(i-1)/(n-1) for i in f])
    return [mrft,f]

def ac(x, lag=20, unbiased=True):
    """ac calculates the auto-covariance for series x out to k lags.
    The result has k+1 elements. The first element is the covariance at lag
    zero; succeeding elements 2:k+1 are the covariances at lags 1 to k."""
    N = size(x)
    if not isinstance(N,int):
        if N[0] < N[1]:
            x=transpose(x)
            N=size(x)
        else:   raise ValueError('Hey! Vectors only!')
    if lag >= N:  raise ValueError('Hey! Too big a lag!')
    #x = x - mean(x)
    c = zeros(lag+1)
    for i in xrange(lag+1):
        #print i, N
        c[i] = sum(x[:(N-i)]*x[i:])
        if unbiased:    c[i] = c[i]/(N-i)
        else:   c[i] = c[i]/N
    return c

def ssaeig(x, M):
    """Syntax: [E,V,C]=ssaeig(x, M);
    This function starts an SSA of series 'x', for embedding dimension 'M'.
    Returns:    E - eigenfunction matrix in standard form
                     (columns are the eigenvectors, or T-EOFs)
                V - vector containing variances (unnormalized eigenvalues)
                C - Covariance Matrix
                E and V are ordered from large to small.
    See section 2 of Vautard, Yiou, and Ghil, Physica D 58, 95-126, 1992."""
    from scipy.linalg import toeplitz, eigh
    N=size(x)
    if not isinstance(N, int):
        if N[0] < N[1]: x = transpose(x); N = size(x)
        else:   raise ValueError('Hey! Vectors only!')
    if M-1 >= N:  raise ValueError('Hey! Too big a lag!')

    acov=ac(x, M-1)            # calculate autocovariance estimates
    Tc = toeplitz(acov)        # create Toeplitz matrix (trajectory matrix)
    C = Tc
    L,E = eigh(Tc)          # calculate eigenvectors, values of T
    V = abs(L)              # create eigenvalue vector
    ind = argsort(V)        # sort eigenvalues
    ind = ind[M::-1]
    V = V[ind]
    E = E[:][:,ind]         # sort eigenvectors
    return [E,V,C]

def pc(x, E):
    #        Syntax: [A]=pc(x, E);
    #  This function calculates the principal components of the series x
    #  from the eigenfunction matrix E, which is output from ssaeig.m
    #  Returns:      A - principal components matrix (N-M+1 x M)
    #  See section 2.4 of Vautard, Yiou, and Ghil, Physica D 58, 95-126, 1992.
    #
    #  Matlab version written by Eric Breitenberger.    Version date 5/20/95
    #  Please send comments and suggestions to eric@gi.alaska.edu

    N = size(x)
    if not isinstance(N,int):
        if N[0] < N[1]:
            x = hstack(x)
            N = size(x)
        else:   raise ValueError('Hey! Vectors only!')
    Y = x-mean(x)
    M, cE = shape(E)
    if M != cE:    raise ValueError('E is improperly dimensioned')

    A = zeros((N-M+1,M))
    for i in range(N-M+1):
        w = Y[i:i+M]
        A[i][:] = dot(w, E)
    return A

def rc(A, E):
    # Syntax: R = rc(A,E)
    # This function calculates the 'reconstructed components' using the
    # eigenvectors (E, from ssaeig.m) and principal components (A, from pc.m).
    # R is N x M, where M is the embedding dimension used in ssaeig.m.
    #
    # See section 2.5 of Vautard, Yiou, and Ghil, Physica D 58, 95-126, 1992.
    #
    # Matlab version written by Eric Breitenberger.   Version date 5/18/95
    # Please send comments and suggestions to eric@gi.alaska.edu
    from scipy.signal import lfilter
    (M, Ec) = shape(E)
    (Ar, Ac) = shape(A)
    if M != Ec: raise ValueError('E is improperly dimensioned.')
    if Ac != M: raise ValueError('A is improperly dimensioned.')
    N = Ar+M-1
    R = zeros((N,M))
    Z = zeros((M-1,M))
    A = concatenate((A, Z), axis=0)
    # Calculate RCs
    for k in range(M):
        R[:,k] = lfilter(E[:,k], M, A[:,k], axis=0)
    # Adjust first M-1 rows and last M-1 rows
    for i in range(M-1):
        R[i][:] = R[i,:]*(M/float(i+1))
        R[(N-i-1)][:] = R[(N-i-1),:]*(M/float(i+1))
    return R

def itc(V,n):
    # Syntax: [kaic,kmdl,aic,mdl]=itc(V,n);
    # Compute signal/noise separation using information-theoretic
    # criteria. Two estimates are returned: the Akaike
    # information-theoretic criterion (AIC), and the minimum
    # description length (MDL). The order for which AIC or MDL is
    # minimum gives the number of significant components in the
    # signal. The two methods often give considerably different
    # results: AIC usually performs better than MDL in low SNR
    # situations, but MDL is a consistent estimator, and thus
    # performs better in large-sample situations.
    #
    # See Wax and Kailath, 1985, Trans. IEEE, ASSP-33, 387-392.
    #
    # Input:   V: an eigenspectrum (sorted in decreasing order)
    #          n: the number of samples used to compute V.
    # Outputs: kaic: the order for which AIC is minimum;
    #          kmdl: the order for which MDL is minimum;
    #          aic: vector containing AIC estimates;
    #          mdl: vector containing MDL estimates.
    #
    # Written by Eric Breitenberger, version 10/4/95, please send
    # any comments and suggestions to eric@gi.alaska.edu

    p=len(V)
    V=V[p-1:1:-1]

    # Calculate log-likelihood function:
    L=zeros((1,p))
    nrm=arange(p,0,-1) #p:-1:1
    sumlog=cumsum(log(V))
    sumlog=sumlog[p-1:1:-1]
    logsum=log(cumsum[V-1])
    logsum=logsum[p-1:1:-1]

    L=n*nrm**((sumlog/nrm)-logsum+log(nrm))

    # Calculate penalty function:
    pen=arange(0,p)
    pen=pen**(2*p-pen)

    # Calculate AIC and MDL, and find minima:
    aic=-L+pen
    mdl=-L+pen*log(n)/2
    kaic=nonzero(aic==min(aic))-1
    kmdl=nonzero(mdl==min(mdl))-1
    return [kaic,kmdl,aic,mdl]

def ssaconf(V, N):
    # Calculates various heuristic confidence limits for SSA.
    # Syntax: [f,g,v]=ssaconf(V,N);
    # Given a singular value vector V and the number of points
    # in the original series N, ssaconf returns vectors containing
    # the 95% confidence interval. These are calculated according
    # to the variance formulas of:
    #     f: Fraedrich 1986
    #     g: Ghil and Mo 1991a
    #     v: Vautard, Yiou, and Ghil 1992.
    # All estimates are for the 95% confidence level.
    # These simple estimates may be adequate for some purposes,
    # but none of them adequately consider the autocorrelation
    # of the time series. The estimates g and v are very similar
    # for N>>M. They are usually fairly conservative, as they
    # correspond to a decorrelation time of M. The Fraedrich
    # estimate is valid only for uncorrelated data, so it tends
    # to give error estimates which are too small.
    #
    # Written by Eric Breitenberger, version date 11/3/95.
    # Please send comments to eric@gi.alaska.edu

    M=len(V)
    f=2.*sqrt(2./N)*V
    g=sqrt(2.*M/(N-M))*V
    v=1.96*sqrt(M/(2.*N))*V
    return [f, g, v]

if __name__ == '__main__':
    """ The following lines were used to test each function. """
    x = loadtxt('test1.txt', delimiter=',')
    mx = mean(x)
    Y = x - mx
    [E, V, C] = ssaeig(array(Y),240)
    A = pc(x, E)
    R = rc(A, E)
    [mE, fE]=eoffreq(E)
    print 'len(x): %i\n'%len(x), 'len(V): %i; V=\n'%len(V), V[:20], V[-20:]
    print 'shape(E): (%i,%i); Ecs=\n'%shape(E), E[-1][:20], E[:20][:,-1]
    print 'shape(C): (%i,%i); C=\n'%shape(C), C[0][:20], C[:20][:,0]
    print 'shape(A): (%i,%i); A=\n'%shape(A), A[0][:20], A[:20][:,0]
    print 'shape(R): (%i,%i); R=\n'%shape(R), R[0][:20], R[:20][:,0]
    rx = [i+mx for i in sum(R, axis=1)]
    rn = [rx[i]-x[i] for i in range(len(x))]
    print 'sum(rn) = {}'.format(sum(rn)), max(rn), min(rn)
