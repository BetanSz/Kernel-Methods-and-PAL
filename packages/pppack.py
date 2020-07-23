# !/usr/bin/env python
# --*-- coding:utf-8 --*--
"""
This module contains the definitions of the objects used to
represent function spaces and functions, in both one- and
multi-dimensional form. A few methods to operate on tensors
as equivalent one-dimensional arrays are also included.

.. todo:: Knot optimization for B-spline approximations is
          not supported yet.

.. todo::
   Chebyshev multi-dimensional representations are not
   provided and should be implemented in the FORTRAN library.
   This development is postponed, as it hasn't been needed
   so far.

.. todo::
   Conversion from B-spline form to PP form in case of
   multi-dimensional representation is currently not
   supported.

.. todo:: Provide a license agreement to the package.
"""
__title__="PPPACK"
__author__ = "D. Tomatis"
__company__ = "CEA/DEN/DANS/DM2S/SERMA/LPEC"
__date__ = "25/04/2017"
__version__ = "v1.1.0"

import numpy as np
import itertools as itt
import sys
#import lib.pppack as ppk
import lib.pppack.pppack as ppk
#import lib.chebyshev_interp_1d.chebyshev_interp_1d as pch
from IPython import embed

def tau_spline(tau, gtau, ntau, gamma):
    """
    This is the las example fo the book for trying to deal with oscilations in splines

    """
    print 'start'

    tau_aux=np.zeros(len(tau), order='F')
    gtau_aux=np.zeros(len(tau), order='F')

    for i in range(len(tau)):
        tau_aux[i]=np.float64(tau[i])
        gtau_aux[i]=np.float64(gtau[i])

    ppk.tautsp(tau=tau_aux,gtau=gtau_aux,ntau=np.int32(ntau),gamma=np.float64(gamma))
    """
    error!
    ValueError: failed to create intent(cache|hide)|optional array-- must have defined dimensions but got (4,-1,)
    """
    print 'end'


class Fd(list):
    r"""
    Class representing the space :math:`Fd \in \mathbb{R}` of
    real-valued functions defined on a *d*-dimensional domain
    obtained by Cartisian products as
    :math:`V^{(d)} = V_{1}\times \ldots \times V_{d}`, with
    :math:`V_i \in \mathbb{R}\, \forall i`.

    *Fd* is also defined by the tensor product of
    one-dimensional function space :math:`F_i(V_i)`.
    """
    _els=0 # nb of elements in V_i used to generate Fd spaces
    def __init__(self, *args, **kwargs):
        """
        Init a multivariate function space as a list of
        univariate function spaces. Uniqueness of the
        components of a Fd instance is ensured.
        """
        #print "entrando a __init__ Fd con",args
        #print "len(args)=",len(args)
        if "name" in kwargs: self.name=kwargs["name"]
        for F in args:
            if isinstance(F,(S,PP)):
                #print "Fd:init:obj en args es S/PP"
                if F.name!='':
                    nm=F.name
                else: # use a standard name
                    nm="F%d"%(self._els); Fd._els+=1
                    F.name=nm
                self.append(F)
            elif isinstance(F,Fd):
                #print "Fd:init:obj en args es Fd"
                for Fi in F: self.append(Fi)
            else: raise TypeError("Unknown function space")
        #print "Fd:init:Chekeos finales"
        Fs=[Fi.name for Fi in self]
        if len(Fs)!=len(set(Fs)):
            raise AttributeError("Fd has non-unique F elements %s", Fs)

    @property
    def dims(self):
        """Return the dimension of the multivariate space."""
        return max(1,len(self))

    @property
    def a(self): return [(F.name,F.a) for F in self]

    @property
    def b(self): return [(F.name,F.b) for F in self]

    def __str__(self):
        l="F^(%d) = "%self.dims
        l+=" x ".join([ V.name for V in self ])
        l+="\n with,\n"
        l+='\n'.join([ ("%s = "%V.name+str(V)) for V in self ])
        return l

    def __repr__(self):
        l="Fd("+", ".join([ repr(F) for F in self ])+')'
        return l

    def __mul__(self,Q):
        #print "Fd: Doing mult"
        return Fd(self,Q)

class PP(Fd):
    r"""
    Class representing a Piecewise Polynomial function space of
    polynomial order :math:`k`, breakpoints :math:`\vec{\xi}` and
    continuities of derivativies :math:`\vec{\nu}`,
    :math:`PP_{k,\vec{\xi}(,\vec{\nu})}`. Simple polynomial space
    :math:`P_k \subseteq PP_{k,\xi}`, with equality in case of
    :math:`\xi = [a,b]`. In this case, the data points :math:`\tau`
    are expected in :math:`[a,b]`.

    :Example:

    >>> P0 = PP(4,[-1., 0., .5, 1.])
    >>> print P0
    PP_{ k=4, xi, nu }, with
     xi = [-1.000,0.000,0.500,1.000],  nu = [-,3,3,-]
    """
    def __init__(self,k,xi,nu=None,name=''):
        r"""
        Init :math:`PP_{k,\xi[,\nu]}`, with default :math:`\nu_i=k-1,
        \,\forall i` (inside breakpoints) if :math:`\nu=\mbox{None}`.
        """
        self.name=name
        if k<0: raise ValueError("invalid negative order k")
        self.k=k
        self.xi=sorted(xi, key=float)
        # check for unique breakpoints in xi
        for x in xi:
            if sum([c==x for c in xi])>1:
                raise ValueError("multiple breakpoint %e detected"%x)
        if nu!=None:
            if len(nu)!=self.l-1:
                 raise ValueError("invalid size of nu/xi")
            if any( vi<0 for vi in nu ):
                 raise ValueError("invalid negative element in nu")
            if any( vi>k for vi in nu ):
                 raise ValueError("invalid element > k in nu")
            self.nu=nu
        else:
            # k-1 null jumps of derivatives up to the (k-1)-th order
            # at breakpoints (default)
            self.nu = [(k-1) for i in range(2,self.l+1)]

    @property
    def l(self): return (len(self.xi)-1)

    @property
    def a(self): return self.xi[0]

    @property
    def b(self): return self.xi[-1]

    def PP2S(self):
        r"""
        Get :math:`S_{k,t}` by the theorem of Curry and Schoenberg (1966)
         :math:`S_{k,t} = PP_{k,\vec{\xi},\vec{\nu}} \mbox{ on }
         [t_k, t_{n+1}]`.
        """

        n=self.k*self.l - sum(self.nu)
        t=[self.xi[0] for i in range(k)]
        if len(self.nu)>0:
            for i in range(1,self.l): t.extend( \
                [self.xi[i] for j in range(k-self.nu[i-1])])
        t.extend([self.xi[-1] for i in range(k)])
        if n+k!=len(t): raise RuntimeError("PP2S: dims mismatch")
        return S(k,t)

    def cmpPPtau(self):
        r"""
        Compute the data points from :math:`P_{k,\vec{\xi},\vec{\nu}}`
        (osculatory interpolation as default).
        """
        self.basefncs="trpowf" # truncated power functions
        tau=[self.xi[0] for i in range(self.k)]
        for i in range(1,self.l):
            tau.extend([self.xi[i] for j in range(self.k-self.nu[i-1])])
        return tau

    def cmptau(self):
        """General interface, here alias to cmpPPtau.""" #!!!Why is this necessary or useful?
        return self.cmpPPtau()

    def __repr__(self):
        return ("PP(k=%d,xi="%self.k+str(self.xi)+ \
                ",nu="+str(self.nu)+")")

    def __str__(self):
        l="PP_{ k=%d, xi, nu }, with\n xi = ["%self.k
        l+=','.join(["%.3f"%c for c in self.xi])+"],  nu = [-,"
        return l+','.join([ "%d"%n for n in self.nu])+',-]'

class S(Fd):
    r"""
    Class representing the function space populated by Bsplines
    of polynomial order *k* and knot sequence :math:`\vec{t}`.

    :Example:

    >>> S0 = S(4,[1., 1., 2., 3., 4., 4.])
    >>> print S0
    S_{ k=4, t }, with
     t = [1.000,1.000,2.000,3.000,4.000,4.000]
    """
    def __init__(self,k,t,name=''):
        r"""Initialize :math:`S(k,\vec{t})`."""
        self.name=name
        if k<0: raise ValueError("invalid negative order k=%i"%k)
        self.k=k
        self.t=t
        if self.n<=0: raise ValueError("invalid knot sequence, n=%i"%self.n)

    @property
    def n(self): return (len(self.t)-self.k)

    @property
    def a(self): return self.t[0]

    @property
    def b(self): return self.t[-1]

    @property
    def tstar(self):
        r"""Compute the Greville control points :math:`t^*`."""
        tstar=np.zeros((self.n+1,), order='F')
        for i in range(self.n+1):
            tstar[i] =np.sum(self.t[i+1:i+self.k])
        tstar/=float(self.k-1)
        return tstar



    def tau4Ik(self,k_is_even=True):
        r"""
        Even order interpolation at knots, say :math:`k=2m`, reproducing
        :math:`I_k`, that is spline interpolation of order *k* with the
        **not-a-knot** end condition.
        """
        k,n,t=self.k,self.n,self.t
        if k%2!=0:
            k_is_even=False
            print("WARNING: k order is not even")
        tau,m=np.zeros((n,), order='F'),k/2
        tau[0],tau[-1]=self.a,self.b
        dt=(t[k]-t[k-1])/float(m)
        for i in range( m ): tau[i+1]=tau[i]+dt
        for i in range(n-k): tau[m+i]=t[k+i]
        j=n-k+m
        if not k_is_even:
            m+=1
            dt=(t[n]-t[n-1])/float(m)
        for i in range(m): tau[j+i]=tau[j+i-1]+dt
        return tau

    def cmptau(self):
        """General interface, here alias to tau4Ik."""
        return self.tau4Ik()

    def __repr__(self):
        return "S(k=%d,t="%self.k+str(self.t)+")"

    def __str__(self):
        l="S_{ k=%d, t }, with\n t = ["%self.k
        return l+','.join([ "%.3f"%v for v in self.t])+']'

class fd(Fd):
    r"""
    Class representing a multivariate function *f* as

    .. math:: f\in Fd: V^{(d)} \rightarrow \mathbb{R}

    with :math:`V^{(d)} = V_1 \times \ldots \times V_d`,
    :math:`V_i=\times[a_i,b_i] \in \mathbb{R},\, \forall\ 1\leq i\leq d`

    In virtue of the tensor product :math:`f(V_1,\ldots,V_d)= F_1(V_1)
    \bigotimes \ldots \bigotimes F_d(V_d)` and
    provided that one-dimensional basis functions are available
     .. math::
              \mathcal{B}(F_i)=\{\psi_{i,j_i}(x_i)\}:=
              \{\psi_{j_i}(x_i)\}\text{ for }
              1\leq i\leq d\text{ and } 1\leq j_i

    Then:
     .. math::
              f&=F_1(V_1) \bigotimes \ldots \bigotimes
              F_d(V_d)= \sum_{j_1} C_{j_1} \psi_{j_1}(x_1)
              \bigotimes\ldots \bigotimes
              \sum_{j_d} C_{j_d} \psi_{j_d}(x_d)
              =\sum_{\vec{j_1}}\ldots\sum_{\vec{j_d}} C_{j_1,\ldots,j_d}
              \prod_{i}\psi_{j_i}(x_i)
              =\sum_{\vec{j}} C_{\vec{j}} \psi_{\vec{j}}(\vec{x}),
              \mbox{ with } \vec{j}=[j_1, \ldots, j_d]

    The coefficients :math:`C_{\vec{j}}` are *d*-dimensional
    FORTRAN array, stored as equivalent one-dimensional array:

    .. math:: C_{\vec{j}} = C(j_1 + n_1(j_2 - 1 + n_2(j_3 - 1
                          + \ldots + n_{d-1}(j_d - 1) \ldots ))).

    The C version with zero-based numbering is:

    .. math:: C_{\vec{j}} = C(j_1 + n_1(j_2 + n_2(j_3
                          + \ldots + n_{d-1}(j_d) \ldots ))).

    The approximant function *f* to *g* is determined by applying
    the linear functional :math:`\lambda_{\vec{i}}` to *g*, and
    so to :math:`\{\psi_{j_i}\}`. The inversion of the associated
    directional Gramian matrices is provided by the specific
    class to which the functions :math:`\{\psi_{j_i}\}` belong.
    Currently, only :math:`\lambda_{i}g_i = [\tau_i]g(x_i,\dot)`
    is available in the module.

    :math:`C_{\vec{j}}` and :math:`g_{\vec{j}}` (gtau in the
    following) are equivalent one-dimensional FORTRAN arrays,
    whereas data points follow as lists :math:`(\tau_{j_i},
    1\leq j_i\leq I_i,\; \forall i)`.

    .. note:: (*tau*, *gtau*) is to be considered as private to
              ensure the correspondence with the fitting
              coefficients *C*.
    """
    coef=None  # C_{\vec{i}}
    _tau=None  # data points
    _gtau=None # value of the function to fit at the data points

    # tau and gtau (plus optional b.c.) are related biunivocally to coef
    # hence, they are treated here as private attributes, only getters are
    # provided explicitly. Setters come with the method cmpcoef.
    @property
    def tau(self): return self._tau

    @property
    def gtau(self): return self._gtau

    @tau.setter
    def tau(self, tau): self._tau=tau

    @gtau.setter
    def gtau(self, g): self._gtau=g

    @classmethod
    def inFd(cls,Fdin):
        """Init from the space *Fd* to which the function belongs."""
        if not isinstance(Fdin, Fd):
            raise AttributeError("non-Fd obj in input")
        return fd(Fdin)

    def __str__(self):
        l=self.__class__.__name__ \
         +" \in "+super(Fd,self).__str__()+'\n'
        if not self.coef is None:
            l+="coef = \n"+str(self.coef)
            l+="\n tau = "+str(self._tau)
            l+="\ngtau = "+str(self._gtau)
        else:
            l+="coef (to be determined)"
        return l+"\n---"

    def cmpcoef(self, tau, gtau):
        r"""
        Compute the tensor coefficients of the *d*-dimensional
        function fd, with input values :math:`gtau=fd(x_\vec{j})`
        as one-dimensional FORTRAN array and data points tau in a
        list as :math:`[ [x_{j_i=1}^{J_i}], \forall j_i ]`. The
        caller must ensure the correspondence between the *d*-tuples
        *x* with the gtau values according to the definition
        of the equivalent one-dimensional array.

        The output array of an inner product applies an index shift
        w.r.t. the input array, according to De Boor's intuition,
        thus always performing efficient operations on contiguous
        data in memory :cite:`deBoor1979TOMS`. Methods for vectorize
        and matricize are here needed because the C-wrapper f2py
        breaks the useful internal FORTRAN array management. This
        drawback will hopefully be removed in future versions, after
        f2py upgrading.
        """
        Ni=[len(t) for t in tau]
        N=np.prod(Ni)
        # interpolation requires as many data points as function values
        if N!=len(gtau): raise ValueError("Mismatch in data sites %i!=%i"%(N,len(gtau)))


        #print 'here'
        #print tau
        #print gtau
        #print sys.exit()

        coef=np.float64(gtau)
        for i in range(self.dims):
            F=self[i]
            ###?print fd.__class__, F.__class__ #why one makes sense and the other no? out: <type 'type'> <class 'packs.pppack.PP'>
            ###? the fact that this is mandatory since s and not S contains the way to find coefficients doesn't undeline that the distinction between S and s is somewhat faible?
            #print "fd:cmpcoef:mi objeto F es",F
            #print "fd:cmpcoef: f instance Fd=",isinstance(F,Fd)
            if not isinstance(F,fd):
                if   isinstance(F,S):  f=s.inS(F)
                elif isinstance(F,PP): f=pp.inPP(F)
                # other types are already prevented by Fd init
            #print "fd:cmpcoef:mi objeto f es",f
            # compute the coefficient along index i
            #print "fd:cmpcoef:tau[i]",tau[i]
            #print "fd:cmpcoef:coef=",coef
            #print "fd:cmpcoef:gtau",gtau
            #print f
            f.cmpcoef(tau=tau[i],gtau= coef )
            #print "fd:cmpcoef:f.coef.shape",f.coef.shape
            #print "fd:cmpcoef:f.coef=",f.coef
            ### comment old part not using FORTRAN ndarray referencing
            ###f.cmpcoef(tau=tau[i],gtau= matricize(coef,Ni) )

            # override the list object with the one-dimensional
            # function type f in F (inheriting F methods) to access
            # the methods of f later on at fd evaluation
            if not isinstance(F,fd): self[i]=f
            coef=np.float64(f.coef)
            #print "fd:cmpcoef:the variable coef takes the value of f.coef",coef
            ###coef=vectorize(f.coef)
            ###Ni=Ni[1:]+Ni[:1] # shift periodically of one position
        #print "fd:cmpcoef:final step the final value is coef=", coef
        # tau, gtau and coef must change at the same time
        self._tau,self._gtau,self.coef=tau,gtau,coef

    def __call__(self, x, j=0):
        r"""
        Evaluate :math:`D^j f^{(d)}` at
        :math:`\vec{x}=[x_1, \ldots, x_d]`.
        """

        try: lx=len(x)
        except TypeError: return "x=%s is not a list/tuple!"%str(x)
        if lx!=self.dims:
            raise ValueError("input x, of len(x)=%i is not a %d-tuple."%(lx,self.dims))
        if self.coef is None: raise AttributeError("missing coef.")
        ### comment to enable FORTRAN ndarray spacial treatment
        ###Ni=[len(t) for t in self._tau]
        ###N=np.prod(Ni)

        val=np.float64(self.coef)
        #val=self.coef
        #print val[0][0],isinstance(val[0][0],float) ###? is this cast of type necessary?
        #print "fd:call: x=",x
        #print "fd:call: (self.coef) val=",val
        for i in range(self.dims):
            f=self[i]
            #print "fd:call: tomo el espacio f=",f
            ###f.coef = matricize(val,Ni)
            #### For the first dimension (dims=0) after the first time and x is evaluated this line is redundant.
            # It could be added an if (self.dims>0 or "this is the first calculation of x")
            #print "fd:call: en f.coef=",f.coef
            #print "fd:call: meto val=",val
            f.coef = np.float64(val)


            val = f(x[i],j)
            #print "fd:call: calculo val=f(x[i],j)",val, "x[i]=",x[i]
            ###del Ni[0] # crop dimension in the val tensor
        #print "fd:call: from the final vector val",val, "i return val[0]", val[0]
        #print val
        #sys.exit()
        return val[0]

class pp(fd,PP):
    r"""
    Class representing a functional element *pp* in *PP* as:

    .. math:: f(x) = \sum_{j=1}^k C_{j,i} P_{j,i}(x),
              \mbox{ with }\xi_i \le x < \xi_{i+1}

    The function and its derivatives are represented (and evaluated)
    as (de Boor 1987, ch. 7 pp. 89):

    .. math:: D^j f(x) = \sum_{m=j}^{k-1} C_{m+1,i}(x-\xi_i)^{m-j}
                       / (m-j)!

    where either

    .. math::
              \begin{align*}
                  i = 1 \mbox{ and }&~~~~~~~~   x < \xi_2,\mbox{ or} \\
                  1<i<l \mbox{ and }&\xi_i \leq x < \xi_{i+1},\mbox{ or}\\
                  i = l \mbox{ and }&\xi_l \leq x.
              \end{align*}

    .. note:: (tau, gtau) is to be considered as private to ensure
              the correspondence with the fitting coefficients.
    """
    # f.coef becomes here COEF(K,L) as C_{i,j} above
    basefncs=None # type of base used in the P_k(,xi,nu) representation
    # inherit the init from PP

    @classmethod
    def inPP(cls,Pkxi):
        """init from the space *PP* the function belongs to."""
        if not isinstance(Pkxi, PP):
            raise AttributeError("non-PP obj in input")
        ### change Pkxi to Pkxinu ?
        return pp(Pkxi.k,Pkxi.xi,Pkxi.nu,Pkxi.name)

    def cmpcoef(self, tau=None, gtau=None, dgtau1=0., dgtauN=0., \
                ibcbeg=0, ibcend=0, method="newton"):
        r"""
        Compute the coefficients of the piecewise polynomial stored
        in coef. Available methods are:

         * CUBSPL for cubic splines + b.c. (see below),
         * NEWTON for simple polynomials in Newton form by divided
           differences (default by method input arg in case of
           :math:`P_k`) or CHEBYSHEV polynomial interpolation.

        Output are the fitting coefficients, -> coef(k,L=N-1).
        ibcbeg/ibcend:

         0. no b.c., so use the not-a-knot condition
            (:math:`jump D^3f(\xi_1)=0`)
         1. complete cubspl with :math:`f'(\tau_1)=\mbox{coef(2,1)}`
         2. :math:`f''(\tau_1)=\mbox{coef(2,1)}` for natural condition

    do j = 1, k
        (replace :math:`\tau_1` with :math:`\tau_N` above for the
        corresponding b.c. of ibcend in the input dgtau1/N).

        .. note:: Input derivatives are expected in dgtau1/N. They
                  will be inserted in coef (a.k.a. C in the fortran
                  subroutine CUBSPL).

        .. note:: The 'private' (tau, gtau) couple is stored after
                  the calculation of the fitting coefficients.
        """
        #print gtau
        ltau,lgtau=len(tau),np.size(gtau) #!!!! why not len?
        k,l=np.int32(self.k),self.l

        # embed()
        flag_nonewton=False
        if l > 1 and flag_nonewton:
            #print "pp:cmpcoef:cubspl"
            print 'here'
            print l
            print k
            # if k!=4: raise ValueError("4th order required!")
            n=np.int32(l+1)
            coef=np.zeros((k,n), order='F')
            if n!=lgtau: raise ValueError("gtau values mismatch")
            coef[0,:]=gtau # remind output coef[0,:]=gtau!
            if ibcbeg>0: coef[1, 0]=dgtau1
            if ibcend>0: coef[1,-1]=dgtauN
            ppk.cubspl( tau=np.float64(self.xi), \
                        c=np.float64(coef), \
                        ibcbeg=np.int32(ibcbeg), \
                        ibcend=np.int32(ibcend), n=n )
            self._tau=self.cmpPPtau(); self._gtau=gtau #!!!! what is that self._tau=self.cmpPPtau()
        else:
            #print "pp:cmpcoef:no cubspl"
            if ltau!=k: raise ValueError("k data points are needed (ltau=%i)!=(k=%i)"%(ltau,k))
            if lgtau%ltau!=0: raise ValueError("tau/gtau mismatch")
            # check sorting order of tau
            # if not all(tau[i] <= tau[i+1] for i in range(k-1)):
            #     raise ValueError("unsorted tau")
            m=np.int32(lgtau/ltau)
            if method=="newton":
                #print "pp:cmpcoef:method newton, m=",m
                if m==1:
                    coef=ppk.newton_coef_1d ( n=k, \
                        x=np.float64(tau), f=np.float64(gtau) )
                else:
                    coef=ppk.newton_coef_nd ( n=k, m=m, \
                        x=np.float64(tau), f=np.float64(gtau) )
            elif method=="chebyshev":
                if m==1:
                    coef,xmin,xmax=pch.chebyshev_coef_1d ( nd=k, \
                        xd=np.float64(tau), yd=np.float64(gtau) )
                else: raise ValueError("ND chebyshev not available yet")
            else: raise ValueError("unknown input argument: method")
            self._tau,self._gtau=tau,gtau
            self.basefncs=method
        # output coef seems not being as fortr<an contiguous, so...
        if not coef.flags['F_CONTIGUOUS']:
            coef=np.asfortranarray(coef, dtype=np.float64)
        self.coef=coef

    def pp2s(self, tau=None, gtau=None):
        r"""
        Converts from piecewise polynomial to B-spline form. Data points
        follow by default the even order interpolation at knots to build
        the interpolant :math:`I_k(g)`.
        """
        s_from_pp=s.inS( self.PP2S() )
        if  tau is None:
          s_from_pp._tau=self.tau4Ik()
          s_from_pp._gtau=self._gtau
        elif gtau is None: raise ValueError("missing input gtau")
        return s_from_pp

    def __call__(self, x, j=0):
        """
        Evaluate :math:`D^j pp` by PPVALU at :math:`x`.
        """

        if self.coef is None: raise AttributeError("missing coef.")
        ltau,lcoef=len(self._tau),np.size(self.coef)
        k,l,j4=np.int32(self.k),np.int32(self.l),np.int32(j)
        flag_nonewton=False
        if self.l>1 and flag_nonewton:
            # embed()
            f = ppk.ppvalu( breaks=np.float64(self.xi), l=l, k=k, \
                              coef=np.float64(self.coef[:,:-1]), \
                            x=np.float64(x), jderiv=j4 )
        else:
            m=np.int32(lcoef/ltau)
            if self.basefncs=="newton":
                tau,coef=np.float64(self._tau),np.float64(self.coef)
                if m==1:
                    ni=np.int32(np.size(x))
                    if isinstance(x,list):
                        if not x.flags['F_CONTIGUOUS']:
                            x=np.asfortranarray(x, dtype=np.float64)
                    elif ni==1: x=[x]
                    else: raise ValueError("Invalid input x="+str(x))
                    #print "value_1d for m=%i"%m
                    f = ppk.newton_value_1d ( n=k, ni=ni, \
                            x=tau, d=coef, xi=np.float64(x) )
                else:
                    #print "value_nd for m=%i"%m
                    f = ppk.newton_value_nd ( n=k, m=m, \
                            x=tau, d=coef, xi=np.float64(x) )
            elif self.basefncs=="chebyshev":
                if m==1:
                    xmin,xmax=np.float64(self.a),np.float64(self.b)
                    f = pch.chebyshev_value_1d ( nd=k, \
                            c=np.float64(self.coef), \
                            xmin=xmin, xmax=xmax, xi=np.float64(x) )
                else:
                    raise ValueError("ND chebyshev not available yet")
            else:
                raise RuntimeError("missing call method for this P_k")
        return f

class s(fd,S):
    r"""
    Class representing a functional element *s* in *S* as:

    .. math:: f(x) = \sum_{i=1}^n \alpha_i B_{i,k,t}(x)

    if :math:`t_j \leq x \leq t_{j+1}` for some
    :math:`j \in [k,n]`, then by the compact support of B-splines it is

    .. math:: f(x) = \sum_{i=j-k+1}^j \alpha_i B_{i,k,t}(x).

    .. note:: (tau, gtau) is to be considered as private to ensure
              the correspondence with the fitting coefficients.
    """
    # f.coef becomes here BCOEF(N) as \alpha_i above
    # inherit the init from S

    @classmethod
    def inS(cls,Skt):
        """Init from the space *S* to which the function belongs."""
        if not isinstance(Skt, S):
            raise AttributeError("non-S obj in input")
        return s(Skt.k,Skt.t,Skt.name)

    def cmpcoef(self, tau, gtau, VDS=False):
        r"""
        Compute the coefficients of the spline of order K with knots
        t(1:n+k), which takes on the value gtau(i) at tau(i), for
        i = 1 to n. The i-th equation of the linear system
        A * BCOEF = B
        for the B-spline coefficients of the interpolant enforces
        interpolation at tau(1:n). Hence, B(i) = gtau(i), for all i,
        and A is a band matrix with 2*k-1 bands, if it is invertible.

        If VDS==True, the Schoenberg's Variational Diminishing Spline
        approximation is used (BCOEF(i) = gtau( tau_i^star ) for all i).
        """
        ltau,lgtau=len(tau),np.size(gtau)
        n,k=np.int32(self.n),np.int32(self.k)

        if ltau!=self.n: raise ValueError("data points (len(tau)=%i) is not n=%i"%(ltau,n))
        if lgtau%ltau!=0: raise ValueError("tau/gtau data mismatch")
        if not VDS:
            # check non-singularity of the matrix B_i(tau_j)
            if self.singularBij(tau):
                raise RuntimeError("Singular B_i(tau_j) matrix detected"+ \
                                   " by the Schoenberg-Whitney's theorem")
            m=np.int32(lgtau/ltau)
            #print "s: Splint collocation of number of dimensions =",m
            if m==1:
                #print "s:performing splint"
                q,self.coef,iflag = ppk.splint( n=n, k=k, \
                                        tau=np.float64(tau), \
                                       gtau=np.float64(gtau), \
                                          t=np.float64(self.t) )
            else:
                #print "s:performing spli2d"
                work, \
                q,self.coef,iflag = ppk.spli2d( n=n, k=k, m=m, \
                                        tau=np.float64(tau), \
                                       gtau=np.float64(gtau), \
                                          t=np.float64(self.t) )
        else:
            #print("WARNING: remind BCOEF(i) <- gtau( tau_i^star ) for all i")
            self.coef=gtau
        if iflag!=1: raise RuntimeError("SPLINT/SPLI2D failed")
        # tau, gtau and coef must change at the same time
        self._tau,self._gtau=tau,gtau

    def singularBij(self,tau):
        r"""
        Check non-singularity of the matrix :math:`B_i(\tau_j), i,j=1,
        \ldots,n` by applying the theorem of Schoenberg and Whitney (1953).
        """
        # check strictly increasing order of tau
        if not all(tau[i] < tau[i+1] for i in range(self.n-1)):
            raise ValueError('the tau sequence is not strictly increasing')
        singularBij,n,k=False,self.n,self.k
        for i in range(n):
            # verify t_i = ... = t_{i+r} = tau_j with r < k
            if self.t[0] < tau[i] < self.t[-1]:
                tm=np.sum(self.t==tau[i])
                if tm>k-1:
                    print("%d internal knots (>k-1) = tau[%d] = %e"%(tm,i,tau[i]))
                    singularBij=True; break
            # verify t_i < tau_i < tau_{i+k}, for all i
            if not( self.t[i] < tau[i] < self.t[i+k] ):
                if i==0   and ( self.t[i]<=tau[i]< self.t[i+k] ): continue
                if i==n-1 and ( self.t[i]< tau[i]<=self.t[i+k] ): continue
                print("?!? i=%d, t[i]=%e < tau[i]=%e < t[i+%d]=%e ?!?"% \
                      (i,self.t[i],tau[i],k,self.t[i+self.k]))
                singularBij=True; break
        return singularBij

    def s2pp(self):
        """
        Converts from B-spline to piecewise polynomial form.
        """
        k,n=np.int32(self.k),np.int32(self.n)
        scrtch, breaks, coef, l = ppk.bsplpp( k=k, n=n, \
                      t=np.float64(self.t), bcoef=np.float64(self.coef) )
        if sum(breaks[l+1:])>0.: raise RuntimeError("BSPLPP failed.")
        pp_from_s=pp(self.k,breaks[:l+1])
        pp_from_s.coef=coef[:,:l]
        # set tau/gtau
        pp_from_s.tau=pp_from_s.cmpPPtau(); pp_from_s.gtau=coef[0,:l]
        return pp_from_s

    def __call__(self, x, j=0):
        """
        Evaluate :math:`D^j s` by BVALUE at :math:`x`.
        """
        if self.coef is None: raise AttributeError("missing coef.")
        ltau,lcoef=len(self._tau),np.size(self.coef)
        n,k,j4=np.int32(self.n),np.int32(self.k),np.int32(j)
        m=np.int32(lcoef/ltau)
        if m==1:
            #print "s: corriendo bvalue con"
            #print "s: n=",n,"k=",k,"x=",x,"coef=",self.coef
            f = ppk.bvalue( n=n, k=k, x=x, jderiv=j4, \
                      t=np.float64(self.t), bcoef=np.float64(self.coef) )
            # put value in numpy.ndarray to agree with fd.__call__
            f=np.array([f,], order='F')
        else:
            #print "s:call: running bvalnd because m=",m
            #print "s:call: using coef=",np.float64(self.coef)
            f = ppk.bvalnd( n=n, m=m, k=k, x=x, jderiv=j4, \
                      t=np.float64(self.t), bcoef=np.float64(self.coef) )
        return f



if __name__ == "__main__":
    import copy as cp
    k,t,xi=4,np.zeros((17,), order='F'),[0.,1.]
    n=len(t)-k
    t[0:k]=0.; t[n:]=10.
    for i in range(1,10): t[k+i-1]=float(i)



    S1=S(k,t)

    P1=Fd(PP(3,xi),name="quadratic")
    P2=   PP(2,xi); P2.name="linear"
    V1=P1*P2; V1.name="poly"
    S2=cp.deepcopy(S1)
    V2=Fd(S1,S2, name="2D-bspline")
    V3=Fd(P1,P2,PP(2,xi,name="linear2"), name="3D-poly")
    V4=Fd(P1,S2,P2, name="3D-mix")


    def cmpgtau(tau,g):
        Ni=[len(t) for t in tau]
        N=np.prod(Ni)
        gtau=np.zeros((N,), order='F')
        for x,i in zip(itt.product(*tau), \
                   itt.product(*[range(nx) for nx in Ni])):
            gtau[getidx(i,Ni)]=g(*x)
        return gtau

    def chkapprx(tau,f,g):
        for x in itt.product(*tau):
            fx,gx=f(x),g(*x)
            print("%+13.6e %+13.6e %+7.2f"%(fx,gx,(1.-fx/gx)*100.))

    g3= lambda x,y,z: np.sin(x)/(x+.1)+(z+1.)*y**2+x*np.exp(z)+.25
    g2= lambda x,y: g3(x,y,0.)
    g1= lambda x: g2(x,0)

    print(" --- test 1D poly ---")
    f1=fd.inFd(P1)
    tau=[]
    tau.append([0., .5, 1.]) # quadratic
    gtau=cmpgtau(tau,g1)
    f1.cmpcoef(tau,gtau)
    chkapprx(tau,f1,g1)

    print(" --- test 2D poly ---")
    f2=fd.inFd(V1)
    tau.append([0.,     1.]) # linear
    gtau=cmpgtau(tau,g2)
    f2.cmpcoef(tau,gtau)
    chkapprx(tau,f2,g2)

    print(" --- test 2D Bspline ---")
    h2=fd.inFd(V2)
    tau=[]
    for V in V2: tau.append( V.cmptau() )
    gtau=cmpgtau(tau,g2)
    h2.cmpcoef(tau,gtau)
    chkapprx(tau,h2,g2)

    print(" --- test 3D poly ---")
    f3=fd.inFd(V3)
    tau=[]
    tau.append([0., .5, 1.]) # quadratic
    tau.append([0.,     1.]) # linear
    tau.append([0.,     1.]) # linear
    gtau=cmpgtau(tau,g3)
    f3.cmpcoef(tau,gtau)
    chkapprx(tau,f3,g3)

    print(" --- test 3D mix ---")
    f4=fd.inFd(V4)
    tau=[]
    tau.append([0., .5, 1.]) # quadratic
    tau.append( S2.cmptau() )
    tau.append([0.,     1.]) # linear
    gtau=cmpgtau(tau,g3)
    f4.cmpcoef(tau,gtau)
    chkapprx(tau,f4,g3)

    """
    # test S - Ex. IX.2
    # the suggested tau sequence doesn't satisfy the W-S. theorem
    # I don't understand...
    def f(x): return (x-3.)*(x-6.)*(x-9.)
    k,t=4,np.zeros((17,), order='F')
    n=len(t)-k
    t[0:k]=0.; t[n:]=10.
    for i in range(1,10): t[k+i-1]=float(i)
    s1=s.inS( S(k,t) )
    tau =np.array([s1.t[i+2] for i in range(s1.n)], order='F')
    # to fulfill S-W theorem:
    tau[1] =.5*(tau[1] +tau[2] )
    tau[11]=.5*(tau[10]+tau[11])
    gtau=np.array([f(i) for i in tau], order='F')
    s1.cmpcoef(tau,gtau)
    print(' I   TAU(I)    TAVE(I)   G(TAU(I))  F(TAVE(I))    BCOEF(I)')
    for i in range(n):
        print('%2d  %8.5f  %8.5f  %10.5f  %10.5f  %10.5f'% \
              (i,tau[i],s1.tstar[i],gtau[i],f(s1.tstar[i]),s1.coef[i]))
    print("Bspline representation:"); print(s1)
    print("PP representation:");      print(s1.s2pp())

    # test PP
    k,xi=4,[0.,1.,3.,4.,6.]
    p1=pp(k,xi)
    p1.cmpcoef(gtau=[-1.,1.,2.5,3.,-4.5])
    print p1
    x=np.linspace(xi[0],xi[-1],7)
    y=[p1(xi) for xi in x]
    print("x="+", ".join(["%+.3f"%xi for xi in x]))
    print("y="+", ".join(["%+.3f"%yi for yi in y]))

    # test P with newton form
    def g(x): return (x**2-.5)
    k,xi=3,[-1.,1.]
    p2=pp(k,xi)
    print p2
    tau=np.linspace(xi[0],xi[1],k)
    gtau=[g(i) for i in tau]
    p2.cmpcoef(tau=tau,gtau=gtau)
    xt=np.linspace(xi[0],xi[1],k*4)
    print(" I   XT    F(XT)   G(XT)")
    for i in range(len(xt)):
        print("%2d  %+.3f  %+5.3f  %+5.3f"% (i,xt[i],g(xt[i]),p2(xt[i])))

    # test by doctest
    import doctest
    doctest.testmod()
    """



