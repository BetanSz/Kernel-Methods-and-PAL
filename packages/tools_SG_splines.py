# --*-- coding:utf-8 --*--
"""
Some auxiliary spline functions.

Includes all the knot vector of the thesis.
"""

from tools_SG_common import *
from packages.pppack import *
import scipy.interpolate as spi


class interp1d_picklable:
    """ class wrapper for piecewise linear function
    """

    def __init__(self, xi, yi, flag):
        self.xi = xi
        self.yi = yi
        self.flag = flag
        self.f = spi.interp1d(xi, yi, flag, bounds_error=False)

    def __call__(self, xnew):
        # for delivering a value in the same format as the other functions
        return self.f(xnew)[0][0]

    def __getstate__(self):
        return self.xi, self.yi, self.flag

    def __setstate__(self, state):
        self.f = spi.interp1d(state[0], state[1], state[2], bounds_error=False)


def knot_vector(tau_vec, k_vec, method_vec):
    """
    Produced the knot vector. In general:

    """

    t = []
    for x, k, mthd in zip(tau_vec, k_vec, method_vec):
        #In genera:
        # available information is n=len(x)
        # so knot vector will len(t)=n+k
        # so actual free knots will be len(t)-2*k=n-k
        # this free knots are taken equal to the breaks, and in the region of high burnup (which is easer) some tau points are not considered
        # len(t)=len(x)+k
        # free_knots=len(t)-2*k
        if mthd == 'linear':
            r"""
            for k=2
            knot vector of the form:
            \vec(t)= a=t_1=t_2=..=t_k<tau_2<...<tau_{n-1}<t_n=t_{n+1}=..=t_{n+k}=b
            For k=2 an asignation t_i=tau_i is possible fullfilling all degrees of freedom available in t.
            for k>2 this will generate an error, as the resulting tau vector is too long for the n avaialable degrees of freedom
            """
            t.append(np.concatenate([[x[0] for r in range(k)],
                                     cp.copy(x[1:-1]), [x[-1] for r in range(k)]]))

        if mthd == 'oddend':
            """
            For K>2
            Knots are lost at the end
            """
            # if simple assignement and the points before 'b' are lost then
            # lost=len(x)-free_knots=k
            # since x[0] and x[-1] are already taken, the vector of interest is x[1:-1], so -2 to proper indexation
            # Loosing the avant-last data sites generates the not-a-knot condition (page 182)
            lost = k - 2
            t.append(np.concatenate([[x[0] for r in range(k)], cp.copy(
                x[1:-1 - lost]), [x[-1] for r in range(k)]]))


        if mthd == 'odd_end_M1_':
            """
            For K>2
            Knots are lost at the end
            """
            lost = k - 2
            t.append(np.concatenate([[x[0] for r in range(k)], cp.copy(
                x[1:-1 - lost]), [x[-1] for r in range(k)]]))
            pos = k
            t[0] = np.delete(t[0], pos)
            t[0] = np.insert(t[0], pos, t[0][pos])

        if mthd == 'odd_end_M2_':
            """
            For K>2
            Knots are lost at the end
            """
            lost = k - 2
            t.append(np.concatenate([[x[0] for r in range(k)], cp.copy(
                x[1:-1 - lost]), [x[-1] for r in range(k)]]))
            pos = k + 1
            t[0] = np.delete(t[0], pos)
            t[0] = np.insert(t[0], pos, t[0][pos])

        if mthd == 'odd_end_M3_':
            """
            For K>2
            Knots are lost at the end
            """
            lost = k - 2
            t.append(np.concatenate([[x[0] for r in range(k)], cp.copy(
                x[1:-1 - lost]), [x[-1] for r in range(k)]]))
            pos = k + 2
            t[0] = np.delete(t[0], pos)
            t[0] = np.insert(t[0], pos, t[0][pos])


        if mthd == 'odd_end_M4_':
            """
            For K>2
            Knots are lost at the end
            """
            lost = k - 2
            t.append(np.concatenate([[x[0] for r in range(k)], cp.copy(
                x[1:-1 - lost]), [x[-1] for r in range(k)]]))
            pos = k + 3
            t[0] = np.delete(t[0], pos)
            t[0] = np.insert(t[0], pos, t[0][pos])

        if mthd == "preplop":
            """
            This is the first guess of the splop program that is supposed to be close to optimal
            """
            n = len(x)
            t_aux = []
            for i in range(n - k):
                # note that a +1 id added due to strange python indexation
                t_aux.append(sum(x[i + 1:i + k - 1 + 1]) / (k - 1))
            t.append(np.concatenate([[x[0] for a in range(k)], t_aux, [x[-1] for b in range(k)]]))

        if mthd == "Ecostum4":
            """
            The problem here is the smoothness condition. That can be played arround with knot multiplicty. Take the pre_plop array and custom the multiplicty to the XS at hand.
            However, I think due to the SW theorem the amount of inner knots with bigger than 1 mulitplicty is limited by the relative shift, so I can only pick one for K=3. If I go to higher K like 5 I did was able to put two, though at K=4 no, idk why. In any case even if indeed this multiplicty helps to deal with smoothens conditions, all the other high order smoothness conditions at the other points fuck up the results anyway.
            """
            n = len(x)
            t_aux = []
            for i in range(n - k):
                # note that a +1 id added due to strange python indexation
                t_aux.append(sum(x[i + 1:i + k - 1 + 1]) / (k - 1))
            t.append(np.concatenate([[x[0] for a in range(k)], t_aux, [x[-1] for b in range(k)]]))
            pos = 4
            t[0] = np.delete(t[0], pos + 1)
            t[0] = np.insert(t[0], pos, t[0][pos])

        if mthd == "Ecostum1":
            """
            The problem here is the smoothness condition. That can be played arround with knot multiplicty. Take the pre_plop array and custom the multiplicty to the XS at hand.
            """
            n = len(x)
            t_aux = []
            for i in range(n - k):
                # note that a +1 id added due to strange python indexation
                t_aux.append(sum(x[i + 1:i + k - 1 + 1]) / (k - 1))
            t.append(np.concatenate([[x[0] for a in range(k)], t_aux, [x[-1] for b in range(k)]]))
            pos = 3
            t[0] = np.delete(t[0], pos + 1)
            t[0] = np.insert(t[0], pos, t[0][pos])

        if mthd == "Ecostum3":
            """
            """
            n = len(x)
            t_aux = []
            for i in range(n - k):
                t_aux.append(sum(x[i + 1:i + k - 1 + 1]) / (k - 1))
            t.append(np.concatenate([[x[0] for a in range(k)], t_aux, [x[-1] for b in range(k)]]))

        if mthd == "Ecostum2":
            """
            This method considers the available data and increases multiplicity of the first points that are known to be problematic
            the idea here is to increase the multiplicty of the first nodes in a stable way as to not violate WS
            """
            lost = k - 2
            t_aux = np.concatenate([[x[0] for r in range(k)], cp.copy(
                x[1:-1 - lost]), [x[-1] for r in range(k)]])
            # I want linear in the first 5 points
            mult = k - 1
            up_to = 5
            # print t_aux
            # multiplicty should be possible.
            count = 0
            while up_to > 0:
                t_aux = np.insert(t_aux, count + k, t_aux[count + k])
                t_aux = np.delete(t_aux, count + k + 2)
                up_to = up_to - 1
                count = count + 2 + 1

            t.append(t_aux)

        if mthd == "splopM1":
            # old name "E_costum3"
            """
            Putting sistematicly K-1 knot in the first interior knot
            """
            n = len(x)
            if n >= k and k > 2:
                t.append(ppk.splopt(tau=x, n=n, k=k)[1])
            else:
                raise ValueError("invalid sequence: n=%i,k=%i,tau_i=%s" % (n, k, str(x)))
            pos = k
            t[0] = np.delete(t[0], pos + 1)
            t[0] = np.insert(t[0], pos, t[0][pos])

        if mthd == "splopM2":
            """
            Skipping the second knot all togheter. This gets complicated with splop as its less clear the values throughout thek not vector an violation of SW tends to happen.
            I will try now with vec_NaK
            """
            n = len(x)
            if n >= k and k > 2:
                t.append(ppk.splopt(tau=x, n=n, k=k)[1])
            else:
                raise ValueError("invalid sequence: n=%i,k=%i,tau_i=%s" % (n, k, str(x)))

            pos = k + 1
            t[0] = np.delete(t[0], pos + 1)
            t[0] = np.insert(t[0], pos, t[0][pos])

        if mthd == "splopM3":
            n = len(x)
            if n >= k and k > 2:
                t.append(ppk.splopt(tau=x, n=n, k=k)[1])
            else:
                raise ValueError("invalid sequence: n=%i,k=%i,tau_i=%s" % (n, k, str(x)))
            pos = k + 2
            t[0] = np.delete(t[0], pos + 1)
            t[0] = np.insert(t[0], pos, t[0][pos])

        if mthd == "splopM4":
            n = len(x)
            if n >= k and k > 2:
                t.append(ppk.splopt(tau=x, n=n, k=k)[1])
            else:
                raise ValueError("invalid sequence: n=%i,k=%i,tau_i=%s" % (n, k, str(x)))
            pos = k + 3
            t[0] = np.delete(t[0], pos + 1)
            t[0] = np.insert(t[0], pos, t[0][pos])

        if mthd == "splopM5":

            n = len(x)
            if n >= k and k > 2:
                t.append(ppk.splopt(tau=x, n=n, k=k)[1])
            else:
                raise ValueError("invalid sequence: n=%i,k=%i,tau_i=%s" % (n, k, str(x)))
            pos = k + 4
            t[0] = np.delete(t[0], pos + 1)
            t[0] = np.insert(t[0], pos, t[0][pos])

        if mthd == "splopM6":
            n = len(x)
            if n >= k and k > 2:
                t.append(ppk.splopt(tau=x, n=n, k=k)[1])
            else:
                raise ValueError("invalid sequence: n=%i,k=%i,tau_i=%s" % (n, k, str(x)))
            pos = k
            t[0] = np.delete(t[0], pos + 1)
            t[0] = np.insert(t[0], pos, t[0][pos])
            pos = k + 3
            t[0] = np.delete(t[0], pos + 1)
            t[0] = np.insert(t[0], pos, t[0][pos])

        if mthd == "splopM7":

            n = len(x)
            if n >= k and k > 2:
                t.append(ppk.splopt(tau=x, n=n, k=k)[1])
            else:
                raise ValueError("invalid sequence: n=%i,k=%i,tau_i=%s" % (n, k, str(x)))
            pos = k
            t[0] = np.delete(t[0], pos)
            t[0] = np.insert(t[0], pos, t[0][pos])

        if mthd == "splopM8":
            n = len(x)
            if n >= k and k > 2:
                t.append(ppk.splopt(tau=x, n=n, k=k)[1])
            else:
                raise ValueError("invalid sequence: n=%i,k=%i,tau_i=%s" % (n, k, str(x)))
            pos = k + 1
            t[0] = np.delete(t[0], pos)
            t[0] = np.insert(t[0], pos, t[0][pos])

        if mthd == "splopM9":
            n = len(x)
            if n >= k and k > 2:
                t.append(ppk.splopt(tau=x, n=n, k=k)[1])
            else:
                raise ValueError("invalid sequence: n=%i,k=%i,tau_i=%s" % (n, k, str(x)))
            pos = k + 2
            t[0] = np.delete(t[0], pos)
            t[0] = np.insert(t[0], pos, t[0][pos])

        if mthd == "splopM10":
            n = len(x)
            if n >= k and k > 2:
                t.append(ppk.splopt(tau=x, n=n, k=k)[1])
            else:
                raise ValueError("invalid sequence: n=%i,k=%i,tau_i=%s" % (n, k, str(x)))
            i = k
            while t[0][i] < 0.1:
                t[0][i] = t[0][i]
                t[0][i + 1] = t[0][i]
                i = i + k - 1

        if mthd == "etoile":
            # old name E_etoile
            """
            This is supposed to be I^* with very low ||I|| and thus low error bound, in theory. The data sites are consider the Greville sequence and the vector t is reconstructed from there.
            """
            t_aux = []
            [t_aux.append(x[0]) for a in range(k)]
            n = len(x)
            for i in range(n - k):
                t_aux.append(3*x[i+1]-sum(t_aux[i+2:i+2+1+1])) #K=4
            [t_aux.append(x[-1]) for a in range(k)]
            t.append(t_aux)

            # this implementation can produce decrecing knots
            for ti in t:
                check = ti[0]
                for tij in ti:
                    if tij < check:
                        raise ValueError('Decresing knot sequence found in etoile')
                    check = tij

        if mthd == "etoilemodif":
            # old name E_costum4 changes to etoile_modif
            """
            This is supposed to be I^* with very low ||I|| and thus low error bound, in theory. The data sites are consider the Greville sequence and the vector t is reconstructed from there.
            """
            t_aux = []
            [t_aux.append(x[0]) for a in range(k)]
            n = len(x)
            for i in range(n - k):
                t_aux.append((k - 1) * x[i + 1] - sum(t_aux[i + 2:i + 2 + (k - 3) + 1]))  # K
            [t_aux.append(x[-1]) for a in range(k)]
            t.append(t_aux)
            pos = 3
            t[0] = np.delete(t[0], 4)
            t[0] = np.insert(t[0], pos, t[0][pos])

        if mthd == "splop":
            """
            This is supposed to be the optimal knot placement detailed in page 193
            """
            n = len(x)
            if n >= k and k > 2:
                t.append(ppk.splopt(tau=np.array(x), n=n, k=k)[1])
            else:
                raise ValueError("invalid sequence: n=%i,k=%i,tau_i=%s" % (n, k, str(x)))

        if mthd == "notaknot":
            # changes from not_a_knot to notaknot
            """
            This is the not a knot spline for k=2m as presented in page 182. Basically the firsts and lasts knots are lost, preserving a,b.
            This is the complete spline approximation presented in page 201  apparently if information of the derivative is avaialbe.
            """
            t.append(vec_NaK(x, k))

        if mthd == "notaknotM1":
            # changes from not_a_knot to notaknot
            """
            This is the not a knot spline for k=2m as presented in page 182. Basically the firsts and lasts knots are lost, preserving a,b.
            This is the complete spline approximation presented in page 201  apparently if information of the derivative is avaialbe.
            """
            t.append(vec_NaK(x, k))

            pos = k
            for distance in [k]:
                for multy in range(1, k - 2):  # for multiplicity, if k=3 range(1,k-2)= 1 repetition
                    t[0] = np.delete(t[0], pos + 1)
                    t[0] = np.insert(t[0], pos, t[0][pos])
                pos = pos + distance


        if mthd == "notaknotM3":
            """
            Put several K-1 knots at the beguining, less spacing would violete WS
            """
            t.append(vec_NaK(x, k))
            pos = k
            for distance in [k, k, k, k, k]:
                if pos + 1 < len(t[0]) - k:

                    for multy in range(1, k - 2):  # for multiplicity, if k=3 range(1,k-2)= 1 repetition
                        t[0] = np.delete(t[0], pos + 1)
                        t[0] = np.insert(t[0], pos, t[0][pos])
                    pos = pos + distance

        if mthd == "notaknotM4":
            """
            The first third is double knot
            """
            t.append(vec_NaK(x, k))
            N_inner = len(t[0]) - 2 * k
            i = 0
            new_t = []
            for aux in t[0][0:k]:
                new_t.append(aux)
            j = 0
            howmuch = 3  # defines region of double knots, the bigger the smaller. The region is 'length of inner knots'/howmuch
            while i < int(N_inner / howmuch):
                new_t.append(t[0][j + k])
                new_t.append(t[0][j + k])
                j = j + 2
                i = i + 1
            for aux in t[0][j + k:len(t[0])]:
                new_t.append(aux)
            t = [np.array(new_t)]

        if mthd == "notaknotM10":
            """
            The first third is double knot
            """
            t.append(vec_NaK(x, k))
            i = k + 1
            while t[0][i] < 0.1:
                t[0][i] = t[0][i]
                t[0][i + 1] = t[0][i]
                t[0][i - 1] = t[0][i]
                i = i + k - 1

        if mthd == "notaknotM11":
            """
            The first third is double knot
            """
            t.append(vec_NaK(x, k))
            i = k + 1
            shift = -1
            while t[0][i] < 0.1:
                t[0][i] = t[0][i + shift]
                t[0][i + 1] = t[0][i + shift]
                t[0][i - 1] = t[0][i + shift]
                i = i + k - 1
                count = 0

        if mthd == "notaknotM7":
            """
            The first third is double knot
            """
            t.append(vec_NaK(x, k))
            N_inner = len(t[0]) - 2 * k
            pos = k
            t[0] = np.delete(t[0], pos + 1)
            t[0] = np.insert(t[0], pos, t[0][pos])
            pos = k + 1
            t[0] = np.delete(t[0], pos + 1)
            t[0] = np.insert(t[0], pos, t[0][pos])

        if mthd == "notaknotM9":
            """
            """
            t.append(vec_NaK(x, k))
            pos = k
            t[0] = np.delete(t[0], pos)
            t[0] = np.insert(t[0], pos, t[0][pos])
            pos = k
            t[0] = np.delete(t[0], pos + 2)
            t[0] = np.insert(t[0], pos, t[0][pos])


        if mthd == "notaknotM8":
            """
            The first third is double knot
            """
            t.append(vec_NaK(x, k))
            pos = k + 1
            t[0] = np.delete(t[0], pos)
            t[0] = np.insert(t[0], pos, t[0][pos])
            pos = k + 1
            t[0] = np.delete(t[0], pos + 2)
            t[0] = np.insert(t[0], pos, t[0][pos])

        if mthd == "notaknotM5":
            t.append(vec_NaK(x, k))
            pos = k + 1
            t[0] = np.delete(t[0], pos)
            pos = len(t[0]) - k - 2
            t[0] = np.insert(t[0], pos, t[0][pos])


        if mthd == "notaknotM6":
            # changes from not_a_knot to notaknot
            """
            """
            t.append(vec_NaK(x, k))
            pos = k + 1
            t[0] = np.delete(t[0], pos)
            t[0] = np.insert(t[0], pos, t[0][pos])
            pos = k + 1
            t[0] = np.delete(t[0], pos + 2)
            t[0] = np.insert(t[0], pos, t[0][pos])

        if mthd == "hermite":
            """
            This is the not a knot spline for k=2m as presented in page 182. Basically the firsts and lasts knots are lost, preserving a,b.
            This is the complete spline approximation presented in page 201  apparently if information of the derivative is avaialbe.
            """
            if k % 2 != 0:
                print("WARNING: k order is not even")
            n = len(x)
            m = int(k / 2)
            t_aux = []
            for i in range(n - k):
                t_aux.append(x[i + m])  # note that a +1 id added due to strange python indexation
                t_aux.append(x[i + m])  # note that a +1 id added due to strange python indexation
            t.append(np.concatenate([[x[0] for a in range(k)], t_aux, [x[-1] for b in range(k)]]))


    if len(t) == 0:
        raise ValueError('Could not find the method %s', mthd)
    return t


def vec_NaK(x, k):
    if k % 2 != 0:
        print("WARNING: k order is not even, exiting...")
        sys.exit()
    n = len(x)
    m = int(k / 2)

    t_aux = []
    for i in range(n - k):
        t_aux.append(x[i + m])  # note that a +1 id added due to strange python indexation

    return np.concatenate([[x[0] for a in range(k)], t_aux, [x[-1] for b in range(k)]])
