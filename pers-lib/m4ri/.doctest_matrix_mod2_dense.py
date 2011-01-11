# -*- coding: utf-8 -*-
from sage.all_cmdline import *; 
import sage.plot.plot; sage.plot.plot.DOCTEST_MODE=True

def warning_function(f):
    import warnings

    def doctest_showwarning(message, category, filename, lineno, file=f):
        try:
            file.write(warnings.formatwarning(message, category, 'doctest', lineno))
        except IOError:
            pass # the file (probably stdout) is invalid
    return doctest_showwarning

def change_warning_output(file):
    import warnings
    warnings.showwarning = warning_function(file)
cython(open('/home/malb/SAGE/devel/sage/sage/matrix/matrix_mod2_dense.pyx').read())

def example_0():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


Dense matrices over GF(2) using the M4RI library.

AUTHOR: Martin Albrecht <malb@informatik.uni-bremen.de>

EXAMPLES:
    >>> a = matrix(GF(Integer(2)),Integer(3),range(Integer(9)),sparse=False); a###line 7:_sage_    >>> a = matrix(GF(2),3,range(9),sparse=False); a
    [0 1 0]
    [1 0 1]
    [0 1 0]
    >>> a.rank()###line 11:_sage_    >>> a.rank()
    2
    >>> type(a)###line 13:_sage_    >>> type(a)
    <type 'sage.matrix.matrix_mod2_dense.Matrix_mod2_dense'>
    >>> a[Integer(0),Integer(0)] = Integer(1)###line 15:_sage_    >>> a[0,0] = 1
    >>> a.rank()###line 16:_sage_    >>> a.rank()
    3
    >>> parent(a)###line 18:_sage_    >>> parent(a)
    Full MatrixSpace of 3 by 3 dense matrices over Finite Field of size 2

    >>> a**Integer(2)###line 21:_sage_    >>> a^2
    [0 1 1]
    [1 0 0]
    [1 0 1]
    >>> a+a###line 25:_sage_    >>> a+a
    [0 0 0]
    [0 0 0]
    [0 0 0]

    >>> b = a.new_matrix(Integer(2),Integer(3),range(Integer(6))); b###line 30:_sage_    >>> b = a.new_matrix(2,3,range(6)); b
    [0 1 0]
    [1 0 1]

    >>> a*b###line 34:_sage_    >>> a*b
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for '*': 'Full MatrixSpace of 3 by 3 dense matrices over Finite Field of size 2' and 'Full MatrixSpace of 2 by 3 dense matrices over Finite Field of size 2'
    >>> b*a###line 38:_sage_    >>> b*a
    [1 0 1]
    [1 0 0]

    >>> a == loads(dumps(a))###line 42:_sage_    >>> a == loads(dumps(a))
    True
    >>> b == loads(dumps(b))###line 44:_sage_    >>> b == loads(dumps(b))
    True

    >>> a.echelonize(); a###line 47:_sage_    >>> a.echelonize(); a
    [1 0 0]
    [0 1 0]
    [0 0 1]
    >>> b.echelonize(); b###line 51:_sage_    >>> b.echelonize(); b
    [1 0 1]
    [0 1 0]

TESTS:
    >>> FF = FiniteField(Integer(2))###line 56:_sage_    >>> FF = FiniteField(2)
    >>> V = VectorSpace(FF,Integer(2))###line 57:_sage_    >>> V = VectorSpace(FF,2)
    >>> v = V([Integer(0),Integer(1)]); v###line 58:_sage_    >>> v = V([0,1]); v
    (0, 1)
    >>> W = V.subspace([v])###line 60:_sage_    >>> W = V.subspace([v])
    >>> W###line 61:_sage_    >>> W
    Vector space of degree 2 and dimension 1 over Finite Field of size 2
    Basis matrix:
    [0 1]
    >>> v in W###line 65:_sage_    >>> v in W
    True

    >>> M = Matrix(GF(Integer(2)), [[Integer(1),Integer(1),Integer(0)],[Integer(0),Integer(1),Integer(0)]])###line 68:_sage_    >>> M = Matrix(GF(2), [[1,1,0],[0,1,0]])
    >>> M.row_space()###line 69:_sage_    >>> M.row_space()
    Vector space of degree 3 and dimension 2 over Finite Field of size 2
    Basis matrix:
    [1 0 0]
    [0 1 0]

    >>> M = Matrix(GF(Integer(2)), [[Integer(1),Integer(1),Integer(0)],[Integer(0),Integer(0),Integer(1)]])###line 75:_sage_    >>> M = Matrix(GF(2), [[1,1,0],[0,0,1]])
    >>> M.row_space()###line 76:_sage_    >>> M.row_space()
    Vector space of degree 3 and dimension 2 over Finite Field of size 2
    Basis matrix:
    [1 1 0]
    [0 0 1]

TODO:
   - make linbox frontend and use it
     - charpoly ?
     - minpoly ?
   - make Matrix_modn_frontend and use it (?)
"""

def example_1():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


    Free global Gray code tables.
    """

def example_2():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


    Dense matrix over GF(2).
    """

def example_3():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


        Creates a new dense matrix over GF(2).

        INPUT:
            parent -- MatrixSpace (always)
            entries -- ignored
            copy -- ignored
            coerce -- ignored
            alloc -- if True a zero matrix is allocated (default:True)
        """

def example_4():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


        Dense matrix over GF(2) constructor.

        INPUT:
            parent -- MatrixSpace.
            entries -- may be list or 0 or 1
            copy -- ignored, elements are always copied
            coerce -- ignored, elements are always coerced to ints % 2

        EXAMPLES:
            >>> type(random_matrix(GF(Integer(2)),Integer(2),Integer(2)))###line 195:_sage_    >>> type(random_matrix(GF(2),2,2))
            <type 'sage.matrix.matrix_mod2_dense.Matrix_mod2_dense'>

            >>> Matrix(GF(Integer(2)),Integer(3),Integer(3),Integer(1))###line 198:_sage_    >>> Matrix(GF(2),3,3,1)
            [1 0 0]
            [0 1 0]
            [0 0 1]

            >>> Matrix(GF(Integer(2)),Integer(2),Integer(2),[Integer(1),Integer(1),Integer(1),Integer(0)])###line 203:_sage_    >>> Matrix(GF(2),2,2,[1,1,1,0])
            [1 1]
            [1 0]

            >>> Matrix(GF(Integer(2)),Integer(2),Integer(2),Integer(4))###line 207:_sage_    >>> Matrix(GF(2),2,2,4)
            [0 0]
            [0 0]

        TESTS:
            >>> Matrix(GF(Integer(2)),Integer(0),Integer(0))###line 212:_sage_    >>> Matrix(GF(2),0,0)
            []
            >>> Matrix(GF(Integer(2)),Integer(2),Integer(0))###line 214:_sage_    >>> Matrix(GF(2),2,0)
            []
            >>> Matrix(GF(Integer(2)),Integer(0),Integer(2))###line 216:_sage_    >>> Matrix(GF(2),0,2)
            []
        """

def example_5():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)



        Compares \code{self} with \code{right}. While equality and
        inequality are clearly defined, $<$ and $>$ are not.  For
        those first the matrix dimensions of \code{self} and
        \code{right} are compared. If these match then $<$ means that
        there is a position smallest (i,j) in \code{self} where
        \code{self[i,j]} is zero but \code{right[i,j]} is one. This
        (i,j) is smaller than the (i,j) if \code{self} and
        \code{right} are exchanged for the comparision.

        INPUT:
            right -- a matrix
            op -- comparison operation

        EXAMPLE:
            >>> A = random_matrix(GF(Integer(2)),Integer(2),Integer(2))###line 260:_sage_    >>> A = random_matrix(GF(2),2,2)
            >>> B = random_matrix(GF(Integer(2)),Integer(3),Integer(3))###line 261:_sage_    >>> B = random_matrix(GF(2),3,3)
            >>> A < B###line 262:_sage_    >>> A < B
            True
            >>> A = MatrixSpace(GF(Integer(2)),Integer(3),Integer(3))(Integer(1))###line 264:_sage_    >>> A = MatrixSpace(GF(2),3,3)(1)
            >>> B = MatrixSpace(GF(Integer(2)),Integer(3),Integer(3))(Integer(1))###line 265:_sage_    >>> B = MatrixSpace(GF(2),3,3)(1)
            >>> B[Integer(0),Integer(1)] = Integer(1)###line 266:_sage_    >>> B[0,1] = 1
            >>> A < B###line 267:_sage_    >>> A < B
            True

        TESTS:
            >>> A = matrix(GF(Integer(2)),Integer(2),Integer(0))###line 271:_sage_    >>> A = matrix(GF(2),2,0)
            >>> B = matrix(GF(Integer(2)),Integer(2),Integer(0))###line 272:_sage_    >>> B = matrix(GF(2),2,0)
            >>> A < B###line 273:_sage_    >>> A < B
            False
        """

def example_6():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


        The has of a matrix is computed as $\oplus i*a_i$ where the $a_i$ are
        the flattened entries in a matrix (by row, then by column).

        EXAMPLE:
            >>> B = random_matrix(GF(Integer(2)),Integer(3),Integer(3))###line 284:_sage_    >>> B = random_matrix(GF(2),3,3)
            >>> B.set_immutable()###line 285:_sage_    >>> B.set_immutable()
            >>> {B:Integer(0)} # indirect doctest###line 286:_sage_    >>> {B:0} # indirect doctest
            {[0 1 0]
            [0 1 1]
            [0 0 0]: 0}
            >>> M = random_matrix(GF(Integer(2)), Integer(123), Integer(321))###line 290:_sage_    >>> M = random_matrix(GF(2), 123, 321)
            >>> M.set_immutable()###line 291:_sage_    >>> M.set_immutable()
            >>> MZ = M.change_ring(ZZ)###line 292:_sage_    >>> MZ = M.change_ring(ZZ)
            >>> MZ.set_immutable()###line 293:_sage_    >>> MZ.set_immutable()
            >>> hash(M) == hash(MZ)###line 294:_sage_    >>> hash(M) == hash(MZ)
            True
            >>> MS = M.sparse_matrix()###line 296:_sage_    >>> MS = M.sparse_matrix()
            >>> MS.set_immutable()###line 297:_sage_    >>> MS.set_immutable()
            >>> hash(M) == hash(MS)###line 298:_sage_    >>> hash(M) == hash(MS)
            True

        TEST:
            >>> A = matrix(GF(Integer(2)),Integer(2),Integer(0))###line 302:_sage_    >>> A = matrix(GF(2),2,0)
            >>> hash(A)###line 303:_sage_    >>> hash(A)
            Traceback (most recent call last):
            ...
            TypeError: mutable matrices are unhashable
            >>> A.set_immutable()###line 307:_sage_    >>> A.set_immutable()
            >>> hash(A)###line 308:_sage_    >>> hash(A)
            0

        """

def example_7():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


        Set the (i,j) entry of self to the int value.
        """

def example_8():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


        EXAMPLE:
            >>> B = random_matrix(GF(Integer(2)),Integer(3),Integer(3))###line 416:_sage_    >>> B = random_matrix(GF(2),3,3)
            >>> B # indirect doctest###line 417:_sage_    >>> B # indirect doctest
            [0 1 0]
            [0 1 1]
            [0 0 0]
            >>> block_matrix([B, Integer(1), Integer(0), B])###line 421:_sage_    >>> block_matrix([B, 1, 0, B])
            [0 1 0|1 0 0]
            [0 1 1|0 1 0]
            [0 0 0|0 0 1]
            [-----+-----]
            [0 0 0|0 1 0]
            [0 0 0|0 1 1]
            [0 0 0|0 0 0]
        """

def example_9():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


        Matrix addition.

        INPUT:
            right -- matrix of dimension self.nrows() x self.ncols()

        EXAMPLES:
            >>> A = random_matrix(GF(Integer(2)),Integer(10),Integer(10))###line 486:_sage_    >>> A = random_matrix(GF(2),10,10)
            >>> A + A == Matrix(GF(Integer(2)),Integer(10),Integer(10),Integer(0))###line 487:_sage_    >>> A + A == Matrix(GF(2),10,10,0)
            True

            >>> A = random_matrix(GF(Integer(2)),Integer(257),Integer(253))###line 490:_sage_    >>> A = random_matrix(GF(2),257,253)
            >>> A + A == Matrix(GF(Integer(2)),Integer(257),Integer(253),Integer(0))###line 491:_sage_    >>> A + A == Matrix(GF(2),257,253,0)
            True

        TESTS:
            >>> A = matrix(GF(Integer(2)),Integer(2),Integer(0))###line 495:_sage_    >>> A = matrix(GF(2),2,0)
            >>> A+A###line 496:_sage_    >>> A+A
            []
            >>> A = matrix(GF(Integer(2)),Integer(0),Integer(2))###line 498:_sage_    >>> A = matrix(GF(2),0,2)
            >>> A+A###line 499:_sage_    >>> A+A
            []
            >>> A = matrix(GF(Integer(2)),Integer(0),Integer(0))###line 501:_sage_    >>> A = matrix(GF(2),0,0)
            >>> A+A###line 502:_sage_    >>> A+A
            []
        """

def example_10():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


        Matrix addition.

        INPUT:
            right -- matrix of dimension self.nrows() x self.ncols()

        EXAMPLES:
            >>> A = random_matrix(GF(Integer(2)),Integer(10),Integer(10))###line 521:_sage_    >>> A = random_matrix(GF(2),10,10)
            >>> A - A == Matrix(GF(Integer(2)),Integer(10),Integer(10),Integer(0))###line 522:_sage_    >>> A - A == Matrix(GF(2),10,10,0)
            True
        """

def example_11():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


        Matrix multiplication.

        ALGORITHM: Uses the 'Method of the Four Russians
        Multiplication', see self._multiply_m4rm.
        """

def example_12():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


        Multiply matrices using the 'Method of the Four Russians
        Multiplication' (M4RM) or Konrod's method.

        The algorithm is based on an algorithm by Arlazarov, Dinic,
        Kronrod, and Faradzev [ADKF70] and appeared in [AHU]. This
        implementation is based on a description given in Gregory
        Bard's 'Method of the Four Russians Inversion' paper [B06].

        INPUT:
            right -- Matrix
            k -- parameter $k$ for the Gray Code table size. If $k=0$ a
                 suitable value is chosen by the function.
                 ($0<= k <= 16$, default: 0)

        EXAMPLE:
              >>> A = Matrix(GF(Integer(2)), Integer(4), Integer(3), [Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(1)] )###line 557:_sage_    >>> A = Matrix(GF(2), 4, 3, [0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1] )
              >>> B = Matrix(GF(Integer(2)), Integer(3), Integer(4), [Integer(0), Integer(0), Integer(1), Integer(0), Integer(1), Integer(0), Integer(0), Integer(1), Integer(1), Integer(1), Integer(0), Integer(0)] )###line 558:_sage_    >>> B = Matrix(GF(2), 3, 4, [0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0] )
              >>> A###line 559:_sage_    >>> A
              [0 0 0]
              [0 1 0]
              [0 1 1]
              [0 0 1]
              >>> B###line 564:_sage_    >>> B
              [0 0 1 0]
              [1 0 0 1]
              [1 1 0 0]
              >>> A._multiply_m4rm(B, Integer(0))###line 568:_sage_    >>> A._multiply_m4rm(B, 0)
              [0 0 0 0]
              [1 0 0 1]
              [0 1 0 1]
              [1 1 0 0]

        TESTS:
            >>> A = random_matrix(GF(Integer(2)),Integer(0),Integer(0))###line 575:_sage_    >>> A = random_matrix(GF(2),0,0)
            >>> B = random_matrix(GF(Integer(2)),Integer(0),Integer(0))###line 576:_sage_    >>> B = random_matrix(GF(2),0,0)
            >>> A._multiply_m4rm(B, Integer(0))###line 577:_sage_    >>> A._multiply_m4rm(B, 0)
            []
            >>> A = random_matrix(GF(Integer(2)),Integer(3),Integer(0))###line 579:_sage_    >>> A = random_matrix(GF(2),3,0)
            >>> B = random_matrix(GF(Integer(2)),Integer(0),Integer(3))###line 580:_sage_    >>> B = random_matrix(GF(2),0,3)
            >>> A._multiply_m4rm(B, Integer(0))###line 581:_sage_    >>> A._multiply_m4rm(B, 0)
            [0 0 0]
            [0 0 0]
            [0 0 0]
            >>> A = random_matrix(GF(Integer(2)),Integer(0),Integer(3))###line 585:_sage_    >>> A = random_matrix(GF(2),0,3)
            >>> B = random_matrix(GF(Integer(2)),Integer(3),Integer(0))###line 586:_sage_    >>> B = random_matrix(GF(2),3,0)
            >>> A._multiply_m4rm(B, Integer(0))###line 587:_sage_    >>> A._multiply_m4rm(B, 0)
            []

        ALGORITHM: Uses the 'Method of the Four Russians'
        multiplication as implemented in the M4RI library.

        REFERENCES:
            [AHU] A. Aho, J. Hopcroft, and J. Ullman. 'Chapter 6:
                     Matrix Multiplication and Related Operations.'
                     The Design and Analysis of Computer
                     Algorithms. Addison-Wesley, 1974.

            [ADKF70] V. Arlazarov, E. Dinic, M. Kronrod, and
                     I. Faradzev. 'On Economical Construction of the
                     Transitive Closure of a Directed Graph.'
                     Dokl. Akad. Nauk. SSSR No. 194 (in Russian),
                     English Translation in Soviet Math Dokl. No. 11,
                     1970.

            [Bard06] G. Bard. 'Accelerating Cryptanalysis with the
                     Method of Four Russians'. Cryptography E-Print
                     Archive (http://eprint.iacr.org/2006/251.pdf),
                     2006.
        """

def example_13():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


        Classical $O(n^3)$ multiplication.

        This can be quite fast for matrix vector multiplication but
        the other routines fall back to this implementation in that
        case anyway.

        EXAMPLE:
              >>> A = Matrix(GF(Integer(2)), Integer(4), Integer(3), [Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(1)] )###line 637:_sage_    >>> A = Matrix(GF(2), 4, 3, [0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1] )
              >>> B = Matrix(GF(Integer(2)), Integer(3), Integer(4), [Integer(0), Integer(0), Integer(1), Integer(0), Integer(1), Integer(0), Integer(0), Integer(1), Integer(1), Integer(1), Integer(0), Integer(0)] )###line 638:_sage_    >>> B = Matrix(GF(2), 3, 4, [0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0] )
              >>> A###line 639:_sage_    >>> A
              [0 0 0]
              [0 1 0]
              [0 1 1]
              [0 0 1]
              >>> B###line 644:_sage_    >>> B
              [0 0 1 0]
              [1 0 0 1]
              [1 1 0 0]
              >>> A._multiply_classical(B)###line 648:_sage_    >>> A._multiply_classical(B)
              [0 0 0 0]
              [1 0 0 1]
              [0 1 0 1]
              [1 1 0 0]

        TESTS:
            >>> A = random_matrix(GF(Integer(2)),Integer(0),Integer(0))###line 655:_sage_    >>> A = random_matrix(GF(2),0,0)
            >>> B = random_matrix(GF(Integer(2)),Integer(0),Integer(0))###line 656:_sage_    >>> B = random_matrix(GF(2),0,0)
            >>> A._multiply_classical(B)###line 657:_sage_    >>> A._multiply_classical(B)
            []
            >>> A = random_matrix(GF(Integer(2)),Integer(3),Integer(0))###line 659:_sage_    >>> A = random_matrix(GF(2),3,0)
            >>> B = random_matrix(GF(Integer(2)),Integer(0),Integer(3))###line 660:_sage_    >>> B = random_matrix(GF(2),0,3)
            >>> A._multiply_classical(B)###line 661:_sage_    >>> A._multiply_classical(B)
            [0 0 0]
            [0 0 0]
            [0 0 0]
            >>> A = random_matrix(GF(Integer(2)),Integer(0),Integer(3))###line 665:_sage_    >>> A = random_matrix(GF(2),0,3)
            >>> B = random_matrix(GF(Integer(2)),Integer(3),Integer(0))###line 666:_sage_    >>> B = random_matrix(GF(2),3,0)
            >>> A._multiply_classical(B)###line 667:_sage_    >>> A._multiply_classical(B)
            []
        """

def example_14():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


        Strassen-Winograd $O(n^{2.807})$ multiplication [Str69].

        This implementation in M4RI is inspired by Sage's generic
        Strassen implementation [BHS08] but uses a more memory
        efficient operation schedule [DP08].

        The performance of this routine depends on the parameter
        cutoff. On many modern machines 2048 should give acceptable
        performance, a good rule of thumb for calculating the optimal
        cutoff would that two matrices of the cutoff size should fit
        in L2 cache, so: $cutoff = \sqrt{L2 * 8 * 1024^2 / 2}$ where
        $L2$ is the size of the L2 cache in MB.

        INPUT:
            right -- a matrix of matching dimensions.
            cutoff -- matrix dimension where M4RM should be used
                      instead of Strassen (default: let M4RI decide)

        EXAMPLE:
              >>> A = Matrix(GF(Integer(2)), Integer(4), Integer(3), [Integer(0), Integer(0), Integer(0), Integer(0), Integer(1), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(1)] )###line 698:_sage_    >>> A = Matrix(GF(2), 4, 3, [0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1] )
              >>> B = Matrix(GF(Integer(2)), Integer(3), Integer(4), [Integer(0), Integer(0), Integer(1), Integer(0), Integer(1), Integer(0), Integer(0), Integer(1), Integer(1), Integer(1), Integer(0), Integer(0)] )###line 699:_sage_    >>> B = Matrix(GF(2), 3, 4, [0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0] )
              >>> A###line 700:_sage_    >>> A
              [0 0 0]
              [0 1 0]
              [0 1 1]
              [0 0 1]
              >>> B###line 705:_sage_    >>> B
              [0 0 1 0]
              [1 0 0 1]
              [1 1 0 0]
              >>> A._multiply_strassen(B, Integer(0))###line 709:_sage_    >>> A._multiply_strassen(B, 0)
              [0 0 0 0]
              [1 0 0 1]
              [0 1 0 1]
              [1 1 0 0]
              >>> A = random_matrix(GF(Integer(2)),Integer(2701),Integer(3000))###line 714:_sage_    >>> A = random_matrix(GF(2),2701,3000)
              >>> B = random_matrix(GF(Integer(2)),Integer(3000),Integer(3172))###line 715:_sage_    >>> B = random_matrix(GF(2),3000,3172)
              >>> A._multiply_strassen(B, Integer(256)) == A._multiply_m4rm(B, Integer(0))###line 716:_sage_    >>> A._multiply_strassen(B, 256) == A._multiply_m4rm(B, 0)
              True

        TESTS:
            >>> A = random_matrix(GF(Integer(2)),Integer(0),Integer(0))###line 720:_sage_    >>> A = random_matrix(GF(2),0,0)
            >>> B = random_matrix(GF(Integer(2)),Integer(0),Integer(0))###line 721:_sage_    >>> B = random_matrix(GF(2),0,0)
            >>> A._multiply_strassen(B, Integer(0))###line 722:_sage_    >>> A._multiply_strassen(B, 0)
            []
            >>> A = random_matrix(GF(Integer(2)),Integer(3),Integer(0))###line 724:_sage_    >>> A = random_matrix(GF(2),3,0)
            >>> B = random_matrix(GF(Integer(2)),Integer(0),Integer(3))###line 725:_sage_    >>> B = random_matrix(GF(2),0,3)
            >>> A._multiply_strassen(B, Integer(0))###line 726:_sage_    >>> A._multiply_strassen(B, 0)
            [0 0 0]
            [0 0 0]
            [0 0 0]
            >>> A = random_matrix(GF(Integer(2)),Integer(0),Integer(3))###line 730:_sage_    >>> A = random_matrix(GF(2),0,3)
            >>> B = random_matrix(GF(Integer(2)),Integer(3),Integer(0))###line 731:_sage_    >>> B = random_matrix(GF(2),3,0)
            >>> A._multiply_strassen(B, Integer(0))###line 732:_sage_    >>> A._multiply_strassen(B, 0)
            []

        ALGORITHM: Uses Strassen-Winograd matrix multiplication with
        M4RM as base case as implemented in the M4RI library.

        REFERENCES:
            [Str69] Volker Strassen. Gaussian elimination is not
                    optimal. Numerische Mathematik, 13:354-356, 1969.

            [BHS08] Robert Bradshaw, David Harvey and William
                    Stein. strassen_window_multiply_c. strassen.pyx,
                    Sage 3.0, 2008. http://www.sagemath.org

            [DP08] Jean-Guillaume Dumas and Clement Pernet. Memory
                   efficient scheduling of Strassen-Winograd's matrix
                   multiplication algorithm. arXiv:0707.2347v1, 2008.
        """

def example_15():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


        EXAMPLES:
            >>> A = random_matrix(GF(Integer(2)),Integer(100),Integer(100))###line 766:_sage_    >>> A = random_matrix(GF(2),100,100)
            >>> A - A == A - -A###line 767:_sage_    >>> A - A == A - -A
            True
        """

def example_16():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


        Inverts self using the 'Method of the Four Russians'
        inversion.

        If \code{self} is not invertible a \code{ZeroDivisionError} is
        raised.

        EXAMPLE:
            >>> A = Matrix(GF(Integer(2)),Integer(3),Integer(3), [Integer(0), Integer(0), Integer(1), Integer(0), Integer(1), Integer(1), Integer(1), Integer(0), Integer(1)])###line 781:_sage_    >>> A = Matrix(GF(2),3,3, [0, 0, 1, 0, 1, 1, 1, 0, 1])
            >>> MS = A.parent()###line 782:_sage_    >>> MS = A.parent()
            >>> A###line 783:_sage_    >>> A
            [0 0 1]
            [0 1 1]
            [1 0 1]
            >>> ~A###line 787:_sage_    >>> ~A
            [1 0 1]
            [1 1 0]
            [1 0 0]
            >>> A * ~A == ~A * A == MS(Integer(1))###line 791:_sage_    >>> A * ~A == ~A * A == MS(1)
            True

        TESTS:
            >>> A = matrix(GF(Integer(2)),Integer(0),Integer(0))###line 795:_sage_    >>> A = matrix(GF(2),0,0)
            >>> A**(-Integer(1))###line 796:_sage_    >>> A^(-1)
            []
        """

def example_17():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


        Returns a copy of \code{self}.

        EXAMPLES:
             >>> MS = MatrixSpace(GF(Integer(2)),Integer(3),Integer(3))###line 828:_sage_    >>> MS = MatrixSpace(GF(2),3,3)
             >>> A = MS(Integer(1))###line 829:_sage_    >>> A = MS(1)
             >>> A.copy() == A###line 830:_sage_    >>> A.copy() == A
             True
             >>> A.copy() is A###line 832:_sage_    >>> A.copy() is A
             False

             >>> A = random_matrix(GF(Integer(2)),Integer(100),Integer(100))###line 835:_sage_    >>> A = random_matrix(GF(2),100,100)
             >>> A.copy() == A###line 836:_sage_    >>> A.copy() == A
             True
             >>> A.copy() is A###line 838:_sage_    >>> A.copy() is A
             False

             >>> A.echelonize()###line 841:_sage_    >>> A.echelonize()
             >>> A.copy() == A###line 842:_sage_    >>> A.copy() == A
             True

        """

def example_18():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


        Returns list of the elements of \code{self} in row major
        ordering.

        EXAMPLE:
            >>> A = Matrix(GF(Integer(2)),Integer(2),Integer(2),[Integer(1),Integer(0),Integer(1),Integer(1)])###line 863:_sage_    >>> A = Matrix(GF(2),2,2,[1,0,1,1])
            >>> A###line 864:_sage_    >>> A
            [1 0]
            [1 1]
            >>> A.list() #indirect doctest###line 867:_sage_    >>> A.list() #indirect doctest
            [1, 0, 1, 1]

        TESTS:
            >>> A = Matrix(GF(Integer(2)),Integer(3),Integer(0))###line 871:_sage_    >>> A = Matrix(GF(2),3,0)
            >>> A.list()###line 872:_sage_    >>> A.list()
            []
        """

def example_19():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


        Puts self in (reduced) row echelon form.

        INPUT:
            self -- a mutable matrix
            algorithm -- 'm4ri' -- uses M4RI (default)
                         'classical' -- uses classical Gaussian elimination
                         'pluq' -- uses PLUQ factorization
            k --  the parameter 'k' of the M4RI algorithm. It MUST be between
                  1 and 16 (inclusive). If it is not specified it will be calculated as
                  3/4 * log_2( min(nrows, ncols) ) as suggested in the M4RI paper.
            reduced -- return reduced row echelon form (default:True)

        EXAMPLE:
             >>> A = random_matrix(GF(Integer(2)), Integer(10), Integer(10))###line 908:_sage_    >>> A = random_matrix(GF(2), 10, 10)
             >>> B = A.copy(); B.echelonize() # fastest###line 909:_sage_    >>> B = A.copy(); B.echelonize() # fastest
             >>> C = A.copy(); C.echelonize(k=Integer(2)) # force k###line 910:_sage_    >>> C = A.copy(); C.echelonize(k=2) # force k
             >>> E = A.copy(); E.echelonize(algorithm='classical') # force Gaussian elimination###line 911:_sage_    >>> E = A.copy(); E.echelonize(algorithm='classical') # force Gaussian elimination
             >>> B == C == E###line 912:_sage_    >>> B == C == E
             True

        TESTS:
             >>> VF2 = VectorSpace(GF(Integer(2)),Integer(2))###line 916:_sage_    >>> VF2 = VectorSpace(GF(2),2)
             >>> WF2 = VF2.submodule([VF2([Integer(1),Integer(1)])])###line 917:_sage_    >>> WF2 = VF2.submodule([VF2([1,1])])
             >>> WF2###line 918:_sage_    >>> WF2
             Vector space of degree 2 and dimension 1 over Finite Field of size 2
             Basis matrix:
             [1 1]

             >>> A2 = matrix(GF(Integer(2)),Integer(2),[Integer(1),Integer(0),Integer(0),Integer(1)])###line 923:_sage_    >>> A2 = matrix(GF(2),2,[1,0,0,1])
             >>> A2.kernel()###line 924:_sage_    >>> A2.kernel()
             Vector space of degree 2 and dimension 0 over Finite Field of size 2
             Basis matrix:
             []

        ALGORITHM: Uses M4RI library

        REFERENCES:
            [Bard06] G. Bard. 'Accelerating Cryptanalysis with the Method of Four Russians'. Cryptography
                     E-Print Archive (http://eprint.iacr.org/2006/251.pdf), 2006.
        """

def example_20():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


        Returns the pivot columns of \code{self} if \code{self} is in
        row echelon form.

        EXAMPLE:
            >>> A = matrix(GF(Integer(2)),Integer(5),Integer(5),[Integer(0),Integer(1),Integer(0),Integer(1),Integer(0),Integer(0),Integer(1),Integer(0),Integer(1),Integer(1),Integer(0),Integer(1),Integer(0),Integer(1),Integer(0),Integer(0),Integer(0),Integer(0),Integer(1),Integer(0),Integer(0),Integer(0),Integer(1),Integer(0),Integer(1)])###line 1005:_sage_    >>> A = matrix(GF(2),5,5,[0,1,0,1,0,0,1,0,1,1,0,1,0,1,0,0,0,0,1,0,0,0,1,0,1])
            >>> E = A.echelon_form()###line 1006:_sage_    >>> E = A.echelon_form()
            >>> E###line 1007:_sage_    >>> E
            [0 1 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 0]
            [0 0 0 0 1]
            [0 0 0 0 0]
            >>> E._pivots()###line 1013:_sage_    >>> E._pivots()
            [1, 2, 3, 4]
        """

def example_21():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


        Randomize density proportion of the entries of this matrix,
        leaving the rest unchanged.

        EXAMPLES:
            >>> A = matrix(GF(Integer(2)), Integer(5), Integer(5), Integer(0))###line 1038:_sage_    >>> A = matrix(GF(2), 5, 5, 0)
            >>> A.randomize(RealNumber('0.5')); A###line 1039:_sage_    >>> A.randomize(0.5); A
            [0 0 0 1 1]
            [0 1 0 0 1]
            [1 0 0 0 0]
            [0 1 0 0 0]
            [0 0 0 1 0]
            >>> A.randomize(); A###line 1045:_sage_    >>> A.randomize(); A
            [0 0 1 1 0]
            [1 1 0 0 1]
            [1 1 1 1 0]
            [1 1 1 1 1]
            [0 0 1 1 0]

        TESTS:
        With the libc random number generator random(), we had problems
        where the ranks of all of these matrices would be the same
        (and they would all be far too low).  This verifies that the
        problem is gone, with Mersenne Twister.
            >>> MS2 = MatrixSpace(GF(Integer(2)), Integer(1000))###line 1057:_sage_    >>> MS2 = MatrixSpace(GF(2), 1000)
            >>> [MS2.random_element().rank() for i in range(Integer(5))]###line 1058:_sage_    >>> [MS2.random_element().rank() for i in range(5)]
            [999, 998, 1000, 999, 999]

        Testing corner case.
            >>> A = random_matrix(GF(Integer(2)),Integer(3),Integer(0))###line 1062:_sage_    >>> A = random_matrix(GF(2),3,0)
            >>> A###line 1063:_sage_    >>> A
            []
        """

def example_22():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


        EXAMPLE:
            >>> A = random_matrix(GF(Integer(2)),Integer(3),Integer(3)); A###line 1107:_sage_    >>> A = random_matrix(GF(2),3,3); A
            [0 1 0]
            [0 1 1]
            [0 0 0]
            >>> A.rescale_row(Integer(0),Integer(0),Integer(0)); A###line 1111:_sage_    >>> A.rescale_row(0,0,0); A
            [0 0 0]
            [0 1 1]
            [0 0 0]
        """

def example_23():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


        EXAMPLE:
            >>> A = random_matrix(GF(Integer(2)),Integer(3),Integer(3)); A###line 1123:_sage_    >>> A = random_matrix(GF(2),3,3); A
            [0 1 0]
            [0 1 1]
            [0 0 0]
            >>> A.add_multiple_of_row(Integer(0),Integer(1),Integer(1),Integer(0)); A###line 1127:_sage_    >>> A.add_multiple_of_row(0,1,1,0); A
            [0 0 1]
            [0 1 1]
            [0 0 0]
        """

def example_24():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


        EXAMPLE:
            >>> A = random_matrix(GF(Integer(2)),Integer(3),Integer(3))###line 1138:_sage_    >>> A = random_matrix(GF(2),3,3)
            >>> A###line 1139:_sage_    >>> A
            [0 1 0]
            [0 1 1]
            [0 0 0]
            >>> A.swap_rows(Integer(0),Integer(1)); A###line 1143:_sage_    >>> A.swap_rows(0,1); A
            [0 1 1]
            [0 1 0]
            [0 0 0]
        """

def example_25():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


        EXAMPLE:
            >>> A = random_matrix(GF(Integer(2)),Integer(3),Integer(3))###line 1153:_sage_    >>> A = random_matrix(GF(2),3,3)
            >>> A###line 1154:_sage_    >>> A
            [0 1 0]
            [0 1 1]
            [0 0 0]
            >>> A.swap_columns(Integer(0),Integer(1)); A###line 1158:_sage_    >>> A.swap_columns(0,1); A
            [1 0 0]
            [1 0 1]
            [0 0 0]

            >>> A = random_matrix(GF(Integer(2)),Integer(3),Integer(65))###line 1163:_sage_    >>> A = random_matrix(GF(2),3,65)

            >>> B = A.copy()###line 1165:_sage_    >>> B = A.copy()
            >>> B.swap_columns(Integer(0),Integer(1))###line 1166:_sage_    >>> B.swap_columns(0,1)
            >>> B.swap_columns(Integer(0),Integer(1))###line 1167:_sage_    >>> B.swap_columns(0,1)
            >>> A == B###line 1168:_sage_    >>> A == B
            True

            >>> A.swap_columns(Integer(0),Integer(64))###line 1171:_sage_    >>> A.swap_columns(0,64)
            >>> A.column(Integer(0)) == B.column(Integer(64))###line 1172:_sage_    >>> A.column(0) == B.column(64)
            True
            >>> A.swap_columns(Integer(63),Integer(64))###line 1174:_sage_    >>> A.swap_columns(63,64)
            >>> A.column(Integer(63)) == B.column(Integer(0))###line 1175:_sage_    >>> A.column(63) == B.column(0)
            True
        """

def example_26():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


        Returns a string of self in \Magma form. Does not return \Magma
        object but string.

        EXAMPLE:
            >>> A = random_matrix(GF(Integer(2)),Integer(3),Integer(3))###line 1188:_sage_    >>> A = random_matrix(GF(2),3,3)


            'Matrix(GF(2),3,3,StringToIntegerSequence("0 1 0 0 1 1 0 0 0"))'
            >>> A = random_matrix(GF(Integer(2)),Integer(100),Integer(100))###line 1191:_sage_    >>> A = random_matrix(GF(2),100,100)
            >>> B = random_matrix(GF(Integer(2)),Integer(100),Integer(100))###line 1192:_sage_    >>> B = random_matrix(GF(2),100,100)


            True

        TESTS:
            >>> A = random_matrix(GF(Integer(2)),Integer(0),Integer(3))###line 1197:_sage_    >>> A = random_matrix(GF(2),0,3)


            Matrix with 0 rows and 3 columns
            >>> A = matrix(GF(Integer(2)),Integer(2),Integer(3),[Integer(0),Integer(1),Integer(1),Integer(1),Integer(0),Integer(0)])###line 1200:_sage_    >>> A = matrix(GF(2),2,3,[0,1,1,1,0,0])


            'Matrix(GF(2),2,3,StringToIntegerSequence("0 1 1 1 0 0"))'


            [0 1 1]
            [1 0 0]
        """

def example_27():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


        Return the determinant of this matrix over GF(2).

        EXAMPLES:
            >>> matrix(GF(Integer(2)),Integer(2),[Integer(1),Integer(1),Integer(0),Integer(1)]).determinant()###line 1216:_sage_    >>> matrix(GF(2),2,[1,1,0,1]).determinant()
            1
            >>> matrix(GF(Integer(2)),Integer(2),[Integer(1),Integer(1),Integer(1),Integer(1)]).determinant()###line 1218:_sage_    >>> matrix(GF(2),2,[1,1,1,1]).determinant()
            0
        """

def example_28():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


        Returns transpose of self and leaves self untouched.

        EXAMPLE:
            >>> A = Matrix(GF(Integer(2)),Integer(3),Integer(5),[Integer(1), Integer(0), Integer(1), Integer(0), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(0), Integer(1), Integer(1), Integer(0), Integer(1), Integer(0)])###line 1230:_sage_    >>> A = Matrix(GF(2),3,5,[1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0])
            >>> A###line 1231:_sage_    >>> A
            [1 0 1 0 0]
            [0 1 1 0 0]
            [1 1 0 1 0]
            >>> B = A.transpose(); B###line 1235:_sage_    >>> B = A.transpose(); B
            [1 0 1]
            [0 1 1]
            [1 1 0]
            [0 0 1]
            [0 0 0]
            >>> B.transpose() == A###line 1241:_sage_    >>> B.transpose() == A
            True

        TESTS:
            >>> A = random_matrix(GF(Integer(2)),Integer(0),Integer(40))###line 1245:_sage_    >>> A = random_matrix(GF(2),0,40)
            >>> A.transpose()###line 1246:_sage_    >>> A.transpose()
            40 x 0 dense matrix over Finite Field of size 2
        """

def example_29():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


        Augments \code{self} with \code{right}.

        EXAMPLE:
            >>> MS = MatrixSpace(GF(Integer(2)),Integer(3),Integer(3))###line 1269:_sage_    >>> MS = MatrixSpace(GF(2),3,3)
            >>> A = MS([Integer(0), Integer(1), Integer(0), Integer(1), Integer(1), Integer(0), Integer(1), Integer(1), Integer(1)]); A###line 1270:_sage_    >>> A = MS([0, 1, 0, 1, 1, 0, 1, 1, 1]); A
            [0 1 0]
            [1 1 0]
            [1 1 1]
            >>> B = A.augment(MS(Integer(1))); B###line 1274:_sage_    >>> B = A.augment(MS(1)); B
            [0 1 0 1 0 0]
            [1 1 0 0 1 0]
            [1 1 1 0 0 1]
            >>> B.echelonize(); B###line 1278:_sage_    >>> B.echelonize(); B
            [1 0 0 1 1 0]
            [0 1 0 1 0 0]
            [0 0 1 0 1 1]
            >>> C = B.matrix_from_columns([Integer(3),Integer(4),Integer(5)]); C###line 1282:_sage_    >>> C = B.matrix_from_columns([3,4,5]); C
            [1 1 0]
            [1 0 0]
            [0 1 1]
            >>> C == ~A###line 1286:_sage_    >>> C == ~A
            True
            >>> C*A == MS(Integer(1))###line 1288:_sage_    >>> C*A == MS(1)
            True

        TESTS:
            >>> A = random_matrix(GF(Integer(2)),Integer(2),Integer(3))###line 1292:_sage_    >>> A = random_matrix(GF(2),2,3)
            >>> B = random_matrix(GF(Integer(2)),Integer(2),Integer(0))###line 1293:_sage_    >>> B = random_matrix(GF(2),2,0)
            >>> A.augment(B)###line 1294:_sage_    >>> A.augment(B)
            [0 1 0]
            [0 1 1]
            >>> B.augment(A)###line 1297:_sage_    >>> B.augment(A)
            [0 1 0]
            [0 1 1]
            >>> M = Matrix(GF(Integer(2)), Integer(0), Integer(0), Integer(0))###line 1300:_sage_    >>> M = Matrix(GF(2), 0, 0, 0)
            >>> N = Matrix(GF(Integer(2)), Integer(0), Integer(19), Integer(0))###line 1301:_sage_    >>> N = Matrix(GF(2), 0, 19, 0)
            >>> W = M.augment(N)###line 1302:_sage_    >>> W = M.augment(N)
            >>> W.ncols()###line 1303:_sage_    >>> W.ncols()
            19
            >>> M = Matrix(GF(Integer(2)), Integer(0), Integer(1), Integer(0))###line 1305:_sage_    >>> M = Matrix(GF(2), 0, 1, 0)
            >>> N = Matrix(GF(Integer(2)), Integer(0), Integer(1), Integer(0))###line 1306:_sage_    >>> N = Matrix(GF(2), 0, 1, 0)
            >>> M.augment(N)###line 1307:_sage_    >>> M.augment(N)
            []
        """

def example_30():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


        Stack \code{self} on top of \code{other}.

        EXAMPLE:
            >>> A = matrix(GF(Integer(2)),Integer(2),Integer(2),[Integer(1),Integer(0),Integer(0),Integer(1)])###line 1331:_sage_    >>> A = matrix(GF(2),2,2,[1,0,0,1])
            >>> B = matrix(GF(Integer(2)),Integer(2),Integer(2),[Integer(0),Integer(1),Integer(1),Integer(0)])###line 1332:_sage_    >>> B = matrix(GF(2),2,2,[0,1,1,0])
            >>> A.stack(B)###line 1333:_sage_    >>> A.stack(B)
            [1 0]
            [0 1]
            [0 1]
            [1 0]
            >>> B.stack(A)###line 1338:_sage_    >>> B.stack(A)
            [0 1]
            [1 0]
            [1 0]
            [0 1]

        TESTS:
            >>> A = random_matrix(GF(Integer(2)),Integer(0),Integer(3))###line 1345:_sage_    >>> A = random_matrix(GF(2),0,3)
            >>> B = random_matrix(GF(Integer(2)),Integer(3),Integer(3))###line 1346:_sage_    >>> B = random_matrix(GF(2),3,3)
            >>> A.stack(B)###line 1347:_sage_    >>> A.stack(B)
            [0 1 0]
            [0 1 1]
            [0 0 0]
            >>> B.stack(A)###line 1351:_sage_    >>> B.stack(A)
            [0 1 0]
            [0 1 1]
            [0 0 0]
            >>> M = Matrix(GF(Integer(2)), Integer(0), Integer(0), Integer(0))###line 1355:_sage_    >>> M = Matrix(GF(2), 0, 0, 0)
            >>> N = Matrix(GF(Integer(2)), Integer(19), Integer(0), Integer(0))###line 1356:_sage_    >>> N = Matrix(GF(2), 19, 0, 0)
            >>> W = M.stack(N)###line 1357:_sage_    >>> W = M.stack(N)
            >>> W.nrows()###line 1358:_sage_    >>> W.nrows()
            19
            >>> M = Matrix(GF(Integer(2)), Integer(1), Integer(0), Integer(0))###line 1360:_sage_    >>> M = Matrix(GF(2), 1, 0, 0)
            >>> N = Matrix(GF(Integer(2)), Integer(1), Integer(0), Integer(0))###line 1361:_sage_    >>> N = Matrix(GF(2), 1, 0, 0)
            >>> M.stack(N)###line 1362:_sage_    >>> M.stack(N)
            []
        """

def example_31():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


        Return submatrix from the index lowr,lowc (inclusive) with
        dimension nrows x ncols.

        INPUT:
            lowr -- index of start row
            lowc -- index of start column
            nrows -- number of rows of submatrix
            ncols -- number of columns of submatrix

        EXAMPLES:
             >>> A = random_matrix(GF(Integer(2)),Integer(200),Integer(200))###line 1392:_sage_    >>> A = random_matrix(GF(2),200,200)
             >>> A[Integer(0):Integer(2),Integer(0):Integer(2)] == A.submatrix(Integer(0),Integer(0),Integer(2),Integer(2))###line 1393:_sage_    >>> A[0:2,0:2] == A.submatrix(0,0,2,2)
             True
             >>> A[Integer(0):Integer(100),Integer(0):Integer(100)] == A.submatrix(Integer(0),Integer(0),Integer(100),Integer(100))###line 1395:_sage_    >>> A[0:100,0:100] == A.submatrix(0,0,100,100)
             True
             >>> A == A.submatrix(Integer(0),Integer(0),Integer(200),Integer(200))###line 1397:_sage_    >>> A == A.submatrix(0,0,200,200)
             True

             >>> A[Integer(1):Integer(3),Integer(1):Integer(3)] == A.submatrix(Integer(1),Integer(1),Integer(2),Integer(2))###line 1400:_sage_    >>> A[1:3,1:3] == A.submatrix(1,1,2,2)
             True
             >>> A[Integer(1):Integer(100),Integer(1):Integer(100)] == A.submatrix(Integer(1),Integer(1),Integer(99),Integer(99))###line 1402:_sage_    >>> A[1:100,1:100] == A.submatrix(1,1,99,99)
             True
             >>> A[Integer(1):Integer(200),Integer(1):Integer(200)] == A.submatrix(Integer(1),Integer(1),Integer(199),Integer(199))###line 1404:_sage_    >>> A[1:200,1:200] == A.submatrix(1,1,199,199)
             True
        """

def example_32():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


        Serialize \code{self}.

        EXAMPLE:
            >>> A = random_matrix(GF(Integer(2)),Integer(10),Integer(10))###line 1440:_sage_    >>> A = random_matrix(GF(2),10,10)
            >>> f,s = A.__reduce__()###line 1441:_sage_    >>> f,s = A.__reduce__()
            >>> f(*s) == A###line 1442:_sage_    >>> f(*s) == A
            True
        """

def example_33():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


        Return space separated string of the entries in this matrix.

        EXAMPLES:
            >>> w = matrix(GF(Integer(2)),Integer(2),Integer(3),[Integer(1),Integer(0),Integer(1),Integer(1),Integer(1),Integer(0)])###line 1474:_sage_    >>> w = matrix(GF(2),2,3,[1,0,1,1,1,0])
            >>> w._export_as_string()###line 1475:_sage_    >>> w._export_as_string()
            '1 0 1 1 1 0'
        """

def example_34():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


        Return the density of this matrix.

        By density we understand the ration of the number of nonzero
        positions and the self.nrows() * self.ncols(), i.e. the number
        of possible nonzero positions.

        INPUT:
            approx -- return floating point approximation (default: False)

        EXAMPLE:
            >>> A = random_matrix(GF(Integer(2)),Integer(1000),Integer(1000))###line 1512:_sage_    >>> A = random_matrix(GF(2),1000,1000)
            >>> d = A.density(); d###line 1513:_sage_    >>> d = A.density(); d
            62483/125000
            >>> float(d)###line 1515:_sage_    >>> float(d)
            0.499863...
            >>> A.density(approx=True)###line 1517:_sage_    >>> A.density(approx=True)
            0.499773...
            >>> float(len(A.nonzero_positions())/Integer(1000)**Integer(2))###line 1519:_sage_    >>> float(len(A.nonzero_positions())/1000^2)
            0.499863...
        """

def example_35():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


        Return the rank of this matrix.

        EXAMPLE:
            >>> A = random_matrix(GF(Integer(2)), Integer(1000), Integer(1000))###line 1533:_sage_    >>> A = random_matrix(GF(2), 1000, 1000)
            >>> A.rank()###line 1534:_sage_    >>> A.rank()
            999

            >>> A = matrix(GF(Integer(2)),Integer(10), Integer(0))###line 1537:_sage_    >>> A = matrix(GF(2),10, 0)
            >>> A.rank()###line 1538:_sage_    >>> A.rank()
            0
        """

def example_36():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


        Return the right kernel of this matrix, as a vector
        space. This is the space of vectors x such that self\*x=0.

        INPUT: ignored.

        By convention if self has 0 rows, the kernel is of dimension
        0, whereas the kernel is whole domain if self has 0 columns.


        EXAMPLE::

            >>> A = random_matrix(GF(Integer(2)), Integer(10), Integer(10))###line 1570:_sage_    >>> A = random_matrix(GF(2), 10, 10)
            >>> A.rank()###line 1571:_sage_    >>> A.rank()
            >>> A.right_kernel()###line 1572:_sage_    >>> A.right_kernel()

        ALGORITHM: PLUQ Factorization.
        """

def example_37():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


    Returns the parity of the number of bits in a.

    EXAMPLES:
        >>> from sage.matrix.matrix_mod2_dense import parity###line 1622:_sage_    >>> from sage.matrix.matrix_mod2_dense import parity
        >>> parity(Integer(1))###line 1623:_sage_    >>> parity(1)
        1L
        >>> parity(Integer(3))###line 1625:_sage_    >>> parity(3)
        0L
        >>> parity(Integer(0x10000101011))###line 1627:_sage_    >>> parity(0x10000101011)
        1L
    """

def example_38():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


    Deserialize a matrix encoded in the string \code{s}.

    INPUT:
        r -- number of rows of matrix
        c -- number of columns of matrix
        s -- a string
        size -- length of the string s

    EXAMPLE:
        >>> A = random_matrix(GF(Integer(2)),Integer(100),Integer(101))###line 1651:_sage_    >>> A = random_matrix(GF(2),100,101)
        >>> _,(r,c,s,s2) = A.__reduce__()###line 1652:_sage_    >>> _,(r,c,s,s2) = A.__reduce__()
        >>> from sage.matrix.matrix_mod2_dense import unpickle_matrix_mod2_dense_v1###line 1653:_sage_    >>> from sage.matrix.matrix_mod2_dense import unpickle_matrix_mod2_dense_v1
        >>> unpickle_matrix_mod2_dense_v1(r,c,s,s2) == A###line 1654:_sage_    >>> unpickle_matrix_mod2_dense_v1(r,c,s,s2) == A
        True
        >>> loads(dumps(A)) == A###line 1656:_sage_    >>> loads(dumps(A)) == A
        True
    """

def example_39():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


    Returns a dense matrix over GF(2) from a 1-bit PNG image read from
    \code{filename}. No attempt is made to verify that the filname string
    actually points to a PNG image.

    INPUT:
        filename -- a string

    EXAMPLE:
        >>> from sage.matrix.matrix_mod2_dense import from_png, to_png###line 1699:_sage_    >>> from sage.matrix.matrix_mod2_dense import from_png, to_png
        >>> A = random_matrix(GF(Integer(2)),Integer(10),Integer(10))###line 1700:_sage_    >>> A = random_matrix(GF(2),10,10)
        >>> fn = tmp_filename()###line 1701:_sage_    >>> fn = tmp_filename()
        >>> to_png(A, fn)###line 1702:_sage_    >>> to_png(A, fn)
        >>> B = from_png(fn)###line 1703:_sage_    >>> B = from_png(fn)
        >>> A == B###line 1704:_sage_    >>> A == B
        True
    """

def example_40():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


    Saves the matrix \code{A} to filename as a 1-bit PNG image.

    INPUT:
        A -- a matrix over GF(2)
        filename -- a string for a file in a writeable position

    EXAMPLE:
        >>> from sage.matrix.matrix_mod2_dense import from_png, to_png###line 1741:_sage_    >>> from sage.matrix.matrix_mod2_dense import from_png, to_png
        >>> A = random_matrix(GF(Integer(2)),Integer(10),Integer(10))###line 1742:_sage_    >>> A = random_matrix(GF(2),10,10)
        >>> fn = tmp_filename()###line 1743:_sage_    >>> fn = tmp_filename()
        >>> to_png(A, fn)###line 1744:_sage_    >>> to_png(A, fn)
        >>> B = from_png(fn)###line 1745:_sage_    >>> B = from_png(fn)
        >>> A == B###line 1746:_sage_    >>> A == B
        True
    """

def example_41():	r""">>> set_random_seed(0L)

>>> change_warning_output(sys.stdout)


    Return PLUQ factorization of A.

    INPUT:
        A -- matrix
        algorithm -- 'standard' asymptotically fast (default)
                     'mmpf' M4RI inspired
                     'naive' naive cubic
        param -- either k for 'mmpf' is chosen or matrix multiplication
                 cutoff for 'standard' (default: 0)

    EXAMPLE:
        >>> from sage.matrix.matrix_mod2_dense import pluq###line 1782:_sage_    >>> from sage.matrix.matrix_mod2_dense import pluq
        >>> A = random_matrix(GF(Integer(2)),Integer(4),Integer(4)); A###line 1783:_sage_    >>> A = random_matrix(GF(2),4,4); A
        [0 1 0 1]
        [0 1 1 1]
        [0 0 0 1]
        [0 1 1 0]

        >>> LU, P, Q = pluq(A)###line 1789:_sage_    >>> LU, P, Q = pluq(A)
        >>> LU###line 1790:_sage_    >>> LU
        [1 0 1 0]
        [1 1 0 0]
        [0 0 1 0]
        [1 1 1 0]

        >>> P###line 1796:_sage_    >>> P
        [0, 1, 2, 3]

        >>> Q###line 1799:_sage_    >>> Q
        [1, 2, 3, 3]
    """


if __name__ ==  '__main__':
    verbose = False
    do_timeit = False
    output_filename = '/home/malb/SAGE/devel/sage/sage/matrix/matrix_mod2_dense.pyx.timeit.sobj'

    import sys
    sys.path = sys.path + ['/usr/local/sage-3.4.1/local/bin']
    import sagedoctest

    # execfile('/home/malb/SAGE/devel/sage/sage/matrix/matrix_mod2_dense.pyx')
    m = sys.modules[__name__]
    m.__file__ = '/home/malb/SAGE/devel/sage/sage/matrix/matrix_mod2_dense.pyx'

    # configure special sage doc test runner
    runner = sagedoctest.SageDocTestRunner(checker=None, verbose=verbose, optionflags=0)
    runner._collect_timeit_stats = do_timeit
    runner._reset_random_seed = True

    runner = sagedoctest.testmod_returning_runner(m,
                   # filename='/home/malb/SAGE/devel/sage/sage/matrix/matrix_mod2_dense.pyx',
                   verbose=verbose,
                   globs=globals(),
                   runner=runner)
    runner.save_timeit_stats_to_file_named(output_filename)
    quit_sage(verbose=False)
    sys.exit(runner.failures)
