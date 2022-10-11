r"""
Tools for manipulating Abelian groups
"""
from sage.misc.misc_c import prod
from sage.misc.cachefunc import cached_function

from sage.arith.functions import lcm

from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.modules.free_module import VectorSpace
from sage.rings.number_field.number_field import CyclotomicField


def abelian_group_element_encoding(invariants, g):
    r"""
    Return the integer in {1, ..., n} corresponding to the encoding of the group
    element ``g``.

    INPUT:

    - ``invariants`` -- the number defining the Abelian group

    - ``g`` -- the tuple defining the group element

    EXAMPLES::

        sage: from hecke_hurwitz.abelian_groups import abelian_group_element_encoding, abelian_group_element_decoding
        sage: abelian_group_element_encoding((10, 2), (0, 0))
        1
        sage: abelian_group_element_encoding((10, 2), (0, 1))
        11
        sage: [abelian_group_element_encoding((10,2), (a,b)) for b in range(2) for a in range(10)]
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
        sage: assert all(abelian_group_element_encoding((10,3,2), abelian_group_element_decoding((10,3,2), i)) == i for i in range(1, 61))
    """
    result = 1
    m = 1
    for a, k in zip(g, invariants):
        result += m * (int(a) % k)
        m *= k
    return result

def abelian_group_element_decoding(invariants, i):
    r"""
    EXAMPLES::

        sage: from hecke_hurwitz.abelian_groups import abelian_group_element_decoding
        sage: abelian_group_element_decoding((10, 2), 1)
        (0, 0)
        sage: abelian_group_element_decoding((10, 2), 11)
        (0, 1)
    """
    result = []
    i -= 1
    for k in invariants:
        result.append(i % k)
        i //= k
    return tuple(result)


def abelian_group_add(invariants, a, b):
    return tuple((aa+bb) % inv for (aa, bb, inv) in zip(a, b, invariants))


def abelian_group_mul(invariants, a, b):
    return tuple((aa*bb) % inv for (aa, bb, inv) in zip(a, b, invariants))


def abelian_group_mul_permutation(invariants, a):
    r"""
    Return the permutation ax + b in the Abelian group with the given
    invariants.

    EXAMPLES::

        sage: from hecke_hurwitz.abelian_groups import abelian_group_mul_permutation
        sage: abelian_group_mul_permutation((10, 2), (-1, 1))
        (2,10)(3,9)(4,8)(5,7)(12,20)(13,19)(14,18)(15,17)
    """
    n = prod(invariants)
    S = SymmetricGroup(n)
    result = []
    for i in range(1, n+1):
        x = abelian_group_element_decoding(invariants, i)
        y = abelian_group_mul(invariants, a, x)
        result.append(abelian_group_element_encoding(invariants, y))
    return S(result)


def abelian_group_add_permutation(invariants, a):
    n = prod(invariants)
    S = SymmetricGroup(n)
    result = []
    for i in range(1, n+1):
        x = abelian_group_element_decoding(invariants, i)
        y = abelian_group_add(invariants, a, x)
        result.append(abelian_group_element_encoding(invariants, y))
    return S(result)


def abelian_group_permutation_gens(invariants):
    gens = []
    k = len(invariants)
    for i in range(k):
        a = [0] * k
        a[i] = 1
        gens.append(abelian_group_add_permutation(invariants, a))
    return gens


def symmetric_group_gens(domain):
    r"""
    EXAMPLES::

        sage: symmetric_group_gens([1, 3])
        sage: symmetric_group_gens([1, 3, 5, 7])
    """
    if len(domain) == 1:
        return []
    transposition = '({}, {})'.format(domain[0], domain[1]) 
    if len(domain) == 2:
        return [transposition]
    cycle = '(' + ','.join(str(i) for i in domain) + ')'
    return [transposition, cycle]

@cached_function
def abelian_group_characters(invariants):
    r"""
    EXAMPLES::

        sage: from hecke_hurwitz.abelian_groups import abelian_group_characters
        sage: ch = matrix(abelian_group_characters((10,2)))
        sage: ch.rank()
        20
        sage: matrix([[s.hermitian_inner_product(t).abs() for s in ch.rows()] for t in ch.rows()]) == 20 * identity_matrix(20)
        True
    """
    n = prod(invariants)
    elements = [abelian_group_element_decoding(invariants, i) for i in range(1, n+1)]
    N = lcm(invariants)
    d = prod(invariants)
    K = CyclotomicField(N)
    V = VectorSpace(K, d)
    z = K.gen()
    zetas = [z**(N/i) for i in invariants]
    char_list = []
    for a in elements:
        ch = [prod(zeta**(aa*bb) for zeta, aa, bb in zip(zetas, a, b)) for b in elements]
        char_list.append(ch)
    return tuple(char_list)


