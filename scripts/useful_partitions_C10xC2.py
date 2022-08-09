r"""
Exploration of the candidate group A = Z/10 x Z/2
"""

import itertools
import operator

from hecke_hurwitz.invariant_partitions import invariant_partitions
from hecke_hurwitz.abelian_groups import *

# TODO: make this a generic methods
def C10xC2_action():
    r"""
    Return the action of A as a permutation group on {1, 2, ..., 20}.

    The bijection used is `Z/10 \times Z/2: (a,b) \mapsto 1 + a + 10b`.

    EXAMPLES::

        sage: C10xC2_action()
        Permutation Group with generators [(2,4,10,8)(3,7,9,5)(12,14,20,18)(13,17,19,15), (2,10)(3,9)(4,8)(5,7)(12,20)(13,19)(14,18)(15,17)]
    """
    p_one = abelian_group_mul_permutation((10, 2), (-1, 1))
    p_three = abelian_group_mul_permutation((10, 2), (3, 1))
    return PermutationGroup([p_one, p_three])

# TODO: implement generic maps
#  group to int -> done in abelian_groups.element_encoding
#  int to group -> done in abelian_groups.element_decoding

def hermitian_inner_product(x, y):
    assert x.parent() is y.parent()
    return sum(xx.conjugate() * yy for xx, yy in zip(x, y))

# TODO: is that really the action we want to consider??
def perm_vector_action(perm, v):
    return v.parent()([v[perm(i+1)-1] for i in range(20)])

def useful_tuples(partition):
    if len(partition) < 5:
        return

    # condition 1: M0 contains only the identity
    if partition[0] != [1]:
        return

    # find atoms that projects onto {i} and are -1 invariant for i=1,2,3,4
    m1_candidates = []
    m2_candidates = []
    m3_candidates = []
    m4_candidates = []
    minus_one = [None] * 21
    three = [None] * 21
    for a in range(10):
        for b in range(2):
            x = 1 + a + 10*b
            y = 1 + (- a)%10 + 10 * b
            minus_one[x] = y

            y = 1 + (3 * a)%10 + 10 * b
            three[x] = y

    for i, atom in enumerate(partition):
        if any(minus_one[x] not in atom for x in atom):
            continue
        if (1 + 1) in atom or (1 + 1 + 10) in atom:
            if len(atom) >= 3:
                m1_candidates.append(i)
        if (1 + 2) in atom or (1 + 2 + 10) in atom:
            m2_candidates.append(i)
        if (1 + 3) in atom or (1 + 3 + 10) in atom:
            if len(atom) >= 3:
                m3_candidates.append(i)
        if (1 + 4) in atom or (1 + 4 + 10) in atom:
            m4_candidates.append(i)
    for i1, i2, i3, i4 in itertools.product(m1_candidates, m2_candidates, m3_candidates, m4_candidates):
        if len(set([i1,i2,i3,i4])) != 4:
            continue
        atom1 = partition[i1]
        atom2 = partition[i2]
        atom3 = partition[i3]
        atom4 = partition[i4]
        if len(atom1) == len(atom3) and all(three[x] in atom3 for x in atom1):
            assert all(three[x] in atom1 for x in atom3)
            yield (i1, i2, i3, i4)

def useful_partitions():
    r"""
    Run through all useful partitions of A

    A = bigcup M_i

    (1) M_0 = {identity} = {(0,0)}
    (2) for i={1,2,3,4} projection of M_ to the Z/10 factor includes i and -i
    (3) for i={1,2,3,4} M_i is invariant under multiplication by -1 on the first factor
    (4) M_1 and M_3 are exchaned by multiplication by 3 (on the first factor) and both
        contain at least three elements

    EXAMPLES::

        sage: partition, i1, i2, i3, i4 = next(useful_partitions(as_int=True))
        sage: partition
        [[1], [2], [3, 9, 12, 20], [4], [5, 7, 14, 18], [6, 11, 16], [8], [10], [13, 19], [15, 17]]
        sage: print(i1, i2, i3, i4)
        2, 8, 4, 9
    """
    for partition in invariant_partitions(C10xC2_action(), range(1, 21)):
        for i1, i2, i3, i4 in useful_tuples(partition):
            yield partition, i1, i2, i3, i4

def str_atom(atom):
    return '[' + ', '.join('({},{})'.format(a,b) for a,b in atom) + ']'

def str_partition(partition):
    return '[' + ', '.join(str_atom(atom) for atom in partition) + ']'

def str_useful_partition(partition, i1, i2, i3, i4, sep='\n'):
    r"""
    EXAMPLES::

        sage: up = [[[1], [2], [3, 9, 12, 20], [4], [5, 7, 14, 18], [6, 11, 16], [8], [10], [13, 19], [15, 17]], 2, 8, 4, 9]
        sage: str_useful_partition(*up)
        'M0=[(0,0)]\nM1=[(2,0), (8,0), (1,1), (9,1)]\nM2=[(2,1), (8,1)]\nM3=[(4,0), (6,0), (3,1), (7,1)]\nM4=[(4,1), (6,1)]\nR=[[(1,0)], [(3,0)], [(5,0), (0,1), (5,1)], [(7,0)], [(9,0)]]\n'
    """
    ppartition = [[((x - 1) % 10, (x - 1)//10) for x in atom] for atom in partition]
    atom0 = str_atom(ppartition[0])
    atom1 = str_atom(ppartition[i1])
    atom2 = str_atom(ppartition[i2])
    atom3 = str_atom(ppartition[i3])
    atom4 = str_atom(ppartition[i4])
    for i in sorted([0, i1, i2, i3, i4], reverse=True):
        del ppartition[i]
    ppartition = str_partition(ppartition)
    return f'M0={atom0}{sep}M1={atom1}{sep}M2={atom2}{sep}M3={atom3}{sep}M4={atom4}{sep}R={ppartition}'

def write_useful_partitions(filename):
    with open(filename, 'w') as output:
        for i, up in enumerate(useful_partitions()):
            output.write('useful partition {}\n'.format(i))
            output.write(str_useful_partition(*up))
            output.write('\n')
            output.write('-' * 80)
            output.write('\n')

def perm_finest_partition(perm):
    r"""
    Return the finest partition preserved by the permutation ``perm``.

    EXAMPLES::

        sage: S = SymmetricGroup(20)
        sage: perm_finest_partition(S('(1,2)'))
        [[1], [2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20], [11]]
    """
    char_list = characters()
    char_list_conjugate = [ch.conjugate() for ch in char_list]
    atoms = []
    for i, v in enumerate(char_list):
        vv = [perm_vector_action(perm, v).dot_product(w) for w in char_list_conjugate]
        atom = [i for i in range(20) if not vv[i].is_zero()]
        assert i in atom
        if i == atom[0]:
            atoms.append(atom)
    return [[i + 1 for i in atom] for atom in atoms]


def cc_stabilizer(cc, partition):
    r"""
    Iterate through elements of the conjugacy class ``cc`` that stabilizes the
    splitting of `\CC[A]` associated to ``partition`` 

    EXAMPLES::

        sage: S = SymmetricGroup(20)
        sage: perm_finest_partition(S('(1,2)'))

    """
    if sum(cc) != 20:
        cc = list(cc) + [1]*(20-sum(cc))
    char_list = characters()
    char_list_conjugate = [v.conjugate() for v in char_list]
    n = sum(len(atom) for atom in partition)
    S = SymmetricGroup(n)
    for perm in S.conjugacy_class(cc):
        stab = True
        for atom in partition:
            for i in atom:
                vv = perm_vector_action(perm, char_list_conjugate[i-1])
                s = sum(char_list[j].dot_product(vv) for j in range(20)) / 20
                assert s.abs().is_one(), s
                dp = sum(char_list[j-1].dot_product(vv) for j in atom) / 20
                if not dp.abs().is_one():
                    stab = False
                    break
            if not stab:
                break
        if stab:
            yield perm
