r"""
A = Z/10 x Z/2
G = (Z/10)^x \simeq Z/2 x Z/2
  = {1,-1,3,-3}
"""

import itertools
import operator

from hecke_hurwitz.invariant_partitions import invariant_partitions

def our_action():
    # Z/10 x Z/2: (a,b) -> 1 + a + 10b
    p_one = [(1 + (-1*a)%10 + 10*b) for b in range(2) for a in range(10)]
    p_three = [(1 + (3*a)%10 + 10*b) for b in range(2) for a in range(10)]
    S = SymmetricGroup(20)
    return PermutationGroup([S(p_one), S(p_three)])
#    return libgap.Group([libgap.PermList(p_one), libgap.PermList(p_three)])

def useful_tuples(partition):
    r"""
    A G-partition is useful if
    - 
    """
    # Z/10 x Z/2: (a,b) -> 1 + a + 10b

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

        sage: next(useful_partitions())
        ([[(0, 0)],
          [(2, 0), (8, 0), (1, 1), (9, 1)],
          [(2, 1), (8, 1)],
          [(4, 0), (6, 0), (3, 1), (7, 1)],
          [(4, 1), (6, 1)]],
         [[(1, 0)], [(3, 0)], [(5, 0), (0, 1), (5, 1)], [(7, 0)], [(9, 0)]])
    """
    for partition in invariant_partitions(our_action(), range(1, 21)):
        for i1, i2, i3, i4 in useful_tuples(partition):
            ppartition = [[((x - 1) % 10, (x - 1)//10) for x in atom] for atom in partition]
            atom0 = ppartition[0]
            atom1 = ppartition[i1]
            atom2 = ppartition[i2]
            atom3 = ppartition[i3]
            atom4 = ppartition[i4]
            for i in sorted([0, i1, i2, i3, i4], reverse=True):
                del ppartition[i]
            yield [atom0, atom1, atom2, atom3, atom4], ppartition

def str_atom(atom):
    return '[' + ', '.join('({},{})'.format(a,b) for a,b in atom) + ']'

def str_partition(partition):
    return '[' + ', '.join(str_atom(atom) for atom in partition) + ']'

def write_useful_partitions(filename):
    with open(filename, 'w') as output:
        for i, (M, R) in enumerate(useful_partitions()):
            output.write('useful partition {}\n'.format(i))
            output.write('M0={}\nM1={}\nM2={}\nM3={}\nM4={}\nR={}\n'.format(str_atom(M[0]), str_atom(M[1]), str_atom(M[2]), str_atom(M[3]), str_atom(M[4]), str_partition(R)))
            output.write('-' * 80)
            output.write('\n')
