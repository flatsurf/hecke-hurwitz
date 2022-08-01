r"""
Module to compute invariant partitions of permutation groups.

Given a subgroup G of the symmetric group on `\{1, 2, \ldots, n\}`
we want to compute the list of set partitions of `\{1,2,\ldots,n\}`
that are stabilized by ``G``. In the case the action is transitive,
these are the so called "blocks". In non-transitive situation, one
decompose the action into a list of transitive actions, compute
the blocks of each of them and merge them when possible.
"""

import itertools
import operator

from sage.misc.flatten import flatten
from sage.libs.gap.libgap import libgap

from sage.combinat.set_partition import SetPartitions

def group_by(elements, test):
    r"""
    Return a triple ``(representatives, equivalence_classes, group_ids)``
    """
    representatives = []
    equivalence_classes = []
    group_ids = []
    for i, x in enumerate(elements):
        equivalent = False
        for j, y in enumerate(representatives):
            if test(x, y):
                equivalent = True
                break
        if equivalent:
            group_ids.append(j)
            equivalence_classes[j].append(i)
        else:
            group_ids.append(len(equivalence_classes))
            equivalence_classes.append([i])
            representatives.append(x)
    return representatives, equivalence_classes, group_ids

def invariant_partitions(G, domain=None):
    r"""
    Iterate through ``G``-invariant partitions of the domain.

    INPUT:

    - ``G`` - a Sage permtation group

    - ``domain`` - (optional) the domain on which ``G`` acts

    EXAMPLES::

        sage: from hecke_hurwitz.invariant_partitions import invariant_partitions
        sage: G = PermutationGroup(['(1,2,3)','(4,5)'])
        sage: list(invariant_partitions(G))
        [[[1], [2], [3], [4], [5]],
         [[1], [2], [3], [4, 5]],
         [[1, 2, 3], [4], [5]],
         [[1, 2, 3, 4, 5]],
         [[1, 2, 3], [4, 5]]]
    """
    G = libgap(G)

    assert G.IsCommutative()
    if domain:
        orbits = libgap.Orbits(G, list(domain))
    else:
        orbits = libgap.Orbits(G)
    ground_set = set().union(*orbits)

    # group orbits by isomorphism types to avoid perform identical block computations
    # NOTE: if non-Abelian, one need to switch to libgap.IsConjugate instead of operator.eq
    orbit_stabilizers = [libgap.Stabilizer(G, orb[0]) for orb in orbits]
    stabilizers_cc, values, orbit_to_cc = group_by(orbit_stabilizers, operator.eq)
    assert len(orbit_to_cc) == len(orbits)
    actions = [libgap.ActionHomomorphism(G, libgap.RightCosets(G, H), libgap.OnRight) for H in stabilizers_cc]
    perm_groups = [libgap.Image(action) for action in actions]

    # recompute orbit so that it matches the canonical ordering of the action
    representatives_cc = [[libgap.RepresentativeAction(g, 1, i) for i in range(1, libgap.Index(G, stab) + 1)] for g, stab in zip(perm_groups, stabilizers_cc)]
    representatives_cc_preimages = [[libgap.PreImagesRepresentative(action, p) for p in rep] for action, rep in zip(actions, representatives_cc)]
    orbit_lifts = []
    for i, orbit in enumerate(orbits):
        j = orbit_to_cc[i]
        assert len(orbit) == len(representatives_cc[j]), (len(orbit), len(representatives_cc[j]))
        orbit_lifts.append([libgap.OnPoints(orbit[0], r) for r in representatives_cc_preimages[j]])
    orbits = orbit_lifts
    del orbit_lifts

    # for each orbit isomorphism class compute blocks, stabilizers and orbits
    blocks_cc = []
    for g, stab in zip(perm_groups, stabilizers_cc):
        n = libgap.Index(G, stab)
        if n == 1:
            blocks_cc.append(libgap.Concatenation([[1]], libgap.AllBlocks(g)))
        else:
            blocks_cc.append(libgap.Concatenation([[1]], libgap.AllBlocks(g), [list(range(1,n+1))]))
    block_cc_orbits = [[libgap.Orbit(g, block, libgap.OnSets) for block in b] for g, b in zip(perm_groups, blocks_cc)]
    block_cc_stabilizers = [[libgap.Stabilizer(g, block, libgap.OnSets) for block in b] for g, b in zip(perm_groups, blocks_cc)]

    # go through all blocks of each orbit
    for prod_id in itertools.product(*[range(len(blocks_cc[j])) for j in orbit_to_cc]):
        # identify stablizers equality for possible fusion
        # NOTE: if non-Abelian this is wrong
        prod_block_orbits = [[sorted(orbit[x-1] for x in b) for b in block_cc_orbits[j][k]] for orbit,j,k in zip(orbits, orbit_to_cc, prod_id)]
        assert set(flatten(prod_block_orbits)) == ground_set
        block_stabilizers = [libgap.Stabilizer(G, orbit[0], libgap.OnSets) for orbit in prod_block_orbits]
        #libgap.PreImage(actions, stab
        representatives, equivalence_classes, group_id = group_by(block_stabilizers, operator.eq)
        assert sum(len(cl) for cl in equivalence_classes) == len(prod_block_orbits) == len(block_stabilizers)

        # choose which orbits to merge
        for xx in itertools.product(*[SetPartitions(cl) for cl in equivalence_classes]):
            xx = tuple(tuple(y) for atom in xx for y in atom)
            u = [[prod_block_orbits[y] for y in atom] for atom in xx]
            assert set().union(*xx) == set(range(sum(len(x) for x in xx)))
            # choose which element in each orbit to merge
            # for each element of the set partition pick an element of the orbit and merge according to the action
            # note that we fix the first
            singletons = [pos for pos in range(len(xx)) if len(xx[pos]) == 1]
            non_singletons = [pos for pos in range(len(xx)) if len(xx[pos]) != 1]
            for yy in itertools.product(*[itertools.product(*[prod_block_orbits[y] for y in xx[pos][1:]]) for pos in non_singletons]):
                partition = [None] * len(xx)
                for pos in singletons:
                    # no merge
                    x = xx[pos]
                    partition[pos] = prod_block_orbits[x[0]][0]
                for pos, y in zip(non_singletons, yy):
                    x = xx[pos]
                    partition[pos] = tuple(sorted(sum((tuple(h) for h in y), tuple(prod_block_orbits[x[0]][0]))))
                # merge xx[i][0] with yy[i]
                partition = libgap.Concatenation([libgap.Orbit(G, p, libgap.OnSets) for p in partition])
                assert set().union(*partition) == ground_set
                partition = list((list(map(int,x)) for x in partition))
                partition.sort()
                yield partition


