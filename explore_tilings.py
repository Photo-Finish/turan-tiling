"""
Turan & Tiling — Computational Exploration
============================================
Explores constraints on polygons / polygon combinations that can tile the plane,
using the contact-graph perspective.

Key question: What does the contact graph of a plane tiling tell us about
the polygons involved?

Author: Generated exploration
Date: 2026-08-02
"""

import itertools
import math
from collections import defaultdict, Counter
from fractions import Fraction

# ============================================================
# PART 1: Fundamental constraint from Euler + contact graph
# ============================================================

def average_face_degree(k):
    """
    For a monohedral edge-to-edge tiling by congruent k-gons:
    On a torus (periodic): V - E + F = 0
    V = n (tiles), E = kn/2, F = n(k/2 - 1)
    Average face degree = 2E/F = 2k/(k-2)
    """
    if k <= 2:
        return None
    return 2 * k / (k - 2)


def check_monohedral_feasibility():
    """Check which k allow monohedral edge-to-edge tiling."""
    print("=" * 60)
    print("PART 1: Monohedral edge-to-edge tiling constraints")
    print("=" * 60)
    print()
    print("For k-gon monohedral tiling (interior, edge-to-edge):")
    print("  V = n tiles, E = kn/2, F = n(k/2-1)")
    print("  avg face degree m_avg = 2k/(k-2)")
    print("  Requirement: m_avg >= 3 (at least 3 tiles meet at each vertex)")
    print()

    for k in range(3, 13):
        m_avg = average_face_degree(k)
        feasible = m_avg >= 3 if m_avg else False
        marker = "[OK] FEASIBLE" if feasible else "[XX] IMPOSSIBLE"
        print(f"  k={k:2d}: m_avg = {m_avg:8.4f}  ->  {marker}")

    print()
    print("CONCLUSION: Only k in {3,4,5,6} can monohedrally tile the plane")
    print("            in an edge-to-edge manner. (Matches known classification!)")
    print()


# ============================================================
# PART 2: Vertex types for regular polygon tilings
# ============================================================

def enumerate_vertex_types(max_k=12, max_m=6):
    """
    Enumerate all solutions to:
    Σ(1/k_i) = (m-2)/2
    with k_i >= 3, m = number of polygons meeting at a vertex.
    
    These are the Archimedean vertex types for regular polygon tilings.
    """
    print("=" * 60)
    print("PART 2: Vertex types for regular polygon tilings")
    print("=" * 60)
    print()
    print("Equation: Σ 1/k_i = (m-2)/2  (angles of regular k_i-gons sum to 360°)")
    print()

    all_solutions = defaultdict(list)

    for m in range(3, max_m + 1):
        target = Fraction(m - 2, 2)
        print(f"--- m = {m} (target Σ 1/k_i = {target}) ---")

        # Generate k_i values (3 to max_k)
        # Use recursion with pruning
        solutions = []

        def search(remaining_sum, current, start_k):
            if len(current) == m:
                if remaining_sum == 0:
                    solutions.append(tuple(sorted(current)))
                return
            # With all remaining slots filled with min value (3),
            # the remaining sum must be achievable
            remaining_slots = m - len(current)
            min_possible = Fraction(remaining_slots, max_k)
            max_possible = Fraction(remaining_slots, 3)
            if remaining_sum < min_possible or remaining_sum > max_possible:
                return

            for k_i in range(start_k, max_k + 1):
                term = Fraction(1, k_i)
                if term > remaining_sum:
                    continue  # skip, try larger k (which gives smaller 1/k)
                search(remaining_sum - term, current + [k_i], k_i)

        search(target, [], 3)

        # Deduplicate
        unique = sorted(set(solutions))
        all_solutions[m] = unique

        for sol in unique:
            # Verify
            total = sum(Fraction(1, k) for k in sol)
            angle_sum = sum((k - 2) * 180 / k for k in sol)
            sol_str = "[" + ",".join(str(x) for x in sol) + "]"
            print(f"    {sol_str:22s}  Σ1/k={float(total):.6f}  Σangles={angle_sum:.1f}°")

        if not unique:
            print("    (none)")

        print()

    return all_solutions


# ============================================================
# PART 3: Monohedral tiling analysis via contact graph
# ============================================================

def analyze_monohedral_contact_graph(k, n=1000):
    """
    Analyze the contact graph of a monohedral tiling by k-gons.
    """
    print("=" * 60)
    print(f"PART 3a: Contact graph analysis for k={k}-gon tiling")
    print("=" * 60)
    print()

    E = k * n / 2
    F = E - n  # for torus: V-E+F=0, V=n
    m_avg = 2 * E / F if F > 0 else float('inf')

    print(f"  Number of tiles (V):      {n}")
    print(f"  Edges in contact graph:   {E:.0f} (= kn/2)")
    print(f"  Faces (tiling vertices):  {F:.0f} (= n(k/2-1))")
    print(f"  Avg tiles per vertex:     {m_avg:.4f}")
    print()

    # Turán bound comparison
    turan_bound = 3 * n * n / 8  # for K_5-free, r=4
    planar_bound = 3 * n - 6
    print(f"  --- Bounds on edge count ---")
    print(f"  Actual edges:           {E:.0f}")
    print(f"  Turán bound (K₅-free):  {turan_bound:.0f}  (≈3n^2/8)")
    print(f"  Planar bound:           {planar_bound:.0f}  (3n-6)")
    print(f"  Turán/actual ratio:     {turan_bound/E:.2f}")
    print(f"  Planar/actual ratio:    {planar_bound/E:.2f}")
    print()

    # Chromatic number considerations
    # For a bipartite contact graph (like square tiling), chi=2
    # Turán: K₃-free -> bipartite -> chi=2
    # For K₅-free and not K₃-free, chi can be 3 or 4
    print(f"  --- Chromatic constraints ---")
    print(f"  k={k}: contact graph is {k}-regular planar")
    if k == 4:
        print(f"  Square tiling -> contact graph is bipartite (chi=2, K₃-free)")
    elif k == 3:
        print(f"  Triangular tiling -> contact graph is tripartite?")
        print(f"  Regular triangle tiling: each triangle touches 3 others")
        print(f"  The contact graph of triangular tiling is the hexagonal")
        print(f"  lattice (3-regular), which is bipartite! (chi=2)")
    elif k == 6:
        print(f"  Hexagonal tiling -> contact graph is triangular lattice (6-regular)")
        print(f"  The contact graph is the dual of the hexagonal lattice,")
        print(f"  which is the triangular lattice. It's 3-colorable but not bipartite.")
    print()

    return {
        'k': k, 'n': n, 'E': E, 'F': F, 'm_avg': m_avg,
        'turan_bound': turan_bound, 'planar_bound': planar_bound
    }


def analyze_all_monohedral():
    """Run contact graph analysis for all feasible k."""
    results = []
    for k in [3, 4, 5, 6]:
        r = analyze_monohedral_contact_graph(k)
        results.append(r)
    return results


# ============================================================
# PART 4: Mixed polygon tilings (Archimedean)
# ============================================================

def analyze_archimedean_tilings():
    """
    Analyze the 11 Archimedean (uniform) tilings.
    Each has a vertex type (a.b.c...) and known proportions of tiles.
    """
    print("=" * 60)
    print("PART 4: Archimedean (uniform) tilings")
    print("=" * 60)
    print()

    # The 11 Archimedean tilings (vertex type, name)
    archimedean = [
        ((3,3,3,3,3,3), "3^6 - Regular triangular"),
        ((4,4,4,4), "4^4 - Square grid"),
        ((6,6,6), "6^3 - Regular hexagonal"),
        ((3,3,3,4,4), "3^3.4^2 - Snub square"),
        ((3,3,4,3,4), "3^2.4.3.4 - Snub square (variant)"),
        ((3,3,3,3,6), "3^4.6 - Snub hexagonal"),
        ((3,6,3,6), "3.6.3.6 - Trihexagonal"),
        ((3,3,6,6), "3^2.6^2 - Truncated hexagonal"),
        ((3,4,6,4), "3.4.6.4 - Rhombitrihexagonal"),
        ((4,8,8), "4.8^2 - Truncated square"),
        ((3,12,12), "3.12^2 - Truncated hexagonal"),
    ]

    for vertex_type, name in archimedean:
        m = len(vertex_type)
        # Count tile types
        counter = Counter(vertex_type)
        total_tiles_per_vertex = sum(1/k for k in vertex_type)  # each tile counted at 1/k of its vertices

        # Each k-gon appears at each of its k vertices
        # So proportion of each type ∝ counter[k]/k
        proportions = {}
        total_prop = sum(counter[k] / k for k in counter)
        for k in counter:
            proportions[k] = (counter[k] / k) / total_prop

        avg_sides = sum(k * proportions[k] for k in proportions)

        print(f"  {name}")
        print(f"    Vertex type: {vertex_type}")
        print(f"    m = {m} tiles per vertex")
        print(f"    Tile proportions: {', '.join(f'{k}-gon: {proportions[k]*100:.1f}%' for k in sorted(proportions))}")
        print(f"    Average # sides: {avg_sides:.3f}")

        # Contact graph parameters
        # avg degree = weighted avg of k values
        # But this is subtle - each k-gon has k contacts
        avg_degree = avg_sides
        print(f"    Avg degree in contact graph: {avg_degree:.3f}")
        m_from_formula = 2 * avg_degree / (avg_degree - 2) if avg_degree > 2 else float('inf')
        print(f"    Predicted avg tiles/vertex: {m_from_formula:.3f} (actual: {m})")
        print()

    print("NOTE: All these satisfy Σ 1/k_i = (m-2)/2.")
    print("      The contact graph approach correctly predicts the relationship")
    print("      between average polygon sides and average vertex valence.")
    print()


# ============================================================
# PART 5: Turán-type constraints on tiling contact graphs
# ============================================================

def turan_analysis():
    """
    Explore what Turán's theorem says about tiling contact graphs.
    """
    print("=" * 60)
    print("PART 5: Turán-type analysis of contact graphs")
    print("=" * 60)
    print()

    print("""
Key insight: The contact graph G of any plane tiling is:
  1. Planar -> contains no K₅ or K₃,₃ subdivision
  2. Edge-to-edge with k-gons -> each vertex degree <= k

Turán's theorem (K_{r+1}-free) gives an UPPER bound on edges.
But this bound is O(n^2), while planar graphs only allow O(n) edges.
So the Turán bound is too loose for direct application.

However, we can ask: what if we only know that G is K_{r+1}-free,
NOT that it's planar? What does Turán tell us?

For a tiling by k-gons:
  - G has n vertices, degree approx k (interior)
  - E = kn/2
  - For G to not exceed Turán bound t(n,r):
      kn/2 <= t(n,r) ≈ (1 - 1/r)n^2/2
      k <= (1 - 1/r)n

  This is automatically satisfied for any fixed k as n grows.
  
But consider: if the contact graph has a K_{r+1} subgraph, that means
r+1 tiles are ALL pairwise adjacent (share boundaries).
This is geometrically impossible for r >= 4 in the plane (planarity).
So G MUST be K₅-free -> r=4 -> Turán says E <= 3n^2/8.

More interesting: what about K₃-free (triangle-free) contact graphs?
  - Square tilings: contact graph IS bipartite = K₃-free
  - For a K₃-free contact graph, Turán (r=2) says E <= n^2/4
  - But for k-gon tilings, E = kn/2
  - So kn/2 <= n^2/4 -> k <= n/2  (trivially true)

The REAL power: if we classify tiles into "types" and forbid
certain adjacency patterns between types, Turán gives density bounds.
""")

    print("--- Turán bounds for different forbidden clique sizes ---")
    print()
    print(f"{'r':>3} {'Forbids':>8} {'t(n,r)/n^2':>12} {'Max k for k-gon tiling':>25}")
    print("-" * 55)
    for r in range(1, 7):
        density = (1 - 1/r) / 2
        # For a k-gon tiling: E = kn/2, Turán bound: (1-1/r)n^2/2
        # So k <= (1-1/r)n. For large n, this is trivial.
        # But in terms of density: k/(2n) <= (1-1/r)/2n... hmm
        print(f"  {r:2d}  {'K_' + str(r+1):>7s}  {density:.6f}    {'no constraint (trivial for fixed k)'}")

    print()
    print("CONCLUSION: For MONOHEDRAL tilings, Turán gives no meaningful constraint")
    print("beyond planarity. But for analyzing adjacency PATTERNS between different")
    print("tile types, Turán can constrain which type-to-type adjacencies are possible.")
    print()


# ============================================================
# PART 6: Tile-type adjacency and Turán
# ============================================================

def type_adjacency_analysis():
    """
    Suppose tiles are colored by type. Consider the adjacency graph restricted
    to edges between same-type tiles vs different-type tiles.
    """
    print("=" * 60)
    print("PART 6: Type-to-type adjacency constraints")
    print("=" * 60)
    print()

    print("""
Consider tiling with multiple polygon types (e.g., squares + triangles).
Split the contact graph into r "parts" by tile type.

Observation: In many tilings, SAME-type tiles don't touch each other
(or touch rarely). For instance:

  - Square grid: all tiles same type, all touch (K₃-free though)
  - Trihexagonal (3.6.3.6): triangles only touch hexagons, hexagons only
    touch triangles -> the contact graph is BIPARTITE!
  - Snub square (3^3.4^2): each triangle touches 2 other triangles + 1 square

If tiles of type i NEVER touch tiles of type i, then the contact graph
is r-partite with EMPTY parts -> it's a complete r-partite graph = Turán graph!

So: the contact graph being exactly T(n,r) means each tile type forms
an independent set, and every tile touches every tile of different types.

Can this happen in a planar tiling?
  - T(n,2) = complete bipartite K_{a,b}: planar only if a<=2 or b<=2
    (contains K_{3,3} otherwise)
  - T(n,3): contains K_{3,3,3}... which contains K_{3,3} -> not planar for n>=9
  - T(n,4): contains K₅ -> not planar

So TRUE Turán graphs CANNOT be realized as planar contact graphs
(except very small cases). This is why planar tilings have
constrained contact patterns!

What CAN be realized? If same-type adjacency is FORBIDDEN:
  - Need an r-partite PLANAR graph
  - Max edges in r-partite planar graph << Turán bound
  - This gives a "planar Turán" type problem!
""")

    print("--- Planar r-partite graphs: max edges ---")
    print()
    print("For an r-partite planar graph on n vertices:")
    print("  r=2 (bipartite planar): E <= 2n-4 (n>=3)")
    print("  r=3 (tripartite planar): E <= ?")
    print("  r=4: E <= 3n-6 (same as planar, since K₅ avoided anyways)")
    print()
    print("This highlights a RICH area: planar Turán numbers!")
    print("Given a forbidden graph H, what's max edges in a PLANAR graph")
    print("that avoids H?")
    print()


# ============================================================
# PART 7: Generate possible vertex configurations
# ============================================================

def generate_angle_partitions(k, target_degrees=360):
    """
    For a k-gon, enumerate ways to assign integer angles such that
    sum = (k-2)*180 and each angle divides the 360° at a tiling vertex.
    """
    print("=" * 60)
    print(f"PART 7: Angle partitions for {k}-gons in edge-to-edge tilings")
    print("=" * 60)
    print()

    angle_sum = (k - 2) * 180
    print(f"  Total interior angle sum: {angle_sum}°")
    print(f"  At each tiling vertex, angles meeting must sum to 360°")
    print(f"  Each angle must be <= 180° (convex) and > 0°")
    print()

    # For monohedral tilings, each angle must divide 360° so that
    # several copies of the SAME polygon can meet at that vertex.

    # Find divisors of 360 between 1 and 180
    divisors = [d for d in range(1, 181) if 360 % d == 0]
    # Also consider angles that can combine to 360 in non-regular tilings
    # For regular polygons: angle = (k-2)*180/k
    regular_angle = (k - 2) * 180 / k
    n_at_vertex = 360 / regular_angle

    print(f"  Regular {k}-gon: each angle = {regular_angle:.1f}°")
    print(f"  Number meeting at a vertex: {n_at_vertex:.3f}")
    print(f"  m_avg = {2*k/(k-2):.4f} (from contact graph formula)")
    print()

    # Enumerate ALL possible k-tuples of divisors that sum to angle_sum
    # This is huge for large k, so only do small k or limited enumeration
    if k <= 6:
        print(f"  Enumerating angle assignments (each angle divides 360°)...")
        solutions_by_multiplicity = defaultdict(list)

        def enum_angles(remaining, current):
            if len(current) == k:
                if remaining == 0:
                    solutions_by_multiplicity[tuple(sorted(current))].append(current)
                return
            if remaining <= 0:
                return
            start = current[-1] if current else 1
            for a in divisors:
                if a < start:
                    continue
                if a > remaining:
                    break
                # Check: 360 % a == 0 (for edge-to-edge, a must divide 360)
                if remaining - a < (k - len(current) - 1) * a:
                    break  # can't fill remaining with values >= a
                enum_angles(remaining - a, current + [a])

        enum_angles(angle_sum, [])

        print(f"  Found {len(solutions_by_multiplicity)} distinct (sorted) angle tuples")
        # Show a few examples
        for i, (angles, _) in enumerate(sorted(solutions_by_multiplicity.items())):
            if i >= 10:
                print(f"  ... and {len(solutions_by_multiplicity)-10} more")
                break
            mults = defaultdict(int)
            for a in angles:
                mults[a] += 1
            mult_str = ", ".join(f"{m}×{a}°" for a, m in sorted(mults.items()))
            print(f"    {mult_str}")
        print()

    return


# ============================================================
# PART 8: Summary - what the contact graph approach gives us
# ============================================================

def summary():
    print("=" * 60)
    print("SUMMARY: Contact Graph -> Tiling Constraints")
    print("=" * 60)
    print()
    print("""
-------------------------------------------------------------
| 1. MONOHEDRAL EDGE-TO-EDGE TILING                           |
|    k-gons -> k-regular planar contact graph                  |
|    Euler on torus: 2k/(k-2) = avg tiles per vertex >= 3     |
|    => k <= 6                                                  |
|    [OK] Triangle, quadrilateral, pentagon, hexagon possible    |
|    [XX] Heptagon+ impossible (edge-to-edge, monohedral)        |
-------------------------------------------------------------
| 2. REGULAR POLYGON VERTEX TYPES                             |
|    Σ 1/k_i = (m-2)/2                                        |
|    Enumerates all Archimedean vertex types                  |
|    Same equation as angle sum = 360°                        |
-------------------------------------------------------------
| 3. TURÁN CONNECTION                                         |
|    Pure Turán bound (O(n^2)) too loose vs planar (O(n))      |
|    BUT: r-partite classification of tiles + Turán           |
|    => Planar Turán numbers: max edges in planar H-free graph |
|    => If same-type adjacency forbidden, contact graph         |
|      is r-partite planar -> much tighter constraints         |
-------------------------------------------------------------
| 4. "重边" (MULTI-EDGES) REVISITED                            |
|    - Convex polygons: no true multi-edges                   |
|    - Weighted edges (contact length) -> isoperimetric const. |
|    - Non-edge-to-edge: same tile pair can touch on          |
|      multiple sides -> multi-edges in contact graph          |
|    - 3D polyhedra: true multi-faces possible                |
-------------------------------------------------------------
""")


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    check_monohedral_feasibility()
    vertex_types = enumerate_vertex_types(max_k=12, max_m=6)
    analyze_all_monohedral()
    analyze_archimedean_tilings()
    turan_analysis()
    type_adjacency_analysis()
    for k in [3, 4, 5, 6]:
        generate_angle_partitions(k)
    summary()



