"""
Cluster Expansion Calculation Example

Demonstrates cluster expansion for configurational thermodynamics:
1. Basic cluster expansion setup
2. Structure generation
3. Training from DFT data
4. Formation energy predictions
5. Monte Carlo simulations
6. Ground state determination

Formula: E = J₀ + Σᵢ Jᵢ ⟨σᵢ⟩ + Σᵢⱼ Jᵢⱼ ⟨σᵢσⱼ⟩ + ...

Reference:
Sanchez, Ducastelle & Gratias, "Generalized cluster description of multicomponent systems,"
Physica A 128, 334 (1984)

Note: Requires icet package: pip install icet
"""

from ase.build import bulk
from ase.calculators.emt import EMT
import numpy as np

try:
    from icet import ClusterSpace, StructureContainer, ClusterExpansion
    from icet.tools import enumerate_structures

    ICET_AVAILABLE = True
except ImportError:
    ICET_AVAILABLE = False
    print("Warning: icet not installed. Install with: pip install icet")


def example_1_basic_cluster_space():
    """Example 1: Create and inspect cluster space."""
    if not ICET_AVAILABLE:
        print("Skipping: icet not available")
        return

    print("=" * 60)
    print("Example 1: Basic Cluster Space")
    print("=" * 60)

    # Create primitive structure
    prim = bulk("Cu", "fcc", a=3.6)

    print(f"\nPrimitive structure:")
    print(f"  Formula: {prim.get_chemical_formula()}")
    print(f"  Atoms: {len(prim)}")

    # Define cluster space for Cu-Au alloy
    # Cutoffs: [2-body, 3-body, 4-body] in Angstroms
    cutoffs = [6.0, 5.0, 4.0]

    cs = ClusterSpace(prim, cutoffs, chemical_symbols=[["Cu", "Au"]])

    print(f"\nCluster space:")
    print(f"  Chemical symbols: Cu, Au")
    print(f"  Cutoffs: {cutoffs} Å")
    print(f"  Number of orbits: {len(cs)}")

    # Inspect clusters
    print(f"\nCluster statistics:")
    n_pairs = sum(1 for orbit in cs if orbit.order == 2)
    n_triplets = sum(1 for orbit in cs if orbit.order == 3)
    n_quadruplets = sum(1 for orbit in cs if orbit.order == 4)

    print(f"  Zerolet (constant): 1")
    print(f"  Pairs: {n_pairs}")
    print(f"  Triplets: {n_triplets}")
    print(f"  Quadruplets: {n_quadruplets}")
    print(f"  Total parameters: {len(cs)}")

    # Print some orbit details
    print(f"\nFirst few orbits:")
    for i, orbit in enumerate(cs[:5]):
        print(
            f"  Orbit {i}: order={orbit.order}, "
            f"radius={orbit.radius:.3f} Å, "
            f"multiplicity={orbit.multiplicity}"
        )

    print()


def example_2_structure_enumeration():
    """Example 2: Generate training structures."""
    if not ICET_AVAILABLE:
        print("Skipping: icet not available")
        return

    print("=" * 60)
    print("Example 2: Structure Enumeration")
    print("=" * 60)

    # Primitive cell
    prim = bulk("Cu", "fcc", a=3.6)

    # Enumerate structures in supercells
    print(f"\nEnumerating structures...")
    print(f"  Base: Cu FCC")
    print(f"  Substituent: Au")
    print(f"  Supercell sizes: 1-4 atoms")

    structures = enumerate_structures(prim, sizes=range(1, 5), chemical_symbols=[["Cu", "Au"]])

    # Convert to list to count
    structure_list = list(structures)

    print(f"\nGenerated structures: {len(structure_list)}")

    # Group by composition
    compositions = {}
    for structure in structure_list:
        formula = structure.get_chemical_formula()
        compositions[formula] = compositions.get(formula, 0) + 1

    print(f"\nBy composition:")
    for formula in sorted(compositions.keys()):
        print(f"  {formula}: {compositions[formula]} structures")

    # Show a few examples
    print(f"\nExample structures:")
    for i, structure in enumerate(structure_list[:5]):
        n_Cu = structure.get_chemical_symbols().count("Cu")
        n_Au = structure.get_chemical_symbols().count("Au")
        x_Au = n_Au / (n_Cu + n_Au)

        print(
            f"  {i + 1}. {structure.get_chemical_formula()}: "
            f"x(Au) = {x_Au:.3f}, "
            f"atoms = {len(structure)}"
        )

    print()


def example_3_training_cluster_expansion():
    """Example 3: Train cluster expansion from synthetic data."""
    if not ICET_AVAILABLE:
        print("Skipping: icet not available")
        return

    print("=" * 60)
    print("Example 3: Training Cluster Expansion")
    print("=" * 60)

    # Setup
    prim = bulk("Cu", "fcc", a=3.6)
    cutoffs = [6.0, 5.0]
    cs = ClusterSpace(prim, cutoffs, chemical_symbols=[["Cu", "Au"]])

    print(f"\nCluster space: {len(cs)} parameters")

    # Generate training structures
    structures = list(
        enumerate_structures(prim, sizes=range(1, 6), chemical_symbols=[["Cu", "Au"]])
    )

    print(f"Training structures: {len(structures)}")

    # Calculate energies (using EMT as fast approximation)
    print(f"\nCalculating formation energies...")

    sc = StructureContainer(cs)

    # Reference energies (pure elements)
    cu_pure = bulk("Cu", "fcc", a=3.6)
    cu_pure.calc = EMT()
    E_Cu = cu_pure.get_potential_energy() / len(cu_pure)

    au_pure = bulk("Au", "fcc", a=4.08)
    au_pure.calc = EMT()
    E_Au = au_pure.get_potential_energy() / len(au_pure)

    print(f"  E(Cu) = {E_Cu:.4f} eV/atom")
    print(f"  E(Au) = {E_Au:.4f} eV/atom")

    # Calculate formation energies
    for structure in structures[:20]:  # Limit for speed
        structure.calc = EMT()
        E_total = structure.get_potential_energy()

        # Calculate formation energy
        n_Cu = structure.get_chemical_symbols().count("Cu")
        n_Au = structure.get_chemical_symbols().count("Au")
        E_ref = n_Cu * E_Cu + n_Au * E_Au
        E_form = (E_total - E_ref) / len(structure)

        # Add to container
        sc.add_structure(structure, properties={"energy": E_form})

    print(f"  Added {len(sc)} structures to training set")

    # Train cluster expansion
    print(f"\nTraining...")

    from trainstation import Optimizer

    opt = Optimizer(sc.get_fit_data())
    opt.train()

    # Get parameters
    parameters = opt.parameters

    print(f"  Training RMSE: {opt.rmse:.6f} eV/atom")
    print(f"  Number of parameters: {len(parameters)}")

    # Create cluster expansion
    ce = ClusterExpansion(cs, parameters)

    print(f"\nCluster expansion created")
    print(f"  Effective cluster interactions (ECIs):")
    for i, eci in enumerate(parameters[:10]):
        print(f"    J{i} = {eci:.6f} eV")

    print()


def example_4_formation_energy_prediction():
    """Example 4: Predict formation energies with trained CE."""
    if not ICET_AVAILABLE:
        print("Skipping: icet not available")
        return

    print("=" * 60)
    print("Example 4: Formation Energy Prediction")
    print("=" * 60)

    # Quick training (simplified from example 3)
    prim = bulk("Cu", "fcc", a=3.6)
    cutoffs = [6.0]
    cs = ClusterSpace(prim, cutoffs, chemical_symbols=[["Cu", "Au"]])

    # Minimal training set
    structures = list(
        enumerate_structures(prim, sizes=range(1, 5), chemical_symbols=[["Cu", "Au"]])
    )[:15]

    sc = StructureContainer(cs)

    cu_pure = bulk("Cu", "fcc", a=3.6)
    cu_pure.calc = EMT()
    E_Cu = cu_pure.get_potential_energy() / len(cu_pure)

    au_pure = bulk("Au", "fcc", a=4.08)
    au_pure.calc = EMT()
    E_Au = au_pure.get_potential_energy() / len(au_pure)

    for structure in structures:
        structure.calc = EMT()
        E_total = structure.get_potential_energy()

        n_Cu = structure.get_chemical_symbols().count("Cu")
        n_Au = structure.get_chemical_symbols().count("Au")
        E_ref = n_Cu * E_Cu + n_Au * E_Au
        E_form = (E_total - E_ref) / len(structure)

        sc.add_structure(structure, properties={"energy": E_form})

    from trainstation import Optimizer

    opt = Optimizer(sc.get_fit_data())
    opt.train()

    ce = ClusterExpansion(cs, opt.parameters)

    print(f"\nTrained cluster expansion (RMSE: {opt.rmse:.6f} eV/atom)")

    # Test predictions
    print(f"\nTesting predictions:")
    print(f"{'Structure':<15} {'DFT (eV/at)':<15} {'CE (eV/at)':<15} {'Error':<15}")
    print("-" * 65)

    test_structures = list(enumerate_structures(prim, sizes=[4], chemical_symbols=[["Cu", "Au"]]))[
        :5
    ]

    for structure in test_structures:
        # DFT
        structure.calc = EMT()
        E_total = structure.get_potential_energy()
        n_Cu = structure.get_chemical_symbols().count("Cu")
        n_Au = structure.get_chemical_symbols().count("Au")
        E_ref = n_Cu * E_Cu + n_Au * E_Au
        E_form_dft = (E_total - E_ref) / len(structure)

        # CE prediction
        E_form_ce = ce.predict(structure)

        error = E_form_ce - E_form_dft
        formula = structure.get_chemical_formula()

        print(f"{formula:<15} {E_form_dft:<15.6f} {E_form_ce:<15.6f} {error:+.6f}")

    print()


def example_5_monte_carlo_simulation():
    """Example 5: Monte Carlo simulation at finite temperature."""
    if not ICET_AVAILABLE:
        print("Skipping: icet not available")
        return

    try:
        from mchammer.ensembles import CanonicalEnsemble
        from mchammer.calculators import ClusterExpansionCalculator

        MCHAMMER_AVAILABLE = True
    except ImportError:
        MCHAMMER_AVAILABLE = False
        print("Skipping: mchammer not available (part of icet)")
        return

    print("=" * 60)
    print("Example 5: Monte Carlo Simulation")
    print("=" * 60)

    # Quick CE setup
    prim = bulk("Cu", "fcc", a=3.6)
    cutoffs = [6.0]
    cs = ClusterSpace(prim, cutoffs, chemical_symbols=[["Cu", "Au"]])

    # Dummy ECIs for demonstration
    parameters = [0.0] + [0.001 * i for i in range(len(cs) - 1)]
    ce = ClusterExpansion(cs, parameters)

    print(f"\nCluster expansion with {len(parameters)} parameters")

    # Create supercell for MC
    from ase.build import make_supercell

    supercell = make_supercell(prim, [[4, 0, 0], [0, 4, 0], [0, 0, 4]])

    # Set composition (50% Cu, 50% Au)
    symbols = ["Cu"] * (len(supercell) // 2) + ["Au"] * (len(supercell) // 2)
    np.random.shuffle(symbols)
    supercell.set_chemical_symbols(symbols)

    print(f"\nSupercell:")
    print(f"  Size: {len(supercell)} atoms")
    print(f"  Composition: {supercell.get_chemical_formula()}")

    # Setup MC
    calc = ClusterExpansionCalculator(supercell, ce)
    temperature = 600  # K

    print(f"\nMonte Carlo settings:")
    print(f"  Temperature: {temperature} K")
    print(f"  Ensemble: Canonical (NVT)")

    mc = CanonicalEnsemble(supercell, calc, temperature=temperature)

    # Run MC
    n_steps = 1000
    print(f"  Steps: {n_steps}")
    print(f"\nRunning...")

    mc.run(n_steps)

    # Results
    print(f"\nResults:")
    print(f"  Acceptance ratio: {mc.acceptance_ratio:.3f}")

    # Get final energy
    final_energy = calc.calculate_total(supercell) / len(supercell)
    print(f"  Final energy: {final_energy:.6f} eV/atom")

    print()


def example_6_ground_state_search():
    """Example 6: Find ground state structures."""
    if not ICET_AVAILABLE:
        print("Skipping: icet not available")
        return

    print("=" * 60)
    print("Example 6: Ground State Search")
    print("=" * 60)

    # Setup CE
    prim = bulk("Cu", "fcc", a=3.6)
    cutoffs = [6.0]
    cs = ClusterSpace(prim, cutoffs, chemical_symbols=[["Cu", "Au"]])

    # Simple ECIs favoring ordering
    parameters = [0.0, -0.05, -0.02] + [0.001] * (len(cs) - 3)
    ce = ClusterExpansion(cs, parameters)

    print(f"\nSearching for ground states...")

    # Enumerate and evaluate
    structures = enumerate_structures(prim, sizes=range(1, 9), chemical_symbols=[["Cu", "Au"]])

    results = []
    for structure in structures:
        E_form = ce.predict(structure)
        n_Cu = structure.get_chemical_symbols().count("Cu")
        n_Au = structure.get_chemical_symbols().count("Au")
        x_Au = n_Au / (n_Cu + n_Au)

        results.append(
            {
                "structure": structure,
                "energy": E_form,
                "x_Au": x_Au,
                "formula": structure.get_chemical_formula(),
            }
        )

    # Sort by composition, then energy
    results.sort(key=lambda x: (x["x_Au"], x["energy"]))

    print(f"\nEvaluated {len(results)} structures")

    # Find ground states (lowest energy at each composition)
    print(f"\nGround state line:")
    print(f"{'x(Au)':<10} {'Formula':<15} {'E_form (eV/at)':<20}")
    print("-" * 50)

    prev_x = -1
    for r in results:
        if abs(r["x_Au"] - prev_x) > 0.01:  # New composition
            print(f"{r['x_Au']:<10.3f} {r['formula']:<15} {r['energy']:<20.6f}")
            prev_x = r["x_Au"]

    print()


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("Cluster Expansion Calculation Examples")
    print("=" * 60 + "\n")

    if not ICET_AVAILABLE:
        print("ERROR: icet package not installed")
        print("Install with: pip install icet")
        print()
    else:
        try:
            example_1_basic_cluster_space()
            example_2_structure_enumeration()
            example_3_training_cluster_expansion()
            example_4_formation_energy_prediction()
            example_5_monte_carlo_simulation()
            example_6_ground_state_search()

            print("=" * 60)
            print("All examples completed successfully!")
            print("=" * 60)

        except Exception as e:
            print(f"\nError occurred: {e}")
            import traceback

            traceback.print_exc()
