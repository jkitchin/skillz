"""
Surface Energy Calculation Example

Demonstrates surface energy calculations using slab models:
1. Basic surface energy calculation
2. Convergence testing (slab thickness, vacuum)
3. Different surface facets
4. Asymmetric slabs

Formula: γ = (E_slab - N × E_bulk) / (2 × A)

Reference:
Fiorentini & Methfessel, "Extracting convergent surface energies from slab
calculations," J. Phys.: Condens. Matter 8, 6525 (1996)
"""

from ase.build import bulk, fcc111, fcc100, surface
from ase.optimize import BFGS
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def get_bulk_energy_per_atom(element="Cu", crystal="fcc", a=3.6):
    """Calculate bulk energy per atom."""
    bulk_atoms = bulk(element, crystal, a=a)
    bulk_atoms.calc = EMT()
    E_bulk = bulk_atoms.get_potential_energy()
    E_per_atom = E_bulk / len(bulk_atoms)
    return E_per_atom


def example_1_basic_surface_energy():
    """Example 1: Basic surface energy calculation for Cu(111)."""
    print("=" * 60)
    print("Example 1: Basic Surface Energy - Cu(111)")
    print("=" * 60)

    # Bulk energy per atom
    E_bulk_per_atom = get_bulk_energy_per_atom("Cu", "fcc", 3.6)
    print(f"\nBulk energy per atom: {E_bulk_per_atom:.4f} eV")

    # Create (111) slab
    slab = fcc111("Cu", size=(3, 3, 7), vacuum=10.0)

    print(f"\nSlab properties:")
    print(f"  Number of atoms: {len(slab)}")
    print(f"  Cell dimensions: {slab.cell.cellpar()[:3]}")

    # Relax slab (fix bottom 2 layers)
    mask = [atom.position[2] < 6.0 for atom in slab]
    slab.set_constraint(FixAtoms(mask=mask))

    slab.calc = EMT()
    opt = BFGS(slab, logfile=None)
    opt.run(fmax=0.02)

    # Calculate surface energy
    E_slab = slab.get_potential_energy()
    N = len(slab)

    # Surface area (cross product of first two cell vectors)
    cell = slab.get_cell()
    A = np.linalg.norm(np.cross(cell[0], cell[1]))

    # Surface energy (factor of 2 for two surfaces)
    gamma = (E_slab - N * E_bulk_per_atom) / (2 * A)

    print(f"\nResults:")
    print(f"  Slab energy: {E_slab:.4f} eV")
    print(f"  Surface area: {A:.3f} ų")
    print(f"  Surface energy: {gamma:.6f} eV/ų")
    print(f"  Surface energy: {gamma * 1000:.2f} meV/ų")
    print(f"  Surface energy: {gamma * 16.0217:.2f} J/m²")
    print()


def example_2_convergence_slab_thickness():
    """Example 2: Convergence test for slab thickness."""
    print("=" * 60)
    print("Example 2: Slab Thickness Convergence")
    print("=" * 60)

    E_bulk_per_atom = get_bulk_energy_per_atom("Cu", "fcc", 3.6)

    # Test different number of layers
    n_layers = [5, 7, 9, 11, 13]
    surface_energies = []

    print(f"\nTesting slab thickness:")
    print(f"{'Layers':<10} {'Atoms':<10} {'γ (meV/ų)':<15}")
    print("-" * 40)

    for n in n_layers:
        slab = fcc111("Cu", size=(3, 3, n), vacuum=10.0)

        # Fix bottom 2 layers
        z_positions = [atom.position[2] for atom in slab]
        z_min = min(z_positions)
        mask = [atom.position[2] < (z_min + 4.0) for atom in slab]
        slab.set_constraint(FixAtoms(mask=mask))

        slab.calc = EMT()
        opt = BFGS(slab, logfile=None)
        opt.run(fmax=0.02)

        E_slab = slab.get_potential_energy()
        N = len(slab)
        cell = slab.get_cell()
        A = np.linalg.norm(np.cross(cell[0], cell[1]))

        gamma = (E_slab - N * E_bulk_per_atom) / (2 * A)
        surface_energies.append(gamma * 1000)  # Convert to meV/ų

        print(f"{n:<10} {N:<10} {gamma * 1000:<15.2f}")

    # Check convergence
    converged = abs(surface_energies[-1] - surface_energies[-2]) < 1.0
    print(f"\nConverged (< 1 meV/ų): {converged}")
    print(f"Recommended: {n_layers[-2]} layers")

    # Plot
    plt.figure(figsize=(8, 6))
    plt.plot(n_layers, surface_energies, "o-", linewidth=2, markersize=8)
    plt.xlabel("Number of Layers", fontsize=12)
    plt.ylabel("Surface Energy (meV/ų)", fontsize=12)
    plt.title("Cu(111) Surface Energy Convergence", fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("surface_convergence_layers.png", dpi=150)
    print(f"\nSaved: surface_convergence_layers.png")
    print()


def example_3_convergence_vacuum():
    """Example 3: Convergence test for vacuum thickness."""
    print("=" * 60)
    print("Example 3: Vacuum Thickness Convergence")
    print("=" * 60)

    E_bulk_per_atom = get_bulk_energy_per_atom("Cu", "fcc", 3.6)

    # Test different vacuum thicknesses
    vacuum_sizes = [6, 8, 10, 12, 15, 20]
    surface_energies = []

    print(f"\nTesting vacuum thickness:")
    print(f"{'Vacuum (Å)':<15} {'γ (meV/ų)':<15}")
    print("-" * 35)

    for vac in vacuum_sizes:
        slab = fcc111("Cu", size=(3, 3, 7), vacuum=vac)

        # Fix bottom 2 layers
        z_positions = [atom.position[2] for atom in slab]
        z_min = min(z_positions)
        mask = [atom.position[2] < (z_min + 4.0) for atom in slab]
        slab.set_constraint(FixAtoms(mask=mask))

        slab.calc = EMT()
        opt = BFGS(slab, logfile=None)
        opt.run(fmax=0.02)

        E_slab = slab.get_potential_energy()
        N = len(slab)
        cell = slab.get_cell()
        A = np.linalg.norm(np.cross(cell[0], cell[1]))

        gamma = (E_slab - N * E_bulk_per_atom) / (2 * A)
        surface_energies.append(gamma * 1000)

        print(f"{vac:<15} {gamma * 1000:<15.2f}")

    # Check convergence
    converged = abs(surface_energies[-1] - surface_energies[-2]) < 0.5
    print(f"\nConverged (< 0.5 meV/ų): {converged}")
    print(f"Recommended: {vacuum_sizes[-3]} Å")

    # Plot
    plt.figure(figsize=(8, 6))
    plt.plot(vacuum_sizes, surface_energies, "o-", linewidth=2, markersize=8)
    plt.xlabel("Vacuum Thickness (Å)", fontsize=12)
    plt.ylabel("Surface Energy (meV/ų)", fontsize=12)
    plt.title("Cu(111) Surface Energy vs Vacuum", fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("surface_convergence_vacuum.png", dpi=150)
    print(f"\nSaved: surface_convergence_vacuum.png")
    print()


def example_4_different_facets():
    """Example 4: Surface energies for different facets."""
    print("=" * 60)
    print("Example 4: Different Surface Facets")
    print("=" * 60)

    E_bulk_per_atom = get_bulk_energy_per_atom("Cu", "fcc", 3.6)

    facets = [
        ("(111)", fcc111("Cu", size=(3, 3, 7), vacuum=10)),
        ("(100)", fcc100("Cu", size=(3, 3, 7), vacuum=10)),
    ]

    results = []

    print(f"\nCalculating surface energies:")
    print(f"{'Facet':<10} {'Atoms':<10} {'γ (meV/ų)':<15} {'γ (J/m²)':<15}")
    print("-" * 55)

    for facet_name, slab in facets:
        # Fix bottom 2 layers
        z_positions = [atom.position[2] for atom in slab]
        z_min = min(z_positions)
        mask = [atom.position[2] < (z_min + 4.0) for atom in slab]
        slab.set_constraint(FixAtoms(mask=mask))

        slab.calc = EMT()
        opt = BFGS(slab, logfile=None)
        opt.run(fmax=0.02)

        E_slab = slab.get_potential_energy()
        N = len(slab)
        cell = slab.get_cell()
        A = np.linalg.norm(np.cross(cell[0], cell[1]))

        gamma = (E_slab - N * E_bulk_per_atom) / (2 * A)
        gamma_meV = gamma * 1000
        gamma_Jm2 = gamma * 16.0217

        results.append((facet_name, gamma_meV))

        print(f"{facet_name:<10} {N:<10} {gamma_meV:<15.2f} {gamma_Jm2:<15.3f}")

    # Surface energy ordering
    print(f"\nStability order (lowest to highest):")
    sorted_results = sorted(results, key=lambda x: x[1])
    for i, (facet, gamma) in enumerate(sorted_results, 1):
        print(f"  {i}. {facet}: {gamma:.2f} meV/ų")

    print(f"\nTypically for FCC metals: γ(111) < γ(100) < γ(110)")
    print()


def example_5_surface_relaxation():
    """Example 5: Analyze surface relaxation."""
    print("=" * 60)
    print("Example 5: Surface Relaxation Analysis")
    print("=" * 60)

    # Create unrelaxed slab
    slab_unrelaxed = fcc111("Cu", size=(3, 3, 7), vacuum=10.0)
    z_positions_init = [atom.position[2] for atom in slab_unrelaxed]

    # Relax slab
    slab = slab_unrelaxed.copy()

    # Fix bottom 2 layers
    z_min = min(z_positions_init)
    mask = [atom.position[2] < (z_min + 4.0) for atom in slab]
    slab.set_constraint(FixAtoms(mask=mask))

    slab.calc = EMT()
    opt = BFGS(slab, trajectory="surface_relax.traj", logfile=None)
    opt.run(fmax=0.01)

    # Analyze layer spacings
    z_positions_final = [atom.position[2] for atom in slab]

    # Group atoms by layers (z-position)
    unique_z = sorted(set(np.round(z_positions_final, 2)))

    print(f"\nLayer analysis:")
    print(f"  Number of layers: {len(unique_z)}")

    # Calculate interlayer spacings
    if len(unique_z) > 1:
        spacings = np.diff(unique_z)
        bulk_spacing = spacings[len(spacings) // 2]  # Middle spacing

        print(f"\n  Interlayer spacings:")
        for i, d in enumerate(spacings):
            rel_change = ((d - bulk_spacing) / bulk_spacing) * 100
            print(f"    d{i + 1}-{i + 2}: {d:.3f} Å ({rel_change:+.1f}%)")

        # Surface relaxation (top layer)
        if len(spacings) > 2:
            d12 = spacings[-1]  # Top spacing
            rel_relax = ((d12 - bulk_spacing) / bulk_spacing) * 100
            print(f"\n  Top layer relaxation: {rel_relax:+.2f}%")

            if rel_relax < 0:
                print(f"  Contraction (inward relaxation)")
            else:
                print(f"  Expansion (outward relaxation)")

    print()


def example_6_surface_energy_formula():
    """Example 6: Demonstrate surface energy formula components."""
    print("=" * 60)
    print("Example 6: Surface Energy Formula Breakdown")
    print("=" * 60)

    E_bulk_per_atom = get_bulk_energy_per_atom("Cu", "fcc", 3.6)

    slab = fcc111("Cu", size=(3, 3, 7), vacuum=10.0)

    # Fix bottom layers
    z_positions = [atom.position[2] for atom in slab]
    z_min = min(z_positions)
    mask = [atom.position[2] < (z_min + 4.0) for atom in slab]
    slab.set_constraint(FixAtoms(mask=mask))

    slab.calc = EMT()
    opt = BFGS(slab, logfile=None)
    opt.run(fmax=0.02)

    # Components
    E_slab = slab.get_potential_energy()
    N = len(slab)
    E_bulk_total = N * E_bulk_per_atom
    E_excess = E_slab - E_bulk_total

    cell = slab.get_cell()
    A = np.linalg.norm(np.cross(cell[0], cell[1]))

    gamma = E_excess / (2 * A)

    print(f"\nFormula: γ = (E_slab - N × E_bulk) / (2 × A)")
    print(f"\nComponents:")
    print(f"  E_slab       = {E_slab:.4f} eV")
    print(f"  N            = {N} atoms")
    print(f"  E_bulk/atom  = {E_bulk_per_atom:.4f} eV")
    print(f"  N × E_bulk   = {E_bulk_total:.4f} eV")
    print(f"  E_excess     = {E_excess:.4f} eV")
    print(f"  A            = {A:.3f} ų")
    print(f"  Factor of 2  = two surfaces")
    print(f"\nResult:")
    print(f"  γ = {gamma:.6f} eV/ų")
    print(f"  γ = {gamma * 1000:.2f} meV/ų")
    print(f"  γ = {gamma * 16.0217:.2f} J/m²")
    print()


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("Surface Energy Calculation Examples")
    print("=" * 60 + "\n")

    try:
        example_1_basic_surface_energy()
        example_2_convergence_slab_thickness()
        example_3_convergence_vacuum()
        example_4_different_facets()
        example_5_surface_relaxation()
        example_6_surface_energy_formula()

        print("=" * 60)
        print("All examples completed successfully!")
        print("=" * 60)

    except Exception as e:
        print(f"\nError occurred: {e}")
        import traceback

        traceback.print_exc()
