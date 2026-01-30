"""
Basic Structure Relaxation Example

Demonstrates fundamental structure optimization workflows:
1. Geometry optimization (fixed cell)
2. Unit cell optimization (variable cell)
3. Lattice constant determination
4. Space group identification

Uses EMT calculator for fast execution.
"""

from ase.build import bulk
from ase.optimize import BFGS
from ase.constraints import ExpCellFilter
from ase.calculators.emt import EMT
import spglib
import numpy as np


def example_1_geometry_optimization():
    """Example 1: Basic geometry optimization with fixed cell."""
    print("=" * 60)
    print("Example 1: Geometry Optimization (Fixed Cell)")
    print("=" * 60)

    # Create slightly distorted Cu structure
    atoms = bulk("Cu", "fcc", a=3.7)  # Intentionally wrong lattice constant
    atoms.rattle(stdev=0.1)  # Add random displacement

    # Set calculator
    atoms.calc = EMT()

    # Initial state
    E_init = atoms.get_potential_energy()
    forces_init = atoms.get_forces()
    max_force_init = np.max(np.linalg.norm(forces_init, axis=1))

    print(f"\nInitial state:")
    print(f"  Energy: {E_init:.4f} eV")
    print(f"  Max force: {max_force_init:.4f} eV/Å")

    # Optimize
    opt = BFGS(atoms, trajectory="geometry_opt.traj")
    opt.run(fmax=0.01)

    # Final state
    E_final = atoms.get_potential_energy()
    forces_final = atoms.get_forces()
    max_force_final = np.max(np.linalg.norm(forces_final, axis=1))

    print(f"\nFinal state:")
    print(f"  Energy: {E_final:.4f} eV")
    print(f"  Max force: {max_force_final:.4f} eV/Å")
    print(f"  Energy change: {E_final - E_init:.4f} eV")
    print(f"  Steps: {opt.get_number_of_steps()}")
    print()


def example_2_cell_optimization():
    """Example 2: Unit cell optimization (positions and cell)."""
    print("=" * 60)
    print("Example 2: Unit Cell Optimization")
    print("=" * 60)

    # Create Cu with wrong lattice constant
    atoms = bulk("Cu", "fcc", a=3.7)
    atoms.calc = EMT()

    # Initial cell
    a_init = atoms.cell.cellpar()[0]
    V_init = atoms.get_volume()
    E_init = atoms.get_potential_energy()

    print(f"\nInitial state:")
    print(f"  Lattice constant: {a_init:.3f} Å")
    print(f"  Volume: {V_init:.3f} ų")
    print(f"  Energy: {E_init:.4f} eV")

    # Optimize both cell and positions
    ecf = ExpCellFilter(atoms)
    opt = BFGS(ecf, trajectory="cell_opt.traj")
    opt.run(fmax=0.01)

    # Final cell
    a_final = atoms.cell.cellpar()[0]
    V_final = atoms.get_volume()
    E_final = atoms.get_potential_energy()

    print(f"\nFinal state:")
    print(f"  Lattice constant: {a_final:.3f} Å")
    print(f"  Volume: {V_final:.3f} ų")
    print(f"  Energy: {E_final:.4f} eV")
    print(f"  Lattice change: {a_final - a_init:.3f} Å")
    print(f"  Energy change: {E_final - E_init:.4f} eV")
    print()


def example_3_lattice_constant():
    """Example 3: Accurate lattice constant via equation of state."""
    print("=" * 60)
    print("Example 3: Lattice Constant from EOS")
    print("=" * 60)

    from ase.eos import EquationOfState

    # Scan over lattice constants
    lattice_constants = np.linspace(3.4, 3.8, 11)
    energies = []
    volumes = []

    print("\nScanning lattice constants...")
    for a in lattice_constants:
        atoms = bulk("Cu", "fcc", a=a)
        atoms.calc = EMT()
        E = atoms.get_potential_energy()
        V = atoms.get_volume()

        energies.append(E)
        volumes.append(V)

        print(f"  a = {a:.3f} Å: E = {E:.4f} eV")

    # Fit equation of state
    eos = EquationOfState(volumes, energies, eos="birchmurnaghan")
    v0, e0, B = eos.fit()

    # Convert volume to lattice constant (for FCC: V = a³/4)
    a0 = (4 * v0) ** (1 / 3)

    print(f"\nFitted parameters:")
    print(f"  Equilibrium lattice constant: {a0:.4f} Å")
    print(f"  Equilibrium energy: {e0:.4f} eV")
    print(f"  Bulk modulus: {B / 1e9:.2f} GPa")

    # Compare to experimental value (Cu: ~3.615 Å)
    print(f"\n  (Experimental Cu: ~3.615 Å)")
    print()


def example_4_space_group():
    """Example 4: Space group and symmetry analysis."""
    print("=" * 60)
    print("Example 4: Space Group Identification")
    print("=" * 60)

    # Create and optimize structure
    atoms = bulk("Cu", "fcc", a=3.6)
    atoms.calc = EMT()

    # Optimize to get clean structure
    ecf = ExpCellFilter(atoms)
    opt = BFGS(ecf, logfile=None)
    opt.run(fmax=0.001)

    # Get space group using spglib
    cell = (atoms.cell, atoms.get_scaled_positions(), atoms.get_atomic_numbers())

    # Space group
    spacegroup = spglib.get_spacegroup(cell, symprec=1e-5)
    print(f"\nSpace group: {spacegroup}")

    # Detailed symmetry dataset
    dataset = spglib.get_symmetry_dataset(cell, symprec=1e-5)

    print(f"\nDetailed symmetry information:")
    print(f"  International symbol: {dataset['international']}")
    print(f"  Hall symbol: {dataset['hall']}")
    print(f"  Space group number: {dataset['number']}")
    print(f"  Point group: {dataset['pointgroup']}")

    # Lattice parameters
    a, b, c, alpha, beta, gamma = atoms.cell.cellpar()
    print(f"\nLattice parameters:")
    print(f"  a = {a:.4f} Å")
    print(f"  b = {b:.4f} Å")
    print(f"  c = {c:.4f} Å")
    print(f"  α = {alpha:.2f}°")
    print(f"  β = {beta:.2f}°")
    print(f"  γ = {gamma:.2f}°")

    # Crystal system
    crystal_systems = [
        "triclinic",
        "monoclinic",
        "orthorhombic",
        "tetragonal",
        "trigonal",
        "hexagonal",
        "cubic",
    ]
    crystal_system = crystal_systems[(dataset["number"] - 1) // 15]
    print(f"  Crystal system: {crystal_system}")

    # Symmetry operations
    rotations, translations = spglib.get_symmetry(cell)
    print(f"\nNumber of symmetry operations: {len(rotations)}")
    print()


def example_5_different_structures():
    """Example 5: Optimize different crystal structures."""
    print("=" * 60)
    print("Example 5: Multiple Crystal Structures")
    print("=" * 60)

    structures = [
        ("Cu", "fcc", 3.6),
        ("Fe", "bcc", 2.87),
        ("Al", "fcc", 4.05),
    ]

    print("\nOptimizing different structures:")
    print("-" * 40)

    for element, structure, a_init in structures:
        atoms = bulk(element, structure, a=a_init)
        atoms.calc = EMT()

        # Optimize cell
        ecf = ExpCellFilter(atoms)
        opt = BFGS(ecf, logfile=None)
        opt.run(fmax=0.01)

        # Results
        a_opt = atoms.cell.cellpar()[0]
        E_per_atom = atoms.get_potential_energy() / len(atoms)

        # Space group
        cell = (atoms.cell, atoms.get_scaled_positions(), atoms.get_atomic_numbers())
        spacegroup = spglib.get_spacegroup(cell, symprec=1e-5)

        print(f"\n{element} ({structure.upper()}):")
        print(f"  Lattice constant: {a_opt:.4f} Å")
        print(f"  Energy per atom: {E_per_atom:.4f} eV")
        print(f"  Space group: {spacegroup}")

    print()


def example_6_trajectory_analysis():
    """Example 6: Analyze optimization trajectory."""
    print("=" * 60)
    print("Example 6: Trajectory Analysis")
    print("=" * 60)

    from ase.io import read
    import matplotlib

    matplotlib.use("Agg")  # Non-interactive backend
    import matplotlib.pyplot as plt

    # Create and optimize
    atoms = bulk("Cu", "fcc", a=3.7)
    atoms.rattle(stdev=0.2)
    atoms.calc = EMT()

    opt = BFGS(atoms, trajectory="analysis.traj")
    opt.run(fmax=0.01)

    # Read trajectory
    traj = read("analysis.traj", ":")

    # Extract properties
    energies = [a.get_potential_energy() for a in traj]
    forces = [a.get_forces() for a in traj]
    max_forces = [np.max(np.linalg.norm(f, axis=1)) for f in forces]
    volumes = [a.get_volume() for a in traj]

    print(f"\nOptimization summary:")
    print(f"  Total steps: {len(traj)}")
    print(f"  Initial energy: {energies[0]:.4f} eV")
    print(f"  Final energy: {energies[-1]:.4f} eV")
    print(f"  Energy change: {energies[-1] - energies[0]:.4f} eV")
    print(f"  Initial max force: {max_forces[0]:.4f} eV/Å")
    print(f"  Final max force: {max_forces[-1]:.4f} eV/Å")

    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

    ax1.plot(energies, "o-")
    ax1.set_xlabel("Step")
    ax1.set_ylabel("Energy (eV)")
    ax1.set_title("Energy vs Optimization Step")
    ax1.grid(True, alpha=0.3)

    ax2.plot(max_forces, "o-")
    ax2.set_xlabel("Step")
    ax2.set_ylabel("Max Force (eV/Å)")
    ax2.set_title("Force Convergence")
    ax2.axhline(y=0.01, color="r", linestyle="--", label="fmax=0.01")
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_yscale("log")

    plt.tight_layout()
    plt.savefig("trajectory_analysis.png", dpi=150)
    print(f"\nSaved plot: trajectory_analysis.png")
    print()


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("Basic Structure Relaxation Examples")
    print("=" * 60 + "\n")

    try:
        example_1_geometry_optimization()
        example_2_cell_optimization()
        example_3_lattice_constant()
        example_4_space_group()
        example_5_different_structures()
        example_6_trajectory_analysis()

        print("=" * 60)
        print("All examples completed successfully!")
        print("=" * 60)

    except Exception as e:
        print(f"\nError occurred: {e}")
        import traceback

        traceback.print_exc()
