"""
Phonon Calculation Example using Phonopy

Demonstrates phonon calculations:
1. Basic phonon calculation
2. Phonon band structure
3. Phonon density of states
4. Thermal properties

Reference:
Togo & Tanaka, "First principles phonon calculations in materials science,"
Scr. Mater. 108, 1 (2015)

Note: Requires phonopy package: pip install phonopy
"""

from ase.build import bulk
from ase.calculators.emt import EMT
import numpy as np

try:
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms

    PHONOPY_AVAILABLE = True
except ImportError:
    PHONOPY_AVAILABLE = False
    print("Warning: phonopy not installed. Install with: pip install phonopy")


def ase_to_phonopy(atoms):
    """Convert ASE atoms to Phonopy atoms."""
    return PhonopyAtoms(
        symbols=atoms.get_chemical_symbols(), cell=atoms.get_cell(), positions=atoms.get_positions()
    )


def example_1_basic_phonon():
    """Example 1: Basic phonon calculation."""
    if not PHONOPY_AVAILABLE:
        print("Skipping: phonopy not available")
        return

    print("=" * 60)
    print("Example 1: Basic Phonon Calculation - Si")
    print("=" * 60)

    # Create structure
    atoms = bulk("Si", "diamond", a=5.43)
    atoms.calc = EMT()

    print(f"\nStructure:")
    print(f"  Formula: {atoms.get_chemical_formula()}")
    print(f"  Atoms: {len(atoms)}")
    print(f"  Lattice: {atoms.cell.cellpar()[0]:.3f} Å")

    # Convert to phonopy
    phonopy_atoms = ase_to_phonopy(atoms)

    # Create phonopy object with 2×2×2 supercell
    phonon = Phonopy(phonopy_atoms, supercell_matrix=[[2, 0, 0], [0, 2, 0], [0, 0, 2]])

    print(f"\nSupercell:")
    print(f"  Matrix: 2×2×2")
    print(f"  Atoms in supercell: {len(phonon.supercell)}")

    # Generate displacements
    phonon.generate_displacements(distance=0.01)
    supercells = phonon.supercells_with_displacements

    print(f"  Displaced configurations: {len(supercells)}")

    # Calculate forces for each displacement
    print(f"\nCalculating forces...")
    forces = []
    for i, scell in enumerate(supercells):
        # Convert back to ASE
        cell = scell.cell
        positions = scell.positions
        numbers = scell.numbers

        from ase import Atoms as ASEAtoms

        ase_scell = ASEAtoms(numbers=numbers, positions=positions, cell=cell, pbc=True)
        ase_scell.calc = EMT()

        f = ase_scell.get_forces()
        forces.append(f)

        if (i + 1) % 10 == 0 or i == len(supercells) - 1:
            print(f"  Progress: {i + 1}/{len(supercells)}")

    # Set forces
    phonon.forces = forces

    # Produce force constants
    phonon.produce_force_constants()

    print(f"\nForce constants calculated")
    print()


def example_2_phonon_band_structure():
    """Example 2: Calculate and plot phonon band structure."""
    if not PHONOPY_AVAILABLE:
        print("Skipping: phonopy not available")
        return

    print("=" * 60)
    print("Example 2: Phonon Band Structure")
    print("=" * 60)

    # Setup (similar to example 1)
    atoms = bulk("Cu", "fcc", a=3.6)
    atoms.calc = EMT()

    phonopy_atoms = ase_to_phonopy(atoms)

    phonon = Phonopy(phonopy_atoms, supercell_matrix=[[2, 0, 0], [0, 2, 0], [0, 0, 2]])

    phonon.generate_displacements(distance=0.01)

    # Calculate forces (abbreviated)
    forces = []
    for scell in phonon.supercells_with_displacements:
        from ase import Atoms as ASEAtoms

        ase_scell = ASEAtoms(
            numbers=scell.numbers, positions=scell.positions, cell=scell.cell, pbc=True
        )
        ase_scell.calc = EMT()
        forces.append(ase_scell.get_forces())

    phonon.forces = forces
    phonon.produce_force_constants()

    # Set band structure path
    # For FCC: Γ-X-W-K-Γ-L
    path = [
        [[0.0, 0.0, 0.0], [0.5, 0.0, 0.5]],  # Γ-X
        [[0.5, 0.0, 0.5], [0.5, 0.25, 0.75]],  # X-W
        [[0.5, 0.25, 0.75], [0.375, 0.375, 0.75]],  # W-K
        [[0.375, 0.375, 0.75], [0.0, 0.0, 0.0]],  # K-Γ
        [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],  # Γ-L
    ]

    labels = ["$\\Gamma$", "X", "W", "K", "$\\Gamma$", "L"]

    # Calculate band structure
    phonon.run_band_structure(path, labels=labels)

    # Plot
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        phonon.plot_band_structure()
        plt.ylabel("Frequency (THz)")
        plt.title("Cu Phonon Band Structure")
        plt.tight_layout()
        plt.savefig("phonon_bands.png", dpi=150)
        print(f"\nSaved: phonon_bands.png")
    except Exception as e:
        print(f"\nPlotting skipped: {e}")

    print()


def example_3_phonon_dos():
    """Example 3: Phonon density of states."""
    if not PHONOPY_AVAILABLE:
        print("Skipping: phonopy not available")
        return

    print("=" * 60)
    print("Example 3: Phonon Density of States")
    print("=" * 60)

    # Setup
    atoms = bulk("Al", "fcc", a=4.05)
    atoms.calc = EMT()

    phonopy_atoms = ase_to_phonopy(atoms)
    phonon = Phonopy(phonopy_atoms, supercell_matrix=[[2, 0, 0], [0, 2, 0], [0, 0, 2]])

    phonon.generate_displacements(distance=0.01)

    # Calculate forces
    forces = []
    for scell in phonon.supercells_with_displacements:
        from ase import Atoms as ASEAtoms

        ase_scell = ASEAtoms(
            numbers=scell.numbers, positions=scell.positions, cell=scell.cell, pbc=True
        )
        ase_scell.calc = EMT()
        forces.append(ase_scell.get_forces())

    phonon.forces = forces
    phonon.produce_force_constants()

    # Calculate DOS
    phonon.run_mesh([20, 20, 20])
    phonon.run_total_dos()

    # Plot
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        phonon.plot_total_dos()
        plt.xlabel("Frequency (THz)")
        plt.ylabel("DOS")
        plt.title("Al Phonon DOS")
        plt.tight_layout()
        plt.savefig("phonon_dos.png", dpi=150)
        print(f"\nSaved: phonon_dos.png")
    except Exception as e:
        print(f"\nPlotting skipped: {e}")

    print()


def example_4_thermal_properties():
    """Example 4: Calculate thermal properties."""
    if not PHONOPY_AVAILABLE:
        print("Skipping: phonopy not available")
        return

    print("=" * 60)
    print("Example 4: Thermal Properties")
    print("=" * 60)

    # Setup
    atoms = bulk("Cu", "fcc", a=3.6)
    atoms.calc = EMT()

    phonopy_atoms = ase_to_phonopy(atoms)
    phonon = Phonopy(phonopy_atoms, supercell_matrix=[[2, 0, 0], [0, 2, 0], [0, 0, 2]])

    phonon.generate_displacements(distance=0.01)

    forces = []
    for scell in phonon.supercells_with_displacements:
        from ase import Atoms as ASEAtoms

        ase_scell = ASEAtoms(
            numbers=scell.numbers, positions=scell.positions, cell=scell.cell, pbc=True
        )
        ase_scell.calc = EMT()
        forces.append(ase_scell.get_forces())

    phonon.forces = forces
    phonon.produce_force_constants()

    # Calculate thermal properties
    phonon.run_mesh([20, 20, 20])
    phonon.run_thermal_properties(t_min=0, t_max=1000, t_step=10)

    # Get thermal properties
    tp_dict = phonon.get_thermal_properties_dict()

    temperatures = tp_dict["temperatures"]
    free_energy = tp_dict["free_energy"]
    entropy = tp_dict["entropy"]
    heat_capacity = tp_dict["heat_capacity"]

    print(f"\nThermal properties:")
    print(f"{'T (K)':<10} {'F (kJ/mol)':<15} {'S (J/K/mol)':<15} {'Cv (J/K/mol)':<15}")
    print("-" * 58)

    # Print selected temperatures
    for i in [0, 10, 20, 30, 50, 100]:
        if i < len(temperatures):
            T = temperatures[i]
            F = free_energy[i]
            S = entropy[i]
            Cv = heat_capacity[i]
            print(f"{T:<10.1f} {F:<15.3f} {S:<15.3f} {Cv:<15.3f}")

    # Plot
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))

        ax1.plot(temperatures, free_energy)
        ax1.set_xlabel("Temperature (K)")
        ax1.set_ylabel("Free Energy (kJ/mol)")
        ax1.grid(True, alpha=0.3)

        ax2.plot(temperatures, entropy)
        ax2.set_xlabel("Temperature (K)")
        ax2.set_ylabel("Entropy (J/K/mol)")
        ax2.grid(True, alpha=0.3)

        ax3.plot(temperatures, heat_capacity)
        ax3.set_xlabel("Temperature (K)")
        ax3.set_ylabel("Heat Capacity (J/K/mol)")
        ax3.grid(True, alpha=0.3)

        # Debye model comparison for Cv
        # Classical limit: Cv → 3R for high T
        R = 8.314  # J/K/mol
        ax3.axhline(y=3 * R, color="r", linestyle="--", label="Classical limit (3R)")
        ax3.legend()

        ax4.axis("off")
        info_text = f"Material: Cu (FCC)\n"
        info_text += f"Supercell: 2×2×2\n"
        info_text += f"Mesh: 20×20×20\n"
        info_text += f"\nAt 300 K:\n"
        idx_300 = np.argmin(np.abs(temperatures - 300))
        info_text += f"  Cv = {heat_capacity[idx_300]:.2f} J/K/mol\n"
        info_text += f"  S = {entropy[idx_300]:.2f} J/K/mol"
        ax4.text(0.1, 0.5, info_text, fontsize=12, family="monospace")

        plt.tight_layout()
        plt.savefig("thermal_properties.png", dpi=150)
        print(f"\nSaved: thermal_properties.png")

    except Exception as e:
        print(f"\nPlotting error: {e}")

    print()


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("Phonon Calculation Examples")
    print("=" * 60 + "\n")

    if not PHONOPY_AVAILABLE:
        print("ERROR: phonopy package not installed")
        print("Install with: pip install phonopy")
        print()
    else:
        try:
            example_1_basic_phonon()
            example_2_phonon_band_structure()
            example_3_phonon_dos()
            example_4_thermal_properties()

            print("=" * 60)
            print("All examples completed successfully!")
            print("=" * 60)

        except Exception as e:
            print(f"\nError occurred: {e}")
            import traceback

            traceback.print_exc()
