#!/usr/bin/env python3
"""
ENCUT convergence test for VASP.

Tests a range of energy cutoffs and reports convergence.
"""

import os
import subprocess
import matplotlib.pyplot as plt

# Test parameters
ENCUT_VALUES = [300, 350, 400, 450, 500, 550, 600]
CONVERGENCE_THRESHOLD = 0.001  # eV/atom
VASP_COMMAND = "mpirun -np 16 vasp_std"


def setup_directory(encut):
    """Create test directory and input files."""
    dirname = f"encut_{encut}"
    os.makedirs(dirname, exist_ok=True)

    # Copy base files
    os.system(f"cp base_POSCAR {dirname}/POSCAR")
    os.system(f"cp POTCAR {dirname}/POTCAR")

    # Create KPOINTS (fixed at 4×4×4 for speed)
    with open(f"{dirname}/KPOINTS", "w") as f:
        f.write("Automatic mesh\n0\nGamma\n  4 4 4\n  0 0 0\n")

    # Create INCAR with modified ENCUT
    with open("base_INCAR", "r") as f:
        lines = f.readlines()

    with open(f"{dirname}/INCAR", "w") as f:
        for line in lines:
            if line.startswith("ENCUT"):
                f.write(f"ENCUT = {encut}\n")
            else:
                f.write(line)


def run_vasp(dirname):
    """Run VASP in directory."""
    print(f"Running VASP in {dirname}...")
    cwd = os.getcwd()
    os.chdir(dirname)

    with open("vasp.out", "w") as f:
        subprocess.run(VASP_COMMAND, shell=True, stdout=f, stderr=f)

    os.chdir(cwd)


def extract_energy(dirname):
    """Extract final energy from OSZICAR."""
    try:
        with open(f"{dirname}/OSZICAR", "r") as f:
            lines = f.readlines()
            # Last line has final energy
            energy = float(lines[-1].split()[2])
        return energy
    except:
        print(f"Warning: Could not extract energy from {dirname}")
        return None


def main():
    """Main convergence test."""
    print("=" * 60)
    print("ENCUT Convergence Test")
    print("=" * 60)

    # Check for POTCAR
    if not os.path.exists("POTCAR"):
        print("\nERROR: POTCAR not found!")
        print("Copy Cu POTCAR: cp $VASP_PP/potpaw_PBE/Cu/POTCAR .")
        return

    # Setup and run calculations
    energies = []
    for encut in ENCUT_VALUES:
        setup_directory(encut)
        run_vasp(f"encut_{encut}")
        energy = extract_energy(f"encut_{encut}")
        energies.append(energy)
        print(f"ENCUT = {encut} eV: E = {energy:.6f} eV")

    # Get number of atoms
    with open("base_POSCAR", "r") as f:
        lines = f.readlines()
        n_atoms = sum([int(x) for x in lines[6].split()])

    # Analyze convergence
    print("\n" + "=" * 60)
    print("Convergence Analysis")
    print("=" * 60)
    print(f"{'ENCUT (eV)':<12} {'ΔE (meV/atom)':<20}")
    print("-" * 32)

    converged_encut = None
    for i in range(len(ENCUT_VALUES) - 1):
        if energies[i] is not None and energies[i + 1] is not None:
            dE = abs(energies[i + 1] - energies[i]) / n_atoms * 1000  # meV/atom
            print(f"{ENCUT_VALUES[i]:<12} {dE:<20.3f}")

            if dE < CONVERGENCE_THRESHOLD * 1000 and converged_encut is None:
                converged_encut = ENCUT_VALUES[i + 1]

    print("\n" + "=" * 60)
    if converged_encut:
        print(f"Converged at ENCUT = {converged_encut} eV")
    else:
        print("Not converged. Increase ENCUT range.")
    print("=" * 60)

    # Plot results
    plt.figure(figsize=(10, 6))

    plt.subplot(1, 2, 1)
    plt.plot(ENCUT_VALUES, energies, "o-")
    plt.xlabel("ENCUT (eV)")
    plt.ylabel("Energy (eV)")
    plt.title("Energy vs ENCUT")
    plt.grid(True)

    plt.subplot(1, 2, 2)
    dE_per_atom = []
    for i in range(len(energies) - 1):
        if energies[i] is not None and energies[i + 1] is not None:
            dE_per_atom.append(abs(energies[i + 1] - energies[i]) / n_atoms * 1000)
    plt.plot(ENCUT_VALUES[:-1], dE_per_atom, "o-")
    plt.axhline(
        y=CONVERGENCE_THRESHOLD * 1000,
        color="r",
        linestyle="--",
        label=f"Threshold ({CONVERGENCE_THRESHOLD * 1000} meV/atom)",
    )
    plt.xlabel("ENCUT (eV)")
    plt.ylabel("|ΔE| (meV/atom)")
    plt.title("Convergence")
    plt.legend()
    plt.grid(True)
    plt.yscale("log")

    plt.tight_layout()
    plt.savefig("encut_convergence.png", dpi=150)
    print(f"\nPlot saved: encut_convergence.png")


if __name__ == "__main__":
    main()
