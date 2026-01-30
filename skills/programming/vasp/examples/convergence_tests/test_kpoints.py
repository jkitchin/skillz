#!/usr/bin/env python3
"""
k-point convergence test for VASP.

Tests a range of k-point meshes and reports convergence.
"""

import os
import subprocess
import matplotlib.pyplot as plt

# Test parameters
K_MESHES = [2, 4, 6, 8, 10, 12]
CONVERGENCE_THRESHOLD = 0.001  # eV/atom
VASP_COMMAND = "mpirun -np 16 vasp_std"
ENCUT = 520  # Use converged ENCUT from test_encut.py


def setup_directory(k):
    """Create test directory and input files."""
    dirname = f"kpoint_{k}x{k}x{k}"
    os.makedirs(dirname, exist_ok=True)

    # Copy base files
    os.system(f"cp base_POSCAR {dirname}/POSCAR")
    os.system(f"cp POTCAR {dirname}/POTCAR")

    # Create KPOINTS
    with open(f"{dirname}/KPOINTS", "w") as f:
        f.write(f"Automatic mesh\n0\nGamma\n  {k} {k} {k}\n  0 0 0\n")

    # Create INCAR with fixed ENCUT
    with open("base_INCAR", "r") as f:
        lines = f.readlines()

    with open(f"{dirname}/INCAR", "w") as f:
        for line in lines:
            if line.startswith("ENCUT"):
                f.write(f"ENCUT = {ENCUT}\n")
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
            energy = float(lines[-1].split()[2])
        return energy
    except:
        print(f"Warning: Could not extract energy from {dirname}")
        return None


def main():
    """Main convergence test."""
    print("=" * 60)
    print("k-Point Convergence Test")
    print("=" * 60)
    print(f"Using ENCUT = {ENCUT} eV (ensure this is converged!)")
    print()

    # Check for POTCAR
    if not os.path.exists("POTCAR"):
        print("\nERROR: POTCAR not found!")
        print("Copy Cu POTCAR: cp $VASP_PP/potpaw_PBE/Cu/POTCAR .")
        return

    # Setup and run calculations
    energies = []
    for k in K_MESHES:
        setup_directory(k)
        run_vasp(f"kpoint_{k}x{k}x{k}")
        energy = extract_energy(f"kpoint_{k}x{k}x{k}")
        energies.append(energy)
        print(f"k-mesh = {k}×{k}×{k}: E = {energy:.6f} eV")

    # Get number of atoms
    with open("base_POSCAR", "r") as f:
        lines = f.readlines()
        n_atoms = sum([int(x) for x in lines[6].split()])

    # Analyze convergence
    print("\n" + "=" * 60)
    print("Convergence Analysis")
    print("=" * 60)
    print(f"{'k-mesh':<12} {'ΔE (meV/atom)':<20}")
    print("-" * 32)

    converged_k = None
    for i in range(len(K_MESHES) - 1):
        if energies[i] is not None and energies[i + 1] is not None:
            dE = abs(energies[i + 1] - energies[i]) / n_atoms * 1000
            k_label = f"{K_MESHES[i]}³"
            print(f"{k_label:<12} {dE:<20.3f}")

            if dE < CONVERGENCE_THRESHOLD * 1000 and converged_k is None:
                converged_k = K_MESHES[i + 1]

    print("\n" + "=" * 60)
    if converged_k:
        print(f"Converged at k-mesh = {converged_k}×{converged_k}×{converged_k}")
    else:
        print("Not converged. Increase k-mesh range.")
    print("=" * 60)

    # Plot results
    plt.figure(figsize=(10, 6))

    plt.subplot(1, 2, 1)
    plt.plot(K_MESHES, energies, "o-")
    plt.xlabel("k-mesh (N×N×N)")
    plt.ylabel("Energy (eV)")
    plt.title("Energy vs k-mesh")
    plt.grid(True)

    plt.subplot(1, 2, 2)
    dE_per_atom = []
    for i in range(len(energies) - 1):
        if energies[i] is not None and energies[i + 1] is not None:
            dE_per_atom.append(abs(energies[i + 1] - energies[i]) / n_atoms * 1000)
    plt.plot(K_MESHES[:-1], dE_per_atom, "o-")
    plt.axhline(
        y=CONVERGENCE_THRESHOLD * 1000,
        color="r",
        linestyle="--",
        label=f"Threshold ({CONVERGENCE_THRESHOLD * 1000} meV/atom)",
    )
    plt.xlabel("k-mesh (N×N×N)")
    plt.ylabel("|ΔE| (meV/atom)")
    plt.title("Convergence")
    plt.legend()
    plt.grid(True)
    plt.yscale("log")

    plt.tight_layout()
    plt.savefig("kpoint_convergence.png", dpi=150)
    print(f"\nPlot saved: kpoint_convergence.png")


if __name__ == "__main__":
    main()
