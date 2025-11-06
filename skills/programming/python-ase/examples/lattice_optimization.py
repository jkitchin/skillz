#!/usr/bin/env python3
"""
Optimize the lattice constant of a bulk material.

This example demonstrates how to find the equilibrium lattice constant
by fitting an equation of state (EOS) to energy-volume data.
"""

import numpy as np
import matplotlib.pyplot as plt
from ase.build import bulk
from ase.calculators.emt import EMT  # Fast calculator for demonstration
from ase.eos import calculate_eos


def optimize_lattice_constant(
    element='Cu',
    structure='fcc',
    a_initial=3.6,
    a_range=0.3,
    n_points=9,
    plot=True,
    filename='eos.png'
):
    """
    Find the equilibrium lattice constant using equation of state fitting.

    Parameters
    ----------
    element : str
        Chemical element (e.g., 'Cu', 'Al', 'Fe')
    structure : str
        Crystal structure ('fcc', 'bcc', 'hcp', 'diamond')
    a_initial : float
        Initial guess for lattice constant in Angstroms
    a_range : float
        Range to scan around initial guess (±)
    n_points : int
        Number of points for the scan
    plot : bool
        Whether to create a plot
    filename : str
        Output filename for plot

    Returns
    -------
    dict
        Optimized lattice parameters and bulk properties
    """

    print("="*60)
    print(f"Optimizing lattice constant for {element} ({structure})")
    print("="*60)

    # Method 1: Manual scan and fitting
    print(f"\nMethod 1: Manual energy-volume scan")
    print(f"Scanning from a = {a_initial - a_range:.3f} to {a_initial + a_range:.3f} Å")

    # Create range of lattice constants to test
    a_values = np.linspace(a_initial - a_range, a_initial + a_range, n_points)

    energies = []
    volumes = []

    for a in a_values:
        # Create bulk structure with this lattice constant
        atoms = bulk(element, structure, a=a)
        atoms.calc = EMT()

        # Calculate energy
        energy = atoms.get_potential_energy()
        volume = atoms.get_volume()

        energies.append(energy)
        volumes.append(volume)

        print(f"  a = {a:.4f} Å, V = {volume:.3f} Å³, E = {energy:.4f} eV")

    # Convert to numpy arrays
    energies = np.array(energies)
    volumes = np.array(volumes)

    # Fit to equation of state
    from ase.eos import EquationOfState
    eos = EquationOfState(volumes, energies, eos='birchmurnaghan')
    v0, e0, B = eos.fit()

    # Convert optimal volume back to lattice constant
    atoms_test = bulk(element, structure, a=1.0)
    volume_per_unit_a3 = atoms_test.get_volume()  # Volume when a=1
    a_optimal = (v0 / volume_per_unit_a3) ** (1/3)

    print(f"\nFitted parameters:")
    print(f"  Optimal lattice constant: {a_optimal:.4f} Å")
    print(f"  Optimal volume per atom: {v0:.3f} Å³")
    print(f"  Minimum energy: {e0:.4f} eV")
    print(f"  Bulk modulus: {B / 1e9:.2f} GPa")

    # Method 2: Using ASE's built-in calculate_eos
    print(f"\nMethod 2: Using ASE calculate_eos function")

    atoms = bulk(element, structure, a=a_initial)
    atoms.calc = EMT()

    eos_ase = calculate_eos(atoms, trajectory='eos.traj')
    v0_ase, e0_ase, B_ase = eos_ase.fit()

    # Convert volume to lattice constant
    a_optimal_ase = (v0_ase / volume_per_unit_a3) ** (1/3)

    print(f"  Optimal lattice constant: {a_optimal_ase:.4f} Å")
    print(f"  Optimal volume per atom: {v0_ase:.3f} Å³")
    print(f"  Minimum energy: {e0_ase:.4f} eV")
    print(f"  Bulk modulus: {B_ase / 1e9:.2f} GPa")

    # Plot equation of state
    if plot:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        # Energy vs volume
        ax1.plot(volumes, energies, 'o', label='Calculated', markersize=8)

        # Fitted curve
        v_fit = np.linspace(volumes.min(), volumes.max(), 100)
        e_fit = eos.func(v_fit)
        ax1.plot(v_fit, e_fit, '-', label='EOS fit', linewidth=2)
        ax1.axvline(v0, color='r', linestyle='--', label=f'v₀ = {v0:.2f} Ų/atom')

        ax1.set_xlabel('Volume (ų/atom)', fontsize=12)
        ax1.set_ylabel('Energy (eV/atom)', fontsize=12)
        ax1.set_title(f'{element} ({structure}) - Energy vs Volume')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # Energy vs lattice constant
        a_for_volumes = (volumes / volume_per_unit_a3) ** (1/3)
        ax2.plot(a_for_volumes, energies, 'o', label='Calculated', markersize=8)

        a_fit = (v_fit / volume_per_unit_a3) ** (1/3)
        ax2.plot(a_fit, e_fit, '-', label='EOS fit', linewidth=2)
        ax2.axvline(a_optimal, color='r', linestyle='--',
                   label=f'a₀ = {a_optimal:.3f} Å')

        ax2.set_xlabel('Lattice constant (Å)', fontsize=12)
        ax2.set_ylabel('Energy (eV/atom)', fontsize=12)
        ax2.set_title(f'{element} ({structure}) - Energy vs Lattice Constant')
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(filename, dpi=300)
        print(f"\n✓ Plot saved to {filename}")

    print("\n" + "="*60)
    print("RESULTS SUMMARY")
    print("="*60)
    print(f"Element: {element} ({structure})")
    print(f"Optimal lattice constant: {a_optimal:.4f} Å")
    print(f"Cohesive energy: {e0:.4f} eV/atom")
    print(f"Bulk modulus: {B / 1e9:.2f} GPa")
    print("="*60)

    return {
        'element': element,
        'structure': structure,
        'lattice_constant': a_optimal,
        'cohesive_energy': e0,
        'bulk_modulus_GPa': B / 1e9,
        'volumes': volumes,
        'energies': energies
    }


def compare_elements(elements, structure='fcc', plot=True):
    """
    Compare lattice constants for multiple elements.

    Parameters
    ----------
    elements : list
        List of element symbols
    structure : str
        Crystal structure
    plot : bool
        Whether to create comparison plot
    """

    print("\n" + "="*60)
    print(f"Comparing lattice constants for {structure} metals")
    print("="*60)

    results = {}
    for element in elements:
        # Use typical lattice constant as initial guess
        initial_guesses = {'Cu': 3.6, 'Al': 4.05, 'Ni': 3.52, 'Ag': 4.09, 'Au': 4.08}
        a_init = initial_guesses.get(element, 4.0)

        result = optimize_lattice_constant(
            element=element,
            structure=structure,
            a_initial=a_init,
            plot=False
        )
        results[element] = result
        print()

    # Print comparison table
    print("\n" + "="*60)
    print("COMPARISON TABLE")
    print("="*60)
    print(f"{'Element':<10} {'a₀ (Å)':<12} {'E₀ (eV)':<12} {'B (GPa)':<10}")
    print("-"*60)
    for element, res in results.items():
        print(f"{element:<10} {res['lattice_constant']:<12.4f} "
              f"{res['cohesive_energy']:<12.4f} {res['bulk_modulus_GPa']:<10.2f}")
    print("="*60)

    return results


if __name__ == '__main__':
    # Example 1: Single element optimization
    print("\nExample 1: Optimizing Cu lattice constant\n")
    result = optimize_lattice_constant(
        element='Cu',
        structure='fcc',
        a_initial=3.6,
        plot=True,
        filename='cu_eos.png'
    )

    # Example 2: Compare multiple elements
    # Uncomment to run:

    # print("\n\nExample 2: Comparing multiple FCC metals\n")
    # results = compare_elements(
    #     elements=['Cu', 'Al', 'Ni', 'Ag'],
    #     structure='fcc',
    #     plot=True
    # )

    print("\n✓ Calculation complete!")
    print("\nNOTE: This example uses EMT calculator for speed.")
    print("      For production, replace with VASP, GPAW, or other DFT calculator.")
