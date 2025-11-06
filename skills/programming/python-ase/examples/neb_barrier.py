#!/usr/bin/env python3
"""
Calculate reaction barrier using Nudged Elastic Band (NEB) method.

This example demonstrates how to find the minimum energy path and
activation barrier for a simple surface diffusion process.
"""

import numpy as np
import matplotlib.pyplot as plt
from ase import Atoms
from ase.build import fcc111, add_adsorbate
from ase.calculators.emt import EMT
from ase.neb import NEB
from ase.optimize import BFGS
from ase.constraints import FixAtoms
from ase.io import read, write


def create_initial_state(size=(4, 4, 4), adsorbate_position='fcc'):
    """
    Create initial state with adsorbate at one site.

    Parameters
    ----------
    size : tuple
        Surface size as (nx, ny, layers)
    adsorbate_position : str
        Initial adsorption site

    Returns
    -------
    Atoms
        Initial state structure
    """
    slab = fcc111('Pt', size=size, vacuum=10.0)

    # Fix bottom two layers
    n_layers = size[2]
    n_atoms_per_layer = size[0] * size[1]
    constraint = FixAtoms(indices=range(n_atoms_per_layer * 2))
    slab.set_constraint(constraint)

    # Add adsorbate
    add_adsorbate(slab, 'O', height=2.0, position=adsorbate_position)

    return slab


def create_final_state(size=(4, 4, 4), adsorbate_position='hcp'):
    """
    Create final state with adsorbate at different site.

    Parameters
    ----------
    size : tuple
        Surface size as (nx, ny, layers)
    adsorbate_position : str
        Final adsorption site

    Returns
    -------
    Atoms
        Final state structure
    """
    slab = fcc111('Pt', size=size, vacuum=10.0)

    # Fix bottom two layers (same as initial)
    n_layers = size[2]
    n_atoms_per_layer = size[0] * size[1]
    constraint = FixAtoms(indices=range(n_atoms_per_layer * 2))
    slab.set_constraint(constraint)

    # Add adsorbate at different position
    add_adsorbate(slab, 'O', height=2.0, position=adsorbate_position)

    return slab


def calculate_reaction_barrier(
    n_images=7,
    spring_constant=1.0,
    fmax=0.05,
    optimize_endpoints=True,
    plot=True
):
    """
    Calculate reaction barrier using NEB.

    Parameters
    ----------
    n_images : int
        Number of intermediate images (not counting endpoints)
    spring_constant : float
        Spring constant between images in eV/Ų
    fmax : float
        Force convergence criterion in eV/Å
    optimize_endpoints : bool
        Whether to optimize initial and final states first
    plot : bool
        Whether to create energy profile plot

    Returns
    -------
    dict
        Barrier heights and energies
    """

    print("="*60)
    print("NEB Calculation: O diffusion on Pt(111)")
    print("="*60)

    # Step 1: Create and optimize initial state
    print("\n1. Creating initial state (O at FCC site)...")
    initial = create_initial_state(adsorbate_position='fcc')
    initial.calc = EMT()

    if optimize_endpoints:
        print("   Optimizing initial state...")
        opt = BFGS(initial, trajectory='initial.traj', logfile='initial.log')
        opt.run(fmax=fmax)

    E_initial = initial.get_potential_energy()
    print(f"   Initial state energy: {E_initial:.3f} eV")

    # Step 2: Create and optimize final state
    print("\n2. Creating final state (O at HCP site)...")
    final = create_final_state(adsorbate_position='hcp')
    final.calc = EMT()

    if optimize_endpoints:
        print("   Optimizing final state...")
        opt = BFGS(final, trajectory='final.traj', logfile='final.log')
        opt.run(fmax=fmax)

    E_final = final.get_potential_energy()
    print(f"   Final state energy: {E_final:.3f} eV")

    # Step 3: Set up NEB calculation
    print(f"\n3. Setting up NEB with {n_images} intermediate images...")

    # Create image list: [initial] + intermediates + [final]
    images = [initial]

    # Create intermediate images
    for i in range(n_images):
        image = initial.copy()
        image.calc = EMT()
        images.append(image)

    images.append(final)

    print(f"   Total images: {len(images)} (including endpoints)")

    # Create NEB object
    neb = NEB(images, k=spring_constant, climb=True)  # climb=True for CI-NEB

    # Interpolate linearly between initial and final
    print("   Interpolating between initial and final states...")
    neb.interpolate()

    # Step 4: Optimize the band
    print("\n4. Optimizing NEB...")
    print(f"   Convergence criterion: max force < {fmax} eV/Å")

    opt = BFGS(neb, trajectory='neb.traj', logfile='neb.log')

    # Run optimization with progress reporting
    def print_progress():
        forces = [np.sqrt((image.get_forces()**2).sum(axis=1).max())
                 for image in images[1:-1]]
        max_force = max(forces)
        energies = [image.get_potential_energy() for image in images]
        barrier = max(energies) - energies[0]
        print(f"      Step {opt.nsteps}: max force = {max_force:.4f} eV/Å, "
              f"barrier = {barrier:.3f} eV")

    opt.attach(print_progress, interval=5)
    opt.run(fmax=fmax)

    # Step 5: Analyze results
    print("\n5. Analyzing results...")

    energies = [image.get_potential_energy() for image in images]
    energies = np.array(energies)

    # Shift energies relative to initial state
    energies_relative = energies - energies[0]

    # Find transition state
    ts_index = np.argmax(energies)
    E_ts = energies[ts_index]

    # Calculate barriers
    forward_barrier = E_ts - E_initial
    reverse_barrier = E_ts - E_final
    reaction_energy = E_final - E_initial

    print("\n" + "="*60)
    print("RESULTS")
    print("="*60)
    print(f"Initial state energy:      {E_initial:.3f} eV")
    print(f"Final state energy:        {E_final:.3f} eV")
    print(f"Transition state energy:   {E_ts:.3f} eV")
    print(f"Transition state at image: {ts_index}")
    print(f"\nForward barrier:           {forward_barrier:.3f} eV")
    print(f"Reverse barrier:           {reverse_barrier:.3f} eV")
    print(f"Reaction energy:           {reaction_energy:.3f} eV")
    print("="*60)

    # Save transition state
    write('transition_state.traj', images[ts_index])
    print(f"\nTransition state saved to: transition_state.traj")

    # Plot energy profile
    if plot:
        fig, ax = plt.subplots(figsize=(10, 6))

        # Reaction coordinate (image index)
        x = np.arange(len(images))

        # Plot energy profile
        ax.plot(x, energies_relative, 'o-', linewidth=2, markersize=10,
               label='Energy profile')

        # Mark initial, transition, and final states
        ax.plot(0, energies_relative[0], 'go', markersize=15, label='Initial')
        ax.plot(ts_index, energies_relative[ts_index], 'ro', markersize=15,
               label='Transition state')
        ax.plot(len(images)-1, energies_relative[-1], 'bo', markersize=15,
               label='Final')

        # Add horizontal lines for reference
        ax.axhline(y=0, color='g', linestyle='--', alpha=0.5)
        ax.axhline(y=forward_barrier, color='r', linestyle='--', alpha=0.5)

        # Add barrier annotation
        ax.annotate(f'E_a = {forward_barrier:.3f} eV',
                   xy=(ts_index, energies_relative[ts_index]),
                   xytext=(ts_index + 1, energies_relative[ts_index] + 0.1),
                   arrowprops=dict(arrowstyle='->', color='red'),
                   fontsize=12, color='red')

        ax.set_xlabel('Reaction coordinate (image index)', fontsize=12)
        ax.set_ylabel('Relative energy (eV)', fontsize=12)
        ax.set_title('NEB Energy Profile: O diffusion FCC → HCP on Pt(111)',
                    fontsize=14)
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig('neb_profile.png', dpi=300)
        print(f"Energy profile saved to: neb_profile.png")

    return {
        'E_initial': E_initial,
        'E_final': E_final,
        'E_transition_state': E_ts,
        'forward_barrier': forward_barrier,
        'reverse_barrier': reverse_barrier,
        'reaction_energy': reaction_energy,
        'energies': energies,
        'ts_index': ts_index
    }


def analyze_neb_trajectory(trajectory_file='neb.traj'):
    """
    Analyze convergence of NEB optimization.

    Parameters
    ----------
    trajectory_file : str
        Path to NEB trajectory file
    """
    print("\n" + "="*60)
    print("Analyzing NEB convergence")
    print("="*60)

    traj = read(trajectory_file, ':')
    n_steps = len(traj) // (len(traj[0]))  # Number of optimization steps

    print(f"Total optimization steps: {n_steps}")

    # Plot convergence
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    # Energy convergence
    barriers = []
    for i in range(0, len(traj), len(traj) // n_steps):
        energies = [atoms.get_potential_energy() for atoms in traj[i:i+7]]
        barrier = max(energies) - energies[0]
        barriers.append(barrier)

    ax1.plot(barriers, 'o-', linewidth=2)
    ax1.set_xlabel('Optimization step')
    ax1.set_ylabel('Activation barrier (eV)')
    ax1.set_title('Barrier convergence')
    ax1.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('neb_convergence.png', dpi=300)
    print("Convergence plot saved to: neb_convergence.png")


if __name__ == '__main__':
    print("\nNEB Example: Oxygen diffusion on Pt(111) surface")
    print("Process: O atom moving from FCC site to HCP site\n")

    # Run NEB calculation
    results = calculate_reaction_barrier(
        n_images=5,  # Start with fewer images for faster calculation
        spring_constant=1.0,
        fmax=0.05,
        optimize_endpoints=True,
        plot=True
    )

    # Optionally analyze convergence
    # analyze_neb_trajectory('neb.traj')

    print("\n✓ NEB calculation complete!")
    print("\nOutput files:")
    print("  - initial.traj: Initial state")
    print("  - final.traj: Final state")
    print("  - neb.traj: Full NEB trajectory")
    print("  - transition_state.traj: Transition state structure")
    print("  - neb_profile.png: Energy profile plot")
    print("\nNOTE: This example uses EMT calculator for speed.")
    print("      For production, replace with VASP, GPAW, or other DFT calculator.")
    print("      Real NEB calculations may need 7-11 images and tighter convergence.")
