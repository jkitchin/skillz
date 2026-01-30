"""
Nudged Elastic Band (NEB) Calculation Example

Demonstrates reaction barrier calculations using NEB:
1. Basic NEB calculation
2. Climbing image NEB
3. Analysis of reaction path
4. Surface diffusion barrier

Reference:
Henkelman, Uberuaga & Jónsson, "A climbing image nudged elastic band method
for finding saddle points and minimum energy paths," J. Chem. Phys. 113,
9901 (2000)
"""

from ase import Atoms
from ase.build import fcc111, molecule, add_adsorbate
from ase.optimize import BFGS
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.neb import NEB, NEBTools
from ase.io import read, write
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def example_1_basic_neb():
    """Example 1: Basic NEB calculation for surface diffusion."""
    print("=" * 60)
    print("Example 1: Basic NEB - Surface Diffusion")
    print("=" * 60)

    # Create clean slab
    slab = fcc111("Cu", size=(3, 3, 4), vacuum=10.0)
    mask = [atom.position[2] < 6.0 for atom in slab]
    slab.set_constraint(FixAtoms(mask=mask))

    # Initial state: adsorbate at fcc site
    initial = slab.copy()
    add_adsorbate(initial, "Au", height=1.5, position="fcc")
    initial.set_constraint(FixAtoms(mask=mask))
    initial.calc = EMT()

    opt = BFGS(initial, logfile=None)
    opt.run(fmax=0.02)
    E_initial = initial.get_potential_energy()

    print(f"\nInitial state optimized:")
    print(f"  Energy: {E_initial:.4f} eV")

    # Final state: adsorbate at neighboring fcc site
    final = slab.copy()
    # Move to adjacent fcc site
    add_adsorbate(final, "Au", height=1.5, position=(2.5, 1.5))
    final.set_constraint(FixAtoms(mask=mask))
    final.calc = EMT()

    opt = BFGS(final, logfile=None)
    opt.run(fmax=0.02)
    E_final = final.get_potential_energy()

    print(f"\nFinal state optimized:")
    print(f"  Energy: {E_final:.4f} eV")
    print(f"  Reaction energy: {E_final - E_initial:.4f} eV")

    # Create NEB with 5 intermediate images
    images = [initial]
    images += [initial.copy() for i in range(5)]
    images += [final]

    # Interpolate
    neb = NEB(images)
    neb.interpolate()

    # Set calculators for intermediate images
    for image in images[1:-1]:
        image.set_constraint(FixAtoms(mask=mask))
        image.calc = EMT()

    print(f"\nRunning NEB with {len(images)} images...")
    print(f"  (This may take a minute)")

    # Optimize NEB
    optimizer = BFGS(neb, trajectory="neb.traj", logfile=None)
    optimizer.run(fmax=0.05)

    print(f"  Optimization steps: {optimizer.get_number_of_steps()}")

    # Analyze results
    nebtools = NEBTools(images)

    # Get barrier
    Ef, dE = nebtools.get_barrier()

    print(f"\nResults:")
    print(f"  Forward barrier: {Ef:.4f} eV")
    print(f"  Reaction energy: {dE:.4f} eV")
    print(f"  Reverse barrier: {Ef - dE:.4f} eV")

    # Plot energy path
    energies = [image.get_potential_energy() - E_initial for image in images]

    plt.figure(figsize=(8, 6))
    plt.plot(energies, "o-", linewidth=2, markersize=8)
    plt.axhline(y=0, color="k", linestyle="--", alpha=0.3)
    plt.axhline(y=Ef, color="r", linestyle="--", alpha=0.5, label=f"Barrier: {Ef:.3f} eV")
    plt.xlabel("Image Number", fontsize=12)
    plt.ylabel("Energy (eV)", fontsize=12)
    plt.title("NEB Energy Path", fontsize=14)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("neb_basic.png", dpi=150)
    print(f"\nSaved: neb_basic.png")
    print()


def example_2_climbing_image_neb():
    """Example 2: Climbing image NEB for accurate barrier."""
    print("=" * 60)
    print("Example 2: Climbing Image NEB")
    print("=" * 60)

    # Use same setup as example 1 but with climbing image
    slab = fcc111("Cu", size=(3, 3, 4), vacuum=10.0)
    mask = [atom.position[2] < 6.0 for atom in slab]

    # Initial state
    initial = slab.copy()
    add_adsorbate(initial, "Au", height=1.5, position="fcc")
    initial.set_constraint(FixAtoms(mask=mask))
    initial.calc = EMT()
    opt = BFGS(initial, logfile=None)
    opt.run(fmax=0.02)
    E_initial = initial.get_potential_energy()

    # Final state
    final = slab.copy()
    add_adsorbate(final, "Au", height=1.5, position=(2.5, 1.5))
    final.set_constraint(FixAtoms(mask=mask))
    final.calc = EMT()
    opt = BFGS(final, logfile=None)
    opt.run(fmax=0.02)

    # Create NEB
    images = [initial]
    images += [initial.copy() for i in range(7)]  # More images for CI-NEB
    images += [final]

    neb = NEB(images, climb=True)  # Enable climbing image
    neb.interpolate()

    for image in images[1:-1]:
        image.set_constraint(FixAtoms(mask=mask))
        image.calc = EMT()

    print(f"\nRunning CI-NEB with {len(images)} images...")

    optimizer = BFGS(neb, trajectory="cineb.traj", logfile=None)
    optimizer.run(fmax=0.05)

    # Analyze
    nebtools = NEBTools(images)
    Ef, dE = nebtools.get_barrier()

    print(f"\nCI-NEB Results:")
    print(f"  Forward barrier: {Ef:.4f} eV")
    print(f"  Reaction energy: {dE:.4f} eV")

    # Find highest energy image (should be transition state)
    energies = [img.get_potential_energy() for img in images]
    max_idx = np.argmax(energies)

    print(f"  Transition state at image {max_idx}")
    print(f"  TS energy: {energies[max_idx]:.4f} eV")

    print(f"\nNote: CI-NEB gives more accurate barrier by pushing")
    print(f"highest image to exact transition state")
    print()


def example_3_reaction_path_analysis():
    """Example 3: Detailed analysis of reaction path."""
    print("=" * 60)
    print("Example 3: Reaction Path Analysis")
    print("=" * 60)

    # Run a simple NEB (reuse setup)
    slab = fcc111("Cu", size=(3, 3, 4), vacuum=10.0)
    mask = [atom.position[2] < 6.0 for atom in slab]

    initial = slab.copy()
    add_adsorbate(initial, "Au", height=1.5, position="fcc")
    initial.set_constraint(FixAtoms(mask=mask))
    initial.calc = EMT()
    opt = BFGS(initial, logfile=None)
    opt.run(fmax=0.02)

    final = slab.copy()
    add_adsorbate(final, "Au", height=1.5, position=(2.5, 1.5))
    final.set_constraint(FixAtoms(mask=mask))
    final.calc = EMT()
    opt = BFGS(final, logfile=None)
    opt.run(fmax=0.02)

    images = [initial]
    images += [initial.copy() for i in range(5)]
    images += [final]

    neb = NEB(images)
    neb.interpolate()

    for image in images[1:-1]:
        image.set_constraint(FixAtoms(mask=mask))
        image.calc = EMT()

    optimizer = BFGS(neb, logfile=None)
    optimizer.run(fmax=0.05)

    # Detailed analysis
    nebtools = NEBTools(images)

    print(f"\nReaction path analysis:")

    # Energies
    energies = [img.get_potential_energy() for img in images]
    E_ref = energies[0]
    rel_energies = [E - E_ref for E in energies]

    print(f"\n  Image energies (relative to initial):")
    for i, E in enumerate(rel_energies):
        marker = " (TS)" if E == max(rel_energies) else ""
        print(f"    Image {i}: {E:.4f} eV{marker}")

    # Forces on images
    print(f"\n  Max forces on images:")
    for i, img in enumerate(images):
        if i == 0 or i == len(images) - 1:
            continue
        forces = img.get_forces()
        max_force = np.max(np.linalg.norm(forces, axis=1))
        print(f"    Image {i}: {max_force:.4f} eV/Å")

    # Adsorbate positions along path
    print(f"\n  Adsorbate position along path:")
    for i, img in enumerate(images):
        ads_pos = img.positions[-1]  # Last atom is adsorbate
        print(f"    Image {i}: ({ads_pos[0]:.2f}, {ads_pos[1]:.2f}, {ads_pos[2]:.2f})")

    # Get fitted barrier
    Ef_fitted, dE_fitted = nebtools.get_barrier(fit=True)
    Ef_raw, dE_raw = nebtools.get_barrier(fit=False)

    print(f"\n  Barrier comparison:")
    print(f"    Raw (max image): {Ef_raw:.4f} eV")
    print(f"    Fitted (spline): {Ef_fitted:.4f} eV")
    print(f"    Difference: {abs(Ef_fitted - Ef_raw):.4f} eV")

    print()


def example_4_multiple_barriers():
    """Example 4: Path with multiple barriers."""
    print("=" * 60)
    print("Example 4: Multiple Transition States")
    print("=" * 60)

    print("\nNote: This is a simplified demonstration.")
    print("Real multi-barrier paths require careful initial/final states.")

    # Create a more complex path (artificial example)
    # Initial: fcc site
    slab = fcc111("Cu", size=(4, 4, 4), vacuum=10.0)
    mask = [atom.position[2] < 6.0 for atom in slab]

    initial = slab.copy()
    add_adsorbate(initial, "Au", height=1.5, position="fcc")
    initial.set_constraint(FixAtoms(mask=mask))
    initial.calc = EMT()
    opt = BFGS(initial, logfile=None)
    opt.run(fmax=0.02)

    # Final: distant fcc site (requires crossing multiple sites)
    final = slab.copy()
    add_adsorbate(final, "Au", height=1.5, position=(4.5, 4.5))
    final.set_constraint(FixAtoms(mask=mask))
    final.calc = EMT()
    opt = BFGS(final, logfile=None)
    opt.run(fmax=0.02)

    # More images needed for longer path
    images = [initial]
    images += [initial.copy() for i in range(11)]
    images += [final]

    neb = NEB(images)
    neb.interpolate()

    for image in images[1:-1]:
        image.set_constraint(FixAtoms(mask=mask))
        image.calc = EMT()

    print(f"\nOptimizing NEB with {len(images)} images...")
    optimizer = BFGS(neb, logfile=None)
    optimizer.run(fmax=0.05)

    # Analyze path
    energies = [img.get_potential_energy() for img in images]
    E_ref = energies[0]
    rel_energies = [E - E_ref for E in energies]

    # Find local maxima (transition states)
    ts_indices = []
    for i in range(1, len(rel_energies) - 1):
        if rel_energies[i] > rel_energies[i - 1] and rel_energies[i] > rel_energies[i + 1]:
            ts_indices.append(i)

    print(f"\nFound {len(ts_indices)} transition state(s):")
    for idx in ts_indices:
        print(f"  Image {idx}: {rel_energies[idx]:.4f} eV")

    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(rel_energies, "o-", linewidth=2, markersize=6)
    for idx in ts_indices:
        plt.plot(
            idx, rel_energies[idx], "r*", markersize=15, label="TS" if idx == ts_indices[0] else ""
        )
    plt.axhline(y=0, color="k", linestyle="--", alpha=0.3)
    plt.xlabel("Image Number", fontsize=12)
    plt.ylabel("Relative Energy (eV)", fontsize=12)
    plt.title("Multi-Barrier Diffusion Path", fontsize=14)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("neb_multiple_barriers.png", dpi=150)
    print(f"\nSaved: neb_multiple_barriers.png")
    print()


def example_5_neb_convergence():
    """Example 5: Test convergence with number of images."""
    print("=" * 60)
    print("Example 5: Image Number Convergence")
    print("=" * 60)

    slab = fcc111("Cu", size=(3, 3, 4), vacuum=10.0)
    mask = [atom.position[2] < 6.0 for atom in slab]

    # Initial and final states
    initial = slab.copy()
    add_adsorbate(initial, "Au", height=1.5, position="fcc")
    initial.set_constraint(FixAtoms(mask=mask))
    initial.calc = EMT()
    opt = BFGS(initial, logfile=None)
    opt.run(fmax=0.02)

    final = slab.copy()
    add_adsorbate(final, "Au", height=1.5, position=(2.5, 1.5))
    final.set_constraint(FixAtoms(mask=mask))
    final.calc = EMT()
    opt = BFGS(final, logfile=None)
    opt.run(fmax=0.02)

    # Test different numbers of intermediate images
    n_images_list = [3, 5, 7, 9]
    barriers = []

    print(f"\nTesting image number convergence:")
    print(f"{'N_images':<12} {'Barrier (eV)':<15} {'Steps':<10}")
    print("-" * 40)

    for n_intermediate in n_images_list:
        images = [initial]
        images += [initial.copy() for i in range(n_intermediate)]
        images += [final]

        neb = NEB(images)
        neb.interpolate()

        for image in images[1:-1]:
            image.set_constraint(FixAtoms(mask=mask))
            image.calc = EMT()

        optimizer = BFGS(neb, logfile=None)
        optimizer.run(fmax=0.05)

        nebtools = NEBTools(images)
        Ef, dE = nebtools.get_barrier()
        barriers.append(Ef)

        print(f"{len(images):<12} {Ef:<15.4f} {optimizer.get_number_of_steps():<10}")

    # Check convergence
    if len(barriers) > 1:
        convergence = abs(barriers[-1] - barriers[-2])
        print(f"\nConvergence: {convergence:.4f} eV")
        print(f"Recommended: {n_images_list[-2]} intermediate images")

    print()


def example_6_save_trajectory():
    """Example 6: Save and visualize NEB trajectory."""
    print("=" * 60)
    print("Example 6: NEB Trajectory Visualization")
    print("=" * 60)

    slab = fcc111("Cu", size=(3, 3, 4), vacuum=10.0)
    mask = [atom.position[2] < 6.0 for atom in slab]

    initial = slab.copy()
    add_adsorbate(initial, "Au", height=1.5, position="fcc")
    initial.set_constraint(FixAtoms(mask=mask))
    initial.calc = EMT()
    opt = BFGS(initial, logfile=None)
    opt.run(fmax=0.02)

    final = slab.copy()
    add_adsorbate(final, "Au", height=1.5, position=(2.5, 1.5))
    final.set_constraint(FixAtoms(mask=mask))
    final.calc = EMT()
    opt = BFGS(final, logfile=None)
    opt.run(fmax=0.02)

    images = [initial]
    images += [initial.copy() for i in range(5)]
    images += [final]

    neb = NEB(images)
    neb.interpolate()

    for image in images[1:-1]:
        image.set_constraint(FixAtoms(mask=mask))
        image.calc = EMT()

    print(f"\nRunning NEB and saving trajectory...")
    optimizer = BFGS(neb, trajectory="neb_full.traj")
    optimizer.run(fmax=0.05)

    # Save final path
    write("neb_path.traj", images)

    print(f"\nTrajectory files created:")
    print(f"  neb_full.traj  - Full optimization trajectory")
    print(f"  neb_path.traj  - Final converged path")

    # Create energy profile plot
    nebtools = NEBTools(images)
    Ef, dE = nebtools.get_barrier(fit=True)

    # Get fitted path
    s, E, Sfit, Efit = nebtools.get_fit()

    plt.figure(figsize=(10, 6))
    plt.plot(s, E, "o", label="Images", markersize=8)
    plt.plot(Sfit, Efit, "-", label="Spline fit", linewidth=2)
    plt.axhline(y=E[0], color="g", linestyle="--", alpha=0.5, label="Initial")
    plt.axhline(y=E[-1], color="b", linestyle="--", alpha=0.5, label="Final")
    plt.xlabel("Reaction Coordinate", fontsize=12)
    plt.ylabel("Energy (eV)", fontsize=12)
    plt.title(f"NEB Energy Profile (Barrier: {Ef:.3f} eV)", fontsize=14)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("neb_profile.png", dpi=150)

    print(f"  neb_profile.png - Energy profile plot")
    print(f"\nBarrier: {Ef:.4f} eV")
    print()


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("Nudged Elastic Band (NEB) Examples")
    print("=" * 60 + "\n")

    try:
        example_1_basic_neb()
        example_2_climbing_image_neb()
        example_3_reaction_path_analysis()
        example_4_multiple_barriers()
        example_5_neb_convergence()
        example_6_save_trajectory()

        print("=" * 60)
        print("All examples completed successfully!")
        print("=" * 60)

    except Exception as e:
        print(f"\nError occurred: {e}")
        import traceback

        traceback.print_exc()
