#!/usr/bin/env python3
"""
Calculate adsorption energy for an adsorbate on a metal surface.

This example calculates the adsorption energy of oxygen on a Pt(111) surface
using the formula:
    E_ads = E(surface+adsorbate) - E(surface) - E(adsorbate)

A negative value indicates favorable (exothermic) adsorption.
"""

from ase.build import fcc111, add_adsorbate
from ase.calculators.emt import EMT  # Fast calculator for demonstration
from ase.optimize import BFGS
from ase.constraints import FixAtoms
from ase import Atoms


def calculate_adsorption_energy(
    metal="Pt",
    surface_indices=(1, 1, 1),
    adsorbate="O",
    size=(4, 4, 4),
    vacuum=10.0,
    adsorption_height=2.0,
    adsorption_site="fcc",
    fmax=0.05,
):
    """
    Calculate the adsorption energy of an adsorbate on a metal surface.

    Parameters
    ----------
    metal : str
        Metal element symbol (e.g., 'Pt', 'Au', 'Ag')
    surface_indices : tuple
        Miller indices of the surface (e.g., (1,1,1), (1,0,0))
    adsorbate : str
        Adsorbate species (e.g., 'O', 'H', 'CO')
    size : tuple
        Surface size as (nx, ny, layers)
    vacuum : float
        Vacuum spacing in Angstroms
    adsorption_height : float
        Initial height of adsorbate above surface in Angstroms
    adsorption_site : str
        Adsorption site ('fcc', 'hcp', 'ontop', 'bridge')
    fmax : float
        Force convergence criterion in eV/Angstrom

    Returns
    -------
    dict
        Dictionary containing energies and adsorption energy
    """

    print("=" * 60)
    print(f"Calculating adsorption of {adsorbate} on {metal}{surface_indices}")
    print("=" * 60)

    # Step 1: Create and optimize clean surface
    print("\n1. Optimizing clean surface...")
    slab = fcc111(metal, size=size, vacuum=vacuum)

    # Fix bottom two layers (common practice for surfaces)
    # This simulates the bulk material below
    n_layers = size[2]
    n_atoms_per_layer = size[0] * size[1]
    fix_indices = list(range(n_atoms_per_layer * 2))  # Fix bottom 2 layers

    slab.set_constraint(FixAtoms(indices=fix_indices))

    # Attach calculator (replace EMT with VASP/GPAW for production)
    slab.calc = EMT()

    # Optimize geometry
    opt = BFGS(slab, trajectory="clean_surface.traj", logfile="clean_surface.log")
    opt.run(fmax=fmax)

    E_slab = slab.get_potential_energy()
    print(f"   Clean surface energy: {E_slab:.3f} eV")

    # Step 2: Add adsorbate and optimize
    print(f"\n2. Optimizing surface with {adsorbate} adsorbate...")
    slab_with_ads = slab.copy()

    # Add adsorbate at specified site and height
    add_adsorbate(slab_with_ads, adsorbate, height=adsorption_height, position=adsorption_site)

    # Keep the same constraints (fix bottom layers)
    slab_with_ads.set_constraint(FixAtoms(indices=fix_indices))
    slab_with_ads.calc = EMT()

    # Optimize
    opt = BFGS(
        slab_with_ads,
        trajectory="surface_with_adsorbate.traj",
        logfile="surface_with_adsorbate.log",
    )
    opt.run(fmax=fmax)

    E_slab_ads = slab_with_ads.get_potential_energy()
    print(f"   Surface + adsorbate energy: {E_slab_ads:.3f} eV")

    # Step 3: Calculate isolated adsorbate energy
    print(f"\n3. Calculating isolated {adsorbate} energy...")

    # Create isolated adsorbate in a box
    if adsorbate == "O":
        # Atomic oxygen
        isolated = Atoms(adsorbate, positions=[(0, 0, 0)])
    elif adsorbate == "H":
        # Atomic hydrogen
        isolated = Atoms(adsorbate, positions=[(0, 0, 0)])
    elif adsorbate == "CO":
        # CO molecule
        from ase.build import molecule

        isolated = molecule(adsorbate)
    else:
        # Try to build from molecule database
        from ase.build import molecule

        try:
            isolated = molecule(adsorbate)
        except:
            # If not in database, create single atom
            isolated = Atoms(adsorbate, positions=[(0, 0, 0)])

    # Center in a box with vacuum
    isolated.center(vacuum=vacuum)
    isolated.calc = EMT()

    # For molecules, optimize; for atoms, just calculate energy
    if len(isolated) > 1:
        opt = BFGS(isolated, trajectory="isolated_adsorbate.traj")
        opt.run(fmax=fmax)

    E_ads_isolated = isolated.get_potential_energy()
    print(f"   Isolated {adsorbate} energy: {E_ads_isolated:.3f} eV")

    # Step 4: Calculate adsorption energy
    print("\n4. Calculating adsorption energy...")
    E_adsorption = E_slab_ads - E_slab - E_ads_isolated

    print("\n" + "=" * 60)
    print("RESULTS")
    print("=" * 60)
    print(f"E(surface + adsorbate) = {E_slab_ads:.3f} eV")
    print(f"E(surface)             = {E_slab:.3f} eV")
    print(f"E(adsorbate)           = {E_ads_isolated:.3f} eV")
    print(f"E_adsorption           = {E_adsorption:.3f} eV")

    if E_adsorption < 0:
        print(f"\nAdsorption is FAVORABLE (exothermic)")
        print(f"Binding strength: {abs(E_adsorption):.3f} eV")
    else:
        print(f"\nAdsorption is UNFAVORABLE (endothermic)")

    print("=" * 60)

    return {
        "E_surface_with_adsorbate": E_slab_ads,
        "E_surface": E_slab,
        "E_adsorbate": E_ads_isolated,
        "E_adsorption": E_adsorption,
    }


if __name__ == "__main__":
    # Example 1: Oxygen on Pt(111)
    results = calculate_adsorption_energy(
        metal="Pt", surface_indices=(1, 1, 1), adsorbate="O", size=(4, 4, 4), adsorption_site="fcc"
    )

    # Example 2: You can easily test different conditions
    # Uncomment to run:

    # results = calculate_adsorption_energy(
    #     metal='Au',
    #     surface_indices=(1, 1, 1),
    #     adsorbate='CO',
    #     size=(3, 3, 4),
    #     adsorption_site='ontop'
    # )

    print("\n✓ Calculation complete!")
    print("  Check output files: clean_surface.traj, surface_with_adsorbate.traj")
    print("\nNOTE: This example uses EMT calculator for speed.")
    print("      For production, replace with VASP, GPAW, or other DFT calculator.")
