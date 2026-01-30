"""
Adsorption Energy Calculation Example

Demonstrates adsorption energy calculations:
1. Basic adsorption energy
2. Different adsorption sites
3. Coverage effects
4. Optimal geometry search

Formula: E_ads = E_slab+ads - E_slab - E_molecule

Reference:
Hammer & Nørskov, "Theoretical surface science and catalysis,"
Adv. Catal. 45, 71 (2000)
"""

from ase.build import fcc111, molecule, add_adsorbate
from ase.optimize import BFGS
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
import numpy as np


def example_1_basic_adsorption():
    """Example 1: Basic CO adsorption on Cu(111)."""
    print("=" * 60)
    print("Example 1: Basic Adsorption Energy - CO on Cu(111)")
    print("=" * 60)

    # 1. Clean slab energy
    print("\n1. Calculating clean slab energy...")
    slab = fcc111("Cu", size=(3, 3, 4), vacuum=10.0)

    # Fix bottom 2 layers
    mask = [atom.position[2] < 6.0 for atom in slab]
    slab.set_constraint(FixAtoms(mask=mask))

    slab.calc = EMT()
    opt = BFGS(slab, logfile=None)
    opt.run(fmax=0.02)

    E_slab = slab.get_potential_energy()
    print(f"   E_slab = {E_slab:.4f} eV")

    # 2. Molecule in gas phase
    print("\n2. Calculating gas phase molecule...")
    mol = molecule("CO")
    mol.center(vacuum=10.0)
    mol.calc = EMT()
    opt_mol = BFGS(mol, logfile=None)
    opt_mol.run(fmax=0.01)

    E_mol = mol.get_potential_energy()
    print(f"   E_CO = {E_mol:.4f} eV")

    # 3. Adsorbed system
    print("\n3. Calculating adsorbed system...")
    slab_ads = slab.copy()
    add_adsorbate(slab_ads, mol, height=2.0, position="ontop")

    # Fix bottom layers
    slab_ads.set_constraint(FixAtoms(mask=mask))

    slab_ads.calc = EMT()
    opt_ads = BFGS(slab_ads, logfile=None)
    opt_ads.run(fmax=0.02)

    E_total = slab_ads.get_potential_energy()
    print(f"   E_slab+CO = {E_total:.4f} eV")

    # 4. Adsorption energy
    E_ads = E_total - E_slab - E_mol

    print(f"\n4. Adsorption energy:")
    print(f"   E_ads = E_total - E_slab - E_mol")
    print(f"   E_ads = {E_total:.4f} - {E_slab:.4f} - {E_mol:.4f}")
    print(f"   E_ads = {E_ads:.4f} eV")

    if E_ads < 0:
        print(f"   Exothermic adsorption (favorable)")
    else:
        print(f"   Endothermic adsorption (unfavorable)")

    print()


def example_2_adsorption_sites():
    """Example 2: Compare different adsorption sites."""
    print("=" * 60)
    print("Example 2: Adsorption at Different Sites")
    print("=" * 60)

    # Clean slab
    slab_clean = fcc111("Cu", size=(3, 3, 4), vacuum=10.0)
    mask = [atom.position[2] < 6.0 for atom in slab_clean]
    slab_clean.set_constraint(FixAtoms(mask=mask))
    slab_clean.calc = EMT()

    opt = BFGS(slab_clean, logfile=None)
    opt.run(fmax=0.02)
    E_slab = slab_clean.get_potential_energy()

    # Gas phase molecule
    mol = molecule("O")
    mol.center(vacuum=10.0)
    mol.calc = EMT()
    E_mol = mol.get_potential_energy()

    # Test different sites
    sites = ["ontop", "bridge", "fcc", "hcp"]
    adsorption_energies = {}

    print(f"\nTesting adsorption sites:")
    print(f"{'Site':<12} {'E_ads (eV)':<15} {'Favorable':<12}")
    print("-" * 42)

    for site in sites:
        slab_ads = slab_clean.copy()

        # Add adsorbate
        add_adsorbate(slab_ads, mol, height=1.5, position=site)
        slab_ads.set_constraint(FixAtoms(mask=mask))
        slab_ads.calc = EMT()

        # Relax
        opt = BFGS(slab_ads, logfile=None)
        opt.run(fmax=0.02)

        # Calculate adsorption energy
        E_total = slab_ads.get_potential_energy()
        E_ads = E_total - E_slab - E_mol
        adsorption_energies[site] = E_ads

        favorable = "Yes" if E_ads < 0 else "No"
        print(f"{site:<12} {E_ads:<15.4f} {favorable:<12}")

    # Find preferred site
    preferred_site = min(adsorption_energies, key=adsorption_energies.get)
    print(f"\nPreferred site: {preferred_site}")
    print(f"Strongest binding: {adsorption_energies[preferred_site]:.4f} eV")

    # Typically for FCC(111): fcc/hcp (hollow) > bridge > ontop
    print(f"\nTypical order: hollow (fcc/hcp) > bridge > ontop")
    print()


def example_3_height_optimization():
    """Example 3: Optimize adsorbate height."""
    print("=" * 60)
    print("Example 3: Adsorbate Height Optimization")
    print("=" * 60)

    # Clean slab
    slab = fcc111("Cu", size=(3, 3, 4), vacuum=10.0)
    mask = [atom.position[2] < 6.0 for atom in slab]
    slab.set_constraint(FixAtoms(mask=mask))
    slab.calc = EMT()
    opt = BFGS(slab, logfile=None)
    opt.run(fmax=0.02)
    E_slab = slab.get_potential_energy()

    # Gas phase
    mol = molecule("O")
    mol.center(vacuum=10.0)
    mol.calc = EMT()
    E_mol = mol.get_potential_energy()

    # Scan heights
    heights = np.linspace(1.0, 3.0, 11)
    energies = []

    print(f"\nScanning initial heights:")
    print(f"{'Height (Å)':<15} {'E_ads (eV)':<15}")
    print("-" * 32)

    for h in heights:
        slab_ads = slab.copy()
        add_adsorbate(slab_ads, mol, height=h, position="fcc")
        slab_ads.set_constraint(FixAtoms(mask=mask))
        slab_ads.calc = EMT()

        # Relax
        opt = BFGS(slab_ads, logfile=None)
        opt.run(fmax=0.02)

        E_total = slab_ads.get_potential_energy()
        E_ads = E_total - E_slab - E_mol
        energies.append(E_ads)

        print(f"{h:<15.2f} {E_ads:<15.4f}")

    # Find optimal
    optimal_idx = np.argmin(energies)
    optimal_height = heights[optimal_idx]
    optimal_energy = energies[optimal_idx]

    print(f"\nOptimal initial height: {optimal_height:.2f} Å")
    print(f"Strongest binding: {optimal_energy:.4f} eV")
    print()


def example_4_coverage_effects():
    """Example 4: Coverage dependence of adsorption energy."""
    print("=" * 60)
    print("Example 4: Coverage Effects")
    print("=" * 60)

    # Define coverage via slab size
    coverages = [
        (2, 2, "1/4 ML"),  # One adsorbate on 2×2 slab = 1/4 ML
        (3, 3, "1/9 ML"),  # One adsorbate on 3×3 slab = 1/9 ML
    ]

    # Gas phase
    mol = molecule("O")
    mol.center(vacuum=10.0)
    mol.calc = EMT()
    E_mol = mol.get_potential_energy()

    results = {}

    print(f"\nTesting coverage dependence:")
    print(f"{'Coverage':<15} {'Slab Size':<15} {'E_ads (eV)':<15}")
    print("-" * 48)

    for nx, ny, coverage_label in coverages:
        # Clean slab
        slab = fcc111("Cu", size=(nx, ny, 4), vacuum=10.0)
        mask = [atom.position[2] < 6.0 for atom in slab]
        slab.set_constraint(FixAtoms(mask=mask))
        slab.calc = EMT()
        opt = BFGS(slab, logfile=None)
        opt.run(fmax=0.02)
        E_slab = slab.get_potential_energy()

        # Adsorbed
        slab_ads = slab.copy()
        add_adsorbate(slab_ads, mol, height=1.5, position="fcc")
        slab_ads.set_constraint(FixAtoms(mask=mask))
        slab_ads.calc = EMT()
        opt = BFGS(slab_ads, logfile=None)
        opt.run(fmax=0.02)

        E_total = slab_ads.get_potential_energy()
        E_ads = E_total - E_slab - E_mol

        results[coverage_label] = E_ads
        print(f"{coverage_label:<15} {nx}×{ny}×4{'':<8} {E_ads:<15.4f}")

    print(f"\nNote: Lower coverage generally gives stronger binding")
    print(f"(Less adsorbate-adsorbate repulsion)")
    print()


def example_5_molecular_adsorption():
    """Example 5: Molecular adsorption (CO)."""
    print("=" * 60)
    print("Example 5: Molecular Adsorption - CO")
    print("=" * 60)

    # Clean slab
    slab = fcc111("Pt", size=(3, 3, 4), vacuum=10.0)
    mask = [atom.position[2] < 6.0 for atom in slab]
    slab.set_constraint(FixAtoms(mask=mask))
    slab.calc = EMT()
    opt = BFGS(slab, logfile=None)
    opt.run(fmax=0.02)
    E_slab = slab.get_potential_energy()

    # CO molecule
    co = molecule("CO")
    co.center(vacuum=10.0)
    co.calc = EMT()
    opt = BFGS(co, logfile=None)
    opt.run(fmax=0.01)
    E_co = co.get_potential_energy()
    d_co = co.get_distance(0, 1)

    print(f"\nGas phase CO:")
    print(f"  Bond length: {d_co:.3f} Å")
    print(f"  Energy: {E_co:.4f} eV")

    # Test orientations
    orientations = ["vertical", "tilted"]

    for orientation in orientations:
        slab_ads = slab.copy()

        if orientation == "vertical":
            # C down, O up
            add_adsorbate(slab_ads, co, height=2.0, position="ontop")
        else:
            # Tilted
            co_tilted = co.copy()
            co_tilted.rotate(45, "y")
            add_adsorbate(slab_ads, co_tilted, height=2.0, position="ontop")

        slab_ads.set_constraint(FixAtoms(mask=mask))
        slab_ads.calc = EMT()
        opt = BFGS(slab_ads, logfile=None)
        opt.run(fmax=0.02)

        E_total = slab_ads.get_potential_energy()
        E_ads = E_total - E_slab - E_co

        # Get final CO bond length
        co_indices = list(range(len(slab), len(slab_ads)))
        d_co_ads = slab_ads.get_distance(co_indices[0], co_indices[1])

        print(f"\n{orientation.capitalize()} orientation:")
        print(f"  E_ads: {E_ads:.4f} eV")
        print(f"  CO bond length: {d_co_ads:.3f} Å")
        print(f"  Bond stretch: {d_co_ads - d_co:.3f} Å")


def example_6_dissociative_adsorption():
    """Example 6: Compare molecular vs dissociative adsorption."""
    print("=" * 60)
    print("Example 6: Molecular vs Dissociative Adsorption - O2")
    print("=" * 60)

    # Clean slab
    slab = fcc111("Cu", size=(3, 3, 4), vacuum=10.0)
    mask = [atom.position[2] < 6.0 for atom in slab]
    slab.set_constraint(FixAtoms(mask=mask))
    slab.calc = EMT()
    opt = BFGS(slab, logfile=None)
    opt.run(fmax=0.02)
    E_slab = slab.get_potential_energy()

    # O2 molecule
    o2 = molecule("O2")
    o2.center(vacuum=10.0)
    o2.calc = EMT()
    E_o2 = o2.get_potential_energy()

    # Single O atom (for dissociation)
    o_atom = molecule("O")
    o_atom.center(vacuum=10.0)
    o_atom.calc = EMT()
    E_o = o_atom.get_potential_energy()

    # 1. Molecular adsorption
    print(f"\n1. Molecular adsorption (O2):")
    slab_mol = slab.copy()
    add_adsorbate(slab_mol, o2, height=2.0, position="bridge")
    slab_mol.set_constraint(FixAtoms(mask=mask))
    slab_mol.calc = EMT()
    opt = BFGS(slab_mol, logfile=None)
    opt.run(fmax=0.02)

    E_mol_ads = slab_mol.get_potential_energy() - E_slab - E_o2
    print(f"   E_ads(O2) = {E_mol_ads:.4f} eV")

    # 2. Dissociative adsorption (2 O atoms)
    print(f"\n2. Dissociative adsorption (2×O):")
    slab_diss = slab.copy()
    add_adsorbate(slab_diss, o_atom, height=1.5, position="fcc")
    # Add second O atom at different site
    add_adsorbate(slab_diss, o_atom, height=1.5, position=(2.5, 2.5))

    slab_diss.set_constraint(FixAtoms(mask=mask))
    slab_diss.calc = EMT()
    opt = BFGS(slab_diss, logfile=None)
    opt.run(fmax=0.02)

    E_diss_ads = slab_diss.get_potential_energy() - E_slab - 2 * E_o
    print(f"   E_ads(2×O) = {E_diss_ads:.4f} eV")

    # Compare
    print(f"\n3. Comparison:")
    print(f"   Dissociation energy = E_ads(2×O) - E_ads(O2)")
    print(f"                       = {E_diss_ads:.4f} - {E_mol_ads:.4f}")
    print(f"                       = {E_diss_ads - E_mol_ads:.4f} eV")

    if E_diss_ads < E_mol_ads:
        print(f"   Dissociative adsorption preferred")
    else:
        print(f"   Molecular adsorption preferred")

    print()


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("Adsorption Energy Calculation Examples")
    print("=" * 60 + "\n")

    try:
        example_1_basic_adsorption()
        example_2_adsorption_sites()
        example_3_height_optimization()
        example_4_coverage_effects()
        example_5_molecular_adsorption()
        example_6_dissociative_adsorption()

        print("=" * 60)
        print("All examples completed successfully!")
        print("=" * 60)

    except Exception as e:
        print(f"\nError occurred: {e}")
        import traceback

        traceback.print_exc()
