"""
Defect Formation Energy Calculations

Demonstrates calculation of various defect formation energies:
1. Vacancy formation energy
2. Interstitial defects
3. Substitutional defects
4. Stacking fault energy
5. Peierls barrier for dislocation motion
6. Solute-dislocation interaction energy

Key formulas:
- Vacancy: E_f^v = E(N-1) - (N-1)/N × E(N) - μ
- Interstitial: E_f^i = E(N+1) - N/(N-1) × E(N) + μ
- Stacking fault: γ_sf = (E_faulted - E_perfect) / A

References:
Freysoldt et al., "First-principles calculations for point defects in solids,"
Rev. Mod. Phys. 86, 253 (2014)

Vitek, "Intrinsic stacking faults in body-centred cubic crystals,"
Phil. Mag. 18, 773 (1968)

Olmsted et al., "Atomistic simulations of dislocation mobility in Al, Ni and Al/Mg alloys,"
Model. Simul. Mater. Sci. Eng. 13, 371 (2005)
"""

from ase.build import bulk, fcc111, add_adsorbate, make_supercell
from ase.optimize import BFGS
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.neb import NEB
import numpy as np


def example_1_vacancy_formation():
    """Example 1: Calculate vacancy formation energy."""
    print("=" * 60)
    print("Example 1: Vacancy Formation Energy")
    print("=" * 60)

    # Perfect supercell
    atoms_perfect = bulk("Cu", "fcc", a=3.6)
    supercell = make_supercell(atoms_perfect, [[4, 0, 0], [0, 4, 0], [0, 0, 4]])
    supercell.calc = EMT()

    opt = BFGS(supercell, logfile=None)
    opt.run(fmax=0.02)

    E_perfect = supercell.get_potential_energy()
    N = len(supercell)

    print(f"\nPerfect supercell:")
    print(f"  Atoms: {N}")
    print(f"  Energy: {E_perfect:.4f} eV")
    print(f"  Energy per atom: {E_perfect / N:.4f} eV")

    # Create vacancy (remove central atom)
    atoms_vacancy = supercell.copy()
    center_idx = N // 2
    del atoms_vacancy[center_idx]

    atoms_vacancy.calc = EMT()
    opt_vac = BFGS(atoms_vacancy, logfile=None)
    opt_vac.run(fmax=0.02)

    E_vacancy = atoms_vacancy.get_potential_energy()

    print(f"\nVacancy supercell:")
    print(f"  Atoms: {len(atoms_vacancy)}")
    print(f"  Energy: {E_vacancy:.4f} eV")

    # Vacancy formation energy
    # E_f = E(N-1) - (N-1)/N * E(N)
    E_form_vacancy = E_vacancy - ((N - 1) / N) * E_perfect

    print(f"\nVacancy formation energy:")
    print(f"  E_f^v = E(N-1) - (N-1)/N × E(N)")
    print(f"  E_f^v = {E_vacancy:.4f} - {(N - 1) / N:.4f} × {E_perfect:.4f}")
    print(f"  E_f^v = {E_form_vacancy:.4f} eV")
    print(f"  (Experimental Cu: ~1.3 eV)")

    # Analyze local relaxation
    perfect_positions = supercell.get_positions()
    vacancy_positions = atoms_vacancy.get_positions()

    # Compare neighbors of vacancy site
    vacancy_site = supercell.get_positions()[center_idx]
    neighbor_relaxations = []

    for i, pos_vac in enumerate(vacancy_positions):
        # Find corresponding atom in perfect cell
        for j, pos_perf in enumerate(perfect_positions):
            if j != center_idx and np.linalg.norm(pos_vac - pos_perf) < 0.5:
                dist_to_vacancy = np.linalg.norm(pos_perf - vacancy_site)
                displacement = np.linalg.norm(pos_vac - pos_perf)
                if dist_to_vacancy < 5.0:
                    neighbor_relaxations.append((dist_to_vacancy, displacement))
                break

    if neighbor_relaxations:
        neighbor_relaxations.sort()
        print(f"\nNeighbor relaxations:")
        for i, (dist, disp) in enumerate(neighbor_relaxations[:5]):
            print(f"  Shell {i + 1} (r={dist:.2f} Å): Δ={disp:.3f} Å")

    print()


def example_2_interstitial_defect():
    """Example 2: Calculate interstitial formation energy."""
    print("=" * 60)
    print("Example 2: Interstitial Defect Formation")
    print("=" * 60)

    # Perfect supercell
    atoms_perfect = bulk("Cu", "fcc", a=3.6)
    supercell = make_supercell(atoms_perfect, [[3, 0, 0], [0, 3, 0], [0, 0, 3]])
    supercell.calc = EMT()

    opt = BFGS(supercell, logfile=None)
    opt.run(fmax=0.02)

    E_perfect = supercell.get_potential_energy()
    N = len(supercell)

    print(f"\nPerfect supercell:")
    print(f"  Atoms: {N}")
    print(f"  Energy: {E_perfect:.4f} eV/atom")

    # Add interstitial at octahedral site
    # For FCC, octahedral site is at (a/2, a/2, 0) in conventional cell
    cell = supercell.get_cell()
    interstitial_site = (cell[0] + cell[1]) / 2

    atoms_interstitial = supercell.copy()
    from ase import Atom

    atoms_interstitial.append(Atom("Cu", position=interstitial_site))

    atoms_interstitial.calc = EMT()
    opt_int = BFGS(atoms_interstitial, logfile=None)
    opt_int.run(fmax=0.02)

    E_interstitial = atoms_interstitial.get_potential_energy()

    print(f"\nInterstitial supercell:")
    print(f"  Atoms: {len(atoms_interstitial)}")
    print(f"  Energy: {E_interstitial:.4f} eV")

    # Interstitial formation energy
    # E_f^i = E(N+1) - (N+1)/N * E(N)
    E_form_interstitial = E_interstitial - ((N + 1) / N) * E_perfect

    print(f"\nInterstitial formation energy:")
    print(f"  E_f^i = E(N+1) - (N+1)/N × E(N)")
    print(f"  E_f^i = {E_interstitial:.4f} - {(N + 1) / N:.4f} × {E_perfect:.4f}")
    print(f"  E_f^i = {E_form_interstitial:.4f} eV")
    print(f"  (Typical for FCC: 3-4 eV)")

    # Final interstitial position
    int_pos_final = atoms_interstitial.get_positions()[-1]
    displacement = np.linalg.norm(int_pos_final - interstitial_site)

    print(f"\nInterstitial relaxation:")
    print(f"  Initial site: octahedral")
    print(f"  Displacement: {displacement:.3f} Å")

    print()


def example_3_substitutional_defect():
    """Example 3: Calculate substitutional defect energy."""
    print("=" * 60)
    print("Example 3: Substitutional Defect - Au in Cu")
    print("=" * 60)

    # Host material
    atoms_host = bulk("Cu", "fcc", a=3.6)
    supercell_host = make_supercell(atoms_host, [[3, 0, 0], [0, 3, 0], [0, 0, 3]])
    supercell_host.calc = EMT()

    opt = BFGS(supercell_host, logfile=None)
    opt.run(fmax=0.02)

    E_host = supercell_host.get_potential_energy()
    N = len(supercell_host)

    print(f"\nPure Cu supercell:")
    print(f"  Atoms: {N}")
    print(f"  Energy: {E_host:.4f} eV")

    # Substitutional defect (replace central Cu with Au)
    atoms_subst = supercell_host.copy()
    center_idx = N // 2
    atoms_subst[center_idx].symbol = "Au"

    atoms_subst.calc = EMT()
    opt_subst = BFGS(atoms_subst, logfile=None)
    opt_subst.run(fmax=0.02)

    E_subst = atoms_subst.get_potential_energy()

    print(f"\nCu with Au substitution:")
    print(f"  Formula: {atoms_subst.get_chemical_formula()}")
    print(f"  Energy: {E_subst:.4f} eV")

    # Chemical potentials (bulk references)
    cu_bulk = bulk("Cu", "fcc", a=3.6)
    cu_bulk.calc = EMT()
    mu_Cu = cu_bulk.get_potential_energy() / len(cu_bulk)

    au_bulk = bulk("Au", "fcc", a=4.08)
    au_bulk.calc = EMT()
    mu_Au = au_bulk.get_potential_energy() / len(au_bulk)

    print(f"\nChemical potentials:")
    print(f"  μ(Cu) = {mu_Cu:.4f} eV")
    print(f"  μ(Au) = {mu_Au:.4f} eV")

    # Substitutional formation energy
    # E_f = E(subst) - E(host) - μ(Au) + μ(Cu)
    E_form_subst = E_subst - E_host - mu_Au + mu_Cu

    print(f"\nSubstitutional formation energy:")
    print(f"  E_f = E(Cu_{N - 1}Au) - E(Cu_N) - μ(Au) + μ(Cu)")
    print(f"  E_f = {E_form_subst:.4f} eV")

    if E_form_subst < 0:
        print(f"  Au substitution in Cu is favorable")
    else:
        print(f"  Au substitution in Cu is unfavorable")

    print()


def example_4_stacking_fault_energy():
    """Example 4: Calculate stacking fault energy."""
    print("=" * 60)
    print("Example 4: Stacking Fault Energy")
    print("=" * 60)

    # Perfect FCC stacking: ...ABCABC...
    # Intrinsic stacking fault: ...ABCABABC... (remove C layer)

    # Create perfect slab
    slab_perfect = fcc111("Cu", size=(3, 3, 12), vacuum=0)
    slab_perfect.calc = EMT()

    opt = BFGS(slab_perfect, logfile=None)
    opt.run(fmax=0.02)

    E_perfect = slab_perfect.get_potential_energy()
    N_perfect = len(slab_perfect)

    print(f"\nPerfect FCC slab:")
    print(f"  Atoms: {N_perfect}")
    print(f"  Energy: {E_perfect:.4f} eV")

    # Create stacking fault
    # Remove one layer and shift upper half
    slab_fault = slab_perfect.copy()

    # Identify layers by z-coordinate
    z_positions = slab_fault.get_positions()[:, 2]
    unique_z = np.unique(np.round(z_positions, 2))
    z_mid = (unique_z.max() + unique_z.min()) / 2

    # Shift upper half to create intrinsic stacking fault
    # Shift by 1/3 of lattice vector in (111) plane
    positions = slab_fault.get_positions()
    cell = slab_fault.get_cell()

    # Shift vector: a/3 * [110] for intrinsic SF
    shift = (cell[0] - cell[1]) / 3

    for i, pos in enumerate(positions):
        if pos[2] > z_mid:
            slab_fault.positions[i] += shift

    slab_fault.calc = EMT()
    opt_fault = BFGS(slab_fault, logfile=None)
    opt_fault.run(fmax=0.02)

    E_fault = slab_fault.get_potential_energy()

    print(f"\nFaulted slab:")
    print(f"  Energy: {E_fault:.4f} eV")

    # Stacking fault energy per unit area
    cell_fault = slab_fault.get_cell()
    A = np.linalg.norm(np.cross(cell_fault[0], cell_fault[1]))

    gamma_sf = (E_fault - E_perfect) / A

    print(f"\nStacking fault energy:")
    print(f"  γ_sf = (E_fault - E_perfect) / A")
    print(f"  γ_sf = ({E_fault:.4f} - {E_perfect:.4f}) / {A:.3f}")
    print(f"  γ_sf = {gamma_sf:.6f} eV/Å²")
    print(f"  γ_sf = {gamma_sf * 1000:.2f} meV/Å²")
    print(f"  γ_sf = {gamma_sf * 16.0217:.2f} mJ/m²")
    print(f"  (Experimental Cu: ~40 mJ/m²)")

    print()


def example_5_peierls_barrier():
    """Example 5: Calculate Peierls barrier for dislocation motion."""
    print("=" * 60)
    print("Example 5: Peierls Barrier for Dislocation")
    print("=" * 60)

    print(f"\nPeierls barrier: Energy barrier for dislocation glide")
    print(f"between equivalent positions in the lattice")

    # Simplified 2D model for demonstration
    # Create slab with edge dislocation

    # Perfect slab as reference
    slab = fcc111("Al", size=(6, 6, 8), vacuum=10.0)

    # Fix bottom layers
    z_positions = slab.get_positions()[:, 2]
    z_min = z_positions.min()
    mask = [z < z_min + 4.0 for z in z_positions]
    slab.set_constraint(FixAtoms(mask=mask))

    slab.calc = EMT()
    opt = BFGS(slab, logfile=None)
    opt.run(fmax=0.02)

    E_ref = slab.get_potential_energy()

    print(f"\nReference slab:")
    print(f"  Atoms: {len(slab)}")
    print(f"  Energy: {E_ref:.4f} eV")

    # Simplified dislocation: displace atoms to simulate dislocation core
    # In reality, use atomistic dislocation tools (e.g., atomman, matscipy)

    # Create displaced configurations
    n_images = 5
    images = []

    for i in range(n_images):
        image = slab.copy()
        positions = image.get_positions()

        # Apply displacement field for dislocation motion
        # Simplified: shift top layers progressively
        fraction = i / (n_images - 1)
        shift = fraction * (image.cell[0] / 2)

        for j, pos in enumerate(positions):
            if not mask[j]:  # Only move unfixed atoms
                # Smooth displacement
                z_frac = (pos[2] - z_min - 4.0) / (z_positions.max() - z_min - 4.0)
                if z_frac > 0:
                    displacement = shift * (1 - np.exp(-3 * z_frac))
                    image.positions[j, 0] += displacement[0]
                    image.positions[j, 1] += displacement[1]

        image.set_constraint(FixAtoms(mask=mask))
        image.calc = EMT()
        images.append(image)

    # Calculate energies
    print(f"\nDislocation displacement path:")
    print(f"{'Step':<8} {'Displacement':<15} {'Energy (eV)':<15} {'Barrier (meV)':<15}")
    print("-" * 60)

    energies = []
    for i, image in enumerate(images):
        E = image.get_potential_energy()
        energies.append(E)

        displacement_fraction = i / (n_images - 1)
        barrier = (E - E_ref) * 1000  # meV

        print(f"{i:<8} {displacement_fraction:<15.3f} {E:<15.4f} {barrier:<15.2f}")

    # Peierls barrier
    E_barrier = max(energies) - min(energies)

    print(f"\nPeierls barrier:")
    print(f"  ΔE = {E_barrier:.4f} eV")
    print(f"  ΔE = {E_barrier * 1000:.2f} meV")
    print(f"  (Typical for Al: 10-50 meV)")

    print(f"\nNote: This is a simplified model. For accurate calculations,")
    print(f"use specialized dislocation tools (atomman, matscipy) with")
    print(f"anisotropic elasticity theory and proper boundary conditions.")

    print()


def example_6_solute_dislocation_interaction():
    """Example 6: Calculate solute-dislocation interaction energy."""
    print("=" * 60)
    print("Example 6: Solute-Screw Dislocation Interaction")
    print("=" * 60)

    print(f"\nSolute-dislocation interaction: Energy change when a")
    print(f"solute atom is placed near a dislocation core")

    # Reference: Bulk with solute
    bulk_host = bulk("Al", "fcc", a=4.05)
    supercell_bulk = make_supercell(bulk_host, [[4, 0, 0], [0, 4, 0], [0, 0, 4]])
    supercell_bulk.calc = EMT()

    opt = BFGS(supercell_bulk, logfile=None)
    opt.run(fmax=0.02)

    E_bulk_host = supercell_bulk.get_potential_energy()

    # Add Mg solute to bulk
    supercell_bulk_solute = supercell_bulk.copy()
    center = len(supercell_bulk_solute) // 2
    supercell_bulk_solute[center].symbol = "Mg"

    supercell_bulk_solute.calc = EMT()
    opt = BFGS(supercell_bulk_solute, logfile=None)
    opt.run(fmax=0.02)

    E_bulk_solute = supercell_bulk_solute.get_potential_energy()

    print(f"\nBulk reference:")
    print(f"  Pure Al: {E_bulk_host:.4f} eV")
    print(f"  Al + Mg: {E_bulk_solute:.4f} eV")

    # Simplified dislocation model
    # In practice, create dislocation using cylindrical geometry
    # Here: approximate with strained slab

    slab = fcc111("Al", size=(5, 5, 8), vacuum=10.0)

    # Apply shear strain to simulate screw dislocation stress field
    cell = slab.get_cell()
    cell[0, 1] += 0.5  # Shear
    slab.set_cell(cell, scale_atoms=True)

    # Fix edges
    positions = slab.get_positions()
    x_coords = positions[:, 0]
    x_min, x_max = x_coords.min(), x_coords.max()

    mask = [(x < x_min + 2.0) or (x > x_max - 2.0) for x in x_coords]
    slab.set_constraint(FixAtoms(mask=mask))

    slab.calc = EMT()
    opt = BFGS(slab, logfile=None)
    opt.run(fmax=0.02)

    E_disloc = slab.get_potential_energy()

    # Add solute near dislocation core
    slab_solute = slab.copy()

    # Find atom near center (dislocation core)
    center_pos = slab.get_cell().sum(axis=0) / 2
    distances = [np.linalg.norm(pos - center_pos) for pos in positions]
    core_idx = np.argmin(distances)

    slab_solute[core_idx].symbol = "Mg"
    slab_solute.set_constraint(FixAtoms(mask=mask))

    slab_solute.calc = EMT()
    opt = BFGS(slab_solute, logfile=None)
    opt.run(fmax=0.02)

    E_disloc_solute = slab_solute.get_potential_energy()

    print(f"\nDislocation model (strained slab):")
    print(f"  Pure Al: {E_disloc:.4f} eV")
    print(f"  Al + Mg: {E_disloc_solute:.4f} eV")

    # Interaction energy
    # E_int = [E(disloc+solute) - E(disloc)] - [E(bulk+solute) - E(bulk)]
    E_int = (E_disloc_solute - E_disloc) - (E_bulk_solute - E_bulk_host)

    print(f"\nSolute-dislocation interaction energy:")
    print(f"  E_int = [E(d+s) - E(d)] - [E(b+s) - E(b)]")
    print(f"  E_int = [{E_disloc_solute:.4f} - {E_disloc:.4f}] - ")
    print(f"          [{E_bulk_solute:.4f} - {E_bulk_host:.4f}]")
    print(f"  E_int = {E_int:.4f} eV")
    print(f"  E_int = {E_int * 1000:.2f} meV")

    if E_int < 0:
        print(f"\nSolute is attracted to dislocation (solute atmosphere)")
    else:
        print(f"\nSolute is repelled from dislocation")

    print(f"\nNote: This is a simplified model. Accurate calculations require:")
    print(f"  - Proper dislocation geometry (Volterra fields)")
    print(f"  - Large supercells (>1000 atoms)")
    print(f"  - Flexible boundary conditions")
    print(f"  - Tools: matscipy, atomman")

    print()


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("Defect Formation Energy Examples")
    print("=" * 60 + "\n")

    try:
        example_1_vacancy_formation()
        example_2_interstitial_defect()
        example_3_substitutional_defect()
        example_4_stacking_fault_energy()
        example_5_peierls_barrier()
        example_6_solute_dislocation_interaction()

        print("=" * 60)
        print("All examples completed successfully!")
        print("=" * 60)
        print("\nNote: Examples 5 and 6 use simplified models.")
        print("For production calculations of dislocations, use:")
        print("  - matscipy: https://github.com/libAtoms/matscipy")
        print("  - atomman: https://github.com/usnistgov/atomman")

    except Exception as e:
        print(f"\nError occurred: {e}")
        import traceback

        traceback.print_exc()
