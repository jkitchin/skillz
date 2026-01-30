"""
Elastic Constants Calculation Example

Demonstrates elastic constant calculations:
1. Basic elastic constant calculation
2. Full elastic tensor (Cij)
3. Derived mechanical properties
4. Sound velocities and Debye temperature
5. Convergence testing

Uses stress-strain relationships: σij = Cijkl εkl

Reference:
Nielsen & Martin, "First-principles calculation of stress,"
Phys. Rev. Lett. 50, 697 (1983)

Note: Requires elastic package for automated calculations
"""

from ase.build import bulk
from ase.calculators.emt import EMT
from ase.optimize import BFGS
from ase.constraints import ExpCellFilter
import numpy as np


def example_1_basic_elastic_moduli():
    """Example 1: Calculate basic elastic moduli from finite differences."""
    print("=" * 60)
    print("Example 1: Basic Elastic Moduli - Manual Calculation")
    print("=" * 60)

    # Create and optimize structure
    atoms = bulk("Cu", "fcc", a=3.6)
    atoms.calc = EMT()

    ecf = ExpCellFilter(atoms)
    opt = BFGS(ecf, logfile=None)
    opt.run(fmax=0.01)

    V0 = atoms.get_volume()
    a0 = atoms.cell.cellpar()[0]

    print(f"\nEquilibrium state:")
    print(f"  Lattice constant: {a0:.4f} Å")
    print(f"  Volume: {V0:.4f} ų")

    # Bulk modulus via volume deformation
    print(f"\nCalculating bulk modulus...")
    strains = np.linspace(-0.02, 0.02, 9)
    volumes = []
    energies = []

    for strain in strains:
        a_strained = atoms.copy()
        a_strained.set_cell(atoms.cell * (1 + strain), scale_atoms=True)
        a_strained.calc = EMT()

        E = a_strained.get_potential_energy()
        V = a_strained.get_volume()

        volumes.append(V)
        energies.append(E)

    # Fit E vs V to quadratic: E = E0 + 1/2 B V0 ((V-V0)/V0)^2
    # Or use second derivative at V0
    from scipy.optimize import curve_fit

    def energy_volume(V, E0, B):
        return E0 + 0.5 * B * V0 * ((V - V0) / V0) ** 2

    popt, _ = curve_fit(energy_volume, volumes, energies, p0=[min(energies), 100])
    E0_fit, B_fit = popt

    # Convert to GPa (eV/Å³ to GPa: multiply by 160.21766208)
    B_GPa = B_fit * 160.21766208

    print(f"\nBulk modulus:")
    print(f"  B = {B_fit:.4f} eV/ų")
    print(f"  B = {B_GPa:.2f} GPa")
    print(f"  (Experimental Cu: ~140 GPa)")
    print()


def example_2_stress_strain_method():
    """Example 2: Calculate elastic constants from stress-strain."""
    print("=" * 60)
    print("Example 2: Elastic Constants from Stress-Strain")
    print("=" * 60)

    # Equilibrium structure
    atoms = bulk("Cu", "fcc", a=3.6)
    atoms.calc = EMT()

    ecf = ExpCellFilter(atoms)
    opt = BFGS(ecf, logfile=None)
    opt.run(fmax=0.01)

    print(f"\nCalculating elastic constants for cubic system...")
    print(f"Cubic symmetry: C11, C12, C44")

    # For cubic: need 3 independent constants
    strain_magnitude = 0.01

    # Strain 1: Isotropic (εxx = εyy = εzz)
    atoms1 = atoms.copy()
    atoms1.set_cell(atoms.cell * (1 + strain_magnitude), scale_atoms=True)
    atoms1.calc = EMT()
    stress1 = atoms1.get_stress(voigt=True)  # [σxx, σyy, σzz, σyz, σxz, σxy]

    # Strain 2: Volume-conserving (εxx = -εyy)
    atoms2 = atoms.copy()
    cell2 = atoms.cell.copy()
    cell2[0] *= 1 + strain_magnitude
    cell2[1] *= 1 - strain_magnitude
    atoms2.set_cell(cell2, scale_atoms=True)
    atoms2.calc = EMT()
    stress2 = atoms2.get_stress(voigt=True)

    # Strain 3: Shear (εxy)
    atoms3 = atoms.copy()
    cell3 = atoms.cell.copy()
    cell3[0, 1] += cell3[0, 0] * strain_magnitude
    atoms3.set_cell(cell3, scale_atoms=True)
    atoms3.calc = EMT()
    stress3 = atoms3.get_stress(voigt=True)

    # Estimate elastic constants (simplified)
    # C11 from εxx: σxx = C11 εxx + C12 (εyy + εzz)
    # For isotropic strain: σxx = (C11 + 2C12) ε
    C_sum = -stress1[0] / strain_magnitude / 160.21766208  # GPa

    # For volume-conserving: σxx - σyy = 2(C11 - C12) εxx
    C_diff = -(stress2[0] - stress2[1]) / (2 * strain_magnitude) / 160.21766208

    C11 = (C_sum + 2 * C_diff) / 3
    C12 = (C_sum - C_diff) / 3

    # C44 from shear: σxy = 2 C44 εxy
    C44 = -stress3[5] / (2 * strain_magnitude) / 160.21766208

    print(f"\nElastic constants (Voigt notation):")
    print(f"  C11 = {C11:.2f} GPa")
    print(f"  C12 = {C12:.2f} GPa")
    print(f"  C44 = {C44:.2f} GPa")

    # Derived properties
    B = (C11 + 2 * C12) / 3  # Bulk modulus
    G = (C11 - C12 + 3 * C44) / 5  # Shear modulus (Voigt average)
    E = 9 * B * G / (3 * B + G)  # Young's modulus
    nu = (3 * B - 2 * G) / (2 * (3 * B + G))  # Poisson's ratio

    print(f"\nDerived properties:")
    print(f"  Bulk modulus (B): {B:.2f} GPa")
    print(f"  Shear modulus (G): {G:.2f} GPa")
    print(f"  Young's modulus (E): {E:.2f} GPa")
    print(f"  Poisson's ratio (ν): {nu:.3f}")
    print()


def example_3_with_elastic_package():
    """Example 3: Use elastic package for automated calculation."""
    try:
        from elastic import get_elastic_tensor

        ELASTIC_AVAILABLE = True
    except ImportError:
        ELASTIC_AVAILABLE = False

    if not ELASTIC_AVAILABLE:
        print("=" * 60)
        print("Example 3: Elastic Package (Not Installed)")
        print("=" * 60)
        print("\nSkipping: Install with: pip install elastic")
        print()
        return

    print("=" * 60)
    print("Example 3: Automated Elastic Tensor Calculation")
    print("=" * 60)

    # Create structure
    atoms = bulk("Al", "fcc", a=4.05)
    atoms.calc = EMT()

    # Optimize first
    ecf = ExpCellFilter(atoms)
    opt = BFGS(ecf, logfile=None)
    opt.run(fmax=0.01)

    print(f"\nCalculating full elastic tensor...")

    # Calculate elastic tensor
    # get_elastic_tensor applies small strains and computes stress response
    Cij, Bij = get_elastic_tensor(atoms, delta=0.01)

    print(f"\nElastic tensor (GPa):")
    print(f"  (Voigt notation: C_ij)")
    for i in range(6):
        row = "  "
        for j in range(6):
            row += f"{Cij[i, j]:8.2f}"
        print(row)

    # For cubic systems, check symmetry
    print(f"\nCubic elastic constants:")
    print(f"  C11 = {Cij[0, 0]:.2f} GPa")
    print(f"  C12 = {Cij[0, 1]:.2f} GPa")
    print(f"  C44 = {Cij[3, 3]:.2f} GPa")

    print(f"\nCompliance tensor (TPa⁻¹):")
    print(f"  (Voigt notation: S_ij)")
    for i in range(6):
        row = "  "
        for j in range(6):
            row += f"{Bij[i, j]:8.4f}"
        print(row)

    print()


def example_4_mechanical_properties():
    """Example 4: Calculate comprehensive mechanical properties."""
    print("=" * 60)
    print("Example 4: Comprehensive Mechanical Properties")
    print("=" * 60)

    # Simplified calculation for demonstration
    atoms = bulk("Cu", "fcc", a=3.6)
    atoms.calc = EMT()

    ecf = ExpCellFilter(atoms)
    opt = BFGS(ecf, logfile=None)
    opt.run(fmax=0.01)

    # Use approximate values for Cu (from EMT)
    C11 = 180.0  # GPa (approximate)
    C12 = 120.0
    C44 = 80.0

    print(f"\nElastic constants:")
    print(f"  C11 = {C11:.1f} GPa")
    print(f"  C12 = {C12:.1f} GPa")
    print(f"  C44 = {C44:.1f} GPa")

    # Mechanical stability criteria for cubic
    print(f"\nMechanical stability criteria:")
    stable1 = C11 - C12 > 0
    stable2 = C11 + 2 * C12 > 0
    stable3 = C44 > 0

    print(f"  C11 - C12 > 0: {stable1} ({C11 - C12:.1f} GPa)")
    print(f"  C11 + 2C12 > 0: {stable2} ({C11 + 2 * C12:.1f} GPa)")
    print(f"  C44 > 0: {stable3} ({C44:.1f} GPa)")
    print(f"  Mechanically stable: {all([stable1, stable2, stable3])}")

    # Aggregate moduli
    B_Voigt = (C11 + 2 * C12) / 3
    G_Voigt = (C11 - C12 + 3 * C44) / 5

    S11 = (C11 + C12) / ((C11 - C12) * (C11 + 2 * C12))
    S12 = -C12 / ((C11 - C12) * (C11 + 2 * C12))
    S44 = 1 / C44

    B_Reuss = 1 / (S11 + 2 * S12)
    G_Reuss = 5 / (4 * S11 - 4 * S12 + 3 * S44)

    # Voigt-Reuss-Hill average
    B = (B_Voigt + B_Reuss) / 2
    G = (G_Voigt + G_Reuss) / 2

    E = 9 * B * G / (3 * B + G)
    nu = (3 * B - 2 * G) / (2 * (3 * B + G))

    print(f"\nAggregate moduli:")
    print(f"  Bulk modulus (B):")
    print(f"    Voigt: {B_Voigt:.2f} GPa")
    print(f"    Reuss: {B_Reuss:.2f} GPa")
    print(f"    Hill:  {B:.2f} GPa")
    print(f"  Shear modulus (G):")
    print(f"    Voigt: {G_Voigt:.2f} GPa")
    print(f"    Reuss: {G_Reuss:.2f} GPa")
    print(f"    Hill:  {G:.2f} GPa")
    print(f"\n  Young's modulus (E): {E:.2f} GPa")
    print(f"  Poisson's ratio (ν): {nu:.3f}")

    # Hardness indicators
    B_over_G = B / G
    print(f"\n  B/G ratio: {B_over_G:.2f}")
    if B_over_G > 1.75:
        print(f"  Material: Ductile")
    else:
        print(f"  Material: Brittle")

    # Anisotropy
    A = 2 * C44 / (C11 - C12)
    print(f"\n  Zener anisotropy (A): {A:.3f}")
    if abs(A - 1.0) < 0.1:
        print(f"  Nearly isotropic")
    else:
        print(f"  Anisotropic")

    print()


def example_5_sound_velocities():
    """Example 5: Calculate sound velocities and Debye temperature."""
    print("=" * 60)
    print("Example 5: Sound Velocities and Debye Temperature")
    print("=" * 60)

    atoms = bulk("Cu", "fcc", a=3.6)
    atoms.calc = EMT()

    ecf = ExpCellFilter(atoms)
    opt = BFGS(ecf, logfile=None)
    opt.run(fmax=0.01)

    # Material properties
    C11 = 180.0  # GPa
    C12 = 120.0
    C44 = 80.0
    rho = 8960  # kg/m³ (density)

    print(f"\nMaterial properties:")
    print(f"  Density: {rho} kg/m³")
    print(f"  C11 = {C11:.1f} GPa")
    print(f"  C44 = {C44:.1f} GPa")

    # Convert GPa to Pa
    C11_Pa = C11 * 1e9
    C44_Pa = C44 * 1e9

    # Sound velocities (m/s)
    v_long = np.sqrt(C11_Pa / rho)  # Longitudinal
    v_trans = np.sqrt(C44_Pa / rho)  # Transverse/Shear

    print(f"\nSound velocities:")
    print(f"  Longitudinal (v_l): {v_long:.0f} m/s")
    print(f"  Transverse (v_t): {v_trans:.0f} m/s")

    # Average sound velocity
    v_m = ((1 / v_long**3 + 2 / v_trans**3) / 3) ** (-1 / 3)
    print(f"  Mean (v_m): {v_m:.0f} m/s")

    # Debye temperature
    # θ_D = (h/k_B) * v_m * (3N/4πV)^(1/3)
    h = 6.62607015e-34  # J·s
    k_B = 1.380649e-23  # J/K
    N_A = 6.02214076e23  # 1/mol

    V_atom = atoms.get_volume() * 1e-30  # m³
    n = len(atoms) / V_atom  # number density (1/m³)

    theta_D = (h / k_B) * v_m * (3 * n / (4 * np.pi)) ** (1 / 3)

    print(f"\nDebye temperature:")
    print(f"  θ_D = {theta_D:.0f} K")
    print(f"  (Experimental Cu: ~343 K)")
    print()


def example_6_strain_convergence():
    """Example 6: Test convergence with respect to strain magnitude."""
    print("=" * 60)
    print("Example 6: Strain Magnitude Convergence")
    print("=" * 60)

    atoms = bulk("Al", "fcc", a=4.05)
    atoms.calc = EMT()

    ecf = ExpCellFilter(atoms)
    opt = BFGS(ecf, logfile=None)
    opt.run(fmax=0.01)

    # Test different strain magnitudes
    strains = [0.001, 0.005, 0.01, 0.02, 0.03]
    bulk_moduli = []

    print(f"\nTesting strain magnitudes:")
    print(f"{'Strain':<12} {'B (GPa)':<15} {'Diff from prev':<15}")
    print("-" * 45)

    for delta in strains:
        # Volumetric strain
        volumes = []
        energies = []

        for eps in [-delta, 0, delta]:
            a_test = atoms.copy()
            a_test.set_cell(atoms.cell * (1 + eps), scale_atoms=True)
            a_test.calc = EMT()

            V = a_test.get_volume()
            E = a_test.get_potential_energy()

            volumes.append(V)
            energies.append(E)

        # Fit quadratic to get B
        V0 = volumes[1]
        from scipy.optimize import curve_fit

        def parabola(V, E0, B):
            return E0 + 0.5 * B * V0 * ((V - V0) / V0) ** 2

        popt, _ = curve_fit(parabola, volumes, energies)
        B = popt[1] * 160.21766208  # Convert to GPa

        bulk_moduli.append(B)

        if len(bulk_moduli) > 1:
            diff = B - bulk_moduli[-2]
            print(f"{delta:<12.4f} {B:<15.2f} {diff:+.2f}")
        else:
            print(f"{delta:<12.4f} {B:<15.2f} {'--':<15}")

    # Check convergence
    if len(bulk_moduli) > 2:
        converged = abs(bulk_moduli[-1] - bulk_moduli[-2]) < 1.0
        print(f"\nConverged (< 1 GPa): {converged}")
        print(f"Recommended strain: 0.005-0.01")

    print()


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("Elastic Constants Calculation Examples")
    print("=" * 60 + "\n")

    try:
        example_1_basic_elastic_moduli()
        example_2_stress_strain_method()
        example_3_with_elastic_package()
        example_4_mechanical_properties()
        example_5_sound_velocities()
        example_6_strain_convergence()

        print("=" * 60)
        print("All examples completed successfully!")
        print("=" * 60)

    except Exception as e:
        print(f"\nError occurred: {e}")
        import traceback

        traceback.print_exc()
