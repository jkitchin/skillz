"""
Complete Materials Properties Workflow

Demonstrates a comprehensive workflow combining multiple property calculations
for a realistic materials science study: heterogeneous catalysis on Cu surface.

Workflow:
1. Structure optimization and lattice constant determination
2. Elastic properties calculation
3. Surface energy for different facets
4. Adsorption energy for reactant
5. Reaction barrier via NEB
6. Summary and analysis

This example shows how different calculations build on each other
in a typical computational materials science project.
"""

from ase.build import bulk, fcc111, fcc100, molecule, add_adsorbate
from ase.optimize import BFGS
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms, ExpCellFilter
from ase.neb import NEB, NEBTools
from ase.eos import EquationOfState
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


class MaterialsStudy:
    """
    Complete materials property workflow for Cu catalyst.
    """

    def __init__(self, element="Cu", calculator=None):
        """Initialize study with element and calculator."""
        self.element = element
        self.calc = calculator if calculator else EMT()

        # Results storage
        self.results = {
            "element": element,
            "lattice_constant": None,
            "bulk_modulus": None,
            "elastic_constants": {},
            "surface_energies": {},
            "adsorption_energies": {},
            "reaction_barrier": None,
        }

        print("=" * 70)
        print(f"Materials Properties Study: {element}")
        print("=" * 70)
        print()

    def step1_optimize_bulk(self):
        """Step 1: Determine optimal lattice constant."""
        print("Step 1: Bulk Structure Optimization")
        print("-" * 70)

        # Scan lattice constants
        if self.element == "Cu":
            a_range = np.linspace(3.4, 3.8, 9)
        elif self.element == "Al":
            a_range = np.linspace(3.9, 4.2, 9)
        else:
            a_range = np.linspace(3.4, 4.2, 9)

        volumes = []
        energies = []

        print(f"Scanning lattice constants for {self.element}...")

        for a in a_range:
            atoms = bulk(self.element, "fcc", a=a)
            atoms.calc = self.calc
            E = atoms.get_potential_energy()
            V = atoms.get_volume()

            energies.append(E)
            volumes.append(V)

        # Fit equation of state
        eos = EquationOfState(volumes, energies, eos="birchmurnaghan")
        v0, e0, B = eos.fit()

        # Convert to lattice constant (FCC: V = a³/4)
        a0 = (4 * v0) ** (1 / 3)

        self.results["lattice_constant"] = a0
        self.results["bulk_modulus"] = B / 1e9  # Convert to GPa
        self.results["cohesive_energy"] = e0

        print(f"\nResults:")
        print(f"  Lattice constant: {a0:.4f} Å")
        print(f"  Bulk modulus: {B / 1e9:.2f} GPa")
        print(f"  Cohesive energy: {e0:.4f} eV")

        # Save plot
        eos.plot("eos_fit.png")
        print(f"  Saved: eos_fit.png")
        print()

        return a0

    def step2_elastic_properties(self, a0):
        """Step 2: Calculate elastic constants."""
        print("Step 2: Elastic Properties")
        print("-" * 70)

        atoms = bulk(self.element, "fcc", a=a0)
        atoms.calc = self.calc

        # Optimize with cell relaxation
        ecf = ExpCellFilter(atoms)
        opt = BFGS(ecf, logfile=None)
        opt.run(fmax=0.01)

        # Calculate bulk modulus from volume variations
        strain_magnitude = 0.01
        volumes = []
        energies = []

        for eps in [-strain_magnitude, 0, strain_magnitude]:
            a_test = atoms.copy()
            a_test.set_cell(atoms.cell * (1 + eps), scale_atoms=True)
            a_test.calc = self.calc
            E = a_test.get_potential_energy()
            V = a_test.get_volume()

            volumes.append(V)
            energies.append(E)

        # Fit quadratic
        V0 = volumes[1]
        from scipy.optimize import curve_fit

        def parabola(V, E0, B):
            return E0 + 0.5 * B * V0 * ((V - V0) / V0) ** 2

        popt, _ = curve_fit(parabola, volumes, energies)
        B_calc = popt[1] * 160.21766208  # eV/Å³ to GPa

        self.results["elastic_constants"]["bulk_modulus"] = B_calc

        print(f"\nResults:")
        print(f"  Bulk modulus: {B_calc:.2f} GPa")
        print()

        return atoms

    def step3_surface_energies(self, a0):
        """Step 3: Calculate surface energies for different facets."""
        print("Step 3: Surface Energies")
        print("-" * 70)

        # Bulk energy per atom
        bulk_atoms = bulk(self.element, "fcc", a=a0)
        bulk_atoms.calc = self.calc
        E_bulk_per_atom = bulk_atoms.get_potential_energy() / len(bulk_atoms)

        print(f"Bulk energy: {E_bulk_per_atom:.4f} eV/atom")

        # Calculate for different facets
        facets = [
            ("111", fcc111(self.element, size=(3, 3, 7), vacuum=10, a=a0)),
            ("100", fcc100(self.element, size=(3, 3, 7), vacuum=10, a=a0)),
        ]

        print(f"\nCalculating surface energies...")

        for facet_name, slab in facets:
            # Fix bottom layers
            z_positions = slab.get_positions()[:, 2]
            z_min = z_positions.min()
            mask = [z < z_min + 4.0 for z in z_positions]
            slab.set_constraint(FixAtoms(mask=mask))

            slab.calc = self.calc
            opt = BFGS(slab, logfile=None)
            opt.run(fmax=0.02)

            E_slab = slab.get_potential_energy()
            N = len(slab)

            # Surface area
            cell = slab.get_cell()
            A = np.linalg.norm(np.cross(cell[0], cell[1]))

            # Surface energy
            gamma = (E_slab - N * E_bulk_per_atom) / (2 * A)
            gamma_mJ = gamma * 16.0217  # Convert to mJ/m²

            self.results["surface_energies"][facet_name] = gamma_mJ

            print(f"  {self.element}({facet_name}): {gamma_mJ:.2f} mJ/m²")

        # Most stable facet
        most_stable = min(
            self.results["surface_energies"], key=self.results["surface_energies"].get
        )
        print(f"\nMost stable facet: ({most_stable})")
        print()

        return most_stable, a0

    def step4_adsorption_study(self, facet, a0):
        """Step 4: Study CO adsorption on most stable surface."""
        print("Step 4: Adsorption Study - CO on Cu")
        print("-" * 70)

        # Clean slab
        if facet == "111":
            slab = fcc111(self.element, size=(3, 3, 4), vacuum=10, a=a0)
        else:
            slab = fcc100(self.element, size=(3, 3, 4), vacuum=10, a=a0)

        z_positions = slab.get_positions()[:, 2]
        z_min = z_positions.min()
        mask = [z < z_min + 3.0 for z in z_positions]
        slab.set_constraint(FixAtoms(mask=mask))

        slab.calc = self.calc
        opt = BFGS(slab, logfile=None)
        opt.run(fmax=0.02)

        E_slab = slab.get_potential_energy()

        # CO molecule
        co = molecule("CO")
        co.center(vacuum=10.0)
        co.calc = self.calc
        opt_co = BFGS(co, logfile=None)
        opt_co.run(fmax=0.01)

        E_co = co.get_potential_energy()

        print(f"Clean slab: {E_slab:.4f} eV")
        print(f"CO molecule: {E_co:.4f} eV")

        # Test different sites
        sites = ["ontop", "bridge", "fcc"]

        print(f"\nAdsorption at different sites:")

        for site in sites:
            slab_ads = slab.copy()
            add_adsorbate(slab_ads, co, height=2.0, position=site)
            slab_ads.set_constraint(FixAtoms(mask=mask))

            slab_ads.calc = self.calc
            opt = BFGS(slab_ads, logfile=None)
            opt.run(fmax=0.02)

            E_total = slab_ads.get_potential_energy()
            E_ads = E_total - E_slab - E_co

            self.results["adsorption_energies"][site] = E_ads

            print(f"  {site:>10}: E_ads = {E_ads:+.4f} eV")

        # Best site
        best_site = min(
            self.results["adsorption_energies"], key=self.results["adsorption_energies"].get
        )
        print(f"\nPreferred site: {best_site}")
        print()

        return slab, co, best_site, a0

    def step5_reaction_barrier(self, slab_template, co, site, a0):
        """Step 5: Calculate CO diffusion barrier via NEB."""
        print("Step 5: Reaction Barrier - CO Surface Diffusion")
        print("-" * 70)

        # Initial state: CO at site 1
        initial = slab_template.copy()
        add_adsorbate(initial, co, height=2.0, position=site)

        z_positions = initial.get_positions()[:, 2]
        z_min = z_positions.min()
        mask = [z < z_min + 3.0 for z in z_positions]
        initial.set_constraint(FixAtoms(mask=mask))

        initial.calc = self.calc
        opt = BFGS(initial, logfile=None)
        opt.run(fmax=0.02)

        # Final state: CO at neighboring site
        final = slab_template.copy()

        # Move to adjacent site
        cell = slab_template.get_cell()
        if site == "fcc":
            neighbor_offset = cell[0] / 3  # Move to nearby fcc
        else:
            neighbor_offset = cell[0] / 2

        add_adsorbate(final, co, height=2.0, position="ontop")
        # Manually adjust position
        final.positions[-2:] += neighbor_offset

        final.set_constraint(FixAtoms(mask=mask))
        final.calc = self.calc
        opt = BFGS(final, logfile=None)
        opt.run(fmax=0.02)

        print(f"Initial state: {initial.get_potential_energy():.4f} eV")
        print(f"Final state: {final.get_potential_energy():.4f} eV")

        # NEB calculation
        print(f"\nRunning NEB with 5 images...")

        images = [initial]
        for i in range(3):
            images.append(initial.copy())
        images.append(final)

        neb = NEB(images, climb=True)
        neb.interpolate()

        opt_neb = BFGS(neb, trajectory="neb.traj", logfile=None)
        opt_neb.run(fmax=0.05)

        # Analyze barrier
        nebtools = NEBTools(images)

        Ef_forward, dE = nebtools.get_barrier()

        self.results["reaction_barrier"] = Ef_forward

        print(f"\nResults:")
        print(f"  Forward barrier: {Ef_forward:.4f} eV")
        print(f"  Reverse barrier: {Ef_forward - dE:.4f} eV")

        # Plot
        try:
            fig = nebtools.plot_band()
            fig.savefig("neb_barrier.png", dpi=150)
            print(f"  Saved: neb_barrier.png")
        except:
            print(f"  (Band plot skipped)")

        print()

    def step6_summary(self):
        """Step 6: Generate summary report."""
        print("Step 6: Summary Report")
        print("=" * 70)

        print(f"\nMaterial: {self.element}")
        print(f"\n1. Bulk Properties:")
        print(f"   Lattice constant: {self.results['lattice_constant']:.4f} Å")
        print(f"   Bulk modulus: {self.results['bulk_modulus']:.2f} GPa")

        print(f"\n2. Surface Energies:")
        for facet, gamma in self.results["surface_energies"].items():
            print(f"   ({facet}): {gamma:.2f} mJ/m²")

        print(f"\n3. CO Adsorption Energies:")
        for site, E_ads in self.results["adsorption_energies"].items():
            print(f"   {site}: {E_ads:+.4f} eV")

        print(f"\n4. Reaction Barrier:")
        print(f"   CO diffusion: {self.results['reaction_barrier']:.4f} eV")

        print(f"\n{'=' * 70}")
        print(f"Study Complete!")
        print(f"{'=' * 70}\n")

        # Create summary plot
        self._create_summary_plot()

    def _create_summary_plot(self):
        """Create comprehensive summary visualization."""
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))

        # Surface energies
        facets = list(self.results["surface_energies"].keys())
        gamma_values = list(self.results["surface_energies"].values())

        ax1.bar(facets, gamma_values, color="steelblue", alpha=0.7)
        ax1.set_ylabel("Surface Energy (mJ/m²)", fontsize=12)
        ax1.set_title(f"{self.element} Surface Energies", fontsize=14)
        ax1.grid(True, alpha=0.3)

        # Adsorption energies
        sites = list(self.results["adsorption_energies"].keys())
        E_ads_values = list(self.results["adsorption_energies"].values())

        colors = ["green" if E < 0 else "red" for E in E_ads_values]
        ax2.bar(sites, E_ads_values, color=colors, alpha=0.7)
        ax2.axhline(y=0, color="black", linestyle="--", linewidth=1)
        ax2.set_ylabel("Adsorption Energy (eV)", fontsize=12)
        ax2.set_title("CO Adsorption Sites", fontsize=14)
        ax2.grid(True, alpha=0.3)

        # Property summary (text)
        ax3.axis("off")
        summary_text = f"Material: {self.element}\n\n"
        summary_text += f"Lattice constant:\n  {self.results['lattice_constant']:.4f} Å\n\n"
        summary_text += f"Bulk modulus:\n  {self.results['bulk_modulus']:.2f} GPa\n\n"
        summary_text += f"CO diffusion barrier:\n  {self.results['reaction_barrier']:.4f} eV\n\n"
        summary_text += f"Most stable surface:\n  ({facets[np.argmin(gamma_values)]})\n\n"
        summary_text += f"Best adsorption site:\n  {sites[np.argmin(E_ads_values)]}"

        ax3.text(
            0.1, 0.5, summary_text, fontsize=12, family="monospace", verticalalignment="center"
        )

        # Workflow diagram (simple)
        ax4.axis("off")
        workflow = [
            "Workflow Steps:",
            "",
            "1. Bulk optimization",
            "   ↓",
            "2. Elastic properties",
            "   ↓",
            "3. Surface energies",
            "   ↓",
            "4. Adsorption study",
            "   ↓",
            "5. Reaction barriers",
            "   ↓",
            "6. Analysis",
        ]
        ax4.text(
            0.2,
            0.5,
            "\n".join(workflow),
            fontsize=11,
            family="monospace",
            verticalalignment="center",
        )

        plt.tight_layout()
        plt.savefig("workflow_summary.png", dpi=150)
        print(f"Saved: workflow_summary.png")

    def run_complete_workflow(self):
        """Execute complete workflow."""
        a0 = self.step1_optimize_bulk()
        bulk_atoms = self.step2_elastic_properties(a0)
        facet, a0 = self.step3_surface_energies(a0)
        slab, co, site, a0 = self.step4_adsorption_study(facet, a0)
        self.step5_reaction_barrier(slab, co, site, a0)
        self.step6_summary()


if __name__ == "__main__":
    print("\n" + "=" * 70)
    print("Complete Materials Properties Workflow Example")
    print("=" * 70 + "\n")

    print("This workflow demonstrates a typical computational study")
    print("for heterogeneous catalysis:\n")
    print("  - Optimize bulk structure")
    print("  - Calculate mechanical properties")
    print("  - Determine surface stability")
    print("  - Study molecular adsorption")
    print("  - Calculate reaction barriers")
    print()

    try:
        # Run study for Cu
        study = MaterialsStudy(element="Cu", calculator=EMT())
        study.run_complete_workflow()

        print("\nWorkflow completed successfully!")
        print("\nGenerated files:")
        print("  - eos_fit.png: Equation of state fit")
        print("  - neb_barrier.png: Reaction barrier profile")
        print("  - workflow_summary.png: Summary visualization")
        print("  - neb.traj: NEB trajectory for visualization")

    except Exception as e:
        print(f"\nError occurred: {e}")
        import traceback

        traceback.print_exc()
