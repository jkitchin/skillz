"""
Materials Project Search Examples

Demonstrates various search patterns and queries using the Materials Project API.
"""

from mp_api.client import MPRester

# Set your API key as environment variable: export MP_API_KEY="your_key"
# Or pass directly: MPRester(api_key="your_key")


def example_1_basic_formula_search():
    """Search by chemical formula."""
    print("Example 1: Basic Formula Search")
    print("=" * 50)

    with MPRester() as mpr:
        # Simple search
        docs = mpr.materials.summary.search(
            formula="GaN", fields=["material_id", "formula_pretty", "band_gap", "energy_per_atom"]
        )

        print(f"Found {len(docs)} GaN structures:")
        for doc in docs:
            print(f"  {doc.material_id}: {doc.formula_pretty}")
            print(f"    Band gap: {doc.band_gap} eV")
            print(f"    Energy: {doc.energy_per_atom:.3f} eV/atom")
        print()


def example_2_property_screening():
    """Screen materials by multiple properties."""
    print("Example 2: Property Screening - Stable Semiconductors")
    print("=" * 50)

    with MPRester() as mpr:
        # Find stable semiconductors in the Ga-N system
        docs = mpr.materials.summary.search(
            elements=["Ga", "N"],
            energy_above_hull=(0, 0.05),  # Stable or nearly stable
            band_gap=(1.0, 4.0),  # Semiconducting range
            fields=[
                "material_id",
                "formula_pretty",
                "band_gap",
                "is_gap_direct",
                "energy_above_hull",
                "formation_energy_per_atom",
            ],
        )

        print(f"Found {len(docs)} stable Ga-N semiconductors:")
        for doc in docs:
            gap_type = "direct" if doc.is_gap_direct else "indirect"
            print(f"  {doc.material_id}: {doc.formula_pretty}")
            print(f"    Band gap: {doc.band_gap:.2f} eV ({gap_type})")
            print(f"    E_hull: {doc.energy_above_hull:.3f} eV/atom")
            print(f"    E_form: {doc.formation_energy_per_atom:.3f} eV/atom")
        print()


def example_3_structure_retrieval():
    """Get and save crystal structures."""
    print("Example 3: Structure Retrieval and Export")
    print("=" * 50)

    with MPRester() as mpr:
        # Get structure by material ID
        material_id = "mp-149"  # Silicon
        structure = mpr.get_structure_by_material_id(material_id)

        print(f"Retrieved structure: {structure.composition.reduced_formula}")
        print(f"Space group: {structure.get_space_group_info()}")
        print(f"Lattice:\n{structure.lattice}")

        # Save in different formats
        structure.to(filename="Si_mp149.cif")
        structure.to(filename="Si_mp149_POSCAR", fmt="poscar")
        structure.to(filename="Si_mp149.json")

        print("\nSaved structure as:")
        print("  - Si_mp149.cif")
        print("  - Si_mp149_POSCAR")
        print("  - Si_mp149.json")
        print()


def example_4_battery_materials():
    """Search for battery cathode materials."""
    print("Example 4: Battery Cathode Materials")
    print("=" * 50)

    with MPRester() as mpr:
        # Search for lithium-containing stable oxides
        docs = mpr.materials.summary.search(
            elements=["Li", "O"],
            energy_above_hull=(0, 0.02),  # Very stable
            num_elements=(2, 4),  # Binary to quaternary
            fields=[
                "material_id",
                "formula_pretty",
                "energy_above_hull",
                "formation_energy_per_atom",
                "density",
            ],
        )

        print(f"Found {len(docs)} stable Li-containing oxides:")
        for doc in docs[:10]:  # Show first 10
            print(f"  {doc.material_id}: {doc.formula_pretty}")
            print(f"    E_hull: {doc.energy_above_hull:.4f} eV/atom")
            print(f"    Density: {doc.density:.2f} g/cm³")
        print()


def example_5_comparison_analysis():
    """Compare different polymorphs of same composition."""
    print("Example 5: Polymorph Comparison - TiO2")
    print("=" * 50)

    with MPRester() as mpr:
        # Get all TiO2 structures
        docs = mpr.materials.summary.search(
            formula="TiO2",
            fields=[
                "material_id",
                "formula_pretty",
                "energy_per_atom",
                "energy_above_hull",
                "band_gap",
                "density",
                "symmetry",
            ],
        )

        print(f"Found {len(docs)} TiO2 polymorphs:")
        # Sort by stability
        docs_sorted = sorted(docs, key=lambda x: x.energy_above_hull)

        for doc in docs_sorted:
            spg = doc.symmetry.symbol if doc.symmetry else "Unknown"
            print(f"  {doc.material_id} ({spg}):")
            print(f"    E_hull: {doc.energy_above_hull:.3f} eV/atom")
            print(f"    Band gap: {doc.band_gap:.2f} eV" if doc.band_gap else "    Band gap: N/A")
            print(f"    Density: {doc.density:.2f} g/cm³")
        print()


def example_6_metal_search():
    """Find metallic elements and compounds."""
    print("Example 6: Metallic Materials Search")
    print("=" * 50)

    with MPRester() as mpr:
        # Search for stable metallic materials
        docs = mpr.materials.summary.search(
            elements=["Cu", "Ag", "Au"],
            is_metal=True,
            energy_above_hull=(0, 0.01),
            fields=["material_id", "formula_pretty", "energy_above_hull", "density", "nelements"],
        )

        print(f"Found {len(docs)} metallic Cu/Ag/Au materials:")
        for doc in docs:
            print(f"  {doc.material_id}: {doc.formula_pretty}")
            print(f"    Elements: {doc.nelements}")
            print(f"    Density: {doc.density:.2f} g/cm³")
        print()


def example_7_wide_bandgap_oxides():
    """Search for wide band gap oxide semiconductors."""
    print("Example 7: Wide Band Gap Oxide Semiconductors")
    print("=" * 50)

    with MPRester() as mpr:
        # Wide band gap oxides for UV applications
        docs = mpr.materials.summary.search(
            elements=["O"],
            band_gap=(3.0, 6.0),  # Wide band gap
            energy_above_hull=(0, 0.05),
            fields=[
                "material_id",
                "formula_pretty",
                "band_gap",
                "is_gap_direct",
                "energy_above_hull",
                "elements",
            ],
        )

        print(f"Found {len(docs)} wide bandgap oxides:")
        for doc in docs[:15]:  # Show first 15
            gap_type = "direct" if doc.is_gap_direct else "indirect"
            elements_str = "-".join(doc.elements)
            print(f"  {doc.material_id}: {doc.formula_pretty} ({elements_str})")
            print(f"    Band gap: {doc.band_gap:.2f} eV ({gap_type})")
            print(f"    E_hull: {doc.energy_above_hull:.3f} eV/atom")
        print()


def example_8_phase_stability():
    """Analyze phase stability in a chemical system."""
    print("Example 8: Phase Stability - Li-Fe-O System")
    print("=" * 50)

    with MPRester() as mpr:
        # Get all phases in Li-Fe-O system
        docs = mpr.materials.summary.search(
            chemsys="Li-Fe-O",
            fields=[
                "material_id",
                "formula_pretty",
                "energy_above_hull",
                "is_stable",
                "formation_energy_per_atom",
            ],
        )

        stable = [d for d in docs if d.is_stable]
        metastable = [d for d in docs if not d.is_stable and d.energy_above_hull < 0.1]

        print(f"Stable phases ({len(stable)}):")
        for doc in stable:
            print(f"  {doc.material_id}: {doc.formula_pretty}")
            print(f"    E_form: {doc.formation_energy_per_atom:.3f} eV/atom")

        print(f"\nMetastable phases (E_hull < 0.1 eV/atom, {len(metastable)}):")
        for doc in metastable[:5]:
            print(f"  {doc.material_id}: {doc.formula_pretty}")
            print(f"    E_hull: {doc.energy_above_hull:.3f} eV/atom")
        print()


def example_9_small_unit_cells():
    """Find materials with small unit cells for fast calculations."""
    print("Example 9: Small Unit Cell Materials")
    print("=" * 50)

    with MPRester() as mpr:
        # Find binary compounds with small unit cells
        docs = mpr.materials.summary.search(
            num_elements=2,
            nsites=(1, 10),  # 10 or fewer atoms
            energy_above_hull=(0, 0.02),
            band_gap=(0.5, None),  # Semiconducting or insulating
            fields=["material_id", "formula_pretty", "nsites", "band_gap", "symmetry"],
        )

        print(f"Found {len(docs)} binary compounds with ≤10 atoms:")
        for doc in sorted(docs, key=lambda x: x.nsites)[:15]:
            spg = doc.symmetry.symbol if doc.symmetry else "Unknown"
            print(f"  {doc.material_id}: {doc.formula_pretty}")
            print(f"    Atoms: {doc.nsites}, Space group: {spg}")
            print(f"    Band gap: {doc.band_gap:.2f} eV")
        print()


def example_10_export_dataset():
    """Export dataset for machine learning."""
    print("Example 10: Export Dataset for ML")
    print("=" * 50)

    import pandas as pd

    with MPRester() as mpr:
        # Get binary oxides with various properties
        docs = mpr.materials.summary.search(
            elements=["O"],
            num_elements=2,
            energy_above_hull=(0, 0.05),
            fields=[
                "material_id",
                "formula_pretty",
                "band_gap",
                "formation_energy_per_atom",
                "energy_above_hull",
                "density",
                "volume",
                "nsites",
            ],
        )

        # Convert to DataFrame
        data = []
        for doc in docs:
            data.append(
                {
                    "material_id": doc.material_id,
                    "formula": doc.formula_pretty,
                    "band_gap_eV": doc.band_gap,
                    "formation_energy_eV_atom": doc.formation_energy_per_atom,
                    "energy_above_hull_eV_atom": doc.energy_above_hull,
                    "density_g_cm3": doc.density,
                    "volume_A3": doc.volume,
                    "n_atoms": doc.nsites,
                }
            )

        df = pd.DataFrame(data)

        # Save to CSV
        df.to_csv("binary_oxides_dataset.csv", index=False)

        print(f"Exported {len(df)} binary oxides to 'binary_oxides_dataset.csv'")
        print("\nDataset statistics:")
        print(df.describe())
        print()


if __name__ == "__main__":
    # Run all examples
    examples = [
        example_1_basic_formula_search,
        example_2_property_screening,
        example_3_structure_retrieval,
        example_4_battery_materials,
        example_5_comparison_analysis,
        example_6_metal_search,
        example_7_wide_bandgap_oxides,
        example_8_phase_stability,
        example_9_small_unit_cells,
        example_10_export_dataset,
    ]

    for example in examples:
        try:
            example()
        except Exception as e:
            print(f"Error in {example.__name__}: {e}")
            print()

    print("All examples completed!")
