"""
Combined Database Workflow Examples

Demonstrates how to use both Materials Project and AFLOW together
for comprehensive materials research.
"""

import requests
from mp_api.client import MPRester


def example_1_cross_database_comparison():
    """Compare same material from both databases."""
    print("Example 1: Cross-Database Comparison - GaN")
    print("=" * 50)

    # Materials Project
    print("Materials Project Data:")
    with MPRester() as mpr:
        docs = mpr.materials.summary.search(
            formula="GaN",
            fields=["material_id", "band_gap", "formation_energy_per_atom", "density"],
        )

        for doc in docs[:3]:
            print(f"  {doc.material_id}:")
            print(f"    Band gap: {doc.band_gap} eV")
            print(f"    E_form: {doc.formation_energy_per_atom:.3f} eV/atom")
            print(f"    Density: {doc.density:.2f} g/cm³")

    # AFLOW
    print("\nAFLOW Data:")
    url = "http://aflowlib.duke.edu/search/API/?species(Ga,N),nspecies(2)"
    response = requests.get(url, timeout=30)
    results = response.json()

    for result in results[:3]:
        compound = result.get("compound", "N/A")
        egap = result.get("Egap", "N/A")
        enthalpy = result.get("enthalpy_formation_atom", "N/A")
        density = result.get("density", "N/A")

        print(f"  {compound}:")
        print(f"    Band gap: {egap} eV")
        print(f"    E_form: {enthalpy} eV/atom")
        print(f"    Density: {density} g/cm³")
    print()


def example_2_complementary_properties():
    """Get different properties from each database."""
    print("Example 2: Complementary Properties - TiO2")
    print("=" * 50)

    formula = "TiO2"

    # Materials Project: Electronic structure
    print(f"Materials Project - Electronic properties of {formula}:")
    with MPRester() as mpr:
        docs = mpr.materials.summary.search(
            formula=formula,
            fields=["material_id", "band_gap", "is_gap_direct", "energy_above_hull"],
        )

        for doc in docs[:2]:
            gap_type = "direct" if doc.is_gap_direct else "indirect"
            print(f"  {doc.material_id}:")
            print(f"    Band gap: {doc.band_gap} eV ({gap_type})")
            print(f"    Stability: {doc.energy_above_hull:.3f} eV/atom")

    # AFLOW: Elastic properties
    print(f"\nAFLOW - Elastic properties of {formula}:")
    url = "http://aflowlib.duke.edu/search/API/?species(Ti,O),nspecies(2)"
    response = requests.get(url, timeout=30)
    results = response.json()

    for result in results[:2]:
        compound = result.get("compound", "N/A")
        bulk_mod = result.get("bulk_modulus_voigt", "N/A")
        shear_mod = result.get("shear_modulus_voigt", "N/A")

        print(f"  {compound}:")
        print(f"    Bulk modulus: {bulk_mod} GPa")
        print(f"    Shear modulus: {shear_mod} GPa")
    print()


def example_3_validate_predictions():
    """Use AFLOW's ICSD data to validate MP predictions."""
    print("Example 3: Validate Against Experimental (ICSD)")
    print("=" * 50)

    # Get MP prediction
    print("Materials Project prediction:")
    with MPRester() as mpr:
        docs = mpr.materials.summary.search(
            formula="Si",
            theoretical=False,  # Get experimental
            fields=["material_id", "density", "volume", "symmetry"],
        )

        mp_doc = docs[0]
        print(f"  {mp_doc.material_id}:")
        print(f"    Density: {mp_doc.density:.3f} g/cm³")
        print(f"    Volume: {mp_doc.volume:.3f} ų")
        print(f"    Space group: {mp_doc.symmetry.symbol if mp_doc.symmetry else 'N/A'}")

    # Compare with AFLOW ICSD
    print("\nAFLOW (ICSD experimental):")
    url = "http://aflowlib.duke.edu/search/API/?species(Si),nspecies(1),catalog(ICSD)"
    response = requests.get(url, timeout=30)
    results = response.json()

    if results:
        result = results[0]
        print(f"  {result.get('compound', 'N/A')}:")
        print(f"    Density: {result.get('density', 'N/A')} g/cm³")
        print(f"    Volume: {result.get('volume_cell', 'N/A')} ų")
        print(f"    Space group: {result.get('spacegroup_relax', 'N/A')}")
    print()


def example_4_comprehensive_screening():
    """Screen both databases for candidate materials."""
    print("Example 4: Comprehensive Materials Screening")
    print("=" * 50)
    print("Finding stable semiconducting binary oxides\n")

    # Materials Project search
    print("Materials Project results:")
    with MPRester() as mpr:
        mp_docs = mpr.materials.summary.search(
            elements=["O"],
            num_elements=2,
            energy_above_hull=(0, 0.05),
            band_gap=(1.0, 3.0),
            fields=["material_id", "formula_pretty", "band_gap", "energy_above_hull"],
        )

        mp_formulas = set()
        for doc in mp_docs[:10]:
            print(
                f"  {doc.formula_pretty}: Egap = {doc.band_gap:.2f} eV, "
                f"E_hull = {doc.energy_above_hull:.3f}"
            )
            mp_formulas.add(doc.formula_pretty)

    # AFLOW search
    print(f"\nAFLOW results:")
    url = (
        "http://aflowlib.duke.edu/search/API/?"
        "species(O),nspecies(2),Egap(1*,3*),enthalpy_formation_atom(*,0*)"
    )
    response = requests.get(url, timeout=30)
    aflow_results = response.json()

    aflow_formulas = set()
    for result in aflow_results[:10]:
        compound = result.get("compound", "N/A")
        egap = result.get("Egap", "N/A")
        print(f"  {compound}: Egap = {egap} eV")
        # Extract simple formula (simplified)
        if compound != "N/A":
            aflow_formulas.add(compound.split("/")[-1] if "/" in compound else compound)

    # Find overlap
    print(f"\nTotal unique formulas:")
    print(f"  Materials Project: {len(mp_formulas)}")
    print(f"  AFLOW: {len(aflow_formulas)}")
    print()


def example_5_structure_comparison():
    """Download and compare structures from both databases."""
    print("Example 5: Structure Comparison")
    print("=" * 50)

    import os

    os.makedirs("structures_comparison", exist_ok=True)

    material = "GaN"

    # Materials Project
    print(f"Downloading {material} from Materials Project...")
    with MPRester() as mpr:
        docs = mpr.materials.summary.search(formula=material, fields=["material_id"])

        if docs:
            mp_id = docs[0].material_id
            structure = mpr.get_structure_by_material_id(mp_id)
            structure.to(filename=f"structures_comparison/{material}_MP.cif")
            print(f"  Saved: {material}_MP.cif ({mp_id})")

    # AFLOW
    print(f"Downloading {material} from AFLOW...")
    search_url = f"http://aflowlib.duke.edu/search/API/?species(Ga,N),nspecies(2)"
    response = requests.get(search_url, timeout=30)
    results = response.json()

    if results:
        aurl = results[0].get("aurl", "")
        structure_url = f"{aurl}?geometry"
        poscar = requests.get(structure_url, timeout=10).text

        with open(f"structures_comparison/{material}_AFLOW.vasp", "w") as f:
            f.write(poscar)
        print(f"  Saved: {material}_AFLOW.vasp")

    print("\nStructures saved for comparison!")
    print()


def example_6_build_combined_dataset():
    """Build a dataset combining data from both databases."""
    print("Example 6: Build Combined Dataset")
    print("=" * 50)

    try:
        import pandas as pd

        # Target system: binary oxides
        data = []

        # Get MP data
        print("Fetching Materials Project data...")
        with MPRester() as mpr:
            mp_docs = mpr.materials.summary.search(
                elements=["O"],
                num_elements=2,
                energy_above_hull=(0, 0.05),
                fields=[
                    "formula_pretty",
                    "band_gap",
                    "formation_energy_per_atom",
                    "energy_above_hull",
                ],
            )

            for doc in mp_docs[:20]:
                data.append(
                    {
                        "formula": doc.formula_pretty,
                        "source": "MP",
                        "band_gap_eV": doc.band_gap,
                        "formation_energy_eV_atom": doc.formation_energy_per_atom,
                        "energy_above_hull_eV_atom": doc.energy_above_hull,
                    }
                )

        # Get AFLOW data
        print("Fetching AFLOW data...")
        url = (
            "http://aflowlib.duke.edu/search/API/?"
            "species(O),nspecies(2),enthalpy_formation_atom(*,0*)"
        )
        response = requests.get(url, timeout=30)
        aflow_results = response.json()

        for result in aflow_results[:20]:
            compound = result.get("compound", "N/A")
            # Simplify compound name
            formula = compound.split("/")[-1] if "/" in compound else compound

            data.append(
                {
                    "formula": formula,
                    "source": "AFLOW",
                    "band_gap_eV": result.get("Egap"),
                    "formation_energy_eV_atom": result.get("enthalpy_formation_atom"),
                    "energy_above_hull_eV_atom": None,  # Not in AFLOW
                }
            )

        # Create DataFrame
        df = pd.DataFrame(data)
        df.to_csv("combined_materials_data.csv", index=False)

        print(f"\nCombined dataset created: {len(df)} entries")
        print(f"  Materials Project: {len(df[df['source'] == 'MP'])}")
        print(f"  AFLOW: {len(df[df['source'] == 'AFLOW'])}")
        print("\nPreview:")
        print(df.head(10))
        print()

    except ImportError:
        print("pandas required for this example")
        print()


def example_7_complementary_search_strategy():
    """Use both databases with different search strategies."""
    print("Example 7: Complementary Search Strategy")
    print("=" * 50)
    print("Strategy: Use MP for targeted search, AFLOW for broad screening\n")

    target_elements = ["Li", "Co", "O"]

    # Materials Project: Detailed targeted search
    print("Materials Project - Detailed search with stability criteria:")
    with MPRester() as mpr:
        mp_docs = mpr.materials.summary.search(
            elements=target_elements,
            energy_above_hull=(0, 0.02),  # Very stable
            fields=["material_id", "formula_pretty", "energy_above_hull", "band_gap"],
        )

        print(f"Found {len(mp_docs)} highly stable {'-'.join(target_elements)} compounds")
        for doc in mp_docs[:5]:
            print(f"  {doc.material_id}: {doc.formula_pretty}")
            print(f"    Stability: {doc.energy_above_hull:.4f} eV/atom")
            print(f"    Band gap: {doc.band_gap} eV")

    # AFLOW: Broad screening
    print(f"\nAFLOW - Broad screening of {'-'.join(target_elements)} system:")
    url = f"http://aflowlib.duke.edu/search/API/?species({','.join(target_elements)})"
    response = requests.get(url, timeout=30)
    aflow_results = response.json()

    print(f"Found {len(aflow_results)} total {'-'.join(target_elements)} entries")
    print(f"(Use AFLOW to identify interesting compositions,")
    print(f" then get detailed properties from MP)")
    print()


def example_8_error_handling():
    """Demonstrate proper error handling for both APIs."""
    print("Example 8: Robust Error Handling")
    print("=" * 50)

    # Materials Project with error handling
    print("Materials Project query with error handling:")
    try:
        with MPRester() as mpr:
            docs = mpr.materials.summary.search(formula="InvalidFormula123", fields=["material_id"])
            print(f"  Found {len(docs)} results")
    except Exception as e:
        print(f"  Handled error: {type(e).__name__}")

    # AFLOW with error handling
    print("\nAFLOW query with error handling:")
    try:
        url = "http://aflowlib.duke.edu/search/API/?species(Xx,Yy)"  # Invalid
        response = requests.get(url, timeout=10)
        results = response.json() if response.ok else []
        print(f"  Found {len(results)} results")
    except requests.Timeout:
        print("  Handled error: Request timeout")
    except Exception as e:
        print(f"  Handled error: {type(e).__name__}")
    print()


if __name__ == "__main__":
    examples = [
        example_1_cross_database_comparison,
        example_2_complementary_properties,
        example_3_validate_predictions,
        example_4_comprehensive_screening,
        example_5_structure_comparison,
        example_6_build_combined_dataset,
        example_7_complementary_search_strategy,
        example_8_error_handling,
    ]

    for example in examples:
        try:
            example()
        except Exception as e:
            print(f"Error in {example.__name__}: {e}")
            print()

    print("All combined workflow examples completed!")
