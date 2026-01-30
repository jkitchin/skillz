"""
AFLOW Database Search Examples

Demonstrates various search patterns using the AFLOW REST API.
No API key required!
"""

import requests
import json


def example_1_basic_property_query():
    """Get a specific property from a known material."""
    print("Example 1: Basic Property Query")
    print("=" * 50)

    # Get formation enthalpy of FCC Silver
    url = "http://aflowlib.duke.edu/AFLOWDATA/ICSD_WEB/FCC/Ag1/?enthalpy_formation_atom"
    response = requests.get(url)
    enthalpy = float(response.text.strip())

    print(f"FCC Silver formation enthalpy: {enthalpy:.4f} eV/atom")

    # Get band gap
    url = "http://aflowlib.duke.edu/AFLOWDATA/ICSD_WEB/FCC/Ag1/?Egap"
    response = requests.get(url)
    egap = float(response.text.strip())

    print(f"FCC Silver band gap: {egap:.4f} eV")
    print()


def example_2_list_properties():
    """List all available properties for an entry."""
    print("Example 2: List Available Properties")
    print("=" * 50)

    url = "http://aflowlib.duke.edu/AFLOWDATA/ICSD_WEB/FCC/Ag1/?keywords"
    response = requests.get(url)
    keywords = response.text.strip().split(",")

    print("Available properties for FCC Ag:")
    for i, kw in enumerate(keywords[:20], 1):  # Show first 20
        print(f"  {i}. {kw}")
    print(f"  ... and {len(keywords) - 20} more")
    print()


def example_3_get_structure():
    """Download crystal structure."""
    print("Example 3: Get Crystal Structure (POSCAR)")
    print("=" * 50)

    # Get POSCAR for FCC Silver
    url = "http://aflowlib.duke.edu/AFLOWDATA/ICSD_WEB/FCC/Ag1/?geometry"
    response = requests.get(url)
    poscar = response.text

    # Save to file
    with open("Ag_FCC.vasp", "w") as f:
        f.write(poscar)

    print("POSCAR content (first 10 lines):")
    print("\n".join(poscar.split("\n")[:10]))
    print("\nSaved to: Ag_FCC.vasp")
    print()


def example_4_aflux_search_by_composition():
    """Search using AFLUX by chemical composition."""
    print("Example 4: AFLUX Search by Composition")
    print("=" * 50)

    # Find all Gold-containing materials
    url = "http://aflowlib.duke.edu/search/API/?species(Au),nspecies(1)"
    response = requests.get(url, timeout=30)
    results = response.json()

    print(f"Found {len(results)} pure Gold structures:")
    for result in results[:5]:  # Show first 5
        print(f"  Compound: {result.get('compound', 'N/A')}")
        print(f"  Formation enthalpy: {result.get('enthalpy_formation_atom', 'N/A')} eV/atom")
        print(f"  Band gap: {result.get('Egap', 'N/A')} eV")
        print(f"  URL: {result.get('aurl', 'N/A')}")
        print()


def example_5_aflux_property_range():
    """Search by property ranges using AFLUX."""
    print("Example 5: AFLUX Property Range Search")
    print("=" * 50)

    # Find materials with band gap between 2-4 eV
    url = "http://aflowlib.duke.edu/search/API/?Egap(2*,4*),nspecies(2)"
    response = requests.get(url, timeout=30)
    results = response.json()

    print(f"Found {len(results)} binary materials with Egap = 2-4 eV:")
    for result in results[:10]:  # Show first 10
        compound = result.get("compound", "N/A")
        egap = result.get("Egap", "N/A")
        species = result.get("species", [])

        print(f"  {compound}: Egap = {egap} eV ({', '.join(species)})")
    print()


def example_6_stable_compounds():
    """Find thermodynamically stable compounds."""
    print("Example 6: Stable Compounds (negative formation enthalpy)")
    print("=" * 50)

    # Search for stable binary Ti-O compounds
    url = (
        "http://aflowlib.duke.edu/search/API/?"
        "species(Ti,O),nspecies(2),enthalpy_formation_atom(*,0*)"
    )
    response = requests.get(url, timeout=30)
    results = response.json()

    print(f"Found {len(results)} stable Ti-O compounds:")
    for result in results:
        compound = result.get("compound", "N/A")
        enthalpy = result.get("enthalpy_formation_atom", "N/A")
        egap = result.get("Egap", "N/A")

        print(f"  {compound}:")
        print(f"    Formation enthalpy: {enthalpy} eV/atom")
        print(f"    Band gap: {egap} eV")
    print()


def example_7_elastic_properties():
    """Search for materials with specific elastic properties."""
    print("Example 7: Materials with High Bulk Modulus")
    print("=" * 50)

    # Find materials with bulk modulus > 200 GPa
    url = "http://aflowlib.duke.edu/search/API/?bulk_modulus_voigt(200*,*),nspecies(1)"
    response = requests.get(url, timeout=30)
    results = response.json()

    print(f"Found {len(results)} elements with bulk modulus > 200 GPa:")
    for result in results[:10]:
        compound = result.get("compound", "N/A")
        bulk_mod = result.get("bulk_modulus_voigt", "N/A")
        shear_mod = result.get("shear_modulus_voigt", "N/A")

        print(f"  {compound}:")
        print(f"    Bulk modulus: {bulk_mod} GPa")
        print(f"    Shear modulus: {shear_mod} GPa")
    print()


def example_8_download_multiple_structures():
    """Download structures for an entire chemical system."""
    print("Example 8: Download Ga-N System Structures")
    print("=" * 50)

    import os

    # Search for Ga-N compounds
    search_url = "http://aflowlib.duke.edu/search/API/?species(Ga,N),nspecies(2)"
    response = requests.get(search_url, timeout=30)
    results = response.json()

    print(f"Found {len(results)} Ga-N structures")

    # Create directory
    os.makedirs("GaN_structures", exist_ok=True)

    # Download first 5 structures
    for i, result in enumerate(results[:5]):
        aurl = result.get("aurl", "")
        compound = result.get("compound", f"structure_{i}").replace("/", "_")

        # Get POSCAR
        structure_url = f"{aurl}?geometry"
        try:
            poscar_response = requests.get(structure_url, timeout=10)
            poscar = poscar_response.text

            # Save
            filename = f"GaN_structures/{compound}.vasp"
            with open(filename, "w") as f:
                f.write(poscar)

            print(f"  Downloaded: {compound}")
        except Exception as e:
            print(f"  Error downloading {compound}: {e}")

    print()


def example_9_complete_entry_data():
    """Get complete entry data as JSON."""
    print("Example 9: Complete Entry Data")
    print("=" * 50)

    url = "http://aflowlib.duke.edu/AFLOWDATA/ICSD_WEB/FCC/Ag1/?aflowlib_entries"
    response = requests.get(url)
    data = response.json()

    print("Complete entry data for FCC Ag:")
    print(f"  Compound: {data.get('compound', 'N/A')}")
    print(f"  Species: {data.get('species', 'N/A')}")
    print(f"  Number of atoms: {data.get('natoms', 'N/A')}")
    print(f"  Volume: {data.get('volume_cell', 'N/A')} Ų")
    print(f"  Density: {data.get('density', 'N/A')} g/cm³")
    print(f"  Space group: {data.get('spacegroup_relax', 'N/A')}")
    print(f"  Bravais lattice: {data.get('Bravais_lattice_relax', 'N/A')}")
    print(f"  Formation enthalpy: {data.get('enthalpy_formation_atom', 'N/A')} eV/atom")
    print(f"  Band gap: {data.get('Egap', 'N/A')} eV")
    print()


def example_10_catalog_specific_search():
    """Search within specific AFLOW catalog."""
    print("Example 10: Catalog-Specific Search (ICSD)")
    print("=" * 50)

    # Search only ICSD catalog
    url = "http://aflowlib.duke.edu/search/API/?species(Fe,O),nspecies(2),catalog(ICSD)"
    response = requests.get(url, timeout=30)
    results = response.json()

    print(f"Found {len(results)} Fe-O compounds in ICSD catalog:")
    for result in results[:10]:
        compound = result.get("compound", "N/A")
        spacegroup = result.get("spacegroup_relax", "N/A")
        egap = result.get("Egap", "N/A")

        print(f"  {compound} (SG: {spacegroup})")
        print(f"    Band gap: {egap} eV")
    print()


def example_11_small_unit_cells():
    """Find materials with small unit cells."""
    print("Example 11: Small Unit Cell Materials")
    print("=" * 50)

    # Binary compounds with < 10 atoms
    url = "http://aflowlib.duke.edu/search/API/?nspecies(2),natoms(*,10*)"
    response = requests.get(url, timeout=30)
    results = response.json()

    print(f"Found {len(results)} binary compounds with ≤10 atoms:")
    for result in results[:15]:
        compound = result.get("compound", "N/A")
        natoms = result.get("natoms", "N/A")
        volume = result.get("volume_cell", "N/A")

        print(f"  {compound}: {natoms} atoms, V = {volume} ų")
    print()


def example_12_using_aflow_package():
    """Use the aflow Python package for easier access."""
    print("Example 12: Using aflow Python Package")
    print("=" * 50)

    try:
        import aflow

        # Search for materials
        results = aflow.search(filter="species(Si),nspecies(1)")

        print("Silicon structures from aflow package:")
        for i, entry in enumerate(results):
            if i >= 5:  # Show only first 5
                break

            print(f"  Compound: {entry.compound}")
            print(f"  Formation enthalpy: {entry.enthalpy_formation_atom} eV/atom")
            print(f"  Band gap: {entry.Egap} eV")
            print(f"  Space group: {entry.spacegroup_relax}")

            # Get ASE Atoms object
            atoms = entry.atoms
            print(f"  Formula from ASE: {atoms.get_chemical_formula()}")
            print()

    except ImportError:
        print("aflow package not installed.")
        print("Install with: pip install aflow")
        print()


def example_13_compare_polymorphs():
    """Compare different polymorphs of the same composition."""
    print("Example 13: Compare Carbon Polymorphs")
    print("=" * 50)

    # Search for pure carbon structures
    url = "http://aflowlib.duke.edu/search/API/?species(C),nspecies(1)"
    response = requests.get(url, timeout=30)
    results = response.json()

    print(f"Found {len(results)} carbon polymorphs:")
    for result in results:
        compound = result.get("compound", "N/A")
        enthalpy = result.get("enthalpy_formation_atom", 0)
        spacegroup = result.get("spacegroup_relax", "N/A")
        density = result.get("density", "N/A")

        print(f"  {compound} (SG: {spacegroup})")
        print(f"    Formation enthalpy: {enthalpy:.4f} eV/atom")
        print(f"    Density: {density} g/cm³")
    print()


def example_14_export_to_dataframe():
    """Export search results to pandas DataFrame."""
    print("Example 14: Export to Pandas DataFrame")
    print("=" * 50)

    try:
        import pandas as pd

        # Search for ternary oxides
        url = "http://aflowlib.duke.edu/search/API/?species(O),nspecies(3)"
        response = requests.get(url, timeout=30)
        results = response.json()

        # Convert to DataFrame
        data = []
        for result in results[:50]:  # First 50
            data.append(
                {
                    "compound": result.get("compound", "N/A"),
                    "enthalpy_formation_atom": result.get("enthalpy_formation_atom"),
                    "Egap_eV": result.get("Egap"),
                    "density_g_cm3": result.get("density"),
                    "volume_A3": result.get("volume_cell"),
                    "natoms": result.get("natoms"),
                    "spacegroup": result.get("spacegroup_relax"),
                }
            )

        df = pd.DataFrame(data)

        # Save to CSV
        df.to_csv("aflow_ternary_oxides.csv", index=False)

        print(f"Exported {len(df)} ternary oxides to CSV")
        print("\nDataset preview:")
        print(df.head())
        print("\nStatistics:")
        print(df.describe())
        print()

    except ImportError:
        print("pandas not installed.")
        print("Install with: pip install pandas")
        print()


def example_15_advanced_filtering():
    """Advanced multi-property filtering."""
    print("Example 15: Advanced Multi-Property Filtering")
    print("=" * 50)

    # Stable semiconducting binary compounds with moderate density
    url = (
        "http://aflowlib.duke.edu/search/API/?"
        "nspecies(2),"
        "enthalpy_formation_atom(*,0*),"  # Stable
        "Egap(1*,3*),"  # Semiconducting
        "density(3*,6*)"
    )  # Moderate density
    response = requests.get(url, timeout=30)
    results = response.json()

    print(f"Found {len(results)} materials matching all criteria:")
    print("(Binary, stable, semiconducting, density 3-6 g/cm³)")
    print()

    for result in results[:10]:
        compound = result.get("compound", "N/A")
        enthalpy = result.get("enthalpy_formation_atom", "N/A")
        egap = result.get("Egap", "N/A")
        density = result.get("density", "N/A")

        print(f"  {compound}:")
        print(f"    E_form: {enthalpy} eV/atom")
        print(f"    E_gap: {egap} eV")
        print(f"    Density: {density} g/cm³")
    print()


if __name__ == "__main__":
    # Run all examples
    examples = [
        example_1_basic_property_query,
        example_2_list_properties,
        example_3_get_structure,
        example_4_aflux_search_by_composition,
        example_5_aflux_property_range,
        example_6_stable_compounds,
        example_7_elastic_properties,
        example_8_download_multiple_structures,
        example_9_complete_entry_data,
        example_10_catalog_specific_search,
        example_11_small_unit_cells,
        example_12_using_aflow_package,
        example_13_compare_polymorphs,
        example_14_export_to_dataframe,
        example_15_advanced_filtering,
    ]

    for example in examples:
        try:
            example()
        except Exception as e:
            print(f"Error in {example.__name__}: {e}")
            print()

    print("All examples completed!")
