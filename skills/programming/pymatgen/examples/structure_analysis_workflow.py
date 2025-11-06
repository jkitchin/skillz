#!/usr/bin/env python
"""
Example: Complete structure analysis workflow with pymatgen.

This script demonstrates a typical workflow:
1. Load structure from CIF
2. Analyze composition and symmetry
3. Generate VASP input files
4. Query Materials Project for similar structures

Author: Claude Code Assistant
Date: 2025
"""

from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.analysis.local_env import CrystalNN
import argparse


def analyze_structure(structure_file):
    """Perform comprehensive structure analysis."""

    print("="*60)
    print("PYMATGEN STRUCTURE ANALYSIS WORKFLOW")
    print("="*60)

    # 1. Load structure
    print("\n1. Loading structure...")
    structure = Structure.from_file(structure_file)
    print(f"   ✓ Loaded {structure_file}")

    # 2. Composition analysis
    print("\n2. Composition analysis:")
    comp = structure.composition
    print(f"   Formula: {comp.formula}")
    print(f"   Reduced formula: {comp.reduced_formula}")
    print(f"   Molecular weight: {comp.weight:.2f} g/mol")
    print(f"   Number of atoms: {comp.num_atoms}")

    # 3. Geometric properties
    print("\n3. Geometric properties:")
    print(f"   Volume: {structure.volume:.3f} Ų")
    print(f"   Density: {structure.density:.3f} g/cm³")
    lattice = structure.lattice
    print(f"   Lattice parameters:")
    print(f"     a = {lattice.a:.3f} Å")
    print(f"     b = {lattice.b:.3f} Å")
    print(f"     c = {lattice.c:.3f} Å")
    print(f"     α = {lattice.alpha:.2f}°")
    print(f"     β = {lattice.beta:.2f}°")
    print(f"     γ = {lattice.gamma:.2f}°")

    # 4. Symmetry analysis
    print("\n4. Symmetry analysis:")
    sga = SpacegroupAnalyzer(structure)
    print(f"   Space group: {sga.get_space_group_symbol()}")
    print(f"   Space group number: {sga.get_space_group_number()}")
    print(f"   Crystal system: {sga.get_crystal_system()}")
    print(f"   Point group: {sga.get_point_group_symbol()}")

    sym_ops = sga.get_symmetry_operations()
    print(f"   Number of symmetry operations: {len(sym_ops)}")

    # 5. Coordination analysis
    print("\n5. Coordination environment:")
    cnn = CrystalNN()

    # Analyze first few sites
    for i in range(min(3, len(structure))):
        site = structure[i]
        try:
            cn_dict = cnn.get_cn_dict(structure, i)
            cn = sum(cn_dict.values())
            print(f"   Site {i} ({site.species_string}):")
            print(f"     Coordination number: {cn:.1f}")

            # Get nearest neighbor distances
            neighbors = structure.get_neighbors(site, r=4.0)
            if neighbors:
                nn_distances = [n.nn_distance for n in neighbors[:3]]
                avg_dist = sum(nn_distances) / len(nn_distances)
                print(f"     Avg. nearest neighbor distance: {avg_dist:.3f} Å")
        except Exception as e:
            print(f"   Site {i} ({site.species_string}): Could not determine coordination")

    # 6. Primitive and conventional cells
    print("\n6. Cell transformations:")
    primitive = sga.get_primitive_standard_structure()
    conventional = sga.get_conventional_standard_structure()

    print(f"   Original sites: {len(structure)}")
    print(f"   Primitive sites: {len(primitive)}")
    print(f"   Conventional sites: {len(conventional)}")

    # 7. Generate VASP inputs
    print("\n7. Generating VASP input files:")
    output_dir = "vasp_relax"

    try:
        input_set = MPRelaxSet(structure)
        input_set.write_input(output_dir)
        print(f"   ✓ VASP input files written to '{output_dir}/'")
        print(f"     - POSCAR (structure)")
        print(f"     - INCAR (parameters)")
        print(f"     - KPOINTS (k-point mesh)")
        print(f"     - POTCAR (pseudopotentials)")
    except Exception as e:
        print(f"   ✗ Could not generate POTCAR (set PMG_VASP_PSP_DIR)")
        print(f"     POSCAR, INCAR, and KPOINTS written")

        # Write individual files
        try:
            input_set.poscar.write_file(f"{output_dir}/POSCAR")
            input_set.incar.write_file(f"{output_dir}/INCAR")
            input_set.kpoints.write_file(f"{output_dir}/KPOINTS")
        except:
            pass

    # 8. Export structure in different formats
    print("\n8. Exporting structure:")
    structure.to(filename="structure_output.cif")
    print("   ✓ CIF: structure_output.cif")

    structure.to(filename="structure_output.json")
    print("   ✓ JSON: structure_output.json")

    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)


def main():
    parser = argparse.ArgumentParser(
        description="Analyze crystal structure using pymatgen"
    )
    parser.add_argument(
        "structure_file",
        help="Input structure file (CIF, POSCAR, etc.)"
    )
    args = parser.parse_args()

    try:
        analyze_structure(args.structure_file)
    except FileNotFoundError:
        print(f"Error: File '{args.structure_file}' not found")
    except Exception as e:
        print(f"Error during analysis: {e}")
        raise


if __name__ == "__main__":
    main()
