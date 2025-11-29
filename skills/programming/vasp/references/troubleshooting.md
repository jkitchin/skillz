# VASP Troubleshooting Guide

Solutions to common VASP errors and problems.

## Table of Contents

1. [Convergence Issues](#convergence-issues)
2. [Geometry/Structure Problems](#geometrystructure-problems)
3. [Memory and Performance](#memory-and-performance)
4. [Input File Errors](#input-file-errors)
5. [Output Interpretation](#output-interpretation)
6. [Advanced Functionals](#advanced-functionals)

---

## Convergence Issues

### SCF Not Converging

**Error:** Electronic iterations reach NELM without converging

**OUTCAR shows:**
```
DAV:  60    -1.23456E+02    0.12345E-01    ...
         1 F= ...  E0= ...
```

**Solutions:**

1. **Try ALGO = All:**
   ```bash
   ALGO = All
   ```
   Tries multiple algorithms automatically.

2. **Increase NELM:**
   ```bash
   NELM = 200
   ```

3. **Adjust mixing:**
   ```bash
   AMIX = 0.2      # Default 0.4
   BMIX = 0.0001   # Default 1.0
   ```

4. **Check structure:**
   - Atoms too close? (< 1 Å separation)
   - Reasonable bonds?
   - Run with `ISYM = 0` to test if symmetry issue

5. **Tighten convergence temporarily:**
   ```bash
   EDIFF = 1E-4    # Looser for testing
   ```

6. **For metals, adjust smearing:**
   ```bash
   ISMEAR = 1
   SIGMA = 0.3     # Larger smearing
   ```

### Ionic Relaxation Not Converging

**Error:** NSW steps exhausted, forces still large

**Solutions:**

1. **Increase NSW:**
   ```bash
   NSW = 200
   ```

2. **Check forces in OUTCAR:**
   ```bash
   grep "TOTAL-FORCE" OUTCAR -A 20
   ```
   If only one atom has large force, may be specific issue.

3. **Try different optimizer:**
   ```bash
   IBRION = 1      # RMM-DIIS instead of CG
   # Or
   IBRION = 3      # Damped MD for rough surfaces
   ```

4. **Reduce step size:**
   ```bash
   POTIM = 0.2     # Smaller steps
   ```

5. **Tighten electronic convergence:**
   ```bash
   EDIFF = 1E-7
   ```

6. **Check for instability:**
   - May be saddle point, not minimum
   - Try starting from different geometry

### ZBRENT Error

**Error Message:**
```
ZBRENT: fatal error in bracketing
```

**Cause:** Line minimization failed (step too large)

**Solutions:**

1. **Reduce POTIM:**
   ```bash
   POTIM = 0.1     # Much smaller
   ```

2. **Use damped MD:**
   ```bash
   IBRION = 3
   POTIM = 0.1
   SMASS = 0.4
   ```

3. **Check initial structure:**
   - Bad geometry?
   - Pre-optimize with cheap method

### EDDDAV Error

**Error Message:**
```
EDDDAV: X eigenvalues not converged
```

**Cause:** Electronic structure not converging in Davidson algorithm

**Solutions:**

1. **Increase NELM:**
   ```bash
   NELM = 200
   ```

2. **Change algorithm:**
   ```bash
   ALGO = All
   ```

3. **Adjust mixing:**
   ```bash
   AMIX = 0.2
   BMIX = 0.0001
   ```

4. **Check system:**
   - Metallic? Use ISMEAR = 1, SIGMA = 0.2
   - Bad structure?

---

## Geometry/Structure Problems

### Atoms Exploding

**Symptom:** Forces become huge, atoms move far, calculation crashes

**Causes:**
- Time step too large (MD)
- Bad initial geometry
- Wrong POTCAR order

**Solutions:**

1. **Check POTCAR order matches POSCAR:**
   ```bash
   grep VRHFIN POTCAR
   # Should match order of elements in POSCAR line 6
   ```

2. **Reduce time step (MD):**
   ```bash
   POTIM = 0.5     # Or smaller
   ```

3. **Check initial structure:**
   ```bash
   # View in VESTA, Ovito, etc.
   # Look for overlapping atoms
   ```

4. **Start with static calculation:**
   ```bash
   IBRION = -1
   NSW = 0
   # Check if forces reasonable before MD/relaxation
   ```

### Cell Shape Becoming Unrealistic (ISIF=3)

**Symptom:** Relaxation with ISIF=3 gives distorted cell

**Solutions:**

1. **Two-step approach:**
   - First: ISIF=2 (relax ions, fixed cell)
   - Then: ISIF=3 (relax cell)

2. **Check symmetry:**
   ```bash
   ISYM = 2        # Use symmetry
   # Or verify symmetry is correct
   ```

3. **Apply pressure if needed:**
   ```bash
   PSTRESS = 0     # Ensure no external pressure
   ```

4. **Use better starting structure:**
   - Get experimental lattice parameters
   - Relax with tighter settings

### Sub-Space Matrix Not Hermitian

**Error Message:**
```
WARNING: Sub-Space-Matrix is not hermitian
```

**Cause:** Numerical instability, often from:
- Too large step size
- Bad structure
- Symmetry issues

**Solutions:**

1. **Reduce POTIM:**
   ```bash
   POTIM = 0.1
   ```

2. **Adjust symmetry precision:**
   ```bash
   SYMPREC = 1E-8
   ```

3. **Check structure quality**

4. **Try without symmetry:**
   ```bash
   ISYM = 0
   ```

---

## Memory and Performance

### Out of Memory

**Error:** Killed, segmentation fault, or "allocation failed"

**Solutions:**

1. **Increase NCORE (reduces memory per core):**
   ```bash
   NCORE = 4       # Or higher
   ```

2. **Use LREAL = Auto:**
   ```bash
   LREAL = Auto    # Real-space projection
   ```

3. **Reduce NBANDS:**
   ```bash
   # VASP auto-sets, but can reduce if too many
   NBANDS = 100    # Lower than default
   ```

4. **Run on more nodes:**
   - Distribute memory across more CPUs

5. **Check system resources:**
   ```bash
   # Monitor during run
   top
   free -h
   ```

### Very Slow Calculation

**Symptom:** VASP runs much slower than expected

**Diagnoses and Solutions:**

1. **Too many electronic steps:**
   ```bash
   # Check OSZICAR
   # If many SCF steps (>30), try:
   ALGO = Fast
   # Or adjust AMIX, BMIX
   ```

2. **Inefficient parallelization:**
   ```bash
   # Set NCORE appropriately
   NCORE = 4       # Typical

   # Use KPAR for many k-points
   KPAR = 4
   ```

3. **LREAL = .FALSE. for large system:**
   ```bash
   LREAL = Auto    # Much faster for >100 atoms
   ```

4. **Dense k-mesh:**
   - Reduce k-points for testing
   - Use KSPACING instead:
     ```bash
     KSPACING = 0.5
     ```

5. **Check I/O:**
   - LWAVE = .FALSE., LCHARG = .FALSE.
   - Use fast storage (SSD, local disk not NFS)

---

## Input File Errors

### POSCAR/POTCAR Mismatch

**Error:** Wrong number of electrons, strange results

**Check:**
```bash
# Elements in POSCAR line 6
head -6 POSCAR | tail -1

# Elements in POTCAR
grep VRHFIN POTCAR

# MUST be same order!
```

**Example:**
- POSCAR has: `Cu O`
- POTCAR must have Cu potential first, then O

### KPOINTS Error

**Error:** Wrong KPOINTS format

**Common mistakes:**
1. Line 1: Comment (anything)
2. Line 2: Must be `0` for automatic
3. Line 3: `Gamma` or `Monkhorst-Pack`
4. Line 4: Three integers (k-mesh)

**Correct format:**
```
Automatic mesh
0
Gamma
  8 8 8
  0 0 0
```

### POTCAR Corrupted

**Error:** VASP crashes immediately, "error reading POTCAR"

**Solutions:**
1. **Re-copy POTCAR:**
   ```bash
   cat ~/vasp/potpaw_PBE/Cu/POTCAR > POTCAR
   ```

2. **Check file size:**
   ```bash
   ls -lh POTCAR
   # Should be >100 kB typically
   ```

3. **Verify not truncated:**
   ```bash
   tail POTCAR
   # Should end with "End of Dataset"
   ```

---

## Output Interpretation

### Negative Frequencies (Phonons)

**Symptom:** Negative eigenvalues in phonon calculation

**Meaning:** Imaginary frequencies (system unstable)

**Causes:**
1. **Structure not relaxed:**
   - Tighten EDIFFG = -0.01
   - Ensure EDIFF = 1E-8 for phonons

2. **Numerical noise:**
   - Small imaginary modes (<10 cm⁻¹) at Γ are often numerical
   - Increase ENCUT
   - Tighter convergence

3. **Real instability:**
   - System wants to distort
   - May indicate phase transition

**Solutions:**
- Re-relax with tight settings
- Check if modes at Γ (acoustic, should be zero)
- If real instability, structure may be wrong

### Entropy Term Too Large

**Check in OUTCAR:**
```
  energy  without entropy=   -10.123  energy(sigma->0) =   -10.125
```

**Entropy term:** Difference between these values

**Should be:** < 1 meV/atom

**If too large:**
```bash
# Reduce SIGMA
SIGMA = 0.1     # Smaller

# Or use more k-points
# Or use ISMEAR = -5 (tetrahedron, for static)
```

### Band Gap is Zero (Should be Insulator)

**Causes:**
1. **DFT underestimates gaps** (expected!)
   - Use HSE06 or GW for accurate gaps

2. **k-points too sparse:**
   - Increase k-mesh

3. **ISMEAR wrong:**
   - Use ISMEAR = 0 or -5, not 1

**Solutions:**
```bash
# For accurate gap
# Option 1: HSE06
LHFCALC = .TRUE.
HFSCREEN = 0.2

# Option 2: Denser k-mesh + tetrahedron
ISMEAR = -5
# Very dense KPOINTS
```

---

## Advanced Functionals

### HSE Not Converging

**Symptom:** Many SCF iterations, no convergence

**Solutions:**

1. **Use ALGO = All or Damped:**
   ```bash
   ALGO = Damped
   TIME = 0.5
   ```

2. **Increase NELM:**
   ```bash
   NELM = 200
   ```

3. **Start from PBE:**
   - Run PBE first, save WAVECAR
   - Use WAVECAR as starting point for HSE
   - May help convergence

4. **Adjust mixing:**
   ```bash
   AMIX = 0.2
   BMIX = 0.0001
   TIME = 0.4
   ```

### GW Results Unrealistic

**Symptom:** Band gap too large or negative

**Causes:**
1. **ENCUTGW not converged**
2. **NBANDS too small**
3. **Bad DFT starting point**

**Solutions:**

1. **Converge ENCUTGW:**
   ```bash
   # Test ENCUTGW = 150, 200, 250, 300
   ```

2. **Increase NBANDS:**
   ```bash
   NBANDS = 300    # Try 3-4× occupied bands
   ```

3. **Check DFT step (step 1):**
   ```bash
   # Ensure EDIFF = 1E-8
   # Verify LOPTICS = .TRUE.
   # Check WAVECAR exists
   ```

### DFT+U Multiple Solutions

**Symptom:** Different energies with different MAGMOM

**Explanation:** DFT+U can have multiple metastable states

**Strategy:**
1. **Try various initial MAGMOM:**
   ```bash
   # Ferromagnetic
   MAGMOM = 5 5 5 5 0 0

   # Antiferromagnetic
   MAGMOM = 5 -5 5 -5 0 0
   ```

2. **Choose lowest energy state**

3. **Document which state you use**

---

## Diagnostic Commands

### Check if calculation finished
```bash
grep "reached required accuracy" OUTCAR
# Or
tail OSZICAR
```

### Check convergence during run
```bash
tail -f OSZICAR
```

### Extract final energy
```bash
grep "energy  without" OUTCAR | tail -1
```

### Check forces
```bash
grep "TOTAL-FORCE" OUTCAR -A 50
```

### Check stress tensor
```bash
grep "in kB" OUTCAR | tail -1
```

### Check timing
```bash
grep "Total CPU time" OUTCAR
grep "Elapsed time" OUTCAR
```

### Check memory usage
```bash
grep "maximum memory" OUTCAR
```

### Validate POTCAR order
```bash
head -6 POSCAR | tail -1          # Elements in POSCAR
grep VRHFIN POTCAR                # Elements in POTCAR
# MUST match order!
```

---

## Prevention Checklist

Before running production calculations:

- [ ] Test with quick calculation first (low ENCUT, sparse k-points)
- [ ] Verify POTCAR order matches POSCAR
- [ ] Check structure visually (VESTA, p4vasp)
- [ ] Confirm k-points and ENCUT are converged
- [ ] Set appropriate ISMEAR for system type
- [ ] Use EDIFF = 1E-6 minimum (1E-8 for phonons)
- [ ] Set NSW and NELM high enough
- [ ] Save OUTCAR, vasprun.xml for all runs
- [ ] Monitor OSZICAR during first few steps

---

## Getting Help

1. **VASP Forum:** https://www.vasp.at/forum/
2. **VASP Wiki:** https://vasp.at/wiki/
3. **Search error message:** Often already solved
4. **Provide details:** INCAR, POSCAR, error message, VASP version

---

## Summary

**Most common issues:**
1. POTCAR order wrong → Check element order
2. SCF not converging → ALGO = All, adjust mixing
3. Relaxation stuck → Reduce POTIM, try IBRION = 1
4. Memory error → Increase NCORE, use LREAL = Auto
5. Slow calculation → Check NCORE, LREAL, reduce k-points for testing

**Golden rules:**
- Test first with fast settings
- Check structure visually
- Monitor OSZICAR during run
- Save all outputs
- Converge ENCUT and k-points
