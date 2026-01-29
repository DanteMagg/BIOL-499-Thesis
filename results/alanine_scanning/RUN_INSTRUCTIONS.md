# Running the Alanine Scanning Pipeline

## Quick Start

```bash
cd /Users/DantesFolder/BIOL-499-Thesis-3/results/alanine_scanning
/Library/Frameworks/Python.framework/Versions/3.12/bin/python3 alanine_scan_ddg.py
```

## Expected Runtime

- **24 structures** (WT + 23 mutants)
- **5 replicates each** = 120 total calculations
- **~5-10 minutes per replicate** (depends on CPU)
- **Total: ~10-20 hours** (overnight run recommended)

## What Gets Generated

1. **relaxed_structures/** - All 120 relaxed PDB files (24 structures × 5 replicates)
2. **energy_log.csv** - Raw MM/GBSA energies for every replicate
3. **alanine_scan_results.csv** - Final ΔΔG values and classifications
4. **alanine_scan_heatmap.png** - Figure 1 visualization
5. **run_log.txt** - Complete execution log with timestamps

## Monitoring Progress

The script logs to both console and `run_log.txt`. You can monitor progress:

```bash
tail -f run_log.txt
```

## If Interrupted

The script processes structures sequentially. If interrupted, you can:
1. Check `run_log.txt` to see which structures completed
2. Modify the script to skip completed structures (add a check for existing relaxed PDBs)
3. Or simply restart - it will overwrite previous results

## Output Files Reference

- **METHODOLOGY.md** - Complete protocol documentation
- **mutant_summary.txt** - List of all mutants generated
- **mutants/** - Original 23 mutant PDB files

