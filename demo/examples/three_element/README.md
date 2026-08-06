# Three-element input-substitution example

This synthetic example is a software-interface test. It deliberately changes
the element labels, count-column names, module IDs, energy-backend IDs, and
row order relative to the manuscript inputs.

Run from the repository root:

```bash
python demo/run_reviewer_demo.py \
  --mode quick \
  --manifest demo/examples/three_element/system_manifest.json
```

The command validates a complete 15-composition, three-element manifold on
four mixing sites and compares two arbitrary matched energy columns. The
numbers are deterministic synthetic fixtures and carry no materials-science
or transferability claim.

To adapt the interface for supplied data:

1. Replace `manifold_energies.csv` and update the count/energy columns in
   `manifold_config.json`.
2. Replace `backend_energies.csv` and update the key, backend IDs, labels, and
   columns in `system_manifest.json`.
3. Keep the bounded claim states unchanged.
4. Treat a passing run as input-contract and internal-validation evidence;
   validate new structures, reference states, energy methods, competing
   phases, interfaces, and experiments separately.
