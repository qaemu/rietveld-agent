# rietveld-agent

Open-source scientific tool that helps crystallographers and laboratory
technicians analyze laboratory powder X-ray diffraction (PXRD) data: validate
inputs against trusted instrument calibrations, retrieve plausible crystalline
phase families from a pinned local catalog, verify hypotheses with bounded
GSAS-II Rietveld refinement, and report evidence-rich, honestly-uncertain
results.

**Status: Phase 0 (foundations). This is a scaffold — no scientific claims yet.**

## Design posture

- **Hybrid scientific workflow, not an LLM over GSAS-II.** GSAS-II is the
  numerical authority for Rietveld refinement; the deterministic core plans and
  executes bounded, policy-approved refinement recipes.
- **Raw data + sample name in; bounded, evidence-rich result out.** The sample
  name is weak context only (it may add `context_only` candidates, never remove
  candidates, never improve scores, never confirm a phase).
- **Defensible uncertainty over false certainty.** Statuses per phase family:
  `supported`, `inconclusive`, `not_selected`, `out_of_domain`, `held`, `failed`.
  No "absent" verdicts without detection-limit studies. No universal accuracy
  numbers. No claims of phase purity, structure solution, publication
  readiness, or robust QPA from Rwp alone.
- **Controlled scientific inputs** (catalog, calibrations, sample-name
  vocabulary, scientific policies, optional models) are hashed, versioned,
  reviewed, replayable, and rollback-able. No update may silently alter the
  meaning of an existing analysis.
- **Local-first, no telemetry, no cloud refinement, Apache-2.0.**
- **AI hosts are interchangeable assistants** (Codex, Claude Code, OpenCode)
  around one deterministic engine, governed by shared versioned contracts.

## Declared initial domain

Crystalline inorganic and mineral powders, constant-wavelength laboratory
X-ray diffraction, Bragg–Brentano geometry, one measured pattern per analysis,
known structures present in the pinned local catalog. Phase-count capability is
a versioned policy field (Phase-0 spike: 2-phase cap; 3 phases experimental).

Organic compounds, strongly disordered/mixed-layer clays, major amorphous
fractions, nanocrystalline-dominated patterns, unknown catalog-absent phases,
and synchrotron/neutron/TOF/capillary/non-Bragg-Brentano data are
out-of-domain or abstention conditions.

## Repository layout

```
core/           deterministic scientific engine (ingest, calibration, catalog,
                retrieval, rerank, hypothesis, verdict, report, governance)
cli/            operator + expert CLI (source of truth for workflow)
admin/          administrator CLI (calibrations, catalog releases, policies)
mcp/            optional deterministic MCP server (same contracts)
skills/         agent skills (setup-rietveld, analyze-raw-pxrd, review-rietveld-run)
                + shared versioned contracts in skills/contracts/
governance/     schemas and policies for controlled scientific inputs
benchmarks/     instrument-aware simulator, evaluation harness, pre-registered thresholds
```

## Phase 0 work (current)

Scaffold, governance schema drafts, GSAS-II environment pinning and
scriptable-API spike, calibration registry spike, COD catalog spike (release
`catalog_0.1.0` in `data/catalog/`, chemistry/SG-validated, manifest-pinned,
rebuilding `data/candidates/library.json`), and an instrument-aware synthetic
simulator spike. Gate G0 checklist: [`docs/G0-CHECKLIST.md`](docs/G0-CHECKLIST.md).

## License

Apache-2.0 for original code. GSAS-II is an external dependency installed via
its official channels (never bundled; see its license). The reference catalog
is a curated, revisioned subset of the Crystallography Open Database
(CC0/public domain, attribution preserved in catalog releases).
