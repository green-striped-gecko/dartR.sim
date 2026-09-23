# Review: gl.sim.offspring (dartR.sim)

- Family mode: analysis (simulates Mendelian offspring from parental genotypes)
- Date: 2026-09-23
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 94a72f0 (`origin/dev`)
- Datasets: `testset.gl` (all-`NA` loci removed, missing data imputed at random); `testset.gs`; synthetic single parents with chosen heterozygous loci
- Baseline: `tests/testthat/test-gl.sim.offspring.R` (snapshot captured pre-review, 13 expectations, all pass)

## Verdict

**Standards: Needs work** — parental inputs are not checked, the output is a bare genlight without locus metadata, and sex is stored outside `ind.metrics`.
**Spec: Needs work** — alleles at different loci are not always inherited independently: random draws are reused across loci, so some pairs of loci are passed on together. `noffpermother` is not the number of offspring per mother.

What works: single-locus Mendelian segregation (heterozygote × heterozygote gives 0.251 / 0.499 / 0.250) and the sex ratio (0.304 for 0.3).

## Findings

**F1 [HIGH, confidence: high] — random draws are recycled, so some loci are inherited together (principle: model correctness)**
`R/gl.sim.offspring.r:55–66` — `ifelse(mmat == 1, sample(c(0, 2), mhet, replace = TRUE), mmat)` draws one value per heterozygous cell, but `ifelse()` recycles that vector over the whole matrix by position. Heterozygous cells whose positions are equal modulo the number of draws get the same draw, so their alleles are inherited together.
Failure scenario:
- A mother heterozygous at loci 1 and 11 only (20 loci) passes on the same allele at both loci in 100% of 2000 offspring; independent inheritance gives 50%.
- A `testset.gl` mother with 13 identifiable heterozygous loci and 1000 offspring: 4 of 78 locus pairs (5.1%) are perfectly linked.
- `dartR.captive::gl.sim.relatedness()` builds full sibs, half sibs and cousins with this function, so their simulated relatedness carries this artificial linkage.
Proposed change: draw only for heterozygous cells (`mmat[het] <- sample(c(0, 2), sum(het), replace = TRUE)`), and the same for fathers. **Consequence: seeded outputs change; loci become independent.**

**F2 [MEDIUM, confidence: high] — `noffpermother` is not the number of offspring per mother (DOC5)**
`R/gl.sim.offspring.r:48–50` — `nInd(mothers) * noffpermother` offspring are made, but each offspring's mother is drawn at random with replacement.
Failure scenario: 5 mothers, `noffpermother = 4`: the mothers get 3, 5, 1, 5 and 6 offspring. A mother can end up with none, so sibship sizes are not what the argument says.
Proposed change: `rep(seq_len(nInd(mothers)), each = noffpermother)`. Fathers stay random per offspring (random mating); document that siblings are full sibs only when one father is given. **Consequence: family sizes are exact; seeded outputs change.**

**F3 [MEDIUM, confidence: high] — parental inputs are not checked (FS4, FS5, DAT1)**
`R/gl.sim.offspring.r:40–50`
- Fathers with the same loci in a different order run silently, and offspring genotypes combine different loci under the mothers' locus names.
- Fathers with fewer loci give "non-conformable arrays".
- SilicoDArT parents give diploid 0/1/2 offspring.
- `noffpermother = 0` gives "subscript out of bounds".
- `sexratio` is not checked to be between 0 and 1.
Proposed change: `utils.check.datatype(accept = "SNP")` for both; stop unless `locNames(fathers)` is identical to `locNames(mothers)`; validate `noffpermother` (whole number ≥ 1) and `sexratio` (0–1).

**F4 [MEDIUM, confidence: high] — offspring lose metadata and parentage (DAT2, DAT5)**
`R/gl.sim.offspring.r:69–84` — the output is a plain `genlight` with empty `loc.all`, no position, chromosome or `loc.metrics`, and no `ind.metrics`. Sex is in `@other$sex`, and the parents of each offspring are not recorded.
Failure scenario: offspring cannot be analysed per chromosome and need `gl.compliance.check()` first; pedigree-based checks (relatedness against true parents) cannot tell which parents produced each offspring.
Proposed change: return a `dartR` object carrying `loc.all`, position, chromosome and `loc.metrics` from the mothers (flags reset), with `ind.metrics` holding `id`, `pop`, `sex` (`"female"`/`"male"`), `mother` and `father` (parent names). Keep `@other$sex` for back-compatibility.

**F5 [LOW, confidence: high] — standard skeleton (FS3, FS8, FS9)**
`R/gl.sim.offspring.r:40–86` — no `utils.flag.start()`, no completion message, no history entry. `if (!is.na(mhet))` is always true (it guards nothing). The missing-data warning has a typo ("you offspring"), and the `fbm` leftovers remain.
Proposed change: add the FS3/FS9 lines and history; drop the dead guard and `fbm` leftovers; fix the typo.

**F6 [LOW, confidence: high] — documentation (DOC1, DOC5, DOC7 proposed)**
`R/gl.sim.offspring.r:1–34` — the mating scheme is unclear: each offspring has a random father, so siblings are half sibs unless one father is given. `@param mothers` says "potential mothers simulated". `@return` says "n individuals". There is no `@family`, and `@author` lacks the Author(s)/Custodian parts.
Proposed change: describe the model (Mendelian segregation, independent loci, one random father per offspring) and the output metadata; fix the tags.

## Proposed changes

1. Independent draws per heterozygous genotype (F1). **Consequence: seeded outputs change; artificial linkage between loci removed.**
2. Exactly `noffpermother` offspring per mother (F2). **Consequence: family sizes change to the stated number; seeded outputs change.**
3. Check inputs: SNP data, identical loci in fathers and mothers, `noffpermother`, `sexratio` (F3). **Consequence: parents with mismatched loci or SilicoDArT data now stop with an error.**
4. Carry metadata; add `ind.metrics` with sex and parents; return a `dartR` object; keep `@other$sex` (F4).
5. Standard skeleton, history, clean-up (F5).
6. Documentation (F6).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run. PLT, DEP: not applicable.
- Spec: Mendelian proportions, sex ratio, independence across loci (synthetic and `testset.gl` parents), offspring per mother, locus mismatch, SilicoDArT, `noffpermother = 0`, missing data, `rbind()` with parents — run.
- Callers: `dartR.captive::gl.sim.relatedness()` (6 calls with one father and one mother; affected by F1, not by F2) and the dartr2shiny `code_args.csv` entry (reads `ind.metrics$sex` of the input) — read.
- FBM path (DAT6): not applicable (the `fbm` argument is disabled).
- Google Group / GitHub issues: not searched.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis |  |
| 2 | approved | Luis |  |
| 3 | approved | Luis |  |
| 4 | approved | Luis |  |
| 5 | approved | Luis |  |
| 6 | approved | Luis |  |

## Addendum notes (not changed)

- `utils.check.datatype()` (dartR.base) treats a SNP object whose genotypes are all 0/1 and that has no `loc.all`/SNP metrics as SilicoDArT. A single parent built by hand without alleles, and with no homozygous-alternative genotypes, is therefore rejected by change 3. Real dartR objects carry `loc.all` and pass; the synthetic test parents now carry alleles. This is a dartR.base heuristic, not changed here.
- `rbind()` of dartR objects drops `@other$ind.metrics` (even `rbind(platypus.gl, platypus.gl[1, ])`), before and after this change; outside this function.

## Outcome

| Change | Evidence | Snapshot result |
|---|---|---|
| 1 | test "loci are inherited independently" | flipped: loci 1 and 11 agree in 50% of offspring (was 100%); 0 perfectly linked pairs for a `testset.gl` mother (was 4 of 78) |
| 2 | test "each mother has noffpermother offspring" | flipped: 4, 4, 4, 4, 4 (was 3, 5, 1, 5, 6) |
| 3 | test "inputs are checked" | flipped: locus mismatch, SilicoDArT, `noffpermother = 0` and `sexratio = 2` stop with named errors |
| 4 | test "metadata and parentage are carried"; `glSim()` parents without names get positions as parent IDs | flipped: `dartR` object, `loc.all`/position/`loc.metrics` from mothers, `ind.metrics` with sex/mother/father; `rbind` with parents works (`dartR.captive` pattern) |
| 5 | history entry; start/end messages at `verbose >= 1` | new |
| 6 | `devtools::document()`; `man/gl.sim.offspring.Rd` regenerated | new `@details` |

Mendelian proportions (0.25 / 0.50 / 0.25) and sex ratio unchanged and passing. Tests: `test-gl.sim.offspring.R` 21 expectations pass; full suite passes. NEWS.md updated.
PR: pending.

```json
{
  "function": "gl.sim.offspring",
  "package": "dartR.sim",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "94a72f0",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "principle:model-correctness", "status": "approved", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "FS4", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DAT2", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "FS3", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 6}
  ],
  "coverage_skipped": ["FBM path: fbm argument disabled", "forum/issues search: not run"],
  "status": "applied",
  "pr": null
}
```
