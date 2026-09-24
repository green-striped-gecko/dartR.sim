# Review: gl.sim.create_dispersal (dartR.sim)

- Family mode: io (writes the dispersal table read by `gl.sim.WF.run()`)
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: cc29b1f (`origin/dev`); dartR.base f9f1be8 (`origin/dev`)
- Datasets: none (the function takes no genlight input); tables for 0–5 populations and each `dispersal_type`, read back with `read.csv()` as `gl.sim.WF.run()` does
- Baseline: `tests/testthat/test-gl.sim.create_dispersal.R` (snapshot captured pre-review, 18 expectations, all pass)

## Verdict

**Standards: Needs work** — no argument is checked, so invalid values give opaque errors or silently write tables that `gl.sim.WF.run()` cannot use; messages go through `message()`; the function returns `NULL`.
**Spec: Ready** — for valid inputs the table matches the documentation and `gl.sim.WF.run()`: each connected pair appears once, with the columns and values it reads. One documentation pointer is wrong.

## Findings

**F1 [MEDIUM, confidence: high] — arguments are not checked (FS5)**
`R/gl.sim.create.dispersal.r:42–97`
- `number_pops = 1`: "arguments imply differing number of rows" for `all_connected` and `line`; `circle` writes population 1 paired with itself.
- `number_pops = 0` writes a pair between populations 1 and 0; `number_pops = 2.5` becomes 2 populations without a message.
- `dispersal_type = "ring"` (or "Line"): "object 'dispersal_pairs' not found".
- `transfer_each_gen = 0` is written; `gl.sim.WF.run()` then evaluates `gen %% 0` (`NA`) in an `if()` and stops mid-simulation.
- Negative or fractional `number_transfers` are written as given.
Proposed change: `match.arg()` for `dispersal_type`; stop unless `number_pops` is a whole number ≥ 2, `number_transfers` a whole number ≥ 0 (0 lets a user switch a pair off), and `transfer_each_gen` a whole number ≥ 1. **Consequence: `number_pops` of 0, 1 or fractional, and `transfer_each_gen = 0`, which write a file today, stop with an error.**

**F2 [LOW, confidence: high] — output folder not checked (FS7)**
`R/gl.sim.create.dispersal.r:106–113` — a folder that does not exist gives "cannot open the connection".
Proposed change: fall back to `tempdir()` with a warning at `verbose >= 1`; the path actually used is reported at `verbose >= 2`.
(Changed in Phase C: the first draft proposed `gl.check.wd(outpath, verbose = verbose)`, which gives the same fallback but prints its own start/end messages inside this function's; the check is written directly instead.)

**F3 [LOW, confidence: high] — messages (VRB2, FS3)**
`R/gl.sim.create.dispersal.r:54–56, 99–104, 117–119` — messages use `message()` rather than `cat()`, so they go to stderr unlike other dartR functions; the saved path is printed with a trailing "/" (`file.path(outpath, outfile, "\n")`); `utils.flag.start()` gets the outdated `build = "Jody"`, printed at `verbose = 5`.
Proposed change: `cat(report(...))`; print the path without the trailing "/"; drop `build`.

**F4 [LOW, confidence: high] — returns `NULL` (FS10)**
`R/gl.sim.create.dispersal.r:121` — the value is the `NULL` returned by `write.table()`.
Proposed change: return the table invisibly, so it can be inspected or edited in R before writing it again. Printing is unchanged (invisible).

**F5 [LOW, confidence: high] — documentation (DOC1, DOC6 proposed, DOC7 proposed)**
`R/gl.sim.create.dispersal.r:1–40` — the vignette pointer names package "dartR" and closes with a curly quote (`"dartR”`); `@return` does not mention a returned value; `@author` lacks the Author(s)/Custodian parts; `@details` is missing and the argument checks of F1 are not stated.
Proposed change: point to `dartR.sim`, ASCII quotes; document the returned table and the valid ranges; fix the tags.

## Proposed changes

1. Check `number_pops`, `dispersal_type`, `number_transfers`, `transfer_each_gen` (F1). **Consequence: invalid values that write a file today stop with an error.**
2. Check the output folder with `gl.check.wd()` (F2).
3. `cat(report())` messages; path printed correctly; `build` dropped (F3).
4. Return the table invisibly (F4).
5. Documentation (F5).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run. DAT, DEP, PLT: not applicable.
- Spec: tables for 0, 1, 2, 3, 5 populations × three dispersal types; column names and values against the code in `gl.sim.WF.run()` that reads the file; invalid type; zero and fractional values; missing folder; return value; messages — run.
- Callers: `gl.sim.WF.run()` reads the file through `file_dispersal` (`read.csv()`; columns `pop1`, `pop2`, `number_transfers`, `transfer_each_gen`); the Shiny app in `utils.sims.r` points users to this function in a tooltip; no calls in sibling `dartR.*` packages or dartr2shiny.
- Input-data checks (FS4, DAT1–DAT4): not applicable; no genlight input.
- Google Group / GitHub issues: not searched.

### Note outside this function (not proposed here)

`gl.sim.WF.run()` decides which sexes migrate from the phase setting (`number_transfers_phase1/2`), not from each row: `number_transfers >= 2` moves both sexes, `1` moves males only (`R/gl.sim.WF.run.r:403–410, 515–520`). The per-row value then sets the count. A table whose rows were raised to 3 while the phase is set to 1 moves only the male share (2 per direction) and no females. This function's documentation says the column "can be modified by hand"; the interaction belongs to `gl.sim.WF.run()` and needs its own change there.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | consequence approved |
| 2 | approved | Luis | implemented without `gl.check.wd()` (see F2) |
| 3 | approved | Luis |  |
| 4 | approved | Luis |  |
| 5 | approved | Luis |  |

## Outcome

| Change | Evidence | Snapshot result |
|---|---|---|
| 1 | test "invalid arguments stop" | flipped: `number_pops` 1, 0, 2.5, `"ring"`, `transfer_each_gen = 0`, `number_transfers` -1 or 1.5 stop with named errors (were opaque errors, a 1-1 self-pair, a pair with population 0, truncation, or a written `0`). `number_transfers = 0` still writes |
| 2 | test "missing outpath falls back to tempdir()" | flipped: warning at `verbose = 1` and the file is in `tempdir()` (was "cannot open the connection") |
| 3 | test "returns the table invisibly; cat messages" | flipped: path printed through `cat()` without a trailing "/"; silent at `verbose = 0`; no "Build = Jody" line |
| 4 | same test | flipped: returns the table invisibly (was `NULL`); equal to the file read back |
| 5 | `devtools::document()`; `man/gl.sim.create_dispersal.Rd` regenerated | docs only |

Unchanged and passing: pairs for `all_connected`, `line` and `circle`; columns and values read by `gl.sim.WF.run()`; the `gl.sim.WF.run()` test that runs a simulation on a table from this function. Tests: `test-gl.sim.create_dispersal.R` 24 expectations pass; full suite 223 pass. `R CMD check` on the tracked files: 0 errors; 1 warning (installed packages built under R 4.4.3, local environment); 1 note (timestamps). NEWS.md updated.
PR: pending.

```json
{
  "function": "gl.sim.create_dispersal",
  "package": "dartR.sim",
  "family": "io",
  "skill_version": "2.0.0",
  "commit": "cc29b1f",
  "verdict_standards": "needs_work",
  "verdict_spec": "ready",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "LOW", "confidence": "high", "rule": "FS7", "status": "approved", "change": 2},
    {"id": "F3", "severity": "LOW", "confidence": "high", "rule": "VRB2", "status": "approved", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "FS10", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "DOC6", "status": "approved", "change": 5}
  ],
  "coverage_skipped": ["input-data checks: no genlight input", "forum/issues search: not run"],
  "status": "pr-open",
  "pr": null
}
```
