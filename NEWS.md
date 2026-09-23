# dartR.sim 1.2.2.9000

## gl.sim.WF.table

* Arguments passed through `...` no longer put quotes around
  `chromosome_name`; before, any `...` override made `real_loc = TRUE` and
  the map/targets chromosome lookups fail unless `chromosome_name` was also
  passed.
* `real_freq = TRUE` with `real_loc = FALSE` now returns loci typed `"real"`.
  Before, it stopped with an error or returned those loci with `NA` values.
* Recombination maps with intervals of any size are now supported. Each
  interval's own `from`/`to` is used, and neutral loci are spread over the
  whole mapped chromosome. **Output changes** for maps whose intervals are not
  aligned to `chunk_bp` from position 1 (including `fly_recom_map.csv`).
  `NA` in `cM` now marks an interval as not recombining instead of breaking
  the call.
* `loci_deleterious >= chunk_number` now returns exactly the number of loci
  requested (before, 149 gave 100 and 150 gave 200). **Output changes** when
  `loci_deleterious` is not a multiple of `chunk_number`.
* Caps on `q` (0.5) and `s` (0.99, -0.5) now apply to each locus class
  according to that class's own distribution setting, are documented, and
  are reported at `verbose >= 1`. **Output changes** when deleterious and
  advantageous settings differ.
* Names passed through `...` that are not simulation variables now stop the
  function with an error. Before, they were ignored.
* Errors now carry their message (`stop(error(...))`).
* `fly_recom_map.csv` and `fly_targets_of_selection.csv` now ship in
  `inst/extdata`.
