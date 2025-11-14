#' Simulate diploid genotypes from per-population allele frequencies
#'
#' @description
#' This function generates a diploid SNP dataset by sampling genotypes for a
#' specified number of individuals per population from user-provided allele
#' frequencies. The result is returned as an `adegenet::genlight` object with
#' population and individual metadata.
#'
#' @details
#' The input `df` must have three columns: population name, locus name, and the
#' frequency of the first allele for that population–locus
#' combination. For each population, the function simulates two haploid
#' chromosomes per individual by independently drawing alleles at each locus
#' according to the provided allele frequency, then merges the two chromosomes
#' into diploid genotypes (0, 1, 2 copies of the first allele). The procedure
#' assumes Hardy–Weinberg proportions and linkage equilibrium (i.e., loci are
#' sampled independently and there is no within-population structure beyond the
#' supplied allele frequencies).  
#'
#' Sex labels are assigned as "Male"/"Female" in alternating blocks (stored as
#' factors `"m"`/`"f"` in the returned object), and a placeholder phenotype is
#' set to `"control"` for all individuals. Locus allele labels are initialized
#' to `"G/C"` as a placeholder. Computation of chromosomes and genotype strings
#' is implemented with `Rcpp` for speed.
#'
#' @param df A `data.frame` with **three** columns: (1) population name,
#'   (2) locus name, and (3) frequency of the first allele (numeric in \[0, 1\]).
#'   The function internally renames these to `popn`, `locus`, and `frequency`.
#' @param pop.sizes A numeric (integer) vector of population sizes, with one
#'   element **per unique population** in `df`, in the same order as
#'   `unique(df$popn)`.
# @param fbm If TRUE, the genlight object will be converted to a filebacked 
# large matrix format, which is faster if the dataset is large 
# [default FALSE, because still in a testing phase].
# If you want to back convert use gl.gen2fbm and gl.fbm2gen.
#'
#' @return
#' A `genlight` object with:
#' \itemize{
#'   \item diploid SNP genotypes encoded as allele counts (0, 1, 2 for copies of
#'         the first allele),
#'   \item `pop()` set to population names,
#'   \item individual IDs of the form `"0_<popIndex>_<i>"`,
#'   \item `other$ind.metrics` containing `sex` (`"m"`/`"f"`) and `phenotype`
#'         (`"control"`).
#' }
#'
#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @importFrom Rcpp cppFunction
#' @importFrom data.table rbindlist
#' @export
#'
#' @examples
# if (isTRUE(getOption("dartR_fbm"))) platypus.gl <- gl.gen2fbm(platypus.gl)
#' t1 <- gl.filter.callrate(platypus.gl,threshold = 1, mono.rm = TRUE)
#' r1 <- gl.allele.freq(t1, by='popxloc' )
#' r2 <- r1[,c("popn",'locus',"frequency")]
#' res <- gl.sim.ind.af(df = r2, pop.sizes= c(50,50,50))

gl.sim.ind.af <- function(df,
                          pop.sizes
                          # ,
                          # fbm = FALSE
                          ) {
  # --- Basic validation --------------------------------------------------------
  if (!is.data.frame(df) || ncol(df) != 3L) {
    stop("'df' must be a data.frame with exactly three columns: ",
         "population, locus, frequency.")
  }
  colnames(df) <- c("popn", "locus", "frequency")
  
  if (!is.numeric(df$frequency) || anyNA(df$frequency)) {
    stop("'df$frequency' must be numeric with no NAs.")
  }
  if (any(df$frequency < 0 | df$frequency > 1)) {
    stop("All allele frequencies must be in [0, 1].")
  }
  
  # Unique populations in the order they appear
  pops_in_df <- unique(df$popn)
  n_pops <- length(pops_in_df)
  
  # Coerce pop.sizes to integer and align names/order if provided
  if (!is.numeric(pop.sizes)) {
    stop("'pop.sizes' must be numeric (integers).")
  }
  pop.sizes <- as.integer(pop.sizes)
  
  if (length(pop.sizes) != n_pops && is.null(names(pop.sizes))) {
    stop("Length of 'pop.sizes' (", length(pop.sizes),
         ") must equal the number of unique populations in 'df' (", n_pops, "). ",
         "Alternatively, provide a *named* vector whose names match df$popn.")
  }
  
  if (!is.null(names(pop.sizes))) {
    # Reorder sizes to match pops_in_df; check names match exactly
    missing_names <- setdiff(pops_in_df, names(pop.sizes))
    extra_names   <- setdiff(names(pop.sizes), pops_in_df)
    if (length(missing_names)) {
      stop("Missing sizes for populations: ", paste(missing_names, collapse = ", "))
    }
    if (length(extra_names)) {
      stop("Unknown population names in 'pop.sizes': ",
           paste(extra_names, collapse = ", "))
    }
    pop.sizes <- pop.sizes[pops_in_df]
  } else {
    # Unnamed: assume order corresponds to unique(df$popn)
    if (length(pop.sizes) != n_pops) {
      stop("Unnamed 'pop.sizes' must be length ", n_pops, ".")
    }
  }
  
  # --- Ensure consistent locus set & order across populations ------------------
  # Define a canonical locus order from the whole df
  loci_all  <- unique(df$locus)
  n_loci    <- length(loci_all)
  
  # Split per population and (optionally) check each has the same locus set
  df_pops <- split(df, df$popn)
  # if (check.loci) {
    bad <- vapply(df_pops, function(d) {
      !setequal(d$locus, loci_all)
    }, logical(1))
    if (any(bad)) {
      stop("All populations must provide the same set of loci. ",
           "Populations failing this check: ",
           paste(names(df_pops)[bad], collapse = ", "))
    }
  # }
  # Reorder each population's table by canonical locus order
  df_pops <- lapply(df_pops, function(d) {
    d[match(loci_all, d$locus), , drop = FALSE]
  })
  
  # --- Rcpp helpers ------------------------------------------------------------
  # Faster haplotype generator:
  # For j haplotypes, and a vector q of locus-wise allele-1 frequencies,
  # returns a string of '0'/'1' per haplotype (length = length(q))
  make_chr <- function(){}  # dummy for R CMD check
  Rcpp::cppFunction(plugins = "cpp11", code = '
  Rcpp::StringVector make_chr(const int j, const Rcpp::NumericVector q) {
    Rcpp::RNGScope scope;
    const int L = q.size();
    Rcpp::StringVector out(j);
    for (int i = 0; i < j; ++i) {
      std::string s; s.reserve(L);
      for (int z = 0; z < L; ++z) {
        // Bernoulli draw with prob q[z]
        double u = unif_rand();
        char allele = (u < q[z]) ? \'1\' : \'0\';
        s.push_back(allele);
      }
      out[i] = s;
    }
    return out;
  }')
  
  # Convert two haplotype strings into per-locus diploid genotype counts (0,1,2)
  make_geno <- function(){}  # dummy for R CMD check
  Rcpp::cppFunction(plugins = "cpp11", code = '
  Rcpp::List make_geno(const Rcpp::StringMatrix mat) {
    const int ind = mat.nrow();
    const int L   = Rf_length(mat(0,0)); // length of first string
    Rcpp::List out(ind);
    for (int i = 0; i < ind; ++i) {
      std::string chr1 = Rcpp::as<std::string>(mat(i,0));
      std::string chr2 = Rcpp::as<std::string>(mat(i,1));
      Rcpp::IntegerVector temp(L);
      for (int j = 0; j < L; ++j) {
        int a = (chr1[j] == \'1\') ? 1 : 0;
        int b = (chr2[j] == \'1\') ? 1 : 0;
        temp[j] = a + b;
      }
      out[i] = temp;
    }
    return out;
  }')
  
  # --- Simulate per-population haplotypes and build individual table -----------
  pops_idx   <- seq_len(n_pops)
  pop_list   <- vector("list", length = n_pops)
  names(pop_list) <- pops_in_df
  
  for (i in pops_idx) {
    pop_name   <- pops_in_df[i]
    n_ind      <- pop.sizes[i]
    if (n_ind < 1L) stop("Population '", pop_name, "' has size < 1.")
    
    # Haplotypes for all individuals in this population (2 per individual)
    q_vec   <- df_pops[[i]][, "frequency"]
    chr_all <- make_chr(j = n_ind * 2L, q = q_vec)
    
    # Split into two haplotype sets (first/second chromosome) by alternating index
    grp     <- rep(1:2, length.out = length(chr_all))
    chr_pops <- split(chr_all, grp)
    
    # Create an individual table
    pop_df <- data.frame(
      SEX  = rep(c("Male", "Female"), length.out = n_ind),   # robust to odd sizes
      POP  = rep(pop_name, n_ind),
      CHR1 = chr_pops[[1]],
      CHR2 = chr_pops[[2]],
      id   = sprintf("0_%s_%d", pop_name, seq_len(n_ind)),
      stringsAsFactors = FALSE
    )
    
    pop_list[[i]] <- pop_df
  }
  
  # Bind all populations
  df_geno <- data.table::rbindlist(pop_list, use.names = TRUE, fill = TRUE)
  
  # --- Build genotype matrix (0/1/2) via Rcpp and convert to genlight ----------
  # Two-string haplotype matrix
  hap_mat <- as.matrix(df_geno[, c("CHR1", "CHR2")])
  
  # Per-individual list of integer vectors (length = n_loci)
  ped_list <- make_geno(hap_mat)
  
  # Convert list of vectors to adegenet::SNPbin objects
  txt <- lapply(ped_list, function(e) suppressWarnings(as.integer(e)))
  snp_list <- lapply(txt, function(e) new("SNPbin", snp = e, ploidy = 2L))
  
  # Create genlight
  gl <- new("genlight", snp_list, ploidy = 2L)
  indNames(gl) <- df_geno$id
  pop(gl)      <- df_geno$POP
  locNames(gl) <- as.character(loci_all)
  
  # --- Individual metrics (sex, phenotype placeholders) -----------------------
  # Store metrics in genlight@other$ind.metrics
  ind_metrics <- data.frame(
    fid       = df_geno$POP,
    iid       = df_geno$id,
    sex       = factor(ifelse(df_geno$SEX == "Male", "m", "f")),
    phenotype = factor(rep("control", nrow(df_geno))),
    stringsAsFactors = FALSE
  )
  gl$other$ind.metrics <- ind_metrics
  
  # Placeholder locus allele labels
  gl$loc.all <- rep("G/C", nLoc(gl))
  
  # Reset internal flags if available in your environment (safe to skip if absent)
  if (exists("utils.reset.flags", mode = "function")) {
    gl <- utils.reset.flags(gl, verbose = 0)
  }
  
  # if (fbm) gl <- gl.gen2fbm(gl)
  
  return(gl)
  

}
