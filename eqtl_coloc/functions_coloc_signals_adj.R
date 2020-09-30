##############################################################
#
# Modified version of `coloc.signal`:
#   - coloc.signal -> finemap.signals -> map.mask works with p-values too
#
##############################################################

# coloc.signals_adj --------------------------------
# Add funcionality to use p-values too for method = mask
# Uses `finemap.signals_adj`

coloc.signals_adj <- function(dataset1, dataset2, MAF = NULL, LD = NULL, method = c("single", "cond", "mask"), mode = c("iterative", "allbutone"), p1 = 1e-04, p2 = 1e-04, p12 = NULL, maxhits = 3, r2thr = 0.01, pthr = 1e-06)
{
  method <- match.arg(method)
  mode <- match.arg(mode)
  if (is.null(p12))
    stop("default value for p12 has been removed. please read ... and choose a value appropriate for your study")
  if (!("MAF" %in% names(dataset1)) & !is.null(MAF))
    dataset1$MAF <- MAF
  if (!("MAF" %in% names(dataset2)) & !is.null(MAF))
    dataset2$MAF <- MAF
  if (!("LD" %in% names(dataset1)) & !is.null(LD))
    dataset1$LD <- LD
  if (!("LD" %in% names(dataset2)) & !is.null(LD))
    dataset2$LD <- LD
  if (!("method" %in% names(dataset1)) & !is.null(method))
    dataset1$method <- method
  if (!("method" %in% names(dataset2)) & !is.null(method))
    dataset2$method <- method
  if (dataset1$method == "cond") {
    coloc::check.dataset(dataset1, 1, req = "MAF")
    coloc:::check.ld(dataset1, dataset1$LD)
  }
  else {
    coloc::check.dataset(dataset1, 1)
  }
  if (dataset2$method == "cond") {
    coloc::check.dataset(dataset2, 2, req = "MAF")
    coloc:::check.ld(dataset2, dataset2$LD)
  }
  else {
    coloc::check.dataset(dataset2, 2)
  }
  zthr = qnorm(pthr/2, lower.tail = FALSE)
  fm1 <- finemap.signals_adj(dataset1, method = dataset1$method, maxhits = maxhits, r2thr = r2thr, pthr = pthr)
  fm2 <- finemap.signals_adj(dataset2, method = dataset2$method, maxhits = maxhits, r2thr = r2thr, pthr = pthr)
  if (is.null(fm1))
    fm1 <- find.best.signal_adj(dataset1)
  if (is.null(fm2))
    fm2 <- find.best.signal_adj(dataset2)
  if (!is.null(fm1) && dataset1$method == "cond") {
    cond1 <- coloc:::est_all_cond(dataset1, fm1, mode = mode)
    X1 <- dataset1[intersect(names(dataset1), c("N", "sdY", "type", "s"))]
  }
  if (!is.null(fm2) && dataset2$method == "cond") {
    cond2 <- coloc:::est_all_cond(dataset2, fm2, mode = mode)
    X2 <- dataset2[intersect(names(dataset2), c("N", "sdY", "type", "s"))]
  }
  if (dataset1$method %in% c("single", "mask") & dataset2$method %in% c("single", "mask")) {
    col <- coloc:::coloc.detail(dataset1, dataset2, p1 = p1, p2 = p2, p12 = p12)
    res <- coloc:::coloc.process(col, hits1 = names(fm1), hits2 = names(fm2),
                         LD1 = dataset1$LD, LD2 = dataset2$LD, r2thr = r2thr,
                         mode = mode, p12 = p12)
  }
  if (dataset1$method == "cond" & dataset2$method == "cond") {
    todo <- expand.grid(i = seq_along(cond1), j = seq_along(cond2))
    res <- vector("list", nrow(todo))
    for (k in 1:nrow(todo)) {
      col <- coloc:::coloc.detail(dataset1 = c(cond1[[todo[k, "i"]]], X1), dataset2 = c(cond2[[todo[k, "j"]]], X2),
                          p1 = p1, p2 = p2, p12 = p12)
      res[[k]] <- coloc:::coloc.process(col, hits1 = names(fm1)[todo[k, "i"]], hits2 = names(fm2)[todo[k, "j"]], LD1 = dataset1$LD,
                                LD2 = dataset2$LD, p12 = p12)
    }
  }
  if (dataset1$method == "cond" && dataset2$method %in% c("single", "mask")) {
    res <- vector("list", length(cond1))
    for (k in 1:length(res)) {
      col <- coloc:::coloc.detail(c(cond1[[k]], X1), dataset2, p1 = p1, p2 = p2, p12 = p12)
      res[[k]] <- coloc:::coloc.process(col, hits1 = names(fm1)[k],
                                hits2 = names(fm2), LD1 = dataset1$LD, LD2 = dataset2$LD,
                                r2thr = r2thr, p12 = p12)
    }
  }
  if (dataset1$method %in% c("single", "mask") && dataset2$method == "cond") {
    res <- vector("list", length(cond2))
    for (k in 1:length(res)) {
      col <- coloc:::coloc.detail(dataset1, c(cond2[[k]], X2), p1 = p1, p2 = p2, p12 = p12)
      res[[k]] <- coloc:::coloc.process(col, hits1 = names(fm1),
                                hits2 = names(fm2)[k], LD1 = dataset1$LD, LD2 = dataset2$LD,
                                r2thr = r2thr, p12 = p12)
    }
  }
  if (dataset1$method == "cond" || dataset2$method == "cond") {
    summary <- data.table::rbindlist(lapply(res, "[[", "summary"))
    rowvars <- c("SNP.PP.H4", "z.df1", "z.df2")
    results <- res[[1]]$results
    data.table::setnames(results, rowvars, paste0(rowvars, ".row1"), skip_absent = TRUE)
    for (i in setdiff(1:length(res), 1)) {
      thisone <- res[[i]]$results[, grep("snp|SNP.PP.H4|z.df1|z.df2", names(res[[i]]$results)), drop = FALSE, with = FALSE]
      data.table::setnames(thisone, rowvars, paste0(rowvars, ".row", i), skip_absent = TRUE)
      results <- merge(results, thisone, by = "snp", all = TRUE)
    }
  } else {
    summary <- res[["summary"]]
    results <- res[["results"]]
  }
  d1 <- data.table::data.table(hit1 = names(fm1), hit1.margz = c(fm1))
  res <- merge(summary, d1, by = "hit1")
  d2 <- data.table::data.table(hit2 = names(fm2), hit2.margz = c(fm2))
  res <- merge(res, d2, by = "hit2")
  ret <- list(summary = res, results = results, priors = c(p1 = p2, p2 = p2, p12 = p12))
  class(ret) <- c("coloc_abf", class(ret))
  ret
}

# finamap.signals_adj --------------------------------
# Add funcionality to use p-values too for method = mask
# Uses `map_mask_adj`

finemap.signals_adj <- function(D, LD = D$LD, method = c("single", "mask", "cond"), r2thr = 0.01, sigsnps = NULL, pthr = 1e-06, maxhits = 3)
{
  method <- match.arg(method)
  if (method == "cond") {
    coloc::check.dataset(D, req = "MAF")
    coloc:::check.ld(D, LD)
  } else {
    coloc::check.dataset(D)
  }
  if (!is.null(sigsnps) && is.null(names(sigsnps)))
    stop("sigsnps should be a named numeric vector, with snp ids as the names")
  if (!all(names(sigsnps) %in% D$snp))
    stop("not all sigsnps found in D$snp")
  zthr = qnorm(pthr/2, lower.tail = FALSE)
  if (D$type == "cc" & method == "cond") {
    message("approximating linear analysis of binary trait")
    D <- coloc:::bin2lin(D)
  }
  YY = if (D$type == "quant") {
    sum(D$N * D$sdY^2)
  } else {
    sum(D$N * D$s * (1 - D$s))
  }
  LD <- LD[D$snp, D$snp]
  hits <- NULL
  while (length(hits) < maxhits) {
    newhit = if (method == "mask") {
      map_mask_adj(D, LD, r2thr, sigsnps = names(hits))
    } else {
      coloc:::map_cond(D = D, LD = LD, YY = YY, sigsnps = names(hits))
    }
    if (is.null(newhit) || !length(newhit) || abs(newhit) < zthr)
      break
    hits <- c(hits, newhit)
    if (method == "single")
      break
  }
  hits
}

# map_mask_adj ------------------------------------------------
# Add funcionality to use p-values too

map_mask_adj <- function(D, LD, r2thr = 0.01, sigsnps = NULL)
{
  x <- data.table::as.data.table(D[c("beta", "varbeta", "pvalues", "snp", "MAF")])
  if ("pvalues" %in% colnames(x)) {
    x[, `:=`(z, qnorm(pvalues/2, lower.tail = FALSE))]
  } else {
    x[, `:=`(z, beta/sqrt(varbeta))]
  }
  use <- rep(TRUE, nrow(x))
  if (!is.null(sigsnps)) {
    expectedz <- rep(0, nrow(x))
    friends <- apply(abs(LD[, sigsnps, drop = FALSE]) > sqrt(r2thr), 1, any, na.rm = TRUE)
    use <- use & !friends
    if (!any(use))
      return(NULL)
    imask <- match(sigsnps, x$snp)
    expectedz <- LD[, sigsnps, drop = FALSE] %*% x$z[imask]
    zdiff <- abs(x$z[use]) - abs(expectedz[use])
  } else {
    zdiff <- abs(x$z)
  }
  wh <- which.max(zdiff)
  structure(x$z[use][wh], names = x$snp[use][wh])
}

# find.best.signal_adj ------------------------------------------------
# Add funcionality to use p-values too

find.best.signal_adj <- function(D)
{
  if (D$type == "cc" & D$method == "cond") {
    message("approximating linear analysis of binary trait")
    D <- coloc:::bin2lin(D)
  }
  YY = if (D$type == "quant") {
    sum(D$N * D$sdY^2)
  } else {
    sum(D$N * D$s * (1 - D$s))/sum(D$N)
  }
  z = if ("pvalues" %in% names(D)) {
    qnorm(D$pvalues/2, lower.tail = FALSE)
  } else {
    D$beta/sqrt(D$varbeta)
  }
  zdiff <- abs(z)
  wh <- which.max(zdiff)
  structure(z[wh], names = D$snp[wh])
}
