#!/usr/bin/env Rscript
# Five external resamples matching bootstrap_abundance_refits.R exactly.
.libPaths(c(file.path(getwd(), "r_library"), .libPaths()))
suppressPackageStartupMessages(library(OneNet))
options(future.globals.maxSize = 4 * 1024^3)
root <- normalizePath(getwd())
systems <- c("butyrate_assembly_2021", "carlstrom_phyllosphere_2019", "schafer_phyllosphere_2022")
helpers <- c("prepare_count_like", "edge_table", "adapt_mean_stability_safe")
for (expr in parse(file.path(root, "scripts/run_onenet_benchmarks.R"))) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      as.character(expr[[2]])[1] %in% helpers) eval(expr)
}
canonical <- function(d) {
  a <- as.character(d$taxon_1); b <- as.character(d$taxon_2)
  d$taxon_1 <- pmin(a,b); d$taxon_2 <- pmax(a,b); d
}
ap <- function(y,s) {
  ix <- order(s, decreasing=TRUE); y <- y[ix]; s <- s[ix]
  ends <- c(which(diff(s)!=0), length(s))
  hit <- cumsum(y)[ends]; recall <- hit/sum(y)
  sum(diff(c(0,recall))*hit/ends)
}
outdir <- file.path(root,"results/onenet_bootstrap_5")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
tasks <- expand.grid(bootstrap=1:5, system=systems, stringsAsFactors=FALSE)
run <- function(i) {
  system <- tasks$system[i]; b <- tasks$bootstrap[i]
  stem <- file.path(outdir,paste0(system,"_",b))
  if (file.exists(paste0(stem,"_metrics.csv"))) {
    previous <- read.csv(paste0(stem,"_metrics.csv"))
    if (all(previous$fit_status == "complete")) return(previous)
    file.copy(paste0(stem,"_metrics.csv"),paste0(stem,"_previous_failure.csv"),overwrite=TRUE)
  }
  started <- Sys.time()
  log <- file(paste0(stem,".log"),open="wt")
  sink(log); sink(log,type="message")
  on.exit({sink(type="message");sink();close(log)}, add=TRUE)
  result <- tryCatch({
    abundance <- read.csv(file.path(root,"cleaned_data",paste0(system,"_abundance.csv")),check.names=FALSE)
    set.seed(20260825L + match(system,systems)*100L+b)
    rows <- sample.int(nrow(abundance),nrow(abundance),replace=TRUE)
    write.csv(data.frame(sampled_row=rows),paste0(stem,"_rows.csv"),row.names=FALSE)
    sampled <- abundance[rows,-1,drop=FALSE]
    variable <- vapply(sampled,function(x)length(unique(x))>1,logical(1))
    counts <- prepare_count_like(sampled[,variable,drop=FALSE])
    fitseed <- 20260821L+match(system,systems)*1000L+b
    set.seed(fitseed)
    inference <- all_inferences_new(data=counts,rep.num=30,
      methods=c("PLNnetwork","SpiecEasi","gCoda","EMtree","Magma","SPRING"),
      parallel=FALSE,cores=1,seed=fitseed)
    saveRDS(inference,paste0(stem,"_inference.rds"))
    adapted <- adapt_mean_stability_safe(inference,mean_stability=0.8)
    scores <- canonical(edge_table(colnames(counts),compute_aggreg_measures(adapted$freqs)$mean))
    allpairs <- as.data.frame(t(combn(names(abundance)[-1],2)),stringsAsFactors=FALSE)
    names(allpairs) <- c("taxon_1","taxon_2")
    scores <- merge(canonical(allpairs),scores,all.x=TRUE)
    scores$onenet_score[is.na(scores$onenet_score)] <- 0
    write.csv(scores,paste0(stem,"_scores.csv"),row.names=FALSE)
    truth <- canonical(read.csv(file.path(root,"cleaned_data",paste0(system,"_tested_pairs.csv"))))
    evaluated <- merge(truth[!is.na(truth$interaction_label),],scores)
    expected <- c(104,989,1524)[match(system,systems)]
    stopifnot(nrow(evaluated)==expected, !anyDuplicated(paste(evaluated$taxon_1,evaluated$taxon_2)))
    data.frame(analysis_set=system,bootstrap=b,method="OneNet",
      auprc=ap(evaluated$interaction_label,evaluated$onenet_score),
      omitted_taxa=sum(!variable),fit_status="complete",detail="",
      elapsed_seconds=as.numeric(difftime(Sys.time(),started,units="secs")))
  },error=function(e)data.frame(analysis_set=system,bootstrap=b,method="OneNet",
      auprc=NA_real_,omitted_taxa=NA_integer_,fit_status="failed",detail=conditionMessage(e),
      elapsed_seconds=as.numeric(difftime(Sys.time(),started,units="secs"))))
  write.csv(result,paste0(stem,"_metrics.csv"),row.names=FALSE)
  result
}
workers <- as.integer(Sys.getenv("ONENET_EXTERNAL_WORKERS","3"))
message("Starting 15 fits with ",workers," workers at ",Sys.time())
results <- do.call(rbind,parallel::mclapply(seq_len(nrow(tasks)),run,mc.cores=workers,mc.preschedule=FALSE))
write.csv(results,file.path(root,"results/onenet_abundance_bootstrap_5_auprc.csv"),row.names=FALSE)
print(results)
if(any(results$fit_status!="complete")) stop("Some fits failed; inspect individual logs before updating figure")
