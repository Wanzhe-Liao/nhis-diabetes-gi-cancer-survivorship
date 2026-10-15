source('scripts/rev1_runtime.R')
source('scripts/rev1_domain_design.R')
suppressPackageStartupMessages(library(survey))
options(survey.lonely.psu='fail')
# A PSU repeated across years remains one variance unit, while redesigns separate it.
x <- data.frame(year=c(1997,2005,2006,2015,2016,2018),design_strata=1L,design_psu=1L)
k <- rev1_period_identifiers(x)$design_psu_prefixed
stopifnot(k[1]==k[2],k[3]==k[4],k[5]==k[6],length(unique(k))==3L)
# Two-PSU domain example. One PSU contains no domain rows, but must be sampled.
f <- data.frame(publicid=letters[1:4],year=c(2006,2007,2006,2007),
                design_strata=1L,design_psu=c(1L,1L,2L,2L),sa_wgt_pool=1,
                gi_any=c(TRUE,TRUE,FALSE,FALSE),y=c(2,4,0,0))
f <- rev1_period_identifiers(f)
b <- f[f$gi_any,];sm <- rev1_prepare_domain_sampler(b,f)
stopifnot(length(sm)==1L,length(sm[[1]]$psu_indices)==2L,
          length(sm[[1]]$psu_indices[[2]])==0L)
set.seed(731)
z <- replicate(5000,{r<-rev1_resample_domain(b,sm);sum(r$y*r$sa_wgt_pool)})
stopifnot(all(z %in% c(0,12)),any(z==0),any(z==12))
d <- rev1_domain_population_design(f)
exact <- as.numeric(vcov(svytotal(~y,d)))
stopifnot(abs(var(z)/exact-1)<0.02)
# Compare to the maintained survey package's formal Rao-Wu replicate weights.
set.seed(731)
ref <- survey::subbootweights(f$design_strata_prefixed,f$design_psu_prefixed,
                              replicates=5000,compress=FALSE)
ref_total <- colSums(ref$repweights * (f$y*f$sa_wgt_pool))
stopifnot(abs(var(z)/var(ref_total)-1)<0.02)
# A multi-PSU stratum must use multiplicity as weights without duplicate rows.
g <- data.frame(publicid=letters[1:6],year=2016L,design_strata=100L,
                design_psu=rep(1:3,each=2),sa_wgt_pool=1,gi_any=TRUE,y=1:6)
g <- rev1_period_identifiers(g);sg <- rev1_prepare_domain_sampler(g,g)
set.seed(37)
for (i in 1:100) {
 r<-rev1_resample_domain(g,sg)
 stopifnot(!anyDuplicated(r$publicid),all(r$sa_wgt_pool %in% c(1.5,3)),
           sum(r$sa_wgt_pool)==6)
}
full <- rev1_period_identifiers(read_rev1_parquet(rev1_fullsample_path()))
stopifnot(nrow(full)==646201L,length(unique(full$design_strata_prefixed))==691L,
          length(unique(full$design_psu_prefixed))==1937L)
cat('PASS: design periods; full-frame domain PSU sampling; exact two-PSU variance; survey::subbootweights comparison; unique-row replicate weights; full population dimensions.\n')
# Independently reconstruct the fitted model covariance from weighted influences.
# No survey variance routine is used here; absent-domain PSUs contribute zeros.
fit <- readRDS(file.path(rev1_output_dir(),'rev1_fits.rds'))$fitA
cox_fit <- fit; class(cox_fit) <- 'coxph'
influence <- residuals(cox_fit,type='dfbeta',weighted=TRUE)
domain <- fit$survey.design
if (!is.null(fit$na.action)) domain <- domain[-as.integer(fit$na.action),]
stopifnot(nrow(influence)==nrow(domain$variables))
units <- unique(full[c('design_strata_prefixed','design_psu_prefixed')])
psu_sum <- rowsum(influence,domain$variables$design_psu_prefixed,reorder=FALSE)
totals <- matrix(0,nrow(units),ncol(influence))
index <- match(rownames(psu_sum),units$design_psu_prefixed)
stopifnot(!anyNA(index)); totals[index,] <- psu_sum
manual <- matrix(0,ncol(totals),ncol(totals))
for (index in split(seq_len(nrow(units)),units$design_strata_prefixed)) {
 n <- length(index)
 centered <- sweep(totals[index,,drop=FALSE],2,colMeans(totals[index,,drop=FALSE]),'-')
 manual <- manual+n/(n-1)*crossprod(centered)
}
stopifnot(max(abs(manual-vcov(fit)))<1e-10)
cat('PASS: independent PSU influence aggregation matches fitted Taylor covariance.\n')
