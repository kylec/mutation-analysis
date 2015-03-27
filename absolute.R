library(ABSOLUTE)
## usage
## Rscript absolute.R [project dir] [project] [outputdir] [maffile] [segfile] [patientid] [build] [review]
## Rscript absolute.R ~/Projects fap absolute Vilar14.maf Vilar14.seg Vilar14 hg19 1
# review = 1 when you run review mode
## Args
args <- commandArgs(TRUE)

genome = args[7]
platform = "Illumina_WES"
primary.disease = args[2]
sample.name = args[6]
seg.dat.fn = args[5]
maf.fn = args[4] 
sigma.p <- 0
max.sigma.h <- 0.02
min.ploidy <- 0.95
max.ploidy <- 6
max.as.seg.count <- 16000
max.non.clonal <- 0
max.neg.genome <- 0
copy_num_type <- "total"
min.mut.af=0.1

## create absolute output folder for a patient
results.dir=file.path(args[1], args[2], args[3], sample.name)
log.dir=file.path(results.dir, "log")
## create absolute review object
# create review object
obj.name <- sample.name
review.dir <- file.path(results.dir)
absolute.files <- file.path(results.dir, paste0(sample.name,".ABSOLUTE.RData"))
calls.path = file.path(review.dir, paste0(sample.name, ".PP-calls_tab.txt"))
modes.path = file.path(review.dir, paste0(sample.name, ".PP-modes.data.RData"))
output.path = file.path(review.dir)
print("Start absolute...")
print(length(args))
# create when it's non-review mode
if (length(args)==7) {
  print("Run absolute in non-review mode")
  if (!file.exists(log.dir)) {
      dir.create(log.dir, recursive=TRUE)
  }
  if (!file.exists(results.dir)) {
      dir.create(results.dir, recursive=TRUE)
  }
  # send r output to file
  sink(file=file.path(log.dir, paste0(sample.name, ".abs.out.txt")))
  RunAbsolute(seg.dat.fn, sigma.p, max.sigma.h, min.ploidy, max.ploidy, primary.disease, platform, sample.name, results.dir, max.as.seg.count, max.non.clonal, max.neg.genome, copy_num_type, maf.fn, min.mut.af, verbose=TRUE)
  sink() 
  CreateReviewObject(obj.name, absolute.files, review.dir, "total", verbose=TRUE)
  ExtractReviewedResults(calls.path, "final", modes.path, output.path, "absolute", "total")
} else {
  print("Running review mode-re-extract absolute results...")
  ## mark up the file output/abs_summary/DRAWS_summary.PP-calls_tab.txt by prepending a column with your desired solution calls
  ExtractReviewedResults(calls.path, "final", modes.path, output.path, "absolute", "total")
}
