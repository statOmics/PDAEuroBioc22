library(tidyverse)
library(limma)
library(QFeatures)
library(msqrob2)


## -----------------------------------------------------------------------------
(basename(f <- msdata::quant(full.names = TRUE)))

grep("Intensity\\.", names(read.delim(f)), value = TRUE)
(ecols <- grep("Intensity\\.", names(read.delim(f))))

qf <- readQFeatures(
    f, fnames = 1, ecol = ecols,
    name = "peptideRaw", sep = "\t")
qf


## -----------------------------------------------------------------------------
qf[[1]]
qf[["peptideRaw"]]

assay(qf[[1]])[1:10, 1:3]

rowData(qf[["peptideRaw"]])[, c("Proteins", "Sequence", "Charges")]

colData(qf)

colnames(qf)[[1]]
(new_names <- sub("Intensity\\.", "", colnames(qf)[[1]]))

qf <- renameColname(qf, i = 1, new_names) |>
  renamePrimary(new_names)

qf$lab <- rep("lab3", 6)
qf$condition <- factor(rep(c("A", "B"), each = 3))
qf$spikeConcentration <- rep(c(A = 0.25, B = 0.74),
                             each = 3)

colData(qf)


## ----nNA----------------------------------------------------------------------
qf <- zeroIsNA(qf, "peptideRaw")
na <- nNA(qf[[1]])
na

table(na$nNArows$nNA)

rowData(qf[[1]])$keepNA <- na$nNArows$nNA <= 4


## ----logTransg----------------------------------------------------------------
qf <- logTransform(qf, base = 2,
                   i = "peptideRaw",
                   name = "peptideLog")
qf



## ----filtering----------------------------------------------------------------
sug <- smallestUniqueGroups(rowData(qf[["peptideRaw"]])$Proteins)
filterFeatures(qf, ~ Proteins %in% sug)

filterFeatures(qf, ~ Reverse != "+")
filterFeatures(qf, ~ Potential.contaminant != "+")

filterFeatures(qf, ~ keepNA)

qf <- qf |>
    filterFeatures(~ Proteins %in% sug) |>
    filterFeatures(~ Reverse != "+") |>
    filterFeatures(~ Potential.contaminant != "+") |>
    filterFeatures(~ keepNA)
qf


## ----normalise----------------------------------------------------------------
qf <- normalize(qf,
                i = "peptideLog",
                name = "peptideNorm",
                method = "center.median")
qf


## ----densityplot--------------------------------------------------------------
as_tibble(longFormat(qf[, , 2:3], colvars = "condition")) %>%
    ggplot(aes(x = value, group = primary, colour = condition)) +
    geom_density() +
    facet_grid(assay ~ .)


## ----mdsplot------------------------------------------------------------------
assay(qf[["peptideNorm"]]) |>
    limma::plotMDS(col = as.numeric(qf$condition))


## -----------------------------------------------------------------------------
qf <- aggregateFeatures(qf,
  i = "peptideNorm",
  fcol = "Proteins",
  na.rm = TRUE,
  name = "proteinMedian",
  fun = matrixStats::colMedians)
qf


## -----------------------------------------------------------------------------
assay(qf[["proteinMedian"]]) %>%
  limma::plotMDS(col = as.numeric(qf$condition))


## ----msqrob, warning=FALSE----------------------------------------------------
qf <- msqrob(object = qf,
             i = "proteinMedian",
             formula = ~condition,
             overwrite = TRUE)


## -----------------------------------------------------------------------------
rowData(qf[["proteinMedian"]])[, c("Proteins", ".n", "msqrobModels")]

## -----------------------------------------------------------------------------
getCoef(rowData(qf[["proteinMedian"]])$msqrobModels[[1]])

## -----------------------------------------------------------------------------
L <- makeContrast("conditionB=0", parameterNames = c("conditionB"))
qf <- hypothesisTest(object = qf, i = "proteinMedian", contrast = L)


## ----warning=FALSE------------------------------------------------------------
tmp <- rowData(qf[["proteinMedian"]])$conditionB[complete.cases(rowData(qf[["proteinMedian"]])$conditionB),]
tmp$shapes <- 16

volcanoMedian<- ggplot(tmp,
                  aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)) +
  geom_point(cex = 2.5, shape = tmp$shapes) +
  scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
  theme_bw() +
    ggtitle(paste0("Median: TP = ",sum(tmp$adjPval<0.05&grepl(rownames(tmp),pattern ="UPS"),na.rm=TRUE),
                   " FP = ", sum(tmp$adjPval<0.05&!grepl(rownames(tmp),pattern ="UPS"),na.rm=TRUE)))
volcanoMedian


## -----------------------------------------------------------------------------
sigNames <- rowData(qf[["proteinMedian"]])$conditionB %>%
  rownames_to_column("proteinMedian") %>%
  filter(adjPval < 0.05) %>%
  pull(proteinMedian)

heatmap(assay(qf[["proteinMedian"]][sigNames, ]),cexRow = 1, cexCol = 1)

sigProteins <- rowData(qf[["proteinMedian"]])$conditionB %>%
  rownames_to_column("proteinMedian") %>%
   filter(grepl("UPS",proteinMedian)) %>%
  pull(proteinMedian)

heatmap(assay(qf[["proteinMedian"]])[sigProteins, ], cexCol = 1)

## -----------------------------------------------------------------------------
rowData(qf[["proteinMedian"]])$conditionB %>%
  rownames_to_column(var = "protein") %>%
  mutate(ups = grepl("UPS",protein)) %>%
  ggplot(aes(x = ups, y = logFC, fill = ups)) +
  geom_boxplot() +
  theme_bw() +
  geom_hline(yintercept = log2(0.74 / .25), color = "#00BFC4") +
    geom_hline(yintercept = 0, color = "#F8766D")



## ----ex, warning=FALSE--------------------------------------------------------
## aggregation
qf <- aggregateFeatures(qf,
  i = "peptideNorm",
  fcol = "Proteins",
  na.rm = TRUE,
  name = "proteinRobust",
  fun = MsCoreUtils::robustSummary)

## estimation
qf <- msqrob(object = qf,
             i = "proteinRobust",
             formula = ~ condition,
             overwrite = TRUE)

## inference
L <- makeContrast("conditionB=0", parameterNames = c("conditionB"))
qf <- hypothesisTest(object = qf, i = "proteinRobust", contrast = L)

## volcano plot
tmp <- rowData(qf[["proteinRobust"]])$conditionB[complete.cases(rowData(qf[["proteinRobust"]])$conditionB),]
tmp$shapes <- 16

volcanoRobust <- ggplot(tmp,
                  aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)) +
  geom_point(cex = 2.5, shape = tmp$shapes) +
  scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
  ggtitle(paste0("Median: TP = ",
                 sum(tmp$adjPval<0.05 & grepl(rownames(tmp), pattern ="UPS"), na.rm=TRUE),
                 " FP = ",
                 sum(tmp$adjPval<0.05 & !grepl(rownames(tmp), pattern ="UPS"), na.rm=TRUE)))
volcanoRobust
