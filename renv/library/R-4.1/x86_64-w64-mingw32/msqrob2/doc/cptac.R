## ---- warning=FALSE, message=FALSE--------------------------------------------
library(tidyverse)
library(limma)
library(QFeatures)
library(msqrob2)
library(plotly)
library(gridExtra)

peptidesFile <- msdata::quant(pattern = "cptac_a_b_peptides", full.names = TRUE)

ecols <- grep("Intensity\\.", names(read.delim(peptidesFile)))

pe <- readQFeatures(
    table = peptidesFile, fnames = 1, ecol = ecols,
    name = "peptideRaw", sep = "\t"
)

## -----------------------------------------------------------------------------
cond <- which(strsplit(colnames(pe)[[1]][1], split = "")[[1]] == "A") # find where condition is stored
colData(pe)$condition <- substr(colnames(pe), cond, cond) %>%
    unlist() %>%
    as.factor()

## -----------------------------------------------------------------------------
rowData(pe[["peptideRaw"]])$nNonZero <- rowSums(assay(pe[["peptideRaw"]]) > 0)

## -----------------------------------------------------------------------------
pe <- zeroIsNA(pe, "peptideRaw") # convert 0 to NA

## ---- cache=TRUE--------------------------------------------------------------
MSnbase::plotNA(assay(pe[["peptideRaw"]])) +
    xlab("Peptide index (ordered by data completeness)")

## -----------------------------------------------------------------------------
pe <- logTransform(pe, base = 2, i = "peptideRaw", name = "peptideLog")
limma::plotDensities(assay(pe[["peptideLog"]]))

## -----------------------------------------------------------------------------
pe <- filterFeatures(pe, ~ Proteins %in% smallestUniqueGroups(rowData(pe[["peptideLog"]])$Proteins))

## -----------------------------------------------------------------------------
pe <- filterFeatures(pe, ~ Reverse != "+")
pe <- filterFeatures(pe, ~ Potential.contaminant != "+")

## -----------------------------------------------------------------------------
pe <- filterFeatures(pe, ~ nNonZero >= 2)
nrow(pe[["peptideLog"]])

## -----------------------------------------------------------------------------
pe <- normalize(pe,
    i = "peptideLog",
    name = "peptideNorm",
    method = "center.median"
)

## -----------------------------------------------------------------------------
limma::plotDensities(assay(pe[["peptideNorm"]]))

## -----------------------------------------------------------------------------
boxplot(assay(pe[["peptideNorm"]]),
    col = palette()[-1],
    main = "Peptide distribtutions after normalisation", ylab = "intensity"
)

## -----------------------------------------------------------------------------
limma::plotMDS(assay(pe[["peptideNorm"]]), col = as.numeric(colData(pe)$condition))

## ----warning=FALSE------------------------------------------------------------
pe <- aggregateFeatures(pe,
    i = "peptideNorm", fcol = "Proteins", na.rm = TRUE,
    name = "protein"
)

## -----------------------------------------------------------------------------
plotMDS(assay(pe[["protein"]]), col = as.numeric(colData(pe)$condition))

## ---- warning=FALSE-----------------------------------------------------------
pe <- msqrob(object = pe, i = "protein", formula = ~condition)

## -----------------------------------------------------------------------------
getCoef(rowData(pe[["protein"]])$msqrobModels[[1]])

## -----------------------------------------------------------------------------
L <- makeContrast("conditionB=0", parameterNames = c("conditionB"))
pe <- hypothesisTest(object = pe, i = "protein", contrast = L)

## ----warning=FALSE------------------------------------------------------------
volcano <- ggplot(
    rowData(pe[["protein"]])$conditionB,
    aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)
) +
    geom_point(cex = 2.5) +
    scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
    theme_minimal() +
    ggtitle("Default workflow")
volcano

## -----------------------------------------------------------------------------
sigNames <- rowData(pe[["protein"]])$conditionB %>%
    rownames_to_column("protein") %>%
    filter(adjPval < 0.05) %>%
    pull(protein)
heatmap(assay(pe[["protein"]])[sigNames, ])

## ---- warning=FALSE, message=FALSE--------------------------------------------
topN <- 5
if (length(sigNames) > topN) {
    for (protName in sigNames[1:topN])
    {
        pePlot <- pe[protName, , c("peptideNorm", "protein")]
        pePlotDf <- data.frame(longFormat(pePlot))
        pePlotDf$assay <- factor(pePlotDf$assay,
            levels = c("peptideNorm", "protein")
        )
        pePlotDf$condition <- as.factor(colData(pePlot)[pePlotDf$colname, "condition"])

        # plotting
        p1 <- ggplot(
            data = pePlotDf,
            aes(x = colname, y = value, group = rowname)
        ) +
            geom_line() +
            geom_point() +
            theme_minimal() +
            facet_grid(~assay) +
            ggtitle(protName)
        print(p1)

        # plotting 2
        p2 <- ggplot(pePlotDf, aes(x = colname, y = value, fill = condition)) +
            geom_boxplot(outlier.shape = NA) +
            geom_point(
                position = position_jitter(width = .1),
                aes(shape = rowname)
            ) +
            scale_shape_manual(values = 1:nrow(pePlotDf)) +
            labs(title = protName, x = "sample", y = "peptide intensity (log2)") +
            theme_minimal()
        facet_grid(~assay)
        print(p2)
    }
}

## ----warning=FALSE------------------------------------------------------------
pe <- aggregateFeatures(pe,
    i = "peptideNorm", fcol = "Proteins", na.rm = TRUE,
    name = "proteinMedian", fun = matrixStats::colMedians
)
pe <- msqrob(object = pe, i = "proteinMedian", formula = ~condition)
pe <- hypothesisTest(object = pe, i = "proteinMedian", contrast = L)

## -----------------------------------------------------------------------------
limma::plotMDS(assay(pe[["proteinMedian"]]),
    col = as.numeric(colData(pe)$condition)
)

## ----warning=FALSE------------------------------------------------------------
volcanoMed <- rowData(pe[["proteinMedian"]])[[colnames(L)]] %>%
    ggplot(aes(x = logFC, y = -log10(pval), color = adjPval < 0.01)) +
    geom_point(cex = 2.5) +
    scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
    theme_minimal() +
    geom_vline(xintercept = log2(0.74 / .25), col = "red") +
    ggtitle("median summarisation")

## -----------------------------------------------------------------------------
try(pe <- msqrob(
    object = pe, i = "protein", formula = ~condition,
    modelColumnName = "ridge", ridge = TRUE
)) # note: intentional error

## ----warning=FALSE, message=FALSE---------------------------------------------
pe <- msqrob(
    object = pe, i = "protein", formula = ~ -1 + condition,
    modelColumnName = "ridge", ridge = TRUE
)
Lridge <- makeContrast(
    "ridgeconditionB - ridgeconditionA = 0",
    c("ridgeconditionB", "ridgeconditionA")
)
pe <- hypothesisTest(
    object = pe, i = "protein", contrast = Lridge,
    modelColumn = "ridge"
)

## ----warning=FALSE------------------------------------------------------------
volcanoRidge <- rowData(pe[["protein"]])[[colnames(Lridge)]] %>%
    ggplot(aes(
        x = logFC, y = -log10(pval),
        color = adjPval < 0.01
    )) +
    geom_point(cex = 2.5) +
    scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
    theme_minimal() +
    geom_vline(xintercept = log2(0.74 / 0.25), col = "red") +
    ggtitle(paste("robust ridge"))

## -----------------------------------------------------------------------------
rowData(pe[["protein"]])$ups <- grepl("UPS", rownames(pe[["protein"]]))
rowData(pe[["proteinMedian"]])$ups <- grepl("UPS", rownames(pe[["proteinMedian"]]))

## -----------------------------------------------------------------------------
grid.arrange(grobs = list(volcano, volcanoMed, volcanoRidge), ncol = 1)

## -----------------------------------------------------------------------------
logFC <- data.frame(
    default = rowData(pe[["protein"]])[[colnames(L)]][, 1],
    median = rowData(pe[["proteinMedian"]])[[colnames(L)]][, 1],
    ridge = rowData(pe[["protein"]])[[colnames(Lridge)]][, 1],
    ups = rowData(pe[["protein"]])$ups
)

logFC <- logFC %>% gather(method, log2FC, c("default", "median", "ridge"))
logFC$ups <- as.factor(logFC$ups)
logFC %>% ggplot(aes(x = method, y = log2FC, fill = ups)) +
    geom_boxplot() +
    geom_hline(yintercept = log2(0.74 / .25), color = "red")

## -----------------------------------------------------------------------------
tprFdp <- function(pval, tp, adjPval) {
    ord <- order(pval)
    return(data.frame(
        pval = pval[ord],
        adjPval = adjPval[ord],
        tpr = cumsum(tp[ord]) / sum(tp),
        fdp = cumsum(!tp[ord]) / 1:length(tp)
    ))
}

## -----------------------------------------------------------------------------
tprFdpDefault <- tprFdp(
    rowData(pe[["protein"]])[[colnames(L)]]$pval,
    rowData(pe[["protein"]])$ups,
    rowData(pe[["protein"]])[[colnames(L)]]$adjPval
)
tprFdpMedian <- tprFdp(
    rowData(pe[["proteinMedian"]])[[colnames(L)]]$pval,
    rowData(pe[["proteinMedian"]])$ups,
    rowData(pe[["proteinMedian"]])[[colnames(L)]]$adjPval
)

tprFdpRidge <- tprFdp(
    rowData(pe[["protein"]])[[colnames(Lridge)]]$pval,
    rowData(pe[["protein"]])$ups,
    rowData(pe[["protein"]])[[colnames(Lridge)]]$adjPval
)

hlp <- rbind(
    cbind(tprFdpDefault, method = "default"),
    cbind(tprFdpMedian, method = "median"),
    cbind(tprFdpRidge, method = "ridge")
)
tprFdpPlot <- hlp %>%
    ggplot(aes(x = fdp, y = tpr, color = method)) +
    geom_path()
tprFdpPlot

