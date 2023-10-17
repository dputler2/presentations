# This script performs a Monte Carlo experiment to look at the effect of
# how increasing sample sizes interplay with correlated predictors.

# Base path for the plots (please provide)
basePath_sc <- ""

# A helper function to determine the count of incorrect coefficient signs
wrong_sign <- function(x) {
  a0_si <- sum(as.numeric(x[,1] < 0))
  a1_si <- sum(as.numeric(x[,2] > 0))
  a2_si <- sum(as.numeric(x[,3] < 0))
  a3_si <- sum(as.numeric(x[,4] < 0))
  c(a0_si, a1_si, a2_si, a3_si)
}

# The functions that runs an experiment for a given sample size
experiment <- function(sampSize, b, cnames, corMat, runs = 100) {
  outMat_mn <- matrix(NA, nrow = runs, ncol = 4)
  for(i in 1:runs) {
    draws_df <- as.data.frame(rmvnorm(sampSize, sigma = corMat))
    names(draws_df) <- c("x1", "x2", "x3")
    draws_df$y <- b[1] + b[2]*draws_df$x1 + b[3]*draws_df$x2 + b[4]*draws_df$x3 + rnorm(sampSize)
    outMat_mn[i,] <- coefficients(lm(y ~ x1 + x2 + x3, data = draws_df))
  }
  mean_vn <- apply(outMat_mn, 2, mean)
  sd_vn <- apply(outMat_mn, 2, sd)
  ws_vn <- wrong_sign(outMat_mn)
  meanSD_df <- data.frame(Mean = mean_vn, SD = sd_vn, Wrong_Sign = ws_vn)
  row.names(meanSD_df) <- cnames
  betas_df <- as.data.frame(outMat_mn)
  names(betas_df) <- cnames
  list(Mean_SD = meanSD_df, Betas = betas_df)
}

# The function for creating the customized histograms
custom_hist <- function(x,
                        breaks = 16,
                        xlim = range(x),
                        main = "Histogram",
                        xlab = "Variable",
                        truth = -1,
                        split = 0,
                        good = "first",
                        lcex = 1.0,
                        pcex = NULL) {
    colors_vc <- c("#1F77B4", "#FF7F0E")
    if (good == "first") {
        c1_sc <- colors_vc[1]
        c2_sc <- colors_vc[2]
    } else {
        c1_sc <- colors_vc[2]
        c2_sc <- colors_vc[1]
    }
    hist1_l <- hist(x, breaks = breaks, plot = FALSE)
    totC1_si <- length(hist1_l$mids[hist1_l$mids <= split])
    totC2_si <- length(hist1_l$mids) - totC1_si
    print(c(totC1_si, totC2_si, length(hist1_l$mids)))
    colors_vc <- c(rep(c1_sc, totC1_si), rep(c2_sc, totC2_si))
    if (is.null(pcex)) {
        hist(x,
             breaks = breaks,
             col = colors_vc,
             xlim = xlim,
             main = main,
             xlab = xlab)
    } else {
        hist(x,
             breaks = breaks,
             col = colors_vc,
             xlim = xlim,
             main = main,
             xlab = xlab,
             cex.lab = pcex,
             cex.main = pcex)

    }
    lines(c(truth, truth), c(0, max(hist1_l$counts)), col = "#F0E442", lwd = 4)
    legend("topright",
           c("Correct Sign", "Incorrect Sign", "True Value"),
           fill = c("#1F77B4", "#FF7F0E", "#F0E442"),
           cex = lcex)
    invisible()
}

# Read in the mvtnorm library
library(mvtnorm)

# The covariance/correlation matrix between the variables
# Note: The standard deviation of each variable is one, so the correlation
# and covariance matrix are one in the same
corMatHigh_mn <- rbind(c(1.00, 0.99, 0.30),
                          c(0.99, 1.00, 0.30),
                          c(0.30, 0.30, 1.00))

corMatMod_mn <- rbind(c(1.00, 0.50, 0.15),
                         c(0.50, 1.00, 0.15),
                         c(0.15, 0.15, 1.00))

# True coefficients
b0 <- 0.5
b1 <- -1.0
b2 <- 1.0
b3 <- 1.5
coefs_vn <- c(b0, b1, b2, b3)

# Names of the coefficients/variables
coefNames_vc <- c("Intercept", "x1", "x2", "x3")

# The non-collinear case
set.seed(123)
sampSz50Ident_l <- experiment(50, coefs_vn, coefNames_vc, corMat = diag(rep(1, 3)), runs = 1000)

# The moderately collinear case
set.seed(123)
sampSz50Mod_l <- experiment(50, coefs_vn, coefNames_vc, corMat = corMatMod_mn, runs = 1000)

# The highly collinear case
set.seed(123)
sampSz50_l <- experiment(50, coefs_vn, coefNames_vc, corMat = corMatHigh_mn, runs = 1000)
sampSz100_l <- experiment(100, coefs_vn, coefNames_vc, corMat = corMatHigh_mn, runs = 1000)
sampSz150_l <- experiment(150, coefs_vn, coefNames_vc, corMat = corMatHigh_mn, runs = 1000)
sampSz200_l <- experiment(200, coefs_vn, coefNames_vc, corMat = corMatHigh_mn, runs = 1000)
sampSz2500_l <- experiment(2500, coefs_vn, coefNames_vc, corMat = corMatHigh_mn, runs = 1000)
smlLims_vn <- range(sampSz50_l$Betas$x1)

# Create the analysis histograms
png(filename = paste0(basePath_sc, "/plot1.png"),
    width = 600,
    height = 600)
layoutMat_mn <- rbind(c(1, 1, 2, 2),
                      c(1, 1, 2, 2),
                      c(0, 3, 3, 0),
                      c(0, 3, 3, 0))
layout(layoutMat_mn)
custom_hist(sampSz50Ident_l$Betas$x1, xlim = smlLims_vn, breaks = 16, xlab = "Estimate of b1", main = "No Collinearity Case", lcex = 1.2, pcex = 1.4)
custom_hist(sampSz50Mod_l$Betas$x1, xlim = smlLims_vn, breaks = 16, xlab = "Estimate of b1", main = "Moderate Collinearity Case", lcex = 1.2, pcex = 1.4)
custom_hist(sampSz50_l$Betas$x1, xlim = smlLims_vn, breaks = 16, xlab = "Estimate of b1", main = "High Collinearity Case", lcex = 1.2, pcex = 1.4)
dev.off()

png(filename = paste0(basePath_sc, "/plot2.png"),
    width = 600,
    height = 600)
layout(matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE))
custom_hist(sampSz50_l$Betas$x1, xlim = smlLims_vn, breaks = 16, xlab = "Estimate of b1", main = "Sample Size 50", lcex = 0.95, pcex = 1.2)
custom_hist(sampSz100_l$Betas$x1, xlim = smlLims_vn, breaks = 16, xlab = "Estimate of b1", main = "Sample Size 100", lcex = 0.95, pcex = 1.2)
custom_hist(sampSz150_l$Betas$x1, xlim = smlLims_vn, breaks = 16, xlab = "Estimate of b1", main = "Sample Size 150", lcex = 0.95, pcex = 1.2)
custom_hist(sampSz200_l$Betas$x1, xlim = smlLims_vn, breaks = 16, xlab = "Estimate of b1", main = "Sample Size 200", lcex = 0.95, pcex = 1.2)
dev.off()

png(filename = paste0(basePath_sc, "/plot3.png"),
    width = 600,
    height = 300)
layout(matrix(1:2, nrow = 1, ncol = 2, byrow = TRUE))
custom_hist(sampSz50Ident_l$Betas$x1, xlim = smlLims_vn, breaks = 16, xlab = "Estimate of b1", main = "No Collinearity, Sample Size 50", lcex = 0.8, pcex = 0.9)
custom_hist(sampSz2500_l$Betas$x1, xlim = smlLims_vn, breaks = 16, xlab = "Estimate of b1", main = "High Collinearity, Sample Size 2500", lcex = 0.8, pcex = 0.9)
dev.off()

# Statistical summaries of the coefficients estimates
cat("\nSample size of 50, no collinearity\n")
print(sampSz50Ident_l$Mean_SD)
cat("\nSample size of 50, moderate collinearity\n")
print(sampSz50Mod_l$Mean_SD)
cat("\nSample size of 50\n")
print(sampSz50_l$Mean_SD)
cat("\n\nSample size of 100\n")
print(sampSz100_l$Mean_SD)
cat("\n\nSample size of 150\n")
print(sampSz150_l$Mean_SD)
cat("\n\nSample size of 200\n")
print(sampSz200_l$Mean_SD)
cat("\n\nSample size of 2500\n")
print(sampSz2500_l$Mean_SD)
