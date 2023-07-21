# Load packages
library(readxl)
library(tidyverse)
library(magrittr)
library(meta)
library(dmetar)

# Load and clean data
pjt_data <- read_excel("/Users/changus/Documents/R/matrix_dec.xlsx") |>
  janitor::clean_names() |>
  dplyr::rename("pjxsi_le" = "p_jx_si_le")
pjt_data_full <- pjt_data[
  !is.na(pjt_data$pj_si) & !is.na(pjt_data$si_le) & !is.na(pjt_data$pj_le),]
pjt_data_h5 <- pjt_data[!is.na(pjt_data$pjxsi_le),]

# Set data as numeric
pjt_data$demin <- as.numeric(pjt_data$demin)
pjt_data$hdi <- as.numeric(pjt_data$hdi)
pjt_data$male_percent <- as.numeric(pjt_data$male_percent)
pjt_data$minority_percent <- as.numeric(pjt_data$minority_percent)
pjt_data$mean_age <- as.numeric(pjt_data$mean_age)

# test correlation between PJ and Legit
h1 <- metacor(
  cor = pj_le,
  n = n,
  studlab = author,
  data = pjt_data,
  prediction = TRUE,
  fixed = FALSE,
  random = TRUE,
  hakn = TRUE,
  title = "h1"
)

# test correlation between SI and Legit

h2 <- metacor(
  cor = si_le,
  n = n,
  studlab = author,
  data = pjt_data,
  prediction = TRUE,
  fixed = FALSE,
  random = TRUE,
  hakn = TRUE,
  title = "h2"
)


# test correlation between PJxSI and Legit

h5 <- metacor(
  cor = pjxsi_le,
  n = n,
  studlab = author,
  data = pjt_data,
  fixed = FALSE,
  random = TRUE,
  hakn = TRUE,
  title = "h5"
)

h5_full <- metacor(
  cor = pjxsi_le,
  n = n,
  studlab = author,
  data = pjt_data_h5,
  fixed = FALSE,
  random = TRUE,
  hakn = TRUE,
  title = "h5_full data"
)

# influential analysis of H1,2,5
# h1
h1_inf <- InfluenceAnalysis(h1, random = TRUE)
h1_removed <- update.meta(h1, exclude = c(38, 40)) |>
  summary()
# h2
h2_inf <- InfluenceAnalysis(h2, random = TRUE)
h2_removed <- update.meta(h2, exclude = c(12)) |>
  summary()
# h5
h5_inf <- InfluenceAnalysis(h5, random = TRUE)
h5_removed <- update.meta(h5, exclude = c(12)) |>
  summary()

# trim and fill----
# h1
h1_tf <- trimfill(h1)
h1_removed_tf <- trimfill(update(h1, 
                                 subset = -c(38, 40)))
#h2
h2_tf <- trimfill(h2)
h2_removed_tf <- trimfill(update(h2, 
                                 subset = -c(12)))
#h5
h5_tf <- trimfill(h5)
h5_removed_tf <- trimfill(update(h5, 
                                 subset = -c(12)))

# metareg----
h1_age <- metareg(h1, ~mean_age)
h1_gender <- metareg(h1, ~male_percent)
h1_race <- metareg(h1, ~minority_percent)
h1_dem <- metareg(h1, ~demin)
h1_hdi <- metareg(h1, ~hdi)

h2_age <- metareg(h2, ~mean_age)
h2_gender <- metareg(h2, ~male_percent)
h2_race <- metareg(h2, ~minority_percent)
h2_dem <- metareg(h2, ~demin)
h2_hdi <- metareg(h2, ~hdi)

h5_age <- metareg(h5, ~mean_age)
h5_gender <- metareg(h5, ~male_percent)
h5_race <- metareg(h5, ~minority_percent)
h5_dem <- metareg(h5, ~demin)
h5_hdi <- metareg(h5, ~hdi)

# subgroup analysis-------------------------------------------------------------
# NA vs Euro vs rest of the world
#h1_location <- update.meta(h1, subgroup = region, tau.common = FALSE)
#h2_location <- update.meta(h2, subgroup = region, tau.common = FALSE)
#h5_location <- update.meta(h5, subgroup = region, tau.common = FALSE)

# study design
pjt_data$design_recoded <- recode(pjt_data$design, "L" = "C", "R" = "E")
h1_design <- update.meta(h1, subgroup = pjt_data$design_recoded, tau.common = FALSE)

# si measure
# create data with only si measure
pjt_data_si <- pjt_data |> 
  filter(!is.na(si_measure)) |> 
  filter(!is.na(si_le))

# run h2 wtih only si measure
h2_si <- metacor(
  cor = si_le,
  n = n,
  studlab = author,
  data = pjt_data_si,
  fixed = FALSE,
  random = TRUE,
  hakn = TRUE,
  title = "h2_si"
)

# sub group analysis for si measure
h2_si_sub_group <- update.meta(h2_si, subgroup = si_measure, tau.common = FALSE)


# h3, h4 & h6-------------------------------------------------------------------
# load packages for metaSEM
library("metaSEM")
library("lavaan")

# load data as matrix
pjt_data_matrix <- as.matrix(pjt_data[, 5:10])
pjt_data_matrix_full <- as.matrix(pjt_data_full[, 5:10]) 


# reshape into 3 x 3
pjt_data_matrix <- lapply(
  split(pjt_data_matrix, 1:160),
  function(x) {
    mat <- matrix(1, ncol = 3, nrow = 3);
    mat[upper.tri(mat, diag = TRUE)] <- x;
    mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)];
    mat
  }
)

pjt_data_matrix_full <- lapply(
  split(pjt_data_matrix_full, 1:19),
  function(x) {
    mat <- matrix(1, ncol = 3, nrow = 3);
    mat[upper.tri(mat, diag = TRUE)] <- x;
    mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)];
    mat
  }
)

# name and match sample size
names(pjt_data_matrix) <- pjt_data$author
pjt <- list(data = pjt_data_matrix, n = pjt_data$n)

names(pjt_data_matrix_full) <- pjt_data_full$author
pjt_full <- list(data = pjt_data_matrix_full, n = pjt_data_full$n)

# rem analysis
rem <- tssem1(
  pjt$data, pjt$n,
  method = "REM",
  RE.type = "Diag")

rem_full <- tssem1(
  pjt_full$data, pjt_full$n,
  method = "REM",
  RE.type = "Diag")

# h34
# # matrix A and S method (error cannot be solved)
A_h34 <- matrix(c(0          , 0           , 0            ,
                 0           , 0           , 0            ,
                 "0.3*PJ_Leigt", "0.3*SI_Legit", 0            ),
               ncol = 3, nrow=3, byrow=TRUE)
dimnames(A_h34)[[1]] <- dimnames(A_h34)[[2]] <- c("PJ", "SI", "Legit")
A_h34 <- as.mxMatrix(A_h34)

S_h34 <- matrix(c(1        , "0.2*PJ_SI", 0                 ,
                  "0.2*PJ_SI", 1        , 0                 ,
                  0        , 0        , "0.1*ErrorVarLegit"),
                ncol = 3, nrow=3, byrow=TRUE)
dimnames(S_h34)[[1]] <- dimnames(S_h34)[[2]] <- c("PJ", "SI", "Legit")
S_h34 <- as.mxMatrix(S_h34)

h34_matrix_method <- tssem2(rem, 
                           Amatrix = A_h34, 
                           Smatrix = S_h34, 
                           intervals.type = "LB", 
                           diag.constraints = TRUE)

multivariate_model <- "## Legit is modeled by PJ and SI
                       Legit ~ PJ2Legit*PJ + SI2Legit*SI
                       ## Variances of predictors are fixed at 1
                       PJ ~~ 1*PJ
                       SI ~~ 1*SI
                       ## Correlation between the predictors
                       PJ ~~ PJCorSI*SI
                       ## Error variance
                       Legit ~~ ErrorVarLegit*Legit"

ram_h34 <- lavaan2RAM(
  multivariate_model,
  obs.variables = c("PJ", "SI", "Legit"))

h34 <- tssem2(
  rem, RAM = ram_h34,
  intervals.type = "LB",
  model.name = "multivariate_model")

# h6
# matrix A and S method (error cannot be solved)
A_h6 <- matrix(c(0            , 0          , 0            ,
              "0.3*PJ_SI", 0            , 0            ,
              "0.3*PJ_Leigt", "0.3*SI_Legit", 0            ),
              ncol = 3, nrow=3, byrow=TRUE)
dimnames(A_h6)[[1]] <- dimnames(A_h6)[[2]] <- c("PJ", "SI", "Legit")
A_h6 <- as.mxMatrix(A_h6)

S_h6 <- Diag(c(1, "0.2*ErrVarSI", "0.2*ErrVarLegit"))
dimnames(S_h6)[[1]] <- dimnames(S_h6)[[2]] <- c("PJ", "SI", "Legit")
S_h6 <- as.mxMatrix(S_h6)

h6_matrix_method <- tssem2(rem_full, 
             Amatrix = A_h6, 
             Smatrix = S_h6, 
             intervals.type = "LB", 
             diag.constraints = TRUE,
             mx.algebras = list(
               indirectEffect = mxAlgebra(PJ_SI*SI_Legit,
                                          name="indirectEffect")))

# ram method
mediation_model <- "## Legit is modeled by PJ and SI
                       Legit ~ PJ2Legit*PJ + SI2Legit*SI
                       ## SI is modelled by PJ
                       SI ~ PJ2SI*PJ
                       ## Variances of predictors are fixed at 1
                       PJ ~~ 1*PJ
                       ## Error variance
                       Legit ~~ ErrorVarLegit*Legit
                       SI ~~ ErrorVarSI*SI"

ram_h6 <- lavaan2RAM(
  mediation_model,
  obs.variables = c("PJ", "SI", "Legit"))

h6 <- tssem2(
  rem_full, RAM = ram_h6,
  intervals.type = "LB",
  model.name = "GEM",
  mx.algebras = list(indirectEffect = mxAlgebra(PJ2SI*SI2Legit,
                                                name="indirectEffect")))

# plot--------------------------------------------------------------------------
library(semPlot)

h34_plot <- meta2semPlot(h34_matrix_method)
labels <- c("Procedural\njustice","Social\nidentity","Legitimacy")
semPaths(h34_plot, 
         whatLabels = "est", 
         edge.color = "black", 
         layout="tree",
         sizeMan = 10,
         edge.label.cex = 1,
         rotation=2,
         nodeLabels = labels)

h6_plot <- meta2semPlot(h6_matrix_method)
semPaths(h6_plot, 
         whatLabels = "est", 
         edge.color = "black",
         sizeMan = 10,
         edge.label.cex = 1,
         layout = "tree2", 
         rotation = 1,
         nodeLabels = labels)



h6_graph <- plot(h6, whatLabels="path", edge.label.cex=0.8)

