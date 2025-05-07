if (!require("tidyverse")) {
  install.packages("tidyverse")
}
library("tidyverse")

if (!require("reticulate")) {
  install.packages("reticulate")
}
library("reticulate")

if (!require("BiocManager")) {
  install.packages("BiocManager")
}
library("BiocManager")

if (!require("BCRANK")) {
  BiocManager::install("BCRANK", ask = FALSE, update = FALSE)
}
library("BCRANK")

if (!require("Biostrings")) {
  BiocManager::install("Biostrings", ask = FALSE, update = FALSE)
}
library("Biostrings")

if (!require("shiny")) {
  install.packages("shiny")
}
library("shiny")

if (!require("shinyjs")) {
  install.packages("shinyjs")
}
library("shinyjs")

if (!require("shinycssloaders")) {
  install.packages("shinycssloaders")
}
library("shinycssloaders")