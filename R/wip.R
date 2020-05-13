
##########
# File: aa_properties.R
##########

# cdrProp <- function (.data, .prop = c("hydro", "polarity+turn"), .region = c("vhj", "h", "v", "j")) {
#   proplist = lapply(.prop, function(.prop) { chooseprop(.prop) })
#
#   if (!has_class(.data, "list")) {
#     .data = list(Data = .data)
#   }
#
#   data("aaproperties")
#
#   res = .data[IMMCOL$cdr3aa]
#
#   for (property in proplist){
#   new = .data %>%
#     select(IMMCOL$cdr3aa) %>%
#     mutate(hydro = aapropeval(IMMCOL$cdr3aa, property)) %>%
#     collect()
#
#   colnames(new) <- c("IMMCOL$cdr3aa", property)
#   res <- dplyr::bind_cols(res, new[property])
#   }
#
#   add_class(res, "immunr_cdr_prop")
# }


# chooseprop <- function(prop) {
#   switch(prop,
#          alpha = "alpha",
#          beta = "beta",
#          charge = "charge",
#          core = "core",
#          hydro = "hydropathy",
#          ph = "pH",
#          polar = "polarity",
#          rim = "rim",
#          surf = "surface",
#          turn = "turn",
#          vol = "volume",
#          str = "strength",
#          dis = "disorder",
#          high = "high_contact",
#          stop("Unknown property name"))
#   }
#
# aapropeval <- function(seq, col){
#   aaproperty <- AA_PROP[,c("amino.acid", col)]
#   seq <- strsplit(x = seq, split = "")
#   aaseqpropvalue <- lapply(seq, function(seq) {
#     sum(aaproperty[seq, ][[col]], na.rm = TRUE) / length(seq) })
#   return(aaseqpropvalue)
# }
#
# cdrPropAnalysis <- function (.data, .method = c("t.test")) {
#
# }


##########
# File: graph.R
##########

# mutationNetwork <- function (.data, .method = c("hamm", "lev"), .err = 2) {
#   require(igraph)
#   UseMethod("mutationNetwork")
# }
#
# mutationNetwork.character <- function (.data, .method = c("hamm", "lev"), .err = 2) {
#   add_class(res, "immunr_mutation_network")
# }
#
# mutationNetwork.immunr_shared_repertoire <- function (.data, .method = c("hamm", "lev"), .err = 2) {
#   add_class(res, "immunr_mutation_network")
# }
#
# mutationNetwork.tbl <- function (.data, .col, .method = c("hamm", "lev"), .err = 2) {
#   select_(.data, .dots = .col)
#   add_class(res, "immunr_mutation_network")
# }
#
# mut.net = mutationNetwork


##########
# File: post_analysis.R
#########

# overlap => distance matrix
# gene usage => N-dimensional vector of values
# diversity => vector of values

# (both) vector of values => distance matrix
# N-dimensional vector of values => clustering
# N-dimensional vector of values => dimensionality reduction
# N-dimensional vector of values => statistical test
# N-dimensional vector of values => grouped statistical test
# vector of values => grouped statistical test

# distance matrix => clustering
# distance matrix => dimensionality reduction
# distance matrix => vis

# dimensionality reduction => clustering
# dimensionality reduction => vis

# clustering => vis

# statistical test => vis
# grouped statistical test => vis

# distance matrix
# - cor
# - js
# - cor
# - cosine

# clustering
# - hclust
# - dbscan
# - kmeans

# dimensionality reduction
# - tsne
# - pca
# - mds

# stat test / grouped stat test
# - kruskall
# - wilcox

# postAnalysis <- function (.data)

# immunr_clustering_preprocessing <- function (...) {
# check for the right input class
# preprocess data somehow
# }


##########
# File: stat_tests.R
#########

#' #' Statistical analysis of groups
#' #'
#' immunr_permut <- function () {
#'   stop(IMMUNR_ERROR_NOT_IMPL)
#' }
#'
#' # groups
#' immunr_mann_whitney <- function () {
#'   stop(IMMUNR_ERROR_NOT_IMPL)
#' }
#'
#' immunr_kruskall <- function (.dunn = T) {
#'   stop(IMMUNR_ERROR_NOT_IMPL)
#' }
#'
#' immunr_logreg <- function () {
#'   stop(IMMUNR_ERROR_NOT_IMPL)
#' }

##########
# File: metadata.R
#########

# read_metadata <- function (.obj) {
#
# }
#
#
# write_metadata <- function (.obj) {
#
# }
#
#
# check_metadata <- function (.data, .meta) {
#   .meta = collect(.meta)
#
#   (length(.data) == length(unique(names(.data)))) &
#     (length(.meta$Sample) == length(unique(.meta$Sample))) &
#     (sum(!(names(.data) %in% .meta$Sample)) == 0)
# }
