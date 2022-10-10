#' @description
#' Description.
#'
#' @importFrom stats cor hclust as.dist cutree
#'
#' @param reaction_consistencies A param.
#' @param cluster_strength A param.
#' @param ... A param.
#'
#' @return An output.
#'
#' @importFrom magrittr %>% %<>%
#'
#' @noRd
get_metareactions <- function(reaction_consistencies, ..., cluster_strength) {
    reaction_id <- metareaction_code <- NULL
    pairwise_reaction_correlations <- stats::cor(t(reaction_consistencies), method = "spearman")
    pairwise_reaction_distances <- stats::as.dist(1 - pairwise_reaction_correlations)
    reaction_hierarchy <- stats::hclust(pairwise_reaction_distances, method = "complete")
    metareactions <-
        stats::cutree(reaction_hierarchy, h = cluster_strength) %>%
        tibble::enframe(name = "reaction_id", value = "metareaction_code") %>%
        dplyr::transmute(
            reaction_id = reaction_id,
            metareaction_id = paste("group", metareaction_code, sep = "_")
        )
    metareactions
}
