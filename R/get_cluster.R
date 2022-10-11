#' @title avg_sil
#' @export
#' @description
#' Given a matrix or data.frame with point coordinates and the number of clusters to
#' get, the average silhouette is computed.
#'
#' @importFrom cluster silhouette
#' 
#' @inheritParams get_cluster
#'
#' @return the average silhouette index.
#'
#' @importFrom magrittr %>% %<>%
#'
#' @examples 
#' set.seed(123)
#' mat <- matrix(c(
#'     rnorm(n = 100, mean = 10),
#'     rnorm(n = 100, mean = 20),
#'     rnorm(n = 100, mean = 1),
#'     rnorm(n = 100, mean = 2)), ncol = 2)
#' avg_sil(mat, 2)

avg_sil <- function(df, k, method = "average") {
    hc_res <- get_cluster(df = df, k = k, method = method, 
        onlyMembership = FALSE)
    hc <- hc_res$hc
    d <- hc_res$d
    ss <- cluster::silhouette(cutree(hc, k), d)
    mean(ss[, 3])
}

#' @title get_optimal_clusters
#' @export
#' @description
#' Compute the optimal number of clusters for a given set of UMAP coordinates.
#'
#' @import ggplot2
#' @importFrom purrr map_dbl
#'
#' @param UMAPcoords a matrix or data.frame with the dimensional reduction coordinates.
#' @param k.values a vector for the number of clusters to test.
#' @param plotIt logical. plot or don't plot the output.
#' @inheritParams get_cluster
#'
#' @return An optional graphical representation of the average silhouette index 
#' for several clustering configurations and the optimal number of clusters.
#' 
#' @examples
#' set.seed(123)
#' mat <- matrix(c(
#'     rnorm(n = 100, mean = 10), 
#'     rnorm(n = 100, mean = 20),
#'     rnorm(n = 100, mean = 1), 
#'     rnorm(n = 100, mean = 2)), ncol = 2)
#' get_optimal_clusters(mat)

get_optimal_clusters <- function (UMAPcoords, k.values = 2:15, 
    method = "average", plotIt = TRUE) {
    # extract avg silhouette for 2-15 clusters
    avg_sil_values <- purrr::map_dbl(k.values, avg_sil, df = UMAPcoords)
    df_to_plot <- data.frame("n_cluster" = k.values, 
                             "avg_silhouette" = avg_sil_values)
    if (plotIt) {
        p <- ggplot(df_to_plot, aes(x = n_cluster, y = avg_silhouette)) + 
            geom_line() + 
            geom_point() + 
            ggtitle("Average Silhouette indexes", 
                    subtitle = "Hierarchical clustering on euclidean distances") +
            xlab("Number of clusters") + ylab("Average Silhouette")
        print(p)
    }
    
    optimalClusterNum <- which(avg_sil_values == max(avg_sil_values))
    
    return (k.values[optimalClusterNum])
}

#' @title get_cluster
#' @export
#' @description
#' Given a matrix or data.frame with point coordinates and the number of clusters to
#' get returns the cluster membership of each point or the cluster object with the
#' distance matrix.
#'
#' @importFrom stats dist hclust cutree
#' 
#' @param df a matrix or data.frame with coordinates.
#' @param k the desired number of clusters.
#' @inheritParams stats::hclust
#' @param onlyMembership logical. Return only the class membership or the cluster
#' object.
#'
#' @return Vector of cluster membership or a list with the cluster object and the
#' distance matrix.
#'
#' @examples 
#' set.seed(123)
#' mat <- matrix(c(
#'     rnorm(n = 100, mean = 10), 
#'     rnorm(n = 100, mean = 20),
#'     rnorm(n = 100, mean = 1), 
#'     rnorm(n = 100, mean = 2)), ncol = 2)
#' get_cluster(mat, 2)

get_cluster <- function(df, k, method = "average", onlyMembership = TRUE){
    df <- na.omit(df)
    d <- stats::dist(df)
    hc <- stats::hclust(d, method = method)
    if(onlyMembership){
        membs <- stats::cutree(hc, k = k)
        return(membs)
    } else {
        return(list("hc" = hc, "d" = d))
    }
}
