#' @title compare_clusters
#' @export
#' @description
#' Given two groups of cells the function returns several objects to study the 
#' differentially consistent reactions.
#'
#' @importFrom dplyr summarise
#' @importFrom DT datatable
#' @importFrom reshape2 melt
#' @importFrom tidytext reorder_within scale_y_reordered
#' @importFrom magrittr %>% 
#' @import ggplot2
#' 
#' @param cell_info the object containing the cell's metadata information
#' @param compass_analyzer the object created by CompassAnalyzer$new
#' @param compass_data the object created by CompassData$new
#' @param variable the name of the variable to test
#' @param cluster_A,cluster_B vector of labels for the chosen variable levels 
#' to compare
#' @param for_metareactions logical to decide wether to perform tests for 
#' reactions or metareactions
#' @param adjusted_p_value_th adjusted p-value threshold to detect significant 
#' reactions (default 0.05)
#' @param cohens_d_th cohen's D statistic threshold to detect effect-size 
#' significant reactions (default 1)
#'
#' @return a list of results.
#' \itemize{
#'     \item{\code{wilcoxon_results_with_metadata}}{ wilcoxon test results}
#'     \item{\code{DC_core_metabolism}}{ table of differentially consistent 
#'     reactions or metareactions}
#'     \item{\code{DC_core_metabolism_table}}{ ready-to-plot table of
#'     differentially consistent reactions or metareactions}
#'     \item{\code{p_DC_reactions}}{ plot for the differentially consistent 
#'     reactions when \code{for_metareactions = TRUE}}
#'     \item{\code{p_DC_reactions_cohensd}}{ plot for the effect sizes of the 
#'     differentially consistent reactions when \code{for_metareactions = TRUE}}}

compare_clusters <- function(
        cell_info,
        compass_analyzer, 
        compass_data, 
        variable,
        cluster_A, cluster_B, 
        for_metareactions = FALSE,
        adjusted_p_value_th = 0.05,
        cohens_d_th = 1)
{
    # Get cell_ids
    group_A_cell_ids <-
        cell_info %>%
        filter(get(variable) %in% cluster_A) %>%
        pull(cell_id)
    group_B_cell_ids <-
        cell_info %>%
        filter(get(variable) %in% cluster_B) %>%
        pull(cell_id)
    
    # Perform wilcox test
    wilcoxon_results <- compass_analyzer$conduct_wilcoxon_test(
        consistencies_matrix =
            if(for_metareactions) {
                compass_data$metareaction_consistencies 
            } else compass_data$reaction_consistencies,
        group_A_cell_ids = group_A_cell_ids,
        group_B_cell_ids = group_B_cell_ids,
        for_metareactions = for_metareactions
    ) %>% # Adding significant informations
        mutate(
            group = ifelse(adjusted_p_value <= adjusted_p_value_th, 
                           yes = ifelse(cohens_d >= cohens_d_th, 
                                        yes = "Significant and Effect Size",
                                        no = "Significant"),
                           no = ifelse(cohens_d >= cohens_d_th, 
                                       yes = "Effect Size",
                                       no = "Not Significant"))
        )
    # Add metadata        
    wilcoxon_results_with_metadata <- wilcoxon_results %>%
        inner_join( # Add reaction partitions
            compass_data$reaction_partitions,
            by = ifelse(test = for_metareactions, 
                        yes = "metareaction_id", 
                        no = "reaction_id")
        ) %>%
        inner_join( # Add reaction metadata
            compass_data$reaction_metadata,
            by = "reaction_no_direction"
        ) %>% 
        mutate( # Create core reaction
            core = confidence %in% c(0,4) & !is.na(EC_number)
        )
    
    # Checking core metareactions
    core_metareaction_df <- wilcoxon_results_with_metadata %>%
        group_by(metareaction_id) %>%
        dplyr::summarise(core_metareaction = ifelse(
            test = sum(core) > 0, 
            yes = "Core Metabolism",
            no = "Not Core Metabolism"))
    
    # Adding core metareactions
    wilcoxon_results_with_metadata <- wilcoxon_results_with_metadata %>%
        inner_join(
            core_metareaction_df,
            by = "metareaction_id"
        )
    
    # Changing core levels
    wilcoxon_results_with_metadata$core = factor(
        wilcoxon_results_with_metadata$core, 
        levels = c(TRUE, FALSE), 
        labels = c("Core Metabolism", "Not Core Metabolism"))
    
    if(for_metareactions){
        # Counting the number of reactions for each metareaction
        # Counting the number of subsystem for each metareaction
        reaction_metareaction_df <- wilcoxon_results_with_metadata %>%
            group_by(metareaction_id) %>%
            dplyr::summarise(
                reaction_number = length(unique(reaction_id)),
                subsystem_number = length(unique(subsystem))) 
        # Extract core metabolism metareactions
        DC_core_metabolism <- wilcoxon_results_with_metadata %>% 
            filter(adjusted_p_value < adjusted_p_value_th & core == "Core Metabolism") %>%
            select(metareaction_id, cohens_d, adjusted_p_value) %>%
            arrange(adjusted_p_value) %>%
            left_join(reaction_metareaction_df, by = "metareaction_id") %>%
            distinct()
    } else {
        # Extract core metabolism reactions
        DC_core_metabolism <- wilcoxon_results_with_metadata %>% 
            filter(adjusted_p_value < adjusted_p_value_th & core == "Core Metabolism") %>%
            select(reaction_id, reaction_name, subsystem,
                   associated_genes, cohens_d, adjusted_p_value) %>%
            arrange(adjusted_p_value)
    }
    table_cap <- paste0(
        "Group: ",
        paste0(cluster_A, collapse = ", "),
        " vs Group: ",
        paste0(cluster_B, collapse = ", "),
        " - Differentially consistent core ",
        ifelse(for_metareactions, yes = "metareactions", no = "reactions")
    )
    DC_core_metabolism_table <- DT::datatable(DC_core_metabolism, 
                                              caption = table_cap, filter = "top", extensions = 'Buttons',
                                              options = list(dom = 'Blfrtip', buttons = c('copy', 'csv', 'excel'),
                                                             lengthMenu = list(c(10,25,50,-1), c(10,25,50,"All"))))
    ##### OUTPUTS #####
    out <- list(
        "wilcoxon_results_with_metadata" = wilcoxon_results_with_metadata, 
        "DC_core_metabolism" = DC_core_metabolism,
        "DC_core_metabolism_table" = DC_core_metabolism_table)
    
    ##### PLOTS #####
    
    if(!for_metareactions){
        sig_subsystem <- wilcoxon_results_with_metadata %>% 
            filter(adjusted_p_value <= adjusted_p_value_th) 
        
        sig_subsystem_summary <- sig_subsystem %>%
            group_by(subsystem, core) %>%
            dplyr::summarise(
                pos = sum(cohens_d > 0),
                neg = -sum(cohens_d < 0))
        
        p_DC_reactions <- ggplot(data = reshape2::melt(sig_subsystem_summary), 
                                 mapping = aes(
                                     x = value, 
                                     y = tidytext::reorder_within(
                                         x = subsystem, 
                                         by = value, 
                                         within = core), 
                                     fill = variable)) +
            facet_wrap(facets = ~ core, nrow = 1, scales = "free") + 
            geom_col(orientation = "y") +
            geom_text(aes(label = ifelse(value > 0, value, "")), hjust = "right") +
            geom_text(aes(label = ifelse(value < 0, abs(value), "")), hjust = "left") +
            scale_x_continuous(trans = ggallin::pseudolog10_trans) +
            scale_y_reordered() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) + 
            theme(legend.position = "top") +
            ylab("Subsystem") + xlab("Number of reactions") +
            ggtitle(
                label = "Number of differentially consistent reactions", 
                subtitle = "Colored by direction and stratified by core metabolism") +
            scale_fill_discrete(name = "More consistent in:", 
                                labels = c(paste0(cluster_A, collapse = ", "), paste0(cluster_B, collapse = ", ")))
        
        p_DC_reactions_cohensd <- ggplot(data = sig_subsystem, 
                                         mapping = aes(
                                             x = cohens_d, 
                                             y = reorder_within(
                                                 x = subsystem, 
                                                 by = cohens_d, 
                                                 within = core, 
                                                 fun = median), 
                                             color = as.factor(sign(cohens_d)))) +
            facet_wrap(~ core, scales = "free") +
            scale_y_reordered() +
            geom_jitter(width = 0, height = 0.3, alpha = 0.5) +
            geom_vline(aes(xintercept = 0), lty = 2, color = "grey") +
            theme(legend.position = "top") +
            ylab("Subsystem") + xlab("Cohen's D Effect Size") +
            ggtitle(
                label = "Effect size of the differentially consistent reactions", 
                subtitle = "Colored by direction and stratified by core metabolism, ordered by median Cohen's D") +
            scale_color_discrete(name = "More consistent in:", 
                                 labels = c(paste0(cluster_A, collapse = ", "), paste0(cluster_B, collapse = ", ")))
        out <- append(out, values = list(
            "p_DC_reactions" = p_DC_reactions,
            "p_DC_reactions_cohensd" = p_DC_reactions_cohensd))
    }
    return(out)
}

#' @title plot_metareactions
#' @export
#' @description
#' Given the wilcoxon results generated by the compare_clusters() function and
#' the name of some metareactions, the alluvium plot of their metabolism and 
#' subsystems is generated.
#'
#' @importFrom dplyr summarise
#' @importFrom ggrepel geom_text_repel 
#' @importFrom magrittr %>% 
#' @import ggplot2 ggalluvial
#' 
#' @param wilcoxon_results_with_metadata the object created by the 
#' compare_clusters() function
#' @param met_id a vector containing the metareaction names
#'
#' @return an alluvium plot.

plot_metareaction <- function(
        met_id, 
        wilcoxon_results_with_metadata){
    
    subset_df_summary <- wilcoxon_results_with_metadata %>%
        filter(metareaction_id %in% met_id) %>%
        group_by(metareaction_id, subsystem, core) %>%
        dplyr::summarise(
            Freq = length(unique(reaction_id))) %>%
        arrange(Freq)
    
    p <- ggplot(subset_df_summary, 
                aes(y = Freq, axis1 = core, axis2 = metareaction_id, axis3 = subsystem)) +
        geom_alluvium(aes(fill = as.factor(metareaction_id))) +
        geom_stratum() +
        geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 4, min.y = 5) +
        geom_text_repel(stat = "stratum", aes(label = after_stat(stratum)), size = 4, max.y = 5, nudge_x = 1) +
        scale_x_discrete(limits = c("Metabolism", "Metareaction ID", "Subsystem"), expand = c(0.5, 0)) +
        theme(legend.position = "none") +
        ggtitle(label = "Alluvium represetation for selected Metareactions", subtitle = "by metabolism type and subsystem")
    
    return(p)
}