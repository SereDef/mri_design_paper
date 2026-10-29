
# ResearchSuite ----------------------------------------
# r_lib <- '/home/r057600/R/x86_64-pc-linux-gnu-library/4.3'

# data_dir <- '/mnt/data/genr/mrdata/GenR_MRI/summary_data/'
# proj_dir <- '/home/r057600/mri_design_paper'

data_dir <- '..'

library(genR.helpR)
library(ggplot2)

core <- file.path(data_dir, "core", "genr_mri_core_data_20260811.rds") |> 
  readRDS()
  
# F05, F09, F13, F17, F17U

dim(core) # >> 9901  111

data_overview(core)

# quick helper: get a set of variables 
gv <- function(var) {
  waves <- paste0('f', c('05', '09', '13', '17'))
  var_names <- paste(var, waves, sep ='_')
  
  return(var_names)
}

lapply(core[, gv('mri_consent')], table) 

summary(core[, gv('age_child_mri')])

names(core)

core_long <- core |>
  dplyr::select(idc, matches("_f(05|09|13|17)$")) |>
  tidyr::pivot_longer(
    cols = -idc,
    names_to  = c(".value", "time"),
    names_pattern = "^(.*)_f(\\d{2})$" # stem = everything before the last _fage
  ) |> 
  dplyr::mutate(time = factor(time, levels = c("05", "09", "13", "17"), 
                              labels = c(' 8 y', '10 y', '14 y', '18 y')))

age_plot <- function(cohort, name, y_max, y_tick = 100, y_label_position = "right",  
                     subplot_space = 0.01, alpha = 0.95) {
  
  y_lims   <- c(0, y_max)
  y_breaks <- seq(0, y_max, by = y_tick)
  y_labs <- c('', y_breaks[-1]) # do not mark 0 count
  
  base_theme <- theme(axis.title = element_text(face = "bold"),
                      legend.position = "left",
                      base_family = "Helvetica",
          panel.grid.major = element_line(linewidth = 0.3, color = "grey90"),
          panel.grid.minor = element_line(linewidth = 0.2, color = "grey95"))
  
  x_lims = c(0, 27)
  x_label = "Age (years)"
  x = expr(age_child_mri)
  binwidth = 0.08 # one mont
  
  xy_labs <- labs(x = x_label, y = name)
    
  theme_adjustment <- base_theme + 
    theme(# axis.title.x = element_blank(), # remove x labs
          # axis.text.x = element_blank(),
          # axis.ticks.x = element_blank(),
          plot.margin = margin(5, 5, subplot_space, 5)) # add space on bottom

  y_aes <- expr(after_stat(count))
    
  
  time_colors <- c(" 8 y" = "salmon", "10 y" = "#58aaa1", "14 y" = "#abc0da", "18 y" = "#1c4b75")
  # "#ba3f8d"
  
  # Build plot
  ggplot(cohort, aes(x = !!(x))) +
    geom_histogram(aes(fill = time, y = !!y_aes), # # "#abc0da",
                   alpha = alpha, binwidth = binwidth, 
                   position = position_stack(reverse = TRUE),
                   show.legend = TRUE) +
    scale_fill_manual(name = "Visit", 
                      limits = names(time_colors), # always display all 3
                      values = time_colors, drop = FALSE) +
    guides(fill = guide_legend(override.aes = list(fill = unname(time_colors)))) + 
    scale_x_continuous(breaks = seq(0, x_lims[2], by = 2), 
                       limits = x_lims, expand = c(0, 0)) + 
    scale_y_continuous(limits = y_lims, breaks = y_breaks, labels = y_labs,
                       expand = c(0.01, 0.01), position = y_label_position) +
    # geom_vline(xintercept = 0, color = "#f7f2b4", linewidth=1.2) +
    xy_labs +
    theme_minimal() + 
    theme_adjustment
}
# Main age plots

p_ages <- age_plot(core_long, 'Counts', y_max = 500)
p_ages

ggsave(filename = "mri_age_counts.png", plot = p_ages,
       width = 12, height = 4.5, dpi = 300)

# ===============================================================================

library(ComplexUpset)

waves <- c("8y", "10y", "14y", "18y")   # f05, f09, f13, f17 – adjust labels if needed

followup_color_map <- c(
  `1` = '#FFC541', 
  `2` = '#FFC3CB',
  `3` = '#F95493',
  `4` = '#B8168A',
  `5` = '#9D4EDD', 
  `6` = '#3A0CA3' )

time_color_map <- c("8y" = "salmon", "10y" = "#58aaa1", "14y" = "#abc0da", "18y" = "#1c4b75")

make_upset <- function(criterium, waves = c("8y", "10y","14y","18y"), yes = 'yes') {
  
  vars <- intersect(gv(criterium), names(core))
  
  df <- as.data.frame(lapply(core[, vars], function(v) v %in% yes))
  names(df) <- waves
  df <- df[rowSums(df) > 0, ] # drop children with no MRI visit at all
  
  intersections <- unique(upset_data(df, 
                                     rev(waves))$plot_intersections_subset)
  
  followup_query <- lapply(intersections, function (i) {
    i_vec <- strsplit(i, '-')[[1]]
    
    i_color <- followup_color_map[[as.character(length(i_vec))]]
    
    upset_query(intersect = i_vec, color= i_color, fill = i_color,
                only_components=c('intersections_matrix', 'Intersection size'))
  })
  
  size = get_size_mode('exclusive_intersection')
  
  upset_mri <- ComplexUpset::upset(
    df, rev(waves),
    name = "Visit",
    sort_sets = FALSE,
    sort_intersections_by = c("degree", "cardinality"),
    # queries = list( # reuse your visit colours on the set-size bars
    #   upset_query(set = "8y",  fill = "salmon"),
    #   upset_query(set = "10y", fill = "#58aaa1"),
    #   upset_query(set = "14y", fill = "#abc0da"),
    #   upset_query(set = "18y", fill = "#1c4b75")
    # ),
    queries = followup_query,
    width_ratio = 0.2, height_ratio = 0.8, wrap=TRUE,
    set_sizes = (
      upset_set_size()
      + geom_text(aes(label = ..count.., y = ..count.. - 1000),
                  hjust = 1.1,
                  color = "white",
                  stat = 'count')
      + expand_limits(y = 4000)
      + theme(axis.text.x = element_blank())),
    base_annotations = list(
      'Intersection size' = intersection_size(
        counts = TRUE,
        text_mapping = aes(
          # label = paste0(!!size,'\n(',!!upset_text_percentage(),')'), 
          label = !!size,
          color = 'on_background',
          y = ifelse(!!size > 500, !!size - (!!size*0.1), !!size)) # ,
        # mapping = aes(fill = Array)) + 
      ) +
        scale_fill_manual(values = time_color_map))) + 
    ggtitle(criterium) + 
    theme(plot.title = element_text(face = "bold", 
                                    size = 14), 
          plot.margin = margin(t = 10, r = 10, b = 10, l = 40))
  ggsave(filename = paste0("mri_repmeas_counts_",criterium,".png"), plot = upset_mri,
         width = 12, height = 7, dpi = 300, bg = "white")
  
}


make_upset('mri_consent')
make_upset('t1_has_nii')
make_upset('dti_has_nii')
make_upset('rsfmri_has_nii')
make_upset('asl_has_nii', waves=c('14y','18y'))
make_upset('usability', yes = 'usable')


# TODO: show trajectory for some global-ish measure 


 