
r_lib <- '/home/r057600/R/x86_64-pc-linux-gnu-library/4.3'

data_dir <- '/mnt/data/genr/mrdata/GenR_MRI/summary_data/'
proj_dir <- '/home/r057600/mri_design_paper'


library(GenR.helpR)
library(ggplot2)

core <- file.path(data_dir, "core", "genr_mri_core_data_20260811.rds") |> 
  readRDS()
  
# F05, F09, F13, F17, F17U

dim(core) # 9901  111

data_overview(core)


age_plot <- function(cohort, name, y_max, y_invert = FALSE, time_period = 'life',
                     y_tick = 100, y_label_position = "right",  
                     subplot_space = 0.01, alpha = 0.95) {
  
  y_lims   <- c(0, y_max)
  y_breaks <- seq(0, y_max, by = y_tick)
  y_labs <- c('', y_breaks[-1]) # do not mark 0 count
  
  base_theme <- theme(axis.title = element_text(face = "bold"),
                      legend.position = "none",
                      base_family = "Helvetica",
          panel.grid.major = element_line(linewidth = 0.3, color = "grey90"),
          panel.grid.minor = element_line(linewidth = 0.2, color = "grey95"))
  
  if (time_period == 'birth') {
    x_lims = c(30, 45)
    x_label = "Gestational age at birth (weeks)"
    x = expr(Gest_age_birth)
    binwidth = 1 # one week
    
    y_label_position = "left"
    
    base_theme <- base_theme + theme(axis.title.y = element_blank())
  } else {
    x_lims = c(0, 27)
    x_label = "Age (years)"
    x = expr(Age)
    binwidth = 0.08 # one month
  }

  if (y_invert) {
    # Note: I use a standard continuous scale but with negative limits and breaks 
    # instead of scale_y_reverse (so stacks are correct)
    
    y_lims <- c(-y_max, 0)
    y_breaks <- -y_breaks
    
    xy_labs <- labs(x = x_label, y = name)
    
    theme_adjustment <- base_theme + 
      theme(plot.margin = margin(subplot_space, 5, 5, 5)) # add space on top
    
    # For inverted plot, map y to negative count so it stacks downwards
    y_aes <- expr(-after_stat(count))
    
  } else {
    
    xy_labs <- labs(y = name)
    
    theme_adjustment <- base_theme + 
      theme(axis.title.x = element_blank(), # remove x labs
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin = margin(5, 5, subplot_space, 5)) # add space on bottom
  
    y_aes <- expr(after_stat(count))
    
  }
  
  array_colors <- c("450K" = "#58aaa1", "EPICV1" = "#abc0da", "EPICV2" = "#1c4b75")
  # "#ba3f8d"
  
  # Build plot
  ggplot(cohort, aes(x = !!(x))) +
    geom_histogram(aes(fill = Array, y = !!y_aes),
                   alpha = alpha, binwidth = binwidth, position = "stack",
                   show.legend = TRUE) +
    scale_fill_manual(name = "Array", 
                      limits = names(array_colors), # always display all 3
                      values = array_colors, drop = FALSE) +
    guides(fill = guide_legend(override.aes = list(fill = unname(array_colors)))) + 
    scale_x_continuous(breaks = seq(0, x_lims[2], by = 2), 
                       limits = x_lims, expand = c(0, 0)) + 
    scale_y_continuous(limits = y_lims, breaks = y_breaks, labels = y_labs,
                       expand = c(0.01, 0.01), position = y_label_position) +
    geom_vline(xintercept = 0, color = "#f7f2b4", linewidth=1.2) +
    xy_labs +
    theme_minimal() + 
    theme_adjustment
}

# Main age plots

genr_max = 240 # 2300
alsp_max = 400 # 900

p_genr <- age_plot(genr |> filter(Period != "Birth"), 'Generation R', genr_max)
p_alsp <- age_plot(alsp |> filter(Period != "Birth"), 'ALSPAC', alsp_max, 
                   y_invert = TRUE)
 