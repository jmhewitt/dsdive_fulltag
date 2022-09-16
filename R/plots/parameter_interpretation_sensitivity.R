parameter_interpretation_sensitivity_plot_script = tar_target(
  name = parameter_interpretation_sensitivity_plot,
  command = {
   
    df = parameter_interpretation_plot
    
    alpha = .2
    
    # visualize parameter effects across multiple starts
    pl = ggplot(df %>% 
                  mutate(rep = factor(rep)), 
                aes(x = time, y = first_deep, group = rep, 
                    col = init_stage, 
                    lty = init_stage, pch = init_stage)) + 
      geom_point(alpha = alpha) +
      geom_line(alpha = alpha) + 
      scale_linetype_discrete('Initial movement', labels = titlecase) + 
      scale_shape_discrete('Initial movement', labels = titlecase) + 
      scale_color_brewer('Initial movement', labels = titlecase, 
                         type = 'qual', palette = 'Dark2') + 
      xlab('Lighting conditions') + 
      ylab('Deep depth hitting time (min)') + 
      facet_grid(init_stage~prototype, labeller = titlecase) + 
      theme_few()
    
    # save plot
    
    f = file.path('output', 'figures', 'parameter_interpretation_sensitivity')
    dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
    
    ggsave(pl, filename = file.path(f, paste(tar_name(), '.pdf', sep = '')),
           width = 8, height = 8, dpi = 'print')
    
    
  }
)