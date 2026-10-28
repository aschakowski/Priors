###############################################################################
#
# Title: Visualize patch leaving results
#
# Description: 
#
# Authors: Schakowski, A.
#
# Last updated: 02/04/2024 (DD/MM/YYYY)
#
################################################################################
visualize_patch_leaving <- function(){
  
  # load model estimates
  patch_leaving_model_fit = readRDS(file = "data/processed_data/model_fit/patch_leaving_models/f_patch_leaving_model_prior_spatial_fit.rds")
  stan_data = readRDS(file = "data/processed_data/stan_data_patch_leaving.rds")
  
  # summarize as table
  #summary = patch_leaving_model_fit$summary(variables = c("beta", "sigma_lake", "sigma_trip", "sigma_id"), ~quantile(.x, probs = c(0.025, 0.975)))
  
  # parameter estimates and convergence diagnostics
  # parameter
  # beta1 = leave baseline
  # beta2 = patch discovery
  # beta3 = time without fish
  # beta4 = success feature
  # beta5 = loss feature
  # beta6 = prior feature
  # beta7 = social feature
  
  # visualize parameter estimates
  main_effects = patch_leaving_model_fit %>% 
    spread_draws(beta[parameter]) %>% 
    pivot_wider(names_from = parameter, values_from = beta)
  gc()
  
  lake_effects = patch_leaving_model_fit %>% 
    spread_draws(v_lake[lake, parameter])
  gc()
  
  lake_effects = left_join(lake_effects, main_effects)
  
  lake_effects = lake_effects %>% 
    mutate(lake_effect = ifelse(parameter == 1, v_lake + `1`, 
                                ifelse(parameter == 2, v_lake + `2`,
                                       ifelse(parameter == 3, v_lake + `3`,
                                              ifelse(parameter == 4, v_lake + `4`,
                                                     ifelse(parameter == 5, v_lake + `5`,
                                                            ifelse(parameter == 6, v_lake + `6`,
                                                                   ifelse(parameter == 7, v_lake + `7`))))))))
  
  # compute standardized predictors
  lake_effects$parameter = lake_effects$parameter %>% recode("Baseline", "Fish\nDiscovery", "Time w/o\nCatch (min)", "Successful", "Unsuccessful", "Prior", "Social")
  
  catch_bias = lake_effects %>% 
    mutate(name = factor(parameter, levels = rev(c("Baseline", "Time w/o\nCatch (min)", "Fish\nDiscovery")))) %>%
    filter(name %in% (c("Baseline", "Time w/o\nCatch (min)", "Fish\nDiscovery"))) %>% 
    mutate(lake_effect = ifelse(parameter == "Time w/o\nCatch (min)", lake_effect * 6, lake_effect)) %>% 
    #mutate(lake_effect = ifelse(parameter == "Baseline", logit(1 - (1-inv_logit(lake_effect))^6), lake_effect)) %>% #probability to not leave within 60 s based on 10 s hazard estimate
    ggplot(aes(x = lake_effect, y = name, group = as.factor(lake))) + 
    geom_vline(xintercept = 0, linetype = "dashed", size = .4) + 
    # .95
    annotate("rect", xmin = (quantile(main_effects$`2`, .025)), xmax =  (quantile(main_effects$`2`, .975)), ymin = 0.65, ymax = 1.35,
             alpha = 1,fill = "#BAE4B3")+
    #.66
    annotate("rect", xmin = (quantile(main_effects$`2`, .17)), xmax =  (quantile(main_effects$`2`, .83)), ymin = 0.65, ymax = 1.35,
             alpha = 1,fill = "#31A354")+
    # .5
    annotate("rect", xmin = (quantile(main_effects$`2`, .5))-0.005, xmax =  (quantile(main_effects$`2`, .5))+0.005, ymin = 0.65, ymax = 1.35,
             alpha = 1,fill = "black")+
    
    annotate("rect", xmin = (quantile(main_effects$`3`*6, .025)), xmax =  (quantile(main_effects$`3`*6, .975)), ymin = 1.65, ymax = 2.35,
             alpha = 1,fill = "#BAE4B3")+
    annotate("rect", xmin = (quantile(main_effects$`3`*6, .17)), xmax =  (quantile(main_effects$`3`*6, .83)), ymin = 1.65, ymax = 2.35,
             alpha = 1,fill = "#31A354")+
    annotate("rect", xmin = (quantile(main_effects$`3`*6, .5))-0.005, xmax =  (quantile(main_effects$`3`*6, .5))+0.005, ymin = 1.65, ymax = 2.35,
             alpha = 1,fill = "black")+
    
    #annotate("rect", xmin = (quantile(logit(1 - (1-inv_logit(main_effects$`1`))^6), .025)), xmax =  (quantile(logit(1 - (1-inv_logit(main_effects$`1`))^6), .925)), ymin = 2.65, ymax = 3.35,
    #         alpha = 1,fill = "#BAE4B3")+
    #annotate("rect", xmin = (quantile(logit(1 - (1-inv_logit(main_effects$`1`))^6), .17)), xmax =  (quantile(logit(1 - (1-inv_logit(main_effects$`1`))^6), .83)), ymin = 2.65, ymax = 3.35,
    #         alpha = 1,fill = "#31A354")+
    #annotate("rect", xmin = (quantile(logit(1 - (1-inv_logit(main_effects$`1`))^6), .49)), xmax =  (quantile(logit(1 - (1-inv_logit(main_effects$`1`))^6), .51)), ymin = 2.65, ymax = 3.35,
    #         alpha = 1,fill = "black")+
    
    annotate("rect", xmin = (quantile(main_effects$`1`, .025)), xmax =  (quantile(main_effects$`1`, .925)), ymin = 2.65, ymax = 3.35,
             alpha = 1,fill = "#BAE4B3")+
    annotate("rect", xmin = (quantile(main_effects$`1`, .17)), xmax =  (quantile(main_effects$`1`, .83)), ymin = 2.65, ymax = 3.35,
             alpha = 1,fill = "#31A354")+
    annotate("rect", xmin = (quantile(main_effects$`1`, .5))-0.005, xmax =  (quantile(main_effects$`1`, .5))+0.005, ymin = 2.65, ymax = 3.35,
             alpha = 1,fill = "black")+
    
    stat_pointinterval(position = position_dodge(width = .8), color = "black", alpha = .10, size = 0) +
    
    #
    #ylab("")+
    labs(tag = "A")+
    ggtitle("Resource Return")+
    xlab("Posterior Est. (95%CrI)\nlogit P(Leave) ~ Predictor")+
    ylab("Predictor")+
    #scale_x_continuous(breaks = seq(from = 0, to = 0.08, by = .02))+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6, angle = 55, hjust = .5),#, hjust = .5
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.9, .2),
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.background = element_rect(colour = "white", fill = "white"),
          text = element_text(family = "sans"))
  catch_bias
  
  spatial_weights = lake_effects %>% 
    filter(parameter %in% c("Successful","Unsuccessful", "Social", "Prior")) %>% 
    mutate(name = factor(parameter, levels = c("Prior","Social", "Unsuccessful", "Successful"))) %>%
    mutate(lake_effect = ifelse(parameter == "Prior", -lake_effect, lake_effect)) %>% 
    ggplot(aes(x = (lake_effect), y = name, group = as.factor(lake))) + 
    #stat_pointinterval(position = position_dodge(width = .6), color = "black", alpha = .25, size = 5) + 
    geom_vline(xintercept = 0, linetype = "dashed", size = .4) + 
    # .95
    annotate("rect", xmin = (quantile(main_effects$`7`, .025)), xmax =  (quantile(main_effects$`7`, .975)), ymin = 1.65, ymax = 2.35,
             alpha = 1,fill = "#BAE4B3")+
    #.66
    annotate("rect", xmin = (quantile(main_effects$`7`, .17)), xmax =  (quantile(main_effects$`7`, .83)), ymin = 1.65, ymax = 2.35,
             alpha = 1,fill = "#31A354")+
    # .5
    annotate("rect", xmin = (quantile(main_effects$`7`, .5))-0.001, xmax =  (quantile(main_effects$`7`, .5))+0.001, ymin = 1.65, ymax = 2.35,
             alpha = 1,fill = "black")+
    
    annotate("rect", xmin = (quantile(-main_effects$`6`, .025)), xmax =  (quantile(-main_effects$`6`, .975)), ymin = 0.65, ymax = 1.35,
             alpha = 1,fill = "#BAE4B3")+
    annotate("rect", xmin = (quantile(-main_effects$`6`, .17)), xmax =  (quantile(-main_effects$`6`, .83)), ymin = 0.65, ymax = 1.35,
             alpha = 1,fill = "#31A354")+
    annotate("rect", xmin = (quantile(-main_effects$`6`, .5))-0.001, xmax =  (quantile(-main_effects$`6`, .5))+0.001, ymin = 0.65, ymax = 1.35,
             alpha = 1,fill = "black")+
    
    annotate("rect", xmin = (quantile(main_effects$`5`, .025)), xmax =  (quantile(main_effects$`5`, .975)), ymin = 2.65, ymax = 3.35,
             alpha = 1,fill = "#BAE4B3")+
    annotate("rect", xmin = (quantile(main_effects$`5`, .17)), xmax =  (quantile(main_effects$`5`, .83)), ymin = 2.65, ymax = 3.35,
             alpha = 1,fill = "#31A354")+
    annotate("rect", xmin = (quantile(main_effects$`5`, .5))-0.001, xmax =  (quantile(main_effects$`5`, .5))+0.001, ymin = 2.65, ymax = 3.35,
             alpha = 1,fill = "black")+
    
    annotate("rect", xmin = (quantile(main_effects$`4`, .025)), xmax =  (quantile(main_effects$`4`, .975)), ymin = 3.65, ymax = 4.35,
             alpha = 1,fill = "#BAE4B3")+
    annotate("rect", xmin = (quantile(main_effects$`4`, .17)), xmax =  (quantile(main_effects$`4`, .83)), ymin = 3.65, ymax = 4.35,
             alpha = 1,fill = "#31A354")+
    annotate("rect", xmin = (quantile(main_effects$`4`, .5))-0.001, xmax =  (quantile(main_effects$`4`, .5))+0.001, ymin = 3.65, ymax = 4.35,
             alpha = 1,fill = "black")+
    
    stat_pointinterval(position = position_dodge(width = .8), color = "black", alpha = .10, size = 0) +
    ylab("Feature")+
    labs(tag = "B")+
    ggtitle("Spatial Features")+
    xlab("Posterior Est. (95%CrI)\nlogit P(Leave) ~ Feature")+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6, angle = 55, hjust = .5),#, hjust = .5
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.9, .2),
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.background = element_rect(colour = "white", fill = "white"),
          text = element_text(family = "sans"))
  spatial_weights
  
  
  
  # make panel
  panel_tmp_results = arrangeGrob(
    grobs = list(catch_bias, 
                 spatial_weights),
    widths = c(1, 1),
    heights = c(1),
    layout_matrix = rbind(c(1, 2))
  )
  
  
  ggsave(panel_tmp_results, file = "6_plots/panel_tmp_results.png", 
         width = 12,height = 6, bg = "white", dpi = 600, units = "cm")
  ggsave(panel_tmp_results, file = "6_plots/panel_tmp_results.svg", 
         width = 12,height = 6, bg = "white", dpi = 600, units = "cm", device = svglite::svglite())
  ggsave(panel_tmp_results, file = "6_plots/panel_tmp_results.pdf", 
         width = 12,height = 6, bg = "white", dpi = 600, units = "cm")
  
  gc()
  ################################################################################
  # visualize model comparison (supplement)
  
  # load comparison data
  comparison = readRDS(file = "data/processed_data/model_fit/patch_leaving_models/comparison/comparison.rds")
  comparison = as.data.frame(comparison)
  comparison$name = c("+Spatial\nSocial", "+Spatial\nPrior","+Spatial\nPersonal", "+Fish\nDiscovery", "+Time w/o\nCatch", "Baseline")
  comparison$model_id = c(2, 1, 3, 4, 5, 6)
  model_comparison = comparison %>% 
    ggplot(aes(x = elpd_diff, y = fct_reorder(name, model_id))) + 
    geom_col(color = "black", width = .6) +
    geom_errorbarh(aes(xmin = elpd_diff - 4*se_diff, xmax = elpd_diff + 4*se_diff), height = .3, size = .8) + 
    ggtitle("Model Comparison Spot-Leaving") + 
    xlab(expression(Delta~"ELPD (rel. to best)")) + 
    ylab("Model") + 
    #labs(tag = "C")+
    scale_x_continuous(labels = scales::comma_format(big.mark = ','))+
    geom_vline(xintercept = 0)+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7, angle = 45, hjust = .5),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.9, .2),
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.background = element_rect(colour = "white", fill = "white"),
          text = element_text(family = "sans"))
  model_comparison
  
  ggsave(model_comparison, file = "6_plots/model_comparison_tmp.png", 
         width = 7,height = 8, bg = "white", dpi = 300, units = "cm")
  ggsave(model_comparison, file = "6_plots/model_comparison_tmp.svg", 
         width = 7,height = 8, bg = "white", dpi = 300, units = "cm")
  ggsave(model_comparison, file = "6_plots/model_comparison_tmp.pdf", 
         width = 7,height = 8, bg = "white", dpi = 300, units = "cm")
  ggsave(model_comparison, file = "6_plots/model_comparison_tmp.tiff", 
         width = 7,height = 8, bg = "white", dpi = 300, units = "cm")
  
  ##############################################################################
  # confidence effects
  patch_leaving_model_fit_confidence = readRDS(file = "data/processed_data/model_fit/patch_leaving_models/g_patch_leaving_model_confidence_fit.rds")
  patch_leaving_model_fit_confidence$code()
  
  plot_s1_b = patch_leaving_model_fit_confidence %>% 
    spread_draws(beta[parameter]) %>%
    filter(parameter %in% c(7, 8, 9)) %>% 
    mutate(beta = -beta) %>% 
    mutate(parameter = ifelse(parameter == 7, "Likely", ifelse(parameter == 8, "Very Likely", "Definitely"))) %>% 
    mutate(parameter = factor(parameter, levels = rev(c("Definitely", "Very Likely", "Likely")))) %>% 
    ggplot(aes(y = parameter, x = beta)) + 
    stat_pointinterval(position = position_dodge(width = .8))+
    labs(tag = "A")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    #ggtitle("Confidence in Prior Beliefs")+
    xlab("Posterior Est. (95%CrI)\n logit P(Leave) ~ Prior | Confidence")+
    #scale_x_continuous(breaks = seq(from = -.3, to = .3, by = .1))+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 6),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5, angle = 65, hjust = .5),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = "none",
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 5, margin=margin(0,0,0,0)),
          legend.text = element_text(size = 5),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white", linewidth = .2),
          text = element_text(family = "sans"),
          legend.box.margin = margin(0, 0, 0, 0),
          legend.key.height = unit(.1, "cm"),
          plot.margin=unit(c(.1, .2, .1, .2), "cm"))
  plot_s1_b
  
  # contrasts of confidence model
  contrasts = patch_leaving_model_fit_confidence %>% 
    spread_draws(beta[parameter]) %>%
    filter(parameter %in% c(7, 8, 9)) %>% 
    mutate(beta = -beta) %>% 
    mutate(parameter = ifelse(parameter == 7, "Likely", ifelse(parameter == 8, "Very Likely", "Definitely"))) %>% 
    pivot_wider(names_from = parameter, values_from = beta) %>% 
    mutate(`Definitely - Very Likely` = Definitely - `Very Likely`,
           `Definitely - Likely` = Definitely - Likely,
           `Very Likely - Likely` = `Very Likely` - Likely)
  
  plot_s1_c = contrasts %>% 
    pivot_longer(cols = 7:9) %>% 
    mutate(name = factor(name, levels = rev(c("Definitely - Very Likely", "Definitely - Likely", "Very Likely - Likely")))) %>% 
    ggplot(aes(y = name, x = value)) + 
    stat_pointinterval() + 
    geom_vline(xintercept = 0, linetype = "dashed")+
    xlab("Posterior Est. (95%CrI)\nContrasts logit P(Leave) ~ Prior | Confidence")+
    labs(tag = "B") +
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 6),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5, angle = 65, hjust = .5),#, hjust = .5
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = "none",
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 5, margin=margin(0,0,0,0)),
          legend.text = element_text(size = 5),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white", linewidth = .2),
          text = element_text(family = "sans"),
          legend.box.margin = margin(0, 0, 0, 0),
          legend.key.height = unit(.1, "cm"),
          plot.margin=unit(c(.1, .2, .1, .2), "cm"))
  plot_s1_c
  
  # make panel
  panel_conf = gridExtra::arrangeGrob(
    grobs = list(plot_s1_b,
                 plot_s1_c),
    widths = c(1, 1),
    layout_matrix = rbind(c(1,2))
  )
  ggsave(panel_conf, file = "6_plots/plot_confidence_tmp.png", 
         width = 12, height = 6, bg = "white", dpi = 300, units = "cm")
  ggsave(panel_conf, file = "6_plots/plot_confidence_tmp.pdf", 
         width = 12, height = 6, bg = "white", dpi = 300, units = "cm")
  ggsave(panel_conf, file = "6_plots/plot_confidence_tmp.svg", 
         width = 12, height = 6, bg = "white", dpi = 300, units = "cm", device = svglite::svglite)
  
}


################################################################################
# END
################################################################################