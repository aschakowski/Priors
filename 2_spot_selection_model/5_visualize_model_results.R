################################################################################
#
# TITLE: Visualize spot selection model
#
# AUTHORS: Schakowski, A. (schakowski@mpib-berlin.mpg.de)
#
# DESCRIPTION: This script visualizes the results of the spot selection model 
# (main text figures + supplement)
#
# LAST UPDATED (DD/MM/YYYY): 10/10/2026
#
################################################################################
visualize_model_results <- function(){
  
  ################################################################################
  # simulate fitted features
  ################################################################################
  # first need to compute the spatial features based on fitted bandwidth parameters
  # sourcing helper function to do so
  source("3_spot_selection_model/scripts/simulate_fitted_features.r")
  simulate_fitted_features()
  
  ################################################################################
  # load data
  ################################################################################
  # load spatial features
  fitted_features = readRDS(file = "data/processed_data/fitted_features.rds")
  
  # load model fit
  spot_selection_model_fit = readRDS("data/processed_data/model_fit/spot_selection_model/i_temporal_fit.rds")
 
  # spot data containing all angling spots
  spot_selection_data = fread(file = "data/spot_selection_data_cleaned.csv")
  
  # idnetifiers
  identifiers_spot_selection = readRDS("data/processed_data/identifiers_spot_selection_model_key.rds")
  
  # load ecological knowledge
  ecological_knowledge = readRDS(file = "data/ecological_knowledge_polygons.rds")
  
  ################################################################################
  # Visualize feature weights
  ################################################################################
  spot_selection_model_fit$code()
  lake_effects = spot_selection_model_fit %>% 
    spread_draws(weights[feature], v_lake[lake, feature]) %>% 
    mutate(lake_effect = weights + v_lake)
  lake_effects$feature = lake_effects$feature %>% recode("Locality","Locality x\nTime", "Successful","Unsuccessful", "Social","Prior", 
                                                         "Successful x\nTime", "Unsuccessful x\nTime", "Social x\nTime", "Prior x\nTime", 
                                                         "Prior x\nSuccessful", "Prior x\nUnsuccessful", "Prior x\nSocial")
  
  # compute standardized predictors
  scale_features = fitted_features %>% 
    group_by(trip, lake) %>% 
    summarize(sd_social = sd(social_feature),
              sd_success = sd(success_feature),
              sd_loss = sd(loss_feature),
              sd_locality = sd(locality_feature),
              sd_prior = sd(knowledge_feature),
              sd_prior_success = sd(knowledge_feature)*sd(success_feature),
              sd_prior_loss = sd(knowledge_feature)*sd(loss_feature),
              sd_prior_social = sd(knowledge_feature)*sd(social_feature)) %>% 
    ungroup() %>% 
    group_by(lake) %>% 
    summarize(sd_social = mean(sd_social),
              sd_prior = mean(sd_prior),
              sd_success = mean(sd_success),
              sd_loss = mean(sd_loss),
              sd_locality = mean(sd_locality),
              sd_prior_success = mean(sd_prior_success),
              sd_prior_loss = mean(sd_prior_loss),
              sd_prior_social = mean(sd_prior_social))
  
  # compute standardized predictors
  lake_effects = left_join(lake_effects, scale_features)
  
  lake_effects = lake_effects %>% 
    mutate(standardized_weight = ifelse(feature == "Social", lake_effect*sd_social, 
                                        ifelse(feature == "Prior", lake_effect*sd_prior, 
                                               ifelse(feature == "Successful", lake_effect*sd_success, 
                                                      ifelse(feature == "Unsuccessful", lake_effect*sd_loss,
                                                             ifelse(feature == "Prior x\nUnsuccessful", lake_effect*sd_prior_loss,
                                                                    ifelse(feature == "Prior x\nSuccessful", lake_effect*sd_prior_success,
                                                                           ifelse(feature == "Prior x\nSocial", lake_effect*sd_prior_social, 
                                                                                  ifelse(feature == "Locality", lake_effect*sd_locality, 
                                                                                         ifelse(feature == "Locality x\nTime", lake_effect*sd_locality,
                                                                                                ifelse(feature == "Social x\nTime", lake_effect*sd_social,
                                                                                                       ifelse(feature == "Prior x\nTime", lake_effect*sd_prior,
                                                                                                              ifelse(feature == "Successful x\nTime", lake_effect*sd_success,
                                                                                                                     ifelse(feature == "Unsuccessful x\nTime", lake_effect*sd_loss, 0))))))))))))),
           standardized_main_effect = ifelse(feature == "Social", (weights)*mean(scale_features$sd_social), 
                                             ifelse(feature == "Prior", (weights)*mean(scale_features$sd_prior), 
                                                    ifelse(feature == "Successful", (weights)*mean(scale_features$sd_success), 
                                                           ifelse(feature == "Unsuccessful", weights*mean(scale_features$sd_loss),
                                                                  ifelse(feature == "Prior x\nUnsuccessful", weights * mean(scale_features$sd_prior_loss),
                                                                         ifelse(feature == "Prior x\nSuccessful", weights * mean(scale_features$sd_prior_success),
                                                                                ifelse(feature == "Prior x\nSocial", weights*mean(scale_features$sd_prior_social), 
                                                                                       ifelse(feature == "Locality", weights*mean(scale_features$sd_locality), 
                                                                                              ifelse(feature == "Locality x\nTime", weights*mean(scale_features$sd_locality),
                                                                                                     ifelse(feature == "Social x\nTime", weights*mean(scale_features$sd_social),
                                                                                                            ifelse(feature == "Prior x\nTime", weights*mean(scale_features$sd_prior),
                                                                                                                   ifelse(feature == "Successful x\nTime", weights*mean(scale_features$sd_success),
                                                                                                                          ifelse(feature == "Unsuccessful x\nTime", weights*mean(scale_features$sd_loss),0))))))))))))))
  # visualize estimates
  model_results = lake_effects %>% 
    filter(feature %in% c("Successful", "Unsuccessful","Prior", "Social")) %>% 
    mutate(name = factor(feature, levels=c("Prior","Social", "Unsuccessful", "Successful"))) %>%
    mutate(standardized_weight = ifelse(feature %in% c("Prior"), -standardized_weight, standardized_weight)) %>% 
    ggplot(aes(x = standardized_weight, y = name, group = as.factor(lake))) + 
    geom_vline(xintercept = 0, linetype = "dashed", size = .4) + 
    # .95
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Successful"], .025), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Successful"], .975), ymin = 3.65, ymax = 4.35,
             alpha = 1,fill = "#BAE4B3")+
    #.66
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Successful"], .17), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Successful"], .83), ymin = 3.65, ymax = 4.35,
             alpha = 1,fill = "#31A354")+
    # .5
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Successful"], .5)-0.005, xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Successful"], .5)+0.005, ymin = 3.65, ymax = 4.35,
             alpha = 1,fill = "black")+
    
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social"], .025), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social"], .975), ymin = 1.65, ymax = 2.35,
             alpha = 1,fill = "#BAE4B3")+
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social"], .17), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social"], .83), ymin = 1.65, ymax = 2.35,
             alpha = 1,fill = "#31A354")+
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social"], .5)-0.005, xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social"], .5)+0.005,ymin = 1.65, ymax = 2.35,
             alpha = 1,fill = "black")+
    
    annotate("rect", xmin = quantile(-lake_effects$standardized_main_effect[lake_effects$feature == "Prior"], .025), xmax =  quantile(-lake_effects$standardized_main_effect[lake_effects$feature == "Prior"], .975), ymin = 0.65, ymax = 1.35,
             alpha = 1,fill = "#BAE4B3")+
    annotate("rect", xmin = quantile(-lake_effects$standardized_main_effect[lake_effects$feature == "Prior"], .17), xmax =  quantile(-lake_effects$standardized_main_effect[lake_effects$feature == "Prior"], .83), ymin = 0.65, ymax = 1.35,
             alpha = 1,fill = "#31A354")+
    annotate("rect", xmin = quantile(-lake_effects$standardized_main_effect[lake_effects$feature == "Prior"], .5)-0.005, xmax =  quantile(-lake_effects$standardized_main_effect[lake_effects$feature == "Prior"], .5)+0.005,ymin = 0.65, ymax = 1.35,
             alpha = 1,fill = "black")+
    
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Unsuccessful"], .025), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Unsuccessful"], .975), ymin = 2.65, ymax = 3.35,
             alpha = 1,fill = "#BAE4B3")+
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Unsuccessful"], .17), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Unsuccessful"], .83), ymin = 2.65, ymax = 3.35,
             alpha = 1,fill = "#31A354")+ 
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Unsuccessful"], .5)-0.005, xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Unsuccessful"], .5)+0.005, ymin = 2.65, ymax = 3.35,
             alpha = 1,fill = "black")+
    ylab("Feature")+
    stat_pointinterval(position = position_dodge(width = .8), color = "black", alpha = .10, size = 0) +
    
    #stat_pointinterval(color = "black", alpha = .15) +
    labs(tag = "D", color = "value")+
    ggtitle("Feature Weights")+
    xlab("Posterior Est. (95%CrI)\nSpot Selection ~ Feature")+
    scale_x_continuous(breaks = seq(from = -3, to = 6, by = 1))+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 6),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5, angle = 55, hjust = .5),#, hjust = .5
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.9, .2),
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.background = element_rect(colour = "white", fill = "white"),
          text = element_text(family = "sans"))
  model_results
  
  ################################################################################
  # Visualize interaction effects
  ################################################################################
  model_results_interaction = lake_effects %>% 
    filter(feature %in% c("Prior x\nUnsuccessful", "Prior x\nSuccessful", "Prior x\nSocial")) %>% 
    mutate(name = factor(feature, levels=c("Prior x\nUnsuccessful", "Prior x\nSuccessful", "Prior x\nSocial"))) %>%
    mutate(standardized_weight = ifelse(feature %in% c("Prior x\nUnsuccessful", "Prior x\nSuccessful", "Prior x\nSocial"), -standardized_weight, standardized_weight)) %>% 
    ggplot(aes(x = standardized_weight, y = name, group = as.factor(lake))) + 
    geom_vline(xintercept = 0, linetype = "dashed", size = .4) + 
    # .95
    annotate("rect", xmin = quantile(-lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nUnsuccessful"], .025), xmax =  quantile(-lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nUnsuccessful"], .975), ymin = .65, ymax = 1.35,
             alpha = 1,fill = "#BAE4B3")+
    #.66
    annotate("rect", xmin = quantile(-lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nUnsuccessful"], .17), xmax =  quantile(-lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nUnsuccessful"], .83), ymin = .65, ymax = 1.35,
             alpha = 1,fill = "#31A354")+
    # .5
    annotate("rect", xmin = quantile(-lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nUnsuccessful"], .5)-.0005, xmax =  quantile(-lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nUnsuccessful"], .5)+.0005, ymin = .65, ymax = 1.35,
             alpha = 1,fill = "black")+
    
    annotate("rect", xmin = quantile(-lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nSuccessful"], .025), xmax =  quantile(-lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nSuccessful"], .975), ymin = 1.65, ymax = 2.35,
             alpha = 1,fill = "#BAE4B3")+
    annotate("rect", xmin = quantile(-lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nSuccessful"], .17), xmax =  quantile(-lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nSuccessful"], .83),  ymin = 1.65, ymax = 2.35,
             alpha = 1,fill = "#31A354")+
    annotate("rect", xmin = quantile(-lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nSuccessful"], .5)-.0005, xmax =  quantile(-lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nSuccessful"], .5)+.0005,  ymin = 1.65, ymax = 2.35,
             alpha = 1,fill = "black")+
    stat_pointinterval(position = position_dodge(width = .8), color = "black", alpha = .10, size = 0)  +
    
    annotate("rect", xmin = quantile(-lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nSocial"], .025), xmax =  quantile(-lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nSocial"], .975), ymin = 2.65, ymax = 3.35,
             alpha = 1,fill = "#BAE4B3")+
    annotate("rect", xmin = quantile(-lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nSocial"], .17), xmax =  quantile(-lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nSocial"], .83),  ymin = 2.65, ymax = 3.35,
             alpha = 1,fill = "#31A354")+
    annotate("rect", xmin = quantile(-lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nSocial"], .5)-.0005, xmax =  quantile(-lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nSocial"], .5)+.0005,  ymin = 2.65, ymax = 3.35,
             alpha = 1,fill = "black")+
    stat_pointinterval(position = position_dodge(width = .8), color = "black", alpha = .10, size = 0)  +
    
    
    ylab("Feature")+
    labs(tag = "E", color = "value")+
    ggtitle("Conditional Weights")+
    xlab("Posterior Est. (95%CrI)\nSpot Selection ~ Feature")+
    scale_x_continuous(breaks = seq(from = -2, to = 2, by = .05))+
    scale_y_discrete(labels = c("Unsuccessful x\nPrior", "Successful x\nPrior", "Social x\nPrior"))+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 6),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5, angle = 55, hjust = .5),#, hjust = .5
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.9, .2),
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.background = element_rect(colour = "white", fill = "white"),
          text = element_text(family = "sans"))
  model_results_interaction
  
  ##############################################################################
  lake_effects %>% 
    ggplot(aes(x = feature, y = weights)) + stat_pointinterval()
  model_results_time = lake_effects %>% 
    filter(feature %in% c("Prior x\nTime", "Social x\nTime", "Successful x\nTime", "Unsuccessful x\nTime")) %>% 
    mutate(name = factor(feature, levels=rev(c("Successful x\nTime", "Unsuccessful x\nTime", "Social x\nTime", "Prior x\nTime")))) %>%
    mutate(standardized_weight = ifelse(feature %in% c("Prior x\nTime"), -standardized_weight, standardized_weight)) %>% 
    ggplot(aes(x = standardized_weight, y = name, group = as.factor(lake))) + 
    geom_vline(xintercept = 0, linetype = "dashed", size = .4) + 
    
    # .95
    annotate("rect", xmin = quantile(-lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nTime"], .025), xmax =  quantile(-lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nTime"], .975), ymin = .65, ymax = 1.35,
             alpha = 1,fill = "#BAE4B3")+
    #.66
    annotate("rect", xmin = quantile(-lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nTime"], .17), xmax =  quantile(-lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nTime"], .83), ymin = .65, ymax = 1.35,
             alpha = 1,fill = "#31A354")+
    # .5
    annotate("rect", xmin = quantile(-lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nTime"], .5)-.005, xmax =  quantile(-lake_effects$standardized_main_effect[lake_effects$feature == "Prior x\nTime"], .5)+.005, ymin = .65, ymax = 1.35,
             alpha = 1,fill = "black")+
    
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social x\nTime"], .025), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social x\nTime"], .975), ymin = 1.65, ymax = 2.35,
             alpha = 1,fill = "#BAE4B3")+
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social x\nTime"], .17), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social x\nTime"], .83),  ymin = 1.65, ymax = 2.35,
             alpha = 1,fill = "#31A354")+
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social x\nTime"], .5)-.005, xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social x\nTime"], .5)+.005,  ymin = 1.65, ymax = 2.35,
             alpha = 1,fill = "black")+
    stat_pointinterval(position = position_dodge(width = .8), color = "black", alpha = .10, size = 0)  +
    
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Unsuccessful x\nTime"], .025), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Unsuccessful x\nTime"], .975), ymin = 2.65, ymax = 3.35,
             alpha = 1,fill = "#BAE4B3")+
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Unsuccessful x\nTime"], .17), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Unsuccessful x\nTime"], .83),  ymin = 2.65, ymax = 3.35,
             alpha = 1,fill = "#31A354")+
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Unsuccessful x\nTime"], .5)-.005, xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Unsuccessful x\nTime"], .5)+.005,  ymin = 2.65, ymax = 3.35,
             alpha = 1,fill = "black")+
    
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Successful x\nTime"], .025), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Successful x\nTime"], .975), ymin = 3.65, ymax = 4.35,
             alpha = 1,fill = "#BAE4B3")+
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Successful x\nTime"], .17), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Successful x\nTime"], .83),  ymin = 3.65, ymax = 4.35,
             alpha = 1,fill = "#31A354")+
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Successful x\nTime"], .5)-.005, xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Successful x\nTime"], .5)+.005,  ymin = 3.65, ymax = 4.35,
             alpha = 1,fill = "black")+
    stat_pointinterval(position = position_dodge(width = .8), color = "black", alpha = .10, size = 0)  +
    ylab("Feature")+
    labs(tag = "F", color = "value")+
    ggtitle("Temporal Weights")+
    xlab("Posterior Est. (95%CrI)\nSpot Selection ~ Feature")+
    scale_x_continuous(breaks = seq(from = -2, to = 4, by = 1))+
    #scale_y_discrete(labels = c("Unsuccessful x\nPrior", "Successful x\nPrior", "Social x\nPrior"))+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 6),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5, angle = 55, hjust = .5),#, hjust = .5
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.9, .2),
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.background = element_rect(colour = "white", fill = "white"),
          text = element_text(family = "sans"))
  model_results_time
  
  ################################################################################
  # visualize model comparison
  comparison = readRDS(file = "data/processed_data/model_fit/spot_selection_model/comparison/comparison.rds")
  comparison = as.data.frame(comparison)
  comparison$name = c("+Temporal",
                      "+Conditional\nSocial",
                      "+Conditional\nPersonal",
                      "+Prior",
                      "+Social",
                      "+Personal",
                      "+Locality",
                      "Random")
  comparison$model_id = 1:8
  comparison_spot_selection = comparison %>%
    ggplot(aes(x = elpd_diff, y = fct_reorder(name, model_id))) + 
    geom_col(color = "black", width = .6) +
    geom_errorbarh(aes(xmin = elpd_diff - 4*se_diff, xmax = elpd_diff + 4*se_diff), height = .3, size = .8) + 
    ggtitle("Model Comparison") + 
    xlab(expression(Delta~"ELPD (rel. to best)")) + 
    ylab("Model") + 
    geom_vline(xintercept = 0)+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 5, angle = 55, hjust = .5),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.9, .2),
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.background = element_rect(colour = "white", fill = "white"),
          text = element_text(family = "sans"))
  comparison_spot_selection
  ggsave(comparison_spot_selection, file = "6_plots/comparison_spot_selection.png", 
         width = 6, height = 7, bg = "white", dpi = 600, units = "cm")
  ggsave(comparison_spot_selection, file = "6_plots/comparison_spot_selection.pdf", 
         width = 6, height = 7, bg = "white", dpi = 600, units = "cm")
  ggsave(comparison_spot_selection, file = "6_plots/comparison_spot_selection.svg", 
         width = 6, height = 7, bg = "white", dpi = 600, units = "cm")
  
  
  ##############################################################################
  # visualize confidence results for supplement
  ##############################################################################
  conf_model = readRDS("data/processed_data/model_fit/spot_selection_model/h_confidence_fit.rds")
  conf_model$code()
  table(stan_data$confidence)
  
  table(ecological_knowledge_polygons$confidence)
  plot_s1_b = conf_model %>% 
    spread_draws(weights[parameter]) %>%
    filter(.chain != 4) %>% 
    filter(parameter %in% c(5, 6, 7)) %>% 
    mutate(weights = -weights) %>% 
    mutate(parameter = ifelse(parameter == 5, "Likely", ifelse(parameter == 6, "Very Likely", "Definitely"))) %>% 
    mutate(parameter = factor(parameter, levels = rev(c("Definitely", "Very Likely", "Likely")))) %>% 
    ggplot(aes(y = parameter, x = weights)) + 
    stat_pointinterval(position = position_dodge(width = .8))+
    labs(tag = "A")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    #ggtitle("Confidence in Prior Beliefs")+
    xlab("Posterior Est. (95%CrI)\nSpot Selection ~ Prior | Confidence")+
    scale_x_continuous(breaks = seq(from = -6, to = 6, by = 1))+
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
  contrasts = conf_model %>% 
    spread_draws(weights[parameter]) %>% 
    filter(.chain != 4) %>% 
    mutate(weights = -weights) %>% 
    filter(parameter %in% c(5, 6, 7)) %>% 
    mutate(parameter = ifelse(parameter == 5, "Likely", ifelse(parameter == 6, "Very Likely", "Definitely"))) %>% 
    pivot_wider(names_from = parameter, values_from = weights) %>% 
    mutate(`Definitely - Very Likely` = Definitely - `Very Likely`,
           `Definitely - Likely` = Definitely - Likely,
           `Very Likely - Likely` = `Very Likely` - Likely)
  
  plot_s1_c = contrasts %>% 
    pivot_longer(cols = 7:9) %>% 
    mutate(name = factor(name, levels = rev(c("Definitely - Very Likely", "Definitely - Likely", "Very Likely - Likely")))) %>% 
    ggplot(aes(y = name, x = value)) + 
    stat_pointinterval() + 
    geom_vline(xintercept = 0, linetype = "dashed")+
    xlab("Posterior Est. (95%CrI)\nContrasts Spot Selection ~ Prior | Confidence")+
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
  ggsave(panel_conf, file = "6_plots/plot_confidence_spatial.png", 
         width = 12, height = 6, bg = "white", dpi = 300, units = "cm")
  ggsave(panel_conf, file = "6_plots/plot_confidence_spatial.pdf", 
         width = 12, height = 6, bg = "white", dpi = 300, units = "cm")
  ggsave(panel_conf, file = "6_plots/plot_confidence_spatial.svg", 
         width = 12, height = 6, bg = "white", dpi = 300, units = "cm", device = svglite::svglite)
  
  ##############################################################################
  # descriptive results of prior use
  ##############################################################################
  # step 1: run random trip rotations
  source("3_spot_selection_model/7_simulate_alternative_trip.R")
  sim_trip()
  
  # extract prior distance for random trips
  source("3_spot_selection_model/8_extract_prior_dist.R")
  extract_prior_dist()
  
  # use model without temporal effects for "marginal" trip estimates (average across time)
  spot_selection_model_fit = readRDS("data/processed_data/model_fit/spot_selection_model/g_cond_social_fit.rds")
  trip_estimates = spot_selection_model_fit %>% 
    spread_draws(v_trip[trip, parameter]) %>% 
    group_by(trip, parameter) %>% 
    summarize(est = mean(v_trip)) %>% 
    pivot_wider(names_from = parameter, names_prefix = "parameter_", values_from = "est")
  
  # add lake level effects and id level effects
  id_estimates = spot_selection_model_fit %>% 
    spread_draws(v_id[id, parameter]) %>% 
    group_by(id, parameter) %>% 
    summarize(est = mean(v_id)) %>% 
    pivot_wider(names_from = parameter, names_prefix = "parameter_", values_from = "est")
  lake_estimates = spot_selection_model_fit %>% 
    spread_draws(v_lake[lake, parameter]) %>% 
    group_by(lake, parameter) %>% 
    summarize(est = mean(v_lake)) %>% 
    pivot_wider(names_from = parameter, names_prefix = "parameter_", values_from = "est")
  
  # join trip estimates to spot data
  identifiers_spot_selection$day = as.integer(identifiers_spot_selection$day)
  identifiers_spot_selection$year = as.integer(identifiers_spot_selection$year)
  spot_selection_data = left_join(spot_selection_data, identifiers_spot_selection %>% group_by(day, camera_id, year) %>% slice(1) %>% dplyr::select(trip, lake, id, day, year, camera_id))
  spot_selection_data = left_join(spot_selection_data, trip_estimates %>% dplyr::select(trip, prior_use_trip = parameter_6, loss_attraction_trip = parameter_4))
  spot_selection_data = left_join(spot_selection_data, lake_estimates %>% dplyr::select(lake, prior_use_lake = parameter_6, loss_attraction_lake = parameter_4))
  spot_selection_data = left_join(spot_selection_data, id_estimates %>% dplyr::select(id, prior_use_id = parameter_6, loss_attraction_id = parameter_4))
  spot_selection_data$prior_use = spot_selection_data$prior_use_id + spot_selection_data$prior_use_lake + spot_selection_data$prior_use_trip
  spot_selection_data$loss_attraction = spot_selection_data$loss_attraction_id + spot_selection_data$loss_attraction_lake + spot_selection_data$loss_attraction_trip
  spot_selection_data$prior_use = -spot_selection_data$prior_use
  
  # compute explored area for each trip
  areas_data = spot_selection_data %>% 
    group_by(day, year, trip, id) %>% 
    mutate(n_spots = n()) %>% 
    filter(n_spots > 2) %>% 
    reframe(area = (cbind(x[chull(x,y)], y[chull(x, y)])),
              prior_use = mean(prior_use),
              loss_attraction = mean(loss_attraction),
            n = n_spots[1]) %>% 
    group_by(day, year, trip, id) %>% 
    mutate(n = n()) %>% 
    filter(n >= 4) %>% 
    summarize(area = Polygon(area)@area,
              prior_use = mean(prior_use),
              loss_attraction = mean(loss_attraction)) %>% 
    mutate(area = area / (1000 * 1000))
  
  # compute lake, and trip estimates for each explored area entry (scaled)
  areas_data$lake = paste(areas_data$day, areas_data$year)
  
  # model based analyses: predict explored area from prior use, loss attraction sequentially
  model_area = brms::brm(data = areas_data, 
                         formula = brmsformula(area ~ 1 + prior_use + loss_attraction + (1 + prior_use + loss_attraction | lake) + (1 + prior_use + loss_attraction | id)
                                               , sigma ~ 1),
                         family = lognormal(),
                         chains = 4,
                         cores = 32,
                         backend = "cmdstanr",
                         threads = threading(8),
                         seed = stan_seed)
  
  model_area_1 = brms::brm(data = areas_data, 
                           formula = brmsformula(area ~ 1 + prior_use + (1 + prior_use | lake) + (1 + prior_use | id)
                                                , sigma ~ 1),
                           family = lognormal(),
                           chains = 4,
                           cores = 32,
                           backend = "cmdstanr",
                           threads = threading(8),
                           seed = stan_seed)
  
  model_area_2 = brms::brm(data = areas_data, 
                           formula = brmsformula(area ~ 1 + loss_attraction + (1 + loss_attraction | lake) + (1 + loss_attraction | id)
                                                 , sigma ~ 1),
                           family = lognormal(),
                           chains = 4,
                           cores = 32,
                           backend = "cmdstanr",
                           threads = threading(8),
                           seed = stan_seed)
  # save model fits
  saveRDS(model_area, file = "data/processed_data/area_model_1.rds")
  saveRDS(model_area_1, file = "data/processed_data/area_model_2.rds")
  saveRDS(model_area_2, file = "data/processed_data/area_model_3.rds")
  
  # Plot posterior effects
  posterior_predictions = rstantools::posterior_epred(model_area_1, newdata = data.frame(prior_use = seq(from = range(areas_data$prior_use)[1], to = range(areas_data$prior_use)[2], by = .05)), re_formula = NA)
  posterior_predictions = as.data.frame(posterior_predictions)
  
  posterior_predictions$draw = 1:4000
  plot_posterior = posterior_predictions %>% 
    pivot_longer(1:(ncol(posterior_predictions)-1)) %>% 
    group_by(draw) %>% 
    mutate(prior_use = seq(from = range(areas_data$prior_use)[1], to = range(areas_data$prior_use)[2], by = .05)) %>% 
    ggplot(aes(x = prior_use, y = value)) + 
    coord_cartesian(xlim = c(range(areas_data$prior_use)[1], range(areas_data$prior_use)[2]), ylim = c(0, 2))+
    stat_lineribbon(.width = .95, fill = "darkgreen", alpha = .3) + 
    xlab("Prior Belief Use") + 
    labs(tag = "A") +
    geom_point(data = areas_data, aes(x = prior_use, y = area), alpha = .4, size = 1, color = "black")+
    ylab("Explored Area [sqkm]") +
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.9, .2),
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.background = element_rect(colour = "white", fill = "white"),
          text = element_text(family = "sans"))
  plot_posterior
  
  ##################
  posterior_predictions2 = posterior_epred(model_area_2, newdata = data.frame(loss_attraction = seq(from = range(areas_data$loss_attraction)[1], to = range(areas_data$loss_attraction)[2], by = .05)), re_formula = NA)
  posterior_predictions2 = as.data.frame(posterior_predictions2)
  
  posterior_predictions2$draw = 1:4000
  plot_posterior2 = posterior_predictions2 %>% 
    pivot_longer(1:(ncol(posterior_predictions2)-1)) %>% 
    group_by(draw) %>% 
    mutate(loss_attraction = seq(from = range(areas_data$loss_attraction)[1], to = range(areas_data$loss_attraction)[2], by = .05)) %>% 
    ggplot(aes(x = -loss_attraction, y = value)) + 
    coord_cartesian(xlim = c(range(-areas_data$loss_attraction)[1], range(-areas_data$loss_attraction)[2]), ylim = c(0, 2))+
    stat_lineribbon(.width = .95, fill = "orange2", alpha = .3) + 
    ylab("Explored Area [sqkm]") + 
    labs(tag = "B")+
    geom_point(data = areas_data, aes(x = -loss_attraction, y = area), alpha = .4, size = 1, color = "black")+
    #stat_summary_bin(data = areas_data, aes(x = -loss_attraction, y = area), alpha = .4, fun.args = list(mult = 2), fun.data = mean_se)+
    xlab("Avoidance of Unsuccessful Spots") +
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.9, .2),
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.background = element_rect(colour = "white", fill = "white"),
          text = element_text(family = "sans"))
  plot_posterior2
  
  panel_area = gridExtra::arrangeGrob(
    grobs = list(plot_posterior,
                 plot_posterior2),
    widths = c(1, 1),
    layout_matrix = rbind(c(1,2))
  )
  
  ggsave(panel_area, file = "6_plots/panel_area.png", 
         width = 12, height = 6, bg = "white", dpi = 600, units = "cm")
  ggsave(panel_area, file = "6_plots/panel_area.pdf", 
         width = 12, height = 6, bg = "white", dpi = 600, units = "cm")
  ggsave(panel_area, file = "6_plots/panel_area.svg", 
         width = 12, height = 6, bg = "white", dpi = 600, units = "cm", device = svglite::svglite)
  
  ################################################################################
  # predict loss attraction from prior use
  parameter_regression_1 = brm(parameter_3 ~ parameter_6,
                               data = trip_estimates, 
                               family = "gaussian")
  parameter_regression_2 = brm(parameter_4 ~ parameter_6,
                               data = trip_estimates, 
                               family = "gaussian")
  
  posterior_predictions_1 = rstantools::posterior_epred(parameter_regression_1, newdata = data.frame(parameter_6 = seq(from = range(trip_estimates$parameter_6)[1], to = range(trip_estimates$parameter_6)[2], by = .05)), re_formula =NA)
  posterior_predictions_1 = as.data.frame(posterior_predictions_1)
  posterior_predictions_1$draw = 1:4000
  posterior_predictions_1 = posterior_predictions_1 %>% 
    pivot_longer(1:(ncol(posterior_predictions_1)-1)) %>% 
    group_by(draw) %>% 
    mutate(prior_use = seq(from = range(trip_estimates$parameter_6)[1], to = range(trip_estimates$parameter_6)[2], by = 0.05))
  
  posterior_predictions_2 = rstantools::posterior_epred(parameter_regression_2, newdata = data.frame(parameter_6 = seq(from = range(trip_estimates$parameter_6)[1], to = range(trip_estimates$parameter_6)[2], by = .05)), re_formula =NA)
  posterior_predictions_2 = as.data.frame(posterior_predictions_2)
  posterior_predictions_2$draw = 1:4000
  posterior_predictions_2 = posterior_predictions_2 %>% 
    pivot_longer(1:(ncol(posterior_predictions_2)-1)) %>% 
    group_by(draw) %>% 
    mutate(prior_use = seq(from = range(trip_estimates$parameter_6)[1], to = range(trip_estimates$parameter_6)[2], by = 0.05))
  
  # select spots for visualization
  trips = spot_selection_data %>% dplyr::select(trip)
  # parameter 6 = prior (higher values mean less attraction to prior)
  # parameter 4 = loss attraction (higher values mean less avoidance of losses)
  # C (strong loss avoidance and low prior use)
  # D (low loss avoidance and high prior use)
  quantile(-trip_estimates$parameter_6, probs = c(0.10, 0.90))
  quantile(-trip_estimates$parameter_4, probs = c(0.10, 0.90))
  quantile(trip_estimates$parameter_3, probs = c(0.10, 0.90))
  
  trips_c = trip_estimates %>% filter(trip %in% pull(trips)) %>% filter(-parameter_6 < 0.6 & -parameter_4 > 0.6) %>% dplyr::select(trip)
  trips_c = spot_selection_data %>% filter(trip %in% pull(trips_c)) %>% group_by(trip) %>% summarize(r = mean(Return),
                                                                                                     day = day[1],
                                                                                                     year = year[1],
                                                                                                     participant_id = participant_id[1],
                                                                                                     camera_id = camera_id[1],
                                                                                                     n_spots = n())
  trips_d = trip_estimates %>% filter(trip %in% pull(trips)) %>% filter(-parameter_6 > 0.6 & -parameter_4 < -0.6) %>% dplyr::select(trip)
  trips_d = spot_selection_data %>% filter(trip %in% pull(trips_d)) %>% group_by(trip) %>% summarize(r = mean(Return),
                                                                                                     day = day[1],
                                                                                                     year = year[1],
                                                                                                     participant_id = participant_id[1],
                                                                                                     camera_id = camera_id[1],
                                                                                                     n_spots = n())
  
  # trips 86 and 65 come from the same day and roughly match in success rate and number of spots
  trip_estimate_plot1 = trip_estimates %>%
    ggplot(aes(y = parameter_3, x = -parameter_6)) + 
    geom_point(size = .7, alpha = .6, color = "lightgrey") + 
    stat_lineribbon(data = posterior_predictions_1, aes(x = -prior_use, y = value), .width= .95, fill = "darkgreen", alpha = .3) + 
    coord_cartesian(xlim = c(range(-trip_estimates$parameter_6)[1], range(-trip_estimates$parameter_6)[2]), ylim = c(range(trip_estimates$parameter_3)[1], range(trip_estimates$parameter_6)[2]))+
    scale_y_continuous(breaks = seq(from = -3, to = 3, by = 1)) +
    scale_x_continuous(breaks = seq(from = -3, to = 3, by = 1)) +
    geom_vline(xintercept = 0, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    xlab("Prior Belief Use") + 
    labs(tag = "A")+
    geom_point(aes(y = -trip_estimates$parameter_3[trip_estimates$trip == 86], x = -trip_estimates$parameter_6[trip_estimates$trip == 86]),  alpha = .6, color = "red", size = 1.5) + 
    geom_point(aes(y = -trip_estimates$parameter_3[trip_estimates$trip == 65], x = -trip_estimates$parameter_6[trip_estimates$trip == 65]),  alpha = .6, color = "red", size = 1.5) + 
    annotate("text", y = -trip_estimates$parameter_3[trip_estimates$trip == 86], x = -trip_estimates$parameter_6[trip_estimates$trip == 86]-0.2, label = "C")+
    annotate("text", y = -trip_estimates$parameter_3[trip_estimates$trip == 65], x = -trip_estimates$parameter_6[trip_estimates$trip == 65]-0.2, label = "D")+
    ylab("Attraction to Successful Spots") + 
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.9, .2),
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.background = element_rect(colour = "white", fill = "white"),
          text = element_text(family = "sans"))
  trip_estimate_plot1
  
  trip_estimate_plot2 = trip_estimates %>%
    ggplot(aes(y = -parameter_4, x = -parameter_6)) + 
    geom_point(size = .7, alpha = .6, color = "lightgrey") + 
    stat_lineribbon(data = posterior_predictions_2, aes(x = -prior_use, y = -value), .width= .95, fill = "darkgreen", alpha = .3) + 
    coord_cartesian(xlim = c(range(-trip_estimates$parameter_6)[1], range(-trip_estimates$parameter_6)[2]), ylim = c(range(-trip_estimates$parameter_4)[1], range(-trip_estimates$parameter_4)[2]))+
    geom_vline(xintercept = 0, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    scale_y_continuous(breaks = seq(from = -3, to = 3, by = 1)) +
    scale_x_continuous(breaks = seq(from = -3, to = 3, by = 1)) +
    xlab("Prior Belief Use") + 
    labs(tag = "B")+
    geom_point(aes(y = -trip_estimates$parameter_4[trip_estimates$trip == 86], x = -trip_estimates$parameter_6[trip_estimates$trip == 86]),  alpha = .6, color = "red", size = 1.5) + 
    geom_point(aes(y = -trip_estimates$parameter_4[trip_estimates$trip == 65], x = -trip_estimates$parameter_6[trip_estimates$trip == 65]),  alpha = .6, color = "red", size = 1.5) + 
    annotate("text", y = -trip_estimates$parameter_4[trip_estimates$trip == 86], x = -trip_estimates$parameter_6[trip_estimates$trip == 86] - .2, label = "C")+
    annotate("text", y = -trip_estimates$parameter_4[trip_estimates$trip == 65], x = -trip_estimates$parameter_6[trip_estimates$trip == 65] - .2, label = "D")+
    ylab("Avoidance of Unsuccessful Spots") + 
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.9, .2),
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.background = element_rect(colour = "white", fill = "white"),
          text = element_text(family = "sans"))
  trip_estimate_plot2
  
  # Step 2: Load competition areas
  competition_areas = readRDS("data/study_lakes.rds")
  depth_profiles = readRDS("data/depth_profiles.rds")
  profile_transformed = st_transform(depth_profiles[[3]], "epsg:3067")
  trips_c
  ex1=ggplot() + 
    geom_sf(data = competition_areas[[3]], alpha = .5) + 
    geom_sf(data = st_crop(profile_transformed, st_buffer(competition_areas[[3]], 100)), color = "lightgrey", alpha = .8) + 
    geom_sf(data = ecological_knowledge %>% filter(camera_id == 48 & day == 3), fill = "lightblue") + 
    geom_path(data = spot_selection_data %>% filter(trip == 86), aes(x = x, y = y), color = "darkgrey") + 
    geom_point(data = spot_selection_data %>% filter(trip == 86) %>% filter(Return == 0), aes(x = x, y = y), color = "lightgreen", size = 1, alpha = .5) + 
    geom_point(data = spot_selection_data %>% filter(trip == 86) %>% filter(Return == 1), aes(x = x, y = y), color = "darkgreen", size = 1) + 
    geom_point(data = spot_selection_data %>% filter(trip == 86) %>% filter(spot_id == 1), aes(x = x, y = y), color = "black", shape = 4, size = 3) + 
    xlab("Longitude") + 
    ylab("Latitude") + 
    labs(tag = "C") +
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.9, .2),
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.background = element_rect(colour = "white", fill = "white"),
          text = element_text(family = "sans"))
  ex1
  trips_d
  ex2=ggplot() + 
    geom_sf(data = competition_areas[[3]], alpha = .5) + 
    geom_sf(data = st_crop(profile_transformed, st_buffer(competition_areas[[3]], 100)), color = "lightgrey", alpha = .8) + 
    geom_sf(data = ecological_knowledge %>% filter(camera_id == 26 & day == 3), fill = "lightblue") + 
    geom_path(data = spot_selection_data %>% filter(trip == 65), aes(x = x, y = y), color = "darkgrey") + 
    geom_point(data = spot_selection_data %>% filter(trip == 65) %>% filter(Return == 0), aes(x = x, y = y), color = "lightgreen", size = 1, alpha = .5) + 
    geom_point(data = spot_selection_data %>% filter(trip == 65) %>% filter(Return == 1), aes(x = x, y = y), color = "darkgreen", size = 1) + 
    geom_point(data = spot_selection_data %>% filter(trip == 65) %>% filter(spot_id == 1), aes(x = x, y = y), color = "black", shape = 4, size = 3) + 
    xlab("Longitude") + 
    ylab("Latitude") + 
    labs(tag = "D") +
    annotate("segment", x = 652605.5, xend = 652605.5 + 900, y = 6917116, yend = 6917116 + 20, colour = "black", linewidth = .1) +
    annotate("segment", x = 652605.5, xend = 652605.5 + 260, y = 6917116, yend = 6917116 + 650, colour = "black", linewidth = .1) +
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.9, .2),
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.background = element_rect(colour = "white", fill = "white"),
          text = element_text(family = "sans"))
  ex2
  
  library(cowplot)
  ex4 = ggplot() + 
    geom_sf(data = competition_areas[[3]], alpha = .5,fill = "lightgrey") + 
    geom_sf(data = st_crop(profile_transformed, st_buffer(competition_areas[[3]], 100)), color = "lightgrey", alpha = .8) + 
    geom_sf(data = ecological_knowledge %>% filter(camera_id == 26 & day == 3), fill = "lightblue") + 
    coord_sf(xlim = c(min((spot_selection_data %>% filter(trip == 65))$x - 5), max((spot_selection_data %>% filter(trip == 65))$x) + 5),
             ylim = c(min((spot_selection_data %>% filter(trip == 65))$y - 5), max((spot_selection_data %>% filter(trip == 65))$y) + 5))+
    geom_path(data = spot_selection_data %>% filter(trip == 65), aes(x = x, y = y), color = "darkgrey") + 
    geom_point(data = spot_selection_data %>% filter(trip == 65) %>% filter(Return == 0), aes(x = x, y = y), color = "lightgreen", size = 1, alpha = 1) + 
    geom_point(data = spot_selection_data %>% filter(trip == 65) %>% filter(Return == 1), aes(x = x, y = y), color = "darkgreen", size = 1) + 
    geom_point(data = spot_selection_data %>% filter(trip == 65) %>% filter(spot_id == 1), aes(x = x, y = y), color = "black", shape = 4, size = 3) + 
    xlab("Longitude") + 
    ylab("Latitude") + 
    labs(tag = "") +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      plot.margin = margin(0, 0, 0, 0, unit = "pt")
    )
  ex4
  plot.with.inset <-
    ggdraw() +
    draw_plot(ex2) +
    draw_plot(ex4, x = 0.5625, y = .625, width = .35, height = .35)
  #ggsave(plot.with.inset, file = "6_plots/panel_illustration_test.png", 
  #       width = 12, height = 12, bg = "white", dpi = 300, units = "cm")
  
  # make panel
  # panel
  panel_illustrations = gridExtra::arrangeGrob(
    grobs = list(trip_estimate_plot1,
                 trip_estimate_plot2,
                 ex1,
                 plot.with.inset),
    widths = c(1, 1),
    layout_matrix = rbind(c(1, 2),
                          c(3, 4))
  )
  
  ggsave(panel_illustrations, file = "6_plots/panel_illustration.png", 
         width = 12, height = 12, bg = "white", dpi = 600, units = "cm")
  ggsave(panel_illustrations, file = "6_plots/panel_illustration.pdf", 
         width = 12, height = 12, bg = "white", dpi = 600, units = "cm")
  ggsave(panel_illustrations, file = "6_plots/panel_illustration.svg", 
         width = 12, height = 12, bg = "white", dpi = 600, units = "cm", device = svglite::svglite)
  
  # posterior correlation
  spot_selection_model_fit %>% 
    spread_draws(Rho_trip[parameter_1, parameter_2]) %>% 
    filter(parameter_2 == 3 & parameter_1 == 6) %>% 
    reframe(quantile(Rho_trip, probs = c(0.025, .975)),
              mean = mean(Rho_trip))
  spot_selection_model_fit %>% 
    spread_draws(Rho_trip[parameter_1, parameter_2]) %>% 
    filter(parameter_2 == 4 & parameter_1 == 6) %>% 
    reframe(quantile(Rho_trip, probs = c(0.025, .975)),
              mean = mean(Rho_trip))
  
  ##############################################################################
  # additional plots from simulations
  samples = unique(pull(spot_selection_data %>% filter(day == 3) %>% dplyr::select(trip)))
  random_trips = sample(samples, 16, replace = F)
  random_trips_sorted = trip_estimates[trip_estimates$trip %in% random_trips, c("trip", "parameter_6", "parameter_3")]
  random_trips_sorted = random_trips_sorted %>% arrange(parameter_6)
  random_ex = list()
  for (i in 1:nrow(random_trips_sorted)){
    random_trip = pull(random_trips_sorted[i,"trip"])
    random_day = pull(spot_selection_data %>% filter(trip == random_trip) %>% dplyr::select(day) %>% slice(1))
    random_camera_id = pull(spot_selection_data %>% filter(trip == random_trip) %>% dplyr::select(camera_id) %>% slice(1))
    
    # get parameters
    prior_use_random = -pull(random_trips_sorted[i, "parameter_6"])
    loss_avoidance_random = -pull(random_trips_sorted[i, "parameter_3"])
    
    # Step 2: Load competition areas
    random_profile_transformed = st_transform(depth_profiles[[random_day]], "epsg:3067")
    
    random_ex[[i]]=ggplot() + 
      geom_sf(data = competition_areas[[random_day]], alpha = .5) + 
      geom_sf(data = st_crop(profile_transformed, st_buffer(competition_areas[[random_day]], 100)), color = "lightgrey", alpha = .8) + 
      geom_sf(data = ecological_knowledge %>% filter(camera_id == random_camera_id & day == random_day), fill = "lightblue") + 
      geom_path(data = spot_selection_data %>% filter(trip == random_trip), aes(x = x, y = y), color = "darkgrey") + 
      geom_point(data = spot_selection_data %>% filter(trip == random_trip) %>% filter(Return == 0), aes(x = x, y = y), color = "lightgreen", size = 1, alpha = .5) + 
      geom_point(data = spot_selection_data %>% filter(trip == random_trip) %>% filter(Return == 1), aes(x = x, y = y), color = "darkgreen", size = 1) + 
      geom_point(data = spot_selection_data %>% filter(trip == random_trip) %>% filter(spot_id == 1), aes(x = x, y = y), color = "black", shape = 4, size = 3) + 
      ggtitle(paste("Loss Avoidance = ", round(loss_avoidance_random, 2), "\nPrior Attraction = ", round(prior_use_random, 2), sep = "")) + 
      xlab("Longitude") + 
      ylab("Latitude") + 
      labs(tag = LETTERS[i]) +
      theme(panel.background = element_rect(fill = NA),
            axis.title.x = element_text(size = 7),
            axis.title.y = element_text(size = 7),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA),
            legend.position = c(.9, .2),
            plot.tag = element_text(size = 7, face = "bold"),
            title = element_text(size = 5),
            legend.text = element_text(size = 7),
            legend.title = element_text(size = 7),
            legend.background = element_rect(colour = "white", fill = "white"),
            text = element_text(family = "sans"))
  }
  
  panel_random_ex = gridExtra::arrangeGrob(
    grobs = random_ex,
    widths = c(1, 1, 1, 1),
    layout_matrix = matrix(1:16, nrow = 4, ncol = 4)
  )
  
  ggsave(panel_random_ex, file = "6_plots/panel_random_ex.png", 
         width = 16, height = 18, bg = "white", dpi = 600, units = "cm")
  ggsave(panel_random_ex, file = "6_plots/panel_random_ex.pdf", 
         width = 16, height = 18, bg = "white", dpi = 600, units = "cm")
  ggsave(panel_random_ex, file = "6_plots/panel_random_ex.svg", 
         width = 16, height = 18, bg = "white", dpi = 600, units = "cm", device = svglite::svglite)
  
  ################################################################################
  # add descriptive prior use to trip estimate data
  Feature_matrix_trip_selection = readRDS("data/processed_data/feature_matrix_trip_selection.rds")
  identifiers_trip_selection = readRDS("data/processed_data/identifiers_trip_selection.rds")
  ecological_knowledge_distances_trip_selection = readRDS(file = "data/processed_data/knowledge_distance_matrix_trip_selection.rds")
   
  # compute prior use for trips
  res = matrix(nrow = 251, ncol = 1)
  eco = vector()
  
  for (i in 1:length(Feature_matrix_trip_selection)){
    
    if (is.null(dim(ecological_knowledge_distances_trip_selection[[i]]))){
      eco[i] = min(ecological_knowledge_distances_trip_selection[[i]][1])
    } else {
      eco[i] = min(ecological_knowledge_distances_trip_selection[[i]][1,])
    }
  }
  
  identifiers_trip_selection$eco_knowledge = eco
  identifiers_trip_selection = identifiers_trip_selection %>% filter(spot_id != 1)
  
  # now add trip id from model to identifiers
  identifiers_trip_selection$day = as.integer(identifiers_trip_selection$day)
  identifiers_trip_selection = left_join(identifiers_trip_selection, identifiers_spot_selection %>% dplyr::select(day, year, spot_id, camera_id, trip))
  
  prior_use_trips = identifiers_trip_selection %>% 
    mutate(use = ifelse(eco_knowledge <= 100, 1, 0)) %>% 
    group_by(day, year, camera_id, trip) %>% 
    summarize(use = mean(use))
  
  # join to trip estimates for later
  trip_estimates = left_join(trip_estimates, prior_use_trips)
  trip_estimates = left_join(trip_estimates, spot_selection_data %>% group_by(trip) %>% slice(1)) # a bit confusing: parameter_6 is the trip estimate, prior_use is the marginal estimate, use is the descriptive estimate.
  
  # predict prior use from trip estimates
  validation_model = brms::brm(bf(use ~ prior_use,
                                  phi ~ prior_use,
                                  zoi ~ prior_use,
                                  coi ~ prior_use),
                               family = zero_one_inflated_beta(),
                               data = trip_estimates, 
                               cores = 4,
                               chains = 4)
  saveRDS(validation_model, file = "data/processed_data/validation_model_results.rds")
  
  posterior_predictions_zoi = posterior_epred(validation_model, newdata = data.frame(prior_use = seq(from = range(trip_estimates$prior_use)[1], to = range(trip_estimates$prior_use)[2], by = .01)))
  posterior_predictions_zoi = as.data.frame(posterior_predictions_zoi)
  posterior_predictions_zoi$draw = 1:4000
  plot_posterior_zoi = posterior_predictions_zoi %>% 
    pivot_longer(1:(ncol(posterior_predictions_zoi)-1)) %>% 
    group_by(draw) %>% 
    mutate(prior_attraction = seq(from = range(trip_estimates$prior_use)[1], to = range(trip_estimates$prior_use)[2], by = .01)) %>% 
    ggplot(aes(x = prior_attraction, y = value)) + 
    stat_lineribbon(.width = .95, fill = "darkolivegreen", alpha = .3) + 
    ylab("Proportion of spots within prior area") + 
    geom_point(data = trip_estimates, aes(x = prior_use, y = use), alpha = .4, size = 1, color = "black")+
    #stat_summary_bin(data = areas_data, aes(x = -loss_attraction, y = area), alpha = .4, fun.args = list(mult = 2), fun.data = mean_se)+
    xlab("Model estimate: Prior Belief Use") +
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.9, .2),
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.background = element_rect(colour = "white", fill = "white"),
          text = element_text(family = "sans"))
  plot_posterior_zoi
  ggsave(plot_posterior_zoi, file = "6_plots/model_validation.png", 
         width = 6, height = 6, bg = "white", dpi = 600, units = "cm")
  ggsave(plot_posterior_zoi, file = "6_plots/model_validation.pdf", 
         width = 6, height = 6, bg = "white", dpi = 600, units = "cm")
  ggsave(plot_posterior_zoi, file = "6_plots/model_validation.svg", 
         width = 6, height = 6, bg = "white", dpi = 600, units = "cm", device = svglite::svglite)
  
  ################################################################################
  # descriptive measure of prior use for plotting
  # load alternative trips
  Feature_matrix = readRDS("data/processed_data/feature_matrix_trip_selection.rds")
  
  # compute distance to prior area
  knowledge_distance_matrix = readRDS(file = "data/processed_data/knowledge_distance_matrix_trip_selection.rds")
  
  # confidence
  confidence_matrix = readRDS(file = "data/processed_data/confidence_matrix_trip_selection.rds")
  
  # compute prior distance for each location
  prior_use = list()
  conf = list()
  for (i in 1:length(knowledge_distance_matrix)){
    if (is.null(dim(knowledge_distance_matrix[[i]]))){
      prior_use[[i]] = rep(NA, length(knowledge_distance_matrix[[i]]))
      conf[[i]] = rep(NA, length(knowledge_distance_matrix[[i]]))
    } else {
      prior_use[[i]] = apply(knowledge_distance_matrix[[i]], 1, min)
      conf[[i]] = confidence_matrix[,i]
    }
  }
  rm(c)
  prior_use = do.call(c, prior_use)
  conf = do.call(c, conf)
  
  Feature_matrix = Feature_matrix %>% dplyr::bind_rows()
  Feature_matrix$prior_distance = prior_use
  Feature_matrix$confidence = conf
  
  # add camera id
  Feature_matrix = Feature_matrix %>% 
    group_by(unique_trip_id) %>% 
    mutate(camera_id = ifelse(is.na(camera_id), camera_id[1], camera_id))
  
  # make distribution plot: prior use
  Feature_matrix$prior_use = ifelse(Feature_matrix$prior_distance <= 100, 1, 0)
  descr_plot = Feature_matrix %>% 
    group_by(unique_trip_id, iteration) %>% 
    summarize(prior_use = mean(prior_use)) %>% 
    ungroup() %>% 
    mutate(quantile_prior_use = cut(prior_use, breaks = seq(from = 0, to = 1, by = .1), include.lowest = T, right = T)) %>% 
    group_by(quantile_prior_use, iteration) %>% 
    summarize(n = n())
  
  descr_plot_confidence = Feature_matrix %>% 
    filter(confidence > -99) %>% 
    group_by(unique_trip_id, iteration, confidence) %>% 
    summarize(prior_use = mean(prior_use)) %>% 
    ungroup() %>% 
    mutate(quantile_prior_use = cut(prior_use, breaks = seq(from = 0, to = 1, by = .1), include.lowest = T, right = T)) %>% 
    group_by(quantile_prior_use, iteration, confidence) %>% 
    summarize(n = n())
  
  colors = c("Random" = "darkolivegreen","Data" = "black")
  
  hists = ggplot() + 
    stat_interval(data = descr_plot %>% filter(iteration > 0), aes(x = as.numeric(quantile_prior_use)/10 - .05, y = n, color = "Random"), .width = c(.66,.95), alpha = .7) + 
    stat_interval(data = descr_plot %>% filter(iteration > 0), aes(x = as.numeric(quantile_prior_use)/10 - .05, y = n), color = "darkgreen", .width = c(.2), alpha = .7) + 
    #stat_summary(data = descr_plot %>% filter(iteration > 0), aes(x = as.numeric(quantile_prior_use), y = n), geom = "line", color = "green3") + 
    geom_point(data = descr_plot %>% filter(iteration == 0), aes(x = as.numeric(quantile_prior_use)/10 - .05, y = n, color = "Data"), shape = 15) + 
    #stat_summary(data = descr_plot %>% filter(iteration == 0), aes(x = as.numeric(quantile_prior_use), y = n), geom = "line", color = "black") +
    labs(tag = "C") +
    annotate("text", x = .8, y = 175, label = expression(epsilon == 100), size = 7/.pt)+
    xlab("Proportion of spots within prior area")+
    ylab("N trips")+
    scale_x_continuous(breaks = seq(from = 0, to = 1, by = .1))+
    scale_color_manual(values = colors)+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 6),
          axis.title.y = element_text(size = 6),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),#, hjust = .5
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = "none",
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 5, margin=margin(0,0,0,0)),
          legend.text = element_text(size = 5),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white", linewidth = .2),
          text = element_text(family = "sans"),
          legend.box.margin = margin(-.2, -.2, -.2, -.2),
          legend.key.height = unit(.01, "cm"),
          plot.margin=unit(c(.1, .2, .1, .2), "cm"))
  hists
  descr_plot %>% 
    filter(iteration == 0) %>% 
    ungroup() %>% 
    filter(quantile_prior_use == "(0.9,1]") %>% 
    summarize(n/length(unique(Feature_matrix$unique_trip_id)))
  descr_plot %>% 
    filter(iteration > 0) %>% 
    ungroup() %>% 
    filter(quantile_prior_use == "(0.9,1]") %>% 
    summarize(n = mean(n)) %>% 
    summarize(n/length(unique(Feature_matrix$unique_trip_id)))
  
  # prior use as function of confidence
  hists_confidence_1 = ggplot() + 
    stat_interval(data = descr_plot_confidence %>% filter(iteration > 0) %>% filter(confidence == 1), aes(x = as.numeric(quantile_prior_use)/10 - .05, y = n, color = "Random"), .width = c(.5,.95), alpha = .7) + 
    stat_interval(data = descr_plot_confidence %>% filter(iteration > 0) %>% filter(confidence == 1), aes(x = as.numeric(quantile_prior_use)/10 - .05, y = n), color = "darkgreen", .width = c(.2), alpha = .7) + 
    #stat_summary(data = descr_plot %>% filter(iteration > 0), aes(x = as.numeric(quantile_prior_use), y = n), geom = "line", color = "green3") + 
    geom_point(data = descr_plot_confidence %>% filter(iteration == 0) %>% filter(confidence == 1), aes(x = as.numeric(quantile_prior_use)/10 - .05, y = n, color = "Data"), shape = 15) + 
    #stat_summary(data = descr_plot %>% filter(iteration == 0), aes(x = as.numeric(quantile_prior_use), y = n), geom = "line", color = "black") +
    labs(tag = "A") +
    annotate("text", x = .8, y = 175, label = expression(epsilon == 100), size = 7/.pt)+
    xlab("Proportion of spots within prior area")+
    ylab("N trips")+
    ggtitle("Confidence: Likely")+
    scale_x_continuous(breaks = seq(from = 0, to = 1, by = .1))+
    scale_color_manual(values = colors)+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 6),
          axis.title.y = element_text(size = 6),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),#, hjust = .5
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = "none",
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 5, margin=margin(0,0,0,0)),
          legend.text = element_text(size = 5),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white", linewidth = .2),
          text = element_text(family = "sans"),
          legend.box.margin = margin(-.2, -.2, -.2, -.2),
          legend.key.height = unit(.01, "cm"),
          plot.margin=unit(c(.1, .2, .1, .2), "cm"))
  
  hists_confidence_2 = ggplot() + 
    stat_interval(data = descr_plot_confidence %>% filter(iteration > 0) %>% filter(confidence == 2), aes(x = as.numeric(quantile_prior_use)/10 - .05, y = n, color = "Random"), .width = c(.5,.95), alpha = .7) + 
    stat_interval(data = descr_plot_confidence %>% filter(iteration > 0) %>% filter(confidence == 2), aes(x = as.numeric(quantile_prior_use)/10 - .05, y = n), color = "darkgreen", .width = c(.2), alpha = .7) +
    #stat_summary(data = descr_plot %>% filter(iteration > 0), aes(x = as.numeric(quantile_prior_use), y = n), geom = "line", color = "green3") + 
    geom_point(data = descr_plot_confidence %>% filter(iteration == 0) %>% filter(confidence == 2), aes(x = as.numeric(quantile_prior_use)/10 - .05, y = n, color = "Data"), shape = 15) + 
    #stat_summary(data = descr_plot %>% filter(iteration == 0), aes(x = as.numeric(quantile_prior_use), y = n), geom = "line", color = "black") +
    labs(tag = "B") +
    annotate("text", x = .8, y = 175, label = expression(epsilon == 100), size = 7/.pt)+
    xlab("Proportion of spots within prior area")+
    ylab("N trips")+
    ggtitle("Confidence: Very Likely")+
    scale_x_continuous(breaks = seq(from = 0, to = 1, by = .1))+
    scale_color_manual(values = colors)+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 6),
          axis.title.y = element_text(size = 6),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),#, hjust = .5
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = "none",
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 5, margin=margin(0,0,0,0)),
          legend.text = element_text(size = 5),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white", linewidth = .2),
          text = element_text(family = "sans"),
          legend.box.margin = margin(-.2, -.2, -.2, -.2),
          legend.key.height = unit(.01, "cm"),
          plot.margin=unit(c(.1, .2, .1, .2), "cm"))
  
  hists_confidence_3 = ggplot() + 
    stat_interval(data = descr_plot_confidence %>% filter(iteration > 0) %>% filter(confidence == 3), aes(x = as.numeric(quantile_prior_use)/10 - .05, y = n, color = "Random"), .width = c(.5,.95), alpha = .7) + 
    stat_interval(data = descr_plot_confidence %>% filter(iteration > 0) %>% filter(confidence == 3), aes(x = as.numeric(quantile_prior_use)/10 - .05, y = n), color = "darkgreen", .width = c(.2), alpha = .7) +
    #stat_summary(data = descr_plot %>% filter(iteration > 0), aes(x = as.numeric(quantile_prior_use), y = n), geom = "line", color = "green3") + 
    geom_point(data = descr_plot_confidence %>% filter(iteration == 0) %>% filter(confidence == 3), aes(x = as.numeric(quantile_prior_use)/10 - .05, y = n, color = "Data"), shape = 15) + 
    #stat_summary(data = descr_plot %>% filter(iteration == 0), aes(x = as.numeric(quantile_prior_use), y = n), geom = "line", color = "black") +
    labs(tag = "C") +
    annotate("text", x = .8, y = 175, label = expression(epsilon == 100), size = 7/.pt)+
    xlab("Proportion of spots within prior area")+
    ylab("N trips")+
    ggtitle("Confidence: Definitely")+
    scale_x_continuous(breaks = seq(from = 0, to = 1, by = .1))+
    scale_color_manual(values = colors)+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 6),
          axis.title.y = element_text(size = 6),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),#, hjust = .5
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = "none",
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 5, margin=margin(0,0,0,0)),
          legend.text = element_text(size = 5),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white", linewidth = .2),
          text = element_text(family = "sans"),
          legend.box.margin = margin(-.2, -.2, -.2, -.2),
          legend.key.height = unit(.01, "cm"),
          plot.margin=unit(c(.1, .2, .1, .2), "cm"))
  
  panel_descriptive_confidence = gridExtra::arrangeGrob(
    grobs = list(hists_confidence_1,
                 hists_confidence_2,
                 hists_confidence_3),
    widths = c(1, 1, 1),
    heights = c(1),
    layout_matrix = rbind(c(1, 2, 3))
  )
  
  ggsave(panel_descriptive_confidence, file = "6_plots/panel_descriptive_confidence.png", 
         width = 18, height = 6, bg = "white", dpi = 600, units = "cm")
  ggsave(panel_descriptive_confidence, file = "6_plots/panel_descriptive_confidence.pdf", 
         width = 18, height = 6, bg = "white", dpi = 600, units = "cm")
  ggsave(panel_descriptive_confidence, file = "6_plots/panel_descriptive_confidence.svg", 
         width = 18, height = 6, bg = "white", dpi = 600, units = "cm", device = svglite::svglite)
  
  # prior depending on navigational uncertainty
  compute_descriptive_measure_eps <- function(Feature_matrix, eps, option){
    # compute binary measure
    Feature_matrix$prior_use = ifelse(Feature_matrix$prior_distance <= eps, 1, 0)
    
    # compute trip level estimates
    if (option == "confidence"){
      prior_use_trip = Feature_matrix %>% 
        group_by(day, camera_id, unique_trip_id, iteration, confidence) %>% 
        summarize(prior_use = mean(prior_use))
      prior_use_trip$eps = eps
    } else {
      prior_use_trip = Feature_matrix %>% 
        group_by(day, camera_id, unique_trip_id, iteration) %>% 
        summarize(prior_use = mean(prior_use))
      prior_use_trip$eps = eps
      
    }
    return(prior_use_trip)
  }
  
  # marginal
  res = list()
  c=1
  for (eps in seq(from = 0, to = 500, by = 10)){
    res[[c]] = compute_descriptive_measure_eps(Feature_matrix, eps, option = "not_confidence")
    c = c + 1
  }
  res = res %>% dplyr::bind_rows()
  res = res %>% ungroup()
  res = res %>% filter(is.finite(prior_use))
  
  # by confidence
  res_confidence = list()
  c=1
  for (eps in seq(from = 0, to = 500, by = 10)){
    res_confidence[[c]] = compute_descriptive_measure_eps(Feature_matrix, eps, option = "confidence")
    c = c + 1
  }
  res_confidence = res_confidence %>% dplyr::bind_rows()
  res_confidence = res_confidence %>% ungroup()
  res_confidence = res_confidence %>% filter(is.finite(prior_use))
  
  colors = c("Random" = "darkolivegreen", "Data" = "black")
  lines = ggplot() + 
    #stat_summary(data = res %>% filter(iteration == 0), aes(x = eps, y = prior_use)) + 
    stat_lineribbon(data = res %>% filter(iteration > 0) %>% group_by(iteration, eps) %>% dplyr::summarize(prior_use = mean(prior_use, na.rm = T)), aes(x = eps, y = prior_use, color = "Random"), .width = c(.95) , size = .5, alpha = .5, fill = "darkgreen") + 
    geom_path(data = res %>% filter(iteration == 0) %>% group_by(iteration, eps) %>% summarize(prior_use = mean(prior_use, na.rm = T)), aes(x = eps, y = prior_use, color = "Data")) +
    #scale_fill_manual(values = colors)+
    scale_color_manual(values = colors)+
    coord_cartesian(ylim = c(0, 1))+
    geom_vline(xintercept = 100, linetype ="dashed")+
    labs(tag = "B")+
    xlab(expression("Navigational error"~epsilon))+
    ylab("Avg. proportion of spots in prior area")+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 6),
          axis.title.y = element_text(size = 6),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),#, hjust = .5
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
  lines
  res %>% 
    filter(iteration == 0 & eps == 100) %>% 
    group_by(unique_trip_id, eps) %>% 
    summarize(prior_use = mean(prior_use, na.rm = T)) %>% 
    ungroup() %>% 
    summarize(mean(prior_use>0))
  res %>% 
    filter(iteration > 0 & eps == 100) %>% 
    group_by(unique_trip_id, iteration, eps) %>% 
    summarize(prior_use = mean(prior_use, na.rm = T)) %>% 
    ungroup() %>% 
    summarize(mean(prior_use>0))
  
  lines2 = ggplot() + 
    #stat_summary(data = res %>% filter(iteration == 0), aes(x = eps, y = prior_use)) + 
    stat_lineribbon(data = res %>% filter(iteration > 0) %>% group_by(iteration, eps) %>% summarize(prior_use = mean(prior_use>0, na.rm = T)), aes(x = eps, y = prior_use, color = "Random"), .width = c(.95) , size = .5,alpha = .5, fill = "darkgreen") + 
    geom_path(data = res %>% filter(iteration == 0) %>% group_by(iteration, eps) %>% summarize(prior_use = mean(prior_use>0, na.rm = T)), aes(x = eps, y = prior_use, color = "Data")) +
    #scale_fill_manual(values = colors)+
    scale_color_manual(values = colors)+
    labs(tag = "A")+
    xlab(expression("Navigational error"~epsilon))+
    ylab("Proportion of trips using prior")+
    coord_cartesian(ylim = c(0, 1))+
    geom_vline(xintercept = 100, linetype ="dashed")+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 6),
          axis.title.y = element_text(size = 6),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),#, hjust = .5
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.7,.2),
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 5, margin=margin(0,0,0,0)),
          legend.text = element_text(size = 5),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white", linewidth = .2),
          text = element_text(family = "sans"),
          legend.box.margin = margin(0, 0, 0, 0),
          legend.key.height = unit(.1, "cm"),
          plot.margin=unit(c(.1, .2, .1, .2), "cm"))
 lines2
 
 res %>% 
   filter(iteration == 0 & eps == 100) %>% 
   group_by(unique_trip_id, eps) %>% 
   summarize(prior_use = mean(prior_use, na.rm = T)) %>% 
   ungroup() %>% 
   summarize(mean(prior_use))
 res %>% 
   filter(iteration > 0 & eps == 100) %>% 
   group_by(unique_trip_id, iteration, eps) %>% 
   summarize(prior_use = mean(prior_use, na.rm = T)) %>% 
   ungroup() %>% 
   summarize(mean(prior_use))
 
 
 # make lines by confidence
 lines_conf1 = ggplot() + 
   #stat_summary(data = res %>% filter(iteration == 0), aes(x = eps, y = prior_use)) + 
   stat_lineribbon(data = res_confidence %>% filter(iteration > 0) %>% filter(confidence == 1) %>% group_by(iteration, eps) %>% dplyr::summarize(prior_use = mean(prior_use, na.rm = T)), aes(x = eps, y = prior_use, color = "Random"), .width = c(.95) , size = .5, alpha = .5, fill = "darkgreen") + 
   geom_path(data = res_confidence %>% filter(iteration == 0) %>% filter(confidence == 1) %>% group_by(iteration, eps) %>% summarize(prior_use = mean(prior_use, na.rm = T)), aes(x = eps, y = prior_use, color = "Data")) +
   #scale_fill_manual(values = colors)+
   scale_color_manual(values = colors)+
   coord_cartesian(ylim = c(0, 1))+
   geom_vline(xintercept = 100, linetype ="dashed")+
   labs(tag = "A", title = "Confidence: Likely")+
   xlab(expression("Navigational error"~epsilon))+
   ylab("Avg. proportion of spots in prior area")+
   theme(panel.background = element_rect(fill = NA),
         axis.title.x = element_text(size = 6),
         axis.title.y = element_text(size = 6),
         axis.text.x = element_text(size = 5),
         axis.text.y = element_text(size = 5),#, hjust = .5
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
 lines_conf2 = ggplot() + 
   #stat_summary(data = res %>% filter(iteration == 0), aes(x = eps, y = prior_use)) + 
   stat_lineribbon(data = res_confidence %>% filter(iteration > 0) %>% filter(confidence == 2) %>% group_by(iteration, eps) %>% dplyr::summarize(prior_use = mean(prior_use, na.rm = T)), aes(x = eps, y = prior_use, color = "Random"), .width = c(.95) , size = .5, alpha = .5, fill = "darkgreen") + 
   geom_path(data = res_confidence %>% filter(iteration == 0) %>% filter(confidence == 2) %>% group_by(iteration, eps) %>% summarize(prior_use = mean(prior_use, na.rm = T)), aes(x = eps, y = prior_use, color = "Data")) +
   #scale_fill_manual(values = colors)+
   scale_color_manual(values = colors)+
   coord_cartesian(ylim = c(0, 1))+
   geom_vline(xintercept = 100, linetype ="dashed")+
   labs(tag = "B", title = "Confidence: Very Likely")+
   xlab(expression("Navigational error"~epsilon))+
   ylab("Avg. proportion of spots in prior area")+
   theme(panel.background = element_rect(fill = NA),
         axis.title.x = element_text(size = 6),
         axis.title.y = element_text(size = 6),
         axis.text.x = element_text(size = 5),
         axis.text.y = element_text(size = 5),#, hjust = .5
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
 lines_conf3 = ggplot() + 
   #stat_summary(data = res %>% filter(iteration == 0), aes(x = eps, y = prior_use)) + 
   stat_lineribbon(data = res_confidence %>% filter(iteration > 0) %>% filter(confidence == 3) %>% group_by(iteration, eps) %>% dplyr::summarize(prior_use = mean(prior_use, na.rm = T)), aes(x = eps, y = prior_use, color = "Random"), .width = c(.95) , size = .5, alpha = .5, fill = "darkgreen") + 
   geom_path(data = res_confidence %>% filter(iteration == 0) %>% filter(confidence == 3) %>% group_by(iteration, eps) %>% summarize(prior_use = mean(prior_use, na.rm = T)), aes(x = eps, y = prior_use, color = "Data")) +
   #scale_fill_manual(values = colors)+
   scale_color_manual(values = colors)+
   coord_cartesian(ylim = c(0, 1))+
   geom_vline(xintercept = 100, linetype ="dashed")+
   labs(tag = "C", title = "Confidence: Definitely")+
   xlab(expression("Navigational error"~epsilon))+
   ylab("Avg. proportion of spots in prior area")+
   theme(panel.background = element_rect(fill = NA),
         axis.title.x = element_text(size = 6),
         axis.title.y = element_text(size = 6),
         axis.text.x = element_text(size = 5),
         axis.text.y = element_text(size = 5),#, hjust = .5
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
 
 lines2_conf1 = ggplot() + 
   #stat_summary(data = res %>% filter(iteration == 0), aes(x = eps, y = prior_use)) + 
   stat_lineribbon(data = res_confidence %>% filter(iteration > 0) %>% filter(confidence == 1) %>% group_by(iteration, eps) %>% summarize(prior_use = mean(prior_use>0, na.rm = T)), aes(x = eps, y = prior_use, color = "Random"), .width = c(.95) , size = .5,alpha = .5, fill = "darkgreen") + 
   geom_path(data = res_confidence %>% filter(iteration == 0) %>% filter(confidence == 1) %>% group_by(iteration, eps) %>% summarize(prior_use = mean(prior_use>0, na.rm = T)), aes(x = eps, y = prior_use, color = "Data")) +
   #scale_fill_manual(values = colors)+
   scale_color_manual(values = colors)+
   labs(tag = "D", title = "Confidence: Likely")+
   xlab(expression("Navigational error"~epsilon))+
   ylab("Proportion of trips using prior")+
   coord_cartesian(ylim = c(0, 1))+
   geom_vline(xintercept = 100, linetype ="dashed")+
   theme(panel.background = element_rect(fill = NA),
         axis.title.x = element_text(size = 6),
         axis.title.y = element_text(size = 6),
         axis.text.x = element_text(size = 5),
         axis.text.y = element_text(size = 5),#, hjust = .5
         panel.border = element_rect(colour = "black", fill=NA),
         legend.position = c(.7,.2),
         plot.tag = element_text(size = 7, face = "bold"),
         title = element_text(size = 5, margin=margin(0,0,0,0)),
         legend.text = element_text(size = 5),
         legend.title = element_blank(),
         legend.background = element_rect(colour = "black", fill = "white", linewidth = .2),
         text = element_text(family = "sans"),
         legend.box.margin = margin(0, 0, 0, 0),
         legend.key.height = unit(.1, "cm"),
         plot.margin=unit(c(.1, .2, .1, .2), "cm"))
 lines2_conf2 = ggplot() + 
   #stat_summary(data = res %>% filter(iteration == 0), aes(x = eps, y = prior_use)) + 
   stat_lineribbon(data = res_confidence %>% filter(iteration > 0) %>% filter(confidence == 2) %>% group_by(iteration, eps) %>% summarize(prior_use = mean(prior_use>0, na.rm = T)), aes(x = eps, y = prior_use, color = "Random"), .width = c(.95) , size = .5,alpha = .5, fill = "darkgreen") + 
   geom_path(data = res_confidence %>% filter(iteration == 0) %>% filter(confidence == 2) %>% group_by(iteration, eps) %>% summarize(prior_use = mean(prior_use>0, na.rm = T)), aes(x = eps, y = prior_use, color = "Data")) +
   #scale_fill_manual(values = colors)+
   scale_color_manual(values = colors)+
   labs(tag = "E", title = "Confidence: Very Likely")+
   xlab(expression("Navigational error"~epsilon))+
   ylab("Proportion of trips using prior")+
   coord_cartesian(ylim = c(0, 1))+
   geom_vline(xintercept = 100, linetype ="dashed")+
   theme(panel.background = element_rect(fill = NA),
         axis.title.x = element_text(size = 6),
         axis.title.y = element_text(size = 6),
         axis.text.x = element_text(size = 5),
         axis.text.y = element_text(size = 5),#, hjust = .5
         panel.border = element_rect(colour = "black", fill=NA),
         legend.position = c(.7,.2),
         plot.tag = element_text(size = 7, face = "bold"),
         title = element_text(size = 5, margin=margin(0,0,0,0)),
         legend.text = element_text(size = 5),
         legend.title = element_blank(),
         legend.background = element_rect(colour = "black", fill = "white", linewidth = .2),
         text = element_text(family = "sans"),
         legend.box.margin = margin(0, 0, 0, 0),
         legend.key.height = unit(.1, "cm"),
         plot.margin=unit(c(.1, .2, .1, .2), "cm"))
 lines2_conf3 = ggplot() + 
   #stat_summary(data = res %>% filter(iteration == 0), aes(x = eps, y = prior_use)) + 
   stat_lineribbon(data = res_confidence %>% filter(iteration > 0) %>% filter(confidence == 3) %>% group_by(iteration, eps) %>% summarize(prior_use = mean(prior_use>0, na.rm = T)), aes(x = eps, y = prior_use, color = "Random"), .width = c(.95) , size = .5,alpha = .5, fill = "darkgreen") + 
   geom_path(data = res_confidence %>% filter(iteration == 0) %>% filter(confidence == 3) %>% group_by(iteration, eps) %>% summarize(prior_use = mean(prior_use>0, na.rm = T)), aes(x = eps, y = prior_use, color = "Data")) +
   #scale_fill_manual(values = colors)+
   scale_color_manual(values = colors)+
   labs(tag = "F", title = "Confidence: Definitely")+
   xlab(expression("Navigational error"~epsilon))+
   ylab("Proportion of trips using prior")+
   coord_cartesian(ylim = c(0, 1))+
   geom_vline(xintercept = 100, linetype ="dashed")+
   theme(panel.background = element_rect(fill = NA),
         axis.title.x = element_text(size = 6),
         axis.title.y = element_text(size = 6),
         axis.text.x = element_text(size = 5),
         axis.text.y = element_text(size = 5),#, hjust = .5
         panel.border = element_rect(colour = "black", fill=NA),
         legend.position = c(.7,.2),
         plot.tag = element_text(size = 7, face = "bold"),
         title = element_text(size = 5, margin=margin(0,0,0,0)),
         legend.text = element_text(size = 5),
         legend.title = element_blank(),
         legend.background = element_rect(colour = "black", fill = "white", linewidth = .2),
         text = element_text(family = "sans"),
         legend.box.margin = margin(0, 0, 0, 0),
         legend.key.height = unit(.1, "cm"),
         plot.margin=unit(c(.1, .2, .1, .2), "cm"))
 
 panel_descriptive_confidence_lines = gridExtra::arrangeGrob(
   grobs = list(lines_conf1,
                lines_conf2,
                lines_conf3,
                lines2_conf1,
                lines2_conf2,
                lines2_conf3),
   widths = c(1, 1, 1),
   heights = c(1, 1),
   layout_matrix = rbind(c(1, 2, 3),
                         c(4, 5, 6))
 )

 ggsave(panel_descriptive_confidence_lines, file = "6_plots/panel_descriptive_confidence_lines.png", 
        width = 18, height = 12, bg = "white", dpi = 600, units = "cm")
 ggsave(panel_descriptive_confidence_lines, file = "6_plots/panel_descriptive_confidence_lines.pdf", 
        width = 18, height = 12, bg = "white", dpi = 600, units = "cm")
 ggsave(panel_descriptive_confidence_lines, file = "6_plots/panel_descriptive_confidence_lines.svg", 
        width = 18, height = 12, bg = "white", dpi = 600, units = "cm", device = svglite::svglite)
 
  ################################################################################ 
  # make full panel figure spatial search
  panel_model_results = gridExtra::arrangeGrob(
    grobs = list(lines2,
                 lines,
                 hists,
                 model_results,
                 model_results_interaction,
                 model_results_time),
    widths = c(1, 1, 1),
    heights = c(1, 1.1),
    layout_matrix = rbind(c(1, 2, 3),
                          c(4, 5, 6))
  )
  
  ggsave(panel_model_results, file = "6_plots/panel_spatial_results.png", 
         width = 18, height = 12, bg = "white", dpi = 600, units = "cm")
  ggsave(panel_model_results, file = "6_plots/panel_spatial_results.pdf", 
         width = 18, height = 12, bg = "white", dpi = 600, units = "cm")
  ggsave(panel_model_results, file = "6_plots/panel_spatial_results.svg", 
         width = 18, height = 12, bg = "white", dpi = 600, units = "cm", device = svglite::svglite)
  
}


