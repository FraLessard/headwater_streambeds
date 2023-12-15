# ------------------------------------------------------------------------------
# Created by Francis Lessard and Naïm Perreault as part of the project
# on the mapping of streambeds by the Laboratoire d'Hydrologie
# Forestière de l'Université Laval
# Created on 2022-09-01
# Updated on 2023-12-15

# This script allows for an analysis of 464.7 km of stream beds
# The analysis database consists of 40,354 positions in the known
# territory, i.e. 20,177 points of presence of stream beds and 20,177 
# points of absence of stream beds. Analysis were performed using
# an ROC/AUC curve approach of logistic regressions and a 
# classification tree approach. The analysis were separated according
# to a three-class hydrological classification.
# ------------------------------------------------------------------------------

# Install the libraries ----
read_library <- function(...) {
  invisible(lapply(substitute(list(...))[-1], function(x) 
    library(deparse(x), character.only = TRUE)))
}

read_library(tidyverse,
             sf,
             magrittr,
             rpart,
             rpart.plot,
             caret,
             corrplot,
             ggrepel,
             ggfortify,
             emdbook,
             pROC,
             ROCR,
             PresenceAbsence,
             nlme,
             colorspace,
             gridExtra)







# Creation of functions ----
# Classification tree
cart <- function(train,
                 test,
                 equation, 
                 MD, # Max depth
                 NB, # Number of max and fixed branches in the classification trees
                 MS) { 
  # Response variable
  equation[2] %>% as.character -> var_rep 
  # To create a model where it is possible to extract the complexity parameter at a position in the table according to the number of max branches
  NBM = NB
  rpart(equation, train, control = rpart.control(maxdepth = MD, xval = 10, minbucket = 5, minsplit = MS, cp = 0), model=TRUE) %>% 
    .$cptable %>% 
    data.frame -> cp_table
  cp_table$nsplit %>% 
    unique %>% 
    subset(. <= NBM) %>% 
    max -> NBM
  cp_table$CP[cp_table$nsplit == NBM] -> cp
  # Creation of the tree
  rpart(equation, train, control = rpart.control(maxdepth = MD, xval = 10, minbucket = 5,  minsplit = MS, cp = cp*1.05), model=TRUE) %>% 
    {model_cart <<- .} %>% 
    rpart.plot
  # Prediction of the tree
  test$predict <- predict(model_cart, test, type = "class", control = rpart.control(maxdepth = MD, xval = 10, minbucket = 5,  minsplit = MS, cp = cp*1.05))
  # Extraction of tree performance statistics
  confusionMatrix(data=test$predict, reference=test %>% pull(var_rep), positive = "1") %>% 
    {scores <<- .}
  
  # To perform a pediction on the entire database
  predict(model_cart, test, type = "class", control = rpart.control(maxdepth = MD, xval = 10, minbucket = 5,  minsplit = MS, cp = cp*1.05))
}

# Hydrological classification
hydrological.classification <- function(data){
  data$chs <- ifelse (data$dep %in% c("1A","1AA","1AAM","1AAR","1AAY","1AB","1AD","1ADY","1AM","1AR","1ASY",
                                      "1AY","1AYR","1M","1Y","M1","M1A","M1AA","M6S","M7T","M8A","M8AP","M8C","M8PY","R",
                                      "R1","R1A","R1AA","R1BD","R2A","R2AK","R2BE","R3AN","R4","R4GA","R4GS","R5A","R5S",
                                      "R6","R6S","R7","R7T","R8A","R8AP","R8C","R8E","R8P","R9S","RS"), "Shallow soil","")
  
  data$chs <- ifelse (data$ter %in% c("EAU","DH","INO","AL") 
                      | data$dep %in% c("4A","4AR","4AY","4GA","4GAM","4GAR","4GAY","5AM","5AR","5AY",
                                        "7","7AN","7E","7L","7R","7T","7TM","7TY"), "Thick soil with low infiltration rate", data$chs)
  
  data$chs <- ifelse (data$ter %in% c("A","AER","ANT","CAR","GR","HAB","MI","PAI","US","VIL","RO") 
                      | data$dep %in% c("1B","1BC","1BD","1BDY","1BF","1BG","1BI","1BIM","1BIY","1BN","1BP",
                                        "1BPY","1BR","1BT","1P","2","2A","2AE","2AK","2AM","2AR","2AT","2AY","2B","2BD",
                                        "2BDY","2BE","2BEM","2BER","2BEY","2BP","2BR","3","3A","3AC","3AE","3AN","3ANY",
                                        "3D","3DD","3DE","4","4GD","4GS","4GSM","4GSY","4P","5G","5GR","5GSR","5L","5R",
                                        "5S","5SM","5SR","5SY","5Y","6","6A","6AM","6AY","6R","6S","6SM","6SR","6SY","8",
                                        "8A","8AC","8AL","8ALM","8ALY","8AM","8AP","8APM","8APY","8AR","8AS","8ASY","8AY",
                                        "8AYP","8C","8CM","8CY","8E","8F","8G","8M","8P","8PM","8PY","8Y","9","9A","9R",
                                        "9S","9SM","9SY","5A"), "Thick soil with high infiltration rate", data$chs)
  
  data$chs <- ifelse(data$chs == "", "Shallow soil", data$chs)
  
  data$chs <- as.factor(data$chs)
  
  return(data)
}





# Import the data ----
rstudioapi::getSourceEditorContext()$path %>%
  dirname %>% 
  paste0(.,"/data.gpkg") %>% 
  st_read(layer = "presence_absence") %>%
  mutate(sbp = factor(sbp),
         type = factor(type),
         dep = factor(dep),
         ter = factor(ter),
         cer = factor(cer)) %>%
  mutate(slope = tan(slope)*100) %>% # Convert slope from radians to percent 
  hydrological.classification -> data
# sbp = Streambed presence : 0 = Absence / 1 = Presence
# type = Type of stream : Absence, Diffuse, Intermittent or Perennial

# Channel head
rstudioapi::getSourceEditorContext()$path %>%
  dirname %>% 
  paste0(.,"/data.gpkg") %>% 
  st_read(layer = "channel_head") %>%
  mutate(type = factor(type)) %>% 
  mutate(dep = factor(dep)) %>%
  mutate(ter = factor(ter)) %>% 
  mutate(mean_corr = mean*9) %>%
  mutate(slope = tan(slope)*100) %>% 
  hydrological.classification %>% 
  filter(type == "Channel head") -> channelhead
# type = Type of theshold : Channel head and stream head
# Get the median of channel head by hydrological classes
channelhead %>% 
  st_drop_geometry %>% 
  group_by(chs) %>% 
  summarise(med = median(mean_corr)) -> channelhead_median

channelhead$chs %>% table




  
# DownSample by ecological reference framework natural region ----
i <- 1
for (c in data$cer %>% unique) {
  set.seed(1) # Set.seed to maintain replicability
  data %>% 
    filter(cer == c) %>% 
    downSample(.$sbp) -> data_temp
  if(i == 1){ # This condition allows to build the dataset in an iterative way
    data1 <- data_temp
  }else{
    data1 %<>% 
      rbind(data_temp)
  }
  i <- i + 1
}
rm(data_temp,i,c)

data$chs %>%  table
data1$sbp %>% table # The number of 0 (absence of a streambed) is equal to the number of 1 (presence of a streambed)
data1 %>% group_by(sbp, cer) %>% summarise(l = length(sbp)) %>% pivot_wider(names_from = sbp, values_from = l) %>% mutate(p = `0`/(`0` + `1`))


# To keep the geometry in the dataset and also have a version without geometry
data1 %>% 
  dplyr::select(-Class) %>% 
  st_as_sf %>% 
  {data1_geom <<-.} %>% 
  st_drop_geometry -> data1

rm(data)





# Correlation matrix ----
data1 %>%
  dplyr::select(D8, PROB, TPI) %>% # Selection of only numeric variables
  cor %>% 
  corrplot(method = "number", type = "upper", addCoef.col = "black", tl.srt = 45)






# Logistic regression of analysis variables ----
# Creation of the results table
matrix(nrow=9,ncol=4) %>% 
  as_tibble() -> results_sens
colnames(results_sens) <- c("chs", "variable", "roc", "auc")

data1 %>% 
  dplyr::select(where(is.numeric), -grhq) %>% 
  names %>% 
  as.list -> var

mod_list <- list()

i <- 1
# Loop that allows to iterate for each hydrological class (chs)
for(c in data1_geom %>% pull(chs) %>% unique){
  # Loop that allows to iterate for each variable of interest
  for(v in seq_along(var)){
    # Creation of a simplified dataset with only the response variable (sbp) and the variable of interest
    data1_geom %>% 
      st_drop_geometry %>% 
      filter(chs == c) %>% 
      dplyr::select(sbp, var[[v]]) -> data1_geom_temp
    
    mod <- glm(sbp ~ ., family = binomial(link = logit), data = data1_geom_temp) # Creation of a logistic regression glm
    mod_list[[i]] <- mod
    predict(mod, data1_geom_temp, type = "response") -> data1_geom_temp$predict # Creation of a logistic regression prediction field
    # ROC-AUC
    # To adjust the width of the graph on the x-axis
    par(pty = "s") 
    # Allows to create the roc curve
    roc(data1_geom_temp$sbp, data1_geom_temp$predict, plot = TRUE, legacy.axes = TRUE, percent = TRUE, print.auc = TRUE,
        xlab = "False Positive Percentage", ylab = "True Positive Percentage", main = paste0(c, " / ", var[[v]])) -> rocc
    
    results_sens$chs[i] <- c # Hydrological classification
    results_sens$variable[i] <- paste(var[[v]],collapse="") # variable of interest
    results_sens$roc[i] <- list(rocc) # Roc curve
    results_sens$auc[i] <- rocc[["auc"]] %>% as.numeric # AUC of the roc curve
    i <- i+1
  } 
}

# To extract the sensitivity/specificity (variables of interest)
results_sens_temp <- list() # Creating an empty list to extract the sensitivity/specificity
i <- 1
for(v in results_sens %>% pull(variable) %>% unique) { # Extraction of variables of interest only
  for(c in results_sens %>% pull(chs) %>% unique) { # For all hydrological classification
    results_sens %>% 
      filter(chs == c) %>% 
      filter(variable == v) %>% 
      pull(roc) %>% # The following 3 lines are used to extract all the sensitivity/specificity values from the roc curve
      .[[1]] %>% 
      {cbind(.$sensitivities, .$specificities)} %>% 
      as_tibble %>%
      rename(sensitivities =  V1, specificities = V2) %>% 
      mutate(chs = c, variable = v) -> results_sens_temp[[i]] # Allows to keep in memory all the tables extracted in a list
    i <- i +1
  } 
}
bind_rows(results_sens_temp) -> results_sens_temp # All tables in the list are merged

rm(var, i, c, v, data1_geom_temp, mod, rocc)





# Creation of models ----
# Model 1 : GRHQ
# The "Near" tool of ArcGIS has previously been used with all segments of the GRHQ, a distance of less than 6 m is considered an accurate model
data1_geom %<>% 
  mutate(mod1 = ifelse(grhq < 6, 1, 0) %>% factor)





# Model 2 : D8 - Median
# Everything above the median of the D8 with a focal statistic of 6 m by hydrological classification is considered an accurate model
data1_geom %>% 
  filter(chs == "Shallow soil") %>% 
  mutate(mod2 = ifelse(D8 > channelhead_median$med[[1]], 1, 0) %>% factor) -> data1_geom_g

data1_geom %>% 
  filter(chs == "Thick soil with high infiltration rate") %>% 
  mutate(mod2 = ifelse(D8 > channelhead_median$med[[2]], 1, 0) %>% factor) -> data1_geom_i

data1_geom %>% 
  filter(chs == "Thick soil with low infiltration rate") %>% 
  mutate(mod2 = ifelse(D8 > channelhead_median$med[[3]], 1, 0) %>% factor) -> data1_geom_pi

data1_geom <- rbind(data1_geom_g, data1_geom_i, data1_geom_pi)

rm(data1_geom_g, data1_geom_i, data1_geom_pi)





# Model 3 : D8 - Max kappa
# The variable D8 was optimized by hydrological classification according to the maximum Kappa
data1_geom %>% 
  filter(chs == "Shallow soil") %>% 
  mutate(ID = 1:nrow(.)) %>% 
  mutate(pred = predict(mod_list[[1]], ., type = "response")) %>% 
  {data1_geom_g <<- .} %>% 
  st_drop_geometry %>% 
  dplyr::select(ID, sbp, pred) %>% 
  mutate(ID = ID %>% as.factor,
         sbp = sbp %>% as.character %>% as.numeric,
         pred = pred %>% as.numeric) %>%
  optimal.thresholds(req.sens = 0.80, req.spec = 0.80)

data1_geom_g %<>% 
  mutate(mod3 = ifelse(pred > 0.23, 1, 0) %>% factor) %>% 
  dplyr::select(-ID, -pred)



data1_geom %>% 
  filter(chs == "Thick soil with high infiltration rate") %>% 
  mutate(ID = 1:nrow(.)) %>% 
  mutate(pred = predict(mod_list[[4]], ., type = "response")) %>% 
  {data1_geom_i <<- .} %>% 
  st_drop_geometry %>% 
  dplyr::select(ID, sbp, pred) %>% 
  mutate(ID = ID %>% as.factor,
         sbp = sbp %>% as.character %>% as.numeric,
         pred = pred %>% as.numeric) %>%
  optimal.thresholds(req.sens = 0.80, req.spec = 0.80)

data1_geom_i %<>% 
  mutate(mod3 = ifelse(pred > 0.5, 1, 0) %>% factor) %>% 
  dplyr::select(-ID, -pred)



data1_geom %>% 
  filter(chs == "Thick soil with low infiltration rate") %>%
  mutate(ID = 1:nrow(.)) %>% 
  mutate(pred = predict(mod_list[[7]], ., type = "response")) %>% 
  {data1_geom_pi <<- .} %>% 
  st_drop_geometry %>% 
  dplyr::select(ID, sbp, pred) %>% 
  mutate(ID = ID %>% as.factor,
         sbp = sbp %>% as.character %>% as.numeric,
         pred = pred %>% as.numeric) %>%
  optimal.thresholds(req.sens = 0.80, req.spec = 0.80)

data1_geom_pi %<>% 
  mutate(mod3 = ifelse(pred > 0.39, 1, 0) %>% factor) %>% 
  dplyr::select(-ID, -pred)



data1_geom <- rbind(data1_geom_g, data1_geom_i, data1_geom_pi)
rm(data1_geom_g, data1_geom_i, data1_geom_pi)





# Model 4 : Classification tree
# Shallow soil
data1_geom %>% 
  filter(chs == "Shallow soil") -> data1_geom_g
cart(data1_geom_g, data1_geom_g, sbp~PROB+TPI, 3, 6, nrow(data1_geom_g)*0.02) -> data1_geom_g$mod4
model_cart -> model_cart_g
scores

data1_geom %>% 
  filter(chs == "Thick soil with high infiltration rate") -> data1_geom_i
cart(data1_geom_i, data1_geom_i, sbp~PROB+TPI, 3, 6, nrow(data1_geom_i)*0.02) -> data1_geom_i$mod4
model_cart -> model_cart_i
scores

data1_geom %>% 
  filter(chs == "Thick soil with low infiltration rate") -> data1_geom_pi
cart(data1_geom_pi, data1_geom_pi, sbp~PROB+TPI, 3, 6, nrow(data1_geom_pi)*0.02) -> data1_geom_pi$mod4
model_cart -> model_cart_pi
scores

# Combining the dataset once each tree is applied to each hydrological class
data1_geom <- rbind(data1_geom_g, data1_geom_i, data1_geom_pi)

rm(data1_geom_g, data1_geom_i, data1_geom_pi, model_cart, scores)





# Creation of the results table, with all the data and by hydrological classes for the figures
results <- matrix(nrow=16,ncol=5)
colnames(results) <- c("chs","model","kappa","sensitivity","specificity")
results %<>% as_tibble()
results$model <- c("GRHQ",
                 "Channel head",
                 "Max Kappa",
                 "CART") %>% rep(4)

# To calculate the performance of models for all hydrological classes combined
for(i in 1:4){
  data1_geom %>%
    st_drop_geometry -> data1_geom_temp
  mod <- paste0("mod",i)
  confusionMatrix(data=data1_geom_temp[,mod], reference=data1_geom_temp$sbp, positive = "1") -> scores
  results$chs[i] <- "All"
  results$kappa[i] <- scores$overall[2]
  results$sensitivity[i] <- scores$byClass[1]
  results$specificity[i] <- scores$byClass[2]
}

# To calculate the performance of the models for each of the hydrological classes
i <- 5
for(c in data1 %>% pull(chs) %>% unique){
  for(id in 1:4){
    data1_geom %>%
      st_drop_geometry %>%
      filter(chs == c) -> data1_geom_temp
    mod <- paste0("mod",id)
    confusionMatrix(data=data1_geom_temp[,mod], reference=data1_geom_temp$sbp, positive = "1") -> scores
    results$chs[i] <- c
    results$kappa[i] <- scores$overall[2]
    results$sensitivity[i] <- scores$byClass[1]
    results$specificity[i] <- scores$byClass[2]
    i <- i + 1
  }
}

# To calculate the omission and commission and arrange models
results %<>% 
  mutate(model = factor(model, levels = c("GRHQ",
                                          "Channel head",
                                          "Max Kappa",
                                          "CART"))) %>% 
  mutate(omission = 1 - sensitivity) %>% 
  mutate(commission = 1 - specificity)

# To remove the performance of all hydrological classes combined
results %<>% 
  filter(chs != "All")

rm(c, i, id, mod, scores, data1_geom_temp)





# Figures ----
# Channel head
channelhead %>% 
  ggplot +
  theme_classic(base_size = 10) +
  theme(panel.grid.major.y = element_line(colour = "darkgrey", linetype = "dashed"),
        text = element_text(size=20),
        title = element_text(size=20),
        axis.text.x = element_text(color="black")) +
  guides(color = "none",
         size = "none") +
  xlab("Hydrological class") +  
  ylab("Drainage area (ha)") +
  scale_x_discrete(labels = c("Shallow soil\n(n = 542)","Thick soil\nwith high\ninfiltration\nrate\n(n = 411)","Thick soil\nwith low\ninfiltration\nrate\n(n = 80)")) + 
  scale_y_continuous(breaks=c(0,1,2,3,4,5,10,20), label = c(0,"", "", "", "", 5, 10, 20)) +
  scale_color_manual(values = c("deepskyblue", "darkolivegreen3", "darkorchid4")) +
  geom_boxplot(aes(chs, mean_corr/10000, color = chs), width = 0.8) +
  coord_cartesian(ylim = c(1,25)) #+
  #geom_text(data = channelhead_median ,mapping = aes(x = c(1.5,2.5,3.5), y = med, label = paste0(round(med, digits = 2), " ha")), size = 8)
rstudioapi::getSourceEditorContext()$path %>%
  dirname %>% 
  paste0(.,"/channel_head_figure_3.jpg") %>% 
  ggsave(width=6, height=8, dpi = 300)


# Classification tree
rstudioapi::getSourceEditorContext()$path %>%
  dirname %>% 
  paste0(.,"/classification_trees_figure_4.jpg") %>% 
  jpeg(width=10, height=4.5, units = 'in', res = 300)

layout(mat = matrix(c(1, 2, 3), nrow = 1, ncol = 3),
       heights = c(1), # Heights of the two rows
       widths = c(1)) # Widths of the two columns

par(mar = c(0, 0, 2, 0), cex = 1, pty = "s")
rpart.plot(model_cart_g, box.palette = "RdYlGn")
mtext("Shallow soil", side=3, cex = 1)
par(mar = c(0, 0, 2, 0), cex = 1.5, pty = "s")
rpart.plot(model_cart_i, box.palette = "RdYlGn")
mtext("Thick soil with high\ninfiltration rate", side=3, cex = 1)
par(mar = c(0, 0, 2, 0), cex = 1.5, pty = "s")
rpart.plot(model_cart_pi, box.palette = "RdYlGn")
mtext("Thick soil with low\ninfiltration rate", side=3, cex = 1)

dev.off()


# ROC-AUC curve and models
results_sens %>% 
  dplyr::select(-roc) %>% # We remove the ROC curve since it makes the table very heavy
  mutate(auc = round(auc/100,3)) -> label_table # Allows you to create a table containing all the information necessary to put the AUC in the text chart
ggplot() +
  geom_line(data = results_sens_temp, 
            mapping =  aes((100-specificities)/100, sensitivities/100, color = variable), size=1) + # Allows you to create a table containing all the information necessary to put the AUC as a label in the figure
  geom_point(data = results,
             mapping = aes(1-specificity, sensitivity, shape = model), size = 5) +
  scale_shape_manual(values = c(16,15,18,17), labels = c("GRHQ", "Channel head", "Max Kappa", "CART")) +
  scale_color_manual(values=c("brown1", "darkgreen", "cadetblue3")) +
  facet_wrap(~chs, scales='free', labeller = label_wrap_gen(multi_line = TRUE)) +
  theme_bw() +
  guides(shape = guide_legend(override.aes = list(size=5), title="Model")) +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background=element_rect(colour="white", fill="white"),
        text = element_text(size=28)) + # Allows to put the figure with a white background
  guides(color = "none") + # Same
  geom_abline(slope = 1, linetype = 2) + # Allows you to set the diagonal of a ROC curve (no effect)
  xlab("False positive rate") +
  ylab("True positive rate") +
  geom_text_repel(x = 0.80, y = 0.35, seed = 5, max.overlaps = 20, segment.color = NA, 
                  mapping = aes(label = paste0(variable ," : ",auc)), 
                  data = label_table , size = 7, direction = "y", 
                  color = c("brown1", "darkgreen", "cadetblue3") %>% rep(3)) + # To put the AUC in text in the figure
  geom_text_repel(results, mapping = aes(x = commission, y = 1-omission, label = round(kappa, digits = 2)), nudge_x = 0.2, size = 7)
rstudioapi::getSourceEditorContext()$path %>%
  dirname %>% 
  paste0(.,"/roc_auc_figure_5.jpg") %>% 
  ggsave(width=20, height=10, dpi = 300)


# Supplementary materials - Slope
for(c in channelhead$chs %>% unique){
  channelhead %>%
    st_drop_geometry %>% 
    mutate(logd8 = log10(mean_corr)) %>% 
    mutate(logd8cat = cut(logd8, seq(3,6,0.2))) %>% 
    group_by(logd8cat, chs) %>% 
    summarise(n = n(),
              slp = mean(slope),
              ld8 = mean(logd8)) %>% 
    filter(chs == c) -> channelhead_group  
  
  data1 %>%
    st_drop_geometry %>%
    filter(chs == c) %>% 
    mutate(logd8 = log10(D8)) %>% 
    mutate(logd8cat = ntile(logd8, ceiling(nrow(.)/100))) %>%
    group_by(logd8cat, chs) %>% 
    summarise(slp = mean(slope),
              ld8 = mean(logd8)) %>% 
    ggplot +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background=element_rect(colour="white", fill="white"),
          text = element_text(size=15)) +
    xlab("Drainage area (ha)") +
    ylab("Slope (%)") +
    labs(size = "Channel head (n)") +
    geom_point(aes(ld8, slp)) +
    geom_smooth(aes(ld8, slp), se = F) +
    geom_point(aes(ld8, slp, size = n), channelhead_group, color = "red") +
    scale_size_continuous(breaks = ceiling(seq(min(channelhead_group$n), max(channelhead_group$n), length.out = 4))) +
    scale_x_continuous(breaks = c(2:7), labels = c("0.01",
                                                   "0.1",
                                                   "1",
                                                   "10",
                                                   "100",
                                                   "1 000")) -> g
  print(g)
  rstudioapi::getSourceEditorContext()$path %>%
    dirname %>% 
    paste0(.,"/slope_area_",c,".jpg") %>% 
    ggsave(width=12, height=6, dpi = 300)
}
