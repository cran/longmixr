## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>",
fig.width = 6
)

## ----setup, message = FALSE---------------------------------------------------
library(longmixr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggalluvial)
library(FactoMineR)
library(factoextra)
library(lme4)
library(purrr)

## -----------------------------------------------------------------------------
data("fake_questionnaire_data")
str(fake_questionnaire_data)

## -----------------------------------------------------------------------------
fake_questionnaire_data %>% 
  filter(visit == 1) %>% 
  ggplot(aes(x = age_visit_1)) +
  geom_histogram() +
  theme_bw()

## -----------------------------------------------------------------------------
fake_questionnaire_data %>% 
  mutate(visit = as.factor(visit)) %>% 
  select(visit, starts_with("questionnaire_A")) %>% 
  pivot_longer(
    cols = -visit,
    names_to = "item",
    values_to = "level"
  ) %>% 
  ggplot(aes(x = visit, fill = level)) +
  geom_bar() +
  theme_bw() +
  facet_wrap(~item)

## -----------------------------------------------------------------------------
fake_questionnaire_data %>% 
  mutate(visit = as.factor(visit)) %>% 
  select(visit, starts_with("questionnaire_B")) %>% 
  pivot_longer(
    cols = -visit,
    names_to = "item",
    values_to = "level"
  ) %>% 
  ggplot(aes(x = visit, fill = level)) +
  geom_bar() +
  theme_bw() +
  facet_wrap(~item)

## -----------------------------------------------------------------------------
fake_questionnaire_data %>% 
  mutate(visit = as.factor(visit)) %>% 
  select(visit, starts_with("questionnaire_C")) %>% 
  pivot_longer(
    cols = -visit,
    names_to = "variable",
    values_to = "value"
  ) %>% 
  ggplot(aes(x = visit, y = value, fill = visit)) +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~variable)

## -----------------------------------------------------------------------------
fake_questionnaire_data %>%
  # only use the values from the first visit as there is only one unique value
  # per subject
  filter(visit == 1) %>% 
  ggplot(aes(y = single_continuous_variable)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

## -----------------------------------------------------------------------------
quest_A_dim <- fake_questionnaire_data %>% 
  select(starts_with("questionnaire_A")) %>% 
  FAMD(ncp = 5, graph = FALSE)

quest_B_dim <- fake_questionnaire_data %>% 
  select(starts_with("questionnaire_B")) %>% 
  FAMD(ncp = 5, graph = FALSE)

quest_C_dim <- fake_questionnaire_data %>% 
  select(starts_with("questionnaire_C")) %>% 
  prcomp(scale = TRUE)

## ---- message = FALSE---------------------------------------------------------
fviz_screeplot(quest_A_dim, main = "Questionnaire A")

## -----------------------------------------------------------------------------
fviz_screeplot(quest_B_dim, main = "Questionnaire B")

## -----------------------------------------------------------------------------
fviz_screeplot(quest_C_dim, main = "Questionnaire C")

## -----------------------------------------------------------------------------
quest_A_comp <- as.data.frame(quest_A_dim$ind$coord[, 1:3])
colnames(quest_A_comp) <- paste0("quest_A_", 1:3)
quest_B_comp <- as.data.frame(quest_B_dim$ind$coord[, 1])
colnames(quest_B_comp) <- paste0("quest_B_", 1)
quest_C_comp <- as.data.frame(quest_C_dim$x[, 1])
colnames(quest_C_comp) <- paste0("quest_C_", 1)

cluster_data <- bind_cols(
  data.frame(
    ID = fake_questionnaire_data$ID,
    visit = fake_questionnaire_data$visit,
    age = fake_questionnaire_data$age_visit_1
  ),
  quest_A_comp,
  quest_B_comp,
  quest_C_comp
)

## -----------------------------------------------------------------------------
fviz_famd_var(quest_A_dim, repel = TRUE)
fviz_contrib(quest_A_dim, "var", axes = 1)
fviz_contrib(quest_A_dim, "var", axes = 2)
fviz_contrib(quest_A_dim, "var", axes = 3)

## -----------------------------------------------------------------------------
fviz_famd_var(quest_B_dim, repel = TRUE)
fviz_contrib(quest_B_dim, "var", axes = 1)

## -----------------------------------------------------------------------------
fviz_pca_var(quest_C_dim, repel = TRUE)
fviz_contrib(quest_C_dim, "var", axes = 1)

## -----------------------------------------------------------------------------
ggplot(cluster_data, aes(x = age, y = quest_A_1)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ylab("Questionnaire A dimension 1") +
  theme_bw()

## -----------------------------------------------------------------------------
# helper function to regress out the effect
generate_residuals <- function(x, age, ID) {
  data <- data.frame(
    x = x,
    age = age,
    ID = ID)
  model <- lmer(x ~ age + (1 | ID), data = data)
  resid <- residuals(model, type = "response")
  names(resid) <- NULL
  resid
}

cluster_data_resid <- cluster_data %>% 
  # apply the function to all variables ending with a number
  mutate(across(matches("[1-9]$"),
                ~generate_residuals(x = .x, age = age, ID = ID),
                .names = "{.col}_resid")) %>% 
  select(ID, visit, ends_with("resid"), age)

## -----------------------------------------------------------------------------
ggplot(cluster_data_resid, aes(x = age, y = quest_A_1_resid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ylab("Residuals of questionnaire A dimension 1") +
  theme_bw()

## -----------------------------------------------------------------------------
response_names <- c(paste0("quest_A_", 1:3, "_resid"),
                    paste0("quest_B_", 1, "_resid"),
                    paste0("quest_C_", 1, "_resid"))

## -----------------------------------------------------------------------------
list_models <- lapply(response_names, function(x) {
  flexmix::FLXMRmgcv(as.formula(paste0(x, " ~ .")))
})

## -----------------------------------------------------------------------------
set.seed(2378)
cluster_model <- longitudinal_consensus_cluster(data = cluster_data_resid,
                                                id_column = "ID",
                                                max_k = 3,
                                                reps = 5,
                                                p_item = 0.8,
                                                model_list = list_models,
                                                flexmix_formula = as.formula("~s(visit, k = 4) | ID"),
                                                final_linkage = "ward.D2")

## -----------------------------------------------------------------------------
plot(cluster_model)

## -----------------------------------------------------------------------------
cluster_assignments <- get_clusters(cluster_model, number_clusters = 2)

## -----------------------------------------------------------------------------
original_data <- fake_questionnaire_data %>% 
  left_join(cluster_assignments, by = "ID")
cluster_data_resid <- cluster_data_resid %>% 
  left_join(cluster_assignments, by = "ID")

## -----------------------------------------------------------------------------
plot_spaghetti <- function(data = cluster_data_resid,
                           num_clus = 2,
                           var_type = c("quest_A", "quest_B", "quest_C"),
                           var_nums = 1:3,
                           scale_arg = "fixed") {
  var_type <- match.arg(var_type)
  clus_assign <- paste0("assignment_num_clus_", num_clus)
  
  plot_data <- data %>% 
    pivot_longer(
      cols = paste0(var_type, "_", var_nums, "_resid"),
      names_to = "variable"
    ) %>% 
    # this step to create a patient_idID per variable is needed, otherwise
    # ggplot can't differentiate between the different variables
    mutate(ID = paste0(ID, "_", variable),
           across(c(visit, variable, ID), as.factor))
  
  additional_data <- plot_data %>% 
    select(visit, value, variable, .data[[clus_assign]]) %>% 
    group_by(visit, variable, .data[[clus_assign]]) %>% 
    summarise(mean = mean(value),
              sd = sd(value)) %>% 
    mutate(ID = 1)
  
  ggplot(data = plot_data) +
    geom_line(mapping = aes(x = visit, y = value, col = variable, group = ID),
              alpha = 0.4) +
    geom_ribbon(data = additional_data,
                mapping = aes(x = visit, y = mean, ymin = mean - sd, ymax = mean + sd,
                              fill = variable, group = variable), alpha = 0.3) +
    geom_line(data = additional_data,
              mapping = aes(x = visit, y = mean, col = variable, group = variable),
              size = 2) +
    facet_wrap(~as.factor(.data[[clus_assign]]), scales = scale_arg) +
    theme_bw()
}

## ---- message = FALSE---------------------------------------------------------
plot_spaghetti(num_clus = 2, var_type = "quest_A", var_nums = 1:3)

## ---- message = FALSE---------------------------------------------------------
plot_spaghetti(num_clus = 2, var_type = "quest_B", var_nums = 1)

## ---- message = FALSE---------------------------------------------------------
plot_spaghetti(num_clus = 2, var_type = "quest_C", var_nums = 1)

## -----------------------------------------------------------------------------
plot_alluvial <- function(data = original_data,
                          num_clus = 2,
                          var_name) {
  clus_assign <- paste0("assignment_num_clus_", num_clus)
  
  p <- data %>% 
    ggplot(aes(x = visit,
               stratum = .data[[var_name]],
               alluvium = ID,
               fill = .data[[var_name]],
               label = .data[[var_name]])) +
    scale_fill_brewer(type = "qual", palette = "Set2") +
    geom_flow(stat = "alluvium", lode.guidance = "frontback",
              color = "darkgray") +geom_stratum() +
    facet_wrap(~as.factor(.data[[clus_assign]])) +
    theme_bw() +
    theme(legend.position = "bottom") +
    ylab("number of subjects") +
    ggtitle(paste0("Distribution of ", var_name, " across clusters"))
  
  print(p)
}

## -----------------------------------------------------------------------------
colnames(original_data)[grepl("^questionnaire_A", colnames(original_data))] %>% 
  walk(~plot_alluvial(num_clus = 2, var_name = .x))

## -----------------------------------------------------------------------------
colnames(original_data)[grepl("^questionnaire_B", colnames(original_data))] %>% 
  walk(~plot_alluvial(num_clus = 2, var_name = .x))

## -----------------------------------------------------------------------------
# first define plot function similar to the first spaghetti plot function but
# suited for the original data
plot_spaghetti_2 <- function(data,
                             var_name,
                             num_clus = 2,
                             scale_arg = "fixed") {
  
  clus_assign <- paste0("assignment_num_clus_", num_clus)
  plot_data <- data %>% 
    pivot_longer(
      cols = all_of(var_name),
      names_to = "variable"
    ) %>% 
    # this step to create a ID per variable is needed, otherwise
    # ggplot can't differentiate between the different variables
    mutate(ID = paste0(ID, "_", variable),
           across(c(visit, variable, ID), as.factor))
  
  additional_data <- plot_data %>% 
    select(visit, value, variable, .data[[clus_assign]]) %>% 
    group_by(visit, variable, .data[[clus_assign]]) %>% 
    summarise(mean = mean(value, na.rm = TRUE),
              sd = sd(value, na.rm = TRUE)) %>% 
    mutate(ID = 1)
  
  p <- ggplot(data = plot_data) +
    geom_line(mapping = aes(x = visit, y = value, col = variable, group = ID),
              alpha = 0.4) +
    geom_ribbon(data = additional_data,
                mapping = aes(x = visit, y = mean, ymin = mean - sd, ymax = mean + sd,
                              fill = variable, group = variable), alpha = 0.3) +
    geom_line(data = additional_data,
              mapping = aes(x = visit, y = mean, col = variable, group = variable),
              size = 2) +
    facet_wrap(~as.factor(.data[[clus_assign]]), scales = scale_arg) +
    theme_bw()
  
  plot(p)
}

## ---- message = FALSE---------------------------------------------------------
colnames(original_data)[grepl("^questionnaire_C", colnames(original_data))] %>% 
  walk(~plot_spaghetti_2(num_clus = 2, data = original_data, var_name = .x))

## ---- message = FALSE---------------------------------------------------------
original_data %>% 
  filter(visit == 1) %>% 
  mutate(cluster = as.factor(assignment_num_clus_2)) %>% 
  ggplot(aes(x = cluster, y = single_continuous_variable,
             fill = cluster)) +
  geom_boxplot() +
  theme_bw()

