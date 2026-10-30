# Directories and Data Import ---------------------------------------------

output_path <- paste0("./output/", results_path)

dir.create(paste0(output_path, "/self-assessment"))


# Part 1 - Protocols ------------------------------------------------------

form1 <- read_excel("./other-data/self-assessment/Dificuldades Experimentais (Responses) - manual edit.xlsx")

inclusion_sets$`Validation - Decision`[inclusion_sets$`Validation - Decision` == "Low sample size"] <- "Excluded"

inclusion_sets_short <- inclusion_sets |>
  filter(UNIT == "BRI") |>
  select(c("EXP", "LAB", "Validation - Mean score", "Validation - Decision", "Validation - Excluded by")) |>
  filter(!is.na(`Validation - Mean score`))

form1_val <- full_join(form1, inclusion_sets_short, by = c("LAB", "EXP"))
form1_val$method <- str_extract(form1_val$EXP, "[A-Za-z]+")

#---
## Plot 1 - Agreement with Assigned Score ----
#---

plot_score_agree <- ggplot(
  form1_val |> filter(!is.na(score_agree)),
  aes(
    x = factor(score_agree, levels = c("Sim", "Nao"), labels = c("yes", "no")),
    y = `Validation - Mean score`, color = method
  )
) +
  geom_jitter(height = 0, width = 0.2, alpha = 0.6) +
  bri_theme +
  xlab("Does the lab agree with the assigned score?") +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = bri_3color)

bri_ggsave(paste0(output_path, "/self-assessment/", "score_agree.png"), plot = plot_score_agree)

#---
## Adjusting validation decisions ----
#---

form1_val <- form1_val |> mutate(Validation_new = if_else(
  condition = `Validation - Decision` == "Included",
  true = `Validation - Decision`,
  false = `Validation - Excluded by`
))

#---
## Plot 2 - Agreement with Validation Decision -----
#---

plot_decision_agree <- ggplot(
  form1_val |> filter(!is.na(decision_agree)),
  aes(
    x = factor(decision_agree, levels = c("Sim", "Nao"), labels = c("yes", "no")),
    y = `Validation_new`, color = method
  )
) +
  geom_jitter(height = 0.1, width = 0.2, alpha = 0.6) +
  bri_theme +
  xlab("Does the lab agree with the decision?") +
  ylab("Reason for exclusion") +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = bri_3color) +
  scale_y_discrete(limits = c(
    "Inadequate experimental unit", "Insufficient sample size", "Lack of documentation",
    "Non-interpretable results", "Protocol", "Protocol and non-interpretable results",
    "", "Included"
  ))

bri_ggsave(paste0(output_path, "/self-assessment/", "decision_agree.png"), plot = plot_decision_agree)

#---
## Summary Statistics -----
#---

form1_val |>
  filter(!is.na(decision_agree)) |>
  group_by(`Validation - Excluded by`) |>
  count(decision_agree)

form1_val |>
  filter(!is.na(decision_agree)) |>
  group_by(`Validation - Excluded by`, method) |>
  count(decision_agree)

#---
## Function to Count Matches  -----
#---

get_n_matches <- function(data, question, option) {
  option_data <- data
  length(grep(option, option_data[[question]]))
}

protocol_options <- c(
  "1a", "1b", "1c", "1d", "1e", "1f",
  "2a", "2b", "2c", "2d", "2e", "2f", "2g", "2h",
  "3a", "3b", "3c", "3d",
  "4a", "4b", "4c", "4d", "4e", "4f",
  "5a", "5b", "5c",
  "6", "7", "8", "9", "9b"
)

cat_protocol <- sapply(data = form1_val, question = "reasons_cat", X = protocol_options, FUN = get_n_matches) |>
  as.data.frame() |>
  rename("category_counts" = c(1)) |>
  rownames_to_column("protocol_options")

#---
## Protocol Labels  -----
#---

protocol_labels <- c(
  "1a" = "(1a) Different cells and animals",
  "1b" = "(1b) Model behaves differently",
  "1c" = "(1c) Other part behaves differently",
  "1d" = "(1d) Validity criteria unfeasible",
  "1e" = "(1e) Changes after pilot",
  "1f" = "(1f) Experimental reasons (general)",
  "2a" = "(2a) Lab's infrastructure",
  "2b" = "(2b) Regulatory requirements",
  "2c" = "(2c) Expiration of reagents",
  "2d" = "(2d) Deadlines",
  "2e" = "(2e) Costs",
  "2f" = "(2f) Pandemic restrictions",
  "2g" = "(2g) Logistics or infrastructure (general)",
  "2h" = "(2*) Restrictions of the animal facility",
  "3a" = "(3a) Team's expertise",
  "3b" = "(3b) Standard lab protocol",
  "3c" = "(3c) Deliberate choice",
  "3d" = "(3*) Recommendations from supplier",
  "4a" = "(4a) Lab's mistake in developing protocol",
  "4b" = "(4b) Lab's mistake in understanding protocol",
  "4c" = "(4c) Lab's mistake in planning",
  "4d" = "(4d) Lab's mistake in executing experiment",
  "4e" = "(4e) Lab's mistake (general)",
  "4f" = "(4*) Protocol not clear",
  "5a" = "(5a) Coordination team's mistake in acquiring materials",
  "5b" = "(5b) Coordination team's mistake in communicating",
  "5c" = "(5c) Coordination team's mistake (general)",
  "6" = "Other",
  "7" = "No deviations noted",
  "8" = "Answer is not relevant to question",
  "9" = "No answer",
  "9b" = "Lists changes without reasons"
)

#---
## Plot 3 - Protocol Deviation Reasons  -----
#---

plot_protocol_reasons <- ggplot(cat_protocol, aes(y = protocol_options, x = category_counts)) +
  geom_bar(stat = "identity") +
  bri_theme +
  labs(y = "What are the main reasons for the protocol deviations observed?", x = "Number of experiments") +
  scale_y_discrete(labels = protocol_labels, limits = rev(protocol_options))

bri_ggsave(paste0(output_path, "/self-assessment/", "protocol_reasons_complete.png"), plot = plot_protocol_reasons, height = 30, width = 25, units = "cm")

#---
## Plot 4 - Simplified Protocol Deviation Reasons ----
#---

plot_protocol_reasons_simplified <- ggplot(cat_protocol, aes(y = protocol_options, x = category_counts)) +
  geom_bar(stat = "identity") +
  bri_theme +
  labs(y = "What are the main reasons for the protocol deviations observed?", x = "Number of experiments") +
  scale_y_discrete(labels = protocol_labels, limits = rev(c(
    "1a", "1b", "1c", "1d", "1e", "1f", " ",
    "2a", "2b", "2c", "2d", "2e", "2f", "2g", "2h", " ",
    "3a", "3b", "3c", "3d", " ",
    "4a", "4c", "4d", "4e", "4f", " ",
    "5a", "5b", " ",
    "6"
  ))) +
  labs(caption = "Responses could include multiple categories. Total number of experiments = 143, of which 87 provided responses.")

bri_ggsave(paste0(output_path, "/self-assessment/", "protocol_reasons.png"), plot = plot_protocol_reasons_simplified, height = 25, width = 20, units = "cm")

#---
## Figure 5 ----
#---

protocol_labels_fig <- c(
  # Experimental Issues
  "1a" = "Different cells or animals",
  "1b" = "Model behaves differently",
  "1c" = "Experiment works differently",
  "1d" = "Unfeasible validity criteria",
  "1e" = "Changes due to pilot",
  "1f" = "Experimental issues (general)",
  # Infrastructure/Logistics
  "2a" = "Lab infrastructure",
  "2b" = "Regulatory requirements",
  "2c" = "Expired reagents",
  "2d" = "Deadlines",
  "2e" = "Costs",
  "2f" = "Pandemic restrictions",
  "2g" = "Infrastructure/logistics (general)",
  "2h" = "Animal facility limitations",
  # Deliberate Choice
  "3a" = "Team expertise",
  "3b" = "Standard lab protocol",
  "3c" = "Deliberate choice (general)",
  "3d" = "Supplier recommendations",
  # Lab Error
  "4a" = "Error in protocol development",
  "4b" = "Lab's mistake in \nunderstanding protocol",
  "4c" = "Error in planning",
  "4d" = "Error in executing experiment",
  "4e" = "Lab error (general)",
  "4f" = "Unclear protocol",
  # Coordinating Team Error
  "5a" = "Error in acquiring materials",
  "5b" = "Error in communication with lab",
  "5c" = "Coord. team's mistake \n(general)",
  "6" = "Other",
  "7" = "No deviations noted",
  "8" = "Answer is not relevant to question",
  "9" = "No answer",
  "9b" = "Lists changes without reasons"
)

protocol_labels_fig[["6"]] <- "**Other**"

y_levels <- rev(c(
  "**Experimental issues**",
  "1b", "1c", "1a", "1e", "1d", "1f",
  "**Infrastructure/logistics**",
  "2a", "2h", "2b", "2c", "2d", "2f", "2e", "2g",
  "**Deliberate choice**",
  "3b", "3d", "3a", "3c",
  "**Lab error**",
  "4a", "4c", "4d", "4f", "4e",
  "**Coordinating team error**",
  "5a", "5b", " ",
  "6"
))

grid_positions <- setdiff(
  seq_along(y_levels),
  which(y_levels %in% c(
    "**Experimental issues**",
    "**Infrastructure/logistics**",
    "**Deliberate choice**",
    "**Lab error**",
    "**Coordinating team error**",
    " "
  ))
)

plot_protocol_reasons_final <- ggplot(cat_protocol, aes(y = protocol_options, x = category_counts)) +
  geom_hline(yintercept = grid_positions, color = "grey90") +
  geom_bar(stat = "identity") +
  bri_theme +
  labs(y = "What are the main reasons for the protocol deviations observed?", x = "Number of replications") +
  scale_y_discrete(labels = protocol_labels_fig, limits = y_levels, expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(expand = FALSE, clip = "off") +
  theme(
    axis.text.y = ggtext::element_markdown(),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.margin = unit(c(5, 5, 5, 5), "mm"),
    axis.ticks.y = element_line(color = c(
      "black", "white", "black", "black", "white",
      "black", "black", "black", "black", "black",
      "white", "black", "black", "black", "black",
      "white", "black", "black", "black", "black",
      "black", "black", "black", "black", "white",
      "black", "black", "black", "black", "black",
      "black", "white"
    ))
  )

protocol_labels_fig[["6"]] <- "Other"

bri_ggsave(
  paste0(output_path, "/self-assessment/", "Figure 5.png"),
  plot = plot_protocol_reasons_final,
  height = 25, width = 20, units = "cm"
)

#---
## Table S23 ----
#---

tbl_s23 <- cat_protocol |>
  mutate(protocol_labels = protocol_labels_fig) |>
  slice(1:5) |>
  arrange(desc(category_counts)) |>
  add_row(protocol_labels = "Experimental issues", .before = 1) |>
  add_row(
    cat_protocol |>
      mutate(protocol_labels = protocol_labels_fig) |>
      slice(6)
  ) |>
  add_row(protocol_labels = "Infrasctructure/logistics") |>
  add_row(
    cat_protocol |>
      mutate(protocol_labels = protocol_labels_fig) |>
      slice(7:12) |>
      add_row(cat_protocol |>
        mutate(protocol_labels = protocol_labels_fig) |>
        slice(14)) |>
      arrange(desc(category_counts)) |>
      add_row(cat_protocol |>
        mutate(protocol_labels = protocol_labels_fig) |>
        slice(13))
  ) |>
  add_row(protocol_labels = "Deliberate choice") |>
  add_row(
    cat_protocol |>
      mutate(protocol_labels = protocol_labels_fig) |>
      slice(15:18) |>
      arrange(desc(category_counts))
  ) |>
  add_row(protocol_labels = "Lab error") |>
  add_row(
    cat_protocol |>
      mutate(protocol_labels = protocol_labels_fig) |>
      slice(19:22) |>
      add_row(cat_protocol |>
        mutate(protocol_labels = protocol_labels_fig) |>
        slice(24)) |>
      arrange(desc(category_counts)) |>
      add_row(cat_protocol |>
        mutate(protocol_labels = protocol_labels_fig) |>
        slice(23))
  ) |>
  add_row(protocol_labels = "Coordinating team error") |>
  add_row(
    cat_protocol |>
      mutate(protocol_labels = protocol_labels_fig) |>
      slice(25:26) |>
      add_row(cat_protocol |>
        mutate(protocol_labels = protocol_labels_fig) |>
        slice(28))
  ) |>
  select(protocol_labels, category_counts) |>
  rename(Count = category_counts) |>
  rename(Label = protocol_labels) |>
  flextable() |>
  bold(i = c(1, 8, 17, 22, 29, 32)) |> 
  set_table_properties(layout = "autofit") |>
  add_name("Table - Examples of reasons for protocol deviations.") |>
  add_footer_lines("Table shows the various categories and subcategories for protocol deviations (in order of general category frequency, as in Figure 5), with illustrative examples of lab responses for each subcategory.")

save_tbl(
  list(tbl_s23),
  paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table S23.docx")
)
cat("\n### Table S23 generated! ###\n")

#---
## Plot 5 - Could the deviations have been prevented? ----
#---

plot_preventable_deviations <- ggplot(
  form1_val |> filter(!is.na(preventable)),
  aes(
    alpha = factor(preventable, levels = c("Nao", "Sim"), labels = c("no", "yes")),
    color = method,
    x = method,
    fill = method
  )
) +
  geom_bar(position = "fill") +
  bri_theme +
  labs(title = "Could the deviations have been prevented?") +
  ylab("Proportion of answers") +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = bri_3color) +
  scale_fill_manual(values = bri_3color)

bri_ggsave(
  paste0(output_path, "/self-assessment/", "mudancas_preveniveis.png"),
  plot = plot_preventable_deviations
)

#---
## Table - Deviations prevented ----
#---

tbl_dev_prev <- reshape2::dcast(
  form1_val |>
    filter(!is.na(preventable)) |>
    group_by(preventable) |>
    count(method),
  formula = method ~ preventable
) |>
  dplyr::rename(Method = method) |>
  dplyr::rename(No = Nao) |>
  dplyr::rename(Yes = Sim) |>
  flextable() |>
  add_name("Could the deviations have been prevented?")

save_tbl(
  list(tbl_dev_prev),
  paste0(output_path, "/self-assessment/Table - Deviation Prevented.docx")
)

#---
## Plot 6 - Would you change anything in the protocol for a replication? ----
#---

plot_would_change_protocol <- ggplot(
  form1_val |> filter(!is.na(would_change)),
  aes(
    alpha = factor(would_change, levels = c("Nao", "Sim"), labels = c("no", "yes")),
    color = method,
    x = method,
    fill = method
  )
) +
  geom_bar(position = "fill") +
  bri_theme +
  labs(title = "If you were to replicate this experiment again, would you change anything in the protocol?") +
  ylab("Proportion of answers") +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = bri_3color) +
  scale_fill_manual(values = bri_3color)

bri_ggsave(
  paste0(output_path, "/self-assessment/", "mudaria_protocolo.png"),
  plot = plot_would_change_protocol
)

# Part 2 - Personal experience ------------------------------------------------------

form2 <- read_excel("./other-data/self-assessment/Avaliacao Geral Individual (Responses) - manual edit.xlsx")

### About the respondents

form2 <- form2 |> mutate(position_EN = if_else(condition = `Qual o seu cargo/função no laboratório durante a realização do projeto?` == "Professor / Pesquisador", true = "Professor/Researcher", false = if_else(
  condition = `Qual o seu cargo/função no laboratório durante a realização do projeto?` == "Supervisora do Laboratorio de Toxicologia-Paola Alejandra Cappelletti", true = "Professor/Researcher", false = if_else(
    condition = `Qual o seu cargo/função no laboratório durante a realização do projeto?` == "Doutorando", true = "PhD Student", false = if_else(
      condition = `Qual o seu cargo/função no laboratório durante a realização do projeto?` == "Pós-doutorando", true = "Post-doctoral researcher", false = if_else(
        condition = `Qual o seu cargo/função no laboratório durante a realização do projeto?` == "Mestrando", true = "MSc Student", false = if_else(
          condition = `Qual o seu cargo/função no laboratório durante a realização do projeto?` == "Técnico de Laboratório", true = "Lab technician/technologist", false = if_else(
            condition = `Qual o seu cargo/função no laboratório durante a realização do projeto?` == "Aluno de Iniciação Científica", true = "Undergraduate student", false = if_else(
              condition = `Qual o seu cargo/função no laboratório durante a realização do projeto?` == "Analista", true = "Lab technician/technologist", false = if_else(
                condition = `Qual o seu cargo/função no laboratório durante a realização do projeto?` == "Biomédica (Técnica nível superior)", true = "Lab technician/technologist", false = "ERRO"
              )
            )
          )
        )
      )
    )
  )
)))

#---
## Plot 6 - Plotting respondents' career levels ----
#---

career_plot <- ggplot(form2, aes(y = fct_rev(fct_infreq(position_EN)), fill = fct_rev(fct_infreq(position_EN)))) +
  geom_bar(stat = "count") +
  bri_theme +
  ylab("Career level") +
  scale_x_continuous(n.breaks = 10) +
  theme(legend.position = "none") +
  scale_fill_manual(values = bri_6color)

# Saving the plot
bri_ggsave(paste0(output_path, "/self-assessment/", "career_level.png"), plot = career_plot)


#---
## Roles ----
#---

table_roles_EN <- data.frame(
  "role" = c("Development of protocols", "Execution of experiments", "Data analysis", "Filling in forms", "Communication with \nthe organizing committee"),
  "count" = c(
    length(grep("Elaboração de protocolos", form2[[8]])),
    length(grep("Realização de experimentos", form2[[8]])),
    length(grep("Análise de dados", form2[[8]])),
    length(grep("Preenchimento de formulários", form2[[8]])),
    length(grep("Comunicação com a comissão organizadora", form2[[8]]))
  )
)

ggplot(table_roles_EN, aes(y = reorder(role, count), x = count)) +
  geom_col() +
  bri_theme +
  ylab("Role(s) in the project") +
  labs(caption = "Respondents could select multiple options. Total number of respondents = 121.")



table_roles_EN_Prof <- data.frame(
  "role" = c("Development of protocols", "Execution of experiments", "Data analysis", "Filling in forms", "Communication with \nthe organizing committee"),
  "count" = c(
    length(grep("Elaboração de protocolos", (form2 |> filter(position_EN == "Professor/Researcher"))[[8]])),
    length(grep("Realização de experimentos", (form2 |> filter(position_EN == "Professor/Researcher"))[[8]])),
    length(grep("Análise de dados", (form2 |> filter(position_EN == "Professor/Researcher"))[[8]])),
    length(grep("Preenchimento de formulários", (form2 |> filter(position_EN == "Professor/Researcher"))[[8]])),
    length(grep("Comunicação com a comissão organizadora", (form2 |> filter(position_EN == "Professor/Researcher"))[[8]]))
  ),
  "position" = "Professor/Researcher"
)

table_roles_EN_PhD <- data.frame(
  "role" = c("Development of protocols", "Execution of experiments", "Data analysis", "Filling in forms", "Communication with \nthe organizing committee"),
  "count" = c(
    length(grep("Elaboração de protocolos", (form2 |> filter(position_EN == "PhD Student"))[[8]])),
    length(grep("Realização de experimentos", (form2 |> filter(position_EN == "PhD Student"))[[8]])),
    length(grep("Análise de dados", (form2 |> filter(position_EN == "PhD Student"))[[8]])),
    length(grep("Preenchimento de formulários", (form2 |> filter(position_EN == "PhD Student"))[[8]])),
    length(grep("Comunicação com a comissão organizadora", (form2 |> filter(position_EN == "PhD Student"))[[8]]))
  ),
  "position" = "PhD Student"
)

table_roles_EN_tec <- data.frame(
  "role" = c("Development of protocols", "Execution of experiments", "Data analysis", "Filling in forms", "Communication with \nthe organizing committee"),
  "count" = c(
    length(grep("Elaboração de protocolos", (form2 |> filter(position_EN == "Lab technician/technologist"))[[8]])),
    length(grep("Realização de experimentos", (form2 |> filter(position_EN == "Lab technician/technologist"))[[8]])),
    length(grep("Análise de dados", (form2 |> filter(position_EN == "Lab technician/technologist"))[[8]])),
    length(grep("Preenchimento de formulários", (form2 |> filter(position_EN == "Lab technician/technologist"))[[8]])),
    length(grep("Comunicação com a comissão organizadora", (form2 |> filter(position_EN == "Lab technician/technologist"))[[8]]))
  ),
  "position" = "Lab technician/technologist"
)

table_roles_EN_pd <- data.frame(
  "role" = c("Development of protocols", "Execution of experiments", "Data analysis", "Filling in forms", "Communication with \nthe organizing committee"),
  "count" = c(
    length(grep("Elaboração de protocolos", (form2 |> filter(position_EN == "Post-doctoral researcher"))[[8]])),
    length(grep("Realização de experimentos", (form2 |> filter(position_EN == "Post-doctoral researcher"))[[8]])),
    length(grep("Análise de dados", (form2 |> filter(position_EN == "Post-doctoral researcher"))[[8]])),
    length(grep("Preenchimento de formulários", (form2 |> filter(position_EN == "Post-doctoral researcher"))[[8]])),
    length(grep("Comunicação com a comissão organizadora", (form2 |> filter(position_EN == "Post-doctoral researcher"))[[8]]))
  ),
  "position" = "Post-doctoral researcher"
)

table_roles_EN_MSc <- data.frame(
  "role" = c("Development of protocols", "Execution of experiments", "Data analysis", "Filling in forms", "Communication with \nthe organizing committee"),
  "count" = c(
    length(grep("Elaboração de protocolos", (form2 |> filter(position_EN == "MSc Student"))[[8]])),
    length(grep("Realização de experimentos", (form2 |> filter(position_EN == "MSc Student"))[[8]])),
    length(grep("Análise de dados", (form2 |> filter(position_EN == "MSc Student"))[[8]])),
    length(grep("Preenchimento de formulários", (form2 |> filter(position_EN == "MSc Student"))[[8]])),
    length(grep("Comunicação com a comissão organizadora", (form2 |> filter(position_EN == "MSc Student"))[[8]]))
  ),
  "position" = "MSc Student"
)

table_roles_EN_und <- data.frame(
  "role" = c("Development of protocols", "Execution of experiments", "Data analysis", "Filling in forms", "Communication with \nthe organizing committee"),
  "count" = c(
    length(grep("Elaboração de protocolos", (form2 |> filter(position_EN == "Undergraduate student"))[[8]])),
    length(grep("Realização de experimentos", (form2 |> filter(position_EN == "Undergraduate student"))[[8]])),
    length(grep("Análise de dados", (form2 |> filter(position_EN == "Undergraduate student"))[[8]])),
    length(grep("Preenchimento de formulários", (form2 |> filter(position_EN == "Undergraduate student"))[[8]])),
    length(grep("Comunicação com a comissão organizadora", (form2 |> filter(position_EN == "Undergraduate student"))[[8]]))
  ),
  "position" = "Undergraduate student"
)

table_roles_EN <- rbind(table_roles_EN_Prof, table_roles_EN_PhD, table_roles_EN_MSc, table_roles_EN_tec, table_roles_EN_und, table_roles_EN_pd)

table_roles_EN$role <- factor(table_roles_EN$role, levels = rev(c("Data analysis", "Filling in forms", "Execution of experiments", "Development of protocols", "Communication with \nthe organizing committee")))

table_roles_EN$position <- factor(table_roles_EN$position, levels = rev(c("Professor/Researcher", "PhD Student", "Lab technician/technologist", "Post-doctoral researcher", "MSc Student", "Undergraduate student")))


#---
### Plot 7 - Plotting respondents' roles in the project ----
#---

roles_plot <- ggplot(table_roles_EN, aes(y = role, x = count, fill = position)) +
  geom_col() +
  bri_theme +
  ylab("Role(s) in the project") +
  labs(caption = "Respondents could select multiple options. Total number of respondents = 121.") +
  scale_fill_manual(values = bri_6color) +
  theme(legend.position = "bottom", legend.title = element_blank())

# Saving the plot
bri_ggsave(paste0(output_path, "/self-assessment/", "career_roles.png"), plot = roles_plot)

table_activities_EN <- data.frame(
  "role" = c("General meetings", "Webinars", "Lab meetings with \nthe orgnizing committee", "Validation committee", "Spin-off projects"),
  "count" = c(
    length(grep("Reuniões gerais", form2[[10]])),
    length(grep("Webinários", form2[[10]])),
    length(grep("Reuniões do laboratório com a comissão organizadora", form2[[10]])),
    length(grep("Comissão de validação", form2[[10]])),
    length(grep("Projetos spin-off", form2[[10]]))
  )
)

ggplot(table_activities_EN |> filter(count != 0), aes(y = reorder(role, count), x = count)) +
  geom_col() +
  bri_theme +
  ylab("Activities engaged with") +
  labs(caption = "Respondents could select multiple options. Total number of respondents = 121.")



table_activities_EN_Prof <- data.frame(
  "role" = c("General meetings", "Webinars", "Lab meetings with \nthe orgnizing committee", "Validation committee", "Spin-off projects"),
  "count" = c(
    length(grep("Reuniões gerais", (form2 |> filter(position_EN == "Professor/Researcher"))[[10]])),
    length(grep("Webinários", (form2 |> filter(position_EN == "Professor/Researcher"))[[10]])),
    length(grep("Reuniões do laboratório com a comissão organizadora", (form2 |> filter(position_EN == "Professor/Researcher"))[[10]])),
    length(grep("Comissão de validação", (form2 |> filter(position_EN == "Professor/Researcher"))[[10]])),
    length(grep("Projetos spin-off", (form2 |> filter(position_EN == "Professor/Researcher"))[[10]]))
  ),
  "position" = "Professor/Researcher"
)

table_activities_EN_PhD <- data.frame(
  "role" = c("General meetings", "Webinars", "Lab meetings with \nthe orgnizing committee", "Validation committee", "Spin-off projects"),
  "count" = c(
    length(grep("Reuniões gerais", (form2 |> filter(position_EN == "PhD Student"))[[10]])),
    length(grep("Webinários", (form2 |> filter(position_EN == "PhD Student"))[[10]])),
    length(grep("Reuniões do laboratório com a comissão organizadora", (form2 |> filter(position_EN == "PhD Student"))[[10]])),
    length(grep("Comissão de validação", (form2 |> filter(position_EN == "PhD Student"))[[10]])),
    length(grep("Projetos spin-off", (form2 |> filter(position_EN == "PhD Student"))[[10]]))
  ),
  "position" = "PhD Student"
)

table_activities_EN_tec <- data.frame(
  "role" = c("General meetings", "Webinars", "Lab meetings with \nthe orgnizing committee", "Validation committee", "Spin-off projects"),
  "count" = c(
    length(grep("Reuniões gerais", (form2 |> filter(position_EN == "Lab technician/technologist"))[[10]])),
    length(grep("Webinários", (form2 |> filter(position_EN == "Lab technician/technologist"))[[10]])),
    length(grep("Reuniões do laboratório com a comissão organizadora", (form2 |> filter(position_EN == "Lab technician/technologist"))[[10]])),
    length(grep("Comissão de validação", (form2 |> filter(position_EN == "Lab technician/technologist"))[[10]])),
    length(grep("Projetos spin-off", (form2 |> filter(position_EN == "Lab technician/technologist"))[[10]]))
  ),
  "position" = "Lab technician/technologist"
)

table_activities_EN_pd <- data.frame(
  "role" = c("General meetings", "Webinars", "Lab meetings with \nthe orgnizing committee", "Validation committee", "Spin-off projects"),
  "count" = c(
    length(grep("Reuniões gerais", (form2 |> filter(position_EN == "Post-doctoral researcher"))[[10]])),
    length(grep("Webinários", (form2 |> filter(position_EN == "Post-doctoral researcher"))[[10]])),
    length(grep("Reuniões do laboratório com a comissão organizadora", (form2 |> filter(position_EN == "Post-doctoral researcher"))[[10]])),
    length(grep("Comissão de validação", (form2 |> filter(position_EN == "Post-doctoral researcher"))[[10]])),
    length(grep("Projetos spin-off", (form2 |> filter(position_EN == "Post-doctoral researcher"))[[10]]))
  ),
  "position" = "Post-doctoral researcher"
)

table_activities_EN_MSc <- data.frame(
  "role" = c("General meetings", "Webinars", "Lab meetings with \nthe orgnizing committee", "Validation committee", "Spin-off projects"),
  "count" = c(
    length(grep("Reuniões gerais", (form2 |> filter(position_EN == "MSc Student"))[[10]])),
    length(grep("Webinários", (form2 |> filter(position_EN == "MSc Student"))[[10]])),
    length(grep("Reuniões do laboratório com a comissão organizadora", (form2 |> filter(position_EN == "MSc Student"))[[10]])),
    length(grep("Comissão de validação", (form2 |> filter(position_EN == "MSc Student"))[[10]])),
    length(grep("Projetos spin-off", (form2 |> filter(position_EN == "MSc Student"))[[10]]))
  ),
  "position" = "MSc Student"
)

table_activities_EN_und <- data.frame(
  "role" = c("General meetings", "Webinars", "Lab meetings with \nthe orgnizing committee", "Validation committee", "Spin-off projects"),
  "count" = c(
    length(grep("Reuniões gerais", (form2 |> filter(position_EN == "Undergraduate student"))[[10]])),
    length(grep("Webinários", (form2 |> filter(position_EN == "Undergraduate student"))[[10]])),
    length(grep("Reuniões do laboratório com a comissão organizadora", (form2 |> filter(position_EN == "Undergraduate student"))[[10]])),
    length(grep("Comissão de validação", (form2 |> filter(position_EN == "Undergraduate student"))[[10]])),
    length(grep("Projetos spin-off", (form2 |> filter(position_EN == "Undergraduate student"))[[10]]))
  ),
  "position" = "Undergraduate student"
)

table_activities_EN <- rbind(table_activities_EN_Prof, table_activities_EN_PhD, table_activities_EN_MSc, table_activities_EN_tec, table_activities_EN_und, table_activities_EN_pd)

table_activities_EN$role <- factor(table_activities_EN$role, levels = rev(c("General meetings", "Webinars", "Lab meetings with \nthe orgnizing committee", "Spin-off projects", "Validation committee")))

table_activities_EN$position <- factor(table_activities_EN$position, levels = rev(c("Professor/Researcher", "PhD Student", "Lab technician/technologist", "Post-doctoral researcher", "MSc Student", "Undergraduate student")))

ggplot(table_activities_EN, aes(y = role, x = count, fill = position)) +
  geom_col() +
  bri_theme +
  ylab("Activities in the project") +
  labs(caption = "Respondents could select multiple options. Total number of respondents = 121.") +
  scale_fill_manual(values = bri_6color) +
  theme(legend.position = "bottom", legend.title = element_blank())

#---
### Plot 8 - Plotting activities engaged by respondents in the project ----
#---

activities_plot <- ggplot(table_activities_EN |> filter(count != 0), aes(y = role, x = count, fill = position)) +
  geom_col() +
  bri_theme +
  ylab("Activities in the project") +
  labs(caption = "Respondents could select multiple options. Total number of respondents = 121.") +
  scale_fill_manual(values = bri_6color) +
  theme(legend.position = "bottom", legend.title = element_blank())

# Saving the plot
bri_ggsave(paste0(output_path, "/self-assessment/", "career_activities_communication.png"), plot = activities_plot)


#---
## Learnings ----
#---

get_n_matches <- function(data, question, option) {
  option_data <- data
  length(grep(option, option_data[[question]]))
}

learning_options <- c(
  "1a", "1b", "1c", "1d", "1e",
  "2a", "2b", "2c", "2d",
  "3a", "3b", "3c", "3d",
  "4a", "4b", "4c",
  "5a", "5b", "5c", "5d",
  "6a", "6b", "6c", "6d", "6e", "6f",
  "7a", "7b", "7c", "7d", "7e", "7f",
  "8", "9"
)

cat_learning <- sapply(data = form2, question = "Aprendizados-cat", X = learning_options, FUN = get_n_matches) |>
  as.data.frame() |>
  rename("category_counts" = c(1)) |>
  rownames_to_column("learning.options")

learning.labels <- c(
  "1a" = "(1a) Protocol development",
  "1b" = "(1b) Protocol pre-registration",
  "1c" = "(1c) Protocol use",
  "1d" = "(1d) Organization of experimental steps",
  "1e" = "(1e) Experimental planning (general)",
  "2a" = "(2a) Use of controls",
  "2b" = "(2b) Use of validity criteria",
  "2c" = "(2c) Bias control measures",
  "2d" = "(2d) Experimental design (general)",
  "3a" = "(3a) Cell culture",
  "3b" = "(3b) Experimental animals",
  "3c" = "(3c) Specific technique",
  "3d" = "(3d) Experimental technique (general)",
  "4a" = "(4a) Recording information during experiment",
  "4b" = "(4b) Data management",
  "4c" = "(4c) Documentation (general)",
  "5a" = "(5a) Statistical analysis",
  "5b" = "(5b) Results interpretation",
  "5c" = "(5c) Reproducibility criteria",
  "5d" = "(5d) Data analysis (general)",
  "6a" = "(6a) Reproducibility importance and concepts",
  "6b" = "(6b) Reproducibility problems",
  "6c" = "(6c) Importance and issues with reporting",
  "6d" = "(6d) Open science",
  "6e" = "(6e) Research integrity and ethics",
  "6f" = "(6f) View of the scientific process (general)",
  "7a" = "(7a) Personal organization",
  "7b" = "(7b) Management",
  "7c" = "(7c) Collaborative work in the lab",
  "7d" = "(7d) Collaborative work in multicentric project",
  "7e" = "(7e) Personal skills (general)",
  "7f" = "(7*) Collaborative work (not specified)",
  "8" = "Other",
  "9" = "No answer"
)

#---
### Plot 9 - Plotting main learnings from respondents ----
#---

learning_plot <- ggplot(cat_learning, aes(y = learning_options, x = category_counts)) +
  geom_bar(stat = "identity") +
  bri_theme +
  labs(y = "What were your main learnings?", x = "Number of answers") +
  scale_y_discrete(labels = learning.labels, limits = rev(c(
    "1a", "1b", "1c", "1d", "1e", " ",
    "2a", "2b", "2c", "2d", " ",
    "3a", "3b", "3c", "3d", " ",
    "4a", "4b", "4c", " ",
    "5a", "5b", "5c", "5d", " ",
    "6a", "6b", "6c", "6d", "6e", "6f", " ",
    "7a", "7b", "7c", "7d", "7e", "7f", " ",
    "8", "9"
  )))

# Saving the plot
bri_ggsave(paste0(output_path, "/self-assessment/", "learning_complete.png"), plot = learning_plot)

#---
### Plot 10 - Plotting main learnings from respondents ----
#---

learning_plot <- ggplot(cat_learning, aes(y = learning_options, x = category_counts)) +
  geom_bar(stat = "identity") +
  bri_theme +
  labs(y = "What were your main learnings?", x = "Number of answers") +
  scale_y_discrete(labels = learning.labels, limits = rev(c(
    "1a", "1b", "1c", "1d", "1e", " ",
    "2b", "2c", "2d", " ",
    "3a", "3b", "3c", "3d", " ",
    "4a", "4b", "4c", " ",
    "5b", "5d", " ",
    "6a", "6b", "6c", "6d", "6e", "6f", " ",
    "7a", "7b", "7c", "7d", "7e", "7f", " ",
    "8"
  ))) +
  labs(caption = "Responses could include multiple categories. Total number of answers to this form = 121, of which 101 provided responses to this question.")

# Saving the plot
bri_ggsave(paste0(output_path, "/self-assessment/", "learning.png"), plot = learning_plot)


#---
### Plot 11: General order of learning options  ----
#---

learning.labels.fig <- c(
  "1a" = "Protocol development",
  "1b" = "Protocol pre-registration",
  "1c" = "Protocol use",
  "1d" = "Organization of experimental steps",
  "1e" = "Experimental planning (general)",
  "2a" = "Use of controls",
  "2b" = "Use of validity criteria",
  "2c" = "Bias control measures",
  "2d" = "Experimental design (general)",
  "3a" = "Cell culture",
  "3b" = "Experimental animals",
  "3c" = "Specific technique",
  "3d" = "Experimental technique (general)",
  "4a" = "Recording information during experiment",
  "4b" = "Data management",
  "4c" = "Documentation (general)",
  "5a" = "Statistical analysis",
  "5b" = "Results interpretation",
  "5c" = "Reproducibility criteria",
  "5d" = "Data analysis (general)",
  "6a" = "Reproducibility importance and concepts",
  "6b" = "Reproducibility problems",
  "6c" = "Importance and issues with reporting",
  "6d" = "Open science",
  "6e" = "Research integrity and ethics",
  "6f" = "View of the scientific process (general)",
  "7a" = "Personal organization",
  "7b" = "Management",
  "7c" = "Collaborative work in the lab",
  "7d" = "Collaborative work in multicentric project",
  "7e" = "Personal skills (general)",
  "7f" = "Collaborative work (not specified)",
  "8" = "Other",
  "9" = "No answer"
)

plot_learning_order_general <- ggplot(cat_learning, aes(y = learning_options, x = category_counts)) +
  geom_bar(stat = "identity", width = 0.7) +
  bri_theme +
  labs(y = "What were your main learnings?", x = "Number of answers") +
  scale_y_discrete(labels = learning.labels.fig, limits = rev(c(
    "6c", "6a", "6b", "1a", "1e", "7d",
    "4a", "3c", "7c", "7e", "1c", "6f", "1d",
    "6d", "3a", "3d", "4c",
    "4b", "7b", "1b", "2c", "2d",
    "5b", "5d", "6e", "7a",
    "3b", "7f", "2b", "8"
  ))) +
  labs(caption = "Responses could include multiple categories. Total number of answers to this form = 121, of which 101 provided responses to this question.")

# Save Plot 11
bri_ggsave(paste0(output_path, "/self-assessment/", "learning-order-general.png"), plot = plot_learning_order_general, height = 15, width = 20, units = "cm")

#---
### Plot 12: Specific order within categories  ----
#---

plot_learning_order_within <- ggplot(cat_learning, aes(y = learning_options, x = category_counts)) +
  geom_bar(stat = "identity") +
  bri_theme +
  labs(y = "What were your main learnings?", x = "Number of answers") +
  scale_y_discrete(labels = learning.labels.fig, limits = rev(c(
    "1a", "1e", "1c", "1d", "1b", " ",
    "2c", "2d", "2b", " ",
    "3c", "3a", "3d", "3b", " ",
    "4a", "4c", "4b", " ",
    "5b", "5d", " ",
    "6c", "6a", "6b", "6f", "6d", "6e", " ",
    "7d", "7e", "7c", "7b", "7a", "7f", " ",
    "8"
  ))) +
  labs(caption = "Responses could include multiple categories. Total number of answers to this form = 121, of which 101 provided responses to this question.")

# Save Plot 12
bri_ggsave(paste0(output_path, "/self-assessment/", "learning-order-within.png"), plot = plot_learning_order_within, height = 15, width = 20, units = "cm")

#---
### Plot 13: Combination of order within and outside categories  ----
#---

plot_learning_order_combined <- ggplot(cat_learning, aes(y = learning_options, x = category_counts)) +
  geom_bar(stat = "identity") +
  bri_theme +
  labs(y = "What were your main learnings?", x = "Number of answers") +
  scale_y_discrete(labels = learning.labels.fig, limits = rev(c(
    "6c", "6a", "6b", "6f", "6d", "6e", " ",
    "1a", "1e", "1c", "1d", "1b", " ",
    "7d", "7e", "7c", "7b", "7a", "7f", " ",
    "4a", "4c", "4b", " ",
    "3c", "3a", "3d", "3b", " ",
    "2c", "2d", "2b", " ",
    "5b", "5d", " ",
    "8"
  ))) +
  labs(caption = "Responses could include multiple categories. Total number of answers to this form = 121, of which 101 provided responses to this question.")

# Save Plot 13
bri_ggsave(paste0(output_path, "/self-assessment/", "learning-order-combined.png"), plot = plot_learning_order_combined, height = 15, width = 20, units = "cm")


#---
## Difficulties ----
#---

dificulties_options <- c(
  "1a", "1b", "1c",
  "2a", "2b", "2c", "2d",
  "3a", "3b", "3c", "3d",
  "4a", "4b", "4c", "4d", "4e", "4f", "4g", "4h",
  "5a", "5b", "5c", "5d",
  "6a", "6b", "6c", "6d", "6e", "6f", "6g",
  "7a", "7b", "7c", "7d", "7e", "7f", "7g",
  "8", "9", "10"
)

cat_dificulties <- sapply(data = form2, question = "Dificuldades-cat", X = dificulties_options, FUN = get_n_matches) |>
  as.data.frame() |>
  rename("category_counts" = c(1)) |>
  rownames_to_column("dificulties.options")



dificulties.labels <- c(
  "1a" = "(1a) CEUAs",
  "1b" = "(1b) Pandemic",
  "1c" = "(1c) External factors (general)",
  "2a" = "(2a) Lab changed location",
  "2b" = "(2b) Equipment availability",
  "2c" = "(2c) Animal facility availability",
  "2d" = "(2d) Infrastructure (general)",
  "3a" = "(3a) Learning new methods",
  "3b" = "(3b) Materials/reagents",
  "3c" = "(3c) Cells/cell culture/animals",
  "3d" = "(3d) Experimental dificulties (general)",
  "4a" = "(4a) Original protocol missing information/clarity",
  "4b" = "(4b) Original protocol inadequate methods",
  "4c" = "(4c) Developing a protocol",
  "4d" = "(4d) Replication protocol missing information/clarity",
  "4e" = "(4e) Required adaptations",
  "4f" = "(4f) Missing information/clarity (unclear protocol)",
  "4g" = "(4g) Protocol dificulties (general)",
  "4h" = "(4*) Terminoogy",
  "5a" = "(5a) Communication with organizing committee",
  "5b" = "(5b) Missing expertise in organizing committee",
  "5c" = "(5c) Data management",
  "5d" = "(5d) Multicentric project (general)",
  "6a" = "(6a) Missing materials/need for purchase",
  "6b" = "(6b) Delayed delivery",
  "6c" = "(6c) Purchase process",
  "6d" = "(6d) Expired/expiring material on arrival",
  "6e" = "(6e) Mistakes in purchases/deliveries",
  "6f" = "(6f) Licenses",
  "6g" = "(6g) Purchases (general)",
  "7a" = "(7a) Missing personnel",
  "7b" = "(7b) Changes in personnel",
  "7c" = "(7c) Low engagement",
  "7d" = "(7d) Fitting in with other activities",
  "7e" = "(7e) High time demand",
  "7f" = "(7f) Lab's team (general)",
  "7g" = "(7*) Deadlines",
  "8" = "Other",
  "9" = "No answer",
  "10" = "Not a dificulty"
)

#---
### Plot 13: General order of difficulties faced ----
#---

plot_difficulties_general <- ggplot(cat_dificulties, aes(y = dificulties_options, x = category_counts)) +
  geom_bar(stat = "identity") +
  bri_theme +
  labs(y = "What were the main difficulties you faced?", x = "Number of answers") +
  scale_y_discrete(labels = dificulties.labels, limits = rev(c(
    "1a", "1b", "1c", " ",
    "2a", "2b", "2c", "2d", " ",
    "3a", "3b", "3c", "3d", " ",
    "4a", "4b", "4c", "4d", "4e", "4f", "4g", "4h", " ",
    "5a", "5b", "5c", "5d", " ",
    "6a", "6b", "6c", "6d", "6e", "6f", "6g", " ",
    "7a", "7b", "7c", "7d", "7e", "7f", "7g", " ",
    "8", "9", "10"
  ))) +
  labs(caption = "Responses could include multiple categories. Total number of answers to this form = 121, of which 103 provided responses to this question.")

# Save Plot 13
bri_ggsave(paste0(output_path, "/self-assessment/", "difficulties-general.png"), plot = plot_difficulties_general)

#---
### Plot 14: Specific difficulties faced ----
#---

plot_difficulties_specific <- ggplot(cat_dificulties, aes(y = dificulties_options, x = category_counts)) +
  geom_bar(stat = "identity") +
  bri_theme +
  labs(y = "What were the main difficulties you faced?", x = "Number of answers") +
  scale_y_discrete(labels = dificulties.labels, limits = rev(c(
    "1a", "1b", " ",
    "2a", "2b", "2c", "2d", " ",
    "3a", "3b", "3c", "3d", " ",
    "4a", "4b", "4c", "4d", "4e", "4f", "4g", "4h", " ",
    "5a", "5b", "5c", "5d", " ",
    "6a", "6b", "6c", "6d", "6e", "6g", " ",
    "7a", "7b", "7c", "7d", "7e", "7f", "7g", " ",
    "8"
  ))) +
  labs(caption = "Responses could include multiple categories. Total number of answers to this form = 121, of which 103 provided responses to this question.")

# Save Plot 14
bri_ggsave(paste0(output_path, "/self-assessment/", "difficulties-specific.png"), plot = plot_difficulties_specific)


difficulties.labels.fig <- c(
  # External factors
  "1a" = "Animal ethics committees",
  "1b" = "Pandemic",
  "1c" = "External factors (general)",
  # Infrastrucure difficulties
  "2a" = "Relocation of facilities",
  "2b" = "Equipment availability",
  "2c" = "Animal facility limitations",
  "2d" = "Infrastructure (general)",
  # Experimental difficulties
  "3a" = "Learning new methods",
  "3b" = "Materials/reagents",
  "3c" = "Cells/cell culture/animals",
  "3d" = "Experimental difficulties (general)",
  # Protocols
  "4a" = "Unclear original protocol",
  "4b" = "Inadequate original protocol",
  "4c" = "Protocol development",
  "4d" = "Unclear replication protocol",
  "4e" = "Required adaptations",
  "4f" = "Unclear protocol (not specified)",
  "4g" = "Protocols (general)",
  "4h" = "Terminology problems",
  # Multicentre project
  "5a" = "Communication with organizing committee",
  "5b" = "Lack of expertise of coordinating team",
  "5c" = "Data management",
  "5d" = "Multicentric project (general)",
  # Reagent acquisition
  "6a" = "Lack of materials/need for acquisition",
  "6b" = "Delayed delivery",
  "6c" = "Acquisition process",
  "6d" = "Delivery close to expiration",
  "6e" = "Errors in acquisition/delivery",
  "6f" = "Licenses",
  "6g" = "Reagent acquisition (general)",
  # Replication team
  "7a" = "Lack of personnel",
  "7b" = "Changes in personnel",
  "7c" = "Low engagement",
  "7d" = "Balance with other activities",
  "7e" = "High workload",
  "7f" = "Lab's team (general)",
  "7g" = "Deadlines",
  #
  "8" = "Other",
  #
  "9" = "No answer",
  #
  "10" = "Not a difficulty"
)

#---
### Plot 14: General order of difficulties faced ----
#---

plot_difficulties_order_general <- ggplot(cat_dificulties, aes(y = dificulties_options, x = category_counts)) +
  geom_bar(stat = "identity", width = 0.7) +
  bri_theme +
  labs(y = "What were the main difficulties you faced?", x = "Number of answers") +
  scale_y_discrete(labels = difficulties.labels.fig, limits = rev(c(
    "1b", "6b", "3c", "7d", "7b",
    "4f", "7a", "4c", "7f",
    "2c", "2d", "3b", "3d", "4e",
    "6d", "2a", "3a", "6a", "6g", "7g",
    "2b", "4b", "4g", "5a",
    "4a", "5b", "5c", "7c", "7e",
    "1a", "4d", "4h", "5d", "6c", "6e",
    "", "8"
  ))) +
  labs(caption = "Responses could include multiple categories. Total number of answers to this form = 121, of which 103 provided responses to this question.")

# Save Plot 14
bri_ggsave(paste0(output_path, "/self-assessment/", "difficulties-order-general.png"), plot = plot_difficulties_order_general, width = 20, height = 15, units = "cm")

#---
### Plot 15: Specific order of difficulties faced ----
#---

plot_difficulties_order_specific <- ggplot(cat_dificulties, aes(y = dificulties_options, x = category_counts)) +
  geom_bar(stat = "identity") +
  bri_theme +
  labs(y = "What were the main difficulties you faced?", x = "Number of answers") +
  scale_y_discrete(labels = difficulties.labels.fig, limits = rev(c(
    "1b", "1a", " ",
    "2c", "2d", "2a", "2b", " ",
    "3c", "3b", "3d", "3a", " ",
    "4f", "4c", "4e", "4b", "4g", "4a", "4d", "4h", " ",
    "5a", "5b", "5c", "5d", " ",
    "6b", "6d", "6a", "6g", "6c", "6e", " ",
    "7d", "7b", "7a", "7f", "7g", "7c", "7e", " ",
    "8"
  ))) +
  labs(caption = "Responses could include multiple categories. Total number of answers to this form = 121, of which 103 provided responses to this question.")

# Save Plot 15
bri_ggsave(paste0(output_path, "/self-assessment/", "difficulties-order-specific.png"), plot = plot_difficulties_order_specific, width = 20, height = 15, units = "cm")

#---
### Plot 16: Combined order of difficulties faced ----
#---

plot_difficulties_order_combined <- ggplot(cat_dificulties, aes(y = dificulties_options, x = category_counts)) +
  geom_bar(stat = "identity") +
  bri_theme +
  labs(y = "What were the main difficulties you faced?", x = "Number of answers") +
  scale_y_discrete(labels = difficulties.labels.fig, limits = rev(c(
    "1b", "1a", " ",
    "6b", "6d", "6a", "6g", "6c", "6e", " ",
    "3c", "3b", "3d", "3a", " ",
    "7d", "7b", "7a", "7f", "7g", "7c", "7e", " ",
    "4f", "4c", "4e", "4b", "4g", "4a", "4d", "4h", " ",
    "2c", "2d", "2a", "2b", " ",
    "5a", "5b", "5c", "5d", " ",
    "8"
  ))) +
  labs(caption = "Responses could include multiple categories. Total number of answers to this form = 121, of which 103 provided responses to this question.")

# Save Plot 16
bri_ggsave(paste0(output_path, "/self-assessment/", "difficulties-order-combined.png"), plot = plot_difficulties_order_combined, width = 20, height = 15, units = "cm")

### Figure S9 ---------------------------------------------------------------

difficulties.labels.fig[["8"]] <- "**Other**" 

y_levels_s9 <- rev(c(
  "**Replication team**",
  "7d", "7b", "7a", "7g", "7c", "7e", "7f",
  "**Protocols**",
  "4f", "4c", "4e", "4b", "4a", "4d", "4h", "4g",
  "**Reagent acquisition**",
  "6b", "6d", "6a", "6c", "6e", "6g",
  "**Infrastructure difficulties**",
  "2c", "2a", "2b", "2d",
  "**Experimental difficulties**",
  "3c", "3b", "3a", "3d",
  "**External factors**",
  "1b", "1a",
  "**Multicentre project**",
  "5a", "5b", "5c", "5d",
  " ",
  "8"
))

grid_positions_s9 <- setdiff(
  seq_along(y_levels_s9),
  which(y_levels_s9 %in% c(
    "**Replication team**",
    "**Protocols**",
    "**Reagent acquisition**",
    "**Infrastructure difficulties**",
    "**Experimental difficulties**",
    "**External factors**",
    "**Multicentre project**",
    " "
  ))
)

fig_s9 <- ggplot(cat_dificulties, aes(y = dificulties_options, x = category_counts)) +
  geom_hline(yintercept = grid_positions_s9, color = "grey90") +
  geom_bar(stat = "identity") +
  # geom_text(aes(label = category_counts, x = 3)) +

  bri_theme +
  labs(y = "What were the main difficulties you faced?", x = "Number of answers") +
  scale_y_discrete(labels = difficulties.labels.fig, limits = y_levels_s9, expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(expand = FALSE, clip = "off") +
  theme(
    axis.text.y = ggtext::element_markdown(), # Increase line height to add space between labels
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.margin = unit(c(5, 5, 5, 5), "mm"),
    axis.ticks.y = element_line(color = c(
      "black", "white", "black", "black", "black",
      "black", "white", "black", "black", "white",
      "black", "black", "black", "black", "white",
      "black", "black", "black", "black", "white",
      "black", "black", "black", "black", "black",
      "black", "white", "black", "black", "black",
      "black", "black", "black", "black", "black",
      "white", "black", "black", "black", "black",
      "black", "black", "black", "white"
    ))
  )

difficulties.labels.fig[["8"]] <- "Other" 

bri_ggsave(paste0(output_path, "/self-assessment/", "Figure S9.png"), plot = fig_s9, width = 20, height = 25, units = "cm")

## Table S24 ---------------------------------------------------------------


# Vector with the desired order. Note that strings that begin and end with "**"
# correspond to header lines and should not be joined with the original data.
order_vec <- c(
  "Replication team",
  "7d", "7b", "7a", "7g", "7c", "7e", "7f",
  "Protocols",
  "4f", "4c", "4e", "4b", "4a", "4d", "4h", "4g",
  "Reagent acquisition",
  "6b", "6d", "6a", "6c", "6e", "6g",
  "Infrastructure difficulties",
  "2c", "2a", "2b", "2d",
  "Experimental difficulties",
  "3c", "3b", "3a", "3d",
  "External factors",
  "1b", "1a",
  "Multicentre project",
  "5a", "5b", "5c", "5d",
  "8"
)

# Create a tibble that defines the order.
# For rows that are headers (i.e., those that start and end with "**")
# or the blank line (" "), the 'code' column is set to NA.
order_tbl <- tibble(
  order = seq_along(order_vec),
  label = order_vec,
  code = if_else(str_detect(order_vec, "^\\*\\*.*\\*\\*$") | order_vec == " ",
    NA_character_,
    order_vec
  )
)

# To facilitate the join, rename the code column in the original data frame
# (here, the 'dificulties.options' column contains values like "1a", "7d", etc.)
data_tbl <- cat_dificulties |>
  mutate(labels = difficulties.labels.fig) |>
  rename(code = dificulties.options)

# Join the order tibble with the data using 'code'.
# Header lines will have no match and will have NA in 'category_counts'.
final_tbl <- order_tbl |>
  left_join(data_tbl, by = "code") |>
  # Create the 'count' column from 'category_counts'.
  mutate(count = category_counts) |>
  # Select the final columns; to keep the original order, you may include 'order'.
  select(label, count, labels)

# If you prefer column names in uppercase, as in the original code:
final_tbl <- final_tbl |>
  rename(
    Label = label,
    Count = count
  )

# Replace NA values in 'labels' with 'Label' and finalize the selection.
tbl_s24 <- final_tbl |>
  mutate(labels = if_else(is.na(labels), Label, labels)) |>
  select(labels, Count) |>
  rename(Label = labels) |>
  flextable() |>
  bold(i = c(1, 9, 18, 25, 30, 35, 38, 43)) |> 
  set_table_properties(layout = "autofit") |>
  add_name("Table - Examples of difficulties faced by labs replicating the experiments.") |>
  add_footer_lines("Table shows the various categories and subcategories for general difficulties (in order of general category frequency, as in Figure S3), with illustrative examples of participant responses for each subcategory.")

save_tbl(
  list(tbl_s24),
  paste0("output/", results_path, "/_manuscript figures and tables", "/tables/Table S24.docx")
)
cat("\n### Table S24 generated! ###\n")

## Self-assessments ----

alpha_mapping <- data.frame(n = seq(0, 121), alpha = seq(0, 1, length.out = 122))


self_engagem_counts <- form2 |>
  count(`Qual a sua avaliação sobre o seu engajamento individual com o projeto?`) |>
  rename("Item" = "Qual a sua avaliação sobre o seu engajamento individual com o projeto?")
self_engagem_counts$Item <- as.character(self_engagem_counts$Item)
self_engagem_counts[1, 1] <- "Low"
self_engagem_counts[2, 1] <- "Moderately Low"
self_engagem_counts[3, 1] <- "Intermediate"
self_engagem_counts[4, 1] <- "Moderately High"
self_engagem_counts[5, 1] <- "High"
# self_engagem_counts = self_engagem_counts |> mutate(Item = if_else(condition = is.na(Item), "NA", Item))
self_engagem_counts <- self_engagem_counts |> mutate(fake_scale = 2)

self_engagem_counts <- merge(self_engagem_counts, alpha_mapping, by = "n", all.x = TRUE)

### Plot 17: Engagement with the project ----

plot_engagement <- ggplot(subset(self_engagem_counts, !is.na(Item)), aes(x = factor(Item, levels = c("Low", "Moderately Low", "Intermediate", "Moderately High", "High")), y = fake_scale)) +
  geom_bar(fill = bri_color["main"], alpha = subset(self_engagem_counts, !is.na(Item))$alpha, stat = "identity") +
  geom_text(aes(x = factor(Item, levels = c("Low", "Moderately Low", "Intermediate", "Moderately High", "High")), y = fake_scale / fake_scale, label = n, fontface = "bold"), color = bri_color["dark"], alpha = 1, size = 5) +
  bri_theme +
  labs(
    title = "How do you assess your own engagement with the project?",
    subtitle = "Categories were presented as 1 = Low, 5 = High"
  ) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.title.x = element_blank())

# Save Plot 17
bri_ggsave(paste0(output_path, "/self-assessment/", "engagement.png"), plot = plot_engagement, width = 20, height = 5, units = "cm")

lab_engagem_counts <- form2 |>
  count(`Como você avalia o engajamento do seu laboratório nas tarefas que lhe foram designadas?`) |>
  rename("Item" = "Como você avalia o engajamento do seu laboratório nas tarefas que lhe foram designadas?")
lab_engagem_counts$Item <- as.character(lab_engagem_counts$Item)
lab_engagem_counts[1, 1] <- "Low"
lab_engagem_counts[2, 1] <- "Moderately Low"
lab_engagem_counts[3, 1] <- "Intermediate"
lab_engagem_counts[4, 1] <- "Moderately High"
lab_engagem_counts[5, 1] <- "High"
lab_engagem_counts <- lab_engagem_counts |> mutate(fake_scale = 2)
lab_engagem_counts <- merge(lab_engagem_counts, alpha_mapping, by = "n", all.x = TRUE)

### Plot 18: Lab engagement with assigned tasks ----

plot_lab_engagement <- ggplot(subset(lab_engagem_counts, !is.na(Item)), aes(x = factor(Item, levels = c("Low", "Moderately Low", "Intermediate", "Moderately High", "High")), y = fake_scale, alpha = n)) +
  geom_bar(fill = bri_color["main"], alpha = subset(lab_engagem_counts, !is.na(Item))$alpha, stat = "identity") +
  geom_text(aes(x = factor(Item, levels = c("Low", "Moderately Low", "Intermediate", "Moderately High", "High")), y = fake_scale / fake_scale, label = n, fontface = "bold"), color = bri_color["dark"], alpha = 1, size = 5) +
  bri_theme +
  labs(
    title = "How do you assess your lab's engagement on assigned tasks?",
    subtitle = "Categories were presented as 1 = Low (tasks incomplete or delayed, unfinished protocols), \n5 = High (tasks completed on time, all protocols completed)"
  ) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.title.x = element_blank())

# Save Plot 18
bri_ggsave(paste0(output_path, "/self-assessment/", "lab-engagement.png"), plot = plot_lab_engagement, width = 20, height = 5, units = "cm")



lab_extra_counts <- form2 |>
  count(`Como você avalia o engajamento do seu laboratório no projeto, para além das tarefas específicas que foram designadas?`) |>
  rename("Item" = "Como você avalia o engajamento do seu laboratório no projeto, para além das tarefas específicas que foram designadas?")
lab_extra_counts$Item <- as.character(lab_extra_counts$Item)
lab_extra_counts[1, 1] <- "Low"
lab_extra_counts[2, 1] <- "Moderately Low"
lab_extra_counts[3, 1] <- "Intermediate"
lab_extra_counts[4, 1] <- "Moderately High"
lab_extra_counts[5, 1] <- "High"
lab_extra_counts <- lab_extra_counts |> mutate(fake_scale = 2)
lab_extra_counts <- merge(lab_extra_counts, alpha_mapping, by = "n", all.x = TRUE)


### Plot 19: Lab engagement beyond assigned tasks ----

plot_lab_engagement_beyond <- ggplot(subset(lab_extra_counts, !is.na(Item)), aes(x = factor(Item, levels = c("Low", "Moderately Low", "Intermediate", "Moderately High", "High")), y = fake_scale, alpha = n)) +
  geom_bar(fill = bri_color["main"], alpha = subset(lab_extra_counts, !is.na(Item))$alpha, stat = "identity") +
  geom_text(aes(x = factor(Item, levels = c("Low", "Moderately Low", "Intermediate", "Moderately High", "High")), y = fake_scale / fake_scale, label = n, fontface = "bold"), color = bri_color["dark"], alpha = 1, size = 5) +
  bri_theme +
  labs(
    title = "How do you assess your lab's engagement with the project, beyond the assigned tasks?",
    subtitle = "Categories were presented as 1 = Low, 5 = High"
  ) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.title.x = element_blank())

# Save Plot 19
bri_ggsave(paste0(output_path, "/self-assessment/", "lab-engagement-beyond.png"), plot = plot_lab_engagement_beyond, width = 20, height = 5, units = "cm")

understanding_counts <- form2 |>
  count(`Como você avalia seu entendimento sobre os objetivos e metodologia do projeto como um todo?`) |>
  rename("Item" = "Como você avalia seu entendimento sobre os objetivos e metodologia do projeto como um todo?")
understanding_counts$Item <- as.character(understanding_counts$Item)

understanding_counts[1, 1] <- "Moderately Low"
understanding_counts[2, 1] <- "Intermediate"
understanding_counts[3, 1] <- "Moderately High"
understanding_counts[4, 1] <- "High"
understanding_counts <- rbind(understanding_counts, c("Low", 0))
understanding_counts$n <- as.numeric(understanding_counts$n)
understanding_counts <- understanding_counts |> mutate(fake_scale = 2)
understanding_counts <- merge(understanding_counts, alpha_mapping, by = "n", all.x = TRUE)

### Plot 20: Understanding of the project ----

plot_understanding <- ggplot(subset(understanding_counts, !is.na(Item)), aes(x = factor(Item, levels = c("Low", "Moderately Low", "Intermediate", "Moderately High", "High")), y = fake_scale, alpha = n)) +
  geom_bar(fill = bri_color["main"], alpha = subset(understanding_counts, !is.na(Item))$alpha, stat = "identity") +
  geom_text(aes(x = factor(Item, levels = c("Low", "Moderately Low", "Intermediate", "Moderately High", "High")), y = fake_scale / fake_scale, label = n, fontface = "bold"), color = bri_color["dark"], alpha = 1, size = 5) +
  bri_theme +
  labs(
    title = "How do you assess your understanding of the project?",
    subtitle = "Categories were presented as 1 = Low, 5 = High"
  ) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.title.x = element_blank())

# Save Plot 20
bri_ggsave(paste0(output_path, "/self-assessment/", "understanding-assessment.png"), plot = plot_understanding, width = 20, height = 5, units = "cm")

## Assessments of the project ----

overall_counts <- form2 |>
  count(`Qual a sua avaliação sobre a experiência da Iniciativa Brasileira de Reprodutibilidade?`) |>
  rename("Item" = "Qual a sua avaliação sobre a experiência da Iniciativa Brasileira de Reprodutibilidade?")
overall_counts$Item <- as.character(overall_counts$Item)

overall_counts[1, 1] <- "Very negative"
overall_counts[2, 1] <- "Negative"
overall_counts[3, 1] <- "Intermediate"
overall_counts[4, 1] <- "Positive"
overall_counts[5, 1] <- "Very positive"
overall_counts <- overall_counts |> mutate(fake_scale = 2)
overall_counts <- merge(overall_counts, alpha_mapping, by = "n", all.x = TRUE)

### Plot 21: Overall project experience assessment ----

plot_overall_experience <- ggplot(subset(overall_counts, !is.na(Item)), aes(x = factor(Item, levels = c("Very negative", "Negative", "Intermediate", "Positive", "Very positive")), y = fake_scale, alpha = n)) +
  geom_bar(fill = bri_color["main"], alpha = subset(overall_counts, !is.na(Item))$alpha, stat = "identity") +
  geom_text(aes(x = factor(Item, levels = c("Very negative", "Negative", "Intermediate", "Positive", "Very positive")), y = fake_scale / fake_scale, label = n, fontface = "bold"), color = bri_color["dark"], alpha = 1, size = 5) +
  bri_theme +
  labs(
    title = "What is your assessment of the overall experience of the project?",
    subtitle = "Categories were presented as 1 = Very negative, 5 = Very positive"
  ) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.title.x = element_blank())

# Save Plot 21
bri_ggsave(paste0(output_path, "/self-assessment/", "overall-project-assessment.png"), plot = plot_overall_experience, width = 20, height = 5, units = "cm")

workload_counts <- form2 |>
  count(`Como você avalia a carga de trabalho do projeto em relação à sua expectativa inicial?`) |>
  rename("Item" = "Como você avalia a carga de trabalho do projeto em relação à sua expectativa inicial?")
workload_counts$Item <- as.character(workload_counts$Item)

workload_counts[1, 1] <- "Below expectations"
workload_counts[2, 1] <- "Expectations met"
workload_counts[3, 1] <- "Above expectations"
workload_counts[4, 1] <- "Far above expectations"
workload_counts <- rbind(workload_counts, c("Far below expectations", 0))
workload_counts$n <- as.numeric(workload_counts$n)
workload_counts <- workload_counts |> mutate(fake_scale = 2)
workload_counts <- merge(workload_counts, alpha_mapping, by = "n", all.x = TRUE)

### Plot 22: Workload assessment compared to expectations ----

plot_workload_comparison <- ggplot(subset(workload_counts, !is.na(Item)), aes(x = factor(Item, levels = c("Far below expectations", "Below expectations", "Expectations met", "Above expectations", "Far above expectations")), y = fake_scale, alpha = n)) +
  geom_bar(fill = bri_color["main"], alpha = subset(workload_counts, !is.na(Item))$alpha, stat = "identity") +
  geom_text(aes(x = factor(Item, levels = c("Far below expectations", "Below expectations", "Expectations met", "Above expectations", "Far above expectations")), y = fake_scale / fake_scale, label = n, fontface = "bold"), color = bri_color["dark"], alpha = 1, size = 5) +
  bri_theme +
  labs(
    title = "What is your assessment of the workload of the project compared to your initial expectations?",
    subtitle = "Categories were presented as 1 = Far below expectations, 5 = Far above expectations"
  ) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.title.x = element_blank())

# Save Plot 22
bri_ggsave(paste0(output_path, "/self-assessment/", "workload-assessment.png"), plot = plot_workload_comparison, width = 20, height = 5, units = "cm")

## Assessments of the coordinating team ----
coord_presentation_counts <- form2 |>
  count(`Qual a sua avaliação sobre o trabalho da equipe coordenadora da Iniciativa em relação a: [Apresentação do projeto]`) |>
  rename("Item" = "Qual a sua avaliação sobre o trabalho da equipe coordenadora da Iniciativa em relação a: [Apresentação do projeto]")
coord_presentation_counts$Item <- as.character(coord_presentation_counts$Item)

coord_presentation_counts[1, 1] <- "Negative"
coord_presentation_counts[2, 1] <- "Intermediate"
coord_presentation_counts[3, 1] <- "Positive"
coord_presentation_counts[4, 1] <- "Very positive"
coord_presentation_counts <- rbind(coord_presentation_counts, c("Very negative", 0))
coord_presentation_counts$n <- as.numeric(coord_presentation_counts$n)
coord_presentation_counts <- coord_presentation_counts |> mutate(fake_scale = 2)
coord_presentation_counts <- merge(coord_presentation_counts, alpha_mapping, by = "n", all.x = TRUE)

### Plot 23: Assessment of the coordinating committee regarding project presentation ----

plot_coord_presentation <- ggplot(subset(coord_presentation_counts, !is.na(Item)), aes(x = factor(Item, levels = c("Very negative", "Negative", "Intermediate", "Positive", "Very positive")), y = fake_scale, alpha = n)) +
  geom_bar(fill = bri_color["main"], alpha = subset(coord_presentation_counts, !is.na(Item))$alpha, stat = "identity") +
  geom_text(aes(x = factor(Item, levels = c("Very negative", "Negative", "Intermediate", "Positive", "Very positive")), y = fake_scale / fake_scale, label = n, fontface = "bold"), color = bri_color["dark"], alpha = 1, size = 5) +
  bri_theme +
  labs(
    title = "What is your assessment of the coordinating committee \nregarding presentation of the project?",
    subtitle = "Categories were presented as 1 = Very negative, 5 = Very positive"
  ) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.title.x = element_blank()) +
  scale_alpha_continuous(breaks = alpha_mapping$variable, labels = alpha_mapping$alpha)

# Save Plot 23
bri_ggsave(paste0(output_path, "/self-assessment/", "lik-coord-presentation.png"), plot = plot_coord_presentation, width = 20, height = 5, units = "cm")

coord_commspeed_counts <- form2 |>
  count(`Qual a sua avaliação sobre o trabalho da equipe coordenadora da Iniciativa em relação a: [Agilidade na comunicação]`) |>
  rename("Item" = "Qual a sua avaliação sobre o trabalho da equipe coordenadora da Iniciativa em relação a: [Agilidade na comunicação]")
coord_commspeed_counts$Item <- as.character(coord_commspeed_counts$Item)

coord_commspeed_counts[1, 1] <- "Negative"
coord_commspeed_counts[2, 1] <- "Intermediate"
coord_commspeed_counts[3, 1] <- "Positive"
coord_commspeed_counts[4, 1] <- "Very positive"
coord_commspeed_counts <- rbind(coord_commspeed_counts, c("Very negative", 0))
coord_commspeed_counts$n <- as.numeric(coord_commspeed_counts$n)
coord_commspeed_counts <- coord_commspeed_counts |> mutate(fake_scale = 2)
coord_commspeed_counts <- merge(coord_commspeed_counts, alpha_mapping, by = "n", all.x = TRUE)

### Plot 24: Assessment of the coordinating committee regarding speed of communication ----

plot_coord_commspeed <- ggplot(subset(coord_commspeed_counts, !is.na(Item)), aes(x = factor(Item, levels = c("Very negative", "Negative", "Intermediate", "Positive", "Very positive")), y = fake_scale, alpha = n)) +
  geom_bar(fill = bri_color["main"], alpha = subset(coord_commspeed_counts, !is.na(Item))$alpha, stat = "identity") +
  geom_text(aes(x = factor(Item, levels = c("Very negative", "Negative", "Intermediate", "Positive", "Very positive")), y = fake_scale / fake_scale, label = n, fontface = "bold"), color = bri_color["dark"], alpha = 1, size = 5) +
  bri_theme +
  labs(
    title = "What is your assessment of the coordinating committee \nregarding speed of communication?",
    subtitle = "Categories were presented as 1 = Very negative, 5 = Very positive"
  ) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.title.x = element_blank()) +
  scale_alpha_continuous(breaks = alpha_mapping$variable, labels = alpha_mapping$alpha)

# Save Plot 24
bri_ggsave(paste0(output_path, "/self-assessment/", "lik-coord-speed.png"), plot = plot_coord_commspeed, width = 20, height = 5, units = "cm")

coord_clarity_counts <- form2 |>
  count(`Qual a sua avaliação sobre o trabalho da equipe coordenadora da Iniciativa em relação a: [Clareza nas orientações]`) |>
  rename("Item" = "Qual a sua avaliação sobre o trabalho da equipe coordenadora da Iniciativa em relação a: [Clareza nas orientações]")
coord_clarity_counts$Item <- as.character(coord_clarity_counts$Item)

coord_clarity_counts[1, 1] <- "Very negative"
coord_clarity_counts[2, 1] <- "Negative"
coord_clarity_counts[3, 1] <- "Intermediate"
coord_clarity_counts[4, 1] <- "Positive"
coord_clarity_counts[5, 1] <- "Very positive"
coord_clarity_counts <- coord_clarity_counts |> mutate(fake_scale = 2)
coord_clarity_counts <- merge(coord_clarity_counts, alpha_mapping, by = "n", all.x = TRUE)

### Plot 25: Assessment of the coordinating committee regarding clarity in instructions ----

plot_coord_clarity <- ggplot(subset(coord_clarity_counts, !is.na(Item)), aes(x = factor(Item, levels = c("Very negative", "Negative", "Intermediate", "Positive", "Very positive")), y = fake_scale, alpha = n)) +
  geom_bar(fill = bri_color["main"], alpha = subset(coord_clarity_counts, !is.na(Item))$alpha, stat = "identity") +
  geom_text(aes(x = factor(Item, levels = c("Very negative", "Negative", "Intermediate", "Positive", "Very positive")), y = fake_scale / fake_scale, label = n, fontface = "bold"), color = bri_color["dark"], alpha = 1, size = 5) +
  bri_theme +
  labs(
    title = "What is your assessment of the coordinating committee \nregarding clarity in instructions?",
    subtitle = "Categories were presented as 1 = Very negative, 5 = Very positive"
  ) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.title.x = element_blank())

# Save Plot 25
bri_ggsave(paste0(output_path, "/self-assessment/", "lik-coord-instructions.png"), plot = plot_coord_clarity, width = 20, height = 5, units = "cm")

coord_support_counts <- form2 |>
  count(`Qual a sua avaliação sobre o trabalho da equipe coordenadora da Iniciativa em relação a: [Material de apoio disponibilizado]`) |>
  rename("Item" = "Qual a sua avaliação sobre o trabalho da equipe coordenadora da Iniciativa em relação a: [Material de apoio disponibilizado]")
coord_support_counts$Item <- as.character(coord_support_counts$Item)

coord_support_counts[1, 1] <- "Very negative"
coord_support_counts[2, 1] <- "Negative"
coord_support_counts[3, 1] <- "Intermediate"
coord_support_counts[4, 1] <- "Positive"
coord_support_counts[5, 1] <- "Very positive"
coord_support_counts <- coord_support_counts |> mutate(fake_scale = 2)
coord_support_counts <- merge(coord_support_counts, alpha_mapping, by = "n", all.x = TRUE)

### Plot 26: Assessment of the coordinating committee regarding support material provided ----

plot_coord_support <- ggplot(subset(coord_support_counts, !is.na(Item)), aes(x = factor(Item, levels = c("Very negative", "Negative", "Intermediate", "Positive", "Very positive")), y = fake_scale, alpha = n)) +
  geom_bar(fill = bri_color["main"], alpha = subset(coord_support_counts, !is.na(Item))$alpha, stat = "identity") +
  geom_text(aes(x = factor(Item, levels = c("Very negative", "Negative", "Intermediate", "Positive", "Very positive")), y = fake_scale / fake_scale, label = n, fontface = "bold"), color = bri_color["dark"], alpha = 1, size = 5) +
  bri_theme +
  labs(
    title = "What is your assessment of the coordinating committee \nregarding the support material provided?",
    subtitle = "Categories were presented as 1 = Very negative, 5 = Very positive"
  ) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.title.x = element_blank())

# Save Plot 26
bri_ggsave(paste0(output_path, "/self-assessment/", "lik-coord-support.png"), plot = plot_coord_support, width = 20, height = 5, units = "cm")

coord_availab_counts <- form2 |>
  count(`Qual a sua avaliação sobre o trabalho da equipe coordenadora da Iniciativa em relação a: [Disponibilidade para auxílio aos laboratórios]`) |>
  rename("Item" = "Qual a sua avaliação sobre o trabalho da equipe coordenadora da Iniciativa em relação a: [Disponibilidade para auxílio aos laboratórios]")
coord_availab_counts$Item <- as.character(coord_availab_counts$Item)

coord_availab_counts[1, 1] <- "Very negative"
coord_availab_counts[2, 1] <- "Negative"
coord_availab_counts[3, 1] <- "Intermediate"
coord_availab_counts[4, 1] <- "Positive"
coord_availab_counts[5, 1] <- "Very positive"
coord_availab_counts <- coord_availab_counts |> mutate(fake_scale = 2)
coord_availab_counts <- merge(coord_availab_counts, alpha_mapping, by = "n", all.x = TRUE)

### Plot 27: Assessment of the coordinating committee regarding availability to support labs ----

plot_coord_availability <- ggplot(subset(coord_availab_counts, !is.na(Item)), aes(x = factor(Item, levels = c("Very negative", "Negative", "Intermediate", "Positive", "Very positive")), y = fake_scale, alpha = n)) +
  geom_bar(fill = bri_color["main"], alpha = subset(coord_availab_counts, !is.na(Item))$alpha, stat = "identity") +
  geom_text(aes(x = factor(Item, levels = c("Very negative", "Negative", "Intermediate", "Positive", "Very positive")), y = fake_scale / fake_scale, label = n, fontface = "bold"), color = bri_color["dark"], alpha = 1, size = 5) +
  bri_theme +
  labs(
    title = "What is your assessment of the coordinating committee \nregarding the availability to support the labs?",
    subtitle = "Categories were presented as 1 = Very negative, 5 = Very positive"
  ) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.title.x = element_blank())

# Save Plot 27
bri_ggsave(paste0(output_path, "/self-assessment/", "lik-coord-availability.png"), plot = plot_coord_availability, width = 20, height = 5, units = "cm")

coord_transp_counts <- form2 |>
  count(`Qual a sua avaliação sobre o trabalho da equipe coordenadora da Iniciativa em relação a: [Transparência nas decisões tomadas]`) |>
  rename("Item" = "Qual a sua avaliação sobre o trabalho da equipe coordenadora da Iniciativa em relação a: [Transparência nas decisões tomadas]")
coord_transp_counts$Item <- as.character(coord_transp_counts$Item)

coord_transp_counts[1, 1] <- "Very negative"
coord_transp_counts[2, 1] <- "Negative"
coord_transp_counts[3, 1] <- "Intermediate"
coord_transp_counts[4, 1] <- "Positive"
coord_transp_counts[5, 1] <- "Very positive"
coord_transp_counts <- coord_transp_counts |> mutate(fake_scale = 2)
coord_transp_counts <- merge(coord_transp_counts, alpha_mapping, by = "n", all.x = TRUE)


### Plot 28: Assessment of the coordinating committee regarding transparency in decision-making ----

plot_coord_transparency <- ggplot(subset(coord_transp_counts, !is.na(Item)), aes(x = factor(Item, levels = c("Very negative", "Negative", "Intermediate", "Positive", "Very positive")), y = fake_scale, alpha = n)) +
  geom_bar(fill = bri_color["main"], alpha = subset(coord_transp_counts, !is.na(Item))$alpha, stat = "identity") +
  geom_text(aes(x = factor(Item, levels = c("Very negative", "Negative", "Intermediate", "Positive", "Very positive")), y = fake_scale / fake_scale, label = n, fontface = "bold"), color = bri_color["dark"], alpha = 1, size = 5) +
  bri_theme +
  labs(
    title = "What is your assessment of the coordinating committee \nregarding the transparency in decision-making?",
    subtitle = "Categories were presented as 1 = Very negative, 5 = Very positive"
  ) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.title.x = element_blank())

# Save Plot 28
bri_ggsave(paste0(output_path, "/self-assessment/", "lik-coord-transparency.png"), plot = plot_coord_transparency, width = 20, height = 5, units = "cm")


## Subjective replication categories -----

Avaliacao_subjetiva_final <- read_excel("./other-data/self-assessment/Avaliacao_subjetiva_final.xlsx", sheet = "resolvidas")

contagens.subj.Y <- full_join(Avaliacao_subjetiva_final |> filter(`replicou` == "Sim") |> count(`Categoria1`), Avaliacao_subjetiva_final |> filter(`replicou` == "Sim") |> count(`Categoria2`), join_by("Categoria1" == "Categoria2"))
contagens.subj.Y

contagens.subj.N <- full_join(Avaliacao_subjetiva_final |> filter(`replicou` == "Nao") |> count(`Categoria1`), Avaliacao_subjetiva_final |> filter(`replicou` == "Nao") |> count(`Categoria2`), join_by("Categoria1" == "Categoria2"))
contagens.subj.N <- full_join(contagens.subj.N, Avaliacao_subjetiva_final |> filter(`replicou` == "Nao") |> count(`Categoria3`), join_by("Categoria1" == "Categoria3"))
contagens.subj.N
