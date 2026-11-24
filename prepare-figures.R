# This script just reorganizes all figures for presentation in one folder, with meaningful names for the manuscript

output_path = paste0("./output/", results_path)

# Restore tables folder if it was previously renamed to "other tables"
tables_path <- paste0(output_path, "/_manuscript figures and tables/tables")
other_tables_path <- paste0(output_path, "/_manuscript figures and tables/other tables")

if (!dir.exists(tables_path) && dir.exists(other_tables_path)) {
  file.rename(other_tables_path, tables_path)
  message("Restored 'other tables' back to 'tables'")
}

# Preprint ----------------------------------------------------------------

figure_correspondence_table = read_excel("./other-data/Manuscript Figure Correspondence.xlsx", sheet = 1)

paper_fig_path = paste0("./output/", results_path, "/_manuscript figures and tables/manuscript/")
suppressWarnings({dir.create(paper_fig_path)})
suppressWarnings({dir.create(paste0(paper_fig_path, "/tables"))})

# Figures not cited in the manuscript go unnumbered into "other figures"
other_fig_path = paste0("./output/", results_path, "/_manuscript figures and tables/other figures/")
suppressWarnings({dir.create(other_fig_path)})

# Remove previously exported figures so renamed or renumbered figures do not leave stale files
file.remove(list.files(path = paper_fig_path, pattern = "^Figure_", full.names = TRUE))


message(paste0("Saving figures and tables to ", results_path))

pwalk(figure_correspondence_table, function (code_generated_filename, manuscript_name, type, ...) {
  fnf = paste0(output_path, code_generated_filename)
  fnt = paste0(paper_fig_path, str_to_title(type), "_", manuscript_name, str_extract(basename(fnf), "\\..+?$"))
  
  if(type == "table"){
    fnt = paste0(paper_fig_path, "/tables/", str_to_title(type), " ", manuscript_name, str_extract(basename(fnf), "\\..+?$"))
  } 

  if(type == "other figure"){
    fnt = paste0(other_fig_path, manuscript_name, str_extract(basename(fnf), "\\..+?$"))
  }

  if (!is.na(code_generated_filename)) {
    if (file.exists(fnf)) {
    file.copy(from = fnf, to = fnt, overwrite = T)
    message(paste0(basename(fnt), " <-- ", fnf))
    } else {
      message(basename(fnt), paste0(" <-- ERROR - file not found: ", fnf))
    }
  } else {
    message(paste0(basename(fnt), " <-- SKIP"))
  }
})

# Deleting tables folder
file.remove(list.files(path = paste0(output_path, "/_manuscript figures and tables/tables"),
                       pattern = "Table[ _](S\\d+|[123])",
                       full.names = TRUE,
                       ignore.case = TRUE))

# Remove destination folder if it exists before renaming
other_tables_path <- paste0(output_path, "/_manuscript figures and tables/other tables")
if (dir.exists(other_tables_path)) {
  unlink(other_tables_path, recursive = TRUE)
}

file.rename(paste0(output_path, "/_manuscript figures and tables/tables"), other_tables_path)
