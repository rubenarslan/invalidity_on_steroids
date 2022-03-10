render_job <- function (input, params = list(), output_file) 
{
  tmp <- tempfile()
  .rcamisc_global_input_file_for_rstudio_job <<- input
  .rcamisc_global_params_for_rstudio_job <<- params
  .rcamisc_global_output_file_for_rstudio_job <<- output_file
  rmarkdown_render <- rmarkdown::render
  cat("rmarkdown::render_site(input = .rcamisc_global_input_file_for_rstudio_job)", 
      file = tmp, sep = "\n")
  rstudioapi::jobRunScript(tmp, name = input, workingDir = getwd(), output_file, 
                           importEnv = TRUE)
}


render_job("1_biocycle_import.Rmd", output_file = NULL)
render_job("2_impute_by_day.Rmd", output_file = NULL)

render_job("3_biocycle_summarize.Rmd", output_file = NULL)
render_job("1_roney_import.Rmd", output_file = NULL)
render_job("1_gocd2_import.Rmd", output_file = NULL)
render_job("1_goettingen_lab2_import.Rmd", output_file = NULL)
render_job("1_ocmate_import.Rmd", output_file = NULL)
render_job("1_grebe_import.Rmd", output_file = NULL)
render_job("1_blake_import.Rmd", output_file = NULL)
render_job("1_marcinkowska_import.Rmd", output_file = NULL)

render_job("index.Rmd", output_file = NULL)
render_job("4_plots_for_paper.Rmd", output_file = NULL)
