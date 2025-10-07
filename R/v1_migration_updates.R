#' @keywords internal
immunarch_v1_update_oct_2025 <- function() {
  cli::cli_h1("{cli::col_green('immunarch')} {cli::col_yellow('0.10.0')} -- Critical Pre-release Notice")

  cli::cli_alert_warning("Update #1 [Oct 2025] -- Major changes are coming in {cli::col_green('immunarch')} {cli::col_yellow('1.0.0')}!")
  cli::cli_text(cli::col_yellow(cli::spark_line(runif(110, 0, 1))))

  cli::cli_par()
  cli::cli_text()
  cli::cli_text(
    "Hi, this is Vadim Nazarov speaking - author of {cli::col_green('immunarch')}. ",
    "{cli::col_green('immunarch')} is finally graduating out of the {cli::col_yellow('0.x.y')} development cycle. ",
    "I'm preparing our {cli::col_yellow('1.0.0')} release, which will remain stable and free of sudden changes until we approach {cli::col_yellow('2.0.0')}. ",
    "A scientific publication will accompany it for proper citation. ",
    "Significant changes are coming, and I want to ensure you have everything you need to migrate to the new version."
  )

  cli::cli_par()
  cli::cli_text()
  cli::cli_text("Here's a preview of what's coming in {cli::col_green('immunarch')} {cli::col_yellow('1.0.0')}:")
  cli::cli_bullets(c(
    "i" = "Some computationally intensive or advanced features (e.g., distance computations, graph-based analyses, dimensionality reduction techniques) will move to separate packages, making {cli::col_green('immunarch')} much lighter to install and manage.",
    "i" = "New functions will replace older ones to make code more readable and maintainable. Legacy functions will remain temporarily, but they won't be updated and will be removed {cli::col_yellow('around 2027')}.",
    "i" = "We will discontinue support for most custom file formats because the AIRR ecosystem is now mature enough; most tools adhere to the AIRR standard.",
    "i" = "The package will transition from data frames to the new {cli::col_blue('ImmunData')} structure, better suited for handling larger, more complex, and multimodal datasets (e.g., single-cell, spatial).",
    "i" = "{cli::col_blue('ImmunData')} is available in the separate {cli::col_blue('immundata')} package, which you can already install via {cli::col_cyan('pak::pkg_install(\"immundata\")')}.",
    "i" = "The {cli::col_blue('ImmunData')}-based computations will be significantly faster, will support datasets larger than RAM, and will fully adhere to AIRR Community standards.",
    "i" = "There are currently only a handful of functions that implement {cli::col_blue('ImmunData')}-based computations. However, if you want to start learning it, or you have large-scale data, now is the best time: tutorials are available at {.url https://github.com/immunomind/immundata} and {.url https://immunomind.github.io/docs/}."
  ))

  cli::cli_par()
  cli::cli_text()
  cli::cli_text(
    "See the dedicated migration guide for what you can do now and how to prepare for the future:"
  )
  cli::cli_text(">> visit {.url https://immunomind.github.io/docs/tutorials/migration}")

  cli::cli_par()
  cli::cli_text()
  cli::cli_text(
    "See the comprehensive tutorial on how to analyse single-cell AIRR data:"
  )
  cli::cli_text(">> visit {.url https://immunomind.github.io/docs/tutorials/single_cell}")

  cli::cli_par()
  cli::cli_text()
  cli::cli_alert_success("Thank you for supporting {cli::col_green('immunarch')} from its early days. Your feedback, contributions, and trust have driven its evolution, and I deeply appreciate it.")

  cli::cli_par()
  cli::cli_text()
  cli::cli_alert_info("Questions, comments, ideas? I'm available via:")
  cli::cli_text(">> Support email: {.url mailto:support@immunomind.com}")
  cli::cli_text(">> GitHub tickets: {.url https://github.com/immunomind/immunarch}")
  cli::cli_text(">> LinkedIn: {.url https://www.linkedin.com/in/vdnaz/}")

  cli::cli_par()
  cli::cli_text()
  cli::cli_text("--")
  cli::cli_text("Vadim I. Nazarov")
}
