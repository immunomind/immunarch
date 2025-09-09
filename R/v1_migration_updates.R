#' @keywords internal
immunarch_v1_update_sep_2025 <- function() {
  cli::cli_h1("{cli::col_green('immunarch')} {cli::col_yellow('0.9.x')} -- Critical Pre-release Notice")

  cli::cli_alert_warning("Update #1 [Sep 2025] -- Major changes are coming in {cli::col_green('immunarch')} {cli::col_yellow('1.0.0')}!")
  cli::cli_text(cli::col_yellow(cli::spark_line(runif(110, 0, 1))))

  cli::cli_par()
  cli::cli_text()
  cli::cli_text(
    "Hi, this is Vadim Nazarov speaking -- author of {cli::col_green('immunarch')}. ",
    "{cli::col_green('immunarch')} is finally graduating from out of the {cli::col_yellow('0.x.y')} development cycle. ",
    "I'm preparing our {cli::col_yellow('1.0.0')} release, which will remain stable and free of sudden changes until we approach {cli::col_yellow('2.0.0')}, along with ",
    "a scientific publication for proper citations. ",
    "Significant changes are coming, and I want to ensure you have everything you need to migrate to the new version."
  )

  cli::cli_par()
  cli::cli_text()
  cli::cli_text("Here's a preview of what's coming in {cli::col_green('immunarch')} {cli::col_yellow('1.0.0')}:")
  cli::cli_bullets(c(
    "i" = "Some computationally intensive or advanced features (e.g., distance computations, graph-based analyses, dimensionality reduction techniques) will move to separate packages, making {cli::col_green('immunarch')} much more lightweight to install and manage;",
    "i" = "New functions will be introduced instead of the left old ones to make code more readable and maintainable. Legacy functions will remain temporarily, but they won't be updated and will be removed by {cli::col_yellow('~2027')};",
    "i" = "We will discontinue support for most custom file formats because the AIRR ecosystem is now mature enough that the majority of tools adhere to the AIRR standard;",
    "i" = "The package will transition from data frames to the new {cli::col_blue('ImmunData')} structure -- better suited for handling modern larger, more complex, and multi-modal datasets (e.g., single-cell, spatial);",
    "i" = "{cli::col_blue('ImmunData')} is available in the separate {cli::col_blue('immundata')} package, which you can already install via {cli::col_cyan('pak::pkg_install(\"immundata\")')};",
    "i" = "The {cli::col_blue('ImmunData')}-based computations will be significantly faster, support datasets larger than RAM, and fully adhere to AIRR Community standards.",
    "i" = "There currently only a handful functions which implement {cli::col_blue('ImmunData')}-based computations. However, if you want to start learning it, or you have a large-scale data, now it is the best time: the tutorials are available on {cli::col_cyan('https://github.com/immunomind/immundata') and {cli::col_cyan('https://immunomind.github.io/docs/')}"
  ))

  cli::cli_par()
  cli::cli_text()
  cli::cli_text(
    "See the dedicated migration guide for migration on what you can do now and how to prepare for the future:"
  )
  cli::cli_text(">> visit {cli::col_cyan('https://immunomind.github.io/docs/tutorials/migration')}")

  cli::cli_par()
  cli::cli_text()
  cli::cli_alert_success("Thank you for supporting {cli::col_green('immunarch')} from its early days. Your feedback, contributions, and trust have driven its evolution, and I deeply appreciate it.")

  cli::cli_par()
  cli::cli_text()
  cli::cli_alert_info("Questions, comments, ideas? I'm available via:")
  cli::cli_text(">> Support email: {cli::col_cyan('support@immunomind.com')}")
  cli::cli_text(">> GitHub tickets: {cli::col_cyan('https://github.com/immunomind/immunarch')}")
  cli::cli_text(">> LinkedIn: {cli::col_cyan('https://www.linkedin.com/in/vdnaz/')}")

  cli::cli_par()
  cli::cli_text()
  cli::cli_text("--")
  cli::cli_text("Vadim I. Nazarov")
}
