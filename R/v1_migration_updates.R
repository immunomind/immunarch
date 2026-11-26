# We pregenerated and hardcoded these values to preserve the user's random seed.
IMMUNARCH_NEWS_SPARKLINE <- c(
  34, 92, 74, 20, 23, 74, 53, 82, 76, 3,
  68, 88, 96, 4, 32, 11, 99, 19, 68, 13,
  35, 49, 47, 62, 41, 40, 48, 41, 94, 69,
  53, 9, 24, 28, 43, 51, 44, 12, 46, 53,
  73, 75, 95, 87, 31, 68, 96, 25, 53, 25,
  67, 73, 86, 16, 58, 100, 20, 9, 7, 16,
  91, 96, 11, 33, 68, 58, 40, 58, 46, 41,
  86, 51, 9, 82, 19, 25, 23, 14, 21, 72,
  25, 88, 51, 15, 59, 35, 91, 30, 49, 29,
  20, 16, 93, 92, 46, 44, 21, 46, 29, 82,
  8, 54, 3, 58, 41, 52, 39, 52, 32, 55
) / 100

#' @keywords internal
immunarch_v1_update_oct_2025 <- function() {
  cli::cli_h1("{cli::col_green('immunarch')} {cli::col_yellow('1.0')} -- Critical Pre-release Notice")

  cli::cli_alert_warning("Update #1 [Oct 2025] -- Major changes are coming in {cli::col_green('immunarch')} {cli::col_yellow('1.0')}!")
  cli::cli_text(cli::col_yellow(cli::spark_line(IMMUNARCH_NEWS_SPARKLINE)))

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

immunarch_v1_update_nov_2026 <- function() {
  cli::cli_h1("{cli::col_green('immunarch')} {cli::col_yellow('1.0')} -- Critical Pre-release Notice")

  cli::cli_alert_warning(
    "Update #2 [Nov 2026] -- Migrate your scripts before the next release in January"
  )
  cli::cli_text(cli::col_yellow(cli::spark_line(IMMUNARCH_NEWS_SPARKLINE)))

  cli::cli_par()
  cli::cli_text()
  cli::cli_text("Hi, this is Vadim Nazarov, the author of {cli::col_green('immunarch')}.")

  cli::cli_par()
  cli::cli_text()
  cli::cli_text(
    "The next {cli::col_green('immunarch')} version {cli::col_yellow('0.12.0')}, in January 2027, will remove legacy functions that currently raise deprecation warnings, ",
    "such as {.fn repLoad}, {.fn repDiversity}, {.fn repClonality}, and {.fn repOverlap}. ",
    "Code that calls these functions {cli::col_red('will stop working')} after you upgrade."
  )

  cli::cli_par()
  cli::cli_text()
  cli::cli_text(
    "Most migrations should be straightforward. ",
    "Start by asking a coding agent to update your code using the migration guide (link below). ",
    "Run the updated code and check that it gives the results you expect."
  )

  cli::cli_par()
  cli::cli_text()
  cli::cli_text("If you still need help, open a {.href [GitHub issue](https://github.com/immunomind/immunarch/issues)}.")
  cli::cli_text("Describe the problem and, ideally, a small example that reproduces it.")

  cli::cli_par()
  cli::cli_text()
  cli::cli_bullets(c(
    "*" = "Comprehensive tutorial on how to work with the new {cli::col_green('immunarch')} {.url https://immunomind.github.io/docs/tutorials/single_cell/}",
    "*" = "Documentation for all ecosystem packages {.url https://immunomind.github.io/docs/}",
    "*" = "Migration guide {.url https://immunomind.github.io/docs/tutorials/migration/}",
    "*" = "New package to streamline data loading {cli::col_yellow('immunread')} {.url https://immunomind.github.io/immunread/}",
    "*" = "New package to call external immune repertoire tools and access databases from R {cli::make_ansi_style('orange')('immundapter')} {.url https://immunomind.github.io/immundapter/}",
    "*" = "Follow updates for all packages in the ecosystem {.url https://github.com/immunomind/immunarch/issues/432}"
  ))

  cli::cli_par()
  cli::cli_text()
  cli::cli_text("--")
  cli::cli_text("Vadim I. Nazarov")
}
