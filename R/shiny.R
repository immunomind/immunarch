if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("data", "carat", "price", "color", "clarity", "diamonds"))
}


#' Manipulate ggplot plots and create publication-ready plots
#'
#' @concept fixvis
#'
#' @aliases fixVis
#'
#' @importFrom graphics plot
#' @importFrom shiny fluidPage titlePanel sidebarLayout sidebarPanel downloadButton actionButton br
#' @importFrom shiny tabsetPanel tabPanel textOutput sliderInput checkboxInput selectInput textInput textAreaInput
#' @importFrom shiny mainPanel uiOutput renderText renderUI renderPlot plotOutput observe updateSliderInput
#' @importFrom shiny observeEvent downloadHandler splitLayout
#' @importFrom shinythemes shinytheme
#'
#' @description The \code{fixVis} is a built-in software tool for the manipulation
#' of plots, such as adjusting title text font and size, axes, and more. It is a powerful
#' tool designed to produce publication-ready plots with minimal amount of coding.
#'
#' @param .plot A ggplot2 plot.
#'
#' @return No return value because it is an application.
#'
#' @examples
#' if (interactive()) {
#'   # Compute gene usage, visualise it and tweak via fixVis
#'   data(immdata) # load test data
#'   gu <- geneUsage(immdata$data)
#'   p <- vis(gu)
#'   fixVis(p)
#' }
#' @export fixVis
fixVis <- function(.plot = NA) {
  make_legend_tab <- function(.label, .name, .is.title) {
    .full <- function(.l) {
      stringr::str_c("legend", .label, .l, sep = "_")
    }

    objs <- list(title = .name)

    objs <- c(objs, list(br()))

    if (.is.title) {
      objs <- c(objs, list(
        splitLayout(
          cellWidths = c("60%", "40%"),
          checkboxInput(.full("remove"), "Remove the legend"),
          checkboxInput(.full("contin"), "Continuous?")
        ),
        sliderInput(.full("ncol"), "Number of columns:",
          min = 1, max = 40, value = 1, step = 1
        ),
        br(),
        textInput(.full("text"), "Title text:", .name, placeholder = "Samples")
      ))
    }

    objs <- c(objs, list(
      sliderInput(.full("size"), "Text size:",
        min = 1, max = 40, value = ifelse(.is.title, 16, 11), step = .5
      ),
      sliderInput(.full("hjust"), "Text horizontal adjustment:",
        min = 0, max = 1, value = 0, step = .05
      ),
      sliderInput(.full("vjust"), "Text vertical adjustment:",
        min = -4, max = 4, value = .5, step = .25
      ),
      sliderInput(.full("angle"), "Text angle:",
        min = 0, max = 90, value = 0, step = 1
      ),
      selectInput(
        .full("face"), "Face:",
        list(Plain = "plain", Bold = "bold", Italic = "italic", "Bold Italic" = "bold.italic")
      )
    ))

    do.call(tabPanel, objs)
  }

  if (has_no_data(.plot)) {
    diamonds <- ggplot2::diamonds
    .plot <- qplot(x = carat, y = price, fill = cut, shape = cut, color = color, size = clarity, data = diamonds[sample.int(nrow(diamonds), 5000), ]) + theme_classic()
  }

  #
  #### UI ####
  #
  ui <- fluidPage(
    theme = shinytheme("cosmo"),
    titlePanel("FixVis: make your plots publication-ready already!"),
    sidebarLayout(
      sidebarPanel(
        downloadButton("save_plot", "Save"),
        actionButton("console_plot", "Plot to R console"),
        br(),
        br(),
        tabsetPanel(
          tabPanel(
            "General",
            br(),
            textOutput("save_text"),
            br(),
            # textOutput("save_text2"),
            # br(),
            sliderInput("plot_width", "Plot width (in):", min = 2, max = 24, value = 8),
            sliderInput("plot_height", "Plot height (in):", min = 2, max = 20, value = 5),
            checkboxInput("coord_flip", "Flip coordinates"),
            # checkboxInput("do_interactive", "Interactive plot"),
            selectInput("ggplot_theme", "Theme",
              selected = "Pubr",
              list(
                "Linedraw",
                "Black-white",
                "Grey / gray",
                "Light",
                "Dark",
                "Minimal",
                "Classic",
                "Pubr",
                "Pubr w/o lines"
              )
            )
          ),
          tabPanel(
            "Title & subtitle",
            br(),
            tabsetPanel(
              tabPanel(
                "Title",
                br(),
                textInput("title_text", "Title text:", "nice title text", placeholder = "Gene usage"),
                sliderInput("title_text_size", "Title text size:",
                  min = 1, max = 40, value = 25, step = .5
                ),
                sliderInput("title_text_hjust", "Title text horizontal adjustment:",
                  min = 0, max = 1, value = 0, step = .05
                ),
                sliderInput("title_text_vjust", "Title text vertical adjustment:",
                  min = -4, max = 4, value = .5, step = .25
                ),
                sliderInput("title_text_angle", "Title text angle:",
                  min = 0, max = 90, value = 0, step = 1
                ),
                selectInput(
                  "title_face", "Face:",
                  list(Plain = "plain", Bold = "bold", Italic = "italic", "Bold Italic" = "bold.italic")
                )
              ),
              tabPanel(
                "Subtitle",
                br(),
                textAreaInput("subtitle_text", "Subtitle text:", "nice subtitle text",
                  placeholder = "Frequency of Variable gene segments presented in the input samples"
                ),
                sliderInput("subtitle_text_size", "Subtitle text size:",
                  min = 1, max = 40, value = 16, step = .5
                ),
                sliderInput("subtitle_text_hjust", "Subtitle text horizontal adjustment:",
                  min = 0, max = 1, value = 0, step = .05
                ),
                sliderInput("subtitle_text_vjust", "Subtitle text vertical adjustment:",
                  min = -4, max = 4, value = .5, step = .25
                ),
                sliderInput("subtitle_text_angle", "Subtitle text angle:",
                  min = 0, max = 90, value = 0, step = 1
                ),
                selectInput(
                  "subtitle_face", "Face:",
                  list(Plain = "plain", Bold = "bold", Italic = "italic", "Bold Italic" = "bold.italic")
                )
              )
            )
          ),
          tabPanel(
            "Legends",
            br(),
            selectInput(
              "legend_position", "Legend position",
              list(
                "right",
                "top",
                "bottom",
                "left"
              )
            ),
            selectInput(
              "legend_box", "Legend arrangement",
              list(
                "vertical",
                "horizontal"
              )
            ),
            tabsetPanel(
              tabPanel(
                "Color",
                tabsetPanel(
                  make_legend_tab("col_title", "Title (color)", TRUE),
                  make_legend_tab("col_text", "Labels (color)", FALSE)
                )
              ),
              tabPanel(
                "Fill",
                tabsetPanel(
                  make_legend_tab("fill_title", "Title (fill)", TRUE),
                  make_legend_tab("fill_text", "Labels (fill)", FALSE)
                )
              ),
              tabPanel(
                "Size",
                tabsetPanel(
                  make_legend_tab("size_title", "Title (size)", TRUE),
                  make_legend_tab("size_text", "Labels (size)", FALSE)
                )
              ),
              tabPanel(
                "Shape",
                tabsetPanel(
                  make_legend_tab("shape_title", "Title (shape)", TRUE),
                  make_legend_tab("shape_text", "Labels (shape)", FALSE)
                )
              ),
              tabPanel(
                "Linetype",
                tabsetPanel(
                  make_legend_tab("linetype_title", "Title (linetype)", TRUE),
                  make_legend_tab("linetype_text", "Labels (linetype)", FALSE)
                )
              )
            )
          ),
          tabPanel(
            "X axis",
            br(),
            tabsetPanel(
              tabPanel(
                "X title",
                br(),
                textInput("x_text", "X axis label:", "x axis text", placeholder = "V genes"),
                checkboxInput("apply_x2y", "Apply X axis settings to Y axis"),
                br(),
                sliderInput("x_title_size", "X axis title text size:",
                  min = 1, max = 40, value = 16, step = .5
                ),
                sliderInput("x_title_hjust", "X axis title text horizontal adjustment:",
                  min = 0, max = 1, value = 0.5, step = .05
                ),
                sliderInput("x_title_vjust", "X axis title text vertical adjustment:",
                  min = -4, max = 4, value = .5, step = .25
                ),
                sliderInput("x_title_angle", "X axis title text angle:",
                  min = 0, max = 90, value = 0, step = 1
                ),
                selectInput(
                  "x_title_face", "Face:",
                  list(Plain = "plain", Bold = "bold", Italic = "italic", "Bold Italic" = "bold.italic")
                )
              ),
              tabPanel(
                "X ticks",
                br(),
                sliderInput("x_text_size", "X axis text size:",
                  min = 1, max = 40, value = 11, step = .5
                ),
                sliderInput("x_text_hjust", "X axis text horizontal adjustment:",
                  min = -2, max = 2, value = .5, step = .1
                ),
                sliderInput("x_text_vjust", "X axis text vertical adjustment:",
                  min = -4, max = 4, value = .5, step = .25
                ),
                sliderInput("x_text_angle", "X axis text angle:",
                  min = 0, max = 90, value = 90, step = 1
                ),
                selectInput(
                  "x_text_face", "Face:",
                  list(Plain = "plain", Bold = "bold", Italic = "italic", "Bold Italic" = "bold.italic")
                )
              )
            )
          ),
          tabPanel(
            "Y axis",
            br(),
            tabsetPanel(
              tabPanel(
                "Y title",
                br(),
                textInput("y_text", "Y axis label:", "y axis text", placeholder = "Gene frequency"),
                checkboxInput("apply_y2x", "Apply Y axis settings to X axis"),
                br(),
                sliderInput("y_title_size", "Y axis title text size:",
                  min = 1, max = 40, value = 16, step = .5
                ),
                sliderInput("y_title_hjust", "Y axis title text horizontal adjustment:",
                  min = 0, max = 1, value = 0.5, step = .05
                ),
                sliderInput("y_title_vjust", "Y axis title text vertical adjustment:",
                  min = -4, max = 4, value = .5, step = .25
                ),
                sliderInput("y_title_angle", "Y axis title text angle:",
                  min = 0, max = 90, value = 90, step = 1
                ),
                selectInput(
                  "y_title_face", "Face:",
                  list(Plain = "plain", Bold = "bold", Italic = "italic", "Bold Italic" = "bold.italic")
                )
              ),
              tabPanel(
                "Y ticks",
                br(),
                sliderInput("y_text_size", "Y axis text size:",
                  min = 1, max = 40, value = 11, step = .5
                ),
                sliderInput("y_text_hjust", "Y axis text horizontal adjustment:",
                  min = -2, max = 2, value = .5, step = .1
                ),
                sliderInput("y_text_vjust", "Y axis text vertical adjustment:",
                  min = -4, max = 4, value = .5, step = .25
                ),
                sliderInput("y_text_angle", "Y axis text angle:",
                  min = 0, max = 90, value = 0, step = 1
                ),
                selectInput(
                  "y_text_face", "Face:",
                  list(Plain = "plain", Bold = "bold", Italic = "italic", "Bold Italic" = "bold.italic")
                )
              )
            )
          )
        )
      ),
      mainPanel(
        uiOutput("main_plot", style = "position:fixed;")
      )
    )
  )


  #
  #### SERVER ####
  #
  server <- function(input, output, session) {
    create_plot <- function(input) {
      # TODO: make automatic detection of available themes from ggplot2 and other packages
      choose_theme <- function(theme_label) {
        switch(theme_label,
          Linedraw = theme_linedraw(),
          `Black-white` = theme_bw(),
          `Grey / gray` = theme_gray(),
          `Light` = theme_light(),
          `Dark` = theme_dark(),
          `Minimal` = theme_minimal(),
          `Classic` = theme_classic(),
          `Pubr` = theme_pubr(),
          `Pubr w/o lines` = theme_pubr()
        )
      }

      check_empty_str <- function(.str) {
        if (.str == "" || .str == "\n" || .str == "\t") {
          NULL
        } else {
          .str
        }
      }

      get_legend_params <- function(.input, .label) {
        .get <- function(.l) {
          .input[[stringr::str_c("legend", .label, .l, sep = "_")]]
        }

        .remove <- .get("title_remove")

        if (.remove) {
          FALSE # return FALSE to pay respects (please remove it in the next release)
        } else {
          guide_fun <- guide_legend
          if (.get("title_contin")) {
            guide_fun <- guide_colorbar
          }

          guide_fun(
            title = .get("title_text"),
            ncol = .get("title_ncol"),
            title.theme = element_text(
              size = .get("title_size"),
              hjust = .get("title_hjust"),
              vjust = .get("title_vjust"),
              angle = .get("title_angle"),
              face = .get("title_face")
            ),
            label.theme = element_text(
              size = .get("text_size"),
              hjust = .get("text_hjust"),
              vjust = .get("text_vjust"),
              angle = .get("text_angle"),
              face = .get("text_face")
            )
          )
        }
      }

      .plot <- .plot +
        choose_theme(input$ggplot_theme) +
        labs(
          x = check_empty_str(input$x_text),
          y = check_empty_str(input$y_text),
          title = check_empty_str(input$title_text),
          subtitle = check_empty_str(input$subtitle_text),
          fill = input$legend_text,
          color = input$legend_text
        ) +
        guides(
          col = get_legend_params(input, "col"),
          fill = get_legend_params(input, "fill"),
          size = get_legend_params(input, "size"),
          shape = get_legend_params(input, "shape"),
          linetype = get_legend_params(input, "linetype")
        ) +
        theme(
          plot.title = element_text(
            size = input$title_text_size,
            hjust = input$title_text_hjust,
            vjust = input$title_text_vjust,
            angle = input$title_text_angle,
            face = input$title_face
          ),
          plot.subtitle = element_text(
            size = input$subtitle_text_size,
            hjust = input$subtitle_text_hjust,
            vjust = input$subtitle_text_vjust,
            angle = input$subtitle_text_angle,
            face = input$subtitle_face
          ),
          legend.position = input$legend_position,
          legend.box = input$legend_box,
          axis.title.x = element_text(
            size = input$x_title_size,
            hjust = input$x_title_hjust,
            vjust = input$x_title_vjust,
            angle = input$x_title_angle,
            face = input$x_title_face
          ),
          axis.title.y = element_text(
            size = input$y_title_size,
            hjust = input$y_title_hjust,
            vjust = input$y_title_vjust,
            angle = input$y_title_angle,
            face = input$y_title_face
          ),
          axis.text.x = element_text(
            size = input$x_text_size,
            hjust = input$x_text_hjust,
            vjust = input$x_text_vjust,
            angle = input$x_text_angle,
            face = input$x_text_face
          ),
          axis.text.y = element_text(
            size = input$y_text_size,
            hjust = input$y_text_hjust,
            vjust = input$y_text_vjust,
            angle = input$y_text_angle,
            face = input$y_text_face
          )
        )
      if (input$ggplot_theme == "Pubr") {
        .plot <- .plot + theme_cleveland2()
      }

      if (input$coord_flip) {
        .plot <- .plot + coord_flip()
      }

      .plot
    }

    output$save_text <- renderText({
      'To save the plot, press the "Save" button above or drag-n-drop
      the plot to your Desktop or into any file manager (Finder, File Explorer, etc.)'
    })
    output$save_text2 <- renderText({
      'Note: saving via the "Save" button will be different from the drag-n-drop method
      due to R\'s peculiar properties.'
    })

    output$main_plot <- renderUI({
      # if (input$do_interactive) {
      # output$main_plot_helper = renderPlotly(ggplotly(create_plot(input)))
      # plotlyOutput("main_plot_helper")
      # } else {
      output$main_plot_helper <- renderPlot(create_plot(input))
      plotOutput("main_plot_helper", width = input$plot_width * 72, height = input$plot_height * 72)
      # }
    })

    #
    # Assign X settings to Y
    #
    observe({
      if (!is.null(input$apply_x2y)) {
        if (input$apply_x2y) {
          updateSliderInput(session, "y_title_size", value = input$x_title_size)
          updateSliderInput(session, "y_title_hjust", value = input$x_title_size)
          updateSliderInput(session, "y_title_vjust", value = input$x_title_size)
          updateSliderInput(session, "y_title_angle", value = input$x_title_size)

          updateSliderInput(session, "y_text_size", value = input$x_text_size)
          updateSliderInput(session, "y_text_hjust", value = input$x_text_size)
          updateSliderInput(session, "y_text_vjust", value = input$x_text_size)
          updateSliderInput(session, "y_text_angle", value = input$x_text_size)
        }
      }
    })

    #
    # Vice versa: assign Y settings to X
    #
    observe({
      if (!is.null(input$apply_y2x)) {
        if (input$apply_y2x) {
          updateSliderInput(session, "x_title_size", value = input$y_title_size)
          updateSliderInput(session, "x_title_hjust", value = input$y_title_hjust)
          updateSliderInput(session, "x_title_vjust", value = input$y_title_vjust)
          updateSliderInput(session, "x_title_angle", value = input$y_title_angle)

          updateSliderInput(session, "x_text_size", value = input$y_text_size)
          updateSliderInput(session, "x_text_hjust", value = input$y_text_hjust)
          updateSliderInput(session, "x_text_vjust", value = input$y_text_vjust)
          updateSliderInput(session, "x_text_angle", value = input$y_text_angle)
        }
      }
    })

    observeEvent(input$console_plot, {
      plot(create_plot(input))
    })

    #
    # Save plots
    #
    output$save_plot <- downloadHandler(
      filename = paste0("plot shiny ", Sys.time(), ".png"),
      content = function(file) {
        ggsave(file, plot = create_plot(input), width = input$plot_width, height = input$plot_height, device = "png")
      }
    )
  }

  shiny::shinyApp(ui = ui, server = server)
}
