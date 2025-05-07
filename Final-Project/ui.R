

library(shiny)
source("required-libraries.R")
source_python("greedy-motif-search.py")
source("example_data_sets.R")
source_python("randomized-motif-search.py")
source_python("gibbs-sampler.py")
source("bcrank-function.R")


fluidPage(
    useShinyjs(),
    # Application title
    titlePanel("Motif Finding"),

    sidebarLayout(
        sidebarPanel(
            
            # Custom Data Import
            fileInput("custom_data", "Upload DNA Sequences (.txt)",
                      multiple = FALSE,
                      accept = c(".txt")
                      ),
            uiOutput("input_info"),
            
            # Algorithm Choice
            radioButtons("method", "Choose Motif Finder:",
                         choices = c("Greedy" = "greedy",
                                     "Randomized" = "randomized",
                                     "Gibbs" = "gibbs",
                                     "Bcrank" = "bcrank")
                         ),
            # Hyperlink
            uiOutput("method_info"),
            radioButtons("example", "Choose Example Dataset:",
                         choices = c("Example 1" = "example_1",
                                     "Example 2" = "example_2",
                                     "Example 3" = "example_3",
                                     "Example 4" = "example_4",
                                     "Example 5" = "example_5",
                                     "Real Data Example" = "real_data",
                                     "Uploaded Data" = "custom_data"
                                     )
                        ),
            uiOutput("real_data_info"),
            # Picking k
            sliderInput("k", "Motif length", min = 1, step = 1, max = 12, value = 3),
            
            # Buttons for running a single algorithm or all
            actionButton("run", "Run Algorithm"),
            actionButton("run_all", "Run All Algorithms")
        ),
        
        # Right Panel with Output
        mainPanel(
            h3("Data"),
            verbatimTextOutput("data"),
            h3("Estimated Motifs"),
            verbatimTextOutput("runtime"),
            verbatimTextOutput("result") |> withSpinner(),
            h3("Comparing Algorithms"),
            tableOutput("all_results") |> withSpinner()
        )
    )
)
