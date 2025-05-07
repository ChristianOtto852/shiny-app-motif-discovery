library(shiny)
source("required-libraries.R")
source_python("greedy-motif-search.py")
source("example_data_sets.R")
source_python("randomized-motif-search.py")
source_python("gibbs-sampler.py")
source("bcrank-function.R")

function(input, output, session) {
    
    observe({
        disable(selector = "input[type='radio'][value='custom_data']")
    })
    
    observeEvent(input$custom_data, {
        enable(selector = "input[type='radio'][value='custom_data']")
    })
    
    
    dna_sequences <- reactive({
        req(input$example)
        
        if (input$example == "custom_data") {
            req(input$custom_data)
            readLines(input$custom_data$datapath)
        } else if (input$example == "real_data") {
            readLines("real_data_set.txt")
        }
        else {
            switch(input$example,
                   "example_1" = example_1,
                   "example_2" = example_2,
                   "example_3" = example_3,
                   "example_4" = example_4,
                   "example_5" = example_5)
        }
    })
    
    observeEvent(dna_sequences(), {
        seqs <- unlist(dna_sequences())
        max_k <- nchar(seqs[1])
        updateSliderInput(session, "k", max = max_k)
    })
    
    output$method_info <- renderUI({
        req(input$method)
        
        info_links <- list(
            greedy = tags$a(href = "greedy-info.pdf", "Learn about Greedy Search", target = "_blank"),
            randomized = tags$a(href = "randomized-info.pdf", "Learn about Randomized Search", target = "_blank"),
            gibbs = tags$a(href = "gibbs-info.pdf", "Learn about Gibbs Sampling", target = "_blank"),
            bcrank = tags$a(href = "bcrank-info.pdf", "Learn about Bcrank", target = "_blank")
        )
        
        info_links[[input$method]]
    })
    
    output$input_info <- renderUI({
        input_info_link = tags$a(href = "input-info.pdf", "Input Data Information", target = "_blank")
    })
    
    output$real_data_info <- renderUI({
        input_info_link = tags$a(href = "real-data-info.pdf", "Real Data Information", target = "_blank")
    })
    
    
    output$data <- renderText({
        sequences <- dna_sequences()
        top_8 <- head(sequences, 8)
        top_8 <- substr(top_8, 1, 50)
        paste(top_8, collapse = "\n")
    })
    
    observeEvent(input$example, {
        default_k <- switch(input$example,
                            "example_1" = 3,
                            "example_2" = 3,
                            "example_3" = 5,
                            "example_4" = 6,
                            "example_5" = 5,
                            "real_data" = 20,
                            "custom_data" = 3)
        updateNumericInput(session, "k", value = default_k)
    })
    
    motif_result <- eventReactive(input$run, {
        start <- Sys.time()
        motifs <- switch(input$method,
                         "greedy" = greedy_motif_search_pseudocounts(dna_sequences(), input$k),
                         "randomized" = randomized_motif_search(dna_sequences(), input$k),
                         "gibbs" = gibbs_sampler(dna_sequences(), input$k, n = as.integer(100)),
                         "bcrank" = Bcrank_motif_finder(dna_sequences(), input$k)
        )
        end <- Sys.time()
        
        time <- as.numeric(difftime(end, start, units = "secs"))
        
        list(motifs = motifs, runtime = time)
    })
    
    output$result <- renderText({
        motifs_results <- motif_result()
        paste(motifs_results$motifs, collapse = "\n")
        
    })
    
    output$runtime <- renderText({
        motifs_results <- motif_result()
        glue::glue("Runtime: ", round(motifs_results$runtime, 3), " Seconds")
    })
    
    
    motif_all_results <- eventReactive(input$run_all, {
        start_greedy <- Sys.time()
        greedy_all <- greedy_motif_search_pseudocounts(dna_sequences(), input$k)
        end_greedy <- Sys.time()
        start_randomized <- Sys.time()
        randomized_all <- randomized_motif_search(dna_sequences(), input$k)
        end_randomized <- Sys.time()
        start_gibbs <- Sys.time()
        gibbs_all <- gibbs_sampler(dna_sequences(), input$k, n = as.integer(100))
        end_gibbs <- Sys.time()
        start_bcrank <- Sys.time()
        bcrank_all <- Bcrank_motif_finder(dna_sequences(), input$k)
        end_bcrank <- Sys.time()
        
        list(
            time_greedy = round(as.numeric(difftime(end_greedy, start_greedy, units = "secs")), 3),
            time_randomized = round(as.numeric(difftime(end_randomized, start_randomized, units = "secs")), 3),
            time_gibbs = round(as.numeric(difftime(end_gibbs, start_gibbs, units = "secs")), 3),
            time_bcrank = round(as.numeric(difftime(end_bcrank, start_bcrank, units = "secs")), 3),
            motif_greedy = greedy_all,
            motif_randomized = randomized_all,
            motif_gibbs = gibbs_all,
            motif_bcrank = bcrank_all
        )
        
    })
    
    output$all_results <- renderTable({
        all_info <- motif_all_results()
        tibble(
            "Algorithm" = c("Greedy", "Randomized", "Gibbs", "Bcrank"),
            "Time" = c(glue::glue(all_info$time_greedy, " Seconds"), glue::glue(all_info$time_randomized, " Seconds"), glue::glue(all_info$time_gibbs, " Seconds"), glue::glue(all_info$time_bcrank, " Seconds")),
            "Motifs" = c(paste(all_info$motif_greedy, collapse = " "), paste(all_info$motif_randomized, collapse = " "), paste(all_info$motif_gibbs, collapse = " "), paste(all_info$motif_bcrank, collapse = " "))
        )
    })
}


























