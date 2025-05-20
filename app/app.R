# Filippo Gastaldello, Fabio Mazza, Alessandro Romanel - 13/05/25

library(shiny)
library(bslib)
library(tidyverse)
library(DT)
library(hrbrthemes)
library(plotly)
library(viridis)
# LOAD DATA ----
data <- read_tsv("../data/fake_data.tsv")

# UI ----
ui <- page_sidebar(
        fillable = FALSE,
        title = "PLM DB",
        ## Sidebar ----
        sidebar = sidebar(
                    h4("Search options"),
                    ### Search type selector ----
                    radioButtons(
                        inputId = "search_type",
                        label = "Search by",
                        choices = list(
                            "Gene" = "gene",
                            "Transcript" = "transcript",
                            "Variant" = "variant"
                        ),
                        selected = "gene"
                    ),
                    ### ID text input ----
                    textInput(
                        inputId = "search_id",
                        label = "ID",
                        placeholder = "Insert ID here"
                    ),
                    ### Transcript warning switch ----
                    input_switch(
                        id = "transcripts_warning",
                        label = "Only use validated transcripts"
                    )
        ),
        ## Main content ----
        
        ### Value box ----
        layout_column_wrap(
            value_box( 
                title = "Selected haplotypes",
                showcase = icon("user", "fa-regular"),
                textOutput(outputId = "selected_haplotypes"), 
                theme = "text-blue" 
            ), 
            value_box( 
                title = "Affected transcripts",
                showcase = icon("dna", "fa-regular"), 
                textOutput(outputId = "affected_transcripts"), 
                theme = "text-blue" 
            ), 
            value_box( 
                title = "Involved variants",
                showcase = icon("disease", "fa-regular"),
                textOutput(outputId = "involved_variants"), 
                theme = "text-blue" 
            ) 
        ),
        
        ### datatable ----
        card(
            card_header("Haplotypes table"),
            dataTableOutput("datatable"),
            full_screen = TRUE
        ),
        
        layout_column_wrap(
            ### Score distributions ----
            card(
                card_header("Score distribution"),
                layout_sidebar(
                    sidebar = sidebar(
                        selectizeInput(
                            inputId = "model",
                            label = "Select the desired model:",
                            choices = list("ESMv2"="score_esm2",
                                           "PoET"="score_poet",
                                           "HAL9000"="score_hal9000",
                                           "Jarvis"="score_jarvis")
                        )
                    ),
                    plotlyOutput(outputId = "score_distribution")
                ),
                full_screen = TRUE,
                fill = TRUE
            ),
            
            ### Score deltas distributions ----
            card(
                card_header("Score deltas distribution"),
                layout_sidebar(
                    sidebar = sidebar(
                        selectizeInput(
                            inputId = "delta",
                            label = "Select the desired model's delta:",
                            choices = list("ESMv2"="delta_esm",
                                           "PoET"="delta_poet",
                                           "HAL9000"="delta_hal9000",
                                           "Jarvis"="delta_jarvis")
                        )
                    ),
                    plotlyOutput(outputId = "score_delta_distribution")
                ),
                full_screen = TRUE,
                fill = TRUE
            )
        ),
        card(
            card_header("Ancestry"),
            layout_sidebar(
                sidebar = sidebar(
                    radioButtons(
                        inputId = "group_by",
                        label = "Group by",
                        choices = list(
                            "Haplotype" = "haplotype",
                            "Ancestry" = "ancestry"
                        ),
                        selected = "haplotype"
                    )
                ),
                plotlyOutput("population")
            ),
            full_screen = TRUE
        )
)

# SERVER LOGIC ----
server <- function(input, output){
    
    ## Update ID text input ----
    # Observe changes in the radioButtons
    observeEvent(input$search_type, {
        selected_type <- input$search_type
        new_label <- ""
        
        if (selected_type == "gene") {
            new_label <- "Ensembl Gene ID:"
        } else if (selected_type == "transcript") {
            new_label <- "Ensembl Transcript ID:"
        } else if (selected_type == "variant") {
            new_label <- "Variant ID (chr:pos):"
        }
        
        # Update the label of the textInput
        updateTextInput(
            inputId = "search_id",
            label = new_label
        )
    })
    
    ## Update value boxes ----
    output$selected_haplotypes <- renderText({
        nrow(df() %>% select(haplo_id) %>% unique())
    })
    output$affected_transcripts <- renderText({
        nrow(df() %>% select(transcript_id) %>% unique())
    })
        
    output$involved_variants <- renderText({
        subset <- df() %>% select(dna_changes) %>% filter(!dna_changes == "wt")
        variants <- c()
        for(row in subset$dna_changes){
            variants <- c(variants, str_split_1(row, pattern = ","))
        }
        length(unique(variants))
    })
    ## Subset dataframe ---- 
    df <- reactive({
        if (input$search_type == "gene") {
            df <- data %>% filter(gene_id == input$search_id)
        }
        if (input$search_type == "transcript") {
            df <- data %>% filter(transcript_id == input$search_id)
        }
        if (input$search_type == "variant") {
            df <- data %>% filter(input$search_id %in% dna_changes)
        }
        if (input$transcripts_warning) {
            df <- df %>% filter(warning == TRUE)
        }
        return(df)
    })
    
    ## Datatable ----
    output$datatable <- renderDataTable({
        df()
    })
    
    ## Score distribution plots ----
    output$score_distribution <- renderPlotly({
        if (nrow(df()>0)) {
            p <- data %>%
                ggplot() +
                geom_density(aes(x = !!sym(input$model)),
                             fill = "#e6ab47",
                             color = "#e39107",
                             alpha = 0.8,
                             adjust = 0.1) +
                geom_vline(data = df(), 
                           aes(xintercept = !!sym(input$model), text = paste0(haplo_id, " ", input$model, ": ", !!sym(input$model))),
                           color = "#e34907") +
                theme(legend.position = "none")
            ggplotly(p, tooltip = "text")
        }else{
            p <- data %>% ggplot(aes(x = !!sym(input$model))) +
                        geom_density(fill = "#e6ab47", color = "#e39107", alpha = 0.8, adjust = 0.1)
            ggplotly(p)
        }
    })
    
    ## Score deltas distribution plots ----
    output$score_delta_distribution <- renderPlotly({
        if (nrow(df()>0)) {
            p <- data %>%
                ggplot() +
                geom_density(aes(x = !!sym(input$delta)),
                             fill = "#e6ab47",
                             color = "#e39107",
                             alpha = 0.8,
                             adjust = 0.1) +
                geom_vline(data = df(), 
                           aes(xintercept = !!sym(input$delta), text = paste0(haplo_id, " ", input$delta, ": ", !!sym(input$delta))),
                           color = "#e34907") +
                theme(legend.position = "none")
            ggplotly(p, tooltip = "text")
        }else{
            p <- data %>% ggplot(aes(x = !!sym(input$delta))) +
                geom_density(fill = "#e6ab47", color = "#e39107", alpha = 0.8, adjust = 0.1)
            ggplotly(p)
        }
    })
    
    ## Population frequencies ----
    output$population <- renderPlotly({
        validate(need(try(nrow(df())>0),"Waiting for a subset of haplotypes to be selected."))
        if (input$group_by == "ancestry") {
            # Group by ancestry
            p <- df() %>%
                    select(haplo_id, AFR_freq, AMR_freq, EAS_freq, EUR_freq, SAS_freq) %>%
                    pivot_longer(cols = c("AFR_freq", "AMR_freq", "EAS_freq", "EUR_freq", "SAS_freq"),
                                 names_to = "ancestry",
                                 values_to = "freq") %>% 
                    ggplot(aes(x = ancestry, y = freq, fill = haplo_id)) +
                        geom_bar(position = "dodge", stat = "identity") +
                        scale_fill_viridis(discrete = TRUE) +
                        theme(legend.position = "none")
            ggplotly(p)
        }else{
            # Group by haplotype
            p <- df() %>%
                    select(haplo_id, AFR_freq, AMR_freq, EAS_freq, EUR_freq, SAS_freq) %>%
                    pivot_longer(cols = c("AFR_freq", "AMR_freq", "EAS_freq", "EUR_freq", "SAS_freq"),
                                 names_to = "ancestry",
                                 values_to = "freq") %>% 
                    ggplot(aes(x = haplo_id, y = freq, fill = ancestry)) +
                    geom_bar(position = "dodge", stat = "identity") +
                    scale_fill_viridis(discrete = TRUE) +
                    theme(legend.position = "none")
            ggplotly(p)
        }
    })
    
}
 

# RUN APP ----
shinyApp(ui, server)
