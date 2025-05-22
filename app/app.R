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
ui <- page_navbar(
        title = "PLM DB",
        nav_spacer(),
        ## Search ----
        nav_panel(
            title = "Search",
            layout_sidebar(
                fillable = FALSE,
                ### Sidebar ----
                sidebar = sidebar(
                            h4("Search options"),
                            #### Search type selector ----
                            radioButtons(
                                inputId = "search_type",
                                label = "Search by",
                                choices = list(
                                    "Ensembl Gene ID" = "ensg",
                                    "Gene Symbol" = "hgnc",
                                    "Ensembl Transcript ID" = "transcript",
                                    "Variant Coordinates" = "coordinates",
                                    "Variant rsid" = "rsid"
                                ),
                                selected = "ensg"
                            ),
                            #### ID text input ----
                            textInput(
                                inputId = "search_id",
                                label = "ID",
                                placeholder = "Insert ID here"
                            ),
                            #### Transcript warning switch ----
                            input_switch(
                                id = "transcripts_warning",
                                label = "Only use validated transcripts"
                            )
                ),
                ### Main content ----
                
                #### Value box ----
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
                
                #### datatable ----
                card(
                    card_header("Haplotypes table"),
                    dataTableOutput("datatable"),
                    full_screen = TRUE
                ),
                
                layout_column_wrap(
                    #### Score distributions ----
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
                                ),
                                selectizeInput(
                                    inputId = "score_distribution_filter",
                                    label = "Restrict to:",
                                    choice = NULL,
                                    options = list(placeholder = "Select an option")
                                )
                            ),
                            plotlyOutput(outputId = "score_distribution")
                        ),
                        full_screen = TRUE,
                        fill = TRUE
                    ),
                    
                    #### Score deltas distributions ----
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
                                ),
                                selectizeInput(
                                    inputId = "delta_distribution_filter",
                                    label = "Restrict to:",
                                    choice = NULL,
                                    options = list(placeholder = "Select an option")
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
        ),
        ## About ----
        nav_panel(
          title = "About"),
        ## Download ----
        nav_panel(
            title = "Download",
        ),
        ## FAQ ----
        nav_panel(
            title = "FAQ"
        ),
        ## Contacts ----
        nav_panel(
            title = "Contacts"
        ),
        nav_spacer()
)

# SERVER LOGIC ----
server <- function(input, output, session){
    
    ## Update ID text input ----
    # Observe changes in the radioButtons
    observeEvent(input$search_type, {
        selected_type <- input$search_type
        new_label <- ""
        
        if (selected_type == "ensg") {
            new_label <- "Ensembl Gene ID:"
        } else if (selected_type == "hgnc") {
            new_label <- "HGNC gene symbol:"
        }else if (selected_type == "transcript") {
            new_label <- "Ensembl Transcript ID:"
        } else if (selected_type == "coordinates") {
            new_label <- "Variant coordinates (chr:pos.REF>ALT):"
        }else if (selected_type == "rsid") {
            new_label <- "Variant rsid:"
        }
        
        # Update the label of the textInput
        updateTextInput(
            inputId = "search_id",
            label = new_label
        )
    })
    
    ## Subset dataframe ---- 
    df <- reactive({
        if (input$search_type == "ensg") {
            df <- data %>% filter(gene_id == input$search_id)
        }
        if (input$search_type == "hgnc") {
            df <- data %>% filter(hgnc_symbol == input$search_id)
        }
        if (input$search_type == "transcript") {
            df <- data %>% filter(transcript_id == input$search_id)
        }
        if (input$search_type == "coordinates") {
            df <- data %>% filter(input$search_id %in% dna_changes)
        }
        if (input$search_type == "rsid") {
            df <- data %>% filter(input$search_id %in% rsid)
        }
        if (input$transcripts_warning) {
            df <- df %>% filter(warning == TRUE)
        }
        return(df)
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
    
    ## Datatable ----
    output$datatable <- renderDataTable({
        df()
    })
    
    ## Update distribution filters, both for scores and deltas----
    observe({
        # Extract haplotype and trancsript ID from subset df 
        current_df_val <- df()
        haplotypes <- c("none", sort(unique(current_df_val$haplo_id)))
        transcripts <- c("none", sort(unique(current_df_val$transcript_id)))
        updateSelectizeInput(session,
                             inputId = "score_distribution_filter",
                             label = "Restrict to:",
                             choices = list(
                                        "Haplotypes" = as.list(haplotypes),
                                        "transcripts" = as.list(transcripts)
                             ),
                             server = FALSE)
        updateSelectizeInput(session,
                             inputId = "delta_distribution_filter",
                             label = "Restrict to:",
                             choices = list(
                                    "Haplotypes" = as.list(haplotypes),
                                    "transcripts" = as.list(transcripts)
                             ),
                             server = FALSE)
        
    })
    
    ## Score distribution plots ----
    output$score_distribution <- renderPlotly({
        if (nrow(df()>0)) {
            plot_df <- df()
            # Check if filtering options have been applied for the plot
            if (str_detect(input$score_distribution_filter, "ENST")) {
                plot_df <- plot_df %>% filter(transcript_id == input$score_distribution_filter)
            }
            if (str_detect(input$score_distribution_filter, "ENSG")) {
                plot_df <- plot_df %>% filter(haplo_id == input$score_distribution_filter)
            }
            p <- data %>%
                ggplot() +
                geom_density(aes(x = !!sym(input$model)),
                             fill = "#e6ab47",
                             color = "#e39107",
                             alpha = 0.8,
                             adjust = 0.1) +
                geom_vline(data = plot_df, 
                           aes(xintercept = !!sym(input$model), text = paste0(haplo_id, " on transcript ",transcript_id, "\n", input$model, ": ", round(!!sym(input$model),2))),
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
            plot_df <- df()
            # Check if filtering options have been applied for the plot
            if (str_detect(input$delta_distribution_filter, "ENST")) {
                plot_df <- plot_df %>% filter(transcript_id == input$delta_distribution_filter)
            }
            if (str_detect(input$delta_distribution_filter, "ENSG")) {
                plot_df <- plot_df %>% filter(haplo_id == input$delta_distribution_filter)
            }
            p <- data %>%
                ggplot() +
                geom_density(aes(x = !!sym(input$delta)),
                             fill = "#e6ab47",
                             color = "#e39107",
                             alpha = 0.8,
                             adjust = 0.1) +
                geom_vline(data = plot_df, 
                           aes(xintercept = !!sym(input$delta), text = paste0(haplo_id, " on transcript ", transcript_id, "\n", input$delta, ": ", round(!!sym(input$delta),2))),
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
                        theme(legend.position = "none",
                              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
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
                    theme(legend.position = "none",
                          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
            ggplotly(p)
        }
    })
    
}
 

# RUN APP ----
shinyApp(ui, server)
