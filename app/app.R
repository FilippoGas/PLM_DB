# Filippo Gastaldello, Fabio Mazza, Alessandro Romanel - 13/05/25

library(shiny)
library(bslib)
library(tidyverse)
library(DT)
library(hrbrthemes)
library(plotly)
library(viridis)
library(ggpubr)
library(scales)

# LOAD DATA ----
data <- read_tsv("../data/interface_data.tsv")

rsid_map <- read_tsv("../data/rsid_map.tsv")
rsid_map <- rsid_map %>% column_to_rownames("rsid") %>% mutate(haplo_ids = strsplit(haplo_ids, ","))

coords_map <- read_tsv("../data/coords_map.tsv")
coords_map <- coords_map %>% column_to_rownames("variant_coord") %>% mutate(haplo_ids = strsplit(haplo_ids, ","))

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
                            width = "20%",
                            h3("Search options"),
                            #### Search type selector ----
                            radioButtons(
                                inputId = "search_type",
                                label = h5("Search by"),
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
                            accordion(
                                accordion_panel(
                                    title = "Advanced search",
                                    #### Variant type filter ----
                                    checkboxGroupInput(
                                        inputId = "variant_type_filter",
                                        label = h5("Exclude haplotypes with:"),
                                        choices = list("Synonymous variants" = "synonymous_variant",
                                                       "Missense variants" = "missense_variant",
                                                       "Nonsense variants" = "stop_gained",
                                                       "Frameshift variants" = "framesift_variants",
                                                       "In frame deletions" = "inframe_deletion",
                                                       "In frame insertions" = "inframe_insertion",
                                                       "Start loss" = "start_lost",
                                                       "Start gain" = "initiator_codon_variant",
                                                       "Stop loss" = "stop_lost")
                                    ),
                                    #### Ensembl tsl filter ----
                                    sliderInput(
                                        inputId = "tsl_filter",
                                        label = h5("Select desired Ensembl transcript support level (tsl):"),
                                        min = 1,
                                        max = 5,
                                        value = 5
                                    ),
                                    checkboxInput(inputId = "tsl_NA",
                                                  label = "Remove NAs",
                                                  value = FALSE
                                    )
                                ),
                                id = "advanced_search",
                                open = FALSE
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
                    ),
                    max_height = "10%"
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
                                    choices = list("ESMv2"=list("PLL" = "esm_PLL",
                                                                "PLLR max freq." = "esm_PLLR_maxfreq",
                                                                "PLLR wt" = "esm_PLLR_ref")
                                                   )
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
                                    choices = list("ESMv2 Delta" = "delta_esm_PLL"
                                                   )
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
                #### Ancestry ----
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
                            ),
                            selectizeInput(
                                inputId = "ancestry_frequency_filter",
                                label = "Restrict to:",
                                choice = NULL,
                                options = list(placeholder = "Select an option")
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
            title = "About",
            ### Value box ----
            layout_column_wrap(
                value_box( 
                    title = "Genes",
                    textOutput(outputId = "analyzed_genes"), 
                    theme = "text-blue" 
                ), 
                value_box( 
                    title = "Transcripts",
                    textOutput(outputId = "analyzed_transcripts"), 
                    theme = "text-blue" 
                ), 
                value_box( 
                    title = "Haplotypes",
                    textOutput(outputId = "analyzed_haplotypes"), 
                    theme = "text-blue" 
                ),
                value_box( 
                    title = "Variants",
                    textOutput(outputId = "analyzed_variants"), 
                    theme = "text-blue" 
                ),
                max_height = "10%"
            ),
            layout_column_wrap(
                ### Haplotypes per gene distributions ----
                card(
                    card_header("Haplotypes per gene distribution"),
                    plotlyOutput(outputId = "haplotype_gene_distribution"),
                    full_screen = TRUE,
                    fill = TRUE
                ),
                ### Variants per haplotype distributions ----
                card(
                    card_header("Variants per haplotype distribution"),
                    plotlyOutput(outputId = "variant_haplotype_distribution"),
                    full_screen = TRUE,
                    fill = TRUE
                ),
                max_height = "40%"
            ),
            ### Models correlation ----
            card(
                card_header("Model scores correlation"),
                plotlyOutput(outputId = "models_correlation"),
                full_screen = TRUE,
                fill = TRUE 
            )
        ),
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
    
    ## Search ----
    
    ### Update ID text input ----
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
    
    ### Subset dataframe ---- 
    df <- reactive({
        if (input$search_type == "ensg") {
            df <- data %>% filter(gene_id == input$search_id)
        }
        if (input$search_type == "hgnc") {
            df <- data %>% filter(gene_symbol == input$search_id)
        }
        if (input$search_type == "transcript") {
            df <- data %>% filter(transcript_id == input$search_id)
        }
        if (input$search_type == "coordinates") {
            haplotypes <- coords_map[input$search_id,"haplo_ids"][[1]]
            df <- data %>% filter(haplo_id %in% haplotypes)
        }
        if (input$search_type == "rsid") {
            haplotypes <- rsid_map[input$search_id,"haplo_ids"][[1]]
            df <- data %>% filter(haplo_id %in% haplotypes)
        }
        # Check advanced filters
        if (length(input$variant_type_filter)>0) {
            try(
                df <- df %>% filter(!grepl(paste0(input$variant_type_filter, collapse = "|"), df$variant_types))
            )
        }
        if (input$tsl_filter) {
            
        }
        
        return(df)
    })
    
    ### Update value boxes ----
    output$selected_haplotypes <- renderText({
        nrow(df() %>% select(haplo_id) %>% unique())
    })
    output$affected_transcripts <- renderText({
        nrow(df() %>% select(transcript_id) %>% unique())
    })
        
    output$involved_variants <- renderText({
        subset <- df() %>% select(variant_coordinates) %>% filter(!variant_coordinates == "wt")
        variants <- c()
        for(row in subset$variant_coordinates){
            variants <- c(variants, str_split_1(row, pattern = ","))
        }
        length(unique(variants))
    })
    
    ### Datatable ----
    output$datatable <- renderDataTable({
        df()
    })
    
    ### Update distribution filters, both for scores and deltas----
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
    
    ### Score distribution plots ----
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
            if (str_detect(input$model, "PLLR")) {
                p <- data %>% filter(!input$model==0) %>% ggplot() +
                                                            geom_density(
                                                                aes(x = !!sym(input$model)),
                                                                fill = "#e6ab47",
                                                                color = "#e39107",
                                                                alpha = 0.8,
                                                                adjust = 0.5,
                                                                stat = "count") +
                                                            geom_vline(
                                                                data = plot_df, 
                                                                aes(xintercept = !!sym(input$model),
                                                                    text = paste0(haplo_id,
                                                                                  " on transcript ",
                                                                                  transcript_id,
                                                                                  "\n",
                                                                                  input$model,
                                                                                  ": ",
                                                                                  round(!!sym(input$model),4)
                                                                                  )
                                                                    ),
                                                                color = "#e34907") +
                                                            theme(legend.position = "none") +
                                                            scale_y_log10(name = "Density (log10 scale)",
                                                                          breaks = trans_breaks("log10", function(x) 10^x),
                                                                          labels = trans_format("log10", math_format(10^.x)
                                                                                                )
                                                                          )
            }else{
                p <- data %>% ggplot() +
                                geom_density(
                                    aes(x = !!sym(input$model)),
                                    fill = "#e6ab47",
                                    color = "#e39107",
                                    alpha = 0.8,
                                    adjust = 0.1) +
                                geom_vline(
                                    data = plot_df, 
                                    aes(xintercept = !!sym(input$model),
                                        text = paste0(haplo_id,
                                                      " on transcript ",
                                                      transcript_id,
                                                      "\n",
                                                      input$model,
                                                      ": ",
                                                      round(!!sym(input$model),4)
                                        )
                                    ),
                                    color = "#e34907") +
                                theme(legend.position = "none")
            }
            ggplotly(p, tooltip = "text")
        }else{
            if (str_detect(input$model, "PLLR")) {
                p <- data %>% filter(!input$model==0) %>% ggplot(aes(x = !!sym(input$model))) +
                                                            geom_density(fill = "#e6ab47",
                                                                         color = "#e39107",
                                                                         alpha = 0.8,
                                                                         adjust = 0.5,
                                                                         stat = "count"
                                                            ) +
                                                            scale_y_log10(name = "Density (log10 scale)",
                                                                          breaks = trans_breaks("log10", function(x) 10^x),
                                                                          labels = trans_format("log10", math_format(10^.x)
                                                                                                )
                                                                          )
            }else{
                p <- data %>% ggplot(aes(x = !!sym(input$model))) +
                                geom_density(fill = "#e6ab47",
                                             color = "#e39107",
                                             alpha = 0.8,
                                             adjust = 0.1
                                )
            }
            ggplotly(p)
        }
    })
    
    ### Score deltas distribution plots ----
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
                           aes(xintercept = !!sym(input$delta),
                               text = paste0(haplo_id, " on transcript ", transcript_id, "\n", input$delta, ": ", round(!!sym(input$delta),4))),
                           color = "#e34907") +
                theme(legend.position = "none")
            ggplotly(p, tooltip = "text")
        }else{
            p <- data %>% ggplot(aes(x = !!sym(input$delta))) +
                geom_density(fill = "#e6ab47", color = "#e39107", alpha = 0.8, adjust = 0.1)
            ggplotly(p)
        }
    })
    
    ### Update ancestry plot filters ----
    observe({
        # Extract haplotype and trancsript ID from subset df 
        current_df_val <- df()
        transcripts <- c("none", sort(unique(current_df_val$transcript_id)))
        updateSelectizeInput(session,
                             inputId = "ancestry_frequency_filter",
                             label = "Restrict to:",
                             choices = list(
                                 "Transcripts" = as.list(transcripts)
                             ),
                             selected = as.list(transcripts)[2], # the first would be "none"
                             server = FALSE)
        
    })
    
    ### Population frequencies ----
    output$population <- renderPlotly({
        validate(need(try(nrow(df())>0),"Waiting for a subset of haplotypes to be selected."))
        plot_df <- df()
        # Check if filtering options have been applied for the plot
        if (str_detect(input$ancestry_frequency_filter, "ENST")) {
            plot_df <- plot_df %>% filter(transcript_id == input$ancestry_frequency_filter)
        }
        if (input$group_by == "ancestry") {
            # Group by ancestry
            p <- plot_df %>%
                    select(frequency, haplo_id, AFR_freq, AMR_freq, EAS_freq, EUR_freq, SAS_freq) %>%
                    dplyr::rename("Global" = "frequency") %>%
                    pivot_longer(cols = c("Global", "AFR_freq", "AMR_freq", "EAS_freq", "EUR_freq", "SAS_freq"),
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
            p <- plot_df %>%
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
    
    ## About ----
    
    ### Update value boxes ----
    output$analyzed_genes <- renderText({
        nrow(data %>% select(gene_id) %>% unique())
    })
    output$analyzed_haplotypes <- renderText({
        nrow(data %>% select(haplo_id) %>% unique())
    })
    output$analyzed_transcripts <- renderText({
        nrow(data %>% select(transcript_id) %>% unique())
    })
    # Read from stats file
    output$analyzed_variants <- renderText({
        subset <- data %>% select(variant_coordinates) %>% filter(!variant_coordinates == "wt") %>% unique()
        variants <- c()
        for(row in subset$variant_coordinates){
            variants <- c(variants, str_split_1(row, pattern = ","))
        }
        length(unique(variants))
    })
    
    ### haplotypes per gene distribution ----
    output$haplotype_gene_distribution <- renderPlotly({
        plot_df <- data %>%
                        select(haplo_id, gene_id) %>%
                        unique() %>%
                        summarise(haplotype_count = n(), .by = gene_id)
        p <- plot_df %>% ggplot(aes(x = haplotype_count)) +
                            geom_histogram(stat = "count",
                                           binwidth = 2,
                                           fill = "#e6ab47",
                                           color = "#e39107",
                                           alpha = 0.8)
        ggplotly(p)
    })
    ### variants per haplotype distribution ----
    output$variant_haplotype_distribution <- renderPlotly({
        p <- data %>% ggplot(aes(x = n_variants)) +
            geom_histogram(stat = "count",
                           binwidth = 2,
                           fill = "#e6ab47",
                           color = "#e39107",
                           alpha = 0.8)
        ggplotly(p)
    })
    
    ### Model correlation ----
    output$models_correlation <- renderPlotly({
        score_columns <- colnames(data)[grep("score", colnames(data))]
        couples <- expand.grid(score_columns, score_columns) %>%
            filter(!Var1==Var2) %>%
            arrange(Var1)
        couples <- couples %>% slice_head(n = nrow(couples)/2) %>% mutate(collapsed = paste0(Var1, "-", Var2))
        plots <- lapply(couples$collapsed, function(x){
            var1 <- str_split_i(x, "-", 1)
            var2 <- str_split_i(x, "-", 2)
            data %>% ggplot(aes(x = !!sym(var1), y = !!sym(var2))) +
                        geom_point(fill = "#e6ab47",
                                   color = "#e39107") +
                        ggtitle(paste0("Correlaion between ", var1, " and ", var2)) +
                        stat_cor(method = "pearson")
        })
        
        # p <- do.call("ggarrange", c(plots, list("ncol"=length(plots), "nrow"=1)))
        do.call("subplot", c(plots, list = ("nrows"=1)))
        
        
    })
    
    ## Download ----
    
    ## FAQ ----
    
    ## Contacts ----
    
}
 

# RUN APP ----
shinyApp(ui, server)
