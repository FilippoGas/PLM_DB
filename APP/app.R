# Filippo Gastaldello, Fabio Mazza, Alessandro Romanel - 13/05/25

library(shiny)
library(bslib)

# UI ----
ui <- page_sidebar(
        title = "NAR DB",
        ### Sidebar ----
        sidebar = sidebar(
                    h4("Search options"),
                    #### Search type selector ----
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
                    ),
                    #### Model selection ----
                    selectizeInput(
                        inputId = "model",
                        label = "Select the desired model:",
                        choices = list("ESMv2"="esmv2",
                                       "PoET"="poet",
                                       "HAL9000"="hal9000",
                                       "Jarvis"="jarvis")
                    )
        ),
        "Main content",
        card(
            card_header("Card"),
        ) 
)

# SERVER LOGIC ----
server <- function(input, output){
    
    #### Update ID text input ----
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
}

# RUN APP ----
shinyApp(ui, server)