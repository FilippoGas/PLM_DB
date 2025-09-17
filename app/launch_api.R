library(plumber)

pr <- plumb("/srv/shiny-server/api.R")
pr$run(port=8007, host="127.0.0.1", swagger = TRUE, debug = TRUE)

