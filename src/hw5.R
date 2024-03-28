library(shiny)
library(sf)
library(tidyverse)
library(tigris)
library(tmap)

leukemia <- read_csv("./mn_leukemia.csv") %>%
  mutate(fips = as.character(fips)) %>%
  mutate(my_rate = (count / population)*100000) %>%
  inner_join(counties("Minnesota"), by = c("fips" = "GEOID")) %>%
  st_as_sf()

ui <- fluidPage(
  titlePanel("Leukemia rates per 100,000 in Minnesota"),
  sidebarLayout(
    sidebarPanel(
      helpText("Create maps for Leukemia rates in Minnesota"),
      selectInput("sex", label = "Sex:", choices = c("Female", "Male"), selected = "Female"),
      selectInput("year", label = "Year Range:", choices = c("2005-2009", "2010-2014", "2015-2019"), selected = ("2015-2019"))
    ),
    mainPanel(plotOutput("map"))
  )
)

server <- function(input, output) {
  output$map <- renderPlot({
    leukemia %>%
      filter(year == input$year & sex == input$sex) %>%
      tm_shape() +
      tm_polygons(col = "my_rate", title = "Rate per 100,000") +
      tm_layout(paste0("Leukemia rates for ", input$sex, "s during ", input$year), inner.margins = c(0.025, 0.025, 0.075, 0.1), legend.position = c("right", "bottom"))
  })
}

shinyApp(ui, server)
