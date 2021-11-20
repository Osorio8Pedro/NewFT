pos=as.numeric(format(Sys.Date(), "%m"))
meses=c("Enero","Febrero","Marzo","Abril","Mayo","Junio","Julio","Agosto",
        "Septiembre","Octubre","Noviembre","Diciembre")
require("shinydashboard")
require("dplyr")
ui <- dashboardPage(
  dashboardHeader(title=paste("Datos UCI ",format(Sys.Date(), "%d"),meses[pos])),
  dashboardSidebar(selectInput("Region", "Seleccionar Region",
                               choices = c("Nacional","Arica","Tarapaca",
                                           "Antofagasta","Atacama","Coquimbo",
                                           "Valparaiso","R.Metropolitana",
                                           "O'Higgins","Maule","Nuble","Bio-Bio",
                                           "Araucania","Los Rios","Los Lagos",
                                           "Aysen","Magallanes","Norte Grande",
                                           "Norte Chico","Zona Central","Zona Sur",
                                           "Zona Austral","Regiones"),
                               selected = "Nacional"),
  sliderInput(inputId = "date",label = "Rango Fechas:",
              min = as.Date("2020/04/15"),
              max = as.Date("2021/11/20"),
              value = c(as.Date("2021/01/02"),as.Date("2021/11/20")))),
  dashboardBody(
    fluidRow(
      valueBoxOutput("vbox"),
      valueBoxOutput("vbox2"),
      valueBoxOutput("vbox3")
    ),
    fluidRow(splitLayout(cellWidths = c("50%", "50%"),
                         cellArgs = list(style = "horizontal-align: right"),
                         plotOutput(outputId = "distPlot", height = "500px",
                                    width= "500px"),
                         plotOutput(outputId = "regplot",height = "500px",width= "500px"))
  ))
)

