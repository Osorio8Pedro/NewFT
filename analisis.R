names(ftregion)
ftregion<-read.csv("ftregion.txt",sep=";",header=T)
ftregion<-read.csv("ftcomuna.txt",sep=";",header=T)
ftis<-read.csv("ftiniciosintomas.txt",sep=";",header=T)
summary(ftis)
table(ftis$fch_publicacion)

ftis2=ftis[which(ftis$fch_publicacion=="2020-10-26"),]
ftis3=ftregion[which(ftregion$id_region==4),]
ftis3
names(ftis3)
View(ftis3)
table(ftis3)
ftis3=ftis2[which(ftis2$id_region==2),]
ftis3
length(ftis3)
t1=aggregate(cant_casosconfirmados~id_semanaepidemiologica,ftis3,sum)
t1
t2=aggregate(cant_casosconfirmadoschowellw5~id_semanaepidemiologica,ftis3,sum)
t2
plot(t1[,1],t1[,2],type="l")
lines(t2[,1],t2[,2],col="red")


############################
########################
##Análisis Descriptivo


###Cuerpo
 install.packages("remotes")
remotes::install_github("jespermaag/gganatogram")
library(gganatogram)
 install.packages("dplyr")
library(dplyr)

hgMale_key %>%
   filter(organ %in% c( "lung")) %>%
   gganatogram(organism = "human", sex = "male",
               fill = "colour") +
   theme_void() + 
   coord_fixed()

hgMale_key %>%
   filter(type %in% "nervous_system") %>%
   gganatogram(organism = "human", sex = "male",
               fill = "colour", outline = FALSE) +
   theme_void() + 
   coord_fixed()

##calendario

install.packages("calendR")
library(calendR)

# Datos
datos <- rnorm(30, 15, 10)

# Crea un vector donde todos los valores son ligeramente
# inferiores que el menor valor de tus datos
dias <- rep(min(datos) - 0.05, 365)

# Rellena los días que quieras con tus datos
dias[30:59] <- datos

calendR(year = 2021,
        special.days = dias,
        low.col = "white",
        special.col = "#FF0000",
        gradient = TRUE,
        legend.pos = "bottom")
###########
night_owlish <- "https://raw.githubusercontent.com/batpigandme/night-owlish/master/rstheme/night-owlish.rstheme"
rstudioapi::addTheme(night_owlish, apply = TRUE)
##########
set.seed(2)
ftregion<-read.csv("newftregion.txt",sep=";",header=T)
names(ftregion)
data_set <- data.frame(Region = ftregion$txt_nombreregion,
   poblacion = ftregion$cant_poblacion,
                       type = sample(1:4, size = 25, replace = TRUE),
                       store = sample(paste("Store", 1:4),
                                      size = 25, replace = TRUE))

head(data_set)

#############