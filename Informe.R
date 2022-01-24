install.packages("remotes")
remotes::install_github("reconhub/trendbreaker")
install.packages("trendbreaker")
require(ggplot2)
require(gridExtra)
require(RColorBrewer)
library(incidence2) 
library(trendbreaker) # for ASMODEE
library(dplyr)        # for data manipulation
library(future)
ftregion<-read.csv("ftregion.txt",sep=";",header=T)
fechas=names(table(ftregion$fch_confirmado))
l1=length(fechas)#218

names(ftregion)
ftregion8=ftregion %>% select(fch_confirmado,txt_nombreregion,cant_casosconfirmadosdiario,cant_poblacion,porc_positividad,cant_uci,cant_fallecidos)

ftregion8$cant_casosconfirmadosdiario<-((ftregion8$cant_casosconfirmadosdiario)/(ftregion8$cant_poblacion)*100000)
#Incidencia
d=as.numeric(ftregion8$cant_casosconfirmadosdiario)
ftregion8$cant_casosconfirmadosdiario=round(d,2)
#Positividad
positivo=(ftregion8$porc_positividad)*100
ftregion8$porc_positividad=round(positivo,2)
#Mortalidad

ftregion8$fch_confirmado<-as.Date(ftregion8$fch_confirmado)
fecha=as.Date("2020-03-03")
covid19_nhs=ftregion8 %>%
   mutate(forma = fch_confirmado-fecha)%>%
   mutate(forma2 = as.character(fch_confirmado, format="%A"))



covid19_nhs$forma=as.numeric(covid19_nhs$forma)
nombres=c("fecha","nhs_region","Incidencia","cant_pobla","Positividad","UCI","Mortalidad","day","weekday")
colnames(covid19_nhs)<-(nombres)
covid19_nhs$Incidencia=as.numeric(covid19_nhs$Incidencia)
covid19_nhs$Positividad=as.numeric(covid19_nhs$Positividad)
covid19_nhs$Mortalidad=as.numeric(covid19_nhs$Mortalidad)
covid19_nhs$UCI=as.numeric(covid19_nhs$UCI)


nombres=c("Tarapaca","Antofagasta","Atacama","Coquimbo","Valparaiso","Ohiggins","Maule","Bio-Bio","Araucania","Los Lagos","Aysen","Magallanes","Metropolitana","Los Rios","Arica","Nuble")

covid19_nhs=covid19_nhs[2:17,]
covid19_nhs[,1]=nombres
df=df[c(15,1:5,13,6:7,16,8:9,14,10:12),]

nhs_pathways_covid19<-covid19_nhs
t=nhs_pathways_covid19%>%
   arrange(nhs_region)

nhs_pathways_covid19$fecha=as.Date(nhs_pathways_covid19$fecha)

nhs_pathways_covid19<-filter(nhs_pathways_covid19, #nhs_region=='00. Nacional') #nhs_region=='17. Norte Grande' |nhs_region=='18. Norte Chico' | nhs_region=='19. Zona Central' |nhs_region=='20. Zona Sur'| nhs_region=='21. Zona Austral'|nhs_region=='22. Regions Sin Metropolitana')
                             
                             #                         nhs_region=='16. Región de Magallanes y Antártica Chilena')
                             #nhs_pathways_covid19<-filter(nhs_pathways_covid19, id_region=='1'|id_region=='2'|id_region=='3'|id_region=='4'|id_region=='5'|id_region=='6'|id_region=='7'|id_region=='8'|id_region=='9'|id_region=='10'|id_region=='11'|id_region=='12'|id_region=='13'|id_region=='14'|id_region=='15'|id_region=='16')
                             
                             nhs_region=='01. Región de Arica y Parinacota'| nhs_region=='02. Región de Tarapacá'| nhs_region=='03. Región de Antofagasta' | nhs_region=='04. Región de Atacama' | nhs_region=='05. Región de Coquimbo' | nhs_region=='06. Región de Valparaíso' | nhs_region=='07. Región Metropolitana de Santiago'  | nhs_region=="08. Región del Libertador Bernardo O'Higgins"|nhs_region=='09. Región del Maule' |
                                nhs_region=='10. Región del Ñuble'| nhs_region=='11. Región del Bío-Bío'| nhs_region=='12. Región de La Araucanía' | nhs_region=='13. Región de Los Ríos' | nhs_region=='14. Región de Los Lagos' | nhs_region=='15. Región de Aysén del General Ibañez del Campo' | nhs_region=='16. Región de Magallanes y Antártica Chilena'  )


first_date <-fechas[l1-42]
last_date<-fechas[l1]

pathways_recent <- nhs_pathways_covid19 %>%
   filter(fecha >= first_date)%>%
   filter(fecha <= last_date)

lookup <- select(pathways_recent, fecha, day, weekday) %>%  distinct()
nhs_pathways_covid19$fecha<-as.Date(nhs_pathways_covid19$fecha)

#Incidencia
dat <-
   pathways_recent %>%
   incidence(date_index = fecha, groups = nhs_region, count = Incidencia) %>%
   left_join(lookup, by = c("date_index" = "fecha"))


models <- list(
   regression = lm_model(Incidencia ~ day),
   poisson_constant = glm_model(Incidencia ~ 1, family = "poisson"),
   negbin_time = glm_nb_model(Incidencia ~ day),
   negbin_time_weekday = glm_nb_model(Incidencia ~ day + weekday)
)
res <- asmodee(dat, models, method = evaluate_aic, alpha=0.05, k = 7)
plot(res)
#p <- p +geom_hline(yintercept=10,color="red")

#Positividad
dat <-
   pathways_recent %>%
   incidence(date_index = fecha, groups = nhs_region, count = Positividad) %>%
   left_join(lookup, by = c("date_index" = "fecha"))


models <- list(
   regression = lm_model(Positividad ~ day),
   poisson_constant = glm_model(Positividad ~ 1, family = "poisson"),
   negbin_time = glm_nb_model(Positividad ~ day),
   negbin_time_weekday = glm_nb_model(Positividad ~ day + weekday)
)
res <- asmodee(dat, models, method = evaluate_aic, alpha=0.05, k = 7)
plot(res)

#UCI
dat <-
   pathways_recent %>%
   incidence(date_index = fecha, groups = nhs_region, count = UCI) %>%
   left_join(lookup, by = c("date_index" = "fecha"))


models <- list(
   regression = lm_model( UCI ~ day),
   poisson_constant = glm_model( UCI ~ 1, family = "poisson"),
   negbin_time = glm_nb_model( UCI ~ day),
   negbin_time_weekday = glm_nb_model( UCI ~ day + weekday)
)
res <- asmodee(dat, models, method = evaluate_aic, alpha=0.05, k = 7)
plot(res)

#Mortalidad
dat <-
   pathways_recent %>%
   incidence(date_index = fecha, groups = nhs_region, count = Mortalidad) %>%
   left_join(lookup, by = c("date_index" = "fecha"))


models <- list(
   regression = lm_model(Mortalidad ~ day),
   poisson_constant = glm_model(Mortalidad ~ 1, family = "poisson"),
   negbin_time = glm_nb_model(Mortalidad ~ day),
   negbin_time_weekday = glm_nb_model(Mortalidad ~ day + weekday)
)
res <- asmodee(dat, models, method = evaluate_aic, alpha=0.05, k = 7)
plot(res)

