require(ggplot2)
require(gridExtra)
require("RColorBrewer")
require("shinydashboard")
server <- function(input, output) {
dataRE<-function(id_region,f1)
{
	ftregion<-read.csv("ftregion.txt",sep=";",header=T)
	regiones=c("Nacional","Arica","Tarapaca","Antofagasta","Atacama","Coquimbo","Valparaiso","R.Metropolitana","O'Higgins","Maule","Nuble","Bio-Bio","Araucania","Los Rios","Los Lagos","Aysen","Magallanes","Norte Grande","Norte Chico","Zona Central","Zona Sur","Zona Austral","Regiones")
	idregion=c(0,15,1:5,13,6:7,16,8:9,14,10:12,50:55)
	idr=idregion[which(regiones==id_region)]
   res<-ftregion[which(ftregion$id_region==idr),]
   res<- res[which(as.Date(res$fch_confirmado) == f1[2]),]
	res<- res$cant_uci  
	return(res)
}
plotUCI<-function(id_region,f1)
{
	ftregion<-read.csv("ftregion.txt",sep=";",header=T)
	regiones=c("Nacional","Arica","Tarapaca","Antofagasta","Atacama","Coquimbo","Valparaiso","R.Metropolitana","O'Higgins","Maule","Nuble","Bio-Bio","Araucania","Los Rios","Los Lagos","Aysen","Magallanes","Norte Grande","Norte Chico","Zona Central","Zona Sur","Zona Austral","Regiones")
	idregion=c(0,15,1:5,13,6:7,16,8:9,14,10:12,50:55)
	idr=idregion[which(regiones==id_region)]
   res<-ftregion[which(ftregion$id_region==idr),]
   res=res[which(as.Date(res$fch_confirmado)<=f1[2] & as.Date(res$fch_confirmado)>=f1[1]),]
   fecha=seq(f1[1],f1[2],by=1)
   ftregion2=data.frame(res,"fecha"=fecha)
   p <- ggplot(ftregion2, aes(x=fecha, group = 1)) 
   p <- p + geom_line(aes(x=fecha, y=cant_uci),stat="identity", colour='steelblue')
   p <- p + geom_point(aes(x=fecha, y=cant_uci),color='dodgerblue4', fill="#69b3a2", size=5, alpha=0.3)
   p <- p + geom_text(aes(label= cant_uci, x=fecha, y=cant_uci), colour="dodgerblue4",position=position_dodge(width = 1), vjust = -1.5, size = 3.2)
   p <- p + scale_x_date(date_breaks="5 day")
   p <- p + xlab('Fecha') + ylab('Numero de Pacientes Covid-19 en UCI')
   p <- p + theme_bw()
   #p <- p + theme(plot.background = element_rect(fill = "gray"))
   p <- p + theme(
   axis.text = element_text(colour="dodgerblue4"), axis.text.x = element_text(angle = 90)
)
grid.arrange(p,ncol=1)
}
porcocu<-function(id_region,f1)
{
uci=read.csv("Camas_UCI_diarias.txt",sep=",",header=T)
#uci=na.omit(uci1)
#uci=uci1[,1:(length(uci1)-1)]
#uci=uci1[,1:325]
#uci=uci1[326]
#View(uci1)
#write.table(uci,"Camas_UCI_diarias22.txt",row.names=T,sep=",")
fecha=seq(as.Date("2020/04/14"),as.Date("2021/12/31"),by=1)
fecha=fecha[2:(dim(uci)[2]-2)]
pos2=which(fecha==f1[2])+3
if(length(pos2)==0){
pos2=dim(uci)[2]
}
regiones=c("Nacional","Arica","Tarapaca","Antofagasta","Atacama","Coquimbo","Valparaiso","R.Metropolitana","O'Higgins","Maule","Nuble","Bio-Bio","Araucania","Los Rios","Los Lagos","Aysen","Magallanes","Norte Grande","Norte Chico","Zona Central","Zona Sur","Zona Austral","Regiones")
	iduci=c(16,3,15,1,4,7,17,12,14,11,13,6,2,9,8,5,10)
idr=iduci[which(regiones==id_region)]
if(id_region=="Norte Grande")
idr=iduci[2:4]
if(id_region=="Norte Chico")
idr=iducii[5:6]
if(id_region=="Zona Central")
idr=iduci[7:12]
if(id_region=="Zona Sur")
idr=iduci[13:15]
if(id_region=="Zona Austral")
idr=iduci[16:17]
if(id_region=="Regiones")
idr=iduci[-c(1,8)]
nombres=names(table(uci$Region))
aux=uci[which(uci[,1] %in% nombres[idr]),]
aux2=aux[,pos2]
k=length(aux2)/4
porc=sum(aux2[(3*k+1):(4*k)])/sum(aux2[(0*k+1):(1*k)])
porc2=sum(aux2[(k+1):(2*k)])/sum(aux2[(3*k+1):(4*k)])
return(c(paste(round(porc*100,2),"%"),paste(round(porc2*100,2),"%")))
}

UCIReg<-function(id_region,f1)
{
   uci=read.csv("Camas_UCI_diarias.txt",sep=",",header=T)
   #as.data.frame(uci1)
   #uci=uci1[,1:(length(uci1)-1)]
   #write.table(uci,"Camas_UCI_diarias22.txt",row.names=T,sep=",")
   #uci=na.omit(uci1)
   #uci=uci1[,1:325]
   fecha=seq(as.Date("2020/04/14"),as.Date("2021/12/31"),by=1)
   fecha=fecha[2:(dim(uci)[2]-2)]
   pos=which(fecha==f1[2])+3

if(length(pos)==0){
pos=dim(uci)[2]
}
nombres=names(table(uci$Region))
porc=rep(0,length(nombres))
for(i in 1:length(nombres))
{
	aux=uci[which(uci[,1]==nombres[i]),]
	aux2=aux[,pos]
	porc[i]=sum(aux2[4])/aux2[1]
}
pos2=c(16,15,1,4,7,17,14,11,6,2,8,5,10,12,9,3,13)

nombres2=nombres[pos2]
porc2=porc[pos2]
regiones=c("Nacional","Tarapaca","Antofagasta","Atacama","Coquimbo","Valparaiso","Ohiggins","Maule","Bio-Bio","Araucania","Los Lagos","Aysen","Magallanes","Metropolitana","Los Rios","Arica","Nuble")
df=data.frame("Region"=regiones,"id_region"=0:16,"Ocupacion"=porc2)
nac=as.numeric(porc2[1])
df2=df[-1,]
df=df2[c(15,1:5,13,6:7,16,8:9,14,10:12),]
barplot(df[,3],main="% Ocupacion Camas UCI por Region",xlab="",ylab="% Ocupacion UCI",ylim=c(0,1), names.arg=df[,1],cex.names=0.7,las=2,col=brewer.pal(9, "Blues"))
abline(h=nac,col="red",lwd=2)
text(15,0.95,paste("Ocupacion Nacional",round(nac*100,2),"%"),col="red",lwd=2)
}
	output$vbox <- renderValueBox({
    valueBox(
      "Pacientes Covid-19 en UCI",
      value=tags$p(dataRE(input$Region,input$date), style = "font-size: 200%;"),
      icon = icon("arrow-up")
    )}) 
    output$vbox2 <- renderValueBox({
    valueBox(
      "% Ocupacion Camas UCI",
      value=tags$p(porcocu(input$Region,input$date)[1], style = "font-size: 200%;"),color="red"
    )})  
    output$vbox3 <- renderValueBox({
    valueBox(
      "% Camas UCI ocupadas por Covid-19",
      value=tags$p(porcocu(input$Region,input$date)[2], style = "font-size: 200%;")
    )})  
    output$distPlot <- renderPlot({
  	plotUCI(input$Region,input$date)
    }, height = 450, width = 600 )
    output$regplot <- renderPlot({
  	UCIReg(input$Region,input$date)
    }, height = 450, width = 600 )
    output$downplot <- downloadHandler(
    filename = function(){
      paste("UCI", input$Region,".png", sep = "")},
    content = function(file){
        png(file, height = 500, width = 800 )
      plotUCI(input$Region,input$date)
      dev.off()}
  )
}

