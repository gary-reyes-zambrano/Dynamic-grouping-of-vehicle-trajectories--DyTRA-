#rm(list = ls())
#options( digits = 15 )
library(xlsx)

folder_result <- ruta
setwd(folder_result)

ipak(packages)

atipicas <- list()
vela <- list()

ev <- 40

length(CopyDyclee)
tdc <- CopyDyclee[[ev]]
length(tdc$oList)
tt <- data.frame()
for(i in 1:length(tdc$oList)){
  print(length(tdc$oList[[i]]$data$vel_prom))
  print(tdc$oList[[i]]$data$vel_prom)
}
#tt <- rbind(tt, tdc$oList[[3]]$data)



velocidad_ponderada <- 43.97

num_new_g <- 1
atipicas[[num_new_g]] <- data.frame()
vela[[num_new_g]] <- velocidad_ponderada
for( x in c(1)){
  numero_atipica <- x
  #print(round(mean(dyclee$oList[[x]]$tramo$velocidad,2)))
  data <- tdc$oList[[x]]$data
  atipicas[[num_new_g]] <- rbind(atipicas[[num_new_g]], data)
}






packages_map <- c("leaflet","mapview","webshot","htmltools","ggrepel")
ipak(packages)
ipak(packages_map)

setwd(folder_result)
print("MapaRangos")

setwd(folder_result)
dir.create(paste("MapaRangosCD",sep=''))
dir_map <- getwd()
setwd(paste(dir_map,'//',"MapaRangosCD",sep=''))
print("Mapas Interactivos por Rangos de Velocidades")

min_v <- 0
max_v <- 200
#colores <- colorNumeric("viridis", domain = c(min_v,max_v))
#colores <- colorNumeric(brewer.pal(4,"RdYlGn"), domain = c(min_v,max_v),reverse = TRUE) #For normal color view
colores <- colorNumeric(brewer.pal(4,"RdYlBu"), domain = c(min_v,max_v),reverse = TRUE) #For colorblind safe
#https://colorbrewer2.org/#type=diverging&scheme=RdYlBu&n=3


#----------------------------------------
def_cants <- NULL
celda_vel <- TRUE  #True: Mostrar el color de la velocidade n la celda. False: Color diferente para cada grupo
#----------------------------------------
list_cant <- list()

get_map <- function(){
  mapa <- leaflet() %>% addTiles(group = "Open Street Map")
  #mapa <- addTiles(mapa, urlTemplate = "http://mt0.google.com/vt/lyrs=m&hl=en&x={x}&y={y}&z={z}&s=Ga", attribution = 'Google Maps', group = "Google Maps")
  #mapa <- addTiles(mapa, urlTemplate = "https://mts1.google.com/vt/lyrs=s&hl=en&src=app&x={x}&y={y}&z={z}&s=G", attribution = 'Google Satellite Maps', group = "Google Satellite Maps")
  #mapa <- mapa %>% addProviderTiles("Esri.WorldImagery", group = "ESRI") %>% addProviderTiles("Stamen.TonerLines", group = "Stamen")
  #leaflet(data = data) %>% addMarkers() %>% addTiles(group = "OSM(default)") %>% addProviderTiles("Esri.WorldImagery", group = "ESRI") %>% addProviderTiles("Stamen.Toner", group = "Stamen") %>% addLayersControl(baseGroup = c("OSM(default", "ESRI", "Stamen"))
  return(mapa)
}

get_centro <- function(cld){
  dfc <- data.frame(lon = (cld$lon_minima + cld$lon_maxima) / 2 ,
                    lat = (cld$lat_minima + cld$lat_maxima) / 2 )
  return(dfc)
}

get_points_prom_group <- function(grupos){
  cp <- data.frame()
  for(b in grupos){
    cp <- rbind(cp,data.frame(grupo = b$label, cp = round(mean(b$data[,"vel_prom"]),2)))
  }
  #print(cp)
  return(round(mean(cp[,"cp"]),0))
}

#' @param CopyDyClee
#' @param CopyBinnacle

dataContext <-  list(c(maximun=max(dataset$data$longitude),
                       minimun=min(dataset$data$longitude)),
                     c(maximun=max(dataset$data$latitude),
                       minimun=min(dataset$data$latitude)))

#---------------------------------------------
print(paste("Mapa Evolucion ",ev, sep=""))
dc <- CopyDyclee[[ev]]
points_prom_group <- 0
mapa <- get_map()
layerBaseGroup <- c("Open Street Map")
capas <- c()

capas_grupos <- c("<b>Todos los rangos<b/>")
grupos <- dc$aList

if(length(grupos)>0){
  col_rb <- rainbow(length(grupos))
  
  for(gl in 1:length(grupos)){
    un_grupo <- grupos[[gl]]
    s_prom <- round(mean(un_grupo$data[,"cant_segmentos"]),2)
    c_prom <- round(mean(un_grupo$data[,"cant_veh"]),2)
    v_prom <- round(mean(un_grupo$data[,"vel_prom"]),2)
    c_puntos <- sum(un_grupo$data[,"cant_segmentos"])
    c_veh <- sum(un_grupo$data[,"cant_veh"])
    v_min_grupo <- min(un_grupo$data[,"vel_prom"])
    v_max_grupo <- max(un_grupo$data[,"vel_prom"])
    
    nombre_capa <- paste(v_prom," Km/h ",sep="")
    ncv <- c(nombre_capa, "<b>Todos los rangos<b/>")
    capas_grupos <- c(capas_grupos, nombre_capa)
    cg <- un_grupo$data
    
    if(celda_vel){
      tcolor <- colores(v_prom)
    }else{
      tcolor <- col_rb[gl]
    }
    
    for(n in 1:nrow(cg)){
      mapa <- addRectangles(map = mapa,
                            lng1 = cg[n,"lon_maxima"], lat1 = cg[n,"lat_maxima"],
                            lng2 = cg[n,"lon_minima"], lat2 = cg[n,"lat_minima"],
                            fillColor = tcolor, color = tcolor, weight = 2,
                            fillOpacity = 0.8, 
                            group = ncv)
    }
    
    content <- paste(paste("<b>Info Grupo</b>"),
                     paste("ID Grupo:",gl),
                     paste("Vel. Min.:",v_min_grupo,"Km/h"),
                     paste("Vel. Max:",v_max_grupo,"Km/h"),
                     paste("Vel.Prom.",v_prom,"Km/h"),
                     paste("Puntos totales:",c_puntos),
                     paste("Puntos promedio:",c_prom),
                     paste("Vehiculos totales:",c_veh),
                     paste("Vehiculos promedio:",c_prom),
                     sep="<br/>")
    centros <- get_centro(cg[nrow(cg),])
    mapa <- addMarkers(map = mapa,
                       lng = centros$lon,
                       lat = centros$lat,
                       popup = content, group = ncv)
    
  }
  
}

if(length(atipicas)>0){
  ol <- atipicas
  
  col_rb <- rainbow(length(ol))
  
  for(gl in 1:length(ol)){
    un_grupo <- ol[[gl]]
    s_prom <- round(mean(un_grupo[,"cant_segmentos"]),2)
    c_prom <- round(mean(un_grupo[,"cant_veh"]),2)
    v_prom <- vela[[gl]]#round(mean(un_grupo[,"vel_prom"]),2)
    c_puntos <- sum(un_grupo[,"cant_segmentos"])
    c_veh <- sum(un_grupo[,"cant_veh"])
    v_min_grupo <- min(un_grupo[,"vel_prom"])
    v_max_grupo <- max(un_grupo[,"vel_prom"])
    
    nombre_capa <- paste(v_prom," Km/h ",sep="")
    ncv <- c(nombre_capa, "<b>Todos los rangos<b/>")
    capas_grupos <- c(capas_grupos, nombre_capa)
    cg <- un_grupo
    
    if(celda_vel){
      tcolor <- colores(v_prom)
    }else{
      tcolor <- col_rb[gl]
    }
    
    for(n in 1:nrow(cg)){
      mapa <- addRectangles(map = mapa,
                            lng1 = cg[n,"lon_maxima"], lat1 = cg[n,"lat_maxima"],
                            lng2 = cg[n,"lon_minima"], lat2 = cg[n,"lat_minima"],
                            fillColor = tcolor, color = tcolor, weight = 2,
                            fillOpacity = 0.8, 
                            group = ncv)
    }
    
    content <- paste(paste("<b>Info Grupo</b>"),
                     paste("ID Grupo:",gl),
                     paste("Vel. Min.:",v_min_grupo,"Km/h"),
                     paste("Vel. Max:",v_max_grupo,"Km/h"),
                     paste("Vel.Prom.",v_prom,"Km/h"),
                     paste("Puntos totales:",c_puntos),
                     paste("Puntos promedio:",c_prom),
                     paste("Vehiculos totales:",c_veh),
                     paste("Vehiculos promedio:",c_prom),
                     sep="<br/>")
    centros <- get_centro(cg[nrow(cg),])
    mapa <- addMarkers(map = mapa,
                       lng = centros$lon,
                       lat = centros$lat,
                       popup = content, group = ncv)
    
  }
  
}

scaleBarOptions(
  maxWidth = 100,
  metric = TRUE,
  imperial = TRUE,
  updateWhenIdle = TRUE
)

capas <- c(capas, capas_grupos)
mapa <- addLegend(map = mapa,pal = colores,values = c(min_v,max_v), opacity = 1, title = "Velocidades (Km/h)",position = "bottomleft")
mapa <- addScaleBar(mapa, position = "bottomleft", options = scaleBarOptions()) #bottomright
mapa <- hideGroup(mapa, c(capas))
mapa <- addLayersControl(map = mapa, baseGroups = layerBaseGroup,overlayGroups = capas,options = layersControlOptions(collapsed = FALSE)) #options = layersControlOptions(collapsed = FALSE)

mapshot(mapa, url = paste("Evolucion_",ev,".html",sep=""))


#setwd(folder_result)
save.image(file="RDataMapsRangoPonderado.RData")
#=====================================================================================


