
#Cargar la librería con la que se van a trabajar los datos raster
library(terra)
N<-rast("G:/Descarga/N.tif") #Adaptar a la ruta en la que se haya hecho la descarga
G<-rast("G:/Descarga/G.tif") #Adaptar a la ruta en la que se haya hecho la descarga
H<-rast("G:/Descarga/H.tif") #Adaptar a la ruta en la que se haya hecho la descarga
ori<-rast("G:/Descarga/orientaciones.tif") #Adaptar a la ruta en la que se haya hecho la descarga

#Visualizar las capas raster
plot(N)
plot(G)
plot(H)
plot(ori)

#Rango de valores de cada raster sin tener en cuenta los valores nulos (na.rm=TRUE)
range(values(N),na.rm=TRUE)
range(values(G),na.rm=TRUE)
range(values(H),na.rm=TRUE)
range(values(ori),na.rm=TRUE)

#Función de reescalado
reescalado = function(r){
  # obtener los valores máximos mínimos de la capa
  minmax_r = range(values(r), na.rm=TRUE) 
  # rescalar 
  return( (r-minmax_r[1]) / (diff(minmax_r)))
}

#Aplicación de la función de reescalado
N_r<-reescalado(N)
G_r<-reescalado(G)
H_r<-reescalado(H)
ori_r<-reescalado(ori)

#Crear raster multibanda
s <- stack(N_r, G_r, H_r,ori_r)

#Guardar raster
writeRaster(s,"G:/Descarga/s.tif") #Adaptar a la ruta que se desee emplear

#Comprobación de si está instalado el paquete necesario para la descarga desde GitHub
if(!requireNamespace("remotes")){
  install.packages("remotes")
}

# Instalación del paquete unbalanced necesario en algunos algoritmos de clasificación
remotes::install_github("dalpozz/unbalanced")

# Instalación de SegOptim con las últimas actualizaciones
remotes::install_github("joaofgoncalves/SegOptim")

#Directorio de la carpeta bin del programa
otb_path <- "C:/OTB-8.0.1-Win64/bin" #Adaptar a la ruta en la que se haya instalad OTB

#Directorio de trabajo
setwd("G:/Descarga/")  #Adaptar a la ruta en la que se vayan a configurar las salidas de la segmentación

#Ruta de la imagen a segmentar
imagen.path <- "G:/Descarga/s.tif" #Adaptar a la ruta en la que se haya guardado la imagen generada en el paso 1

#Visualización del histograma de las imágenes reescaladas.
par(mfrow=c(2,2))
hist(s,layer=1,breaks=50)
hist(s,layer=2,breaks=50)
hist(s,layer=3,breaks=50)
hist(s,layer=4,breaks=50)

#Tamaño de pixel de la imagen s
res(s)

#Activación de la librería SegOptim
library(SegOptim)

#Aplicación de la segmentación Mean-Shift a través de Orfeo
out_segm_obj <- segmentation_OTB_LSMS(
  inputRstPath  = imagen.path,               #Ruta en la que se encuentra la imagen a segmentar
  SpectralRange = 0.05,                      #Rango espectral
  SpatialRange  = 5,                         #Rango espacial
  MinSize       = 30,                        #Tamaño mínimo del segmento
  lsms_maxiter  = 50,                        #Número máximo de iteraciones
  outputSegmRst = "./segmRaster.tif",        #Imagen raster resultante
  verbose       = TRUE,                      #Muestra en pantalla la evolución del algoritmo
  otbBinPath    = otb_path)                  #Ruta en la que se encuentra el programa Orfeo Toolbox

# Cargar el raster de segmentos y visualizarlo
library(terra)
segm_rst <- rast(out_segm_obj$segm)
plot(segm_rst)

#Transformación del raster a shapefile
segm_sf<- as.polygons(segm_rst)

#Visualización del shapefile
library(mapview)
mapview(segm_sf,map.types="Esri.WorldImagery",color = 'cyan',legend =FALSE,
        alpha.regions = 0)

#Cálculo de la superficie de cada segmento
segm_sf$area<-expanse(segm_sf,unit="m")

#Tamaño medio de los segmentos en hectáreas
mean(segm_sf$area)/10000

#Tamaño máximo de los segmentos en hectáreas
max(segm_sf$area)/10000

#Tamaño minimo de los segmentos en hectáreas
min(segm_sf$area)/10000

#Número de segmentos que componen la segmentación
length(segm_sf$area)

#Aplicación de la segmentación Mean-Shift a través de Orfeo
out_segm_obj <- segmentation_OTB_LSMS(
  inputRstPath  = imagen.path,               #Ruta en la que se encuentra la imagen a segmentar
  SpectralRange = 0.1,                       #Rango espectral
  SpatialRange  = 10,                        #Rango espacial
  MinSize       = 30,                        #Tamaño mínimo del segmento
  lsms_maxiter  = 50,                        #Número máximo de iteraciones
  outputSegmRst = "./segmRaster.tif",        #Imagen raster resultante
  verbose       = TRUE,                      #Muestra en pantalla la evolución del algoritmo
  otbBinPath    = otb_path)                  #Ruta en la que se encuentra el programa Orfeo Toolbox

# Cargar el raster de segmentos y visualizarlo
segm_rst2 <- rast(out_segm_obj$segm)
plot(segm_rst2)

#Transformación del raster a shapefile
segm_sf2<- as.polygons(segm_rst2)

#Visualización del shapefile
mapview(segm_sf2,map.types="Esri.WorldImagery",color = 'red',legend =FALSE,
        alpha.regions = 0)

#Para la primera segmentacion
medias<-zonal(s,segm_rst,fun="mean")      #Calculo de los valores medios en los segmentos para las imágenes de densidad, altura dominante, área basimétrica y orientaciones reescaladas
desviaciones<-zonal(s,segm_rst,fun="sd")  #Calculo de las desviaciones estándar en los segmentos para las imágenes de densidad, altura dominante, área basimétrica y orientaciones reescaladas

#Varianza
varianzas<-(desviaciones^2)
varianzas<-as.data.frame(varianzas)

#Homogeneidad interna
homo.int.N<-sum(segm_sf$area*varianzas$N)/sum(segm_sf$area)
homo.int.G<-sum(segm_sf$area*varianzas$G)/sum(segm_sf$area)
homo.int.H<-sum(segm_sf$area*varianzas$H,na.rm=TRUE)/sum(segm_sf$area)
homo.int.ori<-sum(segm_sf$area*varianzas$orientaciones,na.rm=TRUE)/sum(segm_sf$area)

homo.int.N;homo.int.G;homo.int.H;homo.int.ori

#Heterogeniedad externa
library(spdep)
library(sf)
#Definir los poligonos vecinos
nb <- poly2nb(st_as_sf(segm_sf), queen=TRUE)
nb[1] #Imprimirá en pantalla los polígonos vecinos al primero

#Asignar pesos a los vecinos
lw <- nb2listw(nb, style="W", zero.policy=TRUE)
lw$weights[1] #Imprimirá en pantalla los pesos asignados a cada vecino del primer polígono

#Cálculo del estadístico del índice de Moran
hete.ext.N<- moran(medias$N, lw, length(nb), Szero(lw),NAOK = TRUE)[1]
hete.ext.G<- moran(medias$G, lw, length(nb), Szero(lw),NAOK = TRUE)[1]
hete.ext.H<- moran(medias$H, lw, length(nb), Szero(lw),NAOK = TRUE)[1]
hete.ext.ori<- moran(medias$orientaciones, lw, length(nb), Szero(lw),NAOK = TRUE)[1]

hete.ext.N;hete.ext.G;hete.ext.H;hete.ext.ori

#Para la segunda segmentacion
#Cálculo de superficies
segm_sf2$area<-expanse(segm_sf2,unit="m")
medias2<-zonal(s,segm_rst2,fun="mean")      #Calculo de los valores medios en los segmentos para las imágenes de densidad, altura dominante, área basimétrica y orientaciones reescaladas
desviaciones2<-zonal(s,segm_rst2,fun="sd")  #Calculo de las desviaciones estándar en los segmentos para las imágenes de densidad, altura dominante, área basimétrica y orientaciones reescaladas

#Varianza
varianzas2<-(desviaciones2^2)
varianzas2<-as.data.frame(varianzas2)

#Homogeneidad interna
homo.int.N2<-sum(segm_sf2$area*varianzas2$N)/sum(segm_sf2$area)
homo.int.G2<-sum(segm_sf2$area*varianzas2$G)/sum(segm_sf2$area)
homo.int.H2<-sum(segm_sf2$area*varianzas2$H,na.rm=TRUE)/sum(segm_sf2$area)
homo.int.ori2<-sum(segm_sf2$area*varianzas2$orientaciones,na.rm=TRUE)/sum(segm_sf2$area)

homo.int.N2;homo.int.G2;homo.int.H2;homo.int.ori2

#Heterogeniedad externa
library(spdep)
library(sf)
#Definir los poligonos vecinos
nb2 <- poly2nb(st_as_sf(segm_sf2), queen=TRUE)
nb2[1]

#Asignar pesos a los vecinos
lw2 <- nb2listw(nb2, style="W", zero.policy=TRUE)
lw2$weights[1]

#Cálculo del estadístico del índice de Moran
hete.ext.N2<- moran(medias2$N, lw2, length(nb2), Szero(lw2),NAOK = TRUE)[1]
hete.ext.G2<- moran(medias2$G, lw2, length(nb2), Szero(lw2),NAOK = TRUE)[1]
hete.ext.H2<- moran(medias2$H, lw2, length(nb2), Szero(lw2),NAOK = TRUE)[1]
hete.ext.ori2<- moran(medias2$orientaciones, lw2, length(nb2), Szero(lw2),NAOK = TRUE)[1]

hete.ext.N2;hete.ext.G2;hete.ext.H2;hete.ext.ori2

#Valor medio de homogeneidad interna de la primera segmentación
mean(homo.int.N[[1]],homo.int.N[[1]],homo.int.N[[1]],homo.int.N[[1]])

#Valor medio de homogeneidad interna de la segunda segmentación
mean(homo.int.N2[[1]],homo.int.N2[[1]],homo.int.N2[[1]],homo.int.N2[[1]])

#Valor medio de heterogeneidad externa de la primera segmentación
mean(hete.ext.N[[1]],hete.ext.G[[1]],hete.ext.H[[1]],hete.ext.ori[[1]])

#Valor medio de heterogeneidad externa de la segunda segmentación
mean(hete.ext.N2[[1]],hete.ext.G2[[1]],hete.ext.H2[[1]],hete.ext.ori2[[1]])