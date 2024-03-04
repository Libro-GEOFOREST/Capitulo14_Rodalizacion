# Capitulo14_Rodalizacion

La segmentación como medio para obtener una rodalización dasocrática uniforme de la estructura forestal es una técnica especialmente interesante cuando se combinan datos dasocráticos obtenidos a partir de la modelización LiDAR que se vio en el capítulo 13 del presente libro, junto con información espectral (por ejemplo, un índice NDVI de una imagen Landsat podría aportar información del estado fisiológico del arbolado) o mapas de orientaciones para que los rodales resultantes sigan las reglas de la rodalización tradicional basada en líneas permanentes del terreno.

En el presente ejemplo se va a partir de la modelización del área basimétrica realizada en el capítulo anterior. El usuario puede realizar si lo desea la modelización de otros datos dasocráticos. Aquí se facilitan la densidad (N) en pies/ha y la altura dominante (H) en metros.

```r
#Cargar la librería con la que se van a trabajar los datos raster
library(raster)
N<-raster("G:/Descarga/N.tif") #Adaptar a la ruta en la que se haya hecho la descarga
G<-raster("G:/Descarga/G.tif") #Adaptar a la ruta en la que se haya hecho la descarga
H<-raster("G:/Descarga/H.tif") #Adaptar a la ruta en la que se haya hecho la descarga

#Visualizar las capas raster
plot(N)
plot(G)
plot(H)
```


