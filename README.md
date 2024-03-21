# Capitulo14_Rodalizacion

La segmentación como medio para obtener una rodalización dasocrática uniforme de la estructura forestal es una técnica especialmente interesante cuando se combinan datos dasocráticos obtenidos a partir de la modelización LiDAR que se vio en el capítulo 13 del presente libro, junto con información espectral (por ejemplo, un índice NDVI de una imagen Landsat podría aportar información del estado fisiológico del arbolado) o mapas de orientaciones para que los rodales resultantes sigan las reglas de la rodalización tradicional basada en líneas permanentes del terreno.

En el presente ejemplo se va a partir de la modelización del área basimétrica realizada en el capítulo anterior. El usuario puede realizar si lo desea la modelización de otros datos dasocráticos. En la carpeta Capas_base se facilitan la densidad (N) en pies/ha y la altura dominante (H) en metros, así como el mapa de orientaciones reclasificado del 1 al 4 según el pixel pertenezca a Norte-Este-Sur-Oeste.

```r
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
```

Cuando se realiza un proceso de segmentación, se debe tener en cuenta que los algoritmos utilizarán los niveles digitales de los que están compuestos los píxeles de cada una de las capas raster. Ésto significa que tendrá mucho más peso la capa cuyo rango de valores sea mayor, en este caso la densidad, descrita en número de pies por hectárea y casi no quedará reflejada en la segmentación las orientaciones, puesto que sólo toma valores entre 1 y 4.

```r
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
```

Una vez reescalados, se unen en un mismo raster con tantas bandas como capas utilizadas y se guarda para su utilización durante el proceso de segmentación.
```r
#Función de reescalado
s <- stack(N_r, G_r, H_r,ori_r)
```
