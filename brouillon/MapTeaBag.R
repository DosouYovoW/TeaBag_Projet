library(tidyverse)
library(leaflet)
library(sp)
data =data.frame(lon=c(-71.744876, -72037064, -71.870602, -71.858327),
                 lat=c(46.324630, 46.237920, 46.277615, 46.276309),
                 Name =c('Site10', 'SiteA9', 'Site9', 'Site45'))

str(data)
map = leaflet() %>%
  addTiles() %>%
  addMarkers(lat=data$lat, lng=data$long)

#coord <- as.factor(coord$Name)
data <- data[complete.cases(data),]
data$lon <- as.numeric(data$lon)
data$lat <- as.numeric(data$lat)
data <- as.data.frame(data)

data.sp <- SpatialPointsDataFrame(data[,c(2,1)], data[,-c(2,1)])
map
 
 addProviderTiles("Stamen.TonerHybrid")
  
  for (i in 1:nrow(coord)) {
  #map = map %>% addMarkers(lng=coord$lon[i], lat=coord$lat[i])
  map = map %>% addCircles(lng=coord$lon[i], lat=coord$lat[i], color="blue") # couleur code hexadécimal
}
map

----------------------------------------------------
require(maps)
require(mapdata) 
library(ggplot2)                        
library(ggrepel)
library(ggmap)
Coord <- data.frame(longitude=c(-71.744876, -72.037064, -71.870602, -71.858327),
                    latitude=c(46.324630, 46.237920, 46.277615, 46.276309),
                    stringsAsFactors = FALSE,
                    Name =c('Site10', 'SiteA9', 'Site9', 'Site45'))
#points <- data.frame(lon = c(-71, -72),
                     #lat = c(46, 46.5))

#map_qc <- get_googlemap(center = c(lon = -71, lat = 46),
                        #zoom = 9, scale = 2, maptype ='terrain')
boite <- c(left = -72.1, bottom = 46.2, right = -71.7, top = 46.4)
map <- get_stamenmap(bbox = boite, zoom = 10, maptype = "toner-lite")
ggmap(map)

#boite <- c(left = -70, bottom = 42, right = -72, top = 47)
#map <- get_stamenmap(bbox = boite, zoom = 10, maptype = "toner-lite")
#ggmap(map)
ggmap(map) +
  geom_point(data = Coord, mapping = aes(x = longitude, y = latitude))
ggmap(map) +
  geom_label(data = Coord, mapping = aes(x = longitude, y = latitude, 
                                         label = Name), label.padding  = unit(0.2, "lines"))
  
ggsave("carte.png", width=10, height=8, dpi=300)
  





global <- map("world", "canada")
global
ggplot() +  geom_point(data=global, aes(Coord$lon, Coord$lat), color="red")
ggplot()+ geom_polygon(data = global, aes(x=lon,y=lat, group = group))+
  coord_fixed(1,3)
ggplot()+
  geom_polygon(data="global", aes(x= long, y= lat, group=group), fill= NA, color= "red")+
  coord_fixed(1,3)
gg1<- ggplot()+
  geom_polygon(data=global, aes(x=long, y= lat, group=group), fill="green", color="blue")+
  coord_fixed(1,3)
gg1

gg1+
  geom_point(data=Coord, aes(lon, lat), colour = "red", size = 1)+
  ggtitle("Quebec Map")+
  geom_text_repel(data= Coord, aes(lon, lat, label = cities))+xlim(0,150)+ ylim(0,100)

