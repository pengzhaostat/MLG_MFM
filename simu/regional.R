centroids <- readRDS("GAcentroids.rds")

centroids <- as.data.frame(centroids)


library(ggmap)
library(maps)
library(data.table)

a <- 1:159
map.county <- map_data('county')
counties   <- unique(map.county[, 5:6])
louiscounties <- counties[counties$region == "georgia", ]
age_map <- data.frame(
  state_names = louiscounties$region,
  county_names = louiscounties$subregion,
  Region = factor(a)
)
age_map <- data.table(age_map)
setkey(age_map, state_names, county_names)

map.county <-
  data.table(map.county[map.county$region == "georgia", ])
setkey(map.county, region, subregion)



map.df <- map.county[age_map]




for (i in 1:nrow(centroids)) {
  if (centroids$x[i] - 2 * centroids$y[i] < -150) {
  map.df[map.df$subregion == row.names(centroids)[i], "Region"] <- "1"
  } else if (centroids$x[i] + centroids$y[i] > -51) {
    map.df[map.df$subregion == row.names(centroids)[i], "Region"] <- "2"
  } else {
    map.df[map.df$subregion == row.names(centroids)[i], "Region"] <- "3"
  }
}


map.df$Region <- factor(map.df$Region)

map.df$Cluster <- map.df$Region
# plot inter

age <- ggplot(map.df, aes(x = long, y = lat, fill = Cluster)) +
  geom_polygon(aes(group = group),  color = "black") + 
  coord_map() +
  xlab("Longitude") + ylab("Latitude") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette = "OrRd")
age1 <- age + theme_bw()
age1
asm<-c()
for (i in 1:nrow(centroids)) {
  if (centroids$x[i] - 2 * centroids$y[i] < -150) {
    asm[i] <- 1
  } else if (centroids$x[i] + centroids$y[i] > -51) {
    asm[i] <- 2
  } else {
    asm[i] <- 3
  }
}

betaMat <- t(matrix(nrow = 159, ncol = 6, byrow = TRUE))

for (i in 1:159) {
  ## cluster 1
  betaMat[,asm == 1] <- c(2,0,1,0,4,2)
  ## cluster 2
  betaMat[,asm == 2] <- c(1,0,3,2,0,3)
  ## cluster 3
  betaMat[,asm == 3] <- c(4,1,0,3,0,1)
}
