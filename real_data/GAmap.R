require("data.table")
require("rgdal") # requires sp, will use proj.4 if installed
require("maptools")
require("ggplot2")
require("plyr")
require("gridExtra")
load("GA_data.Rdata")
load("GA_map.Rdata")
data_use<-GA_data
#data_use$county <- as.character(data_use$county)
#data_use$House_market_price <- data_use$House_market_price/1000
#data_use$Population <- data_use$Population/1000
names(data_use) <- c("county", "Premature.Death", "PM2.5",
                     "Unemployment", "Food.Environment.Index","Health.Care.Cost")
data_use$county <- as.character(seq(1:159))
GAmap <- fortify(mapdata)
    
ggplot(GAmap) +
    aes(long, lat, group = subregion) +
    geom_polygon(color = "black", fill = "white") +
    theme_bw()

map.county <- map_data('county')
counties   <- unique(map.county[, 5:6])
GAcounties <- counties[counties$region == "georgia", ]

map.county <-
    data.table(map.county[map.county$region == "georgia", ])
setkey(map.county, region, subregion)

rent_df <- data.table(data.frame(
    state_names = "georgia",
    county_names = GAcounties$subregion,
    rent = data_use$House_Rent
))
setkey(rent_df, state_names, county_names)
map.df <- map.county[rent_df, allow.cartesian=TRUE]




dfList <- purrr::map(c(3,5),function(x) {
                         data <- data.table(data.frame(
                                state_names = "georgia",
                                county_names = GAcounties$subregion,
                                Value = data_use[,x]))
                        setkey(data, state_names, county_names)
                        map.df <- map.county[data]
                        return(map.df)
}
)

plotList <- purrr::map2(dfList, names(data_use)[c(3,5)], function(i, j){
    ggplot(i, aes(
        x = long,
        y = lat,
        group = group,
        fill = Value
    )) + 
        geom_polygon(col = "black") + coord_map() + 
        xlab("Longitude") + ylab("Latitude") + 
        scale_fill_gradientn(colors=c('white','deepskyblue')) + theme_bw() + 
        ggtitle(j)
})

pdf("covariates_GA.pdf", width = 12)
grid.arrange(grobs = plotList, nrow = 1)
dev.off()



data <- data.table(data.frame(
    state_names = "georgia",
    county_names = GAcounties$subregion,
    Value = data_use[,2]/100))
setkey(data, state_names, county_names)
map.df <- map.county[data]

ggplot(map.df, aes(
    x = long,
    y = lat,
    group = group,
    fill = Value
)) + 
    geom_polygon(col = "black") + coord_map() + 
    xlab("Longitude") + ylab("Latitude") + 
    scale_fill_gradientn(colors=c('white','deepskyblue')) + theme_bw()+ggtitle('Premature Death')
