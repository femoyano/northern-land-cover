# Plot results from LC
# Yearly mean values from Glance Land Cover data are rad and plotted

library(ggplot2)
library(lubridate)
library(dplyr)

# Inputs ---
# file <- 'gee_yearlyAreas_ABoVE_scale300.csv'
# name <- 'ABoVE_scale300'
# file <- 'gee_yearlyAreas_ABoVE_scale900.csv'
# name <- 'ABoVE_scale900'
# file <- 'gee_yearlyAreas_NA_scale900.csv'
# name <- 'NorthAmerica - scale900'
file <- 'gee_yearlyAreas_NA50_scale900.csv'
name <- 'NorthAmerica > 50N - scale900'

lc <- c('Water', 'SnowIce', 'Developed', 'Bare', 'Forest', 'Shrub', 'Herb')
datadir <- file.path('..', 'data')
df <- read.csv(file.path(datadir, file))
df$LandCover <- recode(df$LC, !!!lc)
df$area <- df$area / 1000000  # to million of sqkm
cols <- c('Water'='#4169E1', 'SnowIce'='#DCDCDC', 'Developed'='#708090', 'Bare'='#8B4513',
          'Forest'='#006400', 'Shrub'='#9ACD32', 'Herb'='#BDB76B')
df <- df %>% group_by(LC) %>% mutate(areaPct = (area - first(area)) / first(area)  * 100)
df <- df %>% group_by(LC) %>% mutate(areaDif = area - first(area))


lc1 <- c('Forest', 'Herb')
cols1 <- c('Forest'='#006400', 'Herb'='#BDB76B')
df1 <- df[df$LandCover %in% lc1,]

lc2 <- c('Water', 'Developed', 'Bare')
cols2 <- c('Water'='#4169E1', 'Developed'='#708090', 'Bare'='#8B4513')
df2 <- df[df$LandCover %in% lc2,]

lc3 <- c('Forest', 'Shrub', 'Herb')
cols3 <- c('Forest'='#006400', 'Shrub'='#9ACD32', 'Herb'='#BDB76B')
df3 <- df[df$LandCover %in% lc3,]

# Plotting function
plot_ts <- function(df, cols, var, ylab, name, sel) {
    theme_set(theme_bw())
    brks <- seq(2001, 2019)
    lbls <- brks
    fileout <- paste0('LC_2001-2019_', name, '_', sel, '_', var, '.png')
    pathout <- file.path('..', 'figures', 'time_series', fileout)
    png(filename = pathout, width = 500, height = 400)
    p <- ggplot(df, aes(x=year)) + 
        geom_line(aes(y=get(var), col=LandCover)) + 
        labs(title=paste("Land Cover:", name), 
             y=ylab, color=NULL) +  # title and caption
        scale_x_continuous(labels = lbls, breaks = brks) +
        scale_color_manual(values = cols) +  # line color
        theme(axis.text.x = element_text(angle = 90, vjust=0.5, size = 8),  # rotate x axis text
              panel.grid.minor = element_blank())  # turn off minor grid
    print(p)
    dev.off()
}

plot_ts(df = df, cols = cols, var = 'area', ylab = "Area (millions km^2)", name = name, sel = 'all')
plot_ts(df = df2, cols = cols2, var = 'areaPct', ylab = "Area Change (%)", name = name, sel = 'nonveg')
plot_ts(df = df2, cols = cols2, var = 'areaDif', ylab = "Area Change (million km^2)", name = name, sel = 'nonveg')
plot_ts(df = df3, cols = cols3, var = 'areaPct', ylab = "Area Change (%)", name = name, sel = 'veg')
plot_ts(df = df3, cols = cols3, var = 'areaDif', ylab = "Area Change (million km^2)", name = name, sel = 'veg')
