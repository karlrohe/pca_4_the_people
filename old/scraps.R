# https://github.com/karlrohe/pca_4_the_people

source("pca_4_the_people.R")
library(nycflights13)
embeddings = pca_count(1 ~ (month & day)*(dest), 
                       flights, 
                       k = 6)


airport_dat = embeddings$column_features %>% 
  left_join(airports %>% select(dest=faa, lat,lon)) %>% 
  select(lat, lon, contains("_col")) %>% 
  pivot_longer(contains("pc_"),
               names_to = "pc_dimension", values_to = "loadings") %>% 
  filter(pc_dimension == "pc_3_columns") %>% 
  drop_na()


library(maps)
usa_map <- map_data("state")
p <- ggplot() + 
  geom_polygon(data = usa_map, aes(x = long, y = lat, group = group), 
               fill = "white", color = "black") +
  coord_fixed(1.3, xlim = c(-125, -65), ylim = c(25, 50)) 
# i'm only keeping lower 48 states, dropping Anchorage and Honolulu.


p + geom_point(data = airport_dat, aes(x = lon, y = lat, 
                                       size = abs(loadings), color = loadings)) +
  facet_wrap(~ pc_dimension)  +
  scale_color_gradient2(low = "red", high = "blue", mid = "white")




embeddings = pca_count(1 ~ (month & day)*(dest), flights, k = 6)

embeddings$row_features %>% 
  mutate(date = make_date(day = day, month=month, year = 2013)) %>% 
  select(date, contains("pc_")) %>% 
  pivot_longer(contains("pc_"), names_to = "pc_dimension", values_to = "loadings") %>% 
  filter(pc_dimension == "pc_3_rows") %>% 
  ggplot(aes(x = date, y = loadings)) + geom_line()  
  # facet_wrap(~pc_dimension, scales= "free") + geom_smooth()
