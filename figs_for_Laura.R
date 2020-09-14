library(tidyverse)
library(data.table)
pmppm <- function(mass, ppm=4){
  if(mass<200){
    as.numeric(mass)+(c(-ppm, ppm)*200/1000000)
  } else {
    c(mass*(1-ppm/1000000), mass*(1+ppm/1000000))
  }
}
all_stans <- read.csv("falkor_stans.csv") %>% 
  filter(polarity=="pos") %>% filter(!is.na(mix)) %>%
  arrange(as.numeric(mz), as.numeric(rt)) %>%
  mutate(mz=as.numeric(mz))
stan_MS1_data <- readRDS("stan_data.rds")
addiso_peaks <- read.csv("addiso_peaks.csv")
real_peaks <- read.csv("real_peaks.csv")
all_peaks <- real_peaks %>%
  select(feature, mz, rt, into, file_name, M_area) %>%
  rbind(addiso_peaks %>% select(feature, mz, rt, into, file_name, M_area)) %>%
  arrange("feature", "file_name")
all_features <- all_peaks %>%
  group_by(feature) %>%
  summarise(mzmed=median(mz), rtmed=median(rt), avgarea=mean(M_area))
stan_annotations <- read.csv("stan_assignments.csv")



# Proline betaine peak
mz_i <- all_stans$mz[all_stans$compound_name=="Proline betaine"]
eic <- stan_MS1_data[mz%between%pmppm(mz_i)]
eic_zeroed <- split(eic, eic$fileid) %>%
  lapply(function(x){
    rbind(
      data.frame(rt=min(x$rt), mz=median(x$mz), int=0, fileid=unique(x$fileid)),
      x,
      data.frame(rt=max(x$rt), mz=median(x$mz), int=0, fileid=unique(x$fileid))
    )
  }) %>% 
  do.call(what = rbind) %>%
  mutate(fillcol=str_extract(fileid, "Blk|Smp|Std")) %>%
  filter(fillcol=="Std") %>%
  mutate(stdtype=str_extract(fileid, "Mix[1|2]|H2O"))
gp <- ggplot() + 
  geom_polygon(data = eic_zeroed, aes(x=rt, y=int, fill=stdtype)) +
  geom_line(data = eic_zeroed, aes(x=rt, y=int, group=stdtype)) +
  scale_y_continuous(oob = scales::rescale_none, limits = c(0, max(eic_zeroed$int)/5)) +
  facet_wrap(~fileid) +
  scale_fill_viridis_d(alpha = 0.7, option = "C") +
  xlim(c(200, 700)) +
  theme_bw() +
  ggtitle("Proline betaine")
print(gp)




# Acetylcholine/TMAP
mz_i <- all_stans$mz[all_stans$compound_name=="Acetylcholine"]
eic <- stan_MS1_data[mz%between%pmppm(mz_i)]
eic_zeroed <- split(eic, eic$fileid) %>%
  lapply(function(x){
    rbind(
      data.frame(rt=min(x$rt), mz=median(x$mz), int=0, fileid=unique(x$fileid)),
      x,
      data.frame(rt=max(x$rt), mz=median(x$mz), int=0, fileid=unique(x$fileid))
    )
  }) %>% 
  do.call(what = rbind) %>%
  mutate(fillcol=str_extract(fileid, "Blk|Smp|Std")) %>%
  filter(fillcol=="Std") %>%
  mutate(stdtype=str_extract(fileid, "Mix[1|2]|H2O"))
gp <- ggplot() + 
  geom_polygon(data = eic_zeroed, aes(x=rt, y=int, fill=stdtype)) +
  geom_line(data = eic_zeroed, aes(x=rt, y=int, group=stdtype)) +
  scale_y_continuous(oob = scales::rescale_none, limits = c(0, max(eic_zeroed$int)/5)) +
  facet_wrap(~fileid) +
  scale_fill_viridis_d(alpha = 0.7, option = "C") +
  xlim(c(400, 600)) +
  theme_bw() +
  ggtitle("Acetylcholine/3-carboxy")
print(gp)
