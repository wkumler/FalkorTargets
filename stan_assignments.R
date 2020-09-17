# Script to assign feature numbers to each standard for the Falkor data
# Press Alt-O to collapse headings neatly


# Setup things ----
library(data.table)
library(tidyverse)
library(pbapply)
options(dplyr.summarise.inform=F)


# Functions ----
# stan_guesser()
# Custom function that tries to resolve discrepancies between various metrics ## Can you be more specific here?
# Basically a bunch of if/else statements
# Highest priority: isotope validated (isotope_choice)
# Next highest: if two peaks are found and two standards should be found, use RT
#               ordering to match them up (match_choice)
# Next highest: if peaks are dramatically larger in the correct mix and smaller
#               in the other mix, assume it's the larger one (mix_choice)
#               sometimes produces multiple, thus later area and RT matching
# Next highest: if one peak is larger than the others, assume it's the
#               standard (area_choice)
# Lowest: RT matching based on the number in the standards list
#
# Returns a character vector, usually a single peak but can be multiple or none
# concatenated with a semicolon
stan_guesser <- function(isotope_choice, mix_choice, match_choice, area_choice, rt_choice){
  ## This is usually a good place to define all your variables.
  if(all(is.na(c(isotope_choice, mix_choice, match_choice, area_choice, rt_choice)))){
    return("No peaks found")
  }
  if(!is.na(isotope_choice)){
    return(isotope_choice)
  }
  if(!is.na(match_choice)) {
    return(match_choice)
  }
  ## You could consider writing more documentation on what is happening in each step,
  ## or include it all at the top.
  if(is.na(mix_choice)){
    mix_options <- NA
  } else if(!nchar(mix_choice)){
    mix_options <- NA
  } else {
    mix_options <- unlist(strsplit(mix_choice, split = "; "))
  }
  if(length(mix_options)==1&!all(is.na(mix_options))){
    return(mix_options)
  }
  
  if(length(mix_options)>1&!all(is.na(mix_options))){
    if(rt_choice==area_choice){
      return(rt_choice)
    }
    if(rt_choice%in%mix_options&!area_choice%in%mix_options){
      return(rt_choice)
    }
    if(area_choice%in%mix_options&!rt_choice%in%mix_options){
      return(rt_choice)
    }
    if(area_choice%in%mix_options&rt_choice%in%mix_options){
      return(paste(c(area_choice, rt_choice), collapse = "; "))
    }
  }
  if(rt_choice==area_choice){
    return(rt_choice)
  }
}

# pmppm()
# Takes in a mass and produces a ppm window, lower bound of 200m/z = 0.002 Da
pmppm <- function(mass, ppm=4){
  if(mass<200){
    as.numeric(mass)+(c(-ppm, ppm)*200/1000000)
  } else {
    c(mass*(1-ppm/1000000), mass*(1+ppm/1000000))
  }
}



# Read in data ----
all_stans <- read.csv("falkor_stans.csv") %>% 
  filter(polarity=="pos") %>%
  mutate(mz=as.numeric(mz)) %>%
  add_row(compound_type="Custom", 
          compound_name="Pyroglutamic acid from glutamine",
          formula="C5H7NO3", rt=10, mz=129.042593+1.007276,
          ionization_form="[M+H-NH3]", charge=1, kegg_id=NA,
          polarity="pos", date_added=NA, mix="Mix1") %>%
  add_row(compound_type="Custom", 
          compound_name="Pyroglutamic acid from glutamate",
          formula="C5H7NO3", rt=11, mz=129.042593+1.007276,
          ionization_form="[M+H-NH3]", charge=1, kegg_id=NA,
          polarity="pos", date_added=NA, mix="Mix1")

addiso_peaks <- read.csv("addiso_peaks.csv")
real_peaks <- read.csv("real_peaks.csv")
all_peaks <- real_peaks %>%
  select(feature, mz, rt, into, file_name, M_area) %>%
  rbind(addiso_peaks %>% select(feature, mz, rt, into, file_name, M_area)) %>%
  arrange("feature", "file_name")
all_features <- all_peaks %>%
  group_by(feature) %>%
  summarise(mzmed=median(mz), rtmed=median(rt), avgarea=mean(M_area)) 



# Assign automatically ----
# For each standard
# Grab all the features that it could be based on 5ppm m/z window
# If there aren't any features, return NA for everything

# If there are features, calculate 5 metrics for each one:
# Isotope choice:
#   If the compound has an isotope, use the peak closest to the isotope peak in RT
#     (i.e. another compound exists with the same name + C13, N15, or O18)
#   If the compound IS an isotope (i.e. has C13, N15, or O18 in the name)
#     find the two peaks closest in RT and assume those are the isotopologues
#   Otherwise, NA
# Match choice:
#   If there are two standards expected in a mix and two features found,
#   assume that the two found peaks are the standards and use rank-order RT matching
#   If more or fewer peaks found, NA
# Mix choice:
#   Generally, which peaks are bigger in the mix they've been added to?
#   Calculate average z-statistic between the mixes after controlling for H2O vs Matrix
#   Return all feature numbers with a z-score above 10 (arbitrarily chosen)
# Area choice:
#   Which peak has the largest area?
# RT choice:
#   Which peak is closest in RT to the value in the standards list?

# Then use stan_guesser to resolve discrepancies based on the priority hierarchy

stan_assignments <- all_stans %>%
  mutate(mz=as.numeric(mz)) %>%
  split(.$compound_name) %>%
  pblapply(function(stan_data){
    # dput(stan_data)
    possible_stan_features <- all_features %>%
      filter(mzmed %between% pmppm(as.numeric(stan_data["mz"]))) %>%
      arrange(rtmed)
    if(!nrow(possible_stan_features)){
      return(rep(NA, 5))
    }
    
    # If the compound has an isotopologue, use RT matching
    # If the compound IS an isotopologue, same deal
    isotopologue <- stan_data$compound_name %>%
      paste0(", [15N\\d|13C\\d|2H\\d]") %>%
      grep(x = all_stans$compound_name, value = TRUE)
    if(length(isotopologue)==1){
      isotopo_data <- all_stans %>% 
        filter(compound_name==isotopologue)
      isotopo_features <- all_features %>%
        filter(mzmed %between% pmppm(as.numeric(isotopo_data$mz), 5)) 
      ## Small syntax thing, but I think it's a lot easier to read with the spaces around the 
      ## %between% function.
      if(nrow(isotopo_features)==1){
        isotope_choice <- possible_stan_features %>% 
          mutate(rtdiff=abs(rtmed-isotopo_features$rtmed)) %>%
          arrange(rtdiff) %>%
          slice(1) %>%
          pull(feature)
      } else {
        isotope_choice <- NA
      }
    } else if(grepl(pattern = ", 15N|13C|2H", x = stan_data$compound_name)) {
      isotopologue <- gsub(", 15N.*|, 13C.*|, 2H.*", "", stan_data$compound_name)
      isotopo_data <- all_stans %>% 
        filter(compound_name==isotopologue|
                 compound_name==gsub("^D", "", isotopologue)|
                 compound_name==gsub("^L", "", isotopologue))
      isotopo_features <- all_features %>%
        filter(mzmed %between% pmppm(as.numeric(isotopo_data$mz), 5)) 
      isotope_choice <- possible_stan_features %>% 
        left_join(isotopo_features, by=character()) %>%
        select(-starts_with("mzmed")) %>%
        mutate(rtdiff=abs(rtmed.x-rtmed.y)) %>%
        arrange(rtdiff) %>%
        slice(1) %>%
        pull(feature.x)
    } else {
      isotope_choice <- NA
    }
    
    # If there's a peak that differs between the mixes
    if(is.na(stan_data$mix)){
      mix_choice <- NA
    } else {
      stan_peaks <- all_peaks %>%
        filter(feature%in%possible_stan_features$feature) %>%
        select(feature, mz, rt, into, file_name, M_area) %>%
        filter(grepl("Std", .$file_name)) %>%
        filter(!grepl("H2OinMatrix", .$file_name)) %>%
        mutate(correct_mix=grepl(file_name, pattern = stan_data["mix"])) %>%
        mutate(stan_type=str_extract(file_name, "InH2O|InMatrix"))
      
      # ggplot(stan_peaks) +
      #   geom_boxplot(aes(x=stan_type, y=M_area, color=correct_mix)) +
      #   facet_wrap(~feature, scales = "free_y")
      
      if(!nrow(stan_peaks)){
        return(NA)
      }
      
      mix_peaks <- stan_peaks %>% 
        group_by(feature, correct_mix, stan_type) %>%
        summarise(avgarea=mean(M_area), sdarea=sd(M_area)) %>%
        right_join(expand.grid(
          unique(.$feature),
          unique(.$correct_mix),
          unique(.$stan_type)
        ) %>% `names<-`(c("feature", "correct_mix", "stan_type")),
        by=c("feature", "correct_mix", "stan_type")) %>%
        mutate(avgarea=ifelse(is.na(avgarea), 0, avgarea)) %>%
        mutate(sdarea=ifelse(is.na(sdarea), 0, sdarea))
        
      if(length(unique(mix_peaks$feature))==1){
        mix_choice <- unique(mix_peaks$feature)
      } else {
        mix_choice <- mix_peaks %>%
          split(interaction(.$feature, .$stan_type)) %>%
          lapply(function(v){
            diff <- (v$avgarea[v$correct_mix] - v$avgarea[!v$correct_mix])/mean(v$sdarea)
            data.frame(feature=unique(v$feature), 
                       stan_type=unique(v$stan_type),
                       diff_degree=diff)
          }) %>%
          do.call(what = rbind) %>% `rownames<-`(NULL) %>%
          group_by(feature) %>%
          summarise(correct_mix_peak=mean(diff_degree)) %>%
          filter(correct_mix_peak>10) %>%
          pull(feature) %>%
          paste(collapse = "; ")
      }
    }

    # If there's one peak much closer in RT to expected than the others
    expected_rt <- stan_data$rt*60
    rt_choice <- possible_stan_features %>% 
      mutate(rtdiff=abs(rtmed-expected_rt)) %>%
      arrange(rtdiff) %>%
      slice(1) %>%
      pull(feature)
    
    # If there's same number of features as expected peaks, assume 1:1 and order by RT
    possible_other_stans <- all_stans %>%
      filter(mz %between% pmppm(stan_data$mz, 5)) %>% 
      arrange(rt)
    if(nrow(possible_stan_features)==nrow(possible_other_stans)){
      match_choice <- possible_stan_features %>%
        arrange(rtmed) %>%
        cbind(possible_other_stans) %>%
        select(feature, compound_name) %>%
        filter(compound_name==stan_data$compound_name) %>%
        pull(feature)
    } else {
      match_choice <- NA
    }
    
    stan_peaks <- all_peaks %>%
      filter(feature%in%possible_stan_features$feature) %>%
      select(feature, mz, rt, into, file_name, M_area) %>%
      filter(grepl("Std", .$file_name))
    
    if(!nrow(stan_peaks)){
      area_choice <- NA
    } else if(nrow(possible_stan_features)==nrow(possible_other_stans)&
              nrow(possible_stan_features)>1) {
      area_choice <- NA
    } else {
      area_choice <- stan_peaks %>% 
        group_by(feature) %>%
        summarise(avgarea=mean(M_area)) %>%
        arrange(desc(avgarea)) %>%
        slice(1) %>%
        pull(feature)
    }
    
    best_guess <- stan_guesser(isotope_choice, mix_choice, match_choice, area_choice, rt_choice)
    
    return(data.frame(
      best_guess=best_guess,
      isotope_validated=isotope_choice,
      rt_matchup=match_choice,
      mix_matched=mix_choice,
      closer_rt=rt_choice,
      area_choice=area_choice
    ))
  }) %>%
  do.call(what=rbind) %>%
  mutate(compound_name=rownames(.)) %>%
  select(compound_name, everything())


# Summary checks ----
# Check on internal standard assignments
stan_assignments %>%
  filter(!is.na(.$isotope_validated))


# Check on features that are dual-assigned - these should be resolved by manual assignment
stan_assignments[(duplicated(stan_assignments$best_guess, fromLast=TRUE)|
                   duplicated(stan_assignments$best_guess))&
                   !is.na(stan_assignments$best_guess), ] %>%
  split(.$best_guess)

## Is the code below a direct result from checking for manual assignments? It may be good
## to try and automate that process. Filter by those compounds that come up in the duplicate output
## and make a function that will do the work below with manual inputs. 

# Manual assignments ----
leucine_data <- all_stans %>% filter(compound_name=="L-Leucine")
leucine <- all_features %>%
  filter(mzmed %between% pmppm(leucine_data$mz, 5)) %>%
  arrange(rtmed) %>%
  slice(2)
stan_assignments[stan_assignments$compound_name=="L-Leucine",] <- 
  c("L-Leucine", leucine$feature, rep("Manual", ncol(stan_assignments)-2))

N6_lysine_data <- all_stans %>% filter(compound_name=="N6-Acetyl-L-lysine")
N6_lysine <- all_features %>%
  filter(mzmed %between% pmppm(N6_lysine_data$mz, 5)) %>%
  arrange(desc(avgarea)) %>%
  slice(1)
stan_assignments[stan_assignments$compound_name=="N6-Acetyl-L-lysine",] <- 
  c("N6-Acetyl-L-lysine", N6_lysine$feature, rep("Manual", ncol(stan_assignments)-2))

allopurinol_data <- all_stans %>% filter(compound_name=="Allopurinol")
allopurinol <- all_features %>%
  filter(mzmed %between% pmppm(allopurinol_data$mz, 5)) %>%
  arrange(desc(avgarea)) %>%
  slice(1)
hypoxanthine <- all_features %>%
  filter(mzmed %between% pmppm(allopurinol_data$mz, 5)) %>%
  arrange(desc(avgarea)) %>%
  slice(2)
stan_assignments[stan_assignments$compound_name=="Allopurinol",] <- 
  c("Allopurinol", allopurinol$feature, rep("Manual", ncol(stan_assignments)-2))
stan_assignments[stan_assignments$compound_name=="Hypoxanthine",] <- 
  c("Hypoxanthine", hypoxanthine$feature, rep("Manual", ncol(stan_assignments)-2))


# Write out ----
write.csv(x = stan_assignments, file = "stan_assignments.csv")
