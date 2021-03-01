# Colors 
# https://www.w3schools.com/colors/colors_picker.asp
color_S <- "orange"
color_TP <- "springgreen4"
color_AI <- "#2BB07FFF"
color_SA <- "#38598CFF"
color_SI <- "yellow4" # 8b8b00



#sapply(mLBF,typeof)
#
#x <- mLBF %>% mutate_if(is.double,as.numeric)
#
#sapply(x,typeof)
#max(x$Sequences)
#max(x$TvS)
#max(x$Sequences, na.rm=T)




# grab rows with sig, rank, or per >= 2, lables rows that have 
#
#z <- p %>% filter_at(vars(sig,per,rank), any_vars(. >= 2)) %>%
#  mutate(all=case_when(sig+per+rank >=6 ~ 3,
#                       between((sig+per+rank),4,6) ~ 2,
#                       between((sig+per+rank),2,3) ~ 1,
#                       sig+per+rank <= 1 ~ 0))


# total loci for each hypothesis that have at least one issue 
y <- p %>% group_by(Locus,Hypothesis) %>% summarize(loci=n_distinct(sig))

# quick gut check that gives the same # of loci as above. 
#x <- p %>% mutate(q=sig+per+rank)
#y <- x[x$q>1,] # greater than 1 because some of sig are 1 if b and g agree and are significant
