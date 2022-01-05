targets::tar_load(cee_dive_response_probs)
targets::tar_load(cee_dive_response_summaries)

all_samples = do.call(cbind,lapply(cee_dive_response_probs, function(x) x$samples))

ppd_avg = rowMeans(all_samples[,1:4] * 5 / 60)

colMeans(all_samples * 5 / 60)

observed = mean(c(
  cee_dive_response_summaries$exposure_contexts$`data/sattag/ZcTag093_DUML_series_20191021.csv`$time_to_deep,
  cee_dive_response_summaries$exposure_contexts$`data/sattag/ZcTag095_DUML_series_20200703.csv`$time_to_deep,
  cee_dive_response_summaries$exposure_contexts$`data/sattag/ZcTag096_DUML_series_20200703.csv`$time_to_deep,
  cee_dive_response_summaries$exposure_contexts$`data/sattag/ZcTag096_DUML_series_20200703.csv`$time_to_deep
) * 5 / 60)

plot(density(ppd_avg))
abline(v = observed, lty = 3)

ecdf(ppd_avg)(observed)
# equivalent p-value computation to line 20
mean(ppd_avg <= observed)

mean(ppd_avg)



for(animalId in 1:4) {
  plot(density((all_samples * 5 / 60)[,animalId]), main = animalId)
  abline(v=cee_dive_response_summaries$pvals$time_to_deep[animalId] * 5 / 60, lty=3)
}


# TODO: Open question is how to figure out individual "takes" for each animal 
#  here.  Really, as Rob suggests, the goal is to estimate parameters that can 
#  help build a dose-response curve that can be used in simulations.  The 
#  challenge, of course, is that the model I developed doesn't explicitly 
#  estimate response parameters at the population level because there are so 
#  few replicates of movement under response.  
# 
#  well, if we think about a conditional test, we might get somewhere.  the idea 
#  is similar to how we might decompose the test statistic for a chi-sq test.
#
#  so, what if we 

# TODO: send rob the programming agile commentary youtube channel

# TODO: also maybe send rob other programming/software dev resources