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
