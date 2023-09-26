# file: git_05-summary.R
# author: Derek Wong, Ph.D
# date: October 14th, 2021
# modified: June 28th, 2023 by sprokopec

## Load in healthy controls
load(healthy);

sd_dist <- function(sample, mean, sd) {
	sds <- (sample - mean)/sd
	}

get.correlation <- function(x,y) {
	cor(x,y,method = "pearson", use = "complete.obs")
	}

## Generate correlations
correlations <- df.fr3 %>% ungroup() %>% group_by(sample_id) %>%
	dplyr::summarize(
		ratio = get.correlation(ratio, healthy_median$median_ratio),
		ratio_corrected = get.correlation(ratio_corrected, healthy_median$median_ratio_corrected),
		ratio_centered = get.correlation(ratio_centered, healthy_median$median_ratio_centered),
		coverage = get.correlation(coverage, healthy_median$median_coverage),
		coverage_corrected = get.correlation(coverage_corrected,
			healthy_median$median_coverage_corrected),
		coverage_centered = get.correlation(coverage_centered,
			healthy_median$median_coverage_centered),
		combined = get.correlation(combined, healthy_median$median_combined),
		combined_centered = get.correlation(combined_centered,
			healthy_median$median_combined_centered),
		nfrags = sum(nfrags),
		mode_size = unique(mode_size),
		mean_size = unique(mean_size),
		median_size = unique(median_size),
		q25_size = unique(q25_size),
		q75_size = unique(q75_size),
		hqbases_analyzed = 100*sum(nfrags)*2,
		depth = hqbases_analyzed/(504*5e6)
		);

## Generate distance from median
distance <- df.fr3 %>% 
	dplyr::summarize(seqnames = seqnames,
		arm = arm,
		start = start,
		end = end,
		gc = gc,
		ratio_dist = sd_dist(ratio, healthy_median$median_ratio, healthy_median$sd_ratio),
		ratio_corrected_dist = sd_dist(ratio_corrected, healthy_median$median_ratio_corrected,
			healthy_median$sd_ratio_corrected),
		ratio_centered_dist = sd_dist(ratio_centered, healthy_median$median_ratio_centered,
			healthy_median$sd_ratio_centered),
		coverage_dist = sd_dist(coverage, healthy_median$median_coverage,
			healthy_median$sd_coverage),
		coverage_corrected_dist = sd_dist(coverage_corrected,
			healthy_median$median_coverage_corrected, healthy_median$sd_coverage_corrected),
		coverage_centered_dist = sd_dist(coverage_centered,
			healthy_median$median_coverage_centered, healthy_median$sd_coverage_centered),
		combined_dist = sd_dist(combined, healthy_median$median_combined,
			healthy_median$sd_combined),
		combined_centered_dist = sd_dist(combined_centered,
			healthy_median$median_combined_centered, healthy_median$sd_combined_centered)
		);

## Generate summaries
summary_sd <- distance %>% ungroup() %>% group_by(sample_id) %>%
	dplyr::summarize(ratio = mean(ratio_dist),
		ratio_corrected = mean(ratio_corrected_dist),
		ratio_centered = mean(ratio_centered_dist),
		coverage = mean(coverage_dist),
		coverage_corrected = mean(coverage_corrected_dist),
		coverage_centered = mean(coverage_centered_dist),
		combined = mean(combined_dist),
		combined_centered = mean(combined_centered_dist)
		);

summary_df <- bind_rows(correlations, summary_sd);
summary_df <- as.data.frame(t(summary_df));
colnames(summary_df) <- c("correlation_to_median", "sd_from_median");

## Write summary tables
write.table(distance, file.path(outdir, paste0(id, "_5Mb_dist.txt")), sep = "\t", row.names = FALSE);
write.table(summary_df, file.path(outdir, paste0(id, "_summary.txt")), sep = "\t");

rm(correlations, distance, summary_df, sd_dist);

gc();
