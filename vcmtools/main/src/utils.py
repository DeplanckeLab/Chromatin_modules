import numpy as np
import os
import pandas as pd
import scipy.stats as stats


def process_different_marks_pv(
    path_to_output_directory, dataset, dict_for_chr, max_dist
):
    with open(
        os.path.join(
            path_to_output_directory,
            dataset + "_all_marks_VCM_theoretical_corr_with_p_values.txt",
        ),
        "a+",
    ) as file_with_corr:
        for peak1, peak1_values in dict_for_chr.items():
            peak1_variables = peak1.split(":")
            mark1 = peak1_variables[0]
            peak1_start = int(peak1_variables[-2])
            peak1_end = int(peak1_variables[-1])
            peak1_center = peak1_start + (peak1_end - peak1_start) / 2
            for peak2, peak2_values in dict_for_chr.items():
                peak2_variables = peak2.split(":")
                mark2 = peak2_variables[0]
                peak2_start = int(peak2_variables[-2])
                peak2_end = int(peak2_variables[-1])
                peak2_center = peak2_start + (peak2_end - peak2_start) / 2
                try:
                    if (peak1 != peak2) and (
                        abs(peak1_center - peak2_center) <= max_dist
                    ):
                        x = np.array(peak1_values)
                        y = np.array(peak2_values)
                        mask = ~np.isnan(x) & ~np.isnan(y)
                        x = x[mask]
                        y = y[mask]
                        theoretical_correlation, p_value = stats.pearsonr(x, y)
                        file_with_corr.write(
                            "\t".join(
                                [
                                    peak1,
                                    peak2,
                                    str(theoretical_correlation),
                                    str(p_value),
                                ]
                            )
                            + "\n"
                        )
                except:
                    continue


def empirical_pvalue_for_corr(dataset, input_path, output_path):
    with open(
        os.path.join(output_path, dataset + "_all_marks_empirical_corr_p_values.txt"),
        "a+",
    ) as file_with_corr:
        background_distribution = pd.read_csv(input_path, sep="\t", header=None)
        background_distribution.columns = ["peak1", "peak2", "corr", "pvalue"]
        background_distribution["peak_pair"] = list(
            zip(background_distribution["peak1"], background_distribution["peak2"])
        )

        background_mean = background_distribution["corr"].mean()
        background_std_dev = background_distribution["corr"].std()
        dict_with_corr = dict(
            zip(background_distribution["peak_pair"], background_distribution["corr"])
        )
        del background_distribution
        for peak_pair, corr in dict_with_corr.items():
            peak1, peak2 = peak_pair
            p_value = 1 - stats.norm.cdf(corr, background_mean, background_std_dev)
            file_with_corr.write(
                "\t".join([peak1, peak2, str(corr), str(p_value)]) + "\n"
            )
