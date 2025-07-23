#!/usr/bin/env python3

"""
Created on 21 Jul. 2025
"""

__author__ = "Nicolas JEANNE"
__copyright__ = "GNU General Public License"
__email__ = "jeanne.n@chu-toulouse.fr"
__version__ = "1.0.0"

import argparse
from difflib import SequenceMatcher
import logging
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def create_log(path, level):
    """Create the log as a text file and as a stream.

    :param path: the path of the log.
    :type path: str
    :param level: the level og the log.
    :type level: str
    :return: the logging:
    :rtype: logging
    """
    log_level_dict = {"DEBUG": logging.DEBUG,
                      "INFO": logging.INFO,
                      "WARNING": logging.WARNING,
                      "ERROR": logging.ERROR,
                      "CRITICAL": logging.CRITICAL}

    if level is None:
        log_level = log_level_dict["INFO"]
    else:
        log_level = log_level_dict[level]

    if os.path.exists(path):
        os.remove(path)

    logging.basicConfig(format="%(asctime)s %(levelname)s:\t%(message)s",
                        datefmt="%Y/%m/%d %H:%M:%S",
                        level=log_level,
                        handlers=[logging.FileHandler(path), logging.StreamHandler()])
    return logging


def get_samples_data(directory_annotations, data_inputs, domain, out_dir):
    """
    Extract for each sample the start, end, length, condition,
    path of the RMSF and if the sample is rejected or not because the length is different from the majority.

    :param directory_annotations: the samples annotations file directory.
    :type: str
    :param data_inputs: the input conditions data of the script.
    :type data_inputs: Pandas.Dataframe
    :param domain: the domain name.
    :type domain: str
    :param out_dir: the output directory path.
    :type out_dir: str
    :return: the extracted data of the samples on the domain.
    :rtype: Pandas.Dataframe
    """
    data = {"sample": [], "start": [], "end": [], "length": []}
    for annotation_file in os.listdir(directory_annotations):
        sample = annotation_file.replace("_domains.csv", "")
        sample_annotations = pd.read_csv(os.path.join(directory_annotations, annotation_file), sep=",")
        for _, row_annotation in sample_annotations.iterrows():
            if domain == row_annotation["domain"]:
                data["sample"].append(sample)
                data["start"].append(row_annotation["start"])
                data["end"].append(row_annotation["end"])
                data["length"].append(row_annotation["end"] - row_annotation["start"])
    df = pd.DataFrame.from_dict(data)
    df = df.sort_values(by=["length"])

    # get the conditions, the colors and the RMSF file's paths
    conditions = []
    colors = []
    paths_rmsf = []
    for _, row in df.iterrows():
        condition_found = False
        for _, row_input in data_inputs.iterrows():
            if not condition_found:
                dir_rmsf_condition = row_input["RMSF files directory"]
                for rmsf_file in os.listdir(dir_rmsf_condition):
                    if row["sample"] in rmsf_file:
                        conditions.append(row_input["condition"])
                        colors.append(row_input["color"])
                        paths_rmsf.append(os.path.join(dir_rmsf_condition, rmsf_file))
                        condition_found = True
                        break
    df = df.assign(condition=conditions)
    df = df.assign(path_rmsf=paths_rmsf)
    df = df.assign(color=colors)

    # check if the sample is rejected because the domain has a different length than the majority
    try:
        mode = df["length"].mode().values[0]
    except IndexError as ex:
        logging.error(f"domain \"{domain}\" not found in the column \"domain\" of the annotations CSV files in "
                      f"{directory_annotations}.")
        sys.exit(1)
    rejected = []
    for _, row in df.iterrows():
        if row["length"] != mode:
            logging.debug(f"sample {row['sample']} exhibits a different domain size ({row['length']} AA) for "
                          f"\"{domain}\" compared to the majority of samples ({mode} AA).")
            rejected.append(True)
        else:
            rejected.append(False)
    df = df.assign(rejected=rejected)
    path = os.path.join(out_dir, f"{domain.replace(' ', '-')}_samples_domain_sizes.csv")
    df.to_csv(path, index=False)

    # info on the rejected
    if df["rejected"].any():
        df_rejected = df[df["rejected"]]
        if len(df_rejected) != 0:
            conditions_rejected = []
            for condition in sorted(set(df_rejected["condition"])):
                count_condition_rejected = len(df_rejected[df_rejected["condition"] == condition])
                count_condition = len(df[df["condition"] == condition])
                conditions_rejected.append(f"{count_condition_rejected}/{count_condition} {condition}")
            logging.warning(f"{df['rejected'].sum()} rejected of which {', '.join(conditions_rejected)}. See the "
                            f"Domains length file by sample: {path}")
    else:
        logging.info(f"Domains length file by sample: {path}")

    # add the samples count in the condition column
    counts = {}
    df_not_rejected = df[~df.rejected]
    for condition in df_not_rejected["condition"]:
        counts[condition] = counts.get(condition, 0) + 1
    for condition in counts.keys():
        counts[condition] = f"{condition} ({counts[condition]})"
    df["condition"] = df["condition"].map(counts)

    return df


def extract_rmsf_data(data_samples_domains, set_first_residue_to_one):
    """
    Extract the domain RMSF data for the samples.

    :param data_samples_domains: the extracted data of the samples on the domain.
    :type data_samples_domains: Pandas.Dataframe
    :param set_first_residue_to_one: if the first residues of the samples domains should be set to one.
    :type set_first_residue_to_one: bool
    :return: the RMSF data samples domain.
    :rtype: Pandas.Dataframe
    """
    data = {"condition": [], "sample": [], "residue": [], "RMSF": []}
    # select only the not rejected
    data_samples_domains = data_samples_domains[~data_samples_domains.rejected]
    data_samples_domains = data_samples_domains.reset_index()
    # check if the start and end coordinates are the same, if not force the first residue as residue number 1
    if not set_first_residue_to_one and (
            not (data_samples_domains["start"] == data_samples_domains["start"][0]).all() or not (
            data_samples_domains["end"] == data_samples_domains["end"][0]).all()):
        set_first_residue_to_one = True
        logging.warning(f"Not all the samples have the same start or end coordinates for the domain, the domains "
                        f"coordinates are set to 1.")
    for _, row in data_samples_domains.iterrows():
        rmsf_sample = pd.read_csv(row["path_rmsf"], sep=",")
        rmsf_sample_domain = rmsf_sample[rmsf_sample["residues"].between(row["start"], row["end"])]
        for _, row_smp_domain in rmsf_sample_domain.iterrows():
            data["condition"].append(row["condition"])
            data["sample"].append(row["sample"])
            if set_first_residue_to_one:
                residue_init = int(row_smp_domain["residues"])
                domain_start = int(
                    data_samples_domains[data_samples_domains["sample"] == row["sample"]]["start"].iloc[0])
                residue_position = residue_init - domain_start + 1
                data["residue"].append(residue_position)
            else:
                data["residue"].append(row_smp_domain["residues"])
            data["RMSF"].append(row_smp_domain["RMSF"])

    df_rmsf = pd.DataFrame.from_dict(data)
    return df_rmsf


def create_plot_palette(df_input, data_rmsf):
    """
    Create the plot palette for the legend with samples count by condition.

    :param df_input: the script input conditions data.
    :type df_input: Pandas.Dataframe
    :param data_rmsf: the RMSF dataframe.
    :type data_rmsf: Pandas.Dataframe
    :return: the plot palette.
    :rtype: dict
    """
    conditions_with_counts = set(data_rmsf["condition"])
    palette = dict()
    for condition_with_count in conditions_with_counts:
        updated_condition = ""
        color = ""
        previous_match = ""
        for _, row in df_input.iterrows():
            match = SequenceMatcher(None, condition_with_count, row["condition"]).find_longest_match()
            if match:
                str_match = condition_with_count[match.a:match.a + match.size]
                if len(str_match) > len(previous_match):
                    updated_condition = condition_with_count
                    color = row["color"]
                    previous_match = str_match
        if not updated_condition:
            logging.error(f"no substring match for the conditions in the input file "
                          f"{', '.join(df_input['condition'].to_list())}  and the updated condition with "
                          f"the sample count: {condition_with_count}")
            sys.exit(1)
        palette[updated_condition] = color

    return palette


def plot_rmsf(src, out_dir, md_time, domain, palette, fmt):
    """
    Plot the RMSF samples data.

    :param src: the RMSF dataframe.
    :type src: Pandas.Dataframe
    :param out_dir: the output directory path.
    :type: str
    :param md_time: the molecular dynamics run time.
    :type: int
    :param domain: the domain name.
    :type: str
    :param palette: the plot palette.
    :type palette: dict
    :param fmt: the plot image format.
    :type; str
    """
    rmsf_ax = sns.lineplot(data=src, x="residue", y="RMSF", hue="condition", palette=palette)
    plot = rmsf_ax.get_figure()
    plt.suptitle(f"RMSF {domain} ({md_time} ns)", fontsize="large", fontweight="bold")
    plt.xlabel("residue", fontweight="bold")
    plt.ylabel(f"RMSD (\u212B)", fontweight="bold")
    out_path_plot = os.path.join(out_dir, f"{domain.replace(' ', '-')}_{md_time}-ns_RMSF.{fmt}")
    plot.savefig(out_path_plot)
    logging.info(f"RMSF plot for {domain} domain: {os.path.abspath(out_path_plot)}")


if __name__ == "__main__":
    descr = f"""
    {os.path.basename(__file__)} v. {__version__}

    Created by {__author__}.
    Contact: {__email__}
    {__copyright__}

    Distributed on an "AS IS" basis without warranties or conditions of any kind, either express or implied.

    Aggregate the samples RMSF data in one plot to compare between various conditions using the RMSF CSV files from 
    the RMS analysis (https://github.com/njeanne/rms) and the domains annotations file. The plot will draw the mean and 
    95% confidence interval for the samples belonging to the same condition.
    """
    parser = argparse.ArgumentParser(description=descr, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--out", required=True, type=str, help="the path to the output directory.")
    parser.add_argument("-a", "--annotations", required=True, type=str,
                        help="the path to the annotations CSV files directory.")
    parser.add_argument("-t", "--md-time", required=True, type=int,
                        help="the molecular dynamics duration in nanoseconds.")
    parser.add_argument("-d", "--domain", required=True, type=str,
                        help="Name of one domain of the annotation files.")
    parser.add_argument("-r", "--reset-domain-positions", required=False, action="store_true",
                        help="Set the domain first position to 1.")
    parser.add_argument("-x", "--format", required=False, default="svg",
                        choices=["eps", "jpg", "jpeg", "pdf", "pgf", "png", "ps", "raw", "svg", "svgz", "tif", "tiff"],
                        help="the output plots format: 'eps': 'Encapsulated Postscript', "
                             "'jpg': 'Joint Photographic Experts Group', 'jpeg': 'Joint Photographic Experts Group', "
                             "'pdf': 'Portable Document Format', 'pgf': 'PGF code for LaTeX', "
                             "'png': 'Portable Network Graphics', 'ps': 'Postscript', 'raw': 'Raw RGBA bitmap', "
                             "'rgba': 'Raw RGBA bitmap', 'svg': 'Scalable Vector Graphics', "
                             "'svgz': 'Scalable Vector Graphics', 'tif': 'Tagged Image File Format', "
                             "'tiff': 'Tagged Image File Format'. Default is 'svg'.")
    parser.add_argument("-l", "--log", required=False, type=str,
                        help="the path for the log file. If this option is skipped, the log file is created in the "
                             "output directory.")
    parser.add_argument("--log-level", required=False, type=str,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="set the log level. If the option is skipped, log level is INFO.")
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument("input", type=str,
                        help="the path to the CSV, comma separated with the following header: condition, RMSF, "
                             "annotation and color. The file first column is the condition, the second column the path "
                             "of the directory containing the RMSF analysis files, the third column the domains "
                             "annotation CSV files directory path and the fourth column is the color.")
    args = parser.parse_args()

    # create output directory if necessary
    os.makedirs(args.out, exist_ok=True)
    # create the logger
    if args.log:
        log_path = args.log
    else:
        log_path = os.path.join(args.out, f"{os.path.splitext(os.path.basename(__file__))[0]}.log")
    create_log(log_path, args.log_level)

    logging.info(f"version: {__version__}")
    logging.info(f"CMD: {' '.join(sys.argv)}")

    input_df = pd.read_csv(args.input, sep=",")
    samples_domains = get_samples_data(args.annotations, input_df, args.domain, args.out)
    rmsf_extracted = extract_rmsf_data(samples_domains, args.reset_domain_positions)
    plot_palette = create_plot_palette(input_df, rmsf_extracted)
    plot_rmsf(rmsf_extracted, args.out, args.md_time, args.domain, plot_palette, args.format)