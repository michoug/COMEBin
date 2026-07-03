from scripts.unitem_profile import Profile, make_sure_path_exists
import copy
import os

from collections import defaultdict
import pandas as pd

from scripts.unitem_common import read_bins
from scripts.unitem_markers import Markers
from scripts.gen_bins_from_tsv import gen_bins_with_cluster_ids
from scripts.unitem_defaults import CHECKM2_DIR, CHECKM2_QUALITY_REPORT
from filter_small_bins import filter_small_bins

from typing import List, Optional, Union, Dict


def read_bins_nosequences(bin_dirs):
    """Read contig-to-cluster assignments from TSV result files."""

    bins = defaultdict(lambda: defaultdict(set))

    contigs_in_bins = defaultdict(lambda: {})
    for method_id, bin_dir in bin_dirs.items():
        Header = pd.read_csv(bin_dir, sep='\t', nrows=1)
        cluster_ids = pd.read_csv(bin_dir, sep='\t',header=None, usecols=range(1, Header.shape[1])).values[:, 0]
        namelist = pd.read_csv(bin_dir, sep='\t', header=None, usecols=range(1)).values[:, 0]

        for i in range(len(namelist)):
            bins[method_id][cluster_ids[i]].add(namelist[i])
            contigs_in_bins[namelist[i]][method_id] = cluster_ids[i]

    return bins, contigs_in_bins


def get_bin_quality(quality_by_method: Dict[str, Dict[str, tuple]],
                    methods_sorted: List[str]):
    """
    Count high-quality bins for each clustering method using CheckM2 quality scores.

    :param quality_by_method: {method_id: {bin_name: (completeness, contamination)}}
    :param methods_sorted: List of method IDs sorted in a specific order.
    :return: A tuple of (bin_quality_dict, best_method).
    """
    bin_quality_dict = defaultdict(lambda: {})
    sum_list = []
    sumcont5_list = []

    for method_id in methods_sorted:
        quality = quality_by_method.get(method_id, {})

        num_5010 = 0
        num_7010 = 0
        num_9010 = 0
        num_505 = 0
        num_705 = 0
        num_905 = 0

        for comp, cont in quality.values():
            if comp > 50 and cont < 10:
                num_5010 += 1
            if comp > 70 and cont < 10:
                num_7010 += 1
            if comp > 90 and cont < 10:
                num_9010 += 1
            if comp > 50 and cont < 5:
                num_505 += 1
            if comp > 70 and cont < 5:
                num_705 += 1
            if comp > 90 and cont < 5:
                num_905 += 1

        bin_quality_dict[method_id]['num_5010'] = num_5010
        bin_quality_dict[method_id]['num_7010'] = num_7010
        bin_quality_dict[method_id]['num_9010'] = num_9010
        bin_quality_dict[method_id]['num_505'] = num_505
        bin_quality_dict[method_id]['num_705'] = num_705
        bin_quality_dict[method_id]['num_905'] = num_905
        bin_quality_dict[method_id]['sum'] = num_5010 + num_7010 + num_9010 + num_505 + num_705 + num_905
        bin_quality_dict[method_id]['sum_cont5'] = num_505 + num_705 + num_905

        sum_list.append(bin_quality_dict[method_id]['sum'])
        sumcont5_list.append(bin_quality_dict[method_id]['sum_cont5'])

    sum_max = max(sum_list)
    sum_max_method = []

    sumcont5_remain = []

    for i in range(len(methods_sorted)):
        if sum_list[i] == sum_max:
            sum_max_method.append(methods_sorted[i])
            sumcont5_remain.append(sumcont5_list[i])
    if len(sum_max_method) == 1:
        best_method = sum_max_method[0]
    else:
        best_method = sum_max_method[sumcont5_remain.index(max(sumcont5_remain))]

    return bin_quality_dict, best_method


def savecontigs_with_high_bin_quality(orig_bins: Dict[str, Dict], quality_by_bin: Dict[str, tuple],
                                      best_method: str, outpath: str):
    """
    Save contigs belonging to high-quality bins to text files.

    :param orig_bins: {method_id: {cluster_id: set(contig_names)}}
    :param quality_by_bin: {bin_name: (completeness, contamination)} for best_method
    :param best_method: The best clustering method.
    :param outpath: Directory where output text files are written.
    """
    bin_count_5010 = 0
    bin_count_5005 = 0
    with open(outpath + '/' + best_method + '5010_res.txt', 'w') as f1:
        with open(outpath + '/' + best_method + '5005_res.txt', 'w') as f2:
            for bin_id in orig_bins[best_method]:
                comp, cont = quality_by_bin.get(str(bin_id), (0.0, 0.0))
                if comp > 50 and cont < 10:
                    for key in orig_bins[best_method][bin_id]:
                        f1.write(key + '\t' + str(bin_count_5010) + '\n')
                    bin_count_5010 += 1

                if comp > 50 and cont < 5:
                    for key in orig_bins[best_method][bin_id]:
                        f2.write(key + '\t' + str(bin_count_5005) + '\n')
                    bin_count_5005 += 1


def write_estimated_bin_quality(bin_quality_dict, output_file):
    fout = open(os.path.join(output_file), 'w')
    fout.write('Binning_method\tnum_5010\tnum_7010\tnum_9010\tnum_505\tnum_705\tnum_905\tsum\tsum_cont5\n')
    for method_id in bin_quality_dict:
        fout.write(method_id + '\t' + str(bin_quality_dict[method_id]['num_5010']) + '\t'
                   + str(bin_quality_dict[method_id]['num_7010']) + '\t'
                   + str(bin_quality_dict[method_id]['num_9010']) + '\t'
                   + str(bin_quality_dict[method_id]['num_505']) + '\t'
                   + str(bin_quality_dict[method_id]['num_705']) + '\t'
                   + str(bin_quality_dict[method_id]['num_905']) + '\t'
                   + str(bin_quality_dict[method_id]['sum']) + '\t'
                   + str(bin_quality_dict[method_id]['sum_cont5']) + '\n')

    fout.close()


def run_checkm2_on_bins(bins_dir: str, checkm2_out_dir: str, num_threads: int) -> bool:
    """Run CheckM2 predict on a directory of bin FASTA files.

    :param bins_dir: Directory containing bin FASTA files (*.fa).
    :param checkm2_out_dir: Output directory for CheckM2.
    :param num_threads: Number of threads to use.
    :return: True if the quality_report.tsv was produced, False otherwise.
    """
    quality_report = os.path.join(checkm2_out_dir, CHECKM2_QUALITY_REPORT)
    if os.path.exists(quality_report):
        return True

    if not os.path.exists(bins_dir) or not os.listdir(bins_dir):
        return False

    make_sure_path_exists(checkm2_out_dir)
    cmd = ('checkm2 predict --threads %d --input %s --extension fa '
           '--output-directory %s' % (num_threads, bins_dir, checkm2_out_dir))
    os.system(cmd)
    return os.path.exists(quality_report)


def estimate_bins_quality_nobins(contig_file: str, res_path: str, num_threads: int,
                                 ignore_kmeans_res: bool = False) -> str:
    """
    Estimate the quality of bins for each clustering result using CheckM2.

    For each TSV clustering result in res_path, this function:
      1. Creates a per-cluster-ID bin directory (``<result>_checkm2_bins/``) if absent.
      2. Runs CheckM2 on that directory (``<result>_checkm2/``) if not already done.
      3. Counts high-quality bins (completeness/contamination thresholds).

    :param contig_file: Path to the assembly FASTA used as input to COMEBin.
    :param res_path: Directory containing clustering result TSV files.
    :param num_threads: Number of threads for CheckM2.
    :param ignore_kmeans_res: If True, skip TSV files whose names start with 'weight'.
    :return: Filename of the best clustering result TSV.
    """
    markers = Markers()

    filenames = os.listdir(res_path)
    namelist = []
    for filename in filenames:
        if filename.endswith('.tsv'):
            if ignore_kmeans_res:
                if not filename.startswith('weight'):
                    namelist.append(filename)
            else:
                namelist.append(filename)

    namelist.sort()

    bin_dirs = {}
    for res in namelist:
        bin_dirs[res] = (res_path + res)

    bins, contigs_in_bins = read_bins_nosequences(bin_dirs)
    methods_sorted = sorted(bins.keys())
    orig_bins = copy.deepcopy(bins)

    # For each clustering result, create bins with cluster IDs as filenames
    # then run CheckM2 to obtain per-bin quality scores.
    quality_by_method = {}
    for res in namelist:
        tsv_path = res_path + res
        checkm2_bins_dir = tsv_path + '_checkm2_bins'
        checkm2_out_dir = tsv_path + '_checkm2'
        quality_report = os.path.join(checkm2_out_dir, CHECKM2_QUALITY_REPORT)

        if not os.path.exists(checkm2_bins_dir):
            gen_bins_with_cluster_ids(contig_file, tsv_path, checkm2_bins_dir)

        run_checkm2_on_bins(checkm2_bins_dir, checkm2_out_dir, num_threads)

        quality_by_method[res] = markers.read_quality_report(quality_report)

    bin_quality_dict, best_method = get_bin_quality(quality_by_method, methods_sorted)

    savecontigs_with_high_bin_quality(orig_bins, quality_by_method.get(best_method, {}),
                                      best_method, res_path)

    output_file = res_path + 'estimate_res.txt'
    write_estimated_bin_quality(bin_quality_dict, output_file)
    return best_method


def run_get_final_result(logger, args, seed_num: int, num_threads: int = 40,
                         res_name: Optional[str] = None, ignore_kmeans_res: bool = True):
    """
    Run the final step to select the best clustering result based on CheckM2 quality.

    :param seed_num: The seed number.
    :param num_threads: The number of threads (default: 40).
    :param res_name: Unused; kept for API compatibility.
    :param ignore_kmeans_res: Whether to ignore K-means results (default: True).
    """
    logger.info("Seed_num:\t" + str(seed_num))

    best_method = estimate_bins_quality_nobins(
        args.contig_file,
        args.output_path + '/cluster_res/',
        num_threads,
        ignore_kmeans_res=ignore_kmeans_res,
    )

    logger.info('Final result:\t' + args.output_path + '/cluster_res/' + best_method)
    filter_small_bins(logger, args.contig_file, args.output_path + '/cluster_res/' + best_method, args)

