import logging
import os
import sys
import time

# import pkg_resources # replacing deprecated package
if sys.version_info >= (3, 7):
    import importlib.resources as importlib_resources
else:
    import importlib_resources
from contextlib import ExitStack
import atexit

import pandas as pd

from deepec import utils
from deepec import ec_prediction_dl
from deepec import ec_prediction_seq
from deepec import __version__


def main():
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

    start = time.time()
    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    file_manager = ExitStack()
    atexit.register(file_manager.close)

    deep_ec4d_oh_ref = importlib_resources.files('deepec.data') / 'DeepEC_4d.h5'
    deep_ec4d_oh = file_manager.enter_context(importlib_resources.as_file(deep_ec4d_oh_ref))
    mb4d_oh_ref = importlib_resources.files('deepec.data') / 'multilabelbinarizer_dl_4digit.pkl'
    mb4d_oh = file_manager.enter_context(importlib_resources.as_file(mb4d_oh_ref))

    deep_ec3d_oh_ref = importlib_resources.files('deepec.data') / 'DeepEC_3d.h5'
    deep_ec3d_oh = file_manager.enter_context(importlib_resources.as_file(deep_ec3d_oh_ref))
    mb3d_oh_ref = importlib_resources.files('deepec.data') / 'multilabelbinarizer_dl_3digit.pkl'
    mb3d_oh = file_manager.enter_context(importlib_resources.as_file(mb3d_oh_ref))

    deep_ec_enzyme_ref = importlib_resources.files('deepec.data') / 'Binary_class.h5'
    deep_ec_enzyme_oh = file_manager.enter_context(importlib_resources.as_file(deep_ec_enzyme_ref))
    ref_db_file_ref = importlib_resources.files('deepec.data') / 'ref_low_seq_db.dmnd'
    ref_db_file = file_manager.enter_context(importlib_resources.as_file(ref_db_file_ref))

    parser = utils.argument_parser(version=__version__)

    options = parser.parse_args()
    fasta_file = options.fasta_file
    output_dir = options.output_dir
    threshold = options.threshold
    threshold = float(threshold)

    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    # One-hot encoding
    temp_file = '%s/SwissProt_input.csv.gz' % output_dir
    temp_fasta_file = '%s/temp_fasta.fa' % output_dir

    ec_prediction_dl.preprocessing(fasta_file, temp_fasta_file)
    ec_prediction_dl.run_one_hot_encoding(temp_fasta_file, temp_file)

    temp_df = pd.read_csv(temp_file, compression='gzip', index_col=0)

    enzyme_output_file = '%s/Enzyme_prediction.txt' % output_dir
    ec_prediction_dl.predict_dl(temp_df, enzyme_output_file, deep_ec_enzyme_oh)

    dl4_output_file = '%s/4digit_EC_prediction.txt' % output_dir
    ec_prediction_dl.predict_dl(temp_df, dl4_output_file, deep_ec4d_oh, mb4d_oh, threshold)

    dl3_output_file = '%s/3digit_EC_prediction.txt' % output_dir
    ec_prediction_dl.predict_dl(temp_df, dl3_output_file, deep_ec3d_oh, mb3d_oh, threshold)

    target_fasta_file = '%s/low_seq_candidates.txt' % output_dir
    utils.merge_dl_prediction_results(output_dir, dl4_output_file, dl3_output_file, enzyme_output_file, fasta_file,
                                      target_fasta_file)

    output_file = '%s/Blastp_result.txt' % output_dir
    ec_prediction_seq.predict_ec(ref_db_file, target_fasta_file, output_file, output_dir)

    utils.merge_dl_seq_results(output_dir)
    utils.remove_files(output_dir)
    utils.move_files(output_dir)
    logging.info(time.strftime("Elapsed time %H:%M:%S", time.gmtime(time.time() - start)))


if __name__ == '__main__':
    main()
