# Add UI here!
# TODO : add UI to select Base_edit_2 or Indel_searcher_2
MODE = 1  # INDEL SEARCHER

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        prog="CRISPR_toolkit",
        description="Indel search program for CRISPR CAS9 & CPF1",
        epilog="SKKUGE_DEV, 2023-02-09 ~",
    )

    parser.add_argument(
        "-t",
        "--thread",
        default="0",
        type=int,
        dest="multicore",
        help="multiprocessing number, recommendation:t<16",
    )
    parser.add_argument(
        "-c",
        "--chunk_size",
        default="100000",
        type=int,
        dest="chunk_size",
        help="split FASTQ, indicates how many reads will be in a splitted file. \
            file size < 1G recommendation:10000, size > 1G recommendation:100000",
    )
    parser.add_argument(
        "-q",
        "--base_quality",
        default="20",
        dest="base_quality",
        help="NGS read base quality",
    )
    parser.add_argument(
        "--gap_open",
        default="-10",
        type=float,
        dest="gap_open",
        help="gap open: -100~0",
    )
    parser.add_argument(
        "--gap_extend",
        default="1",
        type=float,
        dest="gap_extend",
        help="gap extend: 1~100",
    )
    parser.add_argument(
        "-i",
        "--insertion_window",
        default="4",
        type=int,
        dest="insertion_window",
        help="a window size for insertions",
    )
    parser.add_argument(
        "-d",
        "--deletion_window",
        default="4",
        type=int,
        dest="deletion_window",
        help="a window size for deletions",
    )
    parser.add_argument("--pam_type", dest="pam_type", help="PAM type: Cas9 Cpf1")
    parser.add_argument(
        "--pam_pos", dest="pam_pos", help="PAM position: Forward Reverse"
    )
    parser.add_argument(
        "--python", dest="python", help="The python path including the CRISPResso2"
    )
    parser.add_argument("--user", dest="user_name", help="The user name with no space")
    parser.add_argument(
        "--project", dest="project_name", help="The project name with no space"
    )
    parser.add_argument(
        "--pickle",
        dest="pickle",
        default="False",
        help="Dont remove the pickles in the tmp folder : True, False",
    )
    parser.add_argument(
        "--split",
        dest="split",
        default="False",
        help="Dont remove the split files in the input folder : True, False",
    )
    parser.add_argument(
        "--classfied_FASTQ",
        dest="class_fastq",
        default="True",
        help="Dont remove the ClassfiedFASTQ in the tmp folder : True, False",
    )
    parser.add_argument(
        "--ednafull", dest="ednafull", help="The nucleotide alignment matrix"
    )

    args = parser.parse_args()

    if MODE == 1:
        from Indel_searcher_2 import Run_indel_searcher

        Run_indel_searcher.indel_searcher_runner(args)
