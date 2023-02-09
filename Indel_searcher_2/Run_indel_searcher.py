#!/usr/bin/env python

import logging
import os

# import cPickle as pickle
import pickle
import subprocess as sp
import sys

sys.path.insert(0, os.path.dirname(os.getcwd()))
from Core.CoreSystem import (
    InitialFolder,
    UserFolderAdmin,
    Helper,
    RunMulticore,
    CheckProcessedFiles,
)


class clsIndelSearcherRunner(UserFolderAdmin):
    """
    self.strOutputDir is inherited variable.

    """

    def __init__(self, strSample, strRef, args, InstInitFolder):
        UserFolderAdmin.__init__(
            self, strSample, strRef, args, InstInitFolder.strLogPath
        )
        self.struct_output_dir()

        self.strProjectFile = InstInitFolder.strProjectFile
        self.intChunkSize = args.chunk_size
        self.strQualCutoff = args.base_quality
        self.intInsertionWin = args.insertion_window  # Insertion window 0,1,2,3,4
        self.intDeletionWin = args.deletion_window  # Deletion window 0,1,2,3,4
        self.strPamType = (
            args.pam_type
        )  # CRISPR type : Cpf1(2 cleavages), Cas9(1 cleavage) # TODO: add more PAM types
        self.strPamPos = (
            args.pam_pos
        )  # Barcode target position : Forward (barcode + target), Reverse (target + barcode)
        self.strPickle = args.pickle
        self.strClassFASTQ = args.class_fastq
        self.strSplit = args.split
        self.strLogPath = InstInitFolder.strLogPath
        self.strRefFile = os.path.join(self.strRefDir, 'Reference.fa')

        ## Files needed in the FASTQ directory
        self.strFastqDir = "./Input/{user}/FASTQ/{project}".format(
            user=self.strUser, project=self.strProject
        )
        ## './Input/JaeWoo/FASTQ/JaeWoo_test_samples/Test_sample'
        self.strSampleDir = os.path.join(self.strFastqDir, self.strSample)

        self.strFastq_name = ""
        for strFile in os.listdir(self.strSampleDir):
            if (
                os.path.isfile(self.strSampleDir + "/" + strFile)
                and strFile.split(".")[-1] == "fastq"
            ):
                self.strFastq_name = ".".join(strFile.split(".")[:-1])
        logging.info("File name : %s" % self.strFastq_name)

        ## './Input/JaeWoo/FASTQ/JaeWoo_test_samples/Test_sample/eCas9_rep1_D4.ext.fastq'
        self.strInputFile = os.path.join(
            self.strSampleDir, self.strFastq_name + ".fastq"
        )
        ## './Input/JaeWoo/FASTQ/JaeWoo_test1_samples/Test_sample/Split_files'
        self.strSplitPath = os.path.join(self.strSampleDir, "Split_files")
        Helper.MakeFolderIfNot(self.strSplitPath)

    def _split_into_chunks(self):
        import shlex
        import pathlib

        sp.run(
            shlex.split(
                shlex.quote(
                    f"split {self.strInputFile} -l {4 * self.intChunkSize} -d -a 6 --additional-suffix=.fq {self.strSplitPath}/split_"
                )
            ),
            shell=True,
            check=True,
        )
        logging.info(
            f"The number of split files: {len(list(pathlib.Path(self.strSplitPath).glob('*')))}"
        )

    def MakeReference(self):
        with open(self.strRefFile, "w") as output:
            import pandas as pd

            refs = pd.read_csv(f"{self.strRefDir}/reference.csv")

            ## defensive
            assert refs.isna().sum().sum() == 0, "reference.csv has some blanks."

            listBarcode = Helper.preprocess_ref_file(refs["Barcode"])
            listTarget = Helper.preprocess_ref_file(refs["Target_region"])
            listRef = Helper.preprocess_ref_file(refs["Reference_sequence"])

            # String pre-processing

            for strBar, strTar, strRef in zip(listBarcode, listTarget, listRef):
                output.write(">" + f'{strBar + ":" + strTar}' + "\n" + strRef + "\n")

        # The strings should be composed of only A,T,C,G,N
        Helper.CheckIntegrity(self.strRefFile)  ## defensive

    def RunIndelFreqCalculator(self):
        sp.call(
            "{python} Indel_frequency_calculator.py {outdir} {sample} {logpath}".format(
                python=self.strPython,
                outdir=self.strOutSampleDir,
                sample=self.strSample,
                logpath=self.strLogPath,
            ),
            shell=True,
        )
        sp.call(
            "{python} Summary_all_trim.py {outdir} {sample} {logpath}".format(
                python=self.strPython,
                outdir=self.strOutSampleDir,
                sample=self.strSample,
                logpath=self.strLogPath,
            ),
            shell=True,
        )
        sp.call(
            'cp $(find ./Output/{user}/{project} -name "*.tsv") ./Output/{user}/{project}/All_results'.format(
                user=self.strUser, project=self.strProject
            ),
            shell=True,
        )

    def IndelNormalization(self):

        sp.call(
            "{python} Indel_normalization.py {project_file} {user} {project}".format(
                python=self.strPython,
                project_file=self.strProjectFile,
                user=self.strUser,
                project=self.strProject,
            ),
            shell=True,
        )

    def MakeOutput(self):
        """
        dictResult
        {'TTTGTAGTCATACATCGCAATGTCAA': [0, 0, 0, 0, [], [], [], [], []]}
        dictResultIndelFreq
        {'TTTGCTCAGTCACACGTCACGAGCTG': [['TCATCGACTTGCAGGACATTAGGCGAAGCTTGGCGTAACTAGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTC',
        ['TCATCGACTTGCAGGACGAAGCTTGGCGAAGCTTGGCGTAACTAGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAATA'], '19M3I', 1.0,
        'TCATCGACTTGCAGGACATTAGGCGA', ['TCATCGACTTGCAGGACAT---TAGGCGAAGCTTGGCGTAACTAGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTC---------'], ['TCATCGACTTGCAGGACGAAGCTTGGCGAAGCTTGGCGTAACTAGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAATA']]]}
        strBarcodePamPos
        Foward
        """
        # index name, constant variable.
        intTotal = 0
        intNumIns = 1
        intNumDel = 2
        intNumCom = 3
        intTotalFastq = 4
        intInsFastq = 5
        intDelFastq = 6
        intComFastq = 7
        intIndelInfo = 8

        with open(
            "{outdir}/Tmp/{sample}_Summary.txt".format(
                outdir=self.strOutSampleDir, sample=self.strSample
            ),
            "w",
        ) as Summary, open(
            "{outdir}/Tmp/{sample}_Classified_Indel_barcode.fastq".format(
                outdir=self.strOutSampleDir, sample=self.strSample
            ),
            "w",
        ) as FastqOut, open(
            "{outdir}/Tmp/{sample}_Indel_freq.txt".format(
                outdir=self.strOutSampleDir, sample=self.strSample
            ),
            "w",
        ) as FreqOut:

            for binPickle in os.listdir(
                "{outdir}/Tmp/Pickle".format(outdir=self.strOutSampleDir)
            ):
                with open(
                    "{outdir}/Tmp/Pickle/{pickle}".format(
                        outdir=self.strOutSampleDir, pickle=binPickle
                    ),
                    "rb",
                ) as PickleResult:

                    dictPickleResult = pickle.load(PickleResult)
                    dictResult = dictPickleResult["dictResult"]
                    dictResultIndelFreq = dictPickleResult["dictResultIndelFreq"]
                    strBarcodePamPos = dictPickleResult["strBarcodePamPos"]

                    for strBarcode, listValue in dictResult.items():
                        if strBarcodePamPos == "Reverse":
                            strBarcode = strBarcode[::-1]

                        Summary.write(
                            "{Bar}\t{NumTot}\t{NumIns}\t{NumDel}\t{NumCom}\n".format(
                                Bar=strBarcode,
                                NumTot=listValue[intTotal],
                                NumIns=listValue[intNumIns],
                                NumDel=listValue[intNumDel],
                                NumCom=listValue[intNumCom],
                            )
                        )

                        if self.strClassFASTQ == "True":
                            for strJudge, intFastqKind in [
                                ("total", intTotalFastq),
                                ("insertion", intInsFastq),
                                ("deletion", intDelFastq),
                                ("complex", intComFastq),
                            ]:
                                for listFastq in listValue[intFastqKind]:  ## category
                                    listFastqAddClass = [
                                        listFastq[0]
                                        + ":Barcode_%s:%s" % (strBarcode, strJudge)
                                    ]
                                    FastqOut.write(
                                        "\n".join(listFastqAddClass + listFastq[1:])
                                        + "\n"
                                    )

                    for (
                        strBarcode
                    ) in (
                        dictResultIndelFreq
                    ):  # dictResultIndelFreq [sRef_seq, lQuery, float(iFreq)/iTotal, sTarget_region]

                        if strBarcodePamPos == "Reverse":
                            strBarcode = strBarcode[::-1]

                        for (
                            strRefSeq,
                            listQuery,
                            strINDEL,
                            floFreq,
                            strTargetRegion,
                            listRefNeedle,
                            listQueryNeedle,
                        ) in sorted(
                            dictResultIndelFreq[strBarcode],
                            key=lambda x: x[3],
                            reverse=True,
                        ):
                            for strQuery, strRefNeedle, strQueryNeedle in zip(
                                listQuery, listRefNeedle, listQueryNeedle
                            ):

                                if strBarcodePamPos == "Reverse":
                                    strQuery = strQuery[::-1]
                                    strRefNeedle = strRefNeedle[::-1]
                                    strQueryNeedle = strQueryNeedle[::-1]

                                FreqOut.write(
                                    "\t".join(
                                        [
                                            strBarcode,
                                            strQuery,
                                            strINDEL,
                                            str(round(floFreq, 4)),
                                            strRefNeedle,
                                            strQueryNeedle,
                                        ]
                                    )
                                    + "\n"
                                )

            if self.strPickle == "False":
                logging.info("Delete tmp pickles")
                sp.call(
                    "rm {outdir}/Tmp/Pickle/*.pickle".format(
                        outdir=self.strOutSampleDir
                    ),
                    shell=True,
                )

            elif self.strSplit == "False":
                logging.info("Delete splited input files")
                sp.call(
                    "rm {split_path}/*.fq".format(split_path=self.strSplitPath),
                    shell=True,
                )


@CheckProcessedFiles
def RunPipeline(**kwargs):
    InstInitFolder = kwargs["InstInitFolder"]
    listSamples = kwargs["listSamples"]
    args = kwargs["args"]
    logging = kwargs["logger"]

    setGroup = set()

    for tupSampleInfo in listSamples:
        # tupSampleInfo == (strSample, strRef, strExpCtrl)

        strSample, strRef, strExpCtrl = tupSampleInfo
        setGroup.add(strExpCtrl)

        # Locating input files in the Reference and FASTQ directory
        InstRunner = clsIndelSearcherRunner(strSample, strRef, args, InstInitFolder)

        # Chunking
        logging.info("SplitFile")
        InstRunner._split_into_chunks()

        # If there is no manually created "Reference.fa", make it for future processing
        # Need proper "Barcode.txt", "Reference_sequence.txt", "Target_region.txt"
        logging.info("MakeReference")
        InstRunner.MakeReference()

        logging.info("RunMulticore")
        RunMulticore(InstRunner, args.multicore)  ## from CoreSystem.py
        logging.info("MakeOutput")
        InstRunner.MakeOutput()
        logging.info("RunIndelFreqCalculator")
        InstRunner.RunIndelFreqCalculator()

    if setGroup == {"EXP", "CTRL"}:
        InstRunner.IndelNormalization()
    elif setGroup in [set(), set([]), set([""]), set([" "])]:
        pass
    else:
        logging.error("The group category is not appropriate. : %s" % setGroup)
        logging.error("Please make sure your project file is correct.")
        logging.error("The group category must be Exp or Ctrl")
        raise Exception


def indel_searcher_runner(args):
    """
    This function is the main function that runs the pipeline

    :param args: This is the arguments that are passed to the program
    """
    import datetime

    InstInitFolder = InitialFolder(
        args.user_name, args.project_name, os.path.basename(__file__)
    )
    InstInitFolder.MakeDefaultFolder()

    logging.basicConfig(
        format="%(process)d %(levelname)s %(asctime)s : %(message)s",
        level=logging.DEBUG,
        filename=InstInitFolder.strLogPath,
        filemode="a",
    )
    logger = logging.getLogger()
    logger.addHandler(logging.StreamHandler(sys.stdout))

    logger.info(f"Program start at {datetime.datetime.now()}")
    logger.info(str(vars(args)))

    with open(InstInitFolder.strProjectFile, "r") as f:
        Sample_list = f.readlines()
        listSamples = Helper.preprocess_user_file(Sample_list)

        strInputProject = "./Input/{user}/FASTQ/{project}".format(
            user=args.user_name, project=args.project_name
        )

        RunPipeline(
            InstInitFolder=InstInitFolder,
            strInputProject=strInputProject,
            listSamples=listSamples,
            args=args,
            logger=logger,
        )

    logger.info("Program end")
