#! /usr/bin/env python3

import sys
import os
import subprocess
import shutil
import argparse
import datetime
import string
import sys
import time
import traceback
import types
import tempfile
import shlex
import xml.etree.ElementTree as ET
from Bio import SeqIO
import psutil
import math
import json
import threading

import workflow

from progress import Progress

class NonZeroReturnValueException(Exception):
    def __init__(self, returnvalue, msg):
        self.msg = msg
        self.returnvalue = returnvalue
        return

def logline(log_fh, text):
    log_fh.write(text+"\n")
    log_fh.flush()
    return

big_lock = threading.Lock()
notification_lock = threading.Lock()
notifications = []

if __name__=="__main__":

    parser = argparse.ArgumentParser(description='glaDIAtor')

    parser.add_argument('--project-name',
                        action='store',
                        dest='project_name',
                        required=True,
                        default=None,
                        help='Project name. A result folder with given name will be created.')

    parser.add_argument('--library-data',
                        nargs='+',
                        action='store',
                        dest='library_DDA',
                        required=False,
                        default=None,
                        help='MSMS files for DDA library in raw or mzXML')

    parser.add_argument('--sample-data',
                        nargs='+',
                        action='store',
                        dest='sample_DIA',
                        required=True,
                        default=None,
                        help='MSMS DIA sample files in raw or mzML format')

    parser.add_argument('--databases',
                        nargs='+',
                        action='store',
                        dest='databases',
                        required=True,
                        default=None,
                        help='Fasta files of the database of peptide search space')

    parser.add_argument('--comet-cfg-template', 
                        action='store',
                        dest='comet_cfg_template',
                        default="/opt/gladiator/comet.params.template",
                        required=False,
                        help='Comet config template file.')

    parser.add_argument('--xtandem-cfg-template', 
                        action='store',
                        dest='xtandem_cfg_template',
                        default="/opt/gladiator/xtandem_settings.xml",
                        required=False,
                        help='XTandem config template file.')

    parser.add_argument('--library-FDR', 
                        action='store',
                        dest='library_FDR',
                        required=False,
                        default="0.01",
                        help='Set FDR used in spectral library build. [default: 0.01]')

    parser.add_argument('--feature-alignment-FDR', 
                        nargs='+',
                        action='store',
                        dest='feature_alignment_FDR',
                        required=False,
                        default=["0.01", "0.05"],
                        help='Set target and max FDR used in TRIC alignment in diRT mode. [default: [0.01, 0.05]]')

    parser.add_argument('--threads', 
                        action='store',
                        dest='threads',
                        required=False,
                        default=None,
                        help='Set amount of threads. [default: auto]')

    parser.add_argument('--design-file', 
                        action='store',
                        dest='design_file',
                        required=False,
                        default=None,
                        help='Design file to be used in the differential expression analysis.')

    parser.add_argument('--retain-tmp-files', 
                        action='store_true',
                        dest='retain_tmp_files',
                        required=False,
                        default=False,
                        help='Dont delete temp files. Used for debug purposes.')

    parser.add_argument('--search-engines',
                        nargs='+',
                        action='store',
                        dest='search_engines',
                        required=False,
                        default=["comet", "xtandem"],
                        help='Specify search engines to be used. Default: comet and xtandem')

    args = parser.parse_args()

    cwd = os.getcwd()

    delete_temp_files_flag = not args.retain_tmp_files
    
    max_threads = os.cpu_count()
    if not max_threads:
        max_threads = 1

    if args.threads:
        max_threads = args.threads

    result_root = "/run-files"

    project = args.project_name

    cwd = os.path.join(result_root, project)

    if not os.path.exists(cwd):
        os.mkdir(cwd)

    analysis_name = project
    sample_files = args.sample_DIA
    library_files = []
    if args.library_DDA:
        library_files = args.library_DDA
  
    database_files = args.databases
    pvalue = args.library_FDR
    trig_target_pvalue = args.feature_alignment_FDR[0]
    trig_max_pvalue = args.feature_alignment_FDR[1]

    options = []

    #if not args.library_DDA:
    #    options.append("use_speudospectra_flag")

    if "comet" in args.search_engines:
        options.append("use_comet_flag")
        
    if "xtandem" in args.search_engines:
        options.append("use_xtandem_flag")
        
    delete_tmp_files_flag = not args.retain_tmp_files

    cfgfile = os.path.join(cwd, "config.txt")

    scan = None
    data = None

    if not os.path.isfile(cfgfile):

        cfg = {
            "analysis_name": analysis_name, 
            "files": {
                "samples": sample_files, 
                "library": library_files, 
                "database": database_files
            }, 
            "pvalue": pvalue,
            "trig_target_pvalue": trig_target_pvalue, 
            "trig_max_pvalue": trig_max_pvalue, 
            "options": options
        }

        with open(cfgfile, "w") as fh:
            json.dump(cfg, fh)

    else:
        scan = workflow.scan_project_phases(project, result_root)
        if not "config" in scan:
            print ("Config not found from existing project.")
            sys.exit(1)
        data = scan["config"]

    # Sanity check here

    # 1) At least one DIA file
    
    analysis_state = workflow.State(project)
    event  = workflow.Event()

    if scan: # TODO
        log_name = "rerun-log.txt"
    else:
        log_name = "log.txt"

    with open(os.path.join(cwd, log_name), "w") as log_fh:

        logline(log_fh, "Pipeline command line: " + " ".join(sys.argv))
        logline(log_fh, "Option list: " + "; ".join(options))
        
        # Convert DIA data to open format

        extensions = set()
        for filename in sample_files:
            extension = os.path.splitext(filename)[1].lower()
            extensions.add(extension)

        extensions = list(extensions)

        # TODO: Return the actual error
        if len(extensions) > 1:
            print("A single input type allowed, multile file types found (" + ", ".join(extensions) + ")")
            sys.exit(1)

        # TODO: Return the actual error
        if extension not in [".mzml", ".mzxml", ".raw"]:
            print("Unknown extension (" + extension + ") found.")
            sys.exit(1)

        if (extension == ".raw"):

            if not os.path.exists("/wineprefix64"):
                print("This version of gladiator does not support raw files.")
                sys.exit(1)

            if scan and "conv_DIA" in scan and scan["conv_DIA"]:
                print ("Detected existing converted DIA files.")
                None

            else:
                phaseID = analysis_state.createPhase("Converting RAW DIA files")
                workflow.convertRAW(\
                    event, \
                    log_fh,\
                    analysis_state.getPhaseData(phaseID), \
                    cwd, \
                    sample_files, \
                    True, \
                    "mzml",
                    "DIA")

            if event.kill_flag:
                sys.exit(1)

            # replace sample files with converted versions
            sample_files = [os.path.join(cwd, "DIA", x) for x in os.listdir(os.path.join(cwd, "DIA"))]


        # Convert DDA data to open format

        if library_files:

            extensions = set()
            extension = None

            for filename in library_files:
                extension = os.path.splitext(filename)[1].lower()
                extensions.add(extension)

            extensions = list(extensions)

            # TODO: Return the actual error
            if len(extensions) > 1:
                print("A single input type allowed, multile file types found (" + ", ".join(extensions) + ")")
                sys.exit(1)

            # TODO: Return the actual error
            if extension not in [".mzml", ".mzxml", ".raw"]:
                print("Unknown extension (" + extension + ") found.")
                sys.exit(1)

            if (extension == ".raw"):

                if not os.path.exists("/wineprefix64"):
                    print("This version of gladiator does not support raw files.")
                    sys.exit(1)

                if scan and "conv_DDA" in scan and scan["conv_DDA"]:
                    print ("Detected existing converted DDA files.")
                    None

                else:
                    phaseID = analysis_state.createPhase("Converting RAW DDA files (+picking peaks)")
                    workflow.convertRAWqtofpeakpicker(\
                        event, \
                        log_fh, \
                        analysis_state.getPhaseData(phaseID), \
                        cwd, \
                        library_files, \
                        "DDA")

                if event.kill_flag:
                    sys.exit(1)
                            
                # replace library raw files with converted versions
                library_files = [os.path.join(cwd, "DDA", x) for x in os.listdir(os.path.join(cwd, "DDA"))]


        swaths_min = 0
        swaths_max = 0

        swaths, tswaths = workflow.create_swath_window_files(cwd, sample_files[0])
        swaths_min = swaths[0][0]
        swaths_max = swaths[-1][1]

        pseudospectrafiles = []
        if not args.library_DDA:

            if scan and "pseudospectra" in scan and scan["pseudospectra"]:
                print ("Detected existing converted DDA files.")
                pseudospectrafiles = scan["pseudospectrafiles"]
                None

            else:
                phaseID = analysis_state.createPhase("Building pseudospectra")
                pseudospectrafiles = workflow.runDiaumpire(\
                    event, \
                    log_fh, \
                    analysis_state.getPhaseData(phaseID), \
                    cwd, \
                    sample_files, \
                    "libfree",
                    max_threads)

            if event.kill_flag:
                sys.exit(1)

        # BUILD SEQUENCE DATABASE

        db_filename = "DB.fasta"
        decoy_db_file = "DB_with_decoys.fasta"

        if scan and "DB" in scan and scan["DB"]:
            print ("Detected existing sequence database")
        else:
            phaseID = analysis_state.createPhase("Building database")
            workflow.build_database(\
                event, \
                log_fh, \
                analysis_state.getPhaseData(phaseID), \
                cwd, \
                database_files, \
                db_filename, \
                decoy_db_file)

            if event.kill_flag:
                sys.exit(1)


        # BUILD lib & pseudolib peptide files

        libmethod_pepXMLs = []
        libfreemethod_pepXMLs = []

        if "use_comet_flag" in options:

            comet_cfg = "/opt/gladiator/comet.params.template"

            if library_files:

                if scan and "comet_peptides" in scan and scan["comet_peptides"]:
                    print ("Detected existing Comet spectra search results")

                else:

                    comet_cfg = os.path.join(cwd, "comet_settings.xml")

                    cfg_txt = None
                    with open("/opt/gladiator/comet_settings_template.xml", "r") as fh:
                        cfg_txt = fh.read()

                    param_map = {
                        "PRECURSOR_MASS_TOLERANCE":"20",
                        "FRAGMENT_MASS_TOLERANCE": "{:.2f}".format(float("0.02")/2),
                        "DATABASE_FASTA_FILE": decoy_db_file
                    }    

                    for k in param_map:
                        value = param_map[k]
                        cfg_txt = cfg_txt.replace(k, value)

                    with open(comet_cfg,"w") as fh:
                        fh.write(cfg_txt)


                    # Cleanup

                    filenames = \
                    ["interact_comet_pep-MODELS.html", \
                    "interact_comet_pep.xml", \
                    "interact_comet_pep.xml.index", \
                    "interact_comet_pep.xml.RTcoeff", \
                    "interact_comet_pep.xml.RTstats"]

                    for filename in filenames:
                        filepath = os.path.join(cwd, filename)
                        if os.path.isfile(filepath):
                            os.remove(filepath)

                    # Run Comet

                    existing_pep_xmls = None
                    if scan and "samples_generated_by_comet" in scan:
                        existing_pep_xmls = scan["samples_generated_by_comet"]

                    phaseID = analysis_state.createPhase("Speclib - Matching sequences [Comet]")
                    workflow.runComet(\
                        event, \
                        log_fh, \
                        analysis_state.getPhaseData(phaseID), \
                        cwd, \
                        comet_cfg, \
                        library_files, \
                        existing_pep_xmls, \
                        "interact_comet_pep.xml", \
                        cwd)

                if event.kill_flag:
                    sys.exit(1)

                libmethod_pepXMLs.append("interact_comet_pep.xml")

            if pseudospectrafiles:

                if scan and "comet_pseudo_peptides" in scan and scan["comet_pseudo_peptides"]:
                    print ("Detected existing Comet pseudo spectra search results")

                else:

                    # Cleanup

                    filenames = \
                    ["interact_comet_pseudo_pep-MODELS.html", \
                    "interact_comet_pseudo_pep.xml", \
                    "interact_comet_pseudo_pep.xml.index", \
                    "interact_comet_pseudo_pep.xml.RTcoeff", \
                    "interact_comet_pseudo_pep.xml.RTstats"]

                    for filename in filenames:
                        filepath = os.path.join(cwd, filename)
                        if os.path.isfile(filepath):
                            os.remove(filepath)

                    existing_pep_xmls = None
                    if scan and "pseudo_samples_generated_by_comet" in scan:
                        existing_pep_xmls = scan["pseudo_samples_generated_by_comet"]

                    # Run Comet for speudospectrafiles

                    #decoy_db_file

                    phaseID = analysis_state.createPhase("Pseudospeclib - Matching sequences [Comet]")
                    workflow.runComet(\
                        event, \
                        log_fh, \
                        analysis_state.getPhaseData(phaseID), \
                        cwd, \
                        comet_cfg, \
                        pseudospectrafiles, \
                        existing_pep_xmls, \
                        "interact_comet_pseudo_pep.xml", \
                        cwd)

                if event.kill_flag:
                    sys.exit(1)

                libfreemethod_pepXMLs.append("interact_comet_pseudo_pep.xml")            

        if "use_xtandem_flag" in options:
            xtandem_cfg = "/opt/gladiator/xtandem_settings.xml"
            if library_files:

                if scan and "xtandem_peptides" in scan and scan["xtandem_peptides"]:
                    print ("Detected existing X!Tandem spectra search results")
                else:

                    xtandem_cfg = os.path.join(cwd, "xtandem_settings.xml")

                    cfg_txt = None
                    with open("/opt/gladiator/xtandem_settings_template.xml", "r") as fh:
                        cfg_txt = fh.read()

                    param_map = {
                        "PRECURSOR_MASS_TOLERANCE":"20",
                        "FRAGMENT_MASS_TOLERANCE": "{:.2f}".format(float("0.02")/2),
                        #"DATABASE_FASTA_FILE": decoy_db_file
                    }    

                    for k in param_map:
                        value = param_map[k]
                        cfg_txt = cfg_txt.replace(k, value)

                    with open(xtandem_cfg,"w") as fh:
                        fh.write(cfg_txt)


                    # Cleanup

                    filenames = [\
                    "interact_xtandem_pep-MODELS.html   ", \
                    "interact_xtandem_pep.xml", \
                    "interact_xtandem_pep.xml.index", \
                    "interact_xtandem_pep.xml.RTcoeff   ", \
                    "interact_xtandem_pep.xml.RTstats   "]

                    for filename in filenames:
                        filepath = os.path.join(cwd, filename)
                        if os.path.isfile(filepath):
                            os.remove(filepath)

                    # Run X!Tandem

                    existing_pep_xmls = None
                    if scan and "samples_generated_by_xtandem" in scan:
                        existing_pep_xmls = scan["samples_generated_by_xtandem"]

                    phaseID = analysis_state.createPhase("Speclib - Matching sequences [X!Tandem]")
                    workflow.runXTandem(\
                        event, \
                        log_fh, \
                        analysis_state.getPhaseData(phaseID), \
                        cwd, \
                        xtandem_cfg, \
                        decoy_db_file, \
                        library_files, \
                        existing_pep_xmls, \
                        "interact_xtandem_pep.xml", \
                        cwd, \
                        delete_tmp_files_flag)

                if event.kill_flag:
                    sys.exit(1)

                libmethod_pepXMLs.append("interact_xtandem_pep.xml")

            if pseudospectrafiles:

                if scan and "xtandem_pseudo_peptides" in scan and scan["xtandem_pseudo_peptides"]:
                    print ("Detected existing X!Tandem pseudo spectra search results")

                else:

                    # Cleanup

                    filenames = [\
                    "interact_xtandem_pseudo_pep-MODELS.html", \
                    "interact_xtandem_pseudo_pep.xml", \
                    "interact_xtandem_pseudo_pep.xml.index", \
                    "interact_xtandem_pseudo_pep.xml.RTcoeff", \
                    "interact_xtandem_pseudo_pep.xml.RTstats"]

                    for filename in filenames:
                        filepath = os.path.join(cwd, filename)
                        if os.path.isfile(filepath):
                            os.remove(filepath)

                    # Run X!Tandem for speudospectrafiles

                    existing_pep_xmls = None
                    if scan and "pseudo_samples_generated_by_xtandem" in scan:
                        existing_pep_xmls = scan["pseudo_samples_generated_by_xtandem"]

                    phaseID = analysis_state.createPhase("Pseudospeclib - Matching sequences [X!Tandem]")
                    workflow.runXTandem(\
                        event, \
                        log_fh, \
                        analysis_state.getPhaseData(phaseID), \
                        cwd, \
                        xtandem_cfg, \
                        decoy_db_file, \
                        pseudospectrafiles, \
                        existing_pep_xmls, \
                        "interact_xtandem_pseudo_pep.xml",
                        cwd, \
                        delete_tmp_files_flag)

                if event.kill_flag:
                    sys.exit(1)

                libfreemethod_pepXMLs.append("interact_xtandem_pseudo_pep.xml")         

        if scan and "speclib_cons" in scan and scan["speclib_cons"]:

            print ("Detected existing Spectrum library.")

        else:

            phaseID = analysis_state.createPhase("Building Library")

            # Cleanup

            filenames = \
            ["SpecLib_cons_decoy.TraML",\
            "SpecLib_cons_openswath.tsv",\
            "SpecLib_cons.pepidx",\
            "SpecLib_cons.spidx",\
            "SpecLib_cons.splib",\
            "SpecLib_cons.sptxt",\
            "SpecLib_cons.TraML",\
            "SpecLib_libfree.pepidx",\
            "SpecLib_libfree.spidx",\
            "SpecLib_libfree.splib",\
            "SpecLib_libfree.sptxt",\
            "SpecLib_lib.pepidx",\
            "SpecLib_lib.spidx",\
            "SpecLib_lib.splib",\
            "SpecLib_lib.sptxt",\
            "SpecLib_merged.pepidx",\
            "SpecLib_merged.spidx",\
            "SpecLib_merged.splib",\
            "SpecLib_merged.sptxt",\
            "spectrast.log"]

            for filename in filenames:
                filepath = os.path.join(cwd, filename)
                if os.path.isfile(filepath):
                    os.remove(filepath)

            for method in ["lib", "libfree"]:
                mayu_folder = method+"_mayu"
                directory = os.path.join(cwd, mayu_folder)
                # print ("DEBUG: CHECKING TO REMOVE DIR " + directory)
                if os.path.isdir(directory):
                    shutil.rmtree(directory)
                    # print ("DEBUG: REMOVED " + directory)

            workflow.buildlib(\
                event,\
                log_fh,\
                analysis_state.getPhaseData(phaseID),\
                cwd, \
                "DECOY_", \
                decoy_db_file, \
                str(max_threads), \
                libmethod_pepXMLs, \
                libfreemethod_pepXMLs, \
                "/opt/gladiator/iRT.txt", \
                "swath-windows.txt", \
                pvalue, \
                pvalue, \
                swaths_min, \
                swaths_max)

        if event.kill_flag:
            sys.exit(1)

        if scan and "matrices" in scan and scan["matrices"]:

            print ("Detected existing matrices.")
        
        else:

            # Cleanup

            postfixes = ["-DIA_cutoffs.txt",\
            "-DIA_dscores_top_decoy_peaks.txt",\
            "-DIA_dscores_top_target_peaks.txt",\
            "-DIA_full_stat.csv",\
            "-DIA_mayu.csv",\
            "-DIA_mayu.cutoff",\
            "-DIA_mayu.fasta",\
            "-DIA_qvalues.txt",\
            "-DIA_report.pdf",\
            "-DIA_scorer.bin",\
            "-DIA_summary_stat.csv",\
            "-DIA_svalues.txt",\
            "-DIA.tsv",\
            "-DIA_weights.txt",\
            "-DIA_with_dscore-0_0-None.tr",\
            "-DIA_with_dscore.csv",\
            "-DIA_with_dscore_filtered.csv"]

            for filename in sample_files:
                filepath = os.path.join(cwd, os.path.basename(filename))
                for postfix in postfixes:
                    if filepath.endswith(postfix) and os.path.isfile(filepath):
                        os.remove(filepath)

            filenames = \
                ["DIA-analysis-result.csv",\
                "DIA-peptide-matrix.tsv",\
                "DIA-protein-matrix.tsv"]

            for filename in filenames:
                filepath = os.path.join(cwd, filename)
                if os.path.isfile(filepath):
                    os.remove(filepath)

            # Run OpenSWATH and building the DIA matrices
            phaseID = analysis_state.createPhase("Searching peptides from DIA spectrum")

            workflow.buildDIAMatrix(\
                event, \
                log_fh, \
                analysis_state.getPhaseData(phaseID), \
                cwd, \
                sample_files, \
                "truncated-swath-windows.txt", \
                trig_target_pvalue, \
                trig_max_pvalue, \
                str(max_threads), \
                "/opt/gladiator/iRTAssayLibrary.TraML",\
                None) # TODO: insert design file here

        if event.kill_flag:
            sys.exit(1)

        with notification_lock:
            notifications.append(
                {
                    "title": "Analysis completed",
                    "text": "Analysis of project " + project + " was completed successfully.", 
                }
            )









    sys.exit(0)



