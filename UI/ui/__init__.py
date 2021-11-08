from pyramid.config import Configurator
from pyramid.response import Response
from pyramid.view import view_config
from datetime import datetime
import threading
import workflow
import json
import os
import shutil
from pathlib import Path
import pandas as pd
import time
import plotly.offline as py
import tempfile

import plotly.graph_objects as go
import numpy as np

import annotation
import annotateSwath2stats

result_root = "/run-files/"

big_lock = threading.Lock()

analysis_threads = {}
annotation_threads = {}


class Ticker(threading.Thread):

    def __init__(self):
        threading.Thread.__init__(self)
        self.active = True

    def run(self):

        # Clean progress state if analysis thread has terminated
        while self.active:
            with big_lock:
                for project in state:
                    if "analysis" in state[project]:
                        if project in analysis_threads and analysis_threads[project].isAlive():
                            None
                        else:
                            del state[project]["analysis"]

            time.sleep(5)


ticker = Ticker()

notification_lock = threading.Lock()
notifications = []

analysis_events = {}


class AnalysisThread(threading.Thread):

    def __init__(self, data, folder):
        threading.Thread.__init__(self)
        self.data = data
        self.project_folder = folder
        self.scan = None

    def get_project_folder(self):
        return self.project_folder

    def get_raw_converter(self, event, log_fh, analysis_state):
        if os.path.exists("/wineprefix64"):
            return ("msconvert")
        if os.path.exists("/opt/ThermoRawFileParser/ThermoRawFileParser.exe"):
            return ("thermorawparser")

        phaseID = analysis_state.createPhase("Installing ThermoRawFileParser")
        workflow.install_ThermoRawFileParser(event, log_fh, analysis_state.getPhaseData(phaseID))
        if os.path.exists("/opt/ThermoRawFileParser/ThermoRawFileParser.exe"):
            return ("thermorawparser")
        else:
            return None
    
    def run(self):

        max_threads = os.cpu_count()
        if not max_threads:
            max_threads = 1

        project = self.project_folder

        cwd = os.path.join(result_root, project)

        if self.data:
            with open(os.path.join(cwd, "config.txt"), "w") as fh:
                json.dump(self.data, fh)
        else:
            self.scan = workflow.scan_project_phases(project, result_root)
            if not "config" in self.scan:
                print("Config not found from existing project.")
                return
            self.data = self.scan["config"]

        analysis_name = self.data['analysis_name']
        sample_files = self.data['files']['samples']
        library_files = self.data['files']['library']
        database_files = self.data['files']['database']
        pvalue = self.data['pvalue']
        trig_target_pvalue = self.data['trig_target_pvalue']
        trig_max_pvalue = self.data['trig_max_pvalue']
        options = self.data['options']

        if 'force_threads' in options:
            max_threads = self.data['threads']
            print ("Forcing max " + str(max_threads) + " threads.")

        delete_tmp_files_flag = True

        # Sanity check here

        # 1) At least one DIA file

        global state
        event = None
        with big_lock:
            analysis_state = state[project]["analysis"]
            event = analysis_events[project]

        if self.scan:  # TODO
            log_name = "rerun-log.txt"
        else:
            log_name = "log.txt"

        with open(os.path.join(cwd, log_name), "w") as log_fh:

            # Convert DIA data to open format

            extensions = set()
            for filename in sample_files:
                extension = os.path.splitext(filename)[1].lower()
                extensions.add(extension)

            extensions = list(extensions)

            # TODO: Return the actual error
            if len(extensions) > 1:
                print("A single input type allowed, multile file types found (" +
                      ", ".join(extensions) + ")")
                return

            # TODO: Return the actual error
            if extension not in [".mzml", ".mzxml", ".raw"]:
                print("Unknown extension (" + extension + ") found.")
                return

            if (extension == ".raw"):
                converter = self.get_raw_converter(event, log_fh, analysis_state)

                if self.scan and "conv_DIA" in self.scan and self.scan["conv_DIA"]:
                    print("Detected existing converted DIA files.")
                    None

                else:

                    phaseID = analysis_state.createPhase("Converting RAW DIA files")

                    if converter == "msconvert":
                        workflow.convertRAW(
                            event,
                            log_fh,
                            analysis_state.getPhaseData(phaseID),
                            cwd,
                            sample_files,
                            False,
                            "mzml",
                            "DIA")

                    elif converter == "thermorawparser":
                        workflow.convertRAW_ThermoRawFileParser(
                            event,
                            log_fh,
                            analysis_state.getPhaseData(phaseID),
                            cwd,
                            sample_files,
                            False,
                            "mzml",
                            "DIA")
                    else:
                        print ("Error: No, RAW file converter available")
                        return
                        

                if event.kill_flag:
                    return

                # replace sample files with converted versions
                sample_files = [os.path.join(cwd, "DIA", x)
                                for x in os.listdir(os.path.join(cwd, "DIA")) if Path(x).suffix in [".mzXML", ".mzML"] ]

            swaths_min = 0
            swaths_max = 0

            swaths, tswaths = workflow.create_swath_window_files(
                cwd, sample_files[0])
            swaths_min = swaths[0][0]
            swaths_max = swaths[-1][1]

            pseudospectrafiles = []

            if library_files and ("dda_library" in options):

                extensions = set()
                extension = None

                for filename in library_files:
                    extension = os.path.splitext(filename)[1].lower()
                    extensions.add(extension)

                extensions = list(extensions)

                # TODO: Return the actual error
                if len(extensions) > 1:
                    print("A single input type allowed, multile file types found (" + ", ".join(extensions) + ")")
                    return

                # TODO: Return the actual error
                if extension not in [".mzml", ".mzxml", ".raw"]:
                    print("Unknown extension (" + extension + ") found.")
                    return

                if (extension == ".raw"):

                    converter = self.get_raw_converter(event, log_fh, analysis_state)

                    if self.scan and "conv_DDA" in self.scan and self.scan["conv_DDA"]:
                        print("Detected existing converted DDA files.")
                        None

                    else:
                        phaseID = analysis_state.createPhase(
                            "Converting RAW DDA files (+picking peaks)")

                        if converter == "msconvert":
                        
                            workflow.convertRAWqtofpeakpicker(
                                event,
                                log_fh,
                                analysis_state.getPhaseData(phaseID),
                                cwd,
                                library_files,
                                "DDA")

                        elif converter == "thermorawparser":

                            workflow.convertRAW_ThermoRawFileParser(
                                event,
                                log_fh,
                                analysis_state.getPhaseData(phaseID),
                                cwd,
                                library_files,
                                True,
                                "mzml",
                                "DDA")

                        else:
                            print ("Error: No, RAW file converter available")
                            return


                    if event.kill_flag:
                        return

                    # replace library raw files with converted versions
                    library_files = [os.path.join(
                        cwd, "DDA", x) for x in os.listdir(os.path.join(cwd, "DDA")) if Path(x).suffix in [".mzXML", ".mzML"]]


            else : # BUILD pseudospectral library instead of spectral library

                if self.scan and "pseudospectra" in self.scan and self.scan["pseudospectra"]:
                    # print("Detected existing converted DDA files.")
                    pseudospectrafiles = self.scan["pseudospectrafiles"]
                    None

                else:
                    phaseID = analysis_state.createPhase(
                        "Building pseudospectra")
                    pseudospectrafiles = workflow.runDiaumpire(
                        event,
                        log_fh,
                        analysis_state.getPhaseData(phaseID),
                        cwd,
                        sample_files,
                        "libfree",
                        max_threads)

                if event.kill_flag:
                    return

            # BUILD SEQUENCE DATABASE

            db_filename = "DB.fasta"
            decoy_db_file = "DB_with_decoys.fasta"

            if self.scan and "DB" in self.scan and self.scan["DB"]:
                None
                #print("Detected existing sequence database")
            else:
                phaseID = analysis_state.createPhase("Building database")
                workflow.build_database(
                    event,
                    log_fh,
                    analysis_state.getPhaseData(phaseID),
                    cwd,
                    database_files,
                    db_filename,
                    decoy_db_file)

                if event.kill_flag:
                    return

            # BUILD lib & pseudolib peptide files

            libmethod_pepXMLs = []
            libfreemethod_pepXMLs = []

            param_map = {
                "precursor_tolerance" : "PRECURSOR_MASS_TOLERANCE",
                "fragment_tolerance": "FRAGMENT_MASS_TOLERANCE"
            }    

            if "use_comet_flag" in options:

                comet_cfg = os.path.join(cwd, "comet_settings.xml")

                cfg_txt = None
                with open("/opt/gladiator/comet_settings_template.xml", "r") as fh:
                    cfg_txt = fh.read()

                for k in param_map:
                    placeholder = param_map[k]
                    target = str(self.data[k])
                    if k == "precursor_tolerance":
                        target = "{:.2f}".format(float(target))
                    cfg_txt = cfg_txt.replace(placeholder, target)

                cfg_txt = cfg_txt.replace("DATABASE_FASTA_FILE", decoy_db_file)

                with open(comet_cfg,"w") as fh:
                    fh.write(cfg_txt)

                if library_files:

                    if self.scan and "comet_peptides" in self.scan and self.scan["comet_peptides"]:
                        None
                        #print("Detected existing Comet spectra search results")

                    else:

                        # Cleanup

                        filenames = \
                            ["interact_comet_pep-MODELS.html",
                             "interact_comet_pep.xml",
                             "interact_comet_pep.xml.index",
                             "interact_comet_pep.xml.RTcoeff",
                             "interact_comet_pep.xml.RTstats"]

                        for filename in filenames:
                            filepath = os.path.join(cwd, filename)
                            if os.path.isfile(filepath):
                                os.remove(filepath)

                        # Run Comet

                        existing_pep_xmls = None
                        if self.scan and "samples_generated_by_comet" in self.scan:
                            existing_pep_xmls = self.scan["samples_generated_by_comet"]

                        phaseID = analysis_state.createPhase(
                            "Speclib - Matching sequences [Comet]")
                        workflow.runComet(
                            event,
                            log_fh,
                            analysis_state.getPhaseData(phaseID),
                            cwd,
                            comet_cfg,
                            library_files,
                            existing_pep_xmls,
                            "interact_comet_pep.xml",
                            cwd)

                    if event.kill_flag:
                        return

                    libmethod_pepXMLs.append("interact_comet_pep.xml")

                if pseudospectrafiles:

                    if self.scan and "comet_pseudo_peptides" in self.scan and self.scan["comet_pseudo_peptides"]:
                        None
                        # print("Detected existing Comet pseudo spectra search results")

                    else:

                        # Cleanup

                        filenames = \
                            ["interact_comet_pseudo_pep-MODELS.html",
                             "interact_comet_pseudo_pep.xml",
                             "interact_comet_pseudo_pep.xml.index",
                             "interact_comet_pseudo_pep.xml.RTcoeff",
                             "interact_comet_pseudo_pep.xml.RTstats"]

                        for filename in filenames:
                            filepath = os.path.join(cwd, filename)
                            if os.path.isfile(filepath):
                                os.remove(filepath)

                        existing_pep_xmls = None
                        if self.scan and "pseudo_samples_generated_by_comet" in self.scan:
                            existing_pep_xmls = self.scan["pseudo_samples_generated_by_comet"]

                        # Run Comet for speudospectrafiles

                        phaseID = analysis_state.createPhase(
                            "Pseudospeclib - Matching sequences [Comet]")
                        workflow.runComet(
                            event,
                            log_fh,
                            analysis_state.getPhaseData(phaseID),
                            cwd,
                            comet_cfg,
                            pseudospectrafiles,
                            existing_pep_xmls,
                            "interact_comet_pseudo_pep.xml",
                            cwd)

                    if event.kill_flag:
                        return

                    libfreemethod_pepXMLs.append(
                        "interact_comet_pseudo_pep.xml")

            if "use_xtandem_flag" in options:
                xtandem_cfg = os.path.join(cwd, "xtandem_settings.xml")

                cfg_txt = None
                with open("/opt/gladiator/xtandem_settings_template.xml", "r") as fh:
                    cfg_txt = fh.read()

                for k in param_map:
                    placeholder = param_map[k]
                    target = str(self.data[k])
                    if k == "precursor_tolerance":
                        target = "{:.2f}".format(float(target))
                    cfg_txt = cfg_txt.replace(placeholder, target)

                with open(xtandem_cfg,"w") as fh:
                    fh.write(cfg_txt)

                if library_files:

                    if self.scan and "xtandem_peptides" in self.scan and self.scan["xtandem_peptides"]:
                        None
                        # print("Detected existing X!Tandem spectra search results")
                    else:

                        # Cleanup

                        filenames = [
                            "interact_xtandem_pep-MODELS.html   ",
                            "interact_xtandem_pep.xml",
                            "interact_xtandem_pep.xml.index",
                            "interact_xtandem_pep.xml.RTcoeff   ",
                            "interact_xtandem_pep.xml.RTstats   "]

                        for filename in filenames:
                            filepath = os.path.join(cwd, filename)
                            if os.path.isfile(filepath):
                                os.remove(filepath)

                        # Run X!Tandem

                        existing_pep_xmls = None
                        if self.scan and "samples_generated_by_xtandem" in self.scan:
                            existing_pep_xmls = self.scan["samples_generated_by_xtandem"]

                        phaseID = analysis_state.createPhase(
                            "Speclib - Matching sequences [X!Tandem]")
                        workflow.runXTandem(
                            event,
                            log_fh,
                            analysis_state.getPhaseData(phaseID),
                            cwd,
                            xtandem_cfg,
                            decoy_db_file,
                            library_files,
                            existing_pep_xmls,
                            "interact_xtandem_pep.xml",
                            cwd,
                            delete_tmp_files_flag)

                    if event.kill_flag:
                        return

                    libmethod_pepXMLs.append("interact_xtandem_pep.xml")

                if pseudospectrafiles:

                    if self.scan and "xtandem_pseudo_peptides" in self.scan and self.scan["xtandem_pseudo_peptides"]:
                        None
                        # print("Detected existing X!Tandem pseudo spectra search results")

                    else:

                        # Cleanup

                        filenames = [
                            "interact_xtandem_pseudo_pep-MODELS.html",
                            "interact_xtandem_pseudo_pep.xml",
                            "interact_xtandem_pseudo_pep.xml.index",
                            "interact_xtandem_pseudo_pep.xml.RTcoeff",
                            "interact_xtandem_pseudo_pep.xml.RTstats"]

                        for filename in filenames:
                            filepath = os.path.join(cwd, filename)
                            if os.path.isfile(filepath):
                                os.remove(filepath)

                        # Run X!Tandem for speudospectrafiles

                        existing_pep_xmls = None
                        if self.scan and "pseudo_samples_generated_by_xtandem" in self.scan:
                            existing_pep_xmls = self.scan["pseudo_samples_generated_by_xtandem"]

                        phaseID = analysis_state.createPhase(
                            "Pseudospeclib - Matching sequences [X!Tandem]")
                        workflow.runXTandem(
                            event,
                            log_fh,
                            analysis_state.getPhaseData(phaseID),
                            cwd,
                            xtandem_cfg,
                            decoy_db_file,
                            pseudospectrafiles,
                            existing_pep_xmls,
                            "interact_xtandem_pseudo_pep.xml",
                            cwd,
                            delete_tmp_files_flag)

                    if event.kill_flag:
                        return

                    libfreemethod_pepXMLs.append(
                        "interact_xtandem_pseudo_pep.xml")

            if self.scan and "speclib_cons" in self.scan and self.scan["speclib_cons"]:
                None
                #print("Detected existing Spectrum library.")

            else:

                phaseID = analysis_state.createPhase("Building Library")

                # Cleanup

                filenames = \
                    ["SpecLib_cons_decoy.TraML",
                     "SpecLib_cons_openswath.tsv",
                     "SpecLib_cons.pepidx",
                     "SpecLib_cons.spidx",
                     "SpecLib_cons.splib",
                     "SpecLib_cons.sptxt",
                     "SpecLib_cons.TraML",
                     "SpecLib_libfree.pepidx",
                     "SpecLib_libfree.spidx",
                     "SpecLib_libfree.splib",
                     "SpecLib_libfree.sptxt",
                     "SpecLib_lib.pepidx",
                     "SpecLib_lib.spidx",
                     "SpecLib_lib.splib",
                     "SpecLib_lib.sptxt",
                     "SpecLib_merged.pepidx",
                     "SpecLib_merged.spidx",
                     "SpecLib_merged.splib",
                     "SpecLib_merged.sptxt",
                     "spectrast.log"]

                for filename in filenames:
                    filepath = os.path.join(cwd, filename)
                    if os.path.isfile(filepath):
                        os.remove(filepath)

                for method in ["lib", "libfree"]:
                    mayu_folder = method+"_mayu"
                    directory = os.path.join(cwd, mayu_folder)
                    if os.path.isdir(directory):
                        shutil.rmtree(directory)

                workflow.buildlib(
                    event,
                    log_fh,
                    analysis_state.getPhaseData(phaseID),
                    cwd,
                    "DECOY_",
                    decoy_db_file,
                    str(max_threads),
                    libmethod_pepXMLs,
                    libfreemethod_pepXMLs,
                    "/opt/gladiator/iRT.txt",
                    "swath-windows.txt",
                    pvalue,
                    pvalue,
                    swaths_min,
                    swaths_max)

            if event.kill_flag:
                return

            if self.scan and "matrices" in self.scan and self.scan["matrices"]:
                None
                #print("Detected existing matrices.")

            else:

                # Cleanup

                postfixes = ["-DIA_cutoffs.txt",
                             "-DIA_dscores_top_decoy_peaks.txt",
                             "-DIA_dscores_top_target_peaks.txt",
                             "-DIA_full_stat.csv",
                             "-DIA_mayu.csv",
                             "-DIA_mayu.cutoff",
                             "-DIA_mayu.fasta",
                             "-DIA_qvalues.txt",
                             "-DIA_report.pdf",
                             "-DIA_scorer.bin",
                             "-DIA_summary_stat.csv",
                             "-DIA_svalues.txt",
                             "-DIA.tsv",
                             "-DIA_weights.txt",
                             "-DIA_with_dscore-0_0-None.tr",
                             "-DIA_with_dscore.csv",
                             "-DIA_with_dscore_filtered.csv"]

                for filename in sample_files:
                    filepath = os.path.join(cwd, os.path.basename(filename))
                    for postfix in postfixes:
                        if filepath.endswith(postfix) and os.path.isfile(filepath):
                            os.remove(filepath)

                filenames = \
                    ["DIA-analysis-result.csv",
                     "DIA-peptide-matrix.tsv",
                     "DIA-protein-matrix.tsv"]

                for filename in filenames:
                    filepath = os.path.join(cwd, filename)
                    if os.path.isfile(filepath):
                        os.remove(filepath)

                # Run OpenSWATH and building the DIA matrices
                phaseID = analysis_state.createPhase(
                    "Searching peptides from DIA spectrum")

                workflow.buildDIAMatrix(
                    event,
                    log_fh,
                    analysis_state.getPhaseData(phaseID),
                    cwd,
                    sample_files,
                    "truncated-swath-windows.txt",
                    trig_target_pvalue,
                    trig_max_pvalue,
                    str(max_threads),
                    "/opt/gladiator/iRTAssayLibrary.TraML",
                    None)  # TODO: insert design file here

            if event.kill_flag:
                return

            with notification_lock:
                notifications.append(
                    {
                        "title": "Analysis completed",
                        "text": "Analysis of project " + self.project_folder + " was completed successfully.",
                    }
                )

            # Create matrix.tsv and annotations.tsv

            if 'annotate_peptides' in options:
                if not (self.scan and "annotations" in self.scan and self.scan["annotations"]):

                    cwd = os.path.join(result_root, project)

                    assign_ambiguous = self.data["annotation_assign_ambiguous"]
                    annotation_files = [self.data["annotation_filename"]]
                    merge_unimods = False
                    if "annotation_merge_unimods" in self.data:
                        merge_unimods = self.data["annotation_merge_unimods"]

                    with open(os.path.join(cwd, "annotation_log.txt"), "w") as log_fh:

                        phaseID = analysis_state.createPhase(
                            "Loading annotation file")

                        dicts = None
                        for filename in annotation_files:
                            id_key = self.data['annotation_id_column']
                            dicts = annotation.load_annotations(
                                event,
                                log_fh,
                                analysis_state.getPhaseData(phaseID),
                                id_key,
                                os.path.join(result_root, filename),
                                dicts)

                            if event.kill_flag:
                                return

                        phaseID = analysis_state.createPhase(
                            "Loading detected peptides")

                        df = None

                        df = annotateSwath2stats.read_swaths2stats_peptide_table(
                            event,
                            log_fh,
                            analysis_state.getPhaseData(phaseID),
                            os.path.join(result_root, project, "DIA-peptide-matrix.tsv"))

                        # df = annotation.read_openswathfile(
                        #     event,
                        #     log_fh,
                        #     analysis_state.getPhaseData(phaseID),
                        #     os.path.join(result_root, project, "DIA-analysis-result.csv"))

                        if event.kill_flag:
                            return

                        # if "sample_rename_map_file" in self.data and self.data["sample_rename_map_file"]:
                        #     filename = self.data["sample_rename_map_file"]
                        #     rename_df = pd.read_csv(os.path.join(filename), sep='\t')
                        #     rename_df = rename_df[rename_df["ORIGINAL"].isin(df.columns.values)]
                        #     rename_map = dict(rename_df.values.tolist())
                        #     sorted_orig_samplenames = rename_df["ORIGINAL"].tolist()
                        #     sorted_conv_samplenames = rename_df["RENAMED"].tolist()
                        #     other_cols = [x for x in df.columns.values if x not in sorted_orig_samplenames]
                        #     df = df.reindex(columns = other_cols + sorted_orig_samplenames)
                        #     df.rename(columns=rename_map, inplace=True)

                        phaseID = analysis_state.createPhase("Annotating peptides")

                        ambiguous_threshold = 0
                        if assign_ambiguous:
                            ambiguous_threshold = 2

                        annotateSwath2stats.annotate_df(
                            event,
                            log_fh,
                            analysis_state.getPhaseData(phaseID),
                            df,
                            dicts,
                            os.path.join(result_root, project, "matrix.tsv"),
                            os.path.join(result_root, project, "annotations.tsv"),
                            ambiguous_threshold,
                            merge_unimods,
                            None) # Contaminants

                        # annotation.annotate_df(
                        #     event,
                        #     log_fh,
                        #     analysis_state.getPhaseData(phaseID),
                        #     df,
                        #     dicts,
                        #     os.path.join(result_root, project, "matrix.tsv"),
                        #     os.path.join(result_root, project, "annotations.tsv"),
                        #     ambiguous_threshold,
                        #     merge_unimods)

                        if event.kill_flag:
                            return

                # Draw graphs

                if not (self.scan and "figures" in self.scan and self.scan["figures"]):

                    figures_path = os.path.join(result_root, project, "figures")

                    if not os.path.exists(figures_path):
                        os.mkdir(figures_path)

                    annotation_df = pd.read_csv(os.path.join(
                        result_root, project, "annotations.tsv"), sep="\t", index_col=0)

                    with open(os.path.join(cwd, "figures", "figcfg.tsv"), "w") as fh:
                        dict_variable = {
                            key: None for key in annotation_df.columns.tolist()}
                        json.dump(dict_variable, fh)

                with notification_lock:
                    notifications.append(
                        {
                            "title": "Annotation completed",
                            "text": "Annotation of project " + project + " was completed successfully.",
                        }
                    )

        return


state = {}


def cancel_analysis(request):
    requestData = json.loads(request.body)
    project = requestData['project']
    with big_lock:
        if project in analysis_events:
            analysis_events[project].kill()


def get_status(request):

    status = {}

    with big_lock:

        for project in annotation_threads:

            if annotation_threads[project].isAlive():
                if project not in status:
                    status[project] = {}
                    status[project]["Annotation Running"] = "No"
                    status[project]["Analysis Running"] = "No"

                status[project]["Annotation Running"] = "Yes"
                folder = annotation_threads[project].get_project_folder()
                status[project]["Annotation Folder"] = folder

        for project in analysis_threads:

            if analysis_threads[project].isAlive():
                if project not in status:
                    status[project] = {}
                    status[project]["Annotation Running"] = "No"
                    status[project]["Analysis Running"] = "No"

                status[project]["Analysis Running"] = "Yes"
                folder = analysis_threads[project].get_project_folder()
                status[project]["Analysis Folder"] = folder

    status["Status"] = "Success"
    return status


def remove_analysisrun(request):
    data = json.loads(request.body)
    if "folder" in data:
        folder = data["folder"]
        fullname = os.path.join(result_root, folder)
        if folder and os.path.isdir(fullname):
            shutil.rmtree(fullname)
            return ({"Status": "Success"})
    return ({"Status": "Failed"})


def getProgress(request):

    with big_lock:

        progress = {}

        for project in state:
            progress[project] = {}
            for k in ["analysis"]:
                if k in state[project]:
                    progress[project][k] = state[project][k].load()

        return progress


def load_config(request):
    runData = json.loads(request.body)
    analysis_name = runData['analysis_name']
    project_path = os.path.basename(analysis_name)

    cwd = os.path.join(result_root, project_path)
    cfg = None

    cfg_filename = os.path.join(cwd, "config.txt")
    cfg_exists = os.path.isfile(cfg_filename)
    if cfg_exists:
        with open(cfg_filename, "r") as fh:
            cfg = json.load(fh)
            return ({"Status": "Success", "config": cfg})
    else:
        return ({"Status": "Failed", "config": cfg})


def rerun(request):
    runData = json.loads(request.body)
    analysis_name = runData['analysis_name']

    analysis_events[analysis_name] = workflow.Event()

    global state

    if analysis_name not in state:
        state[analysis_name] = {}

    if "analysis" in state[analysis_name]:
        state[analysis_name]["analysis"].clear()
    else:
        state[analysis_name]["analysis"] = workflow.State(analysis_name)

    # START ANALYSIS THREAD
    global analysis_threads
    analysis_threads[analysis_name] = AnalysisThread(None, analysis_name)
    analysis_threads[analysis_name].start()
    return ({"Status": "Success"})


def run(request):

    runData = json.loads(request.body)

    analysis_name = runData['analysis_name']

    with big_lock:

        # CREATE ANALYSIS FOLDER

        now = datetime.now()
        timestamp = now.strftime("%Y%m%d-%H%M%S")

        if analysis_name:
            result_path = os.path.join(result_root, analysis_name)
            if os.path.exists(result_path):
                analysis_name += "-" + timestamp
                result_path = os.path.join(result_root, analysis_name)

            if os.path.exists(result_path):
                #print("Analysis folder already exists " + result_path)
                return

            os.mkdir(result_path)

        else:
            analysis_name = "analysis-" + timestamp
            result_path = os.path.join(result_root, analysis_name)

            if os.path.exists(result_path):
                #print("Analysis folder already exists " + result_path)
                return

            os.mkdir(result_path)

        analysis_events[analysis_name] = workflow.Event()

        global state

        if analysis_name not in state:
            state[analysis_name] = {}

        if "analysis" in state[analysis_name]:
            state[analysis_name]["analysis"].clear()
        else:
            state[analysis_name]["analysis"] = workflow.State(analysis_name)

        # START ANALYSIS THREAD
        global analysis_threads
        analysis_threads[analysis_name] = AnalysisThread(
            runData, analysis_name)
        analysis_threads[analysis_name].start()
        return ({"Status": "Success"})


def list_projects(request):

    #with big_lock:

    global result_root

    if not os.path.exists(result_root):
        return []

    # Scan folders
    projects = []
    for basename in os.listdir(result_root):
        fullname = os.path.join(result_root, basename)
        if os.path.isdir(fullname):

            cfg_filename = os.path.join(fullname, "config.txt")
            cfg = None
            cfg_exists = os.path.isfile(cfg_filename)

            if cfg_exists:
                with open(cfg_filename, "r") as fh:
                    cfg = json.load(fh)

            is_pepmatrix = os.path.isfile(
                os.path.join(fullname, "DIA-peptide-matrix.tsv"))
            is_protmatrix = os.path.isfile(
                os.path.join(fullname, "DIA-protein-matrix.tsv"))

            exists_pepmatrixfile = os.path.isfile(
                os.path.join(fullname, "matrix.tsv"))
            exists_annotationfile = os.path.isfile(
                os.path.join(fullname, "annotations.tsv"))

            is_running = False
            if basename in analysis_threads:
                is_running = analysis_threads[basename].isAlive()

            figcfg = None
            fig_path = os.path.join(fullname, "figures")
            if os.path.isdir(fig_path):

                figcfg_filename = os.path.join(fig_path, "figcfg.tsv")
                figcfg_exists = os.path.isfile(figcfg_filename)
                if figcfg_exists:
                    with open(figcfg_filename, "r") as fh:
                        figcfg = json.load(fh)
                        
            projects.append({
                "fullname": fullname,
                "folder": basename,
                "complete": is_pepmatrix and is_protmatrix and cfg_exists,
                "annotated": exists_pepmatrixfile and exists_annotationfile,
                "figures": figcfg,
                "config": cfg,
                "selected": False,
                "running": is_running,
            })

    return projects


def listDirectory(request):

    path = json.loads(request.body)

    root = "/" + "/".join(["data"] + [x.strip() for x in path.split("/") if x])

    if not os.path.exists(root):
        return {"files": []}

    l = []

    basenames = os.listdir(root)

    for basename in basenames:
        fullname = os.path.join(root, basename)
        filetype = "unknown"
        filesize = None
        if os.path.isfile(fullname):
            filetype = "file"
            filesize = os.path.getsize(fullname)

        if os.path.islink(fullname):
            filetype = "file"

        if os.path.isdir(fullname):
            filetype = "folder"

        l.append(
            {
                "root": root,
                "basename": basename,
                "type": filetype,
                "fullpath": fullname,
                "filesize": filesize,
                "selected": False
            }
        )

    return {"files": l}


def get_tsv_headers(request):
    data = json.loads(request.body)
    df = pd.read_csv(data["file"], sep='\t', nrows=1)
    headers = df.columns.tolist()
    return (headers)


def get_notifications(request):
    result = None
    with notification_lock:
        global notifications
        result = notifications
        notifications = []
    return (result)


def scan(request):
    result = workflow.scan_project_phases("PD92-libfree", result_root)
    return result


def render_figure(request):

    global result_root
    
    data = json.loads(request.body)

    with big_lock:

        project = data['project']
        annotation_field = data['field']

        figcfg_filename = os.path.join(result_root, project, "figures", "figcfg.tsv")
        figcfg_exists = os.path.isfile(figcfg_filename)
        figcfg = {}
        if figcfg_exists:
            with open(figcfg_filename, "r") as fh:
                figcfg = json.load(fh)
        else:
            print ("Whoops, figcfg.tsv file not found.")
            return None

        if not annotation_field in figcfg:
            figcfg[annotation_field] = {"barplot": None, "pieplot":None}
        elif not figcfg[annotation_field]:
            figcfg[annotation_field] = {"barplot": None, "pieplot":None}

#        if figcfg[annotation_field]:
#            if os.path.isfile(figcfg[annotation_field]):
#                return figcfg[annotation_field]
#            else:
#                print ("Whoops, figure " + figcfg[annotation_field] + " does not exist. Trying to regenerate.")

        matrix_df = pd.read_csv(os.path.join(
            result_root, project, "matrix.tsv"), sep="\t", index_col=0)
        annotation_df = pd.read_csv(os.path.join(
            result_root, project, "annotations.tsv"), sep="\t", index_col=0)

        df = matrix_df.merge(annotation_df[[annotation_field]],
                            left_index=True,
                            right_index=True,
                            how='inner')

        mp_df = df.groupby(annotation_field, as_index=True).agg('sum')
        mp_df = mp_df.drop(index=['unknown', 'ambiguous'])

        # mp_df = mp_df.loc[:, natsorted(list(mp_columns & mg_columns))]

        count_df = mp_df.sum(axis=1)
        count_df = count_df.div(count_df.sum(axis=0), axis=0)
        mp_df = mp_df.loc[count_df[count_df > 0.01].index, :]

        mp_df = mp_df.div(mp_df.sum(axis=0), axis=1)

        count_df = count_df[mp_df.index]

        # pieplot

        entry = {
            "type": "pie",
            "values": count_df.iloc[0:].tolist(),
            "labels": count_df.index.tolist(),
            "texttemplate": "%{label} (%{percent})",
            "textposition": "outside",
            "hole": .4,
            "marker" : {"line":{"width":0, "color":"rgb(0,0,0)"}},
        }

        fig = {
            "layout": {
                "template": "plotly_white",
                # "title": "MP (>= 0.1%)",
                #"barmode": "stack",
                "showlegend":False,
                "font": {
                    "family": 'Arial, sans-serif',
                    "size": 16,
                    "color": '#7f7f7f',
                },
                "margin": {
                    "b": 400,
                    "pad": 4,
                },

                "yaxis": {
                    # "tickformat": ',.0%',
                },

                "xaxis": {
                    "dtick": 1,
                    "type": "category",
                    "tickangle": -90,
                },
                "legend": {
                    "orientation": "v",
                    "font": {
                        "size" :18,
                    },
                },

            },
            "data": [entry]
        }
        fig["layout"]["yaxis"]["tickformat"] = ',.0%'

        tmp_fd1 = tempfile.NamedTemporaryFile(dir=os.path.join(result_root, project, "figures"),
                                            mode="r+",
                                            suffix='.html',
                                            delete=False)
        tmp_fd1.close()

        py.plot(fig,
                filename=os.path.join(result_root, project, "figures", tmp_fd1.name),
                auto_open=False)        

        figcfg[annotation_field]["pieplot"] = tmp_fd1.name

        #barplot

        data = []

        for genus in mp_df.index.tolist():

            entry = {
                "type": "bar",
                "name": genus,
                "x": mp_df.loc[genus].index.tolist(),
                "y": mp_df.loc[genus].tolist(),
            }

            data.append(entry)

        fig = {
            "layout": {
                "template": "plotly_white",
                # "title": "MP (>= 0.1%)",
                "barmode": "stack",
                "font": {
                    "family": 'Arial, sans-serif',
                    "size": 12,
                    "color": '#7f7f7f',
                },
                "margin": {
                    "b": 400,
                    "pad": 4,
                },

                "yaxis": {
                    # "tickformat": ',.0%',
                },

                "xaxis": {
                    "dtick": 1,
                    "type": "category",
                    "tickangle": -90,
                },
                "legend": {
                    "orientation": "v",
                    "font": {
                        "size" :18,
                    },
                },

            },
            "data": data
        }
        fig["layout"]["yaxis"]["tickformat"] = ',.0%'

        tmp_fd2 = tempfile.NamedTemporaryFile(dir=os.path.join(result_root, project, "figures"),
                                            mode="r+",
                                            suffix='.html',
                                            delete=False)
        tmp_fd2.close()

        py.plot(fig,
                filename=os.path.join(result_root, project, "figures", tmp_fd2.name),
                auto_open=False)        

        figcfg[annotation_field]["barplot"] = tmp_fd2.name

        with open(figcfg_filename, "w") as fh:
            json.dump(figcfg, fh)

        return None


def main(global_config, **settings):
    config = Configurator(settings=settings)

    config.add_route('list', '/list')
    config.add_view(listDirectory, route_name='list', renderer="json")

    config.add_route('get_progress', '/get_progress')
    config.add_view(getProgress, route_name='get_progress', renderer="json")

    config.add_route('get_status', '/get_status')
    config.add_view(get_status, route_name='get_status', renderer="json")

    config.add_route('list_projects', '/list_projects')
    config.add_view(list_projects, route_name='list_projects', renderer="json")

    config.add_route('remove_analysisrun', '/remove_analysisrun')
    config.add_view(remove_analysisrun,
                    route_name='remove_analysisrun', renderer="json")

    config.add_route('cancel_analysis', '/cancel_analysis')
    config.add_view(cancel_analysis,
                    route_name='cancel_analysis', renderer="json")

    config.add_route('run', '/run')
    config.add_view(run, route_name='run', renderer="json")

    config.add_route('rerun', '/rerun')
    config.add_view(rerun, route_name='rerun', renderer="json")

    config.add_route('load_config', '/load_config')
    config.add_view(load_config, route_name='load_config', renderer="json")

    config.add_route('get_notifications', '/get_notifications')
    config.add_view(get_notifications,
                    route_name='get_notifications', renderer="json")

    config.add_route('get_tsv_headers', '/get_tsv_headers')
    config.add_view(get_tsv_headers,
                    route_name='get_tsv_headers', renderer="json")

    config.add_route('scan', '/scan')
    config.add_view(scan, route_name='scan', renderer="json")

    config.add_route('render_figure', '/render_figure')
    config.add_view(render_figure, route_name='render_figure', renderer="json")

    global result_root
    config.add_static_view(name=result_root, path=result_root)

    config.add_static_view(name='/', path="assets/")

    ticker.start()

    os.chdir(result_root)

    return config.make_wsgi_app()
