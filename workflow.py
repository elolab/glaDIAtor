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

from progress import Progress

class NonZeroReturnValueException(Exception):
    def __init__(self, returnvalue, msg):
        self.msg = msg
        self.returnvalue = returnvalue
        return


class Event:
    def __init__(self):
        self.kill_flag = False
    
    def kill(self):
        self.kill_flag = True

    def reset(self):
        self.kill_flag = False


class State:
    
    def __init__(self, project):
        self.i = 0
        self.phases = {}
        self.complete = False
        self.project = project
        self.lock = threading.Lock()

    def mark_complete(self):
        self.complete = True

    def clear(self):
        with self.lock:
            self.i = 0
            self.phases = {}

    def createPhase(self, name):
        with self.lock:
            ID = str(self.i)
            self.phases[ID] = {
                "name" : name,
                "type" : "percentage-indicator", 
                "percentage" : 0,
            }
            self.i += 1
            return ID


    def getPhaseData(self, ID):
        with self.lock:
            return self.phases[ID]

    def updatePhase(self, ID, data):
        with self.lock:
            self.phases[ID] = data
            return

    def save(self, data):
        with self.lock:
            self.phases = data

    def load(self):
        with self.lock:
            return ({
                "phases": self.phases,
                "complete":self.complete,
                })


def logline(log_fh, text):
    log_fh.write(text+"\n")
    log_fh.flush()
    return


def install_ThermoRawFileParser(event, log_fh, progressData):

    progress = Progress(progressData, "percentage-indicator")

    step = 0
    n_steps = 3

    cmd = [
        "wget",
        "https://github.com/compomics/ThermoRawFileParser/releases/download/v1.3.4/ThermoRawFileParser.zip"

    ]

    logline(log_fh, "Running pipeline command: "+" ".join(cmd))
    proc = subprocess.Popen(cmd, cwd="/root", stdout=log_fh, stderr=log_fh)

    while True:
        try:
            proc.wait(1)
            break
        except subprocess.TimeoutExpired:
            if event and event.kill_flag:
                proc.kill()
                progress.fail({
                    "cmd": "Killed",
                    "returncode": -1, 
                })

                return

    log_fh.flush()
    if proc.returncode != 0:
        progress.fail({
            "cmd": "Wget returned error",
            "returncode": proc.returncode, 
        })
        raise NonZeroReturnValueException(proc.returncode, 'wget')

    step += 1
    progress.update_n_of_m(step, n_steps)


    if not os.path.exists("/root/ThermoRawFileParser.zip"):
        print ("/root/ThermoRawFileParser.zip does not exist.")
        progress.fail({
            "cmd": "Only mzml is available with hermoRawFileParser conversion method",
            "returncode": -1, 
        })
        return

    os.mkdir("/opt/ThermoRawFileParser")

    step += 1
    progress.update_n_of_m(step, n_steps)
    
    cmd = [
        "unzip",
        "/root/ThermoRawFileParser.zip"

    ]

    logline(log_fh, "Running pipeline command: "+" ".join(cmd))
    proc = subprocess.Popen(cmd, cwd="/opt/ThermoRawFileParser", stdout=log_fh, stderr=log_fh)

    while True:
        try:
            proc.wait(1)
            break
        except subprocess.TimeoutExpired:
            if event and event.kill_flag:
                proc.kill()
                progress.fail({
                    "cmd": "Killed",
                    "returncode": -1, 
                })
                return

    log_fh.flush()
    if proc.returncode != 0:
        progress.fail({
            "cmd": "ThermoRawFileParser returned error",
            "returncode": proc.returncode, 
        })
        raise NonZeroReturnValueException(proc.returncode, '/ThermoRawFileParser/ThermoRawFileParser.exe')

    step += 1
    progress.update_n_of_m(step, n_steps)
    
    progress.update_n_of_m(n_steps, n_steps)
    progress.ready()

    return
    

def convertRAW_ThermoRawFileParser(event, log_fh, progressData, cwd, filenames, peak_picking, to_format, folder):

    progress = Progress(progressData, "percentage-indicator")

    if to_format != "mzml":
        progress.fail({
            "cmd": "Only mzml is available with hermoRawFileParser conversion method",
            "returncode": -1, 
        })
        return False

    
    step = 0
    n_steps = len(filenames)

    if not os.path.exists(os.path.join(cwd, folder)):
        os.mkdir(os.path.join(cwd, folder))

    for filename in filenames:

        outputfilename = os.path.splitext(os.path.basename(filename))[0] + ".mzML"
        
        cmd = [
            "mono",
            "/opt/ThermoRawFileParser/ThermoRawFileParser.exe",
            "-f=1",
            "-m=0",
            "--noZlibCompression",
            "-i="+filename,
            "-b="+os.path.join(cwd, folder, outputfilename)
        ]

        if not peak_picking:
            cmd.append("-p")

        logline(log_fh, "Running pipeline command: "+" ".join(cmd))
        proc = subprocess.Popen(cmd, cwd=cwd, stdout=log_fh, stderr=log_fh)

        while True:
            try:
                proc.wait(1)
                break
            except subprocess.TimeoutExpired:
                if event and event.kill_flag:
                    proc.kill()
                    return

        log_fh.flush()
        if proc.returncode != 0:
            raise NonZeroReturnValueException(proc.returncode, '/ThermoRawFileParser/ThermoRawFileParser.exe')

        step += 1
        progress.update_n_of_m(step, n_steps)
    
    progress.update_n_of_m(n_steps, n_steps)
    progress.ready()

    return True



def convertRAWqtofpeakpicker(event, log_fh, progressData, cwd, filenames, folder):

    progress = Progress(progressData, "percentage-indicator")

    step = 0
    n_steps = len(filenames)

    if not os.path.exists(os.path.join(cwd, folder)):
        os.mkdir(os.path.join(cwd, folder))

    for filename in filenames:

        outputfilename = os.path.splitext(os.path.basename(filename))[0] + ".mzXML"

        cmd = [
            "wine",
            "qtofpeakpicker.exe",
            "--resolution=2000",
            "--area=1",
            "--threshold=1",
            "--smoothwidth=1.1",
            "--in", filename,
            "--out", os.path.join(cwd, folder, outputfilename)
        ]

        run_env = os.environ.copy()
        run_env["WINEPREFIX"] = "/wineprefix64"
        run_env["WINEDEBUG"]="-all,err+all"
        run_env["WINEPATH"]="C:\pwiz"

        logline(log_fh, "Running pipeline command: "+" ".join(cmd))
        proc = subprocess.Popen(cmd, env=run_env, cwd=cwd, stdout=log_fh, stderr=log_fh)

        while True:
            try:
                proc.wait(1)
                break
            except subprocess.TimeoutExpired:
                if event and event.kill_flag:
                    proc.kill()
                    return

        log_fh.flush()
        if proc.returncode != 0:
            raise NonZeroReturnValueException(proc.returncode, 'qtofpeakpicker.exe')

        step += 1
        progress.update_n_of_m(step, n_steps)
    
    progress.update_n_of_m(n_steps, n_steps)
    progress.ready()

    return True

def convertRAW(event, log_fh, progressData, cwd, filenames, peak_picking, to_format, folder):

    progress = Progress(progressData, "percentage-indicator")

    step = 0
    n_steps = len(filenames)

    for filename in filenames:

        cmd = [
            "wine64",
            "msconvert.exe",
            filename,
        ]

        if peak_picking:
            cmd.extend(["--filter", "peakPicking<cwt> true 1-"])

        if to_format == "mzxml":
            cmd.append("--mzXML")
        elif to_format == "mzml":
            cmd.append("--mzML")
        else:
            print ("Unknown format " + to_format)
            return False

        cmd.extend(["-o", folder])
        
        run_env = os.environ.copy()
        run_env["WINEPREFIX"] = "/wineprefix64"
        run_env["WINEDEBUG"]="-all,err+all"
        run_env["WINEPATH"]="C:\pwiz"

        logline(log_fh, "Running pipeline command: "+" ".join(cmd))
        proc = subprocess.Popen(cmd, env=run_env, cwd=cwd, stdout=log_fh, stderr=log_fh)

        while True:
            try:
                proc.wait(1)
                break
            except subprocess.TimeoutExpired:
                if event and event.kill_flag:
                    proc.kill()
                    return


        log_fh.flush()
        if proc.returncode != 0:
            raise NonZeroReturnValueException(proc.returncode, 'msconvert.exe')

        step += 1
        progress.update_n_of_m(step, n_steps)
    
    progress.update_n_of_m(n_steps, n_steps)
    progress.ready()

    return True


def runDiaumpire(event, log_fh, progressData, cwd, DIAfiles, outDir, threads):

    if not os.path.isdir(os.path.join(cwd, "libfree")):
        os.mkdir(os.path.join(cwd, "libfree"))

    if not os.path.isdir(os.path.join(cwd, "libfree-pseudospectra")):
        os.mkdir(os.path.join(cwd, "libfree-pseudospectra"))

    diaumpire_cfg = os.path.join(cwd, "diaumpire-params.txt")

    cfg_txt = None
    with open("/opt/gladiator/diaumpire-params-template.txt", "r") as fh:
        cfg_txt = fh.read()

    cfg_txt = cfg_txt.replace("THREADS", str(threads))

    with open(diaumpire_cfg,"w") as fh:
        fh.write(cfg_txt)

    phase_1_n = len(DIAfiles)
    phase_1_steps = 0
    n_steps = phase_1_n * (2+3) # estimate

    progress = Progress(progressData, "percentage-indicator")
    for i, DIAfile in enumerate(DIAfiles):
        
        mzXMLfile = os.path.join(cwd, "libfree/"+os.path.basename(os.path.splitext(DIAfile)[0]+".mzXML"))

        if not os.path.exists(mzXMLfile):
        
            cmd = [
                "/opt/tpp/bin/msconvert",
            ]
            
            params = [
                #"--filter", '"titleMaker <RunId>.<ScanNumber>.<ScanNumber>.<ChargeState>"',
                "--32",
                "--zlib",
                '--filter', '"peakPicking false 1-"',
                "--mzXML",
                "-o", os.path.join(cwd, "libfree"),
            ]
            
            cmd.append(DIAfile)
            cmd.extend(params)
            
            logline(log_fh, "Running pipeline command: "+" ".join(cmd))
            proc = subprocess.Popen(cmd, cwd=cwd, stdout=log_fh, stderr=log_fh)

            while True:
                try:
                    proc.wait(1)
                    break
                except subprocess.TimeoutExpired:
                    if event and event.kill_flag:
                        proc.kill()
                        return

            log_fh.flush()
            if proc.returncode != 0:
                raise NonZeroReturnValueException(proc.returncode, 'msconvert')

        else:
            print (mzXMLfile + " already exists, skipping msconvert")

        phase_1_steps += 1
        progress.update_n_of_m(phase_1_steps, n_steps)

    filenames = []
    for root, dirs, files in os.walk(os.path.join(cwd, "libfree")):
        for filename in files:
            if filename.endswith(".mzXML"):
                filenames.append(filename)

    phase_2_n = len(filenames)
    phase_2_steps = 0
    n_steps = phase_1_steps + phase_2_n + 3 * phase_2_n # estimate

    for i, filename in enumerate(filenames):
        basename = os.path.basename(filename)
        mgf1_file = os.path.join(cwd, "libfree/"+os.path.splitext(basename)[0]+"_Q1.mgf")
        mgf2_file = os.path.join(cwd, "libfree/"+os.path.splitext(basename)[0]+"_Q2.mgf")
        mgf3_file = os.path.join(cwd, "libfree/"+os.path.splitext(basename)[0]+"_Q3.mgf")

        if not (os.path.exists(mgf1_file) \
            and os.path.exists(mgf2_file) \
            and os.path.exists(mgf3_file)):

            free_gigs = psutil.virtual_memory().available / (1024.0 ** 3)
            memory_str = str(math.trunc(free_gigs))
            
            cmd = [
                "java",
                "-Xms"+memory_str+"g",
                "-Xmx"+memory_str+"g",
                "-jar", "/opt/dia-umpire/DIA_Umpire_SE.jar",
                

            ]
            
            params = [
                os.path.join(cwd, "libfree", filename),
                diaumpire_cfg
            ]

            cmd.extend(params)
            logline(log_fh, "Running pipeline command: "+" ".join(cmd))
            proc = subprocess.Popen(cmd, cwd=cwd, stdout=log_fh, stderr=log_fh)

            while True:
                try:
                    proc.wait(1)
                    break
                except subprocess.TimeoutExpired:
                    if event and event.kill_flag:
                        proc.kill()
                        return

            log_fh.flush()
            if proc.returncode != 0:
                raise NonZeroReturnValueException(proc.returncode, 'dia-umpire')

        else:
            print (mgf1_file + ", " + mgf2_file + "," + mgf3_file + " exists, skipping diaumpire")

        phase_2_steps += 1
        progress.update_n_of_m(phase_1_steps + phase_2_steps, n_steps)

    filenames = []
    for root, dirs, files in os.walk(os.path.join(cwd, "libfree")):
        for filename in files:
            if filename.endswith(".mgf"):
                filenames.append(filename)

    phase_3_n = len(filenames)
    phase_3_steps = 0
    n_steps = phase_1_steps + phase_2_steps + phase_3_n # estimate

    for i, filename in enumerate(filenames):
        basename = os.path.basename(filename)
        pseudospectra_mzXML_file = os.path.join(cwd, "libfree-pseudospectra/"+os.path.splitext(basename)[0]+".mzXML")
        if not os.path.exists(pseudospectra_mzXML_file):
        
            cmd = [
                "/opt/tpp/bin/msconvert",
                "libfree/"+filename,
                "--mzXML",
                "-o", os.path.join(cwd, "libfree-pseudospectra")
            ]

            logline(log_fh, "Running pipeline command: "+" ".join(cmd))
            proc = subprocess.Popen(cmd, cwd=cwd, stdout=log_fh, stderr=log_fh)

            while True:
                try:
                    proc.wait(1)
                    break
                except subprocess.TimeoutExpired:
                    if event and event.kill_flag:
                        proc.kill()
                        return

            log_fh.flush()
            if proc.returncode != 0:
                raise NonZeroReturnValueException(proc.returncode, 'msconvert')

        else:
            print (pseudospectra_mzXML_file + " already exists, skipping conversion.")

        phase_3_steps += 1
        progress.update_n_of_m(phase_1_steps + phase_2_steps + phase_3_steps, n_steps)

    progress.update_n_of_m(n_steps, n_steps)
    progress.ready()

    pseudospectrafiles = []
    for root, dirs, files in os.walk(os.path.join(cwd, "libfree-pseudospectra" )):
        for filename in files:
            if filename.endswith(".mzXML"):
                pseudospectrafiles.append("libfree-pseudospectra/"+filename)
    return pseudospectrafiles
    

# def runMSGFPlus(DDA_DB_filename, DDA_filenames, outputfilename, threads):

#     DDA_pep_xmls = []
    
#     for DDA_filename in DDA_filenames:
        
#         basename = os.path.basename(DDA_filename)
        
#         if os.path.exists(os.path.splitext(basename)[0]+".pepXML"):
#             print (os.path.split(basename)[0]+".pepXML exists, skipping MSGF+")
#             DDA_pep_xmls.append(os.path.splitext(basename)[0]+".pepXML")
            
#         else:
            
#             msgf_cmd = [
#                 "java",
#                 "-Xmx120g",
#                 "-jar", "/opt/msgfplus/MSGFPlus.jar",
#                 "-d", DDA_DB_filename,
#                 "-s", DDA_filename,
#                 "-t", "10ppm",
#                 "-thread", threads,
#                 "-m", "3",
#                 "-tda", "0",
#                 "-o", os.path.splitext(basename)[0]+".mzid"
#                 ]

            
#             logline(log_fh, "Running pipeline command: "+" ".join(msgf_cmd))
#             msgf_proc = subprocess.Popen(msgf_cmd, stdout=log_fh, stderr=log_fh)
#             msgf_proc.wait()
#             log_fh.flush()
#             if msgf_proc.returncode != 0:
#                 raise NonZeroReturnValueException(msgf_proc.returncode, 'MSGF+')

#             idconvert_cmd = [
#                 "/opt/tpp/bin/idconvert", os.path.splitext(basename)[0]+".mzid",
#                 "--pepXML"
#             ]

#             logline(log_fh, "Running pipeline command: "+" ".join(idconvert_cmd))
#             idconvert_proc = subprocess.Popen(idconvert_cmd, stdout=log_fh, stderr=log_fh)
#             idconvert_proc.wait()
#             log_fh.flush()
#             if idconvert_proc.returncode != 0:
#                 raise NonZeroReturnValueException(idconvert_proc.returncode, 'idconvert')

#             DDA_pep_xmls.append(os.path.splitext(basename)[0]+".pepXML")
           
#     if os.path.exists(outputfilename):
#         print (outputfilename + " file already exist, skipping merging MSGF files.")        
#     else:
#         xinteract_cmd = [
#             "/opt/tpp/bin/xinteract",
#             "-OARPd",
#             "-dDECOY_",
#             "-N"+outputfilename,
#             ]
#         xinteract_cmd.extend(DDA_pep_xmls)
#         logline(log_fh, "Running pipeline command: "+" ".join(xinteract_cmd))
#         xinteract_proc = subprocess.Popen(xinteract_cmd, stdout=log_fh, stderr=log_fh)
#         xinteract_proc.wait()
#         log_fh.flush()
#         if xinteract_proc.returncode != 0:
#             raise NonZeroReturnValueException(xinteract_proc.returncode, 'xinteract')
       
#     return
   
        
def runComet(\
    event, \
    log_fh, \
    progressData, \
    cwd,\
    comet_cfg, \
    DDA_filenames, \
    existing_pep_xmls, \
    outputfilename, \
    worktempdir):

    progress = Progress(progressData, "percentage-indicator")
    n_steps = len(DDA_filenames) + 1

    DDA_basenames = []
    DDA_pep_xmls = []
    DDA_basename_to_pep_xml = {}
    
    for DDA_filename in DDA_filenames:
        basename = os.path.basename(DDA_filename)
        slink = os.path.join(cwd, basename)
        if os.path.exists(slink):
            if not os.path.islink(slink):
                progress.fail({
                    "cmd": "link  " + DDA_filename + " " + slink,
                    "returncode": -1, 
                })
        else:
            os.symlink(DDA_filename, slink)

        DDA_basenames.append(basename)
        pepXML = os.path.splitext(basename)[0]+".pep.xml"
        DDA_pep_xmls.append(pepXML)
        DDA_basename_to_pep_xml[basename] = pepXML
        

    existing_pep_xmls_by_samplenames = []
    if existing_pep_xmls:
        for existing_pep_xml in existing_pep_xmls:
            existing_basename = os.path.splitext(existing_pep_xml)[0] # remove xml
            existing_basename = os.path.splitext(existing_basename)[0] # remove pep
            existing_pep_xmls_by_samplenames.append(existing_basename)

    for i, DDA_basename in enumerate(DDA_basenames):

        samplename = os.path.splitext(DDA_basename)[0]

        if samplename in existing_pep_xmls_by_samplenames:
            # print ("DEBUG: Comet skipping " + basename)
            None
        else:

            # DDA_pep_xml = DDA_basename_to_pep_xml[DDA_basename]

            cmd = [
                "/opt/comet/comet-ms",
                "-P"+comet_cfg,
                ]
            cmd.append(DDA_basename)
            logline(log_fh, "Running pipeline command: "+" ".join(cmd))
            comet_proc = subprocess.Popen(cmd, cwd=cwd, stdout=log_fh, stderr=log_fh)
            
            while True:
                try:
                    comet_proc.wait(1)
                    break
                except subprocess.TimeoutExpired:
                    if event and event.kill_flag:
                        comet_proc.kill()
                        return

            log_fh.flush()
            if comet_proc.returncode != 0:
                progress.fail({
                    "cmd": " ".join(cmd),
                    "returncode": comet_proc.returncode, 
                })
                raise NonZeroReturnValueException(comet_proc.returncode, 'Comet')

        progress.update_n_of_m(i+1, n_steps)


    xinteract_cmd = [
        "/opt/tpp/bin/xinteract",
        "-OARPd",
        "-dDECOY_",
        "-N" + outputfilename,
        ]
    xinteract_cmd.extend(DDA_pep_xmls)
    logline(log_fh, "Running pipeline command: "+" ".join(xinteract_cmd))
    xinteract_proc = subprocess.Popen(xinteract_cmd, cwd=cwd, stdout=log_fh, stderr=log_fh)
    while True:
        try:
            xinteract_proc.wait(1)
            break
        except subprocess.TimeoutExpired:
            if event and event.kill_flag:
                xinteract_proc.kill()
                return

    log_fh.flush()
    if xinteract_proc.returncode != 0:
        progress.fail({
            "cmd": " ".join(xinteract_cmd),
            "returncode": xinteract_proc.returncode, 
        })
        raise NonZeroReturnValueException(xinteract_proc.returncode, 'xinteract')

    progress.update_n_of_m(n_steps, n_steps)
    progress.ready()

    return



def runXTandem(\
    event, \
    log_fh, \
    progressData, \
    cwd, \
    xtandem_default_input_filename, \
    DDA_DB_filename, \
    DDA_filenames, \
    existing_pep_xmls, \
    outputfilename, \
    worktempdir, \
    delete_temp_files_flag):

    progress = Progress(progressData, "percentage-indicator")
    n_steps = len(DDA_filenames) + 1

    taxonomy_tmp_fd = \
                        tempfile.NamedTemporaryFile(dir=worktempdir, \
                                                    mode="r+", \
                                                    delete=delete_temp_files_flag)
    taxonomy_xml = '<?xml version="1.0"?>\n'
    taxonomy_xml += '<bioml label="x! taxon-to-file matching list">\n'
    taxonomy_xml += '<taxon label="DB">\n'
    taxonomy_xml += '<file format="peptide" URL="' + DDA_DB_filename + '" />\n'
    taxonomy_xml += '</taxon>\n'
    taxonomy_xml += '</bioml>\n'
    taxonomy_tmp_fd.write(taxonomy_xml)
    taxonomy_tmp_fd.flush()

    DDA_pep_xmls = {}
    for DDA_filename in DDA_filenames:
        basename = os.path.basename(DDA_filename)
        DDA_pep_xmls[DDA_filename]=os.path.splitext(basename)[0]+".tandem.pep.xml"

    tandem_outs = {}

    existing_pep_xmls_by_samplenames = []
    if existing_pep_xmls:
        for existing_pep_xml in existing_pep_xmls:
            existing_basename = os.path.splitext(existing_pep_xml)[0] # remove xml
            existing_basename = os.path.splitext(existing_basename)[0] # remove pep
            existing_basename = os.path.splitext(existing_basename)[0] # remove tandem
            existing_pep_xmls_by_samplenames.append(existing_basename)

    for i, DDA_filename in enumerate(DDA_filenames):

        samplename = os.path.splitext(os.path.basename(DDA_filename))[0]
        tandem_outs[DDA_filename] = samplename+".TANDEM.OUTPUT.xml"

        if samplename in existing_pep_xmls_by_samplenames:
            # print ("DEBUG: X!Tandem skipping " + samplename)
            None
        else:
            input_xml = '<?xml version="1.0"?>\n'
            input_xml += '<bioml>\n'
            input_xml += '<note type="input" label="list path, default parameters">' + xtandem_default_input_filename + '</note>\n'
            input_xml += '<note type="input" label="list path, taxonomy information">'+taxonomy_tmp_fd.name+'</note>\n'
            input_xml += '<note type="input" label="protein, taxon">DB</note>\n'
            input_xml += '<note type="input" label="spectrum, path">' + DDA_filename + '</note>\n'
            input_xml += '<note type="input" label="output, path">' + tandem_outs[DDA_filename] + '</note>\n'
            input_xml += '</bioml>'
            input_tmp_fd = \
                        tempfile.NamedTemporaryFile(dir=worktempdir, \
                                                    mode="r+", \
                                                    delete=delete_temp_files_flag)
            input_tmp_fd.write(input_xml)
            input_tmp_fd.flush()
            tandem_cmd = [
                "/opt/tandem/tandem",
                input_tmp_fd.name,
                tandem_outs[DDA_filename]
                ]

            logline(log_fh, "Running pipeline command: "+" ".join(tandem_cmd))
            tandem_proc = subprocess.Popen(tandem_cmd, cwd=cwd, stdout=log_fh, stderr=log_fh)
            while True:
                try:
                    tandem_proc.wait(1)
                    break
                except subprocess.TimeoutExpired:
                    if event and event.kill_flag:
                        tandem_proc.kill()
                        return

            log_fh.flush()
            if tandem_proc.returncode != 0:
                progress.fail({
                    "cmd": " ".join(tandem_cmd),
                    "returncode": tandem_proc.returncode, 
                })
                raise NonZeroReturnValueException(tandem_proc.returncode, 'Tandem')

            tandemxml_cmd = [
                "/opt/tpp/bin/Tandem2XML",
                tandem_outs[DDA_filename],
                DDA_pep_xmls[DDA_filename]
                ]
            logline(log_fh, "Running pipeline command: "+" ".join(tandemxml_cmd))
            tandemxml_proc = subprocess.Popen(tandemxml_cmd, cwd=cwd, stdout=log_fh, stderr=log_fh)
            
            while True:
                try:
                    tandemxml_proc.wait(1)
                    break
                except subprocess.TimeoutExpired:
                    if event and event.kill_flag:
                        tandemxml_proc.kill()
                        return

            
            log_fh.flush()
            if tandem_proc.returncode != 0:
                progress.fail({
                    "cmd": " ".join(tandemxml_cmd),
                    "returncode": tandemxml_proc.returncode, 
                })
                raise NonZeroReturnValueException(tandemxml_proc.returncode, 'Tandem2XML')

        progress.update_n_of_m(i+1, n_steps)

       
    xinteract_cmd = [
        "/opt/tpp/bin/xinteract",
        "-OARPd",
        "-dDECOY_",
        "-N"+outputfilename,
        ]
    xinteract_cmd.extend([DDA_pep_xmls[x] for x in DDA_pep_xmls])

    logline(log_fh, "Running pipeline command: "+" ".join(xinteract_cmd))
    xinteract_proc = subprocess.Popen(xinteract_cmd, cwd=cwd, stdout=log_fh, stderr=log_fh)

    while True:
        try:
            xinteract_proc.wait(1)
            break
        except subprocess.TimeoutExpired:
            if event and event.kill_flag:
                xinteract_proc.kill()
                return
    
    log_fh.flush()
    if xinteract_proc.returncode != 0:
        progress.fail({
            "cmd": " ".join(xinteract_cmd),
            "returncode": xinteract_proc.returncode, 
        })
        raise NonZeroReturnValueException(xinteract_proc.returncode, 'xinteract')

    progress.update_n_of_m(n_steps, n_steps)
    progress.ready()

    return



def combine_search_engine_results(\
                                  event, \
                                  log_fh, \
                                  progressData,\
                                  cwd,
                                  decoyPrefix, \
                                  outputFilename, \
                                  DDA_DB_filename, \
                                  threads, \
                                  pepXMLs):
    progress = Progress(progressData, "percentage-indicator")
    progress.update_n_of_m(0, 100)

    interprophetparser_cmd = [
        "/opt/tpp/bin/InterProphetParser",
        "DECOY="+decoyPrefix,
        "THREADS="+threads,
    ]
    interprophetparser_cmd.extend(pepXMLs)
    interprophetparser_cmd.append(outputFilename)

    logline(log_fh, "Running pipeline command: "+" ".join(interprophetparser_cmd))
    interprophetparser_proc = subprocess.Popen(interprophetparser_cmd, cwd=cwd, stdout=log_fh, stderr=log_fh)
    
    while True:
        try:
            interprophetparser_proc.wait(1)
            break
        except subprocess.TimeoutExpired:
            if event and event.kill_flag:
                interprophetparser_proc.kill()
                return


    log_fh.flush()
    if interprophetparser_proc.returncode != 0:
        progress.fail({
            "cmd": " ".join(interprophetparser_cmd),
            "returncode": interprophetparser_proc.returncode, 
        })
        raise NonZeroReturnValueException(interprophetparser_proc.returncode, 'InterProphetParser')

    progress.update_n_of_m(100, 100)
    progress.ready()
    return



def buildlib(event, \
             log_fh, \
             progressData,\
             cwd, \
             decoyPrefix, \
             DDA_DB_filename, \
             threads, \
             libmethod_pepXMLs, \
             libfreemethod_pepXMLs, \
             iRT_filename, \
             swaths_filename, \
             protFDR, \
             gFDR, \
             swaths_min, \
             swaths_max):


    method_pepXMLs = {"lib":libmethod_pepXMLs, "libfree":libfreemethod_pepXMLs}
    speclibs = []

    progress = Progress(progressData, "percentage-indicator")
    n_steps = len(method_pepXMLs) +1
    
    for i, method in enumerate(method_pepXMLs):

        if len(method_pepXMLs[method]) > 0:
            
            if len(method_pepXMLs[method]) > 1:
                pepXMLFile = method+"_iprofet.peps.xml"
                
                combine_search_engine_results(event,
                                                log_fh,
                                                None,
                                                cwd,
                                                "DECOY_", \
                                                pepXMLFile, \
                                                DDA_DB_filename, \
                                                threads, \
                                                method_pepXMLs[method])
                
                if event and event.kill_flag:
                    return
                    
            else:
                assert(len(method_pepXMLs[method]) == 1)
                pepXMLFile = method_pepXMLs[method][0]


            mayu_folder = method+"_mayu"

            os.makedirs(os.path.join(cwd, mayu_folder))

            mayu_cmd = [
                "/opt/tpp/bin/Mayu.pl",
                "-A", "../"+pepXMLFile,
                "-C", "../"+DDA_DB_filename,
                "-E", decoyPrefix,
                "-G", str(gFDR),
                "-H", "51",
                "-I", "2",
                "-P", "protFDR="+str(protFDR)+":t",
            ]

            logline(log_fh, "Running pipeline command: "+" ".join(mayu_cmd))
            mayu_proc = subprocess.Popen(mayu_cmd, stdout=log_fh, stderr=log_fh, cwd=os.path.join(cwd, mayu_folder))
        
            while True:
                try:
                    mayu_proc.wait(1)
                    break
                except subprocess.TimeoutExpired:
                    if event and event.kill_flag:
                        mayu_proc.kill()
                        return

            log_fh.flush()
            if mayu_proc.returncode != 0:
                progress.fail({
                    "cmd": " ".join(mayu_cmd),
                    "returncode": mayu_proc.returncode, 
                })
                raise NonZeroReturnValueException(mayu_proc.returncode, 'Mayu')


            shell_cmd = "cat *_psm_protFDR0*.csv |cut -f 5 -d ',' |tail -n+2 |sort -u |head -n1"
            logline(log_fh, "Running pipeline command: "+ shell_cmd)
            shell_proc = subprocess.Popen(shell_cmd, stdout=subprocess.PIPE, stderr=log_fh, cwd=os.path.join(cwd, mayu_folder), shell=True)
            shell_proc.wait()
            log_fh.flush()
            if shell_proc.returncode != 0:
                progress.fail({
                    "cmd": " ".join(shell_cmd),
                    "returncode": shell_proc.returncode, 
                })
                raise NonZeroReturnValueException(shell_proc.returncode, 'shell command')

            cutoff = shell_proc.stdout.readline().decode("utf8").rstrip('\n')

            speclibs.append("SpecLib_"+method+".splib")

            spectrast_cmd1 = [
                "/opt/tpp/bin/spectrast",
                "-cNSpecLib_"+method,
                "-cIHCD",
                "-cf", "\"Protein! ~ "+decoyPrefix+"\"",
                "-cP"+cutoff,
                "-c_IRT"+iRT_filename, 
                "-c_IRR", pepXMLFile 
            ]

            logline(log_fh, "Running pipeline command: "+" ".join(spectrast_cmd1))
            spectrast_cmd1_proc = subprocess.Popen(spectrast_cmd1, cwd=cwd, stdout=log_fh, stderr=log_fh)
            
            while True:
                try:
                    spectrast_cmd1_proc.wait(1)
                    break
                except subprocess.TimeoutExpired:
                    if event and event.kill_flag:
                        spectrast_cmd1_proc.kill()
                        return

            log_fh.flush()
            if spectrast_cmd1_proc.returncode != 0:
                progress.fail({
                    "cmd": " ".join(spectrast_cmd1),
                    "returncode": spectrast_cmd1_proc.returncode, 
                })
                raise NonZeroReturnValueException(spectrast_cmd1_proc.returncode, 'spectrast')

    progress.update_n_of_m(i+1, n_steps)


    # Merge lib and libree method libraries

    if len(speclibs) == 0:
        print("No spectral libraries. Library must be obtained by either DDA lib or libfree method.")
        sys.exit(1)
        
    elif len(speclibs) == 1:
        speclib = speclibs[0]
        
    elif len(speclibs) == 2:
    
        spectrast_concatlib_cm = [
            "/opt/tpp/bin/spectrast",
            "-cNSpecLib_merged",
            "-cJA"
        ]
        assert (len(speclibs)==2)
        spectrast_concatlib_cm.extend(speclibs)

        logline(log_fh, "Running pipeline command: "+" ".join(spectrast_concatlib_cm))
        spectrast_concatlib_cm_proc = subprocess.Popen(spectrast_concatlib_cm, cwd=cwd, stdout=log_fh, stderr=log_fh)

        while True:
            try:
                spectrast_concatlib_cm_proc.wait(1)
                break
            except subprocess.TimeoutExpired:
                if event and event.kill_flag:
                    spectrast_concatlib_cm_proc.kill()
                    return

        log_fh.flush()
        if spectrast_concatlib_cm_proc.returncode != 0:
            progress.fail({
                "cmd": " ".join(spectrast_concatlib_cm),
                "returncode": spectrast_concatlib_cm_proc.returncode, 
            })
            raise NonZeroReturnValueException(spectrast_concatlib_cm_proc.returncode, 'spectrast merge lib')

        speclib = "SpecLib_merged.splib"
        
    else:
        print ("Found " + str(len(speclibs)) + " spectral libraries. "+",".join(speclibs)+" The amount was different than expected." )
        sys.exit(1)
        
   # Build consensus lib 
    
    spectrast_cmd2 = [
        "/opt/tpp/bin/spectrast",
        "-cNSpecLib_cons",
        "-cD"+DDA_DB_filename,
        "-cIHCD",
        "-cAC", speclib,
    ]

    logline(log_fh, "Running pipeline command: "+" ".join(spectrast_cmd2))
    spectrast_cmd2_proc = subprocess.Popen(spectrast_cmd2, cwd=cwd, stdout=log_fh, stderr=log_fh)
    
    while True:
        try:
            spectrast_cmd2_proc.wait(1)
            break
        except subprocess.TimeoutExpired:
            if event and event.kill_flag:
                spectrast_cmd2_proc.kill()
                return

    log_fh.flush()
    if spectrast_cmd2_proc.returncode != 0:
        progress.fail({
            "cmd": " ".join(spectrast_cmd2),
            "returncode": spectrast_cmd2_proc.returncode, 
        })
        raise NonZeroReturnValueException(spectrast_cmd2_proc.returncode, 'spectrast')

    spectrast2tsv_cmd = [
        "spectrast2tsv.py",
        "-l", str(swaths_min)+","+str(swaths_max),
        "-s", "y,b",
        "-d",
        "-e",
        "-o", "6",
        "-n", "6", 
        "-w", swaths_filename, 
        "-k", "openswath",
        "-a", "SpecLib_cons_openswath.tsv",
        "SpecLib_cons.sptxt"
    ]    

    logline(log_fh, "Running pipeline command: "+" ".join(spectrast2tsv_cmd))
    spectrast2tsv_cmd_proc = subprocess.Popen(spectrast2tsv_cmd, cwd=cwd, stdout=log_fh, stderr=log_fh)
    
    while True:
        try:
            spectrast2tsv_cmd_proc.wait(1)
            break
        except subprocess.TimeoutExpired:
            if event and event.kill_flag:
                spectrast2tsv_cmd_proc.kill()
                return


    log_fh.flush()
    if spectrast2tsv_cmd_proc.returncode != 0:
        progress.fail({
            "cmd": " ".join(spectrast2tsv_cmd),
            "returncode": spectrast2tsv_cmd_proc.returncode, 
        })
        raise NonZeroReturnValueException(spectrast2tsv_cmd_proc.returncode, 'spectrast2tsv')


    ConvertTSVToTraML_cmd = [
        "TargetedFileConverter",
        "-in", "SpecLib_cons_openswath.tsv",
        "-out", "SpecLib_cons.TraML"

    ]

    logline(log_fh, "Running pipeline command: "+" ".join(ConvertTSVToTraML_cmd))
    ConvertTSVToTraML_proc = subprocess.Popen(ConvertTSVToTraML_cmd, cwd=cwd, stdout=log_fh, stderr=log_fh)
    
    while True:
        try:
            ConvertTSVToTraML_proc.wait(1)
            break
        except subprocess.TimeoutExpired:
            if event and event.kill_flag:
                ConvertTSVToTraML_proc.kill()
                return

    log_fh.flush()
    if ConvertTSVToTraML_proc.returncode != 0:
        progress.fail({
            "cmd": " ".join(ConvertTSVToTraML_cmd),
            "returncode": ConvertTSVToTraML_proc.returncode, 
        })
        raise NonZeroReturnValueException(ConvertTSVToTraML_proc.returncode, 'TargetedFileConverter')

    OpenSwathDecoyGenerator_cmd = [
        "OpenSwathDecoyGenerator",
        "-in", "SpecLib_cons.TraML",
        "-out", "SpecLib_cons_decoy.TraML",
        "-method", "shuffle",
#        "-append",
#        "-exclude_similar",
#        "-remove_unannotated"
        ]

    logline(log_fh, "Running pipeline command: "+" ".join(OpenSwathDecoyGenerator_cmd))
    OpenSwathDecoyGenerator_proc = subprocess.Popen(OpenSwathDecoyGenerator_cmd, cwd=cwd, stdout=log_fh, stderr=log_fh)
    
    while True:
        try:
            OpenSwathDecoyGenerator_proc.wait(1)
            break
        except subprocess.TimeoutExpired:
            if event and event.kill_flag:
                OpenSwathDecoyGenerator_proc.kill()
                return

    log_fh.flush()
    if OpenSwathDecoyGenerator_proc.returncode != 0:
        progress.fail({
            "cmd": " ".join(OpenSwathDecoyGenerator_cmd),
            "returncode": OpenSwathDecoyGenerator_proc.returncode, 
        })
        raise NonZeroReturnValueException(OpenSwathDecoyGenerator_proc.returncode, 'OpenSwathDecoyGenerator')    

    progress.update_n_of_m(n_steps, n_steps)
    progress.ready()

    return


def buildDIAMatrix(\
    event,\
    log_fh, \
    progressData,
    cwd,
    DIA_filenames, \
    fixed_swaths_filename, \
    target_FDR, \
    max_FDR, \
    threads, \
    irt_assay_library_traml,\
    design_file):    

    successfull_DIA_filenames = []

    progress = Progress(progressData, "percentage-indicator")
    n_steps = len(DIA_filenames) * 2 +1

    # REMOVE
    # for i, DIA_filename in enumerate(DIA_filenames):
    #     successfull_DIA_filenames.append(DIA_filename)
    #     progress.update_n_of_m(i+1, n_steps)

    for i, DIA_filename in enumerate(DIA_filenames):
       
        OpenSwathWorkflow_cmd = [
            "OpenSwathWorkflow",
            "-in", DIA_filename,
            "-tr", "SpecLib_cons_decoy.TraML", 
            "-tr_irt", irt_assay_library_traml, 
            "-out_tsv", os.path.basename(DIA_filename)+"-DIA.tsv", 
            "-min_upper_edge_dist", "1",
            "-sort_swath_maps",
            "-swath_windows_file", fixed_swaths_filename,
            "-force",
            "-threads", threads

        ]

        logline(log_fh, "Running pipeline command: "+" ".join(OpenSwathWorkflow_cmd))
        OpenSwathWorkflow_proc = subprocess.Popen(OpenSwathWorkflow_cmd, cwd=cwd, stdout=log_fh, stderr=log_fh)
        
        while True:
            try:
                OpenSwathWorkflow_proc.wait(1)
                break
            except subprocess.TimeoutExpired:
                if event and event.kill_flag:
                    OpenSwathWorkflow_proc.kill()
                    return

        log_fh.flush()
        if OpenSwathWorkflow_proc.returncode != 0:
            #raise NonZeroReturnValueException(OpenSwathWorkflow_proc.returncode, 'OpenSwathWorkflow')
            logline(log_fh, "DIA sample "+ DIA_filename +" failed. Skipping the sample.")
            
        else:
            successfull_DIA_filenames.append(DIA_filename)

        progress.update_n_of_m(i+1, n_steps)

    successfull_DIA_filenames_after_pyprophet = []

    j = i
    n_steps = len(DIA_filenames) + len(successfull_DIA_filenames) +1
    progress.update_n_of_m(j+1, n_steps)

    for i, DIA_filename in enumerate(successfull_DIA_filenames):

        pyprophet_cmd = shlex.split("pyprophet --delim=tab --export.mayu " + os.path.basename(DIA_filename)+"-DIA.tsv" + " --ignore.invalid_score_columns")
        
        logline(log_fh, "Running pipeline command: "+" ".join(pyprophet_cmd))
        curr_env = os.environ.copy()
        # print(curr_env)
        if "PYTHONPATH" in curr_env:
            del curr_env['PYTHONPATH']
        pyprophet_proc = subprocess.Popen(pyprophet_cmd, cwd=cwd, env=curr_env, stdout=log_fh, stderr=log_fh)

        while True:
            try:
                pyprophet_proc.wait(5)
                break
            except subprocess.TimeoutExpired:
                if event and event.kill_flag:
                    pyprophet_proc.kill()
                    return


        log_fh.flush()
        if pyprophet_proc.returncode != 0:
            logline(log_fh, "DIA sample "+ DIA_filename +" failed. Skipping the sample.")              #raise NonZeroReturnValueException(pyprophet_proc.returncode, 'pyprophet')
        else:
            successfull_DIA_filenames_after_pyprophet.append(DIA_filename)

        progress.update_n_of_m(i+j+1, n_steps)

    feature_alignment_cmd = [
        "feature_alignment.py",
        "--method", "best_overall",
        "--realign_method", "diRT",
        "--max_rt_diff", "90",
        "--target_fdr", str(target_FDR),
        "--max_fdr_quality", str(max_FDR),
        "--out", "DIA-analysis-result.csv", 
        "--in" 
    ]
    feature_alignment_cmd.extend([ os.path.basename(x) + "-DIA_with_dscore.csv" for x in successfull_DIA_filenames_after_pyprophet])

    logline(log_fh, "Running pipeline command: "+" ".join(feature_alignment_cmd))
    feature_alignment_proc = subprocess.Popen(feature_alignment_cmd, cwd=cwd, stdout=log_fh, stderr=log_fh)
    
    while True:
        try:
            feature_alignment_proc.wait(1)
            break
        except subprocess.TimeoutExpired:
            if event and event.kill_flag:
                feature_alignment_proc.kill()
                return

    log_fh.flush()
    if feature_alignment_proc.returncode != 0:
        progress.fail({
            "cmd": " ".join(feature_alignment_cmd),
            "returncode": feature_alignment_proc.returncode, 
        })
        raise NonZeroReturnValueException(feature_alignment_proc.returncode, 'feature_alignment')


    swaths2stats_cmd = [
        "/opt/gladiator/swaths2stats.R",
        "--input", "DIA-analysis-result.csv"
    ]

    if design_file:
        swaths2stats_cmd.extend(["--design-file", design_file])

    logline(log_fh, "Running pipeline command: "+" ".join(swaths2stats_cmd))
    swaths2stats_proc = subprocess.Popen(swaths2stats_cmd, cwd=cwd, stdout=log_fh, stderr=log_fh)
    
    while True:
        try:
            swaths2stats_proc.wait(1)
            break
        except subprocess.TimeoutExpired:
            if event and event.kill_flag:
                swaths2stats_proc.kill()
                return

    
    log_fh.flush()
    if swaths2stats_proc.returncode != 0:
        progress.fail({
            "cmd": " ".join(swaths2stats_cmd),
            "returncode": swaths2stats_proc.returncode, 
        })
        raise NonZeroReturnValueException(swaths2stats_proc.returncode, 'swaths2stats')

    progress.update_n_of_m(n_steps, n_steps)
    progress.ready()

    return


def read_swath_windows(dia_mzML):

    print ("DEBUG: reading_swath_windows: ", dia_mzML)
    
    context = ET.iterparse(dia_mzML, events=("start", "end"))

    windows = {}
    for event, elem in context:

        if event == "end" and elem.tag == '{http://psi.hupo.org/ms/mzml}precursor':
            il_target = None
            il_lower = None
            il_upper = None

            isolationwindow = elem.find('{http://psi.hupo.org/ms/mzml}isolationWindow')
            for cvParam in isolationwindow.findall('{http://psi.hupo.org/ms/mzml}cvParam'):
                name = cvParam.get('name')
                value = cvParam.get('value')

                if (name == 'isolation window target m/z'):
                    il_target = value
                elif (name == 'isolation window lower offset'):
                    il_lower = value
                elif (name == 'isolation window upper offset'):
                    il_upper = value

            ionList = elem.find('{http://psi.hupo.org/ms/mzml}selectedIonList')
           
            selectedion = ionList.find('{http://psi.hupo.org/ms/mzml}selectedIon')

            if selectedion:
            
                for cvParam in selectedion.findall('{http://psi.hupo.org/ms/mzml}cvParam'):
                    name = cvParam.get('name')
                    value = cvParam.get('value')

                    if (name == 'selected ion m/z'):
                        if not il_target:
                            il_target = value
                
            if not il_target in windows:
                windows[il_target] = (il_lower, il_upper)
            else:
                lower, upper = windows[il_target]

#                print ("Target: " + str(il_target))
#                print ("Lower:" + str(lower))
#                print ("Upper:" + str(upper))
#                print ("IL Lower:" + str(il_lower))
#                print ("IL Upper:" + str(il_upper))
                
                assert (il_lower == lower)
                assert (il_upper == upper)
                return windows

    return windows


def create_swath_window_files(cwd, dia_mzML):

    windows = read_swath_windows(dia_mzML)

    swaths = []
    for x in windows:
        target_str = x
        lower_str, upper_str = windows[x]
        target = float(target_str)
        lower = float(lower_str)
        upper = float(upper_str)
        assert (lower > 0)
        assert (upper > 0)
        swaths.append((target - lower, target + upper))
        
    swaths.sort(key=lambda tup: tup[0])

    tswaths = []
    tswaths.append(swaths[0])
    for i in range(1, len(swaths)):
        if swaths[i-1][1] > swaths[i][0]:
            lower_prev, upper_prev = swaths[i-1]
            lower, upper = swaths[i]
            assert (upper_prev < upper)
            tswaths.append((upper_prev, upper))
        else:
            tswaths.append(swaths[i])

    assert (len(swaths) == len(tswaths))
            
    with open(os.path.join(cwd, "swath-windows.txt"), "w") as fh_swaths, open(os.path.join(cwd, "truncated-swath-windows.txt"), "w") as fh_tswaths:
        fh_tswaths.write("LowerOffset\tHigherOffset\n")

        for i in range(len(swaths)):
            fh_swaths.write(str(swaths[i][0]) + "\t" + str(swaths[i][1])  + "\n")
            fh_tswaths.write(str(tswaths[i][0]) + "\t" + str(tswaths[i][1])  + "\n")

    return swaths, tswaths
            
def build_database(\
    event, \
    log_fh, \
    progressData, \
    cwd, \
    fasta_filenames, \
    database_filename, \
    database_decoy_filename):

    progress = Progress(progressData, "percentage-indicator")

    # Combine sequence filenames

    IDs = set()
    seqRecords = []
    n_steps = len(fasta_filenames) +2 # + writing output + making decoydb

    for i, filename in enumerate(fasta_filenames):
        records = SeqIO.index(filename, "fasta")

        for ID in records:
            if ID not in IDs:
                seqRecords.append(records[ID])
                IDs.add(ID)
            else:
                print("Found duplicated sequence ID " + str(ID) + ", skipping this sequence from file " + filename)

        progress.update_n_of_m(i+1, n_steps)

    SeqIO.write(seqRecords, os.path.join(cwd, database_filename), "fasta")

    progress.update_n_of_m(i+2, n_steps)

    # Build Decoy database

    cmd = [
            "DecoyDatabase",
            "-in", database_filename,
            "-out", database_decoy_filename,
            ]
    logline(log_fh, "Running pipeline command: "+" ".join(cmd))
    proc = subprocess.Popen(cmd, cwd=cwd, stdout=log_fh, stderr=log_fh)
    
    while True:
        try:
            proc.wait(1)
            break
        except subprocess.TimeoutExpired:
            if event and event.kill_flag:
                proc.kill()
                return

    log_fh.flush()
    if proc.returncode != 0:
        progress.fail({
                "cmd": " ".join(cmd),
                "returncode": proc.returncode, 
            }
        )
        raise NonZeroReturnValueException(proc.returncode, 'DecoyDatabase')

    progress.update_n_of_m(n_steps, n_steps)
    progress.ready()

    return















def combine_DB(progressData, fasta_filenames, combined_fasta_filename):
    print ("OBSOLETE function: combine_DB")

    progress = Progress(progressData, "percentage-indicator")

    IDs = set()
    seqRecords = []
    n_steps = len(fasta_filenames) +1 # inpute files + writing output

    for i, filename in enumerate(fasta_filenames):
        records = SeqIO.index(filename, "fasta")
        # n_records = len(records)

        for ID in records:
            if ID not in IDs:
                seqRecords.append(records[ID])
                IDs.add(ID)
            else:
                print("Found duplicated sequence ID " + str(ID) + ", skipping this sequence from file " + filename)

        progress.update_n_of_m(i+1, n_steps)

    SeqIO.write(seqRecords, combined_fasta_filename, "fasta")

    progress.update_n_of_m(n_steps, n_steps)

    progress.ready()
    return
    

def scan_project_phases(projectname, result_root):
    
    phases = {}
    phases["projectname"] = projectname

    dir = os.path.join(result_root, projectname)

    phases["folder"] = False
    if not (os.path.isdir(result_root) and os.path.exists(dir)):
        return phases
    phases["folder"] = True

    cfg_file = os.path.join(dir, "config.txt")
    if not os.path.isfile(cfg_file):
        return phases

    cfg = None

    with open(cfg_file) as fh:
        cfg = json.load(fh)
        phases["config"] = cfg


    if os.path.isfile(os.path.join(dir, "DB_with_decoys.fasta")):
        phases["DB"] = True
    else:
        phases["DB"] = False

    samples = cfg["files"]["samples"]
    samples_wo_extension = [ os.path.splitext(os.path.basename(x))[0] for x in samples]

    library_samples = cfg["files"]["library"]
    library_samples_wo_extension = [ os.path.splitext(os.path.basename(x))[0] for x in library_samples]

    dda_path = os.path.join(dir, "DDA")
    if os.path.isdir(dda_path):
        phases["conv_DDA"] = True
        for basename in library_samples_wo_extension:
            if not os.path.isfile(os.path.join(dda_path, basename + ".mzXML")):
                phases["conv_DDA"] = False

    dia_path = os.path.join(dir, "DIA")
    if os.path.isdir(dia_path):
        phases["conv_DIA"] = True
        for basename in samples_wo_extension:
            if not os.path.isfile(os.path.join(dia_path, basename + ".mzML")):
                phases["conv_DIA"] = False

    if "use_speudospectra_flag" in cfg["options"]:

        phases["pseudospectrafiles"] = []
        pseudospectra_path = "libfree-pseudospectra"
        if os.path.isdir(os.path.join(dir, pseudospectra_path) ):

            phases["pseudospectra"] = True
            for basename in samples_wo_extension:
                for stage in ["Q1", "Q2", "Q3"]:
                    rel_filepath = os.path.join(pseudospectra_path, basename + "_" + stage + ".mzXML")
                    if os.path.isfile(os.path.join(dir, rel_filepath)):
                        phases["pseudospectrafiles"].append(rel_filepath)
                    else:
                        phases["pseudospectra"] = False


    if "use_comet_flag" in cfg["options"]:
        
        engine = "comet"

        samples_missed = []
        pseudo_samples_missed = []
        samples_existing = []
        pseudo_samples_existing = []

        for basename in library_samples_wo_extension:

            xml = os.path.join(dir, basename + ".pep.xml")
            if not os.path.isfile(xml):
                samples_missed.append(basename + ".pep.xml")
            else:
                samples_existing.append(basename + ".pep.xml")

        if "use_speudospectra_flag" in cfg["options"]:

            for basename in samples_wo_extension:

                for stage in ["Q1", "Q2", "Q3"]:
                    pseudo_basename = basename + "_" + stage
                    speudo_xml = os.path.join(dir, pseudo_basename + ".pep.xml")
                    if not os.path.isfile(speudo_xml):
                        pseudo_samples_missed.append(pseudo_basename + ".pep.xml")
                    else:
                        pseudo_samples_existing.append(pseudo_basename + ".pep.xml")

        if samples_missed:
            phases["samples_missed_by_" + engine] = samples_missed

        if samples_existing:
            phases["samples_generated_by_" + engine] = samples_existing

        if pseudo_samples_missed:
            phases["pseudo_samples_missed_by_" + engine] = pseudo_samples_missed
        
        if pseudo_samples_existing:
            phases["pseudo_samples_generated_by_" + engine] = pseudo_samples_existing

        pep_xml = os.path.join(dir, "interact_"+engine+"_pep.xml")
        phases[engine+"_peptides"] = os.path.isfile(pep_xml)

        pseudo_pep_xml = os.path.join(dir, "interact_"+engine+"_pseudo_pep.xml")
        phases[engine+"_pseudo_peptides"] = os.path.isfile(pseudo_pep_xml)


    if "use_xtandem_flag" in cfg["options"]:
        
        engine = "xtandem"

        samples_missed = []
        pseudo_samples_missed = []
        samples_existing = []
        pseudo_samples_existing = []

        for basename in library_samples_wo_extension:

            xml = os.path.join(dir, basename + ".tandem.pep.xml")
            if not os.path.isfile(xml):
                samples_missed.append(basename + ".tandem.pep.xml")
            else:
                samples_existing.append(basename + ".tandem.pep.xml")

        if "use_speudospectra_flag" in cfg["options"]:

            for basename in samples_wo_extension:

                for stage in ["Q1", "Q2", "Q3"]:
                    pseudo_basename = basename + "_" + stage
                    speudo_xml = os.path.join(dir, pseudo_basename + ".tandem.pep.xml")
                    if not os.path.isfile(speudo_xml):
                        pseudo_samples_missed.append(pseudo_basename + ".tandem.pep.xml")
                    else:
                        pseudo_samples_existing.append(pseudo_basename + ".tandem.pep.xml")

        if samples_missed:
            phases["samples_missed_by_" + engine] = samples_missed

        if samples_existing:
            phases["samples_generated_by_" + engine] = samples_existing

        if pseudo_samples_missed:
            phases["pseudo_samples_missed_by_" + engine] = pseudo_samples_missed
        
        if pseudo_samples_existing:
            phases["pseudo_samples_generated_by_" + engine] = pseudo_samples_existing
        
        pep_xml = os.path.join(dir, "interact_"+engine+"_pep.xml")
        phases[engine+"_peptides"] = os.path.isfile(pep_xml)

        pseudo_pep_xml = os.path.join(dir, "interact_"+engine+"_pseudo_pep.xml")
        phases[engine+"_pseudo_peptides"] = os.path.isfile(pseudo_pep_xml)

        if os.path.isfile(os.path.join(dir, "SpecLib_cons_decoy.TraML")):
            phases["speclib_cons"] = True
        else:
            phases["speclib_cons"] = False
        
        if os.path.isfile(os.path.join(dir, "DIA-analysis-result.csv")) \
            and os.path.isfile(os.path.join(dir, "DIA-peptide-matrix.tsv")) \
            and os.path.isfile(os.path.join(dir, "DIA-protein-matrix.tsv")):
            phases["matrices"] = True
        else:
            phases["matrices"] = False

        annotation_config_found = False
        if "annotation_filename" in cfg and cfg["annotation_filename"]:
            annotation_config_found = True

        if annotation_config_found and os.path.isfile(os.path.join(dir, "matrix.tsv")) \
            and os.path.isfile(os.path.join(dir, "annotations.tsv")):
            phases["annotations"] = True
        else:
            phases["annotations"] = False

        # TODO: check if actual plot files also exist in the folder
        if annotation_config_found and os.path.isdir(os.path.join(dir, "figures")):
            phases["figures"] = True
        else:
            phases["figures"] = False

    return phases


if __name__=="__main__":
    sys.exit(0)



