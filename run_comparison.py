from genericpath import exists
import subprocess
import sys
import os
import re
import argparse
from configparser import ConfigParser
from multiprocessing import Pool

class compare():
    def __init__(self, inputdir, outputdir):
        self.inputdir=inputdir
        self.outputdir=outputdir
        self.ecod_db="ecod_pdb"
        self.tmpdir="tmp"
        self.ecodfile="ecod.latest.domains.txt"

        os.makedirs(self.outputdir, exist_ok=True)
        os.makedirs(self.tmpdir, exist_ok=True)

    
    
    def process_model(self, model):
        sys.stdout.flush()
        pattern=r'(AF-[A-Z0-9]+-.*)\.cif\.gz'
        found = re.search(pattern,model)
        if found:
            print(found.group(1))
            fpath=os.path.join(self.inputdir,model)
            outputfile=f"{self.outputdir}/{found.group(1)}.m8"

            if not os.path.isfile(outputfile) or (os.path.isfile(outputfile) and os.path.getsize(outputfile) < 0):
                # print(outputfile)
                cmd = ["foldseek", "easy-search", "--format-output", "query,target,qlen,tlen,fident,alnlen,mismatch,evalue,bits", fpath, self.ecod_db, outputfile, self.tmpdir]
                self.run(cmd)

            if os.path.isfile(outputfile) and os.path.getsize(outputfile) > 0:
                with open (outputfile, "r") as fout:
                    count = 0
                    keepline = ""
                    keepaln=0
                    for line in fout.readlines():
                        if count <5:
                            line = line.strip()
                            count +=1
                            # print(line)
                            evalue=float(line.split()[7])
                            qlen=int(line.split()[2])
                            alnlen=int(line.split()[5])

                            aln=(qlen*100)/alnlen
                            if evalue < 1.0e-04 and aln >= 60:
                                if aln > keepaln:
                                    keepaln = aln
                                    keepline = line
                        else:
                            break
                            
                    if keepline:
                        print(keepline)
                        ecod_dom=keepline.split()[1].split(".")[0]
                        cmdgrep=["grep", ecod_dom, self.ecodfile]
                        output=self.run(cmdgrep)
                        print(output)

def run(command_list):
    stdout = ''
    out_log = subprocess.run(command_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT )
    stdout = str(out_log.stdout, 'utf-8')
    # print (stdout)
    return stdout          

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("config", metavar="FILE", help="configuration file")
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"Cannot open '{args.config}': " f"no such file or directory")

    config = ConfigParser()
    config.read(args.config)
    #foldseek easy-search alpha_fold_ecoli/AF-P0AAZ0-F1-model_v2.cif.gz ecod_pdb ecoli_ecod_stats.m8 tmp
    os.environ['PATH'] += os.pathsep + f"{os.getcwd()}/bin/"

    outputdir=config["misc"]["outputdir"]
    inputdir=config["misc"]["afdir"]
    ecod_db=config["misc"]["ecod_db"]
    tmpdir=config["misc"]["tmpdir"]
    ecodfile=config["misc"]["ecod_file"]
    outputfile=config["misc"]["output_file"]
    os.makedirs(outputdir, exist_ok=True)
    os.makedirs(tmpdir, exist_ok=True)

    cmd = ["foldseek", "easy-search", "--format-output", "query,target,qlen,tlen,fident,alnlen,mismatch,evalue,bits", inputdir, ecod_db, outputfile, tmpdir]
    run(cmd)

    # comp=compare(inputdir, outputdir)
    # files_to_process=os.listdir(inputdir)
    # files_to_process=['AF-A0A385XJ53-F1-model_v2.cif.gz', 'AF-U3PVA8-F1-model_v2.cif.gz', 'AF-U3PVA8-F1-model_v2.pdb.gz', 'AF-V9HVX0-F1-model_v2.cif.gz', 'AF-V9HVX0-F1-model_v2.pdb.gz']

    # with Pool(10) as p:
    #     results = p.map(comp.process_model, files_to_process)
    # for f in files_to_process:
    #     comp.process_model(f)