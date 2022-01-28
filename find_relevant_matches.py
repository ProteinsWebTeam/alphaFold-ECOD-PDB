import subprocess
import os
import re
import argparse
from configparser import ConfigParser


def run(command_list):
    stdout = ''
    out_log = subprocess.run(command_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT )
    stdout = str(out_log.stdout, 'utf-8')
    # print (stdout)
    return stdout      

def run_tm_align(afdir, af, ecoddir, ecodpdb):
    af_file = os.path.join(afdir, f"{af}.pdb")
    pdb_file = os.path.join(ecoddir, ecodpdb)

    cmd = ["./TMalign", af_file, pdb_file]
    result = run(cmd)

    tm_align=""
    score = 0.0
    pattern=r'^TM-score'
    for line in result.split('\n'):
        if re.match(pattern, line):
            tm_align+=f"{line}\n"
            if "Chain_2" in line:
                score = float(line.split()[1])

    return tm_align, score

def get_best_match(list_matches):

    keep=dict()
    pattern=r'[0-9]+\.pdbnum\.pdb_[A-Z0-9]+'
    for line in list_matches:
        if re.search(pattern, line[1]): #non existing pdb file => ignore match
            continue
        qstart=int(line[7])
        qend=int(line[8])
        midpoint=(qend+qstart)/2

        if not keep:
            keep[f"{qstart}-{qend}"]=[line]
        else:
            new=True
            for key in keep.keys():
                saved_s,saved_e=key.split('-')
                if midpoint>int(saved_s)-10 and midpoint<int(saved_e)+10: #same domain
                    keep[key].append(line)
                    new=False
            if new:
                keep[f"{qstart}-{qend}"]=[line]
                new=False
    
    final_doms=dict()
    for key,lines in keep.items():
        # print(key)
        keepline = ""
        keepident=0.0
        maxi=5 if len(lines)>5 else len(lines)
        for i in range(0, maxi):
            evalue=float(lines[i][11])
            fident=float(lines[i][4])

            # print(fident, keepident)
            if evalue < 1.0e-04 and fident > keepident:
                    keepident = fident
                    keepline = lines[i]
        if keepline:
            final_doms[key]=keepline
    
    # print(final_doms)
    return final_doms

def get_extra_info(final_doms, afdir, af, ecoddir, ecod_file):   
    final_text=[]
    for info in final_doms.values():
        tm_align_score, score = run_tm_align(afdir, af, ecoddir, info[1])
        if score > 0.5:
            final_text.append(af)
            dom_s=int(info[7])
            dom_e=int(info[8])
            final_text.append(f"Domain: {dom_s}-{dom_e}")
            final_text.append(get_pfam_info(af2pfam, af, dom_s, dom_e))

            text="\t".join(info)
            final_text.append("# Most significant PDB match:")
            final_text.append("#query   target  qlen    tlen    fident  alnlen  mismatch    qstart  qend tstart tend    evalue  bits")
            final_text.append(text)

            ecod_dom=info[1].split(".")[0]
            cmdgrep=["grep", ecod_dom, ecod_file]
            output=run(cmdgrep)

            final_text.append("# ECOD domain information:")
            final_text.append("#uid	ecod_domain_id	manual_rep	f_id	pdb	chain	pdb_range	seqid_range	unp_acc	arch_name	x_name	h_name	t_name	f_name	asm_status	ligand")
            final_text.append(output.split("\n")[0])

            final_text.append("# TMalign scores")
            final_text.append(tm_align_score)
        elif score == 0.0:
            print(af, info)

    if final_text:
        final_text.append("######################\n")

    return final_text

def get_pfam_info(af2pfam, af, af_start, af_end):
    count = 0
    text=""
    pattern=r'AF-([A-Z0-9]+)-.*'
    found = re.search(pattern, af)
    af_prot=found.group(1)
    midpoint1 = (af_start+af_end)/2

    for item in af2pfam[af_prot]:
        pf_start = int(item[-2])
        pf_end = int(item[-1])
        midpoint2 = (pf_start+pf_end)/2
        count2 = 0
        
        if midpoint2 > af_start and midpoint2 < af_end:
            count2+=1
        elif midpoint1 > pf_start and midpoint1 < pf_end:
            count2+=1

        if count2:
            count+=1
            if count > 1:
                text+="\n"
            text+="\t".join(item)
            
    final_text=""          
    if count:
        final_text+="# Pfam entry(ies):\n"
        final_text+="#pfam_acc  pfam_id    clan   pfam_desc  struct_info   af_start    af_end\n"
        final_text+=text
    else:
        final_text+="# No Pfam without PDB structure found for this domain"
    
    return final_text


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("config", metavar="FILE", help="configuration file")
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"Cannot open '{args.config}': " f"no such file or directory")

    config = ConfigParser()
    config.read(args.config)

    afdir=config["misc"]["afdir"]
    ecoddir=config["misc"]["ecoddir"]
    af2pfam_file=config["misc"]["af2pfam_file"]
    ecod_file=config["misc"]["ecod_file"]
    matches=config["misc"]["output_file"]
    output_file=config["misc"]["matches_file"]

    af2pfam=dict()
    if os.path.isfile(af2pfam_file) and os.path.getsize(af2pfam_file) > 0:
        with open(af2pfam_file, "r") as f:
            for line in f.readlines():
                line=line.strip()
                values=line.split(",")
                af=values[0]
                if af not in af2pfam:
                    af2pfam[af]=[values[1:]]
                else:
                    af2pfam[af].append(values[1:])

    if os.path.isfile(matches) and os.path.getsize(matches) > 0:
        with open(matches, "r") as f, open(output_file, "w") as fout:
            af=""
            list_matches=[]
            
            for line in f.readlines():
                line = line.strip()
                values = line.split()
                aftmp=values[0][:-4]
                # print(aftmp)
                
                if not af:
                    af = aftmp
                elif af and aftmp!=af:
                    doms = get_best_match(list_matches)
                    text = []
                    if doms:
                        text = get_extra_info(doms, afdir, af, ecoddir, ecod_file)
                    # else:
                    #     text.append(af)
                    #     text.append("# No significant PDB match found\n")
                    #     text.append("######################\n")
                    fout.write("\n".join(text))

                    af=aftmp
                    list_matches=[]
                    
                    
                list_matches.append(values)

            doms = get_best_match(list_matches)
            text = []
            if doms:
                text = get_extra_info(doms, afdir, af, ecoddir, ecod_file)
            # else:
            #     text.append(af)
            #     text.append("# No significant PDB match found\n")
            #     text.append("######################\n")
            fout.write("\n".join(text))