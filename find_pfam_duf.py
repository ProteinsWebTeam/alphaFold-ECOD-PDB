import argparse
import os
from configparser import ConfigParser

import mysql.connector
import requests


def db_connection(conf):
    mydb = mysql.connector.connect(
        host=conf["host"],
        user=conf["user"],
        password=conf["password"],
        port=conf["port"],
        database=conf["db"],
    )
    mycursor = mydb.cursor()
    return mydb, mycursor


def get_duf(cur):
    list_duf = dict()
    request = """select pfamA_acc, pfamA_id, description 
                from pfamA where pfamA_id like 'DUF%'
            """
    cur.execute(request)
    for row in cur.fetchall():
        list_duf[row[0]] = {"id": row[1], "desc": row[2]}
    return list_duf


def get_pfam_AF(conf):
    myslcon, mycursor = db_connection(conf)
    list_pfam_af = dict()
    request = """select distinct pfamA_acc, pfamseq_acc, seq_start, seq_end
                from AF2 
                order by pfamseq_acc, pfamA_acc, seq_start
            """
    mycursor.execute(request)
    for row in mycursor.fetchall():
        try:
            list_pfam_af[row[0]][row[1]] += f";{row[2]}-{row[3]}"
        except KeyError:
            list_pfam_af[row[0]] = {row[1]: f"{row[2]}-{row[3]}"}

    mycursor.close()
    myslcon.close()

    print(len(list_pfam_af))
    return list_pfam_af


def get_af_file(outputdir, protein):
    url = f"https://alphafold.ebi.ac.uk/files/AF-{protein}-F1-model_v2.pdb"
    filepath = os.path.join(outputdir, f"AF-{protein}-F1-model_v2.pdb")

    if not os.path.isfile(filepath) or os.path.getsize(filepath) == 0:
        r = requests.get(url)
        with open(filepath, "wb") as f:
            f.write(r.content)


def get_pfam_with_pdb(conf):
    myslcon, mycursor = db_connection(conf)
    pfam_with_struct = list()
    request = """select distinct pfamA_acc
                from pdb_pfamA_reg
        """
    mycursor.execute(request)

    for row in mycursor.fetchall():
        pfam_with_struct.append(row[0])

    mycursor.close()
    myslcon.close()

    print(len(pfam_with_struct))
    return pfam_with_struct


def get_pfam(conf, pf_struct):
    myslcon, mycursor = db_connection(conf)
    pfam_dict = dict()

    request = """select pfamA.pfamA_acc, pfamA_id, description, clan_acc
                from pfamA 
                left join clan_membership on clan_membership.pfamA_acc = pfamA.pfamA_acc
                """
    mycursor.execute(request)

    for row in mycursor.fetchall():
        pfam_acc = row[0]

        cl = "NO_CLAN"
        if row[3]:
            cl = row[3]

        struct = "NO_STRUCT"
        if pfam_acc in pf_struct:
            struct = "HAS_STRUCT"

        pfam_dict[pfam_acc] = {"id": row[1], "desc": row[2], "clan": cl, "struct": struct}

    mycursor.close()
    myslcon.close()

    print(len(pfam_dict))
    return pfam_dict


def write_in_file(list_pfam_af, pf_dict):
    af2process = dict()
    with open(af2pfam_file, "w") as af2pf:
        af2pf.write("af_prot\tpfam_acc\tpfam_id\tclan\tpfam_desc\thas_struct\taf_start\taf_end\n")
        for pfam, af_dict in list_pfam_af.items():
            count = 0
            for af, match in af_dict.items():
                # only get file if pfam doesn't have a known structure for first AF encountered
                if count == 0 and pf_dict[pfam]["struct"] == "NO_STRUCT":
                    get_af_file(outputdir, af)
                    try:
                        af2process[af].add(pfam)
                    except KeyError:
                        af2process[af] = set(pfam)
                count += 1
                af2pf.write(
                    f"{af}\t{pfam}\t{pf_dict[pfam]['id']}\t{pf_dict[pfam]['clan']}\t{pf_dict[pfam]['desc']}\t{pf_dict[pfam]['struct']}\t{match}\n"
                )

    print(len(af2process))
    return af2process


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("config", metavar="FILE", help="configuration file")
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"Cannot open '{args.config}': " f"no such file or directory")

    config = ConfigParser()
    config.read(args.config)
    outputdir = config["misc"]["afdir"]
    af2pfam_file = config["misc"]["af2pfam_file"]
    os.makedirs(outputdir, exist_ok=True)

    print("Searching Pfam/AlphaFold matches")
    list_pfam_af = get_pfam_AF(config["pfam_live"])

    print("Searching Pfam with PDB structures")
    pfam_struct = get_pfam_with_pdb(config["pfam_rel"])

    print("Searching pfam info")
    pfam_dict = get_pfam(config["pfam_rel"], pfam_struct)

    print("Searching for Pfam with AlphaFold models")
    list_af_to_process = write_in_file(list_pfam_af, pfam_dict)
