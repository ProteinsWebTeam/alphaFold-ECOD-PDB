#!/bin/bash
if [ ! -z $1 ]; then
    config_file=$1
    . $config_file
else
    echo "missing config file"
    exit
fi

python find_pfam_duf.py $config_file

export PATH=$(pwd)/bin/:$PATH
if [[ ! -f $ecod_db ]];then
    if [[ ! -d $ecoddir ]];then
        echo "Downloading PDB matches from ECOD"
        wget http://prodata.swmed.edu/ecod/distributions/ecod.latest.F70.pdb.tar.gz
        tar xvzf $ecod_pdb_file
        mkdir $ecoddir
        cp -fr $ecodtmpdir $ecoddir
        rm -r $ecodtmpdir
    fi
    echo "Generating ecodDB for foldseek"
    foldseek createdb $ecoddir/ $ecod_db
fi


echo "running foldseek"
foldseek easy-search --format-output "query,target,qlen,tlen,fident,alnlen,mismatch,qstart,qend,tstart,tend,evalue,bits" $afdir $ecod_db $output_file $tmpdir

if [[ ! -f $ecod_file ]];then
    echo "download $ecod_file"
    wget http://prodata.swmed.edu/ecod/distributions/ecod.latest.domains.txt
fi

echo "search relevant matches"
python find_relevant_matches.py $config_file > $matches_file