#!/bin/bash
if [ ! -z $1 ]; then
    config_file=$1
    . $config_file
else
    echo "Missing config file"
    exit
fi

CWD=$(pwd)

export PATH=$(pwd)/bin/:$PATH
if [[ ! -f $ecod_db ]] || [[ ! -s $ecod_db ]];then
    echo "$ecod_db is empty"
    if [[ ! -d $ecoddir ]] || [ ! "$(ls -A $ecoddir)" ];then
        echo "Downloading PDB matches from ECOD"
        wget -O $ecod_pdb_file http://prodata.swmed.edu/ecod/distributions/ecod.latest.F70.pdb.tar.gz
        cd $basedir
	    tar xvzf $ecod_pdb_file
        mkdir $ecoddir
        cd $ecodtmpdir
	# cp -fr $ecodtmpdir $ecoddir
        for dir in ./*; do
            cd "$dir"
            for dir2 in *; do
                echo "$dir2"
                cd "$dir2"
                for file in *; do
                #echo "moving $file" 
                mv $file $ecoddir/    
                done
                cd ..
            done
            cd ..
            done
        cd $CWD
        rm -r "$basedir/data"
    fi
    echo $pwd
    echo $PATH
    echo "Generating ecodDB for foldseek"
    foldseek createdb $ecoddir/ $ecod_db
fi

export PFAM_CONFIG=$pfam_config

python3 find_pfam_duf.py $config_file || exit 2

echo "Running foldseek"
foldseek easy-search --format-output "query,target,qlen,tlen,fident,alnlen,mismatch,qstart,qend,tstart,tend,evalue,bits" $afdir $ecod_db $output_file $tmpdir

if [[ ! -f $ecod_file ]];then
    echo "Download $ecod_file"
    wget -O $ecod_file http://prodata.swmed.edu/ecod/distributions/ecod.latest.domains.txt
fi

echo "Search relevant matches"
python3 find_relevant_matches.py $config_file || exit 2
