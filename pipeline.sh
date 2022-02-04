#!/bin/bash
if [ ! -z $1 ]; then
    config_file=$1
    . $config_file
else
    echo "Missing config file"
    exit
fi

CWD=$(pwd)

python3 find_pfam_duf.py $config_file || exit $1

export PATH=$(pwd)/bin/:$PATH
if [[ ! -f $ecod_db ]];then
    if [[ ! -d $ecoddir ]];then
        echo "Downloading PDB matches from ECOD"
        wget http://prodata.swmed.edu/ecod/distributions/ecod.latest.F70.pdb.tar.gz
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
                mv $file $targetdir/    
                done
                cd ..
            done
            cd ..
            done
        cd $CWD
        rm -r "$basedir/data"
    fi
    echo "Generating ecodDB for foldseek"
    foldseek createdb $ecoddir/ $ecod_db
fi


echo "Running foldseek"
foldseek easy-search --format-output "query,target,qlen,tlen,fident,alnlen,mismatch,qstart,qend,tstart,tend,evalue,bits" $afdir $ecod_db $output_file $tmpdir

if [[ ! -f $ecod_file ]];then
    echo "Download $ecod_file"
    wget http://prodata.swmed.edu/ecod/distributions/ecod.latest.domains.txt
fi

echo "Search relevant matches"
python3 find_relevant_matches.py $config_file || exit $1