#!/bin/bash
# A bash script to run ZEPPI
# ----------------------------------------------------------------------
# Copyright (C) Haiqing Zhao
# Honig Group at Columbia University
# Last Update: 01/12/2023
# ----------------------------------------------------------------------


if [ "$#" -lt 3 ]; then
    printf "\nUsage: bash Run_ZEPPI.sh PPI_list.csv Output.csv -option
where:
    PPI_list  A csv file containing your query PPIs.
    -m  calculate ZEPPI on mutual information & conservation (recommended for heterodimers).
    -md calculate ZEPPI on mutual information & conservation and direct coupling analysis (recommended for homodimers)."
    printf "\nOutput files will be named after the input file with different endings.\n\n"
    exit 1
fi

DCA=false
case $3 in
    -m|-mc|--m|--mc|-default|--default)
    printf "\nZEPPI is computed based on mutual information & conservation.\n\n"
    ;;
    -md|--md)
    printf "\nZEPPI is computed based on mutual information & conservation and direct coupling analysis.\n\n"
    DCA=true
    ;;
    -*|--*)
    echo "Unknown option $3; Use default"
    printf "\nZEPPI is computed based on mutual information & conservation.\n\n"
    ;;
esac

inputfile=$1
outputfile=$2
IFS='/' read -r -a arrIN <<< "$inputfile"
name=${arrIN[-1]::-4}

# Configure path; change to your file path
ZEPPI_base=/Your_path/ZEPPI-main
Method_dir=$ZEPPI_base/Methods
Project_dir=$ZEPPI_base/Demo
Seqmap_dir=$Project_dir/Seqmap
ASA_dir=$Project_dir/ASA
IFC_dir=$Project_dir/IFC
MSA_dir=$Project_dir/MSA

# Setup for running with SGE 
#$ -l mem=16G

# Setup for running with SLURM
#SBATCH --mem=16G

# Calculate coevolution and conservation signals
# The needed input files are: 
# *.msa: the MSA files for each protein
# *.ifc_seq: the interface contact file for the PPI; indices must be based on full-length sequences as used in the MSA file
# *.asr_seq: the surface residue file for the two proteins where the first line is for protein 1 and the second line for protein 2; indices must be based on full-length sequences as usd in the MSA file

Metric_dir=$ZEPPI_base/Scratch/Metrics
cd $ZEPPI_base/Scratch
mkdir Metrics

sed 1d $inputfile | while IFS=$',' read -r f1 f2 f3 f4 f5 #f6 f7 f8 f9 f10 f11 f12
do
    pdb=`echo "$f1" | tr '[:upper:]' '[:lower:]'`
    chain1=$f2
    chain2=$f3
    uni1=$f4
    uni2=$f5
    echo "##### Running ZEPPI for" $pdb $chain1 $chain2 "#####"

    if [[ -s $MSA_dir/$uni1".msa" && -s $MSA_dir/$uni2".msa" ]]; then 
        if [[ -s $IFC_dir/$pdb"_"$chain1$chain2".ifc_seq" ]]; then 
            if [[ -s $ASA_dir/$pdb"_b1_"$chain1$chain2".asr_seq" ]]; then
                echo "## Calculating MI and Conservation score.."
                python $Method_dir/CalcPPI_MI_Con_Zscore.py $MSA_dir/$uni1".msa" $MSA_dir/$uni2".msa" $IFC_dir/$pdb"_"$chain1$chain2".ifc_seq" $ASA_dir/$pdb"_b1_"$chain1$chain2".asr_seq" $Metric_dir/$pdb"_"$chain1$chain2".miZ"
                if $DCA ; then
                    echo "## Calculating DCA score.."
                    python $Method_dir/CalcPPI_DCA_Zscore.py $MSA_dir/$uni1".msa" $MSA_dir/$uni2".msa" $IFC_dir/$pdb"_"$chain1$chain2".ifc_seq" $ASA_dir/$pdb"_b1_"$chain1$chain2".asr_seq" $Metric_dir/$pdb"_"$chain1$chain2".dca"
                fi
            else echo "Surface reisude file is missing."
            fi

        else echo "Interface residue file is missing."
        fi

    else echo "MSA files are missing."
    fi 

done

<<'##END4'
##END4

# Write ZEPPI outputs

printf "\nFinalizing ZEPPI based on the above calculated metrics.\n\n"

if $DCA; then
    python $Method_dir/ZEPPI_collectMIZ.py $inputfile $name"_MIZ.csv" $Metric_dir
    python $Method_dir/ZEPPI_collectDCA.py $inputfile $name"_DCA.csv" $Metric_dir
    python $Method_dir/ZEPPI_final.py $name"_MIZ.csv" $name"_DCA.csv" $outputfile
else
    python $Method_dir/ZEPPI_collectMIZ.py $inputfile $name"_MIZ.csv" $Metric_dir
    python $Method_dir/ZEPPI_final.py $name"_MIZ.csv" $outputfile
fi 


