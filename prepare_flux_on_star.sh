#!/bin/bash

##################
# OPTION PARSING
##################

#getopts w:d:r: VARNAME
#
#while getopts ":a:" opt; do
#  case $opt in
#    a)
#      echo "-a was triggered, Parameter: $OPTARG" >&2
#      ;;
#    \?)
#      echo "Invalid option: -$OPTARG" >&2
#      exit 1
#      ;;
#    :)
#      echo "Option -$OPTARG requires an argument." >&2
#      exit 1
#      ;;
#  esac
#done


###########
# VARIABLES
###########

# dashboard selected lines
dashboard_selection=$1
# the directory where to store all the outputs
working_dir='/users/rg/abreschi/Documents/blueprint/pilot/Flux/antisense'
# the Flux to use
flux_exec="/users/rg/epalumbo/bin/flux-capacitor-1.2.4-SNAPSHOT/bin/flux-capacitor"
# annotation file
annotation_file="/users/rg/projects/encode/scaling_up/whole_genome/Gencode/version15/gencode.v15.annotation.gtf"
# cluster directory
#cluster_dir="/scratch/local"
# Flux read_strand
read_strand="MATE1_SENSE"
# Flux annotation mapping
annotation_mapping="PAIRED_STRANDED" 
# Flux memory
fmem=6G
# Flux memory on cluster
cluster_fmem=8G
# File with the nodes that have permissions for cluster_dir
#nodes_ok='/users/rg/abreschi/nodes.ok'

# #-------#
# | BEGIN |
# #-------#


# I want to create a folder with the name of my sample 
# Each file, instead, will have the LID as identifier
# I can always go back to the metadata from the dashboard selection or from
# the dashboard itself

# Read the file with the good nodes
#nodes=( $(cat $nodes_ok | tr '\n' ' ') )
#index=0
#n_nodes=${#nodes[*]}

# Create the folders into the working directory, if not already existing
samples=($(grep -oP "labExpId=[^;]+" $dashboard_selection | cut -d= -f2 | sort -u))
for sample in ${samples[@]}; do mkdir -p $working_dir/$sample; done  # -p suppresses error when folder exists

# Create a new file for launching the flux jobs
echo "" > launch_flux_on_star.sh
chmod u+x launch_flux_on_star.sh

cat $dashboard_selection | while IFS=$'\t' read file metadata; do

############################################
# ITERATE OVER THE SAMPLES I WANT TO PROCESS
############################################

# For each sample a PARAMETER FILE is created

#sample=$(echo $metadata | grep -oP "cell=[^;]+" | cut -d= -f2)
LID=$(echo $metadata | grep -ioP "labexpid=[^;]+" | cut -d= -f2)
flux_param_file=$working_dir/$LID/${LID}_flux_params.txt
JOB_ID=${LID}


# Remove previous standard error and standard output
rm -f $working_dir/$LID/${LID}.launch_flux_on_star.err
rm -f $working_dir/$LID/${LID}.launch_flux_on_star.out

# ------------------------------------------------
# Write the jobs for cluster and give x permission
# ------------------------------------------------

echo "
#!/bin/bash
. /etc/profile

#$ -e $working_dir/$LID/${LID}.launch_flux_on_star.err
#$ -o $working_dir/$LID/${LID}.launch_flux_on_star.out

#$ -S /bin/bash
#$ -q rg-el6,short
#$ -pe smp 4    #number of cpus
#$ -l virtual_free=$cluster_fmem
#$ -m e 
#$ -M ale.breschi@gmail.com

# create all the directories that the flux will need
# mkdir -p $cluster_dir/${JOB_ID}/TMP_FLUX
# mkdir -p $cluster_dir/${JOB_ID}/TMP_SORT

# Write parameter file
echo \"ANNOTATION_FILE \$TMPDIR/`basename $annotation_file`
MAPPING_FILE $file
STDOUT_FILE \$TMPDIR/${LID}.flux.out.gtf
TMP_DIR \$TMPDIR
READ_STRAND $read_strand
ANNOTATION_MAPPING $annotation_mapping
SORT_IN_RAM true
\" > $flux_param_file


# make sure previous files are deleted
rm -f $working_dir/$LID/${LID}.flux.out.gtf

# go in the directory where the output will be
cd \$TMPDIR

# copy the annotation and mapping file
cp $annotation_file \$TMPDIR

# export important variables
export FLUX_MEM=$fmem

echo 'Running on' \$HOSTNAME
# launch the Flux
date
echo 'running flux'
$flux_exec -p $flux_param_file

# copy the output in my folder
cp \$TMPDIR/${LID}.flux.out.gtf $working_dir/${LID}/

echo 'done'
date

cd $working_dir

" > $working_dir/$LID/${LID}.launch_flux_on_star.sh
chmod u+x $working_dir/$LID/${LID}.launch_flux_on_star.sh

# ---------------------------------
# Write a script to launch the jobs 
# ---------------------------------


echo "
# -N defines the name of the job 
qsub -N $LID $working_dir/$LID/${LID}.launch_flux_on_star.sh
" >> launch_flux_on_star.sh


done


