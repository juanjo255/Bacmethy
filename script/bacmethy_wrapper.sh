#!/bin/bash -login
## The original Bacmethyl was designed for one contig bacteria.
## This is a lazy wrapper on Bacmethy to process multi contig bacteria
## It will split fasta into its contigs run Bacmethy and concat outputs

dir=$(cd $(dirname $0);pwd) 
source $dir"/log_message_function.sh"

usage() {
    echo "
  Bacmethy V1.0.1 01/01/2022
  <liujihong2021@126.com>
  usage:  bacmethy.sh  -m <motifs.gff> -s <motifs.csv> -g <genome.fa> -p <prefix> -t <m6A, m4C or m5C>  [options] -T  -d -i -c -b -a -r -G
  requirement: locally installed PROKKA,  bedtools and MEME  softwares


  options
      -m    FILE    a motif detected file in gff format from PacBio portal (required)
      -s    FILE    a motifs summary file (required)
      -g      FILE    a genome file which complete sequenced (required)
      -p      STRING  prefix of output (required; usually a strains name or a sample name)
      -t      STRING  a type of methylation type (required; one of m6A, m4C or m5C)
      -T      FOLDER  a floder contains TFs files in meme format (required; users can get from Calcute_PPM_console_final.py)
      -d      FLOAT   a undermethylation thresholds of fraction (default: 0.75)
      -i      INT     a identificationQV thresholds (default: 40)
      -c      INT     a coverage thresholds (default: 30)
      -b      INT     number of bps before TSS (default: 500)
      -a      INT     number of bps after TSS (default: 100)
      -r      FILE    FASTA or GBK file to use as 1st priority (default '')
      -G      NULL    Scan the DNA methylation sites on gene Coding region (default only scan Regulation Region)
    "
}

while getopts ":hm:s:g:p:t:T:d:b:a:r:i:n:c:G" opt
do
    case $opt in
        h)  usage
            exit 0
            ;;
        m)    METHYLATIONGFF=$OPTARG ;;
        s)    METHYLATIONCSV=$OPTARG ;;
        g)  FASTA=$OPTARG ;;
        p)  sample=$OPTARG ;;
        t)  METHYLATIONTYPE=$OPTARG ;;
        T) MEME=$OPTARG ;;
        d)  undermethylation_threshold=$OPTARG ;;
        b)  promoter=$OPTARG ;;
        a)  cds=$OPTARG ;;
        r)  proteins=$OPTARG ;;
        i)  identification=$OPTARG ;;
        n)  npu=$OPTARG ;;
        c) coverage=$OPTARG ;;
	G) Gene="true" ;;
        :)
        echo "Option -$OPTARG requires an argument."
	    usage
        exit 1
        ;;
	    ?) 
	        echo "Invalid option: -$OPTARG"
	        usage      
            exit 1
	    ;;
    esac    
done

#If options of 'd', 'c'and 'a' are not defined, setting to default values
undermethylation_threshold=${undermethylation_threshold:-0.75}
promoter=${promoter:-500}
cds=${cds:-100}
identification=${identification:-40}
coverage=${coverage:-30}
npu=${npu:-5}

## Working Dir
wd=$PWD


### START PIPELINE ####
log_warning "This is just a wrapper every error is handle by the original Bacmethy pipeline"
log_warning "Use absolute paths!!"


## Check if FASTA has more than one contig
num_contigs=$(grep -c ">" $FASTA)
if [ "$num_contigs" -le 1 ]; then
  log_info "One Contig Fasta: Running Bacmethy"
  bacmethy.sh "$@"
  exit 1
fi

## Multi-fasta processing
log_info "Multi-Fasta: Splitting by contig before Running Bacmethy"
seqkit_out_dir=$wd"/fasta_split"
seqkit split --quiet --by-id --out-dir "$seqkit_out_dir" $FASTA

wrapper() {

  args=("$@")
  for contig in $(ls "$seqkit_out_dir/" | grep "part_"); 
  do
    log_info "Setting New Options For Contig: $contig"
    contig_id=$(seqkit fx2tab -n -i "$seqkit_out_dir/$contig")

    new_options=()
    options=$(getopt -o :hm:s:g:p:t:T:d:b:a:r:i:n:c:G -- "${args[@]}")
    eval set -- "$options"

    while true; do
        case "$1" in
            -p)
                # $2 contains the argument for -p
                new_options+=("-p" "${2}_${contig_id}")
                shift 2
                ;;
            -g)
                new_options+=("-g" "$contig")
                shift 2 
                ;;
            --)
                shift
                break
                ;;
            *)
                # For all other options, preserve them
                new_options+=("$1")
                shift
                ;;
        esac
    done
    log_info "Options Edited for Each Contig: Check cmd file"
    echo "${new_options[@]}" > $wd"/"$contig_id".cmd"

  done
}



# wrapper "$@"



run_bacmethy() {
  source $dir"/log_message_function.sh"


  local contig="$1"
  shift
  local args=("$@")
  
  contig_id=$(seqkit fx2tab -n -i "$contig")
  
  bacmethy_out_dir=$wd"/bacmethy_"$contig_id
  log_info "Create Dir for Output Files created at: $bacmethy_out_dir"
  mkdir $bacmethy_out_dir
  cd $bacmethy_out_dir
  wd=$bacmethy_out_dir

  log_info "Setting New Options For Contig: $contig"

  new_options=()
  options=$(getopt -o :hm:s:g:p:t:T:d:b:a:r:i:n:c:G -- "${args[@]}")
  eval set -- "$options"

  while true; do
      case "$1" in
          -p)
              # $2 contains the argument for -p
              new_options+=("-p" "${2}_${contig_id}")
              shift 2
              ;;
          -g)
              new_options+=("-g" "$contig")
              shift 2 
              ;;
          --)
              shift
              break
              ;;
          *)
              # For all other options, preserve them
              new_options+=("$1")
              shift
              ;;
      esac
  done
  log_warning "Options Edited for Specific Contig: Check job.txt"
  echo "${new_options[@]}" > $wd"/"$contig_id".log"

  log_info "Running Bacmethy on Contig: $(basename $contig)"
  bacmethy.sh ${new_options[@]} > $wd"/bacmethy_log.txt"
  log_info "See you later alligator!"
  cd ..

}


export -f run_bacmethy
export dir
export wd

# Get contigs and run parallel
# contigs=($(ls "$seqkit_out_dir/" | grep "part_"))

contigs=($(find $wd"/fasta_split/" -name "*part_*"))
#echo ${contigs[@]}
parallel --joblog job.txt --progress run_bacmethy {} "$@" ::: "${contigs[@]}"

rm -r $dir/tmp

