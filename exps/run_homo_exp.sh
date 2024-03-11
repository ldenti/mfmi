#!/usr/bin/env sh

set -e

function log {
    >&2 echo -e "\n##################################"
    >&2 echo "### $(date) # $1"
    >&2 echo -e "##################################\n"
}

mfmi=0
rope=0
rlcsa=0

threads=1
WD="."
while getopts "@:w:123" flag; do
    case "${flag}" in
	h) $(>&2 echo "HELP")
	   exit 0
	   ;;
	@) threads=${OPTARG}
	   ;;
	w) WD=${OPTARG}
	   ;;
	1) mfmi=1
	   ;;
	2) rope=1
	   ;;
	3) rlcsa=1
	   ;;
    esac
done

if [[ $# -lt $((${OPTIND} + 1)) ]]
then
    (>&2 echo "ERROR")
    exit 1
fi

REF=${@:$OPTIND:1}
Q=${@:$OPTIND+1:1}

mkdir -p $WD/times

if [ "$mfmi" -eq "1" ]
then
    INDEX=$WD/mfmi-index
    log "MFMI indexing"
    \time -vo $WD/times/mfmi-index.time ../mfmi index -r -@$threads $REF > $INDEX
    # log "MFMI pingpong"
    # \time -vo $WD/times/mfmi-pingpong.time ../mfmi pingpong $INDEX $Q > $WD/mfmi.sfs
fi

if [ "$rope" -eq "1" ]
then
    INDEX=$WD/rope-index
    log "ROPEBWT2 indexing (buffered with no threads)"
    \time -vo $WD/times/rope-index.time ./ext/pp-ropebwt2/rt2 index $REF > $INDEX
    log "ROPEBWT2 indexing (buffered with threads)"
    \time -vo $WD/times/rope-index.wthreads.time ./ext/pp-ropebwt2/rt2 index -@ $REF > $INDEX.wthreads
    # log "ROPEBWT2 pingpong"
    # \time -vo $WD/times/rope-pingpong.time ./ext/pp-ropebwt2/rt2 pingpong $INDEX $Q > $WD/rope.sfs
fi

if [ "$rlcsa" -eq "1" ]
then
    INDEX=$WD/rlcsa-index
    log "RLCSA indexing (buffered with threads)"
    \time -vo $WD/times/rlcsa-index.time ./ext/pp-rlcsa/rl index -@$threads -i $INDEX $REF
    # log "RLCSA pingpong"
    # \time -vo $WD/times/rlcsa-pingpong.time ./ext/pp-rlcsa/rl pingpong $INDEX $Q > $WD/rlcsa.sfs
fi

log "DONE"
