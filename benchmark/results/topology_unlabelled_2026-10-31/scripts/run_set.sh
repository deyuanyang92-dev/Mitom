#!/bin/bash
# run_set.sh MITOM_DIR SET_DIR LABEL [extra mitom options]: mitom run on every set_dir/in_<phylum>.fasta (no topology
# labels) with the phylum's genetic code and clade -> SET_DIR/LABEL/<phylum>/
# env: THREADS per record (default 2), JOBS records at a time (default 24), PHYLA="A B" only these
M=$1 D=$2 L=$3; shift 3
export PATH=$HOME/miniforge3/envs/mitotools/bin:$PATH OMP_NUM_THREADS=1
declare -A CODE=([Annelida]=5 [Cnidaria]=4 [Chordata]=2 [Porifera]=4 [Mollusca]=5 [Nematoda]=5 [Placozoa]=4)
declare -A CLADE=([Annelida]=Annelida-segmented-worms [Cnidaria]=Cnidaria [Chordata]=Chordata [Porifera]=Porifera-sponges \
                  [Mollusca]=Mollusca [Nematoda]=Nematoda [Placozoa]=Animal)
for f in $D/in_*.fasta; do
  p=$(basename $f .fasta); p=${p#in_}
  [ -n "${PHYLA:-}" ] && [[ " $PHYLA " != *" $p "* ]] && continue
  echo "== $p $(date)"
  python $M/mitom run -i $f -o $D/$L/$p --clade ${CLADE[$p]} --genetic_code ${CODE[$p]} --threads ${THREADS:-2} --jobs ${JOBS:-24} \
      --mitos_refdir $HOME/.wgs_skim/db/mitom/mitos_refs --plot no --keep_intermediate yes "$@" > $D/$L/$p.log 2>&1 \
      || echo "FAILED $p"
done
echo "== done $(date)"
