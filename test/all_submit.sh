set -e
for t in TCGA-CRC  ; do
    bsub -q queue_name -n 24 -W 1440 -o %J.out -e %J.err "source /home/your_dir/.bashrc; sh submit_test.sh $t"
done
