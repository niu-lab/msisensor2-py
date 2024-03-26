set -e

model=./model_dir/
cat $model/sites_list.txt |grep -v "chrX"|grep -v "chrY"|grep -v "chrUn"|grep -v "random" > ./sites_list.txt

test_type=$1
for filename in /your_dis/$test_type/*; do
    echo $filename
    date
    echo -e "${filename##*/} \c" >> $test_type.txt
    python /your_dir/msisensor2.py $filename ./sites_list.txt $model/models/ >> $test_type.txt
done
