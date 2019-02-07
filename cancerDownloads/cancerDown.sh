for i in {34,50,76,500,55,57,402,20,21,31,58,16,17,75,56,8,38,9,83,110,72,46,90,91,92,93,95,96,97,100,418,5,35,47,53,111,409,89,86,3,12,61,40,66,7,19,51,18,67,80,6,71,62,63}
     do wget -O $i\_data.csv  "https://seer.cancer.gov/explorer/download_data.php?site=$i&data_type=1&graph_type=3&compareBy=race&chk_sex_1=1&chk_race_1=1&hdn_data_type=1&advopt_precision=1&showDataFor=sex_1" 
done

