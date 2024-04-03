#dx login

#Imputation file path
imp_file_dir="/Bulk/Imputation/UKB imputation from genotype"

#Output folder
output="/Priyatam/PRScs/SCZ"


for i in {1..22}
do
	dx run swiss-army-knife -icmd="bgenix -g ukb22828_c${i}_b0_v3.bgen -incl-rsids rsid.txt > chr${i}.bgen" -iin="$imp_file_dir/ukb22828_c${i}_b0_v3.bgen" -iin="$imp_file_dir/ukb22828_c${i}_b0_v3.bgen.bgi" -iin="$output/rsid.txt" --destination "$output/" --instance-type "mem1_ssd1_v2_x16" --priority "normal" 
done
