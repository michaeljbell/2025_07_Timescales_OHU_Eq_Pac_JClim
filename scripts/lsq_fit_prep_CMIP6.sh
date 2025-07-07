F=$1

G=${F%.txt}    # removes the .txt string from the end of F

for cvar in 0 1
do 
  grep x_soln_subset.ifit."$cvar" "$G".txt > "$G"_reg_coeff_"$cvar".tmp
  cat ../../scripts/header_mth.txt "$G"_reg_coeff_"$cvar".tmp > Summaries/"$G"_reg_coeff_"$cvar".txt
  rm "$G"_reg_coeff_"$cvar".tmp 
done 
