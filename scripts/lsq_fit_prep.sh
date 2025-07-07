F=$1

G=${F%.txt}    # removes the .txt string from the end of F

echo $G

for cvar in rsq std
do 
   echo "$cvar"_subset "$G".txt > "$G"_"$cvar".tmp
   grep "$cvar"_subset "$G".txt > "$G"_"$cvar".tmp
   echo ../../scripts/header_mth.txt "$G"_"$cvar".tmp > "$G"_"$cvar".txt  # was > /Summaries/"$G"_..
   cat ../../scripts/header_mth.txt "$G"_"$cvar".tmp > "$G"_"$cvar".txt  # was > /Summaries/"$G"_..
#   rm "$G"_"$cvar".tmp 
done 
