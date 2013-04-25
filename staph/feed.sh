while read p; do
  wget $p
done < links.txt
