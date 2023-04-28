set term x11
set yrange [-3:3]
do for [i=1:1000] {p '"$1"'.xyz' index i w lp; pr i; pause -1}
