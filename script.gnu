set logscale x; set logscale y
set term png
set output "graph.png"
set xlabel  "k"
set ylabel  "|| B-B(k) ||"
plot  "dataCross.dat" using 1:2 with linespoints, "dataEigen.dat" using 1:2 with linespoints

