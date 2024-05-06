set terminal postscript eps color "Helvetica" 42
set output "xtplot_case3.eps"
set size 2,1
set hidden3d
set pm3d
set view map
set xrange [0:40]
set yrange [0:10]
set cbrange [0:800]
set cbtics 0,800 offset -1,0
set cblabel "{/Helvetica [m*]}+{/Helvetica [m_{glassy}]}+{/Helvetica [m_{relax}]} ({/Symbol m}M)" offset-8.3,4.5 rotate by 0
unset ztics
unset xmtics
unset ymtics
unset surface
set xlabel "{/Helvetica-Italic t} / min" offset 0,0.2
set ylabel "{/Helvetica-Italic x} / {/Symbol m}m" offset 3.2,0
set xtics offset 0,0
set ytics offset 2.4,0
set palette model RGB defined (0 'black',1 'green')
splot [] [] "res" u ($1/60):($2*1e6):(($3+$5+$6)*1000) notitle
pause -1
