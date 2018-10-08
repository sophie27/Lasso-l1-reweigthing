set terminal pdf dashed enhanced color  
set key font ",19" spacing 0.9 
set grid
set xlabel 'Sparsity' font ",22"
unset label
unset border

set grid xtics ytics mytics




set output 'rwl_1.pdf'
set ylabel 'Probability of recovery' font ",22"
#set key bottom left
set key samplen 2
plot 'RWL.dat' index 0 using 1:3 title "Alg.1-a"  with lp pt 4 ps 0.7  lw 1.5 lc rgb "#FF2000",\
'' index 1 using 1:3 title "Alg.2-a" with lp pt 5 ps 0.7 dashtype 4 lw 1.5 lc rgb "#120A8F",\
'' index 8 using 1:3 title "Alg.1-b" with lp pt 8 ps 0.7  lw 1.5 lc rgb "#FF2000",\
'' index 9 using 1:3 title "Alg.2-b" with lp pt 9 ps 0.7 dashtype 4 lw 1.5 lc rgb "#120A8F",\
'' index 4 using 1:3 title "Alg.1-c" with lp pt 6 ps 0.7  lw 1.5 lc rgb "#FF2000",\
'' index 5 using 1:3 title "Alg.2-c" with lp pt 7 ps 0.7 dashtype 4 lw 1.5 lc rgb "#120A8F",\
'' index 12 using 1:3 title "Lasso" with l lw 1.7 lc rgb "#000000",\


set output 'rwl_3.pdf'
set yrange [1e4:4e5]
set logscale y
set format y "10^{%L}"

set key bottom right
set ylabel 'Number of iterations' font ",22"
plot 'RWL.dat' index 0 using 1:4 title "Alg.1-a"  with lp pt 4 ps 0.7  lw 1.5 lc rgb "#FF2000",\
'' index 1 using 1:4 title "Alg.2-a" with lp pt 5 ps 0.7 dashtype 4 lw 1.5 lc rgb "#120A8F",\
'' index 8 using 1:4 title "Alg.1-b" with lp pt 8 ps 0.7  lw 1.5 lc rgb "#FF2000",\
'' index 9 using 1:4 title "Alg.2-b" with lp pt 9 ps 0.7 dashtype 4 lw 1.5 lc rgb "#120A8F",\
'' index 4 using 1:4 title "Alg.1-c" with lp pt 6 ps 0.7  lw 1.5 lc rgb "#FF2000",\
'' index 5 using 1:4 title "Alg.2-c" with lp pt 7 ps 0.7 dashtype 4 lw 1.5 lc rgb "#120A8F",\
'' index 12 using 1:4 title "Lasso" with l  lw 1.7 lc rgb "#000000",\


set output 'rwl_5.pdf'
set yrange [2e3:2e5]
set ylabel 'Number of iterations' font ",22"
plot 'RWL_NEW.dat' index 2 using 1:4 title "Alg.1-a"  with lp pt 4 ps 0.7  lw 1.5 lc rgb "#FF2000",\
'' index 3 using 1:4 title "Alg.2-a" with lp pt 5 ps 0.7 dashtype 4 lw 1.5 lc rgb "#120A8F",\
'' index 10 using 1:4 title "Alg.1-b" with lp pt 8 ps 0.7  lw 1.5 lc rgb "#FF2000",\
'' index 11 using 1:4 title "Alg.2-b" with lp pt 9 ps 0.7 dashtype 4 lw 1.5 lc rgb "#120A8F",\
'' index 6 using 1:4 title "Alg.1-c" with lp pt 6 ps 0.7  lw 1.5 lc rgb "#FF2000",\
'' index 7 using 1:4 title "Alg.2-c" with lp pt 7 ps 0.7 dashtype 4 lw 1.5 lc rgb "#120A8F",\
'' index 13 using 1:4 title "Lasso" with l  lw 1.7 lc rgb "#000000",\

unset yrange

set output 'rwl_2.pdf'
set logscale y
set format y "10^{%L}"
set key bottom right
set ylabel 'Relative square error' font ",22"
plot 'RWL.dat' index 0 using 1:2 title "Alg.1-a"  with lp pt 4 ps 0.7  lw 1.5 lc rgb "#FF2000",\
'' index 1 using 1:2 title "Alg.2-a" with lp pt 5 ps 0.7 dashtype 4 lw 1.5 lc rgb "#120A8F",\
'' index 8 using 1:2 title "Alg.1-b" with lp pt 8 ps 0.7  lw 1.5 lc rgb "#FF2000",\
'' index 9 using 1:2 title "Alg.2-b" with lp pt 9 ps 0.7 dashtype 4 lw 1.5 lc rgb "#120A8F",\
'' index 4 using 1:2 title "Alg.1-c" with lp pt 6 ps 0.7  lw 1.5 lc rgb "#FF2000",\
'' index 5 using 1:2 title "Alg.2-c" with lp pt 7 ps 0.7 dashtype 4 lw 1.5 lc rgb "#120A8F",\
'' index 12 using 1:2 title "Lasso" with l  lw 1.7 lc rgb "#000000",\




set output 'rwl_4.pdf'
set logscale y
set format y "10^{%L}"
set yrange [1e-3:5e-1]
set key bottom right
set ylabel 'Relative square error' font ",22"
plot 'RWL.dat' index 2 using 1:2 title "Alg.1-a"  with lp pt 4 ps 0.7  lw 1.5 lc rgb "#FF2000",\
'' index 3 using 1:2 title "Alg.2-a" with lp pt 5 ps 0.7 dashtype 4 lw 1.5 lc rgb "#120A8F",\
'' index 10 using 1:2 title "Alg.1-b" with lp pt 8 ps 0.7  lw 1.5 lc rgb "#FF2000",\
'' index 11 using 1:2 title "Alg.2-b" with lp pt 9 ps 0.7 dashtype 4 lw 1.5 lc rgb "#120A8F",\
'' index 6 using 1:2 title "Alg.1-c" with lp pt 6 ps 0.7  lw 1.5 lc rgb "#FF2000",\
'' index 7 using 1:2 title "Alg.2-c" with lp pt 7 ps 0.7 dashtype 4 lw 1.5 lc rgb "#120A8F",\
'' index 13 using 1:2 title "Lasso" with l  lw 1.7 lc rgb "#000000",\








