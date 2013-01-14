set multiplot title "Time domain data" layout 2,1
set xlabel "Time"
set ylabel "Psi real"
plot "timedomain" index 0 w l
set xlabel "E (meV)"
set ylabel "FFT (Psi)"
set xrange[0:50]
plot "timedomain" index 1 w l 
