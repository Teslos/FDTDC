#set terminal ... enhanced
set xlabel "x (nm)"
set ylabel "y (nm)"
set zlabel "{/Symbol F}^2"
set title  "First eigenstate of the infinite well"
set hidden3d
set isosamples 90
splot "eigen_plot" u 1:2:3 w l
set terminal postscript eps enhanced
set output 'eigen_19.5meV.eps'
replot 
