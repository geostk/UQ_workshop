# make sure to put blank lines between each "Scanline" to make pm3d work
set output "conc.jpg"
set terminal jpeg
# set pm3d at s
set pm3d map
unset key
# set grid
set key below
set border 4095
# set surface
# set samples 25
# set isosamples 20
set ticslevel 0
set title "Conc. 64 to 32 micron pyroclasts, kg/cu m"
set zrange [1E-05:1]
set logscale cb
set xlabel "Times, 14"
set ylabel "Times, 14"
set xlabel 'Distance, m'
set ylabel 'Height, m'
splot 'FIELD.TXT' u 1:2:3
