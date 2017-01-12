set multiplot title "Simulation of 8192 objects: r = 1, v₀ = 0, t∊[0,1], Δt = 0.1" layout 2, 2
set logscale
set ytics format "%T"
set xtics format "%T"
set title "Distance vs. Velocity"
set xlabel "log(Distance r)"
set ylabel "log(Velocity v)"
set key bottom right
plot 'Weber.txt' skip 3 u 5:9 title 'Weber', 'Newton.txt' skip 3 u 5:9 title 'Newton'
set title "Distance vs. Acceleration"
set xlabel "log(Distance r)"
set ylabel "log(Acceleration a)"
set key bottom left
plot 'Weber.txt' skip 3 u 5:13 title 'Weber', 'Newton.txt' skip 3 u 5:13 title 'Newton'
set title "Distance vs. Angular Momentum"
set xlabel "log(Distance r)"
set ylabel "log(Angular Momentum ℓ)"
set key top left
plot 'Weber.txt' skip 3 u 5:14 title 'Weber', 'Newton.txt' skip 3 u 5:14 title 'Newton'
set title "Distance vs. Torque"
set xlabel "log(Distance r)"
set ylabel "log(Torque τ)"
set key bottom left
plot 'Weber.txt' skip 3 u 5:15 title 'Weber', 'Newton.txt' skip 3 u 5:15 title 'Newton'
unset multiplot

