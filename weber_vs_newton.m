set multiplot title "Simulation of 8192 objects initially in sphere of r = 1 at t = 0\nobjects' v₀ = 0, m = G¯²   simulation time: t∊[0,1], Δt = 0.1   units: c = 1, ξ = 1" layout 2, 2
set logscale
set ytics format "%T"
set xtics format "%T"
set title "Distance vs. Velocity"
set xlabel "log(Distance r)"
set ylabel "log(Velocity v)"
set key bottom right
plot 'Newton.txt' skip 3 u 5:9 title 'Newton', 'Weber.txt' skip 3 u 5:9 title 'Weber'
set title "Distance vs. Acceleration"
set xlabel "log(Distance r)"
set ylabel "log(Acceleration a)"
set key bottom left
plot 'Newton.txt' skip 3 u 5:13 title 'Newton', 'Weber.txt' skip 3 u 5:13 title 'Weber'
set title "Distance vs. Angular Momentum"
set xlabel "log(Distance r)"
set ylabel "log(Angular Momentum ℓ)"
set key top left
plot 'Newton.txt' skip 3 u 5:14 title 'Newton', 'Weber.txt' skip 3 u 5:14 title 'Weber'
set title "Distance vs. Torque"
set xlabel "log(Distance r)"
set ylabel "log(Torque τ)"
set key bottom left
plot 'Newton.txt' skip 3 u 5:15 title 'Newton', 'Weber.txt' skip 3 u 5:15 title 'Weber'
unset multiplot

