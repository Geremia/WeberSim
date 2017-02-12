set title "Weberian Force\nSimulation of 1024 objects initially on spherical shell of r = 1 at t = 0\nobjects' v₀ = 0, m = 1   simulation time: t∊[0,1], Δt = 0.005   units: c = 1, ξ = 1"
unset ytics
set xlabel "Time"
set ylabel "Energy (linear scale in arbitrary units)"
set key top right
plot 'stats.txt' skip 3 u 1:2 title "Kinetic Energy", 'stats.txt' skip 3 u 1:3 title "Weberian Potential Energy", 'stats.txt' skip 3 u 1:($2+$3) title "Total Energy"
