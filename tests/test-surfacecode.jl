
surfacecode_simulate(6,6,0.1)

(errorrate,abortrate) = surfacecode_simulate_many(8,8,0.1,maxtrials=100000,maxfailure=1000)

prange = 0.01:0.01:0.15
(errorrate,abortrate) = surfacecode_errorcurve(6,6,prange;maxtrials = 10000,maxfailure = 100)
plot(prange,errorrate, color = :blue)
plot!(prange,abortrate, color = :red)
plot!(yaxis=:log10, ylims=(1e-4, 1), xlims=(0, prange[end]))

prange = 0.01:0.01:0.15
(errorrates,abortrates) = surfacecode_thresholdcurve(4:2:10,prange;maxtrials=100000,maxfailure=1000)

plot(prange,errorrates', line = :solid, label=["L=4" "L=6" "L=8" "L=10"])
plot!(prange,abortrates', line = :dash, label=["L=4" "L=6" "L=8" "L=10"])
plot!(yaxis=:log10, ylims=(1e-7, 1), xlims=(0, prange[end]), label=["L=4","L=6","L=8","L=10"])
