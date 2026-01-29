surfacecode_simulate(6,6,0.1)
surfacecode_simulate_many(8,8,0.1,maxtrials=100000,maxfailure=1000)
prange = 0.01:0.01:0.15
prob = surfacecode_errorcurve(6,6,0.01:0.01:0.3;maxtrials = 10000,maxfailure = 100)
plot(prange,prob)
plot!(yaxis=:log10, ylims=(1e-4, 1), xlims=(0, prange[end]))

proball = surfacecode_thresholdcurve(4:2:10,prange;maxtrials=100000,maxfailure=1000)
plot(prange,proball')
plot!(yaxis=:log10, ylims=(1e-7, 1), xlims=(0, prange[end]), labels=["L=4" "L=6" "L=8" "L=10"])





