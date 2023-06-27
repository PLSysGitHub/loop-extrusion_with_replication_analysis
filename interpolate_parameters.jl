"""
This script was used to find the parameters to set the ori separation.

See Supplementary Information, and Fig. S1
"""

using CurveFit, DelimitedFiles, DSP, StatsBase, Statistics, Plots
include("src/helpers.jl")

pyplot(grid=false,label="",framestyle=:box,colorbar=true, linewidth=2,
    guidefontsize=15, tickfontsize=15,colorbar_tickfontsize=12,legend=:outerright, markersize=7,
    colorbar_titlefontsize=15, legendfontsize=15)

"""
Assume a constant deceleration until t_2, and then a constant velocity. Return the expected distance between origins of replication
"""
function slowing_oris(t, v_f; t_2=10, x_2=1.739)
    v_0=(2*x_2/t_2-v_f)
    a=(v_f-v_0)/t_2
    if t<t_2
        return v_0*t+1/2*a*t^2
    else
        return v_0*t_2+a/2*t_2^2+v_f*(t-t_2)
    end
end

"""
Given the separation between origins at given times, return parameters for constant deceleration upto t_2, and constant velocity after t_2

"""
function slowing_ori_params(times, ori_separations; t_2=10)
    v_f=linear_fit(times[2:end], ori_separations[2:end])[2]
    x_2=ori_separations[2]
    v_0=(2*x_2/t_2-v_f)
    a=(v_f-v_0)/t_2
    return v_f, v_0, a
end


#Path where to save figures
fig_folder="Plots/Interpolation/"
if !isdir(fig_folder)
    mkpath(fig_folder)
end

#data on ori separations extracted from microscopy experiments
times=[0,10,30,45,60, 75] #min
ori_separations=[0.0,1.739,2.088,2.399,2.671,3.000] #μm

#do fits to data
L_0=2.3 #μm, from Messelink et.al. 2021
r_growth=log(2)/126 #min^-1; doubling time 126 min
v_f, v_0, a=slowing_ori_params(times, ori_separations)

#make plots
plot(0:0.1:120, t->L_0*exp.(r_growth*t), label="Exponential growth", color=:gray, xlims=(-3,120),ylabel="Cell length [μm]", xlabel="Time [min]")
png(fig_folder*"growth_rate")

plot(0:0.1:80, t-> slowing_oris(t,v_f), label="Decelerating fit", color=:gray, xlims=(-3,80), ylabel="Origin separation [μm]",xlabel="Time [min]")
scatter!(times, ori_separations, color=:black)
png(fig_folder*"separation_rate")
