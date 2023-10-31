### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ f4c44654-5802-11ee-194a-0fc1aed8aabc
begin
	import Gaston,PGFPlotsX, Contour
	using LinearAlgebra,  Statistics
	using Plots, LaTeXStrings
	using HDF5, FileIO
	using PlutoUI
end

# ╔═╡ 50b414b9-bda7-410c-8624-9ce7a0fd08a1
md"
# Plot figures for 2Q paper
"

# ╔═╡ 0eff32b2-2594-41f4-98d0-9d9d650204ff
# Set up path variables
begin
	PATH_DATA = "../Data/"
	PATH_FIGS = "../Figs/"
end

# ╔═╡ c944c2bc-40fe-4405-a8fc-6d81434662d0
# Load data function
"""
    load_dict(data_type::String; N_exc::Int = 2, direction::String = "L", a_inc::Real = 0.03, delta::Real = -0.15)

Load a dictionary from an HDF5 file based on the specified data type and parameters.

# Arguments

- `data_type::String`: The type of data to load. Supported values are "time," "power," "power_Nexc," "power_t," "meanfield," and "power_Lj."

# Keyword Arguments

- `N_exc::Int = 2`: The number of excitations (default: 2).
- `direction::String = "L"`: The direction (default: "L").
- `a_inc::Real = 0.03`: The increment value (default: 0.03).
- `delta::Real = -0.15`: The delta value (default: -0.15).
- `order::Int = 2`: The delta value (default: 2).
- `version::String = ""`: The version value (default: "").

# Returns

- A dictionary loaded from the corresponding HDF5 file.

The function constructs the file name based on the `data_type` and the provided parameters, and then loads the dictionary from the file located at `PATH_FIGS * FILE_NAME`.

# Examples

```julia
dict = load_dict("time", N_exc=3, direction="R", a_inc=0.03, delta=-0.15,
				         order=2, version="")
```
"""
function load_dict(data_type::String; N_exc::Int = 2, 
				   direction::String = "L", a_inc::Real = 0.03, delta::Real = -0.15,
				   order::Int = 2,
				   version::String = "")
	if data_type == "time"
		FILE_NAME = data_type*"_dependence_Nexc"*string(N_exc)*"_"*direction*
					"_ainc_"*string(a_inc)*"_delta_"*string(delta)*version*".h5"
	elseif data_type == "power"
		FILE_NAME = data_type*"_dependence_Nexc"*string(N_exc)*
					"_delta_"*string(delta)*version*".h5"
	elseif data_type == "power_t"
		FILE_NAME = "power_dependence_t_Nexc"*string(N_exc)*
					"_delta_"*string(delta)*version*".h5"
	elseif data_type == "power_Lj"
		FILE_NAME = data_type*"_dependence_Nexc"*string(N_exc)*version*".h5"
	elseif data_type == "meanfield"
		FILE_NAME = data_type*"_sigma_order"*string(order)*"_ainc_"*string(a_inc)*
					"_delta_"*string(delta)*".h5"
	elseif data_type == "power_Nexc"
		FILE_NAME = "power_dependence_Nexc"*"_delta_"*string(delta)*"_v1.1.h5"
	elseif data_type == "states"
		FILE_NAME = "states_evo_Nexc"*string(3)*"_"*direction*
					"_ainc_"*string(a_inc)*"_delta_"*string(-0.015)*version*".h5"
	end
	return load(PATH_DATA*FILE_NAME)
end

# ╔═╡ 7f3041e6-25d2-40d4-98ad-236ecaf4ef40
# Sets of avaliable parameters
begin
	data_types = ["time", "power", "power_t", "power_Nexc", "power_Lj", "meanfield",
				  "states"]
	directions = ["L", "R"]
	a_incs = [0.03]
	deltas = Dict("time" => [-0.15, -0.015], "power" => [-0.15, -0.015, -0.0015],
				  "power_t" => [-0.15], "meanfield" => [-0.15, -0.015, -0.0015])
	orders = [1, 2]
	N_excs = [2, 3, 5]
end

# ╔═╡ 63847c84-bd9e-4bf7-be5d-5e90505b9386
# Load data
begin
	time_2_15_data = load_dict("time")
	power_2_15_data = load_dict("power")
	data_time = [load_dict("time"; N_exc=i, direction=j, a_inc=k,
							delta=m, version="_v1.1") 
		for i in N_excs, j in directions, k in a_incs, m in deltas["time"]]
	data_power = [load_dict("power"; N_exc=i, delta=m, version="_v1.1") 
		for i in N_excs, m in deltas["power"]]
	data_power_t = [load_dict("power_t"; N_exc=i, delta=m, version="_v1.1") 
		for i in N_excs, m in deltas["power_t"]]
	data_power_Lj = [load_dict("power_Lj"; N_exc=i, version="") 
		for i in N_excs]
	data_power_Nexc = [load_dict("power_Nexc", delta=i) for i in deltas["power"]] 
	data_meanfield = [load_dict("meanfield"; order=i, delta=m, version="") 
		for i in orders, m in deltas["meanfield"]]
	data_meanfield_small = [load_dict("meanfield"; order=i, delta=-0.15, a_inc=m, version="") 
		for i in orders, m in [0.0003, 3e-6]]
	data_states = load_dict("states")
	print("data loaded")
end

# ╔═╡ 3ba29c8b-8c1c-4967-b3d7-9568bccd232c
md"
## Power and $L_j$ dependence
"

# ╔═╡ c8321163-99e6-4e14-8f93-45f308370940
data_power_t[3]

# ╔═╡ 1c744290-a859-4a2d-8e4e-cd951331a476
# Fig 3: transmission, efficiency, populations vs power
let
	# pgfplotsx()
	gr()
	ind = 1
	x = data_power[1, 1]["a_inc_list"]
	x_t = data_power_t[1]["a_inc_list"]
	y = [data_power[i, ind]["transmission_L"] for i in eachindex(N_excs)]
	yt = [data_power_t[i]["transmission_L"] for i in eachindex(N_excs)] 
	y_1 = [data_power[i, ind]["transmission_R"] for i in eachindex(N_excs)]
	y_1t = [data_power_t[i]["transmission_R"] for i in eachindex(N_excs)] 
	y_2 = [data_power[i, ind]["efficiency"] for i in eachindex(N_excs)]
	y_2t = [data_power_t[i]["efficiency"] for i in eachindex(N_excs)]
	y_31 = [data_power[i, ind]["population_D_L"] for i in eachindex(N_excs)]
	y_32 = [data_power[i, ind]["population_D_R"] for i in eachindex(N_excs)]
	y_31t = [data_power_t[i]["population_D_L"] for i in eachindex(N_excs)]
	y_32t = [data_power_t[i]["population_D_R"] for i in eachindex(N_excs)]
	y_41 = [data_power[i, ind]["population_B_L"] for i in eachindex(N_excs)]
	y_42 = [data_power[i, ind]["population_B_R"] for i in eachindex(N_excs)]
	y_41t = [data_power_t[i]["population_B_L"] for i in eachindex(N_excs)]
	y_42t = [data_power_t[i]["population_B_R"] for i in eachindex(N_excs)]
	γ = mean(power_2_15_data["gammas"])
	cs_1 = cgrad(:diverging_bwr_20_95_c54_n256)
	cs_2 = cgrad(:diverging_bwr_20_95_c54_n256, rev=true)
	annotation_pos = [[[[7e-2 for i in N_excs],[0.66, 0.51, 0.38]],
					   [[6e-3, 4e-2, 4e-2], [0.48, 0.44, 0.32]],
					   [[7e1 for i in N_excs], [0.28, 0.14, 0.06]],
					   [[7e1 for i in N_excs], [0.27, 0.16, 0.02]]],
					  [[[7e-2 for i in N_excs],[0.7, 0.55, 0.38]],
					   [[7e-5, 2e-2, 2e-2], [0.55, 0.54, 0.45]],
					   [[6e1 for i in N_excs], [0.28, 0.14, 0.06]],
					   [[7e1 for i in N_excs], [0.27, 0.16, 0.02]]],
					  [[[7e-2 for i in N_excs],[0.7, 0.55, 0.38]],
					   [[1e-6, 2e-2, 2e-2], [0.6, 0.55, 0.44]],
					   [[5e1 for i in N_excs], [0.28, 0.14, 0.06]],
					   [[7e1 for i in N_excs], [0.27, 0.16, 0.02]]]
					]
	x_lims = [(1e-5,(x.^2.0./γ)[end]),
			  (1e-7,(x.^2.0./γ)[end]),
			  (1e-8,(x.^2.0./γ)[end])]

	l = @layout grid(2,2)
	# Transmissions
	p1 = plot(x.^2.0./γ, y,
		framestyle=:box,
		xscale=:log10,
		minorgrid=true,
		xlims=x_lims[ind],
		lw=1.5,
		palette = cs_1[20:length(cs_1)÷8:length(cs_1)÷3],
		title="(a)",titleloc = :left, titlefont = font(14),
		xlabel=L"|a_\mathrm{inc}|^2 / \bar{\gamma}",
		ylabel="Transmission",
		label=reshape([(i == 1) ? L"\mathrm{L}" : :none for i in eachindex(N_excs)], (1, length(N_excs))),
		tickfontsize=12,
		labelfontsize=14,
		legend_font=12,
		#legend=:inside,
	)
	scatter!(x_t.^2.0./γ, yt, markershape=:circle, ms=3, label=:none)
	plot!(x.^2.0./γ, y_1,
		lw=1.5,
		ls=:dashdot,
		palette = cs_2[20:length(cs_1)÷8:length(cs_1)÷3],
		label=reshape([(i == 1) ? L"\mathrm{R}" : :none for i in eachindex(N_excs)], (1, length(N_excs))),
	)
	scatter!(x_t.^2.0./γ, y_1t, markershape=:circle, ms=3, label=:none)
	annotate!(annotation_pos[ind][1][1], annotation_pos[ind][1][2], [text("Nₗ = "*"$(N_excs[i])", :center, 13, cs_1[(20:length(cs_1)÷8:length(cs_1)÷3)[i]]) for i in eachindex(N_excs)])
	# Efficiencies
	p2 = plot(x.^2.0./γ, y_2,
		framestyle=:box,
		xscale=:log10,
		minorgrid=true,
		xlims=x_lims[ind],
		lw=1.5,
		palette = cs_1[20:length(cs_1)÷8:length(cs_1)÷3],
		title="(b)",titleloc = :left, titlefont = font(14),
		xlabel=L"|a_\mathrm{inc}|^2 / \bar{\gamma}",
		ylabel="Efficiency",
		tickfontsize=12,
		labelfontsize=14,
		legend=:none,
	)
	scatter!(x_t.^2.0./γ, y_2t, markershape=:circle, ms=3)
	annotate!(annotation_pos[ind][2][1], annotation_pos[ind][2][2], [text("Nₗ = "*"$(N_excs[i])", :center, 13, cs_1[(20:length(cs_1)÷8:length(cs_1)÷3)[i]], rotation = (i==1) ? 60 : 0) for i in eachindex(N_excs)])
	# Populations D
	p3 = plot(x.^2.0./γ, y_31,
		framestyle=:box,
		xscale=:log10,
		minorgrid=true,
		xlims=x_lims[ind],
		lw=1.5,
		palette = cs_1[20:length(cs_1)÷8:length(cs_1)÷3],
		title="(c)",titleloc = :left, titlefont = font(14),
		xlabel=L"|a_\mathrm{inc}|^2 / \bar{\gamma}",
		ylabel=L"\mathrm{Population} \; |D\rangle",
		tickfontsize=12,
		labelfontsize=14,
		legend=:none,
	)
	scatter!(x_t.^2.0./γ, y_31t, markershape=:circle, ms=3)
	plot!(x.^2.0./γ, y_32,
		lw=1.5,
		ls=:dashdot,
		palette = cs_2[20:length(cs_1)÷8:length(cs_1)÷3],
		label=reshape([(i == 1) ? L"\mathrm{R}" : :none for i in eachindex(N_excs)], (1, length(N_excs))),
	)
	scatter!(x_t.^2.0./γ, y_32t, markershape=:circle, ms=3)
	annotate!(annotation_pos[ind][3][1], annotation_pos[ind][3][2], [text("Nₗ = "*"$(N_excs[i])", :center, 13, cs_1[(20:length(cs_1)÷8:length(cs_1)÷3)[i]]) for i in eachindex(N_excs)])
	# Populations B
	p4 = plot(x.^2.0./γ, y_41,
		framestyle=:box,
		xscale=:log10,
		minorgrid=true,
		xlims=(1e-5,(x.^2.0./γ)[end]),
		lw=1.5,
		palette = cs_1[20:length(cs_1)÷8:length(cs_1)÷3],
		title="(d)",titleloc = :left, titlefont = font(14),
		xlabel=L"|a_\mathrm{inc}|^2 / \bar{\gamma}",
		ylabel=L"\mathrm{Population} \; |B\rangle",
		tickfontsize=12,
		labelfontsize=14,
		legend=:none,
	)
	scatter!(x_t.^2.0./γ, y_41t, markershape=:circle, ms=3)
	plot!(x.^2.0./γ, y_42,
		lw=1.5,
		ls=:dashdot,
		palette = cs_2[20:length(cs_1)÷8:length(cs_1)÷3],
		label=reshape([(i == 1) ? L"\mathrm{R}" : :none for i in eachindex(N_excs)], (1, length(N_excs))),
	)
	scatter!(x_t.^2.0./γ, y_42t, markershape=:circle, ms=3)
	annotate!(annotation_pos[ind][4][1], annotation_pos[ind][4][2], [text("Nₗ = "*"$(N_excs[i])", :center, 13, cs_1[(20:length(cs_1)÷8:length(cs_1)÷3)[i]]) for i in eachindex(N_excs)])
	p = plot(p1, p3, p2, p4, layout=l, size=(1000,700))
	# savefig(p, PATH_FIGS*"fig_3_delta"*string(deltas["power"][ind])*".pdf")
end

# ╔═╡ 1e7e0136-a3de-4f49-8015-bb16166aeb8f
data_power_Lj[1]

# ╔═╡ 047c1b14-28e9-4237-b50f-ce2d9c3fd491
# Fig 5
let
	# pgfplotsx()
	# gr()
	gaston()

	γ = mean(data_power[1]["gammas"])
	
	x = [data_power_Lj[i]["L_j_list"] for i in eachindex(N_excs)] ./ 1e-9
	y = [data_power_Lj[i]["a_inc_list"].^2.0./γ for i in eachindex(N_excs)]
	z_1 = [data_power_Lj[i]["efficiency"] for i in eachindex(N_excs)]
	z_2 = [data_power_Lj[i]["transmission_L"] for i in eachindex(N_excs)]
	z_3 = [data_power_Lj[i]["transmission_R"] for i in eachindex(N_excs)]

	kw = (
		yscale=:log10,
		xlims=(3.22, 3.38),
		ylims=(1e-6, 1e1),
		grid = true, 
		minorgrid = false,
		tickfontfamily=:sanserif,
		fontfamily=:sanserif,
		tickfontsize=12,
		labelfontsize=14,
		legend_font=12,
	)
	kw_1 = (
		xlabel = "Lⱼ (nH)", 
		ylabel = L"|a_\mathrm{inc}|^2 / \bar{\gamma}", 
		zlabel = "Efficiency", 
		legend = false,
		title = "(a)",titleloc = :left, titlefont = font(14),
		right_margin=(10, :mm),
		camera=(30,40),
	)
	kw_2 = (
		xlabel = "Lⱼ (nH)", 
		ylabel = L"|a_\mathrm{inc}|^2 / \bar{\gamma}", 
		cbar_title = "Transmission", 
		colorbar_titlefontsize=14, 
		colorbar_tickfontsize=12,
		legend = false,
		title = "(b) L direction, Nₗ=2",titleloc = :left, titlefont = font(14),
		bottom_margin=(5, :mm),
	)
	kw_3 = (
		xlabel = "Lⱼ (nH)", 
		# ylabel = L"|a_\mathrm{inc}|^2 / \bar{\gamma}", 
		cbar_title = "Transmission", 
		colorbar_titlefontsize=14, 
		colorbar_tickfontsize=12,
		legend = true,
		title = "(c) R direction, Nₗ=2",titleloc = :left, titlefont = font(14),
		bottom_margin=(5, :mm),
	)

	l = @layout [a{0.5w} grid(2,2)]
	
	p_1 = plot(0,0,0)
	## we add to the graphic p, then plot
	for cl in Contour.levels(Contour.contours(x[1], y[1], z_1[1]', 60))
	    lvl = Contour.level(cl) # the z-value of this contour level
	    for line in Contour.lines(cl)
	        _xs, _ys = Contour.coordinates(line) # coordinates of this line segment
	        _zs = 0 * _xs
			plot!(p_1, _xs, _ys, lvl .+ _zs, alpha=0.5, color=:red, legend=false; kw..., kw_1...) # add on surface
			plot!(p_1, _xs, _ys.*0 .+ 1e1, lvl .+ _zs, lw = 0.2, alpha=0.85, color=:gray, legend=false) # y projection
			plot!(p_1, _xs*0 .+ 3.22, _ys, lvl .+ _zs, lw = 0.2, alpha=0.85, color=:gray, legend=false) # x projection
	    end
	end
	for cl in Contour.levels(Contour.contours(x[3], y[3], z_1[3]', 60))
	    lvl = Contour.level(cl) # the z-value of this contour level
	    for line in Contour.lines(cl)
	        _xs, _ys = Contour.coordinates(line) # coordinates of this line segment
	        _zs = 0 * _xs
	        plot!(p_1, _xs, _ys, lvl .+ _zs, alpha=0.5, color=:blue, legend=false) # add on surface
			plot!(p_1, _xs, _ys.*0 .+ 1e1, lvl .+ _zs, lw = 0.2, alpha=0.85, color=:black, legend=false) # y projection
			plot!(p_1, _xs*0 .+ 3.22, _ys, lvl .+ _zs, lw = 0.2, alpha=0.85, color=:black, legend=false) # x projection
	    end
	end
	# savefig(p_1, PATH_FIGS*"fig_5_a.pdf")

	p_2 = contourf(x[1], y[1], z_2[1], color=:devon ; kw..., kw_2...)
	for cl in Contour.levels(Contour.contours(x[1], y[1], z_2[1]', 10))
	    lvl = Contour.level(cl) # the z-value of this contour level
	    for line in Contour.lines(cl)
	        _xs, _ys = Contour.coordinates(line) # coordinates of this line segment
	        _zs = 0 * _xs
			plot!(p_2, _xs, _ys, _zs, 
				lw = 1, alpha=0.85, color=:white, legend=false)
	    end
	end

	p_3 = contourf(x[1], y[1], z_3[1], color=:devon ; kw..., kw_3...)
	for cl in Contour.levels(Contour.contours(x[1], y[1], z_3[1]', 10))
	    lvl = Contour.level(cl) # the z-value of this contour level
	    for line in Contour.lines(cl)
	        _xs, _ys = Contour.coordinates(line) # coordinates of this line segment
	        _zs = 0 * _xs
			plot!(p_3, _xs, _ys, _zs, 
				lw = 1, alpha=0.85, color=:white, legend=false)
	    end
	end

	p_4 = contourf(x[3], y[3], z_2[3], color=:devon ; kw..., kw_2..., title="(d) L direction, Nₗ=5")
	for cl in Contour.levels(Contour.contours(x[3], y[3], z_2[3]', 10))
	    lvl = Contour.level(cl) # the z-value of this contour level
	    for line in Contour.lines(cl)
	        _xs, _ys = Contour.coordinates(line) # coordinates of this line segment
	        _zs = 0 * _xs
			plot!(p_4, _xs, _ys, _zs, 
				lw = 1, alpha=0.85, color=:white, legend=false)
	    end
	end

	p_5 = contourf(x[3], y[3], z_3[3], color=:devon ; kw..., kw_3..., title="(e) R direction, Nₗ=5")
	for cl in Contour.levels(Contour.contours(x[3], y[3], z_3[3]', 10))
	    lvl = Contour.level(cl) # the z-value of this contour level
	    for line in Contour.lines(cl)
	        _xs, _ys = Contour.coordinates(line) # coordinates of this line segment
	        _zs = 0 * _xs
			plot!(p_5, _xs, _ys, _zs, 
				lw = 1, alpha=0.85, color=:white, legend=false)
	    end
	end
	p = plot(p_1, p_2, p_3, p_4, p_5, layout=l, size=(1000,500))
	# savefig(p, PATH_FIGS*"fig_5.pdf")
	# Plots.svg(p, PATH_FIGS*"fig_5.svg")
	p
end

# ╔═╡ 68dbe732-236d-4af2-b98d-0a2f8d2afcfa
md"
## Correlation functions and spectra
"

# ╔═╡ ee86760b-9308-4a25-b8b8-60adb442b7a9
data_time[1,1,1,1]

# ╔═╡ e3eeecd4-0c17-45c8-b83e-cd071a02bad5
# Fig 4: g1, g2, spectra
let
	# pgfplotsx()
	gr()
	ind = [1, 1, 1] # N_exc; DIRECTION: 1-L, 2-R; delta
	g2_start_tau = 2
	
	x = [data_time[i, ind[2], 1, ind[3]]["taulist"] for i in eachindex(N_excs)]
	y_1 = [data_time[i, ind[2], 1, ind[3]]["g1_L"] for i in eachindex(N_excs)]
	y_11 = [data_time[i, ind[2], 1, ind[3]]["g1_R"] for i in eachindex(N_excs)]
	x_1 = [data_time[i, ind[2], 1, ind[3]]["taulist"][g2_start_tau:end] for i in eachindex(N_excs)]
	y_2 = [data_time[i, ind[2], 1, ind[3]]["g2_L"][g2_start_tau:end] for i in eachindex(N_excs)]
	y_21 = [data_time[i, ind[2], 1, ind[3]]["g2_R"][g2_start_tau:end] for i in eachindex(N_excs)]
	x_3 = [data_time[i, ind[2], 1, ind[3]]["w_spec_t"] for i in eachindex(N_excs)]
	y_3 = [data_time[i, ind[2], 1, ind[3]]["spec_L_t"] for i in eachindex(N_excs)]
	y_31 = [data_time[i, ind[2], 1, ind[3]]["spec_R_t"] for i in eachindex(N_excs)]
	y_3_max = maximum(maximum([data_time[i, ind[2], 1, ind[3]]["spec_R"] for i in eachindex(N_excs)]))
	
	if ind[3] == 1 # delta = -0.15
		γ_D = mean(data_time[ind[1], ind[2], 1, ind[3]]["g_D(L_j)"][255:256])
	elseif ind[3] == 2 # delta = -0.015
		γ_D = data_time[ind[1], ind[2], 1, ind[3]]["g_D(L_j)"][251]
	end

	annotation_pos = [[
					   [[4.65, 4.65, 4.65], [0.62, 0.53, 0.45]],
					   [[4, 4, 4], [0.12, 0.32, 0.4]],
					   [[0.4, 0.9, 0.9], [1.1, 2.2, 2.08]],
					   [[1.0, 1.87, 1.87], [0.75, 2.4, 2.53]],
					   [[0.0, 0.0, 0.0], [0.3, 0.42, 0.48]],
					   [[-3.4, -1.7, -0.5], [0.9, 0.75, 0.6]]
					  ],
					  [
					   [[1.7, 1.7, 1.7], [0.69, 0.96, 0.9]],
					   [[1.7, 1.7, 1.7], [0.37, 0.43, 0.48]],
					   [[0.9, 0.9, 0.9], [1.95, 1.65, 1.8]],
					   [[0.5, 0.9, 0.9], [1.35, 2.3, 2.15]],
					   [[1.4, 1.4, 1.4], [0.25, 0.38, 0.44]],
					   [[-0.1, -0.05, -0.05], [0.28, 0.37, 0.19]]  
					  ]
					 ]
	lens_pos = [
		[
			[[1.2, 3.5], [0.484, 0.497],0.58, 0.02],
			[[1, 3.5], [0.35, 0.37],0.16, 0.02],
			[[0.45, 0.55], [0.83, 0.84],0.12, 0.05],
			[[0.45, 0.55], [0.73, 0.74],0.54, 0.05]
		],
		[
			[[1.5, 1.7], [0.64, 0.645],0.38, 0.05],
			[[1.5, 1.7], [0.62, 0.63],0.38, 0.52],
			[[0.6, 0.7], [1.005, 1.02],0.25, 0.05],
			[[0.45, 0.55], [0.73, 0.74],0.54, 0.05]
		]
	]
					
	x_lims = [[(-0.1,5),
			  (-0.05,2),
			  (-20,20)],
			  [(-0.05,2),
			  (-0.025,1),
			  (-0.2,0.2)]]
	
	cs_1 = cgrad(:diverging_bwr_20_95_c54_n256)
	cs_2 = cgrad(:diverging_bwr_20_95_c54_n256, rev=true)
	l = @layout [grid(2,1) b{0.5w}] 
	# g1: L
	p1 = plot(x.*γ_D, real(y_1),
		framestyle=:box,
		#xscale=:linear,
		# minorgrid=true,
		xlims=x_lims[ind[3]][1],
		lw=1.5,
		palette = cs_1[20:length(cs_1)÷8:length(cs_1)÷3],
		title="(a)",titleloc = :left, titlefont = font(14),
		xlabel=L"\gamma_D \tau",
		ylabel=L"g^{(1)}(\tau)",
		label=reshape([(i == 1) ? L"\mathrm{Tran}" : :none for i in eachindex(N_excs)], (1, length(N_excs))),
		tickfontsize=12,
		labelfontsize=14,
		legend_font=12,
		legend=:none,
	)
	plot!(x.*γ_D, real(y_11),
		lw=1.5,
		ls=:dashdot,
		palette = cs_2[20:length(cs_1)÷8:length(cs_1)÷3],
		label=reshape([(i == 1) ? L"\mathrm{Ref}" : :none for i in eachindex(N_excs)], (1, length(N_excs))),
	)
	annotate!(annotation_pos[ind[3]][1][1], annotation_pos[ind[3]][1][2], [text("Nₗ = "*"$(N_excs[i])", :center, 13, cs_1[(20:length(cs_1)÷8:length(cs_1)÷3)[i]]) for i in eachindex(N_excs)])
	annotate!(annotation_pos[ind[3]][2][1], annotation_pos[ind[3]][2][2], [text("Nₗ = "*"$(N_excs[i])", :center, 13, cs_2[(20:length(cs_1)÷8:length(cs_1)÷3)[i]]) for i in eachindex(N_excs)])
	lens!(lens_pos[ind[3]][1][1], lens_pos[ind[3]][1][2], inset = (1, bbox(lens_pos[ind[3]][1][3], lens_pos[ind[3]][1][4], 0.4, 0.4)))
	lens!(lens_pos[ind[3]][2][1], lens_pos[ind[3]][2][2], inset = (1, bbox(lens_pos[ind[3]][2][3], lens_pos[ind[3]][2][4], 0.4, 0.4)))
	# g2: L
	p2 = plot(x_1.*γ_D, real(y_2),
		framestyle=:box,
		#xscale=:linear,
		# minorgrid=true,
		xlims=x_lims[ind[3]][2],
		lw=1.5,
		palette = cs_1[20:length(cs_1)÷8:length(cs_1)÷3],
		title="(b)",titleloc = :left, titlefont = font(14),
		xlabel=L"\gamma_D \tau",
		ylabel=L"g^{(2)}(\tau)",
		label=reshape([(i == 1) ? L"\mathrm{Tran}" : :none for i in eachindex(N_excs)], (1, length(N_excs))),
		tickfontsize=12,
		labelfontsize=14,
		legend_font=12,
		legend=:none,
	)
	plot!(x_1.*γ_D, real(y_21),
		lw=1.5,
		ls=:dashdot,
		palette = cs_2[20:length(cs_1)÷8:length(cs_1)÷3],
		label=reshape([(i == 1) ? L"\mathrm{Ref}" : :none for i in eachindex(N_excs)], (1, length(N_excs))),
	)
	annotate!(annotation_pos[ind[3]][3][1], annotation_pos[ind[3]][3][2], [text("Nₗ = "*"$(N_excs[i])", :center, 13, cs_1[(20:length(cs_1)÷8:length(cs_1)÷3)[i]]) for i in eachindex(N_excs)])
	annotate!(annotation_pos[ind[3]][4][1], annotation_pos[ind[3]][4][2], [text("Nₗ = "*"$(N_excs[i])", :center, 13, cs_2[(20:length(cs_1)÷8:length(cs_1)÷3)[i]]) for i in eachindex(N_excs)])
	lens!(lens_pos[ind[3]][3][1], lens_pos[ind[3]][3][2], inset = (1, bbox(lens_pos[ind[3]][3][3], lens_pos[ind[3]][3][4], (ind[3] == 1) ? 0.35 : 0.6, (ind[3] == 1) ? 0.4 : 0.6)))
	if ind[3] == 1
		lens!(lens_pos[ind[3]][4][1], lens_pos[ind[3]][4][2], inset = (1, bbox(lens_pos[ind[3]][4][3], lens_pos[ind[3]][4][4], 0.35, 0.4)))
	end
	# Spectrum: L
	p3 = plot([x_3[i][5000-200:5000+200] for i in eachindex(x_3)]./(2*pi*1e6), real([y_3[i][5000-200:5000+200] for i in eachindex(x_3)]) ./ y_3_max,
		framestyle=:box,
		#xscale=:linear,
		# minorgrid=true,
		xlims=x_lims[ind[3]][3],
		lw=1.5,
		palette = cs_1[20:length(cs_1)÷8:length(cs_1)÷3],
		title="(c)",titleloc = :left, titlefont = font(14),
		xlabel=L"\omega / 2\pi \; (\mathrm{MHz})",
		ylabel="Spectrum (normalized)",
		label=reshape([(i == 1) ? L"\mathrm{Transmitted}" : :none for i in eachindex(N_excs)], (1, length(N_excs))),
		tickfontsize=12,
		labelfontsize=14,
		legend_font=12,
		legend=:topright,
	)
	plot!([x_3[i][5000-200:5000+200] for i in eachindex(x_3)]./(2*pi*1e6), real([y_31[i][5000-200:5000+200] for i in eachindex(x_3)]) ./ y_3_max,
		lw=1.5,
		ls=:dashdot,
		palette = cs_2[20:length(cs_1)÷8:length(cs_1)÷3],
		label=reshape([(i == 1) ? L"\mathrm{Reflected}" : :none for i in eachindex(N_excs)], (1, length(N_excs))),
	)
	annotate!(annotation_pos[ind[3]][5][1], annotation_pos[ind[3]][5][2], [text("Nₗ = "*"$(N_excs[i])", :center, 13, cs_1[(20:length(cs_1)÷8:length(cs_1)÷3)[i]]) for i in eachindex(N_excs)])
	annotate!(annotation_pos[ind[3]][6][1], annotation_pos[ind[3]][6][2], [text("Nₗ = "*"$(N_excs[i])", :center, 13, cs_2[(20:length(cs_1)÷8:length(cs_1)÷3)[i]], rotation=(ind[3]==2) ? 70 : 80) for i in eachindex(N_excs)])
	p = plot(p1, p2, p3, layout=l, size=(1000,700))
	# savefig(p, PATH_FIGS*"fig_4_delta"*string(deltas["time"][ind[3]])*".pdf")
end

# ╔═╡ 2fcfdee4-6c7a-4c2d-957c-42759a63e153
# Fig 4_1: g2 for 
let
	# pgfplotsx()
	gr()
	ind = [1, 1, 2] # N_exc; DIRECTION: 1-L, 2-R; delta
	g2_start_tau = 2
	
	
	x_1 = [data_time[i, ind[2], 1, ind[3]]["taulist"][g2_start_tau:end] for i in eachindex(N_excs)]
	y_2 = [data_time[i, ind[2], 1, ind[3]]["g2_L"][g2_start_tau:end] for i in eachindex(N_excs)]
	y_21 = [data_time[i, ind[2], 1, ind[3]]["g2_R"][g2_start_tau:end] for i in eachindex(N_excs)]
	
	if ind[3] == 1 # delta = -0.15
		γ_D = mean(data_time[ind[1], ind[2], 1, ind[3]]["g_D(L_j)"][255:256])
	elseif ind[3] == 2 # delta = -0.015
		γ_D = data_time[ind[1], ind[2], 1, ind[3]]["g_D(L_j)"][251]
	end

	annotation_pos = [[
					   [[4.65, 4.65, 4.65], [0.62, 0.53, 0.45]],
					   [[4, 4, 4], [0.12, 0.32, 0.4]],
					   [[0.4, 0.9, 0.9], [1.1, 2.2, 2.08]],
					   [[1.0, 1.87, 1.87], [0.75, 2.4, 2.53]],
					   [[0.0, 0.0, 0.0], [0.3, 0.42, 0.48]],
					   [[-3.4, -1.7, -0.5], [0.9, 0.75, 0.6]]
					  ],
					  [
					   [[1.7, 1.7, 1.7], [0.69, 0.96, 0.9]],
					   [[1.7, 1.7, 1.7], [0.37, 0.43, 0.48]],
					   [[0.9, 0.9, 0.9], [1.95, 1.65, 1.8]],
					   [[0.5, 0.9, 0.9], [1.35, 2.3, 2.15]],
					   [[1.4, 1.4, 1.4], [0.25, 0.38, 0.44]],
					   [[-0.1, -0.05, -0.05], [0.28, 0.37, 0.19]]  
					  ]
					 ]
	lens_pos = [
		[
			[[1.2, 3.5], [0.484, 0.497],0.58, 0.02],
			[[1, 3.5], [0.35, 0.37],0.16, 0.02],
			[[0.45, 0.55], [0.83, 0.84],0.12, 0.05],
			[[0.45, 0.55], [0.73, 0.74],0.54, 0.05]
		],
		[
			[[1.5, 1.7], [0.64, 0.645],0.38, 0.05],
			[[1.5, 1.7], [0.62, 0.63],0.38, 0.52],
			[[0.6, 0.7], [1.005, 1.02],0.25, 0.05],
			[[0.45, 0.55], [0.73, 0.74],0.54, 0.05]
		]
	]
					
	x_lims = [[(-0.1,5),
			  (-0.05,2),
			  (-20,20)],
			  [(-0.05,2),
			  (-0.025,1),
			  (-0.2,0.2)]]
	
	cs_1 = cgrad(:diverging_bwr_20_95_c54_n256)
	cs_2 = cgrad(:diverging_bwr_20_95_c54_n256, rev=true)
	l = @layout [grid(2,1) b{0.5w}] 
	# g2: L
	p2 = plot(x_1.*γ_D, real(y_2),
		framestyle=:box,
		#xscale=:linear,
		# minorgrid=true,
		xlims=x_lims[ind[3]][2],
		lw=1.5,
		palette = cs_1[20:length(cs_1)÷8:length(cs_1)÷3],
		title=L"\Delta L_J = 10^{-3}",titleloc = :center, titlefont = font(14),
		xlabel=L"\gamma_D \tau",
		ylabel=L"g^{(2)}(\tau)",
		label=reshape([(i == 1) ? L"\mathrm{Tran}" : :none for i in eachindex(N_excs)], (1, length(N_excs))),
		tickfontsize=12,
		labelfontsize=14,
		legend_font=12,
		legend=:none,
	)
	plot!(x_1.*γ_D, real(y_21),
		lw=1.5,
		ls=:dashdot,
		palette = cs_2[20:length(cs_1)÷8:length(cs_1)÷3],
		label=reshape([(i == 1) ? L"\mathrm{Ref}" : :none for i in eachindex(N_excs)], (1, length(N_excs))),
	)
	annotate!(annotation_pos[ind[3]][3][1], annotation_pos[ind[3]][3][2], [text("Nₗ = "*"$(N_excs[i])", :center, 13, cs_1[(20:length(cs_1)÷8:length(cs_1)÷3)[i]]) for i in eachindex(N_excs)])
	annotate!(annotation_pos[ind[3]][4][1], annotation_pos[ind[3]][4][2], [text("Nₗ = "*"$(N_excs[i])", :center, 13, cs_2[(20:length(cs_1)÷8:length(cs_1)÷3)[i]]) for i in eachindex(N_excs)])
	lens!(lens_pos[ind[3]][3][1], lens_pos[ind[3]][3][2], inset = (1, bbox(lens_pos[ind[3]][3][3], lens_pos[ind[3]][3][4], (ind[3] == 1) ? 0.35 : 0.6, (ind[3] == 1) ? 0.4 : 0.6)))
	if ind[3] == 1
		lens!(lens_pos[ind[3]][4][1], lens_pos[ind[3]][4][2], inset = (1, bbox(lens_pos[ind[3]][4][3], lens_pos[ind[3]][4][4], 0.35, 0.4)))
	end
	p = plot(p2, size=(500,350))
	# savefig(p, PATH_FIGS*"fig_4_1.pdf")
end

# ╔═╡ f38934e1-d100-4ce4-9488-f9cda003d24f
# Fig 4: g1, g2, spectra for right excitation
let
	# pgfplotsx()
	gr()
	ind = [1, 2, 1] # N_exc; DIRECTION: 1-L, 2-R; delta
	x = data_time[ind[1], ind[2], 1, ind[3]]["taulist"]
	y_2 = [data_time[i, ind[2], 1, ind[3]]["g1_R"] for i in eachindex(N_excs)]
	y_21 = [data_time[i, ind[2], 1, ind[3]]["g1_L"] for i in eachindex(N_excs)]
	
	if ind[3] == 1 # delta = -0.15
		γ_D = mean(data_time[ind[1], ind[2], 1, ind[3]]["g_D(L_j)"][255:256])
	elseif ind[3] == 2 # delta = -0.015
		γ_D = data_time[ind[1], ind[2], 1, ind[3]]["g_D(L_j)"][251]
	end
	
	cs_1 = cgrad(:diverging_bwr_20_95_c54_n256)
	cs_2 = cgrad(:diverging_bwr_20_95_c54_n256, rev=true)
	l = @layout grid(2,2)
	# g1: R
	p2 = plot(x.*γ_D, real(y_2),
		framestyle=:box,
		#xscale=:linear,
		# minorgrid=true,
		xlims=(-0.1,5),
		lw=1.5,
		palette = cs_1[20:length(cs_1)÷8:length(cs_1)÷3],
		title="(b)",titleloc = :left, titlefont = font(14),
		xlabel=L"\gamma_D \tau",
		ylabel=L"g^{(1)}(\tau)",
		label=reshape([(i == 1) ? L"\mathrm{Tran}" : :none for i in eachindex(N_excs)], (1, length(N_excs))),
		tickfontsize=12,
		labelfontsize=14,
		legend_font=12,
		#legend=:inside,
	)
	plot!(x.*γ_D, real(y_21),
		lw=1.5,
		ls=:dashdot,
		palette = cs_2[20:length(cs_1)÷8:length(cs_1)÷3],
		label=reshape([(i == 1) ? L"\mathrm{Ref}" : :none for i in eachindex(N_excs)], (1, length(N_excs))),
	)
	annotate!([3.5, 2.3, 2.3], [0.46, 0.82, 0.78], [text("Nₗ = "*"$(N_excs[i])", :center, 13, cs_1[(20:length(cs_1)÷8:length(cs_1)÷3)[i]]) for i in eachindex(N_excs)])
	annotate!([4.7, 4.7, 4.7], [0.73, 0.89, 0.85], [text("Nₗ = "*"$(N_excs[i])", :center, 13, cs_2[(20:length(cs_1)÷8:length(cs_1)÷3)[i]]) for i in eachindex(N_excs)])
	lens!([3, 4], [0.9888, 0.9912], inset = (1, bbox(0.58, 0.2, 0.3, 0.4)))
	lens!([3, 4], [0.478, 0.481], inset = (1, bbox(0.18, 0.2, 0.3, 0.4)))
end

# ╔═╡ 0b808b21-23da-4f45-8fee-3af0f549cf70
cgrad(:diverging_bwr_20_95_c54_n256, rev=true)[20:length(cgrad(:diverging_bwr_20_95_c54_n256, rev=true))÷8:length(cgrad(:diverging_bwr_20_95_c54_n256, rev=true))÷2]

# ╔═╡ 18d23ed4-5eff-45e7-be12-e83024dd9bf0
md"
## $\rho_{ss}$ depending on power and $N_\mathrm{exc}$
"

# ╔═╡ 17b3b2ae-5040-4eb3-9978-3f97bdeb7701
data_power_Nexc[1]["a_inc_list"]

# ╔═╡ 5c4f9766-9aba-4907-bde9-2af4c25b4d39
transition_rates = [0.000000e+00 	0.000000e+00 	0.000000e+00 	0.000000e+00 	0.000000e+00 	0.000000e+00;
5.521305e+04 	5.521305e+04 	0.000000e+00 	0.000000e+00 	0.000000e+00 	0.000000e+00;
1.023906e+09 	2.526576e-09 	1.023906e+09 	0.000000e+00 	0.000000e+00 	0.000000e+00;
0.000000e+00 	5.523656e+04 	1.023906e+09 	1.023961e+09 	0.000000e+00 	0.000000e+00;
0.000000e+00 	1.023906e+09 	5.526594e+04 	0.000000e+00 	1.023962e+09 	0.000000e+00;
0.000000e+00 	5.521893e+04 	1.023906e+09 	0.000000e+00 	5.960464e-08 	1.023961e+09]

# ╔═╡ 211c9bdc-31bb-45cd-b120-ff9a754afccc
@bind a_slider Slider(1:length(data_power_Nexc[1]["a_inc_list"]))

# ╔═╡ 308435df-b47b-416b-b141-81050a1a9f6b
@bind N_exc Slider(2:8)

# ╔═╡ 59726ecd-f3e2-4a94-8a70-8765f3d7e63a
@bind delt_slider Slider(1:length(deltas["power"]))

# ╔═╡ 714de8ab-1422-4f6e-9b2d-a92262ac91d7
print("a_slider = ", a_slider, "; a_in^2/γ = ", round(data_power_Nexc[delt_slider]["a_inc_list"][a_slider].^2 ./ mean(data_power_Nexc[delt_slider]["gammas"]), digits=2), "; N_exc = ", N_exc, "; δ = ", deltas["power"][delt_slider])

# ╔═╡ 5765f3c8-7b6d-465a-bd5e-58bb314edfeb
# Plot reference fig
let
	gr()
	x = data_power_Nexc[delt_slider]["a_inc_list"].^2 ./ mean(data_power_Nexc[delt_slider]["gammas"])
	y_1 = data_power_Nexc[delt_slider]["transmission_R_"*string(N_exc)]
	y_2 = data_power_Nexc[delt_slider]["transmission_L_"*string(N_exc)]

	y_3 = diag(real(data_power_Nexc[delt_slider]["rho_ss_DB_L_"*string(N_exc)][:,:,a_slider]))

	cs_1 = cgrad(:diverging_bwr_20_95_c54_n256)
	cs_2 = cgrad(:diverging_bwr_20_95_c54_n256, rev=true)

	kw = (
		grid = true, 
		minorgrid = true,
		#tickfontfamily=:sanserif,
		#fontfamily=:sanserif,
		tickfontsize=12,
		labelfontsize=14,
		legend_font=12,
		framestyle=:box,
	)
	
	p_1 = plot(x, y_1,
			 lw=2,
			 xlims=(1e-5, x[end]),
		     xscale=:log10,
			 xlabel=L"|a_\mathrm{in}|^2 / \bar{\gamma}",
			 ylabel=L"\mathrm{Transmission}",
			 label=L"R",
			 palette = cs_2[20:length(cs_1)÷8:length(cs_1)÷3];
			 kw...)
	plot!(x, y_2,
		  lw=2,
		  xlabel=L"|a_\mathrm{in}|^2 / \bar{\gamma}",
		  ylabel=L"\mathrm{Transmission}",
		  label=L"L",
		  palette = cs_1[20:length(cs_1)÷8:length(cs_1)÷3];
		  kw...)
	vline!([x[a_slider]], lw=1, linestyle=:dash, color=:black, label=L"\mathrm{slice}")
	p_2 = plot(bar(1:N_exc^2, y_3),
			   ylims=(0,0.01),
			   legend=false,
			   xlabel="State",
			   ylabel="Population";
			   kw...)
	p = plot(p_1, p_2, layout=(1,2))
end

# ╔═╡ abcf719d-04d5-4ac5-aedf-10555d2687be
# Plot reference fig
let
	gr()
	# pgfplotsx()
	
	a_in = 371
	N_exc_list = [2, 3, 5]
	delt_i = 2

	y_1 = diag(real(data_power_Nexc[delt_i]["rho_ss_DB_L_"*string(N_exc_list[1])][:,:,a_in]))
	y_2 = diag(real(data_power_Nexc[delt_i]["rho_ss_DB_L_"*string(N_exc_list[2])][:,:,a_in]))
	y_3 = diag(real(data_power_Nexc[delt_i]["rho_ss_DB_L_"*string(N_exc_list[3])][:,:,a_in]))

	cs_1 = cgrad(:diverging_bwr_20_95_c54_n256)
	cs_2 = cgrad(:diverging_bwr_20_95_c54_n256, rev=true)
	
	kw = (
		grid = true, 
		minorgrid = true,
		#tickfontfamily=:sanserif,
		#fontfamily=:sanserif,
		tickfontsize=12,
		labelfontsize=14,
		legend_font=12,
		framestyle=:box,
	)
	
	p_1 = plot(bar(1:N_exc_list[1]^2, y_1,
			   	   color=cs_1[1]),
			   xticks=([1, 2, N_exc_list[1] + 1, N_exc_list[1]^2], 
			   		   [L"|G\rangle", L"|D\rangle", L"|B\rangle", L"|E\rangle"]),
			   ylims=(0,0.61),
			   legend=false,
			   xlabel="State",
			   ylabel="Population",
			   title=L"N_l = 2";
			   kw...)
	lens!([2.5, 4.5], [0.0, 0.025], inset = (1, bbox(0.5, 0.05, 0.45, 0.6)))
	p_2 = plot(bar(1:N_exc_list[2]^2, y_2,
			   	   color=cs_1[20],),
			   xticks=([1, 2, N_exc_list[2] + 1, N_exc_list[2]+2, N_exc_list[2]*2+1], 
			   		   [L"|G\rangle", L"|D\rangle", L"|B\rangle", L"|E\rangle", L"|A\rangle"]),
			   ylims=(0,0.61),
			   legend=false,
			   xlabel="State",
			   ylabel="Population",
			   title=L"N_l = 3";
			   kw...)
	lens!([2.5, 9.5], [0.0, 0.025], inset = (1, bbox(0.4, 0.05, 0.55, 0.6)))
	p_3 = plot(bar(1:N_exc_list[3]^2, y_3,
			   	   color=cs_1[80],),
			   xticks=([1, 2, N_exc_list[3] + 1, N_exc_list[3]+2, N_exc_list[3]*2+1], 
			   		   [L"|G\rangle", L"|D\rangle", L"|B\rangle", L"|E\rangle", L"|A\rangle"]),
		       xlims=(0, 12),
			   ylims=(0,0.61),
			   legend=false,
			   xlabel="State",
			   ylabel="Population",
			   title=L"N_l = 5";
			   kw...)
	lens!([2.5, 12.5], [0.0, 0.025], inset = (1, bbox(0.4, 0.05, 0.55, 0.6)))
	p = plot(p_1, p_2, p_3, layout=(1,3), size=(1000, 300))
	# savefig(p, PATH_FIGS*"fig_7.pdf")
end

# ╔═╡ b8451879-de3f-4b3a-8402-358543c2597c
data_states

# ╔═╡ 7c2d7a30-9937-4594-af63-5994935a696d
let
	gr()
	# pgfplotsx()

	tran_rates = data_states["tran_rates_L"]
	γ_D = data_states["gamma_D"]
	γ_B = data_states["gamma_B"]
	γ = mean(data_power_Nexc[2]["gammas"])
	t_list = data_states["tlist"]
	gg = data_states["gg"] .+ 1e-6
	dd = data_states["dd"] .+ 1e-6
	bb = data_states["bb"] .+ 1e-6
	ee = data_states["ee"] .+ 1e-6
	bb2 = data_states["bb2"] .+ 1e-6
	dd2 = data_states["dd2"] .+ 1e-6
	
	labels = [L"|G\rangle" L"|D\rangle" L"|B\rangle" L"|E\rangle" L"|A\rangle" L"|S\rangle"]

	l = @layout [a{0.4w} b{0.6w}]
	
	kw = (
		grid = true, 
		minorgrid = false,
		#tickfontfamily=:sanserif,
		#fontfamily=:sanserif,
		tickfontsize=12,
		labelfontsize=14,
		legend_font=12,
		# framestyle=:box,
	)
	kw_1 = (
		cbar_title = L"\log_{10}(\gamma_{ij} / \bar{\gamma})", 
		colorbar_titlefontsize=14, 
		colorbar_tickfontsize=12,
	)
	
	p_1 = heatmap(log10.(tran_rates / γ),
        xticks=(1:6, labels),  # Set the x-axis tick labels
        yticks=(1:6, labels),  # Set the y-axis tick labels
        aspect_ratio=:equal,  # Make sure the aspect ratio is equal for square cells
		title="(a)    Transition rates",titleloc = :left, titlefont = font(14),
        c=:vik;
		kw..., kw_1...
	) 

	p_2 = plot(t_list*γ_D, [gg, dd, bb, ee, bb2, dd2],
		lw=2,
		linestyle=:solid,
		xlims=(-2e-2, 1e0),
		yscale=:log10,
		yticks=[1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0],
		title="(b)",titleloc = :left, titlefont = font(14),
		xlabel=L"\gamma_D t",
		ylabel=L"\mathrm{Population}",
		framestyle=:box,
		left_margin=(25,:mm),
		palette=:darktest,
		label=labels;
		# color = cs_2[20];
		kw...
	)
	p = plot(p_1, p_2, layout=l, size=(800, 300))
	# savefig(p, PATH_FIGS*"fig_8.pdf")
end

# ╔═╡ 348bd3e6-c6f1-46b9-b162-58557ab9e224
cgrad(:diverging_bwr_20_95_c54_n256)[1]

# ╔═╡ 1ba194c4-d4c1-4873-857d-12b1c07c81bc
md"
## Meanfield approximation
"

# ╔═╡ 730897a2-cbe7-480e-a6ae-65261541aa06
data_meanfield[2, 1]

# ╔═╡ 62803a99-78e3-4663-986a-77c04cc99eba
data_meanfield_small[1]["gammas"]

# ╔═╡ 6fa23358-8fe7-47bf-90c3-2fa4321a25a9
data_time[1,1,1,1]["L_j_list"]

# ╔═╡ c81f01ce-f8f5-4309-b977-0dc734eec48e
mean(data_time[1,1,1,1]["g_D(L_j)"][255:256])

# ╔═╡ 8f5fcba4-fcb4-4507-9560-7f7f551978c8
deltas["meanfield"]

# ╔═╡ 9a8c9810-3bf4-4a14-a53c-d4bc0600bfd8
# Fig 6: g1, g2, spectra
let
	# pgfplotsx()
	gr()
	ind = 1 # delta: 1 -> -0.15, 2 -> -0.015, 3 -> -0.0015

	# parameters
	run(`bash run_circ_parameters.sh 3.31e-9 3.3e-9`)  # choose the right L_j1
	params = load("../Data/params.h5")
	
	DIRECTION = "L"
	a_mag = 0.1*sqrt(1.5e9)
	a_inc = Dict("L" => DIRECTION == "L" ? a_mag : 0.0,
				 "R" => DIRECTION == "R" ? a_mag : 0.0)
	a_mag_s = 0.01*sqrt(1.5e9)
	a_inc_s = Dict("L" => DIRECTION == "L" ? a_mag_s : 0.0,
				   "R" => DIRECTION == "R" ? a_mag_s : 0.0)

	g_nr_ = params["gamma_nr"]
	g_phi_ = params["gamma_phi"]
	w_ = params["w_m"]
	A_ = params["A_m"]
	Γ_ = params["Gamma"]
	Ω_ = params["Omega"]
	w_ext_ = params["w_ext"]
	t_j_ = params["t_j"]
	g_ = 1.0im*sqrt.(0.5*params["gammas"]).*(a_inc["R"] * exp.(-im*w_ext_*t_j_) .+
											 a_inc["L"] * exp.(im*w_ext_*t_j_))
	γ = mean(data_meanfield[2, ind]["gammas"])
	γ_D = minimum(eigen(Γ_).values)

	
	# variables
	x_1 = data_meanfield[1, ind]["a_in_list"].^2 / γ
	x_2 = data_meanfield[1, ind]["t_list"] * γ_D
	x_2_s = data_meanfield_small[1, 1]["t_list"] * γ_D
	
	y_11 = [data_meanfield[i, ind]["transmission_R"] for i in eachindex(orders)]
	y_12 = [data_meanfield[i, ind]["transmission_L"] for i in eachindex(orders)]

	a_1 = [data_meanfield[i, ind]["a(1)"] for i in eachindex(orders)]
	a_2 = [data_meanfield[i, ind]["a(2)"] for i in eachindex(orders)]
	a_1_s = [data_meanfield_small[i]["a(1)"] for i in eachindex(orders)]
	a_2_s = [data_meanfield_small[i]["a(2)"] for i in eachindex(orders)]

	y_2R = [abs.(a_inc["R"] .+ exp(-im*w_[1]*t_j_[1]) * 
			 sqrt(0.5*params["gammas"][1]) * a_1[i] .+ exp(-im*w_[2]*t_j_[2]) * 
			 sqrt(0.5*params["gammas"][2]) * a_2[i]) / a_mag for i in eachindex(orders)]
	y_2L = [abs.(a_inc["L"] .+ exp(im*w_[1]*t_j_[1]) * 
			 sqrt(0.5*params["gammas"][1]) * a_1[i] .+ exp(im*w_[2]*t_j_[2]) * 
			 sqrt(0.5*params["gammas"][2]) * a_2[i]) / a_mag for i in eachindex(orders)]

	y_3 = abs.(data_meanfield[2, ind]["a(1)'a(2)"])
	y_3_s = abs.(data_meanfield_small[2, 1]["a(1)'a(2)"])
	y_3_s0 = abs.(data_meanfield_small[2, 2]["a(1)'a(2)"])

	# plot params
	cs_1 = cgrad(:diverging_bwr_20_95_c54_n256)
	cs_2 = cgrad(:diverging_bwr_20_95_c54_n256, rev=true)
	l = @layout [grid(2,1) b{0.5w}] 

	kw = (
		grid = true, 
		minorgrid = true,
		#tickfontfamily=:sanserif,
		#fontfamily=:sanserif,
		tickfontsize=12,
		labelfontsize=14,
		legend_font=12,
		framestyle=:box,
	)
	
	# Power dependence
	p_1 = plot(x_1, y_11,
			 lw=[1.5 2],
			 linestyle=[:dashdot :solid],
			 xlims=(1e-5, x_1[end]),
		     xscale=:log10,
		     title="(a)",titleloc = :left, titlefont = font(14),
			 xlabel=L"|a_\mathrm{in}|^2 / \bar{\gamma}",
			 ylabel=L"\mathrm{Transmission}",
			 label=[:none L"R"],
			 color = cs_2[20];
			 kw...)
	plot!(x_1, y_12,
		  lw=[1.5 2],
		  linestyle=[:dashdot :solid],
		  label=[:none L"L"],
		  color = cs_1[20];
		  kw...)
	vline!([x_1[32]], lw=1, linestyle=:dash, color=:black, label=:none)#L"\mathrm{slice}")
	annotate!([2e-3, 1.1e1], [0.65, 0.05], [text("with correlations", :center, 13), text("without correlations", :center, 13)])

	# Time dependence
	p_2 = plot(x_2, y_2R,
			 lw=[1.5 2],
			 linestyle=[:dashdot :solid],
			 xlims=(0, 2),
		     title="(b)",titleloc = :left, titlefont = font(14),
			 xlabel=L"\gamma_D t",
			 ylabel=L"|a_\mathrm{out} / a_\mathrm{in}|",
			 label=[:none L"R"],
			 legend=false,
			 color = cs_2[20];
			 kw...)
	plot!(x_2, y_2L,
		 lw=[1.5 2],
		 linestyle=[:dashdot :solid],
		 # xlims=(1e-5, x_1[end]),
		 label=[:none L"L"],
		 color = cs_1[20];
		 kw...)

	p_3 = plot(x_2, y_3,
		 lw=2,
		 linestyle=:solid,
		 # yscale=:log10,
		 # xlims=(3, 5),
		 # ylims=(-0.01, 0.05),
		 title="(c)",titleloc = :left, titlefont = font(14),
		 xlabel=L"\gamma_D t",
		 ylabel=L"|\langle \sigma_1^+ \sigma_2 \rangle|",
		 label=[:none L"R"],
		 color = :black;
		 kw...)
	plot!(x_2_s, y_3_s,
		 lw=2,
		 linestyle=:dash,
		 label=[:none L"L"],
		 color = :black;
		 kw...)
	plot!(x_2_s, y_3_s0,
		 lw=2,
		 linestyle=:dot,
		 label=[:none L"L"],
		 color = :black;
		 kw...)
	annotate!([3.5, 3.5, 3.5], [0.27, 0.04, 0.01], [text(L"|a_\mathrm{in}|^2 / \bar{\gamma} = 3\cdot 10^{-2}", :center, 13), text(L"|a_\mathrm{in}|^2 / \bar{\gamma} = 3\cdot 10^{-4}", :center, 13), text(L"|a_\mathrm{in}|^2 / \bar{\gamma} = 3\cdot 10^{-6}", :center, 13)])
	p = plot(p_1, p_2, p_3, layout=l, size=(1000, 500))
	# savefig(p, PATH_FIGS*"fig_6.pdf")
end

# ╔═╡ 563d4760-e924-41c3-abfb-08e984815600
md"
## Fig 2: gammas and eigenlevels
"

# ╔═╡ 2b93634c-315f-42c5-ba3c-9c4083b604ca
data_time[1,1,1,1]

# ╔═╡ be108252-cdc2-4a6d-ad44-3c07051842cb
cgrad(:linear_blue_5_95_c73_n256)[20:9:256]

# ╔═╡ 970e10c5-5894-490b-8972-b69c3365c59f
# fig 2
let
	# pgfplotsx()
	gr()
	# N_exc; DIRECTION: 1-L, 2-R; a_inc: 0.03; delta
	ind = 1 # delta: 1 -> -0.15, 2 -> -0.015, 3 -> -0.0015

	# parameters
	γ = mean(data_meanfield[2, ind]["gammas"])
	
	# variables
	x_1 = data_time[1,1,1,1]["L_j_list"] / 1e-9
	x_2 = range(-pi, pi, 200)
	y_11 = data_time[1,1,1,1]["g_B(L_j)"] / γ
	y_12 = data_time[1,1,1,1]["g_D(L_j)"] / γ

	y_2 = data_time[1,1,1,1]["H_eig(L_j)"] / (2*pi*1e9)
	y_3 = data_time[2,1,1,1]["H_eig(L_j)"] / (2*pi*1e9)
	y_4 = real.(data_time[3,1,1,1]["H_eig(L_j)"] / (2*pi*1e9))
	

	# plot
	cs_1 = cgrad(:linear_blue_5_95_c73_n256)
	cs_2 = cgrad(:linear_blue_5_95_c73_n256, rev=true)
	l = @layout([[a; b] grid(3, 1)])
	
	kw = (
		grid = true, 
		minorgrid = true,
		tickfontsize=12,
		labelfontsize=14,
		legend_font=12,
		framestyle=:box,
	)

	p_1 = plot(x_1, y_11,
		lw=2,
		linestyle=:solid,
		xlims=(2.95, 3.65),
		# xscale=:log10,
		title="(b)",titleloc = :left, titlefont = font(14),
		xlabel=L"L_j \;\; \mathrm{(nH)}",
		ylabel=L"\gamma_{B,D} / \bar{\gamma}",
		label=L"\gamma_B",
		color = cs_2[20],
		right_margin=(6,:mm),
		legend=:topright;
		kw...)
	plot!(x_1, y_12,
		lw=2,
		linestyle=:solid,
		label=L"\gamma_D",
		color = cs_1[20],
	)

	p_2 = plot(x_1, y_2',
		lw=2,
		xlims=(3.15, 3.45),
		ylims=(-0.25, 0.25),
		# xscale=:log10,
		title=L"\mathrm{(c)} \;\; N_{l} = 2",titleloc = :left, titlefont = font(14),
		xlabel=L"L_j \;\; \mathrm{(nH)}",
		ylabel=L"\omega / 2\pi \;\; \mathrm{(GHz)}",
		label=:none,
		palette=cs_1[20:70:256];
		kw...)
	p_3 = plot(x_1, y_3',
		lw=2,
		xlims=(3.15, 3.45),
		ylims=(-0.5, 0.25),
		# xscale=:log10,
		title=L"\mathrm{(d)} \;\; N_{l} = 3",titleloc = :left, titlefont = font(14),
		xlabel=L"L_j \;\; \mathrm{(nH)}",
		ylabel=L"\omega / 2\pi \;\; \mathrm{(GHz)}",
		label=:none,
		palette=cs_1[20:27:256];
		kw...)
	p_4 = plot(x_1, y_4',
		lw=2,
		xlims=(3.15, 3.45),
		ylims=(-2.0, 0.25),
		# xscale=:log10,
		title=L"\mathrm{(e)} \;\; N_{l} = 5",titleloc = :left, titlefont = font(14),
		xlabel=L"L_j \;\; \mathrm{(nH)}",
		ylabel=L"\omega / 2\pi \;\; \mathrm{(GHz)}",
		label=:none,
		palette=cs_1[20:9:256];
		kw...)

	f = x -> 4.5*sin(x/2)^2
	g = x -> 4.5*(x/2)^2
	p_5 = plot(x_2, f.(x_2),
		lw=2,
		# xlims=(3.15, 3.45),
		ylims=(-0.1, 5.1),
		title="(a)",titleloc = :left, titlefont = font(14),
		xlabel=L"\mathrm{Superconducting} \;\; \mathrm{phase}, \; \phi",
		ylabel=L"\mathrm{Energy}, \;\; (\hbar \omega_{01})",
		label=:none,
		xticks=([-pi, -pi/2, 0, pi/2, pi], [L"-\pi", L"-\pi/2", L"0", L"\pi/2", L"\pi"]),
		palette=cs_1[20:27:256];
		kw...)
	plot!(x_2, g.(x_2),
		lw=1,
		linestyle=:dash,
		color=:grey,
		label=:none)
	levels = [0.5, 1.5, 2.4, 3.1, 3.7]
	plot!([[-2*asin(sqrt(2/9*i)), 2*asin(sqrt(2/9*i))] for i in levels], 
		[[i, i] for i in levels],
		label=:none)
	plot!([0.5, 0.5], [[levels[1], levels[2]], [levels[2], levels[3]]],
		label=:none, color=:black, linestyle=:dot)
	scatter!([0.5, 0.5], [levels[1]+0.07, levels[2]+0.07],
		markershape=:dtriangle, color=:black, label=:none)
	scatter!([0.5, 0.5], [levels[2]-0.07, levels[3]-0.07],
		markershape=:utriangle, color=:black, label=:none)
	annotate!([1.1, 1.6, 2.1], levels[1:3], [text(L"|0\rangle", :center, 18), text(L"|1\rangle", :center, 18), text(L"|2\rangle", :center, 18)])
	annotate!([0, 0], levels[1:2] .+ 0.5, [text(L"\hbar \omega_{01}", :center, 18), text(L"\hbar \omega_{12}", :center, 18)])
	
	p = plot(p_5, p_1, p_2, p_3, p_4, layout=l, size=(800, 1000))
	# savefig(p, PATH_FIGS*"fig_2.pdf")
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Contour = "d38c429a-6771-53c6-b99e-75d170b6e991"
FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
Gaston = "4b11ee91-296f-5714-9832-002c20994614"
HDF5 = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PGFPlotsX = "8314cec4-20b6-5062-9cdb-752b83310925"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
Contour = "~0.6.2"
FileIO = "~1.16.1"
Gaston = "~1.1.0"
HDF5 = "~0.17.0"
LaTeXStrings = "~1.3.0"
PGFPlotsX = "~1.6.0"
Plots = "~1.39.0"
PlutoUI = "~0.7.52"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.3"
manifest_format = "2.0"
project_hash = "132bbda24298e5ca25e96995b155be5cc11d7a8f"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "91bd53c39b9cbfb5ef4b015e8b582d344532bd0a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.0"

[[deps.ArgCheck]]
git-tree-sha1 = "a3a402a35a2f7e0b87828ccabbd5ebfbebe356b4"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.3.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "02aa26a4cf76381be7f66e020a3eddeb27b0a092"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "d9a8f86737b665e15a9641ecbac64deef9ce6724"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.23.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "e460f044ca8b99be31d35fe54fc33a5c33dd8ed7"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.9.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "5372dbbf8f0bdb8c700db5367132925c0771ef7e"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.2.1"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.DataAPI]]
git-tree-sha1 = "8da84edb865b0b5b0100c0666a9bc9a0b71c553c"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.15.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3dbd312d370723b6bb43ba9d02fc36abade4518d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.15"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DefaultApplication]]
deps = ["InteractiveUtils"]
git-tree-sha1 = "c0dfa5a35710a193d83f03124356eef3386688fc"
uuid = "3f0dd361-4fe0-5fc6-8523-80b14ec94d85"
version = "1.1.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "e90caa41f5a86296e014e148ee061bd6c3edec96"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.9"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4558ab818dcceaab612d1bb8c19cee87eda2b83c"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.5.0+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "299dc33549f68299137e51e6d49a13b5b1da9673"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d8db6a5a2fe1381c1ea4ef2cab7c69c2de7f9ea0"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.1+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "d73afa4a2bb9de56077242d98cf763074ab9a970"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.72.9"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "1596bab77f4f073a14c62424283e7ebff3072eca"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.72.9+1"

[[deps.Gaston]]
deps = ["ColorSchemes", "DelimitedFiles", "Random"]
git-tree-sha1 = "843b5df546c02aa77880d8f9cfb65063a8938b0f"
uuid = "4b11ee91-296f-5714-9832-002c20994614"
version = "1.1.0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "e94c92c7bf4819685eb80186d51c43e71d4afa17"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.76.5+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HDF5]]
deps = ["Compat", "HDF5_jll", "Libdl", "MPIPreferences", "Mmap", "Preferences", "Printf", "Random", "Requires", "UUIDs"]
git-tree-sha1 = "ec7df74b7b2022e8252a8bfd4ec23411491adc3b"
uuid = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
version = "0.17.0"

    [deps.HDF5.extensions]
    MPIExt = "MPI"

    [deps.HDF5.weakdeps]
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"

[[deps.HDF5_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "Libdl", "OpenSSL_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "4cc2bb72df6ff40b055295fdef6d92955f9dede8"
uuid = "0234f1f7-429e-5d53-9886-15a909be8d59"
version = "1.12.2+2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "5eab648309e2e060198b45820af1a37182de3cce"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "d75853a0bdbfb1ac815478bacd89cd27b550ace6"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.3"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "f377670cda23b6b7c1c0b3893e37451c5c1a2185"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.5"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6f2675ef130a300a112286de91973805fcc5ffbc"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.91+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f689897ccbe049adb19a065c495e75f372ecd42b"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.4+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "f428ae552340899a935973270b8d98e5a31c49fe"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.1"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "7d6dd4e9212aebaeed356de34ccf262a3cd415aa"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.26"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "0d097476b6c381ab7906460ef1ef1638fbce1d91"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.2"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MPIPreferences]]
deps = ["Libdl", "Preferences"]
git-tree-sha1 = "781916a2ebf2841467cda03b6f1af43e23839d85"
uuid = "3da0fdf6-3ccc-4f1b-acd9-58baa6c99267"
version = "0.1.9"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9ee1618cbf5240e6d4e0371d6f24065083f60c48"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.11"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "51901a49222b09e3743c65b8847687ae5fc78eb2"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a12e56c72edee3ce6b96667745e6cbbe5498f200"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.23+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "2e73fe17cac3c62ad1aebe70d44c963c3cfdc3e3"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.2"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+0"

[[deps.PGFPlotsX]]
deps = ["ArgCheck", "Dates", "DefaultApplication", "DocStringExtensions", "MacroTools", "OrderedCollections", "Parameters", "Requires", "Tables"]
git-tree-sha1 = "3e7a0345b9f37da2cd770a5d47bb5cb6e62c7a81"
uuid = "8314cec4-20b6-5062-9cdb-752b83310925"
version = "1.6.0"

    [deps.PGFPlotsX.extensions]
    ColorsExt = "Colors"
    ContourExt = "Contour"
    MeasurementsExt = "Measurements"
    StatsBaseExt = "StatsBase"

    [deps.PGFPlotsX.weakdeps]
    Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
    Contour = "d38c429a-6771-53c6-b99e-75d170b6e991"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "716e24b21538abc91f6205fd1d8363f39b442851"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.7.2"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "64779bc4c9784fee475689a1752ef4d5747c5e87"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.42.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "f92e1315dadf8c46561fb9396e525f7200cdc227"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.5"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "ccee59c6e48e6f2edf8a5b64dc817b6729f99eb5"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.39.0"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e47cd150dbe0443c3a3651bc5b9cbd5576ab75b7"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.52"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "364898e8f13f7eaaceec55fd3d08680498c0aa6e"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.4.2+3"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "30449ee12237627992a99d5e30ae63e4d78cd24a"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "c60ec5c62180f27efea3ba2908480f8055e17cee"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "a1f34829d5ac0ef499f6d84428bd6b4c71f02ead"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.11.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "9a6ae7ed916312b41236fcef7e0af564ef934769"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.13"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.URIs]]
git-tree-sha1 = "b7a5e99f24892b6824a954199a45e9ffcc1c70f0"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "a72d22c7e13fe2de562feda8645aa134712a87ee"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.17.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "e2d817cc500e960fdbafcf988ac8436ba3208bfd"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.3"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "04a51d15436a572301b5abbb9d099713327e9fc4"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.4+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cf2c7de82431ca6f39250d2fc4aacd0daa1675c0"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.4+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "868e669ccb12ba16eaf50cb2957ee2ff61261c56"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.29.0+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╟─50b414b9-bda7-410c-8624-9ce7a0fd08a1
# ╠═f4c44654-5802-11ee-194a-0fc1aed8aabc
# ╠═0eff32b2-2594-41f4-98d0-9d9d650204ff
# ╟─c944c2bc-40fe-4405-a8fc-6d81434662d0
# ╠═7f3041e6-25d2-40d4-98ad-236ecaf4ef40
# ╠═63847c84-bd9e-4bf7-be5d-5e90505b9386
# ╟─3ba29c8b-8c1c-4967-b3d7-9568bccd232c
# ╠═c8321163-99e6-4e14-8f93-45f308370940
# ╟─1c744290-a859-4a2d-8e4e-cd951331a476
# ╠═1e7e0136-a3de-4f49-8015-bb16166aeb8f
# ╟─047c1b14-28e9-4237-b50f-ce2d9c3fd491
# ╟─68dbe732-236d-4af2-b98d-0a2f8d2afcfa
# ╠═ee86760b-9308-4a25-b8b8-60adb442b7a9
# ╟─e3eeecd4-0c17-45c8-b83e-cd071a02bad5
# ╟─2fcfdee4-6c7a-4c2d-957c-42759a63e153
# ╟─f38934e1-d100-4ce4-9488-f9cda003d24f
# ╠═0b808b21-23da-4f45-8fee-3af0f549cf70
# ╟─18d23ed4-5eff-45e7-be12-e83024dd9bf0
# ╠═17b3b2ae-5040-4eb3-9978-3f97bdeb7701
# ╟─5c4f9766-9aba-4907-bde9-2af4c25b4d39
# ╠═211c9bdc-31bb-45cd-b120-ff9a754afccc
# ╠═308435df-b47b-416b-b141-81050a1a9f6b
# ╠═59726ecd-f3e2-4a94-8a70-8765f3d7e63a
# ╟─714de8ab-1422-4f6e-9b2d-a92262ac91d7
# ╠═5765f3c8-7b6d-465a-bd5e-58bb314edfeb
# ╠═abcf719d-04d5-4ac5-aedf-10555d2687be
# ╠═b8451879-de3f-4b3a-8402-358543c2597c
# ╟─7c2d7a30-9937-4594-af63-5994935a696d
# ╠═348bd3e6-c6f1-46b9-b162-58557ab9e224
# ╟─1ba194c4-d4c1-4873-857d-12b1c07c81bc
# ╠═730897a2-cbe7-480e-a6ae-65261541aa06
# ╠═62803a99-78e3-4663-986a-77c04cc99eba
# ╠═6fa23358-8fe7-47bf-90c3-2fa4321a25a9
# ╠═c81f01ce-f8f5-4309-b977-0dc734eec48e
# ╠═8f5fcba4-fcb4-4507-9560-7f7f551978c8
# ╟─9a8c9810-3bf4-4a14-a53c-d4bc0600bfd8
# ╟─563d4760-e924-41c3-abfb-08e984815600
# ╠═2b93634c-315f-42c5-ba3c-9c4083b604ca
# ╠═be108252-cdc2-4a6d-ad44-3c07051842cb
# ╟─970e10c5-5894-490b-8972-b69c3365c59f
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
