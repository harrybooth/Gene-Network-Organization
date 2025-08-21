# Plotting Recipies

using CairoMakie

custom_neglog_formatter(values) = map(
    v -> "-10" * Makie.UnicodeFun.to_superscript(round(Int64, v)),
    values
	)

custom_poslog_formatter(values) = map(
    v -> "10" * Makie.UnicodeFun.to_superscript(round(Int64, v)),
    values
    )

custom_log_formatter(values) = map(
    v -> "10" * Makie.UnicodeFun.to_superscript(round(Int64, v)),
    values
    )

function compute_x(x, width, gap, dodge, dodge_gap)
    scale_width(dodge_gap, n_dodge) = (1 - (n_dodge - 1) * dodge_gap) / n_dodge
    function shift_dodge(i, dodge_width, dodge_gap)
        (dodge_width - 1) / 2 + (i - 1) * (dodge_width + dodge_gap)
    end
    width *= 1 - gap
    n_dodge = maximum(dodge)
    dodge_width = scale_width(dodge_gap, n_dodge)
    shifts = shift_dodge.(dodge, dodge_width, dodge_gap)
    return x .+ width .* shifts
end

binomialp(k,p,n) = (factorial(n)/(factorial(k)*factorial(n-k)))*(p^k)*((1-p)^(n-k))

function prob_k_mutations(k,p,n)

    binomialp(k,p,n) / (1 - binomialp(0,p,n))

end


mutable struct dynamical_summary_config

    fontsize
    embed_markersize
    color_scheme
    node_colors
    draw_config
    fitness_linewidth
    fitness_markersize
    pheno_linewidth
    color_fade
    caption_padding
    caption_fontsize

end

mutable struct evo_summary_config

    fontsize
    wait_markersize
    wait_linewidth
    color_scheme
    node_colors
    draw_config
    color_fade

    pie_radius
    pie_inner_radius
    pie_colors
    pie_strokewidth
end

mutable struct pred_config
    fontsize
    perf_linewidth
    entropy_markersize
    entropy_linewidth

    entropy_hist_scale 
end

###############

seaborn_palette = palette(:seaborn_colorblind)
ext_palette = palette(:Set1_8)
epi_colors = [reverse(palette(:tab10)[1:4])[1],reverse(palette(:tab10)[1:4])[3],reverse(palette(:tab10)[1:4])[2],reverse(palette(:tab10)[1:4])[4]]

const top_n_colors = seaborn_palette[[1,2,3,5,end]]

const node_colors = seaborn_palette[[8,7,4,9]];

#################################################################

weight_indices = Tuple.(findall(viable_mutations.> 0));

vertex_names = Dict(1=>"A",2=> "B", 3=> "C", 4=> "M")

weight_names = [string(vertex_names[last(t)]) * "=>" * string(vertex_names[first(t)]) for t in weight_indices];

const weight_names_latex = [L"W_{aa}",L"W_{ab}",L"W_{ac}",L"W_{ba}",L"W_{bb}",L"W_{bc}",L"W_{ca}",L"W_{cb}",L"W_{cc}",L"W_{ma}"];

const weight_names_latex_m = reshape(vcat(weight_names_latex,[L"$W_{mb}$", L"$W_{mc}$"]),(3,4));

draw_config_18 = fs18_default();
draw_config_12 = fs12_default();

fontsize_pub = 12.;

################

function plot_figure_1!(fig,trajectories,embedding,top_n,sorted_uep,sorted_counts_uep,end_parents,vertex_top_map,example_mst,tr_choice,ds_config)

    mo_umap = fig[1:4, 1:4] = GridLayout()

    ex1 = fig[5:7, 1:4] = GridLayout()
    rmh0 = fig[8,1:4] = GridLayout()

    top_n_dict = Dict(v_id=>pos for (pos,v_id) in enumerate(sorted_uep[1:top_n]))

    #### Motif Distribution

    color_sorted_counts_uep = [i <= top_n ? ds_config.color_scheme[i] : :grey for i in 1:length(sorted_counts_uep)]

    view_sorted_uep_id = [i <= top_n for i in 1:length(sorted_counts_uep)]

    n_norm = sum(sorted_counts_uep)

    sorted_uep_proportions = sorted_counts_uep[view_sorted_uep_id] ./ n_norm 

    view_color_sorted_uep = color_sorted_counts_uep[view_sorted_uep_id]

    ##############

    ax1 = Axis(mo_umap[3:4,1:top_n], xlabel = L"\text{Dyn: UMAP 1}", ylabel = L"\text{Dyn: UMAP 2}")

    CairoMakie.scatter!(ax1,embedding, color = [haskey(top_n_dict,i) ? (ds_config.color_scheme[top_n_dict[i]],0.5) : (:grey,0.5) for i in end_parents],markersize = ds_config.embed_markersize)

    hidedecorations!(ax1, label = false)

    for i in 1:top_n

        ax_geno = Axis(mo_umap[2,i], backgroundcolor = (ds_config.color_scheme[i],ds_config.color_fade),aspect = DataAspect())

        draw_grn!(ax_geno,vertex_top_map[sorted_uep[i]],ds_config.draw_config,ds_config.node_colors,ds_config.fontsize,false,false)
    end

    #######################

    ax_mo = Axis(mo_umap[1,1:top_n],ylabel  = L"\text{Probabilty}", xlabel = L"M^{(i)}_{N_i}")

    CairoMakie.barplot!(ax_mo,sorted_uep_proportions,color = view_color_sorted_uep)

    ax_mo.xticks = (1:length(sorted_uep_proportions),string.(1:length(sorted_uep_proportions[1:end])))

    ax_mo.yticks = ([0.,0.1,0.2,0.3],[L"0",L"0.1",L"0.2",L"0.3"])

    CairoMakie.ylims!(ax_mo,0.,0.3)

    CairoMakie.hidedecorations!(ax_mo,label = false,ticklabels = false,ticks = false,minorticks = false)

    ###################

    tr_data = filter(tr->tr.inc_metagraph_vertices[tr.H0] ∈ sorted_uep[example_mst],trajectories);
    tr_data_id = findall(tr->tr.inc_metagraph_vertices[tr.H0] ∈ sorted_uep[example_mst],trajectories)

    tr_traj_id = uniqueid(tr_data[tr_choice].topologies)
    tr_traj_id[end] = length(tr_data[tr_choice].topologies)

    tr_traj = tr_data[tr_choice].topologies[tr_traj_id]
    tr_networks = tr_data[tr_choice].geno_traj[tr_traj_id]

    development = DefaultGRNSolver()

    tr_phenotypes = [Individual(reshape(net,(3,4)),grn_parameters,development).phenotype.u[end] for net in tr_networks]

    tr_fd = create_full_fitness_traj(tr_data[tr_choice].fitness_traj_tuple,tr_data[tr_choice].wait_times)

    tr_top_fitness_id = uniqueid(tr_fd)[tr_traj_id]
    
    tr_fd_coarse = map(x->x[1]+1,tr_fd)
    tr_fd_refine = map(x->x[2],tr_fd)

    tr_top_fitness_rf = [(x,tr_fd_refine[x]) for x in tr_top_fitness_id]

    tr_top_stripe_id = [tr_fd_coarse[x] for x in tr_top_fitness_id]

    ###################

    progression_cs = palette(:haline,length(tr_top_fitness_rf))

    ax_fitness = Axis(ex1[1:2,1:length(tr_phenotypes)],xlabel = L"\text{Generation}",ylabel = L"\text{Fitness}")

    hideydecorations!(ax_fitness,label = false,ticklabels = false,ticks = false,minorticks = false)

    rline = CairoMakie.lines!(ax_fitness,tr_fd_refine, color = :grey, linewidth = ds_config.fitness_linewidth)
    cline = CairoMakie.lines!(ax_fitness,tr_fd_coarse, linestyle = "--", color = :blue,linewidth = ds_config.fitness_linewidth)

    CairoMakie.scatter!(ax_fitness,tr_top_fitness_rf, color = [progression_cs[i] for i in 1:length(tr_top_fitness_rf)], markersize = ds_config.fitness_markersize, marker = '★')

    h0 = tr_top_fitness_id[minimum(findall(tr_top_stripe_id .== 1))]
    Ni = tr_top_fitness_id[end]

    if h0 != Ni
        v = Int.(floor((h0+Ni)/2))
        ax_fitness.xticks = ([1,h0,v,Ni],[L"1",L"S_0",L"%$v",L"N_i"])
    else
        ax_fitness.xticks = ([1,h0],[L"1",L"S_0 = N_i"])
    end


    Legend(ex1[1:2,1:length(tr_phenotypes)],  [rline, cline], [L"\mathcal{F}_R(\phi)", L"\mathcal{F}_S(\phi)"], framevisible=false,orientation = :vertical,patchsize = (10, 10),rowgap = 2,halign = :right, valign = :bottom)

    ####################

    ax_pheno_list = []

    for i in 1:length(tr_phenotypes)

        if tr_top_stripe_id[i] == 1
            ax_geno = Axis(ex1[3:4,i], backgroundcolor = (ds_config.color_scheme[example_mst],ds_config.color_fade),aspect = DataAspect())
        else
            ax_geno = Axis(ex1[3:4,i], backgroundcolor = RGBf(0.98, 0.98, 0.98),aspect = DataAspect())
        end

        ax_pheno = Axis(ex1[5,i],alignmode=Mixed(bottom=0))

        for g in 1:3
            CairoMakie.lines!(ax_pheno,tr_phenotypes[i][g,:],linewidth = ds_config.pheno_linewidth, color = ds_config.node_colors[g])
        end

        CairoMakie.scatter!(ax_pheno,[(90,0.8*tr_phenotypes[end][3,50])],color = progression_cs[i], markersize = ds_config.fitness_markersize,marker = '★')

        CairoMakie.hidedecorations!(ax_pheno)

        draw_grn!(ax_geno,tr_traj[i],ds_config.draw_config,ds_config.node_colors,ds_config.fontsize,false,false)

        push!(ax_pheno_list,ax_pheno)
    end

    linkyaxes!(ax_pheno_list...)

    ##########################

    all_prop  = []
    all_dodge  = []
    all_x  = []

    for n in 1:top_n+1

        if n == top_n+1
            pop = filter(tr->!(tr.inc_metagraph_vertices[end] ∈  sorted_uep[1:top_n]),trajectories)
        else
            pop = filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)
        end

        pop_equal = filter(tr->tr.minimal_stripe_subgraphs[tr.H0] == tr.minimal_stripe_subgraphs[end], pop)

        pop_H0_incl_N = filter(tr->Bool(test_inclusion(tr.minimal_stripe_subgraphs[end],tr.minimal_stripe_subgraphs[tr.H0])) & !(tr.minimal_stripe_subgraphs[end] == tr.minimal_stripe_subgraphs[tr.H0]),pop)

        pop_N_incl_H0 = filter(tr->Bool(test_inclusion(tr.minimal_stripe_subgraphs[tr.H0],tr.minimal_stripe_subgraphs[end])) & !(tr.minimal_stripe_subgraphs[end] == tr.minimal_stripe_subgraphs[tr.H0]),pop)

        n_pop = length(pop)

        proportions = [length(pop_equal),length(pop_H0_incl_N),length(pop_N_incl_H0),length(pop) - length(pop_equal) - length(pop_N_incl_H0) - length(pop_H0_incl_N)]

        @assert sum(proportions) == n_pop

        x = [1,2,3,4]

        dodge = [n,n,n,n]

        push!(all_prop,proportions ./ n_pop)
        push!(all_dodge,dodge)
        push!(all_x,x)

    end

    ax_rh0 = Axis(rmh0[1,1],alignmode=Mixed(top=0), ylabel = L"\text{% of trajectories}")

    x = reduce(vcat,all_x)
    dodge = reduce(vcat,all_dodge)
    proportions = reduce(vcat,all_prop)

    CairoMakie.barplot!(ax_rh0,x,proportions,color = [n==top_n+1 ? :grey : ds_config.color_scheme[n] for n in dodge],dodge = dodge)

    CairoMakie.hidedecorations!(ax_rh0,label = false,ticklabels = false,ticks = false,minorticks = false)

    ax_rh0.xticks = (1:4,[L"M^{(i)}_{S_{0}} = M^{(i)}_{N_i}",L"M^{(i)}_{S_{0}} \subset M^{(i)}_{N_i}",L"M^{(i)}_{N_i} \subset M^{(i)}_{S_{0}}",L"\text{MST change}"])

    CairoMakie.ylims!(ax_rh0,0.,1.)

    ax_rh0.yticks = ([0.,0.5,1.],[L"\text{0}",L"\text{50}",L"\text{100}"])

    for (label, layout) in zip(["A", "C", "E"], [mo_umap, ex1, rmh0])
        Label(layout[1, 1, TopLeft()], label,
            fontsize = ds_config.caption_fontsize,
            font = :bold,
            padding = (0,ds_config.caption_padding, ds_config.caption_padding, 0),
            halign = :right)
    end

    Label(mo_umap[3, 1, TopLeft()], "B",
    fontsize = ds_config.caption_fontsize,
    font = :bold,
    padding = (0,ds_config.caption_padding, ds_config.caption_padding, 0),
    halign = :right)

    Label(ex1[3, 1, TopLeft()], "D",
    fontsize = ds_config.caption_fontsize,
    font = :bold,
    padding = (0,ds_config.caption_padding, ds_config.caption_padding, 0),
    halign = :right)
    
    colgap = 5
    rowgap = 10

    colgap!(mo_umap,colgap)
    rowgap!(mo_umap, rowgap)

    colgap!(ex1, colgap)
    rowgap!(ex1, rowgap)

    colgap!(rmh0,colgap)
    rowgap!(rmh0, rowgap)

    rowgap!(fig.layout, Relative(0.01))
    colgap!(fig.layout, Relative(0.01))

end

function create_pheno_transition_bars!(fig,trajectories,all_pheno_mut_chf,evo_config)

    pheno_prop_subplot = fig[1,1] = GridLayout()

    ax_prop = Axis(pheno_prop_subplot[1,1],yticklabelsize = 0.8*evo_config.fontsize,yaxisposition = :right, yticksize= 0.25*evo_config.fontsize,xticklabelrotation = pi/2,xlabelpadding = 0.,xgridvisible = false,ygridvisible = false)

    ax_prop.xticks = (1:4,[L"n=1",L"1 < n < S_0",L"n=S_{0}", L"n>S_{0}"])
    ax_prop.yticks = ([0.,0.5,1.],[L"0",L"50", L"100"])

    ############################### Pheno

    Random.seed!(1234567)

    arrow_colors = shuffle([i for i in palette(:Dark2_8)])

    arrow_colors = vcat([:black,palette(:seaborn_dark)[end]],arrow_colors)

    ph_id = findall(tr->(tr.H0-2 > 0),trajectories)

    mtp_1  = [mc[1] for (mc,tr) in zip(all_pheno_mut_chf[ph_id],trajectories[ph_id])]
    mtp_2  = [mc[2:tr.H0-2] for (mc,tr) in zip(all_pheno_mut_chf[ph_id],trajectories[ph_id])]
    mtp_3  = [mc[tr.H0-1] for (mc,tr) in zip(all_pheno_mut_chf[ph_id],trajectories[ph_id])]
    mtp_4  = [mc[tr.H0:end] for (mc,tr) in zip(all_pheno_mut_chf[ph_id],trajectories[ph_id])];

    mtp_1 = map(x->x[1],mtp_1)
    mtp_2 = map(x->x[1],reduce(vcat,mtp_2))
    mtp_3 = map(x->x[1],mtp_3)
    mtp_4 = map(x->x[1],reduce(vcat,mtp_4));

    pheno_edge_order = [[:neutral],[:other_wnb],[:clb],[:crb],[:cbb],[:mlb_wlb,:other_wlb],[:crb_wlb],[:mrb_wrb,:other_wrb],[:clb_wrb],[:mbb_wbb,:other_wbb]]

    mtp_1_prop = calculate_mut_type_proportion_list(mtp_1,pheno_edge_order)
    mtp_2_prop = calculate_mut_type_proportion_list(mtp_2,pheno_edge_order)
    mtp_3_prop = calculate_mut_type_proportion_list(mtp_3,pheno_edge_order)
    mtp_4_prop = calculate_mut_type_proportion_list(mtp_4,pheno_edge_order);

    m_prop_all = []
    m_pheno_labels = []
    m_time_labels = []

    for (n,m_prop) in enumerate([mtp_1_prop,mtp_2_prop,mtp_3_prop,mtp_4_prop])
        push!(m_prop_all,m_prop)
        push!(m_pheno_labels, [i for i in 1:length(pheno_edge_order)])
        push!(m_time_labels,[n for _ in 1:length(m_prop)])
    end

    m_prop_all = reduce(vcat,m_prop_all)
    m_pheno_labels = reduce(vcat,m_pheno_labels)
    m_time_labels = reduce(vcat,m_time_labels);

    CairoMakie.barplot!(ax_prop,m_time_labels,m_prop_all,stack = m_pheno_labels,color = [arrow_colors[i] for i in m_pheno_labels])

end

function create_wait_times!(fig,trajectories,evo_config,wait_time_summary)

    ax_wait_2 = Axis(fig[1,1], yticklabelcolor = :black,yscale = log10,yticklabelsize = 0.8*evo_config.fontsize,yticksize= 0.25*evo_config.fontsize,xgridvisible = false,ygridvisible = false)

    all_wait_times_1 = reduce(vcat,[tr.wait_times[2] for tr in trajectories])
    all_wait_times_ls0 = reduce(vcat,[tr.wait_times[3:tr.H0-1] for tr in trajectories])
    all_wait_times_s0 = reduce(vcat,[tr.wait_times[tr.H0] for tr in trajectories])
    all_wait_times_hs0 = reduce(vcat,[tr.wait_times[tr.H0+1:end] for tr in trajectories])

    wait_color = :black

    if wait_time_summary == :mean

        mean_wait = [mean(x) for x in [all_wait_times_1,all_wait_times_ls0,all_wait_times_s0,all_wait_times_hs0]]

        std_error_wait = [std(x) ./ sqrt(length(x)) for x in [all_wait_times_1,all_wait_times_ls0,all_wait_times_s0,all_wait_times_hs0]]

        mean_wait_type_labels = [1.,2.,3.,4.]

        wt_l = CairoMakie.lines!(ax_wait_2,mean_wait_type_labels,mean_wait,color = wait_color,linewidth = evo_config.wait_linewidth)
        wt_s = CairoMakie.scatter!(ax_wait_2,mean_wait_type_labels,mean_wait,color = wait_color,markersize = evo_config.wait_markersize)

        CairoMakie.errorbars!(ax_wait_2,mean_wait_type_labels,mean_wait,5 * std_error_wait,color = wait_color,whiskerwidth = evo_config.wait_markersize/2)

    else

        median_wait_time = [median(x) for x in [all_wait_times_1,all_wait_times_ls0,all_wait_times_s0,all_wait_times_hs0]]
        lq_wait_time = [quantile(x, [0.25])[1] for x in [all_wait_times_1,all_wait_times_ls0,all_wait_times_s0,all_wait_times_hs0]]
        uq_wait_time = [quantile(x, [0.75])[1] for x in [all_wait_times_1,all_wait_times_ls0,all_wait_times_s0,all_wait_times_hs0]]

        median_wait_type_labels = [1.,2.,3.,4.]

        wt_l = CairoMakie.lines!(ax_wait_2,median_wait_type_labels,median_wait_time,color = wait_color,linewidth = evo_config.wait_linewidth)
        wt_s = CairoMakie.scatter!(ax_wait_2,median_wait_type_labels,median_wait_time,color = wait_color,markersize = evo_config.wait_markersize)

        CairoMakie.rangebars!(ax_wait_2,median_wait_type_labels,lq_wait_time,uq_wait_time,color = wait_color,whiskerwidth = evo_config.wait_markersize/2)
    end

    ax_wait_2.xticks = (1:4,[L"n=1",L"1<n<S_0",L"n=S_0",L"n>S_0"])
    ax_wait_2.ylabel = L"\text{Mutant wait time}"

end

function create_nweight!(fig,mut_prob,all_count_prop,epi_mutn_prob_1,epi_mutn_prob_lS0,epi_mutn_prob_S0,epi_mutn_prob_hS0)

    null_mut_n = [prob_k_mutations(k,mut_prob,10) for k in 1:4];
    edges_x = 0.5:1:5
    edges_y = 1:4

    titles = [L"n=1",L"1<n<S_0",L"n=S_0",L"n>S_0"]

    all_ax_u = []
    all_hm = []

    for (n,mut_prob_v) in enumerate([epi_mutn_prob_1,epi_mutn_prob_lS0,epi_mutn_prob_S0,epi_mutn_prob_hS0])

        count_prop = all_count_prop[n][1]

        d = count_prop[1:4,1:4] 

        for i in 1:4
            for j in 1:4
                if i<j 
                    d[i,j] = NaN
                end
            end
        end

        fig_u = fig[1,n] = GridLayout()

        ax_u = Axis(fig_u[1:2,1:2],title = titles[n], xgridvisible = false, ygridvisible = false)
        ax_h = Axis(fig_u[3:4,1:2],yreversed = true,ylabel = L"k")

        hidespines!(ax_h)

        if n == 1
            hidexdecorations!(ax_h)
            hideydecorations!(ax_h,ticklabels = false,label = false)
            ax_h.yticks = (1:4,[L"1",L"2",L"3",L"4"])
        else
            hidedecorations!(ax_h)
            hideydecorations!(ax_u)
        end

        hm = CairoMakie.heatmap!(ax_h,edges_y,edges_x,d,colorrange = (0.,1.),colormap = :thermal)

        CairoMakie.lines!(ax_u,null_mut_n, color = :black,linestyle = :dash)
        CairoMakie.scatter!(ax_u,null_mut_n, color = :black,marker = 'x')

        CairoMakie.barplot!(ax_u,mut_prob_v,color = (:grey,0.5))

        ax_u.xticks = (1:4,[L"1",L"2",L"3",L"4"])

        CairoMakie.ylims!(ax_u,0.,0.7)

        ax_u.yticks = ([0.,0.35,0.7], [L"0",L"35",L"70"])

        linkxaxes!(ax_u,ax_h)

        push!(all_ax_u,ax_u)
        push!(all_hm,hm)

        rowgap!(fig_u,Relative(0.01))
        colgap!(fig_u,Relative(0.02))
    end

    linkyaxes!(all_ax_u...)

    rowgap!(fig.layout,Relative(0.01))
    colgap!(fig.layout,Relative(0.025))
end

function create_epi_type_counts_neutral_square!(fig,trajectories,check_path,ignore_neutral,top_n_colors,fitness_eps)

    epi_counts_1 = map(tr->tr.epistasis[1],filter(tr->(tr.H0 > 2),trajectories))

    epi_counts_lS0 = map(tr->tr.epistasis[2:tr.H0-2],filter(tr->(tr.H0-2 > 0),trajectories))

    epi_counts_S0 = map(tr->tr.epistasis[tr.H0-1],trajectories)

    epi_counts_hS0 = map(tr->tr.epistasis[tr.H0:end],trajectories);

    if check_path
        n_ind_1 = map(tr->all(check_reg_path(tr.mutant_info[1].weight_id,tr.mutant_info[1].start_network)),filter(tr->(tr.H0 > 2),trajectories));
        n_ind_lS0 = reduce(vcat,map(tr->[all(check_reg_path(trm.weight_id,trm.start_network)) for trm in tr.mutant_info[2:tr.H0-2]],filter(tr->(tr.H0-2 > 0),trajectories)))
        n_ind_S0 = map(tr->all(check_reg_path(tr.mutant_info[tr.H0-1].weight_id,tr.mutant_info[tr.H0-1].start_network)),trajectories);
        n_ind_hS0 = reduce(vcat,map(tr->[all(check_reg_path(trm.weight_id,trm.start_network)) for trm in tr.mutant_info[tr.H0:end]],trajectories));

        epi_1 = epi_counts_1[map(x->(length(x[1]) == 4) & x[2] ,zip(epi_counts_1,n_ind_1))]
        epi_lS0 = reduce(vcat, epi_counts_lS0)[map(x->(length(x[1]) == 4) & x[2] ,zip(reduce(vcat, epi_counts_lS0),n_ind_lS0))]
        epi_S0 = epi_counts_S0[map(x->(length(x[1]) == 4) & x[2] ,zip(epi_counts_S0,n_ind_S0))]
        epi_hS0 = reduce(vcat, epi_counts_hS0)[map(x->(length(x[1]) == 4) & x[2] ,zip(reduce(vcat, epi_counts_hS0),n_ind_hS0))];
    else
        epi_1 = epi_counts_1[map(x->length(x) == 4, epi_counts_1)]
        epi_lS0 = reduce(vcat, epi_counts_lS0)[map(x->length(x) == 4,reduce(vcat, epi_counts_lS0))]
        epi_S0 = epi_counts_S0[map(x->length(x) == 4 ,epi_counts_S0)]
        epi_hS0 = reduce(vcat, epi_counts_hS0)[map(x->length(x) == 4,reduce(vcat, epi_counts_hS0))];
    end

    if ignore_neutral
        epi_1ch = filter(x->x!=:Neutral,map(x->type_epi(x[1],x[2],x[3],x[4],fitness_eps),epi_1))
        epi_lS0ch = filter(x->x!=:Neutral,map(x->type_epi(x[1],x[2],x[3],x[4],fitness_eps),epi_lS0))
        epi_S0ch = filter(x->x!=:Neutral,map(x->type_epi(x[1],x[2],x[3],x[4],fitness_eps),epi_S0))
        epi_hS0ch = filter(x->x!=:Neutral,map(x->type_epi(x[1],x[2],x[3],x[4],fitness_eps),epi_hS0));

        epi_types = [:ne,:me,:se,:rse]

        print(length(epi_1ch))
        print(":")
        # c = countmap(epi_1ch)
        # epi_1_prop = [haskey(c,p) ? c[p]/length(epi_1ch) : 0. for p in epi_types]

        print(length(epi_lS0ch))
        print(":")
        c = countmap(epi_lS0ch)
        epi_lS0_prop = [haskey(c,p) ? c[p]/length(epi_lS0ch) : 0. for p in epi_types]

        print(length(epi_S0ch))
        print(":")
        c = countmap(epi_S0ch)
        epi_S0_prop = [haskey(c,p) ? c[p]/length(epi_S0ch) : 0. for p in epi_types]

        print(length(epi_hS0ch))
        print(":")
        c = countmap(epi_hS0ch)
        epi_hS0_prop = [haskey(c,p) ? c[p]/length(epi_hS0ch) : 0. for p in epi_types]

        for (i,l) in enumerate([epi_lS0_prop,epi_S0_prop,epi_hS0_prop])
            ax = Axis(fig[1,i],aspect = DataAspect())
            CairoMakie.pie!(ax,l,color = [top_n_colors[3],top_n_colors[end],top_n_colors[2],top_n_colors[4]])
        
            hidexdecorations!(ax)
            hideydecorations!(ax)
            hidespines!(ax)
        end
    else
        epi_1ch = map(x->type_epi(x[1],x[2],x[3],x[4],fitness_eps),epi_1)
        epi_lS0ch = map(x->type_epi(x[1],x[2],x[3],x[4],fitness_eps),epi_lS0)
        epi_S0ch = map(x->type_epi(x[1],x[2],x[3],x[4],fitness_eps),epi_S0)
        epi_hS0ch = map(x->type_epi(x[1],x[2],x[3],x[4],fitness_eps),epi_hS0);
    
        epi_types = [:Neutral,:ne,:me,:se,:rse]
    
        c = countmap(epi_1ch)
        epi_1_prop = [haskey(c,p) ? c[p]/length(epi_1ch) : 0. for p in epi_types]
    
        c = countmap(epi_lS0ch)
        epi_lS0_prop = [haskey(c,p) ? c[p]/length(epi_lS0ch) : 0. for p in epi_types]
    
        c = countmap(epi_S0ch)
        epi_S0_prop = [haskey(c,p) ? c[p]/length(epi_S0ch) : 0. for p in epi_types]
    
        c = countmap(epi_hS0ch)
        epi_hS0_prop = [haskey(c,p) ? c[p]/length(epi_hS0ch) : 0. for p in epi_types]

        grid_positions = [(1,1),(1,2),(2,1),(2,2)]

        for (i,l) in enumerate([epi_1_prop,epi_lS0_prop,epi_S0_prop,epi_hS0_prop])
            ax = Axis(fig[grid_positions[i]...],aspect = DataAspect())
            CairoMakie.pie!(ax,l,color = [top_n_colors[1],top_n_colors[3],top_n_colors[end],top_n_colors[2],top_n_colors[4]])
        
            hidexdecorations!(ax)
            hideydecorations!(ax)
            hidespines!(ax)
        end
    end

end

function create_epi_correlation_summary!(fig,all_2mut_epi_char,all_2mut_epi_char_wcorr,mst_corr,rse_ex,ne_ex,top_n_colors)

    order_ = [([0, 0], :Neutral),([0, 0], :ne),([0, 0], :me),([0, 0], :se),([0, 0], :rse),([0, 1], :Neutral),([0, 1], :ne),([0, 1], :me),([0, 1], :se),([0, 1], :rse),([1, 1], :Neutral),([1, 1], :ne),([1, 1], :me),([1, 1], :se),([1, 1], :rse)]

    n_dict = Dict(:Neutral=>1,:ne=>2,:me=>3,:se=>4,:rse=>5)
    n_dictr = Dict(v=>k for (k,v) in n_dict)

    epi_colors = Dict(:Neutral=>top_n_colors[1],:ne=>top_n_colors[3],:me=>top_n_colors[end],:se=>top_n_colors[2],:rse=>top_n_colors[4])

    fig_top = fig[1,1] = GridLayout()

    n = 5

    ax = Axis(fig_top[1,2],title = string(round(corspearman(rse_ex[n][1],rse_ex[n][2]),digits = 2)) * " _ " * string(round(corspearman(ne_ex[n][1],ne_ex[n][2]),digits = 2)),xgridvisible = false,ygridvisible = false, xreversed = true, yreversed = true)
    
    CairoMakie.scatter!(ax,sign.(rse_ex[n][1]) .* log10.(abs.(rse_ex[n][1])),sign.(rse_ex[n][2]) .* log.(abs.(rse_ex[n][2])),color = (top_n_colors[4],0.75),markersize = 7.)
    
    CairoMakie.scatter!(ax,sign.(ne_ex[n][1]) .* log10.(abs.(ne_ex[n][1])),sign.(ne_ex[n][2]) .* log.(abs.(ne_ex[n][2])),color = (top_n_colors[3],0.75),markersize = 7.)
    
    ax.yticks = ([4,2,0,-2,-4],[L"-10^-1",L"-10^-2",L"0",L"10^-2",L"10^-1"])
    
    ax.xticks = ([2,1,0,-1,-2],[L"-10^-1",L"-10^-2",L"0",L"10^-2",L"10^-1"])

    ##################

    all_labels = reduce(vcat,[[i for _ in all_2mut_epi_char_wcorr[i]] for i in 1:5])
    all_data = abs.(Float64.(reduce(vcat,all_2mut_epi_char_wcorr)))
    
    ax = Axis(fig_top[1,3],xgridvisible = false,ygridvisible = false,ylabel = L"\text{Spearman Correlation.}",xlabel = L"\text{Epi. type}")
    
    CairoMakie.boxplot!(ax,all_labels, all_data,
            color = [epi_colors[n_dictr[i]] for i in all_labels],show_outliers = false)
    

    ax.yticks = ([0,0.4,0.8],[L"0",L"0.4",L"0.8"])

    CairoMakie.ylims!(ax,0.,0.95)

    ##########

    fig_bot = fig[2,1] = GridLayout()

    value_x = [1,1,1,1,1,2,2,2,2,2,3,3,3,3,3]
    value_dodge = [1,2,3,4,5,1,2,3,4,5,1,2,3,4,5]

    value_colors = [epi_colors[last(n)] for n in order_]

    d = filter(x->x!= :p,all_2mut_epi_char[2])

    all_2mut_counts = countmap(d)

    all_2mut_types = countmap(first.(d))

    all_2mut_cond_prob = Dict(k=>all_2mut_counts[k]/all_2mut_types[k[1]] for k in keys(all_2mut_counts))

    prob_wr = [all_2mut_types[[0, 0]],all_2mut_types[[0, 1]],all_2mut_types[[1, 1]]] ./ length(d)

    values = [haskey(all_2mut_cond_prob,i) ? all_2mut_cond_prob[i] : 0. for i in order_]

    ##############

    ax = Axis(fig_bot[1,2],xgridvisible = false,ygridvisible = false,ylabel = L"\text{Frequency}",xlabel = L"\text{weight relation}")

    barplot!(ax,value_x, values,
            dodge = value_dodge,
            color = value_colors)

    CairoMakie.lines!(ax, [1,2,3],prob_wr,color = :black)
    CairoMakie.scatter!(ax, [1,2,3],prob_wr,color = :black,markersize = 6.)

    ax.yticks = ([0,0.25,0.5,0.75],[L"0",L"0.25",L"0.5",L"0.75"])

    CairoMakie.ylims!(ax,0.,0.75)

    #####################

    ax = Axis(fig_bot[1,1],xgridvisible = false,ygridvisible = false,ylabel = L"\text{Spearman correlation}")

    ls0_non_mst_corr = mst_corr[1]
    ls0_inter_mst_corr  = mst_corr[2]
    ls0_mst_corr = mst_corr[3]

    hs0_non_mst_corr = mst_corr[4]
    hs0_inter_mst_corr = mst_corr[5]
    hs0_mst_corr = mst_corr[6]

    side_ls0 = reduce(vcat,[[1 for _ in ls0_non_mst_corr],[2 for _ in ls0_inter_mst_corr],[3 for _ in ls0_mst_corr]])
    side_hs0 = reduce(vcat,[[1 for _ in hs0_non_mst_corr],[2 for _ in hs0_inter_mst_corr],[3 for _ in hs0_mst_corr]])

    ls0_all = abs.(reduce(vcat,[ls0_non_mst_corr,ls0_inter_mst_corr, ls0_mst_corr]))
    hs0_all = abs.(reduce(vcat,[hs0_non_mst_corr,hs0_inter_mst_corr, hs0_mst_corr]))

    side_all = reduce(vcat,[side_ls0,side_hs0])
    values_all = reduce(vcat,[ls0_all,hs0_all])
    cat = reduce(vcat,[[1 for _ in ls0_all],[2 for _ in hs0_all]])

    CairoMakie.boxplot!(ax,cat, values_all, dodge= side_all, color = side_all,show_outliers = :false)

    ax.yticks = ([0,0.4,0.8],[L"0",L"0.4",L"0.8"])
    CairoMakie.ylims!(ax,0.,0.95)

    ax.xticks = (1:2,[L"n < S_0",L"n \geq S_0"])

    colgap!(fig.layout,Relative(0.05))
    
end

function create_contingency_accuracy_summary!(fig,vim_trajectories,all_target_prob,null_H0_dist,top_n,predict_label_to_vertex,predict_colors)

    ax_list = []


    for n in 1:top_n+1

        if n == 1
            ax = Axis(fig[1,n],ylabel  = L"\text{Freq.}", xlabel = L"M^{(i)}_{N_i}")
            hidexdecorations!(ax)
        else
            ax = Axis(fig[1,n],xlabel = L"M^{(i)}_{N_i}")
            hideydecorations!(ax)
            hidexdecorations!(ax)
        end

        pop = vim_trajectories[n]

        incl_count_all = []
    
        for k in 1:top_n+1
            incl_count = sum(map(tr->tr.inc_metagraph_vertices == predict_label_to_vertex[k],pop))
            push!(incl_count_all,incl_count)
        end

        p_dist = incl_count_all ./ sum(incl_count_all)

        p_dist_x = [i for i in 1:top_n+1]

        dodge_current = [1 for _ in 1:top_n+1]

        dodge_H0 = [2 for _ in 1:top_n+1]

        x = vcat(p_dist_x,p_dist_x)

        color = vcat([(predict_colors[i],1.) for i in p_dist_x],[:white for i in p_dist_x])

        stroke_color = [(predict_colors[i],1.) for i in x]

        dodge = vcat(dodge_current,dodge_H0)
        y = vcat(p_dist,null_H0_dist)

        vbar = CairoMakie.barplot!(ax,x,y,dodge = dodge,color = color, strokecolor = stroke_color, strokewidth = 1., gap = 0.2)

        vpred = CairoMakie.scatter!(ax,[1,2,3,4,5], all_target_prob[n], marker = 'x', color = predict_colors[n], markersize = 10.)

        CairoMakie.hidedecorations!(ax,label = false,ticklabels = false,ticks = false,minorticks = false)

        axis_label = vcat(collect(string.(1:top_n)),["o"])
        
        ax.xticks = (1:top_n+1,axis_label)

        push!(ax_list,ax)
    end

    colgap!(fig.layout, Relative(0.02))
    rowgap!(fig.layout, Relative(0.075))

    linkyaxes!(ax_list...)
end
