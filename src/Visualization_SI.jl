function plot_closest_hamming_match!(fig,trajectories,top_n,sorted_uep,vertex_top_map,embedding,ds_config)

    is_other = findall(map(tr->!(tr.inc_metagraph_vertices[end] ∈ sorted_uep[1:top_n]),trajectories));

    is_other_h0 = findall(map(tr->!(tr.inc_metagraph_vertices[end] ∈ sorted_uep[1:top_n]) & (tr.inc_metagraph_vertices[tr.H0] ∈ sorted_uep[1:top_n]),trajectories));

    other_mst = reduce(hcat,[tr.minimal_stripe_subgraphs[end] for tr in trajectories[is_other]]);

    mst8 = reduce(hcat,[vertex_top_map[i] for i in sorted_uep[1:n]])

    mst8_other = pairwise(Hamming(),mst8,other_mst,dims = 2);

    hamming_match = [argmin(v) for v in eachcol(mst8_other)];

    hamming_distances = countmap([(minimum(v),n) for (v,n) in zip(eachcol(mst8_other),[i ∈ is_other_h0 for i in is_other])]);

    uep_position_dict = Dict(v=>n for (n,v) in enumerate(sorted_uep));

    end_parent_pos = [uep_position_dict[i] > top_n ? 0 : uep_position_dict[i] for i in end_parents]

    end_parent_pos[is_other] .= hamming_match;

    ax1 = Axis(fig[1,1], xlabel = L"\text{Dyn: UMAP 1}", ylabel = L"\text{Dyn: UMAP 2}",title = L"\text{Original}")

    CairoMakie.scatter!(ax1,embedding,color = [i == 0 ? (:grey,0.5) : (ds_config.color_scheme[i],0.5) for i in end_parent_pos],markersize = ds_config.embed_markersize)

    hidedecorations!(ax1, label = false)

    ax2 = Axis(fig[2,1], xlabel = L"\text{Dyn: UMAP 1}", ylabel = L"\text{Dyn: UMAP 2}", title = L"\text{With closest hamming match}")

    CairoMakie.scatter!(ax2,embedding,color = [i == 0 ? (:grey,0.5) : (ds_config.color_scheme[i],0.5) for i in end_parent_pos],markersize = ds_config.embed_markersize)

    hidedecorations!(ax2, label = false)
end

function plot_hamming_distances!(fig,trajectories,top_n,sorted_uep,vertex_top_map,ds_config)

    is_other = findall(map(tr->!(tr.inc_metagraph_vertices[end] ∈ sorted_uep[1:top_n]),trajectories));

    is_other_h0 = findall(map(tr->!(tr.inc_metagraph_vertices[end] ∈ sorted_uep[1:top_n]) & (tr.inc_metagraph_vertices[tr.H0] ∈ sorted_uep[1:top_n]),trajectories));

    other_mst = reduce(hcat,[tr.minimal_stripe_subgraphs[end] for tr in trajectories[is_other]]);

    mst8 = reduce(hcat,[vertex_top_map[i] for i in sorted_uep[1:n]])

    mst8_other = pairwise(Hamming(),mst8,other_mst,dims = 2);

    hamming_match = [argmin(v) for v in eachcol(mst8_other)];

    hamming_distances = countmap([(minimum(v),n) for (v,n) in zip(eachcol(mst8_other),[i ∈ is_other_h0 for i in is_other])]);

    hamming_keys = reduce(vcat,[[(i,0),(i,1)] for i in 1:1:maximum(first.(keys(hamming_distances)))])

    ax = Axis(fig[1,1], xlabel = L"\text{Min. hamming distance to MST 1-8}", ylabel = L"\text{Frequency}")

    barplot!(ax,first.(hamming_keys),[haskey(hamming_distances,i) ? hamming_distances[i][1]/length(is_other) : 0. for i in hamming_keys], stack = last.(hamming_keys), color = last.(hamming_keys))

    ax.xticks = (1:1:maximum(first.(keys(hamming_distances))),string.(1:1:maximum(first.(keys(hamming_distances)))))
end


function identify_loop_mutations!(fig,trajectories,evo_config,weight_loop_dict)

    wait_subplot = fig[1,1] = GridLayout()

    ax_wait = Axis(wait_subplot[1,1],yticklabelsize = 0.8*evo_config.fontsize,yaxisposition = :right, xticklabelsvisible = false, yticksize= 0.25*evo_config.fontsize,xgridvisible = false,ygridvisible = false)

    wait_color = :black

    ph_id = findall(tr->(tr.H0-2 > 0),trajectories)

    CairoMakie.ylims!(ax_wait,0,1.2)

    function get_mut_weight_full(tr,range_l,range_u,weight_loop_dict)

        d = [[weight_loop_dict[x] for x in mi[:weight_id]] for mi in tr.mutant_info[range_l:range_u]]
    
        if length(d) >= 1
            return d
        else
            return []
        end
    
    end

    ###############################

    mut_type_prop_all = []
    mut_type_time_labels = []
    mut_type_labels = []

    mut_type_prop = reduce(vcat,[get_mut_weight_full(tr,1,1,weight_loop_dict) for tr in trajectories[ph_id]])

    mut_type_prop_av = [count(x->m ∈ x,mut_type_prop) for m in [:morph_loop,:ac_loop,:abc_loop]] ./ length(mut_type_prop)

    push!(mut_type_prop_all,mut_type_prop_av)
    push!(mut_type_labels, [1,2,3])
    push!(mut_type_time_labels,[1,1,1])

    mut_type_prop = reduce(vcat,[get_mut_weight_full(tr,2,tr.H0-2,weight_loop_dict) for tr in trajectories[ph_id]])

    mut_type_prop_av = [count(x->m ∈ x,mut_type_prop) for m in [:morph_loop,:ac_loop,:abc_loop]] ./ length(mut_type_prop)

    push!(mut_type_prop_all,mut_type_prop_av)
    push!(mut_type_labels, [1,2,3])
    push!(mut_type_time_labels,[2,2,2])

    mut_type_prop = reduce(vcat,[get_mut_weight_full(tr,tr.H0-1,tr.H0-1,weight_loop_dict) for tr in trajectories[ph_id]])

    mut_type_prop_av = [count(x->m ∈ x,mut_type_prop) for m in [:morph_loop,:ac_loop,:abc_loop]] ./ length(mut_type_prop)

    push!(mut_type_prop_all,mut_type_prop_av)
    push!(mut_type_labels, [1,2,3])
    push!(mut_type_time_labels,[3,3,3])

    mut_type_prop = reduce(vcat,[get_mut_weight_full(tr,tr.H0,length(tr.topologies)-1,weight_loop_dict) for tr in trajectories[ph_id]])

    mut_type_prop_av = [count(x->m ∈ x,mut_type_prop) for m in [:morph_loop,:ac_loop,:abc_loop]] ./ length(mut_type_prop)

    push!(mut_type_prop_all,mut_type_prop_av)
    push!(mut_type_labels, [1,2,3])
    push!(mut_type_time_labels,[4,4,4])

    mut_type_prop_all = reduce(vcat,mut_type_prop_all)

    mut_type_time_labels = reduce(vcat,mut_type_time_labels)
    mut_type_labels = reduce(vcat,mut_type_labels); 

    CairoMakie.barplot!(ax_wait,mut_type_time_labels,mut_type_prop_all,dodge = mut_type_labels,color = mut_type_labels,bar_labels = :y,colormap = [:red, :green, :blue],label_formatter = x-> "$(round(x*100,digits = 1))%",label_size = 10.)

    # CairoMakie.hidexdecorations!(ax_wait,ticklabels = false)
    CairoMakie.hideydecorations!(ax_wait,label = false)

    ax_wait.xticks = (1:4,[L"n=1",L"1 < n < S_0",L"n=S_{0}", L"n>S_{0}"])
    ###############################\

    labels_mut =  [L"\text{m/a pathway}",L"\text{a/c pathway}",L"\text{a/b/c pathway}"]

    symbol_mut = [PolyElement(color=c) for c in [:red, :green, :blue]]

    legend_row_gap = 3

    Legend(wait_subplot[1, :, Bottom()], symbol_mut, labels_mut, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), colgap = 4, rowgap = 4, padding=(0.,0.,0f0, evo_config.fontsize+legend_row_gap))


end

######## convenience functions 

function edit_sn(weights)

    new = zeros(12)
    old = zeros(12)
    
    if length(weights) > 0 
        for w in weights
            if w[2] ==1
                new[w[1]] = 1
            else
                new[w[1]] = 1
                old[w[1]] = 1
            end
        end
    end

    return (new,old)

    # return n1 .!= 0

end

### prediction Analysis

function create_confusion_matrices_class(pred_scores,bucket_n,normalize)

    # columns eventual MST outcome

    predc = map(p->p[7],pred_scores)[bucket_n]
    predp = map(p->p[8],pred_scores)[bucket_n];

    predcl = map(v->map(x->x==-1 ? 5 : vertex_to_predict_label[x],v),predc)

    predclp = [[p ? l : n==l ? 6 : l for (p,l) in zip(persist,class_labels)] for (n,(persist,class_labels)) in enumerate(zip(predp,predcl))]

    confmat = reduce(hcat,[[count(x->x==n,class_pred) for n in 1:6] for class_pred in predclp]) 

    # confmatp = confmat ./ sum(confmat,dims = 1)

    if normalize
        confmatp = confmat ./ sum(confmat,dims = 1)
    else
        confmat
    end
end

function create_confusion_matrices_pred(pred_scores,bucket_n,normalize)

    # columns prediction

    predc = map(p->p[9],pred_scores)[bucket_n]
    predp = map(p->p[10],pred_scores)[bucket_n];

    # predcl = map(v->map(x->x==-1 ? 5 : vertex_to_predict_label[x],v),predc)

    predclp = [[p ? l : n==l ? 6 : l for (p,l) in zip(persist,class_labels)] for (n,(persist,class_labels)) in enumerate(zip(predp,predcl))]

    confmat = reduce(hcat,[[count(x->x==n,class_pred) for n in 1:6] for class_pred in predclp]) 

    if normalize
        confmatp = confmat ./ sum(confmat,dims = 1)
    else
        confmat
    end
end

function extract_mut_pred(pred_scores,bucket_n,normalize)

    # columns prediction

    predc = map(p->p[9],pred_scores)[bucket_n]
    predp = map(p->p[10],pred_scores)[bucket_n];
    mut = map(p->p[11],pred_scores)[bucket_n];

    # predcl = map(v->map(x->x==-1 ? 5 : vertex_to_predict_label[x],v),predc)

    predclp = [[p ? l : n==l ? 6 : l for (p,l) in zip(persist,class_labels)] for (n,(persist,class_labels)) in enumerate(zip(predp,predcl))]

    mutr = [[m[findall(x->x==n,class_pred)] for n in 1:6] for (m,class_pred) in zip(mut,predclp)]

end