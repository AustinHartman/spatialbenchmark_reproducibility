rule all:
    input:
        ### Objects ###
        "data/linnarsson_fullbrain_reference.rds",
        "data/vizgen.rds",
        "data/vizgen_rep2.rds",
        "data/vizgen_rep3.rds",
        "data/tenx.rds",
        "data/tenx_rep2.rds",
        "data/tenx_rep3.rds",        
        "data/resolve.rds",
        "data/eelfish.rds",
        "data/merfish.rds",
        "data/starmapplus.rds",
        "data/vizgen_baysor_cortex.rds",
        "data/tenx_baysor_cortex.rds",
        "data/resolve_baysor_cortex.rds",
        "data/eelfish_baysor_cortex.rds",
        "data/merfish_baysor_cortex.rds",
        "data/vizgen_baysor_cortex_annotated.rds",
        "data/vizgen_rep2_baysor_cortex_annotated.rds",
        "data/tenx_baysor_cortex_annotated.rds",
        "data/tenx_rep2_baysor_cortex_annotated.rds",
        "data/resolve_baysor_cortex_annotated.rds",
        "data/eelfish_baysor_cortex_annotated.rds",
        "data/merfish_baysor_cortex_annotated.rds",
        "data/merfish_rep2_baysor_cortex_annotated.rds",
            
        ### Plots ###
        "data/replicate_plot_vizgen.png",
        "data/replicate_plot_tenx.png",
        "data/replicate_plot_merfish.png",
        "data/replicate_plot_sagittal_starmapplus.png",
        "data/replicate_plot_coronal_starmapplus.png",
        "data/figure2_exclusive_expression_Gfap_Slc17a7_vizgen.png",
        "data/figure2_exclusive_expression_Gfap_Slc17a7_tenx.png",
        "data/figure2_exclusive_expression_Gfap_Slc17a7_merfish.png",
        "data/figure2_exclusive_expression_Gfap_Slc17a7_resolve.png",
        "data/figure2_exclusive_expression_Gfap_Slc17a7_starmapplus.png",
        "data/figure2_exclusive_expression_Gfap_Slc17a7_eelfish.png",
        "data/figure3_mecr_mean_counts_plot.png",
        "data/figure3_mecr_median_counts_plot.png",
        "data/figure3_mecr_mean_features_plot.png",
        "data/figure3_mecr_median_features_plot.png",
        "data/figure3_mecr_threshold_mean_counts_plot.png",
        "data/exclusive_expression_Gfap_Slc17a7_linnarsson.png",
        "data/exclusive_expression_Gfap_Slc17a7_allenbrain.png",
        "data/expression_histogram_tenx.png",
        "data/expression_histogram_vizgen.png",
        "data/expression_histogram_resolve.png",
        "data/expression_histogram_eelfish.png",
        "data/expression_histogram_starmapplus.png",
        "data/expression_histogram_merfish.png",
        "data/heatmap_tenx.png",
        "data/heatmap_vizgen.png",
        "data/heatmap_resolve.png",
        "data/heatmap_eelfish.png",
        "data/heatmap_starmapplus.png",
        "data/heatmap_merfish.png",
        "data/heatmap_subset_tenx.png",
        "data/heatmap_subset_vizgen.png",
        "data/heatmap_subset_resolve.png",
        "data/heatmap_subset_merfish.png",
        "data/heatmap_subset_eelfish.png",
        "data/heatmap_subset_starmapplus.png",
        "data/heatmap_subset_linnarsson.png",
        "data/average_expression_Astrocyte_tenx_no_color.png",
        "data/average_expression_Astrocyte_vizgen_no_color.png",
        "data/average_expression_Astrocyte_resolve_no_color.png",
        "data/average_expression_Astrocyte_eelfish_no_color.png",
        "data/average_expression_Astrocyte_merfish_no_color.png",
        "data/average_expression_Astrocyte_starmapplus_no_color.png",
        "data/figure1_umaps.png",
        "data/figure1_umaps_taxonomy_group.png",
        "data/figure1_spatial_plot.png",
        "data/supplementary_figure_panel_overlap_upset_plot.png",
        "data/supplementary_figure_panel_size_by_features_per_cell_scatter.png",
        "data/figure3_exclusive_expression_Gad1_Slc17a7_tenx_default_segs_plot.png",
        "data/figure3_exclusive_expression_Gad1_Slc17a7_vizgen_default_segs_plot.png",
        "data/figure3_mecr_vs_frac_mols_in_cells.png",
        "data/counts_per_cell.png",
        "data/features_per_cell.png",
        "data/counts_per_feature_per_cell.png",
        "data/panel_size_vs_detected_mols_per_cell.png",
        "data/features_per_cell_per_gene.png",
        "data/fraction_molecules_in_cells.png",
        "data/figure3_exclusive_expression_Gad1_Slc17a7_vizgen_baysor.png",
        "data/figure3_exclusive_expression_Gad1_Slc17a7_tenx_baysor.png",
        "data/fraction_molecules_in_cells.png",
        "data/figure4_scrna_violin.png",
        "data/figure4_scrna_dotplot.png",
        "data/figure4_tenx_volcanoplot.png",
        "data/figure4_Slc17a6_Aqp4_molecules.png",
        "data/figure4_Slc17a7_Aqp4_molecules.png",
        "data/figure4_tenx_violinplot.png",
        "data/supfig3_rank_counts_per_cell_vizgen.png",
        "data/supfig3_rank_counts_per_cell_tenx.png",
        "data/supfig3_rank_counts_per_cell_resolve.png",
        "data/supfig3_rank_counts_per_cell_eelfish.png",
        "data/supfig3_rank_counts_per_cell_merfish.png",
        "data/supfig3_background_ratio_plot.png",
        "data/supfig_tenx_volcanoplot.png",
        "data/supfig_tenx_violinplot.png",
        "data/supfig_vizgen_volcanoplot.png",
        "data/supfig_vizgen_violinplot.png",
        "data/supfig_vizgen_default_violinplot.png",
        "data/supfig_vizgen_default_volcanoplot.png",
        "data/supfig_vizgen_default_dotplot.png",
        "data/supfig_vizgen_baysor_dotplot.png",
        "data/supfig_tenx_baysor_dotplot.png",
        "data/supp_mecr_mean_counts_plot.png",
        "data/supp_mecr_median_counts_plot.png",
        "data/supp_mecr_mean_features_plot.png",
        "data/supp_mecr_median_features_plot.png",
        "data/supp_mecr_threshold_mean_counts_plot.png",
        "data/supp_mecr_mean_counts_thalamus_plot.png",
        "data/supp_mecr_median_counts_thalamus_plot.png",
        "data/supp_mecr_mean_features_thalamus_plot.png",
        "data/supp_mecr_median_features_thalamus_plot.png",
        "data/supp_mecr_threshold_mean_counts_thalamus_plot.png",

        ### CSVs ###
        "data/vizgen_cortex_molecules.csv",
        "data/tenx_cortex_molecules.csv",
        "data/resolve_cortex_molecules.csv",
        "data/eelfish_cortex_molecules.csv",
        "data/merfish_cortex_molecules.csv",
        "data/figure2_mecr_boxplot.png",
        "data/wilcox_class_markers_linnarsson.csv",
        "data/tenx_baysor_cortex_annotations.csv",
        "data/tenx_baysor_cortex_markers.csv",
        "data/vizgen_baysor_cortex_annotations.csv",
        "data/vizgen_baysor_cortex_markers.csv",
        "data/resolve_baysor_cortex_annotations.csv",
        "data/resolve_baysor_cortex_markers.csv",
        "data/eelfish_baysor_cortex_annotations.csv",
        "data/eelfish_baysor_cortex_markers.csv",
        "data/merfish_baysor_cortex_annotations.csv",
        "data/merfish_baysor_cortex_markers.csv",

        ### Miscellaneous ###
        "data/baysor_vizgen/",
        "data/baysor_tenx/",
        "data/baysor_resolve/",
        "data/baysor_eelfish/",
        "data/baysor_merfish",

rule build_linnarsson_reference:
    input:
        script = "scripts/linnarsson_reference_construction.R",
        loom = "/brahms/hartmana/spatial_sensitivity_comparison/linnarsson_mouse_brain_reference/l5_all.loom",
    output:
        save_path = "data/linnarsson_fullbrain_reference.rds"
    shell:
        """
        Rscript {input.script} {input.loom} {output.save_path}
        """

rule load_seurat_vizgen:
    input:
        script = "scripts/load_seurat.R",
        vizgen = "/brahms/hartmana/spatial_sensitivity_comparison/vizgen_mouse_brain/mouse_brain_v1/s2r1",
    output:
        "data/vizgen.rds"
    shell:
        """
        Rscript {input.script} {input.vizgen} vizgen
        """

rule load_seurat_vizgen_rep2:
    input:
        script = "scripts/load_seurat.R",
        vizgen = "/brahms/hartmana/spatial_sensitivity_comparison/vizgen_mouse_brain/mouse_brain_v1/s2r2",
    output:
        "data/vizgen_rep2.rds"
    shell:
        """
        Rscript {input.script} {input.vizgen} vizgen vizgen_rep2
        """

rule load_seurat_vizgen_rep3:
    input:
        script = "scripts/load_seurat.R",
        vizgen = "/brahms/hartmana/spatial_sensitivity_comparison/vizgen_mouse_brain/mouse_brain_v1/s2r3",
    output:
        "data/vizgen_rep3.rds"
    shell:
        """
        Rscript {input.script} {input.vizgen} vizgen vizgen_rep3
        """

rule load_seurat_tenx:
    input:
        script = "scripts/load_seurat.R",
        tenx = "/brahms/hartmana/spatial_sensitivity_comparison/10x_xenium/fresh_frozen_mouse_brain_rep1",
    output:
        "data/tenx.rds"
    shell:
        """
        Rscript {input.script} {input.tenx} tenx
        """

rule load_seurat_tenx_rep2:
    input:
        script = "scripts/load_seurat.R",
        tenx = "/brahms/hartmana/spatial_sensitivity_comparison/10x_xenium/fresh_frozen_mouse_brain_rep2",
    output:
        "data/tenx_rep2.rds"
    shell:
        """
        Rscript {input.script} {input.tenx} tenx tenx_rep2
        """

rule load_seurat_tenx_rep3:
    input:
        script = "scripts/load_seurat.R",
        tenx = "/brahms/hartmana/spatial_sensitivity_comparison/10x_xenium/fresh_frozen_mouse_brain_rep3",
    output:
        "data/tenx_rep3.rds"
    shell:
        """
        Rscript {input.script} {input.tenx} tenx tenx_rep3
        """

rule load_seurat_resolve:
    input:
        script = "scripts/load_seurat.R",
        resolve = "/brahms/hartmana/spatial_sensitivity_comparison/Resolve_Mm_BHA1",
    output:
        "data/resolve.rds"
    shell:
        """
        Rscript {input.script} {input.resolve} resolve
        """

rule load_seurat_eelfish:
    input:
        script = "scripts/load_seurat.R",
        eelfish = "/brahms/hartmana/spatial_sensitivity_comparison/eel_fish_latest/eel_fish.h5ad",
    output:
        "data/eelfish.rds"
    shell:
        """
        Rscript {input.script} {input.eelfish} eelfish
        """

rule load_seurat_merfish:
    input:
        script = "scripts/load_seurat.R",
        merfish = "/brahms/hartmana/spatial_sensitivity_comparison/2023_merfish_mouse_brain/mouse_coronal_2/raw_counts_mouse2_coronal.h5ad"
    output:
        "data/merfish.rds"
    shell:
        """
        Rscript {input.script} {input.merfish} merfish merfish 21
        """

rule load_seurat_merfish_rep2:
    input:
        script = "scripts/load_seurat.R",
        merfish = "/brahms/hartmana/spatial_sensitivity_comparison/2023_merfish_mouse_brain/mouse_coronal_2/raw_counts_mouse2_coronal.h5ad"
    output:
        "data/merfish_rep2.rds"
    shell:
        """
        Rscript {input.script} {input.merfish} merfish merfish_rep2 22
        """

rule load_seurat_starmapplus:
    input:
        script = "scripts/load_seurat.R",
        starmapplus = "/brahms/hartmana/spatial_sensitivity_comparison/STARmap_plus_data/sagittal2"
    output:
        "data/starmapplus.rds"
    shell:
        """
        Rscript {input.script} {input.starmapplus} starmapplus
        """

rule load_seurat_starmapplus_rep2:
    input:
        script = "scripts/load_seurat.R",
        starmapplus = "/brahms/hartmana/spatial_sensitivity_comparison/STARmap_plus_data/sagittal3"
    output:
        "data/starmapplus_rep2.rds"
    shell:
        """
        Rscript {input.script} {input.starmapplus} starmapplus starmapplus_rep2
        """

rule load_seurat_starmapplus_rep3:
    input:
        script = "scripts/load_seurat.R",
        starmapplus = "/brahms/hartmana/spatial_sensitivity_comparison/STARmap_plus_data/well05"
    output:
        "data/starmapplus_rep3.rds"
    shell:
        """
        Rscript {input.script} {input.starmapplus} starmapplus starmapplus_rep3
        """

rule load_seurat_starmapplus_rep4:
    input:
        script = "scripts/load_seurat.R",
        starmapplus = "/brahms/hartmana/spatial_sensitivity_comparison/STARmap_plus_data/well11"
    output:
        "data/starmapplus_rep4.rds"
    shell:
        """
        Rscript {input.script} {input.starmapplus} starmapplus starmapplus_rep4
        """

rule convert_molecules_vizgen:
    input:
        script = "scripts/convert_molecules.py",
        vizgen = "/brahms/hartmana/spatial_sensitivity_comparison/vizgen_mouse_brain/mouse_brain_v1/s2r1/detected_transcripts_S2R1.csv",
    output:
        "data/vizgen_molecules.csv"
    conda:
        "base"
    shell:
        """
        python {input.script} {input.vizgen} vizgen
        """

rule convert_molecules_vizgen_rep2:
    input:
        script = "scripts/convert_molecules.py",
        vizgen = "/brahms/hartmana/spatial_sensitivity_comparison/vizgen_mouse_brain/mouse_brain_v1/s2r2/detected_transcripts_S2R2.csv",
    output:
        "data/vizgen_rep2_molecules.csv"
    conda:
        "base"
    shell:
        """
        python {input.script} {input.vizgen} vizgen_rep2
        """

rule convert_molecules_vizgen_rep3:
    input:
        script = "scripts/convert_molecules.py",
        vizgen = "/brahms/hartmana/spatial_sensitivity_comparison/vizgen_mouse_brain/mouse_brain_v1/s2r3/detected_transcripts_S2R3.csv",
    output:
        "data/vizgen_rep3_molecules.csv"
    conda:
        "base"
    shell:
        """
        python {input.script} {input.vizgen} vizgen_rep3
        """

rule convert_molecules_tenx:
    input:
        script = "scripts/convert_molecules.py",
        tenx = "/brahms/hartmana/spatial_sensitivity_comparison/10x_xenium/fresh_frozen_mouse_brain_rep1/transcripts.csv.gz",
    output:
        "data/tenx_molecules.csv"
    conda:
        "base"
    shell:
        """
        python {input.script} {input.tenx} tenx
        """

rule convert_molecules_tenx_rep2:
    input:
        script = "scripts/convert_molecules.py",
        tenx = "/brahms/hartmana/spatial_sensitivity_comparison/10x_xenium/fresh_frozen_mouse_brain_rep2/transcripts.csv.gz",
    output:
        "data/tenx_rep2_molecules.csv"
    conda:
        "base"
    shell:
        """
        python {input.script} {input.tenx} tenx_rep2
        """

rule convert_molecules_tenx_rep3:
    input:
        script = "scripts/convert_molecules.py",
        tenx = "/brahms/hartmana/spatial_sensitivity_comparison/10x_xenium/fresh_frozen_mouse_brain_rep3/transcripts.csv.gz",
    output:
        "data/tenx_rep3_molecules.csv"
    conda:
        "base"
    shell:
        """
        python {input.script} {input.tenx} tenx_rep3
        """

rule convert_molecules_resolve:
    input:
        script = "scripts/convert_molecules.py",
        resolve = "/brahms/hartmana/spatial_sensitivity_comparison/Resolve_Mm_BHA1/Resolve_Mm_BHA1_raw_locations.txt",
    output:
        "data/resolve_molecules.csv"
    conda:
        "base"
    shell:
        """
        python {input.script} {input.resolve} resolve
        """

rule convert_molecules_eelfish:
    input:
        script = "scripts/convert_molecules.py",
        eelfish = "/brahms/hartmana/spatial_sensitivity_comparison/eel_fish_latest/20220821_LBEXP20210718_EEL_Mouse_448_2_RNA.csv",
    output:
        "data/eelfish_molecules.csv"
    conda:
        "base"
    shell:
        """
        python {input.script} {input.eelfish} eelfish
        """

rule convert_molecules_merfish:
    input:
        script = "scripts/convert_molecules.py",
        merfish = "/brahms/hartmana/spatial_sensitivity_comparison/2023_merfish_mouse_brain/mouse_coronal_2/spots_220605_co2_11B_merfish4_adaptor.csv",
    output:
        "data/merfish_molecules.csv"
    conda:
        "base"
    shell:
        """
        python {input.script} {input.merfish} merfish
        """

rule convert_molecules_merfish_rep2:
    input:
        script = "scripts/convert_molecules.py",
        merfish = "/brahms/hartmana/spatial_sensitivity_comparison/2023_merfish_mouse_brain/mouse_coronal_2/spots_220507_wb3_co2_11_5z18R_merfish5.csv",
    output:
        "data/merfish_rep2_molecules.csv"
    conda:
        "base"
    shell:
        """
        python {input.script} {input.merfish} merfish_rep2
        """

rule annotate_seurat_vizgen:
    input:
        script = "scripts/annotate_seurat.R",
        reference = "data/linnarsson_fullbrain_reference.rds",
        vizgen = "data/vizgen.rds",
    output:
        obj = "data/vizgen_annotated.rds"
    shell:
        """
        Rscript {input.script} {input.reference} {input.vizgen} {output.obj}
        """

rule annotate_seurat_vizgen_rep2:
    input:
        script = "scripts/annotate_seurat.R",
        reference = "data/linnarsson_fullbrain_reference.rds",
        vizgen = "data/vizgen_rep2.rds",
    output:
        obj = "data/vizgen_rep2_annotated.rds"
    shell:
        """
        Rscript {input.script} {input.reference} {input.vizgen} {output.obj}
        """

rule annotate_seurat_vizgen_rep3:
    input:
        script = "scripts/annotate_seurat.R",
        reference = "data/linnarsson_fullbrain_reference.rds",
        vizgen = "data/vizgen_rep3.rds",
    output:
        obj = "data/vizgen_rep3_annotated.rds"
    shell:
        """
        Rscript {input.script} {input.reference} {input.vizgen} {output.obj}
        """

rule annotate_seurat_tenx:
    input:
        script = "scripts/annotate_seurat.R",
        reference = "data/linnarsson_fullbrain_reference.rds",
        seurat_object = "data/tenx.rds"
    output:
        obj = "data/tenx_annotated.rds"
    shell:
        """
        Rscript {input.script} {input.reference} {input.seurat_object} {output.obj}
        """

rule annotate_seurat_tenx_rep2:
    input:
        script = "scripts/annotate_seurat.R",
        reference = "data/linnarsson_fullbrain_reference.rds",
        seurat_object = "data/tenx_rep2.rds"
    output:
        obj = "data/tenx_rep2_annotated.rds"
    shell:
        """
        Rscript {input.script} {input.reference} {input.seurat_object} {output.obj}
        """

rule annotate_seurat_tenx_rep3:
    input:
        script = "scripts/annotate_seurat.R",
        reference = "data/linnarsson_fullbrain_reference.rds",
        seurat_object = "data/tenx_rep3.rds"
    output:
        obj = "data/tenx_rep3_annotated.rds"
    shell:
        """
        Rscript {input.script} {input.reference} {input.seurat_object} {output.obj}
        """

rule annotate_seurat_resolve:
    input:
        script = "scripts/annotate_seurat.R",
        reference = "data/linnarsson_fullbrain_reference.rds",
        seurat_object = "data/resolve.rds",
    output:
        obj = "data/resolve_annotated.rds"
    shell:
        """
        Rscript {input.script} {input.reference} {input.seurat_object} {output.obj}
        """

rule annotate_seurat_eelfish:
    input:
        script = "scripts/annotate_seurat.R",
        reference = "data/linnarsson_fullbrain_reference.rds",
        seurat_object = "data/eelfish.rds",
    output:
        obj = "data/eelfish_annotated.rds"
    shell:
        """
        Rscript {input.script} {input.reference} {input.seurat_object} {output.obj}
        """

rule annotate_seurat_merfish:
    input:
        script = "scripts/annotate_seurat.R",
        reference = "data/linnarsson_fullbrain_reference.rds",
        seurat_object = "data/merfish.rds",
    output:
        obj = "data/merfish_annotated.rds"
    shell:
        """
        Rscript {input.script} {input.reference} {input.seurat_object} {output.obj}
        """

rule annotate_seurat_starmapplus:
    input:
        script = "scripts/annotate_seurat.R",
        reference = "data/linnarsson_fullbrain_reference.rds",
        seurat_object = "data/starmapplus.rds",
    output:
        obj = "data/starmapplus_annotated.rds"
    shell:
        """
        Rscript {input.script} {input.reference} {input.seurat_object} {output.obj}
        """

rule replicate_plot_vizgen:
    input:
        script = "scripts/replicate_plot.R",
        rep1 = "data/vizgen.rds",
        rep2 = "data/vizgen_rep2.rds",
    output:
        "data/replicate_plot_vizgen.png"
    shell:
        """
        Rscript {input.script} {input.rep1} {input.rep2} vizgen
        """

rule replicate_plot_tenx:
    input:
        script = "scripts/replicate_plot.R",
        rep1 = "data/tenx.rds",
        rep2 = "data/tenx_rep2.rds",
    output:
        "data/replicate_plot_tenx.png"
    shell:
        """
        Rscript {input.script} {input.rep1} {input.rep2} tenx
        """

rule replicate_plot_merfish:
    input:
        script = "scripts/replicate_plot.R",
        rep1 = "data/merfish.rds",
        rep2 = "data/merfish_rep2.rds",
    output:
        "data/replicate_plot_merfish.png"
    shell:
        """
        Rscript {input.script} {input.rep1} {input.rep2} merfish
        """

rule replicate_plot_sagittal_starmapplus:
    input:
        script = "scripts/replicate_plot.R",
        rep1 = "data/starmapplus.rds",
        rep2 = "data/starmapplus_rep2.rds",
    output:
        "data/replicate_plot_sagittal_starmapplus.png"
    shell:
        """
        Rscript {input.script} {input.rep1} {input.rep2} sagittal_starmapplus
        """

rule replicate_plot_coronal_starmapplus:
    input:
        script = "scripts/replicate_plot.R",
        rep1 = "data/starmapplus_rep3.rds",
        rep2 = "data/starmapplus_rep4.rds",
    output:
        "data/replicate_plot_coronal_starmapplus.png"
    shell:
        """
        Rscript {input.script} {input.rep1} {input.rep2} coronal_starmapplus
        """

rule panel_overlap_upset_plot:
    input:
        script = "scripts/panel_overlap_upset_plot.R",
        vizgen = "data/vizgen_annotated.rds",
        tenx = "data/tenx_annotated.rds",
        resolve = "data/resolve_annotated.rds",
        eelfish = "data/eelfish_annotated.rds",
        merfish = "data/merfish_annotated.rds",
        starmap = "data/starmapplus_annotated.rds",
    output:
        "data/supplementary_figure_panel_overlap_upset_plot.png",
    shell:
        """
        Rscript {input.script} {input.vizgen} {input.tenx} {input.resolve} {input.eelfish} {input.merfish} {input.starmap} {output}
        """

rule figure1_plots:
    input:
        script = "scripts/figure_1_plots.R",
        vizgen = "data/vizgen_annotated.rds",
        tenx = "data/tenx_annotated.rds",
        resolve = "data/resolve_annotated.rds",
        eelfish = "data/eelfish_annotated.rds",
        merfish = "data/merfish_annotated.rds",
        starmapplus = "data/starmapplus_annotated.rds",
    output:
        "data/figure1_umaps.png",
        "data/figure1_umaps_taxonomy_group.png",
    shell:
        """
        Rscript {input.script} {input.vizgen} {input.tenx} {input.resolve} {input.eelfish} {input.merfish} {input.starmapplus} {output}
        """

rule figure1_spatial_plots:
    input:
        script = "scripts/figure_1_spatial_plots.R",
        vizgen = "data/vizgen_annotated.rds",
        vizgen_spatial = "data/vizgen.rds",
        tenx = "data/tenx_annotated.rds",
        tenx_spatial = "data/tenx.rds",
        resolve = "data/resolve_annotated.rds",
        resolve_spatial = "data/resolve.rds",
        eelfish = "data/eelfish_annotated.rds",
        eelfish_spatial = "data/eelfish.rds",
        merfish = "data/merfish_annotated.rds",
        merfish_spatial = "data/merfish.rds",
        starmapplus = "data/starmapplus_annotated.rds",
        starmapplus_spatial = "data/starmapplus.rds",
    output:
        "data/figure1_spatial_plot.png",
    shell:
        """
        Rscript {input.script} {input.vizgen} {input.vizgen_spatial} {input.tenx} {input.tenx_spatial} \
        {input.resolve} {input.resolve_spatial} {input.eelfish} {input.eelfish_spatial} \
        {input.merfish} {input.merfish_spatial} {input.starmapplus} {input.starmapplus_spatial} {output}
        """

rule panel_size_by_features_per_cell_scatter:
    input:
        script = "scripts/panel_size_by_features_per_cell_scatter.R",
        vizgen = "data/vizgen_annotated.rds",
        tenx = "data/tenx_annotated.rds",
        resolve = "data/resolve_annotated.rds",
        eelfish = "data/eelfish_annotated.rds",
        merfish = "data/merfish_annotated.rds",
        starmapplus = "data/starmapplus_annotated.rds",
    output:
        "data/supplementary_figure_panel_size_by_features_per_cell_scatter.png",
    shell:
        """
        Rscript {input.script} {input.vizgen} {input.tenx} {input.resolve} {input.eelfish} {input.merfish} {input.starmapplus} {output}
        """

rule single_cell_expression_histogram_vizgen:
    input:
        script = "scripts/expression_histograms.R",
        ref = "data/linnarsson_fullbrain_reference.rds",
        obj = "data/vizgen_annotated.rds",
    output:
        histogram = "data/expression_histogram_vizgen.png",
        violin = "data/expression_violin_vizgen_stats.png",
    shell:
        """
        Rscript {input.script} {input.ref} {input.obj} {output.histogram} {output.violin}
        """

rule single_cell_expression_histogram_tenx:
    input:
        script = "scripts/expression_histograms.R",
        ref = "data/linnarsson_fullbrain_reference.rds",
        obj = "data/tenx_annotated.rds",
    output:
        histogram = "data/expression_histogram_tenx.png",
        violin = "data/expression_violin_tenx_stats.png",
    shell:
        """
        Rscript {input.script} {input.ref} {input.obj} {output.histogram} {output.violin}
        """

rule single_cell_expression_histogram_resolve:
    input:
        script = "scripts/expression_histograms.R",
        ref = "data/linnarsson_fullbrain_reference.rds",
        obj = "data/resolve_annotated.rds",
    output:
        histogram = "data/expression_histogram_resolve.png",
        violin = "data/expression_violin_resolve_stats.png",
    shell:
        """
        Rscript {input.script} {input.ref} {input.obj} {output.histogram} {output.violin}
        """

rule single_cell_expression_histogram_eelfish:
    input:
        script = "scripts/expression_histograms.R",
        ref = "data/linnarsson_fullbrain_reference.rds",
        obj = "data/eelfish_annotated.rds",
    output:
        histogram = "data/expression_histogram_eelfish.png",
        violin = "data/expression_violin_eelfish_stats.png",
    shell:
        """
        Rscript {input.script} {input.ref} {input.obj} {output.histogram} {output.violin}
        """

rule single_cell_expression_histogram_starmapplus:
    input:
        script = "scripts/expression_histograms.R",
        ref = "data/linnarsson_fullbrain_reference.rds",
        obj = "data/starmapplus_annotated.rds",
    output:
        histogram = "data/expression_histogram_starmapplus.png",
        violin = "data/expression_violin_starmapplus_stats.png",
    shell:
        """
        Rscript {input.script} {input.ref} {input.obj} {output.histogram} {output.violin}
        """

rule single_cell_expression_histogram_merfish:
    input:
        script = "scripts/expression_histograms.R",
        ref = "data/linnarsson_fullbrain_reference.rds",
        obj = "data/merfish_annotated.rds",
    output:
        histogram = "data/expression_histogram_merfish.png",
        violin = "data/expression_violin_merfish_stats.png",
    shell:
        """
        Rscript {input.script} {input.ref} {input.obj} {output.histogram} {output.violin}
        """

rule figure2_exclusive_expression_Gfap_Slc17a7_vizgen:
    input:
        script = "scripts/exclusive_expression_plot.R",
        seurat_object = "data/vizgen_annotated.rds"
    output:
        plot = "data/figure2_exclusive_expression_Gfap_Slc17a7_vizgen.png"
    shell:
        """
        Rscript {input.script} {input.seurat_object} vizgen Gfap Slc17a7 300 300 {output.plot}
        """

rule figure2_exclusive_expression_Gfap_Slc17a7_tenx:
    input:
        script = "scripts/exclusive_expression_plot.R",
        seurat_object = "data/tenx_annotated.rds"
    output:
        plot = "data/figure2_exclusive_expression_Gfap_Slc17a7_tenx.png"
    shell:
        """
        Rscript {input.script} {input.seurat_object} tenx Gfap Slc17a7 300 300 {output.plot}
        """

rule figure2_exclusive_expression_Gfap_Slc17a7_merfish:
    input:
        script = "scripts/exclusive_expression_plot.R",
        seurat_object = "data/merfish_annotated.rds"
    output:
        plot = "data/figure2_exclusive_expression_Gfap_Slc17a7_merfish.png"
    shell:
        """
        Rscript {input.script} {input.seurat_object} merfish Gfap Slc17a7 300 300 {output.plot}
        """

rule figure2_exclusive_expression_Gfap_Slc17a7_resolve:
    input:
        script = "scripts/exclusive_expression_plot.R",
        seurat_object = "data/resolve_annotated.rds"
    output:
        plot = "data/figure2_exclusive_expression_Gfap_Slc17a7_resolve.png"
    shell:
        """
        Rscript {input.script} {input.seurat_object} resolve Gfap Slc17a7 300 300 {output.plot}
        """

rule figure2_exclusive_expression_Gfap_Slc17a7_eelfish:
    input:
        script = "scripts/exclusive_expression_plot.R",
        seurat_object = "data/eelfish_annotated.rds"
    output:
        plot = "data/figure2_exclusive_expression_Gfap_Slc17a7_eelfish.png"
    shell:
        """
        Rscript {input.script} {input.seurat_object} eelfish Gfap Slc17a7 300 300 {output.plot}
        """

rule figure2_exclusive_expression_Gfap_Slc17a7_starmapplus:
    input:
        script = "scripts/exclusive_expression_plot.R",
        seurat_object = "data/starmapplus_annotated.rds"
    output:
        plot = "data/figure2_exclusive_expression_Gfap_Slc17a7_starmapplus.png"
    shell:
        """
        Rscript {input.script} {input.seurat_object} starmapplus Gfap Slc17a7 300 300 {output.plot}
        """

rule plot_exclusive_expression_Gfap_Slc17a7_linnarsson:
    input:
        script = "scripts/exclusive_expression_plot.R",
        seurat_object = "data/linnarsson_fullbrain_reference.rds"
    output:
        plot = "data/exclusive_expression_Gfap_Slc17a7_linnarsson.png"
    shell:
        """
        Rscript {input.script} {input.seurat_object} linnarsson Gfap Slc17a7 100 50 black 2 {output.plot}
        """

rule plot_exclusive_expression_Gfap_Slc17a7_allenbrain:
    input:
        script = "scripts/exclusive_expression_plot.R",
        seurat_object = "/brahms/hartmana/vignette_data/allen_cortex.rds",
        # seurat_object = "/brahms/hartmana/spatial_sensitivity_comparison/allen_brain_mouse_whole_cortex_and_hippo_10x/sketched_allen_brain_10x_v4.rds"
    output:
        plot = "data/exclusive_expression_Gfap_Slc17a7_allenbrain.png"
    shell:
        """
        Rscript {input.script} {input.seurat_object} allenbrain Gfap Slc17a7 1000 5000 black 2 {output.plot}
        """

rule plot_heatmap_tenx:
    input:
        script = "scripts/annotation_heatmap.R",
        seurat_object = "data/tenx_annotated.rds"
    output:
        "data/heatmap_tenx.png"
    shell:
        """
        Rscript {input.script} {input.seurat_object} tenx
        """

rule plot_heatmap_vizgen:
    input:
        script = "scripts/annotation_heatmap.R",
        seurat_object = "data/vizgen_annotated.rds"
    output:
        "data/heatmap_vizgen.png"
    shell:
        """
        Rscript {input.script} {input.seurat_object} vizgen
        """

rule plot_heatmap_resolve:
    input:
        script = "scripts/annotation_heatmap.R",
        seurat_object = "data/resolve_annotated.rds"
    output:
        "data/heatmap_resolve.png"
    shell:
        """
        Rscript {input.script} {input.seurat_object} resolve
        """

rule plot_heatmap_merfish:
    input:
        script = "scripts/annotation_heatmap.R",
        seurat_object = "data/merfish_annotated.rds"
    output:
        "data/heatmap_merfish.png"
    shell:
        """
        Rscript {input.script} {input.seurat_object} merfish
        """

rule plot_heatmap_eelfish:
    input:
        script = "scripts/annotation_heatmap.R",
        seurat_object = "data/eelfish_annotated.rds"
    output:
        "data/heatmap_eelfish.png"
    shell:
        """
        Rscript {input.script} {input.seurat_object} eelfish
        """

rule plot_heatmap_starmapplus:
    input:
        script = "scripts/annotation_heatmap.R",
        seurat_object = "data/starmapplus_annotated.rds"
    output:
        "data/heatmap_starmapplus.png"
    shell:
        """
        Rscript {input.script} {input.seurat_object} starmapplus
        """

rule plot_heatmap_subset_linnarsson:
    input:
        script = "scripts/annotation_heatmap_some_cells.R",
        seurat_object = "data/linnarsson_fullbrain_reference.rds"
    output:
        "data/heatmap_subset_linnarsson.png"
    shell:
        """
        Rscript {input.script} {input.seurat_object} linnarsson
        """

rule plot_heatmap_subset_tenx:
    input:
        script = "scripts/annotation_heatmap_some_cells.R",
        seurat_object = "data/tenx_annotated.rds"
    output:
        "data/heatmap_subset_tenx.png"
    shell:
        """
        Rscript {input.script} {input.seurat_object} tenx
        """

rule plot_heatmap_subset_vizgen:
    input:
        script = "scripts/annotation_heatmap_some_cells.R",
        seurat_object = "data/vizgen_annotated.rds"
    output:
        "data/heatmap_subset_vizgen.png"
    shell:
        """
        Rscript {input.script} {input.seurat_object} vizgen
        """

rule plot_heatmap_subset_resolve:
    input:
        script = "scripts/annotation_heatmap_some_cells.R",
        seurat_object = "data/resolve_annotated.rds"
    output:
        "data/heatmap_subset_resolve.png"
    shell:
        """
        Rscript {input.script} {input.seurat_object} resolve
        """

rule plot_heatmap_subset_merfish:
    input:
        script = "scripts/annotation_heatmap_some_cells.R",
        seurat_object = "data/merfish_annotated.rds"
    output:
        "data/heatmap_subset_merfish.png"
    shell:
        """
        Rscript {input.script} {input.seurat_object} merfish
        """

rule plot_heatmap_subset_eelfish:
    input:
        script = "scripts/annotation_heatmap_some_cells.R",
        seurat_object = "data/eelfish_annotated.rds"
    output:
        "data/heatmap_subset_eelfish.png"
    shell:
        """
        Rscript {input.script} {input.seurat_object} eelfish
        """

rule plot_heatmap_subset_starmapplus:
    input:
        script = "scripts/annotation_heatmap_some_cells.R",
        seurat_object = "data/starmapplus_annotated.rds"
    output:
        "data/heatmap_subset_starmapplus.png"
    shell:
        """
        Rscript {input.script} {input.seurat_object} starmapplus
        """

rule plot_average_expression_astrocyte_vizgen:
    input:
        script = "scripts/average_expression_plot.R",
        seurat_object = "data/vizgen_annotated.rds",
        reference = "data/linnarsson_fullbrain_reference.rds",
    output:
        "data/average_expression_Astrocyte_vizgen_no_color.png"
    shell:
        """
        Rscript {input.script} {input.seurat_object} {input.reference} vizgen Astrocyte FALSE
        """

rule plot_average_expression_astrocyte_tenx:
    input:
        script = "scripts/average_expression_plot.R",
        seurat_object = "data/tenx_annotated.rds",
        reference = "data/linnarsson_fullbrain_reference.rds",
    output:
        "data/average_expression_Astrocyte_tenx_no_color.png"
    shell:
        """
        Rscript {input.script} {input.seurat_object} {input.reference} tenx Astrocyte FALSE
        """

rule plot_average_expression_astrocyte_merfish:
    input:
        script = "scripts/average_expression_plot.R",
        seurat_object = "data/merfish_annotated.rds",
        reference = "data/linnarsson_fullbrain_reference.rds",
    output:
        "data/average_expression_Astrocyte_merfish_no_color.png"
    shell:
        """
        Rscript {input.script} {input.seurat_object} {input.reference} merfish Astrocyte FALSE
        """

rule plot_average_expression_astrocyte_resolve:
    input:
        script = "scripts/average_expression_plot.R",
        seurat_object = "data/resolve_annotated.rds",
        reference = "data/linnarsson_fullbrain_reference.rds",
    output:
        "data/average_expression_Astrocyte_resolve_no_color.png"
    shell:
        """
        Rscript {input.script} {input.seurat_object} {input.reference} resolve Astrocyte FALSE
        """

rule plot_average_expression_astrocyte_eelfish:
    input:
        script = "scripts/average_expression_plot.R",
        seurat_object = "data/eelfish_annotated.rds",
        reference = "data/linnarsson_fullbrain_reference.rds",
    output:
        "data/average_expression_Astrocyte_eelfish_no_color.png"
    shell:
        """
        Rscript {input.script} {input.seurat_object} {input.reference} eelfish Astrocyte FALSE
        """

rule plot_average_expression_astrocyte_starmapplus:
    input:
        script = "scripts/average_expression_plot.R",
        seurat_object = "data/starmapplus_annotated.rds",
        reference = "data/linnarsson_fullbrain_reference.rds",
    output:
        "data/average_expression_Astrocyte_starmapplus_no_color.png"
    shell:
        """
        Rscript {input.script} {input.seurat_object} {input.reference} starmapplus Astrocyte FALSE
        """

rule figure2_mecr_boxplot:
    input:
        script = "scripts/figure2_mecr_plot.R",
        linnarsson = "data/linnarsson_fullbrain_reference.rds",
        allen = "/brahms/hartmana/vignette_data/allen_cortex.rds",
        # allen = "/brahms/hartmana/spatial_sensitivity_comparison/allen_brain_mouse_whole_cortex_and_hippo_10x/sketched_allen_brain_10x_v4.rds",
        vizgen = "data/vizgen_annotated.rds",
        tenx = "data/tenx_annotated.rds",
        resolve = "data/resolve_annotated.rds",
        merfish = "data/merfish_annotated.rds",
        starmapplus = "data/starmapplus_annotated.rds",
        eelfish = "data/eelfish_annotated.rds",
    output:
        plot = "data/figure2_mecr_boxplot.png"
    shell:
        """
        Rscript {input.script} {input.linnarsson} {input.allen} {input.vizgen} \
        {input.tenx} {input.resolve} {input.merfish} {input.starmapplus} \
        {input.eelfish} Linnarsson Allen Vizgen 10x Resolve MERFISH STARmapPLUS \
        EELFISH {output.plot}
        """

rule linnarsson_wilcox_class_markers:
    input:
        script = "scripts/wilcox_class_markers.R",
        linnarsson = "data/linnarsson_fullbrain_reference.rds",
    output:
        class_markers = "data/wilcox_class_markers_linnarsson.csv",
        taxonomy_group_markers = "data/wilcox_taxonomy_group_markers_linnarsson.csv",
    shell:
        """
        Rscript {input.script} {input.linnarsson} {output.class_markers} {output.taxonomy_group_markers}
        """

rule plot_exclusive_expression_Gad1_Slc17a7_vizgen_baysor:
    input:
        script = "scripts/exclusive_expression_plot.R",
        seurat_object = "data/vizgen_baysor_cortex_annotated.rds",
    output:
        plot = "data/figure3_exclusive_expression_Gad1_Slc17a7_vizgen_baysor.png"
    shell:
        """
        Rscript {input.script} {input.seurat_object} vizgen_baysor \
        Gad1 Slc17a7 200 200 {output.plot}
        """

rule plot_exclusive_expression_Gad1_Slc17a7_tenx_baysor:
    input:
        script = "scripts/exclusive_expression_plot.R",
        seurat_object = "data/tenx_baysor_cortex_annotated.rds",
    output:
        plot = "data/figure3_exclusive_expression_Gad1_Slc17a7_tenx_baysor.png"
    shell:
        """
        Rscript {input.script} {input.seurat_object} tenx_baysor \
        Gad1 Slc17a7 200 200 '#3258A5' 1 {output.plot}
        """

rule qc_violin_plots:
    input:
        script = "scripts/qc_violin_plots.R",
        vizgen = "data/vizgen_annotated.rds",
        tenx = "data/tenx_annotated.rds",
        resolve = "data/resolve_annotated.rds",
        eelfish = "data/eelfish_annotated.rds",
        merfish = "data/merfish_annotated.rds",
        starmapplus = "data/starmapplus_annotated.rds",
    output:
        counts_per_cell = "data/counts_per_cell.png",
        features_per_cell = "data/features_per_cell.png",
        counts_per_feature_per_cell = "data/counts_per_feature_per_cell.png",
    shell:
        """
        Rscript {input.script} {input.vizgen} {input.tenx} {input.resolve} {input.eelfish} {input.merfish} {input.starmapplus} {output.counts_per_cell} {output.features_per_cell} {output.counts_per_feature_per_cell}
        """

rule qc_scatter_plots:
    input:
        script = "scripts/qc_scatter_plots.R",
        vizgen = "data/vizgen_annotated.rds",
        tenx = "data/tenx_annotated.rds",
        resolve = "data/resolve_annotated.rds",
        eelfish = "data/eelfish_annotated.rds",
        merfish = "data/merfish_annotated.rds",
        starmapplus = "data/starmapplus_annotated.rds",
    output:
        panel_size_vs_detected_mols_per_cell = "data/panel_size_vs_detected_mols_per_cell.png",
        panel_size_vs_detected_mols_per_gene = "data/features_per_cell_per_gene.png",
    shell:
        """
        Rscript {input.script} {input.vizgen} {input.tenx} {input.resolve} {input.eelfish} {input.merfish} {input.starmapplus} \
        {output.panel_size_vs_detected_mols_per_cell} {output.panel_size_vs_detected_mols_per_gene}
        """

rule fraction_molecules_in_cells:
    input:
        script = "scripts/qc_fraction_molecules_in_cells.R",
        vizgen = "data/vizgen_annotated.rds",
        tenx = "data/tenx_annotated.rds",
        resolve = "data/resolve_annotated.rds",
        eelfish = "data/eelfish_annotated.rds",
        merfish = "data/merfish_annotated.rds",
        vizgen_mols = "data/vizgen_molecules.csv",
        tenx_mols = "data/tenx_molecules.csv",
        resolve_mols = "data/resolve_molecules.csv",
        eelfish_mols = "data/eelfish_molecules.csv",
        merfish_mols = "data/merfish_molecules.csv",
    output:
        frac_mols_in_cells = "data/fraction_molecules_in_cells.png",
    shell:
        """
        Rscript {input.script} {input.vizgen} {input.tenx} {input.resolve} {input.eelfish} {input.merfish} \
        {input.vizgen_mols} {input.tenx_mols} {input.resolve_mols} {input.eelfish_mols} {input.merfish_mols} {output.frac_mols_in_cells}
        """

rule supfig3_qc_rank_plot_vizgen:
    input:
        script = "scripts/qc_rank_plots.R",
        molecules = "data/vizgen_molecules.csv",
    output:
        rank_counts_per_cell = "data/supfig3_rank_counts_per_cell_vizgen.png",
        zoomed = "data/supfig3_rank_counts_per_cell_vizgen_zoomed.png",
    shell:
        """
        Rscript {input.script} {input.molecules} vizgen {output.rank_counts_per_cell} {output.zoomed}
        """

rule supfig3_qc_rank_plot_tenx:
    input:
        script = "scripts/qc_rank_plots.R",
        molecules = "data/tenx_molecules.csv",
    output:
        rank_counts_per_cell = "data/supfig3_rank_counts_per_cell_tenx.png",
         zoomed = "data/supfig3_rank_counts_per_cell_tenx_zoomed.png",
    shell:
        """
        Rscript {input.script} {input.molecules} tenx {output.rank_counts_per_cell} {output.zoomed}
        """

rule supfig3_qc_rank_plot_resolve:
    input:
        script = "scripts/qc_rank_plots.R",
        molecules = "data/resolve_molecules.csv",
    output:
        rank_counts_per_cell = "data/supfig3_rank_counts_per_cell_resolve.png",
        zoomed = "data/supfig3_rank_counts_per_cell_resolve_zoomed.png",
    shell:
        """
        Rscript {input.script} {input.molecules} resolve {output.rank_counts_per_cell} {output.zoomed}
        """

rule supfig3_qc_rank_plot_eelfish:
    input:
        script = "scripts/qc_rank_plots.R",
        molecules = "data/eelfish_molecules.csv",
    output:
        rank_counts_per_cell = "data/supfig3_rank_counts_per_cell_eelfish.png",
        zoomed = "data/supfig3_rank_counts_per_cell_eelfish_zoomed.png",
    shell:
        """
        Rscript {input.script} {input.molecules} eelfish {output.rank_counts_per_cell} {output.zoomed}
        """

rule supfig3_qc_rank_plot_merfish:
    input:
        script = "scripts/qc_rank_plots.R",
        molecules = "data/merfish_molecules.csv",
    output:
        rank_counts_per_cell = "data/supfig3_rank_counts_per_cell_merfish.png",
        zoomed = "data/supfig3_rank_counts_per_cell_merfish_zoomed.png",
    shell:
        """
        Rscript {input.script} {input.molecules} merfish {output.rank_counts_per_cell} {output.zoomed}
        """

rule supfig3_qc_background_probe_stats:
    input:
        script = "scripts/supfig3_qc_background_probe_barplot.R",
        vizgen = "data/vizgen_rep3_molecules.csv",
        tenx = "data/tenx_rep3_molecules.csv",
        resolve = "data/resolve_molecules.csv",
        eelfish = "data/eelfish_molecules.csv",
        merfish = "data/merfish_molecules.csv",
    output:
        bg_ratio_plot = "data/supfig3_background_ratio_plot.png",
    shell:
        """
        Rscript {input.script} {input.vizgen} "Vizgen" {input.tenx} "10x" \
         {input.resolve} "Resolve" {input.eelfish} "EEL FISH" \
         {input.merfish} "MERFISH" {output.bg_ratio_plot}
        """

rule cell_cortex_tenx_seurat:
    input:
        script = "scripts/molecules_to_seurat.R",
        molecules = "data/tenx_cortex_molecules.csv",
    output:
        "data/tenx_cell_cortex.rds",
    shell:
        """
        Rscript {input.script} {input.molecules} tenx_cell
        """

rule plot_molecules_tenx_cortex:
    input:
        script = "scripts/plot_molecules.py",
        seurat = "data/tenx_cortex_molecules.csv",
    output:
        "data/tenx_cortex_Gad1_Slc17a7.png",
    conda:
        "base"
    shell:
        """
        python {input.script} {input.seurat} Gad1 Slc17a7 0.5 {output}
        """
    
rule baysor_cortex_vizgen:
    input:
        prepare_script = "scripts/prepare_run_baysor.py",
        run_script = "scripts/run_baysor.sh",
        molecules = "data/vizgen_molecules.csv",
    output:
        "data/vizgen_cortex_molecules.csv",
        directory("data/baysor_vizgen/"),
        "data/baysor_vizgen/segmentation.csv",
    threads: 8
    conda:
        "base"
    shell:
        """
        python {input.prepare_script} {input.molecules} vizgen cortex
        mkdir -p data/baysor_vizgen
        bash {input.run_script} data/vizgen_cortex_molecules.csv 6 100 data/baysor_vizgen/
        """

rule baysor_cortex_vizgen_rep2:
    input:
        prepare_script = "scripts/prepare_run_baysor.py",
        run_script = "scripts/run_baysor.sh",
        molecules = "data/vizgen_rep2_molecules.csv",
    output:
        "data/vizgen_rep2_cortex_molecules.csv",
        directory("data/baysor_vizgen_rep2/"),
        "data/baysor_vizgen_rep2/segmentation.csv",
    threads: 8
    conda:
        "base"
    shell:
        """
        python {input.prepare_script} {input.molecules} vizgen_rep2 cortex
        mkdir -p data/baysor_vizgen_rep2
        bash {input.run_script} data/vizgen_rep2_cortex_molecules.csv 6 100 data/baysor_vizgen_rep2/
        """

rule baysor_cortex_vizgen_rep3:
    input:
        prepare_script = "scripts/prepare_run_baysor.py",
        run_script = "scripts/run_baysor.sh",
        molecules = "data/vizgen_rep3_molecules.csv",
    output:
        "data/vizgen_rep3_cortex_molecules.csv",
        directory("data/baysor_vizgen_rep3/"),
        "data/baysor_vizgen_rep3/segmentation.csv",
    threads: 8
    conda:
        "base"
    shell:
        """
        python {input.prepare_script} {input.molecules} vizgen_rep3 cortex
        mkdir -p data/baysor_vizgen_rep3
        bash {input.run_script} data/vizgen_rep3_cortex_molecules.csv 6 100 data/baysor_vizgen_rep3/
        """

rule baysor_cortex_tenx:
    input:
        prepare_script = "scripts/prepare_run_baysor.py",
        run_script = "scripts/run_baysor.sh",
        molecules = "data/tenx_molecules.csv",
    output:
        "data/tenx_cortex_molecules.csv",
        directory("data/baysor_tenx/"),
        "data/baysor_tenx/segmentation.csv",
    threads: 8
    conda:
        "base"
    shell:
        """
        python {input.prepare_script} {input.molecules} tenx cortex
        mkdir -p data/baysor_tenx
        bash {input.run_script} data/tenx_cortex_molecules.csv 6 100 data/baysor_tenx/
        """

rule baysor_cortex_tenx_rep2:
    input:
        prepare_script = "scripts/prepare_run_baysor.py",
        run_script = "scripts/run_baysor.sh",
        molecules = "data/tenx_rep2_molecules.csv",
    output:
        "data/tenx_rep2_cortex_molecules.csv",
        directory("data/baysor_tenx_rep2/"),
        "data/baysor_tenx_rep2/segmentation.csv",
    threads: 8
    conda:
        "base"
    shell:
        """
        python {input.prepare_script} {input.molecules} tenx_rep2 cortex
        mkdir -p data/baysor_tenx_rep2
        bash {input.run_script} data/tenx_rep2_cortex_molecules.csv 6 100 data/baysor_tenx_rep2/
        """

rule baysor_cortex_tenx_rep3:
    input:
        prepare_script = "scripts/prepare_run_baysor.py",
        run_script = "scripts/run_baysor.sh",
        molecules = "data/tenx_rep3_molecules.csv",
    output:
        "data/tenx_rep3_cortex_molecules.csv",
        directory("data/baysor_tenx_rep3/"),
        "data/baysor_tenx_rep3/segmentation.csv",
    threads: 8
    conda:
        "base"
    shell:
        """
        python {input.prepare_script} {input.molecules} tenx_rep3 cortex
        mkdir -p data/baysor_tenx_rep3
        bash {input.run_script} data/tenx_rep3_cortex_molecules.csv 6 100 data/baysor_tenx_rep3/
        """

rule baysor_cortex_resolve:
    input:
        prepare_script = "scripts/prepare_run_baysor.py",
        run_script = "scripts/run_baysor.sh",
        molecules = "data/resolve_molecules.csv",
    output:
        "data/resolve_cortex_molecules.csv",
        directory("data/baysor_resolve/"),
        "data/baysor_resolve/segmentation.csv",
    threads: 8
    conda:
        "base"
    shell:
        """
        python {input.prepare_script} {input.molecules} resolve cortex
        mkdir -p data/baysor_resolve/
        bash {input.run_script} data/resolve_cortex_molecules.csv 30 100 data/baysor_resolve/
        """

rule baysor_cortex_eelfish:
    input:
        prepare_script = "scripts/prepare_run_baysor.py",
        run_script = "scripts/run_baysor.sh",
        molecules = "data/eelfish_molecules.csv",
    output:
        "data/eelfish_cortex_molecules.csv",
        directory("data/baysor_eelfish/"),
        "data/baysor_eelfish/segmentation.csv",
    threads: 8
    conda:
        "base"
    shell:
        """
        python {input.prepare_script} {input.molecules} eelfish cortex
        mkdir -p data/baysor_eelfish/
        bash {input.run_script} data/eelfish_cortex_molecules.csv 60 25 data/baysor_eelfish/
        """

rule baysor_cortex_merfish:
    input:
        prepare_script = "scripts/prepare_run_baysor.py",
        run_script = "scripts/run_baysor.sh",
        molecules = "data/merfish_molecules.csv",
    output:
        "data/merfish_cortex_molecules.csv",
        directory("data/baysor_merfish/"),
        "data/baysor_merfish/segmentation.csv",
    threads: 8
    conda:
        "base"
    shell:
        """
        python {input.prepare_script} {input.molecules} merfish cortex
        mkdir -p data/baysor_merfish/
        bash {input.run_script} data/merfish_cortex_molecules.csv 6 100 data/baysor_merfish/
        """

rule baysor_cortex_merfish_rep2:
    input:
        prepare_script = "scripts/prepare_run_baysor.py",
        run_script = "scripts/run_baysor.sh",
        molecules = "data/merfish_rep2_molecules.csv",
    output:
        "data/merfish_rep2_cortex_molecules.csv",
        directory("data/baysor_merfish_rep2/"),
        "data/baysor_merfish_rep2/segmentation.csv",
    threads: 8
    conda:
        "base"
    shell:
        """
        python {input.prepare_script} {input.molecules} merfish_rep2 cortex
        mkdir -p data/baysor_merfish_rep2/
        bash {input.run_script} data/merfish_rep2_cortex_molecules.csv 6 100 data/baysor_merfish_rep2/
        """

rule baysor_thalamus_vizgen:
    input:
        prepare_script = "scripts/prepare_run_baysor.py",
        run_script = "scripts/run_baysor.sh",
        molecules = "data/vizgen_molecules.csv",
    output:
        "data/vizgen_thalamus_molecules.csv",
        directory("data/baysor_vizgen_thalamus/"),
        "data/baysor_vizgen_thalamus/segmentation.csv",
    threads: 8
    conda:
        "base"
    shell:
        """
        python {input.prepare_script} {input.molecules} vizgen thalamus
        mkdir -p data/baysor_vizgen_thalamus/
        bash {input.run_script} data/vizgen_thalamus_molecules.csv 6 100 data/baysor_vizgen_thalamus/
        """

rule baysor_thalamus_tenx:
    input:
        prepare_script = "scripts/prepare_run_baysor.py",
        run_script = "scripts/run_baysor.sh",
        molecules = "data/tenx_molecules.csv",
    output:
        "data/tenx_thalamus_molecules.csv",
        directory("data/baysor_tenx_thalamus/"),
        "data/baysor_tenx_thalamus/segmentation.csv",
    threads: 8
    conda:
        "base"
    shell:
        """
        python {input.prepare_script} {input.molecules} tenx thalamus
        mkdir -p data/baysor_tenx_thalamus/
        bash {input.run_script} data/tenx_thalamus_molecules.csv 6 100 data/baysor_tenx_thalamus/
        """

rule baysor_cortex_tenx_seurat:
    input:
        script = "scripts/baysor_to_seurat.R",
        tenx_baysor = "data/baysor_tenx/segmentation.csv",
    output:
        "data/tenx_baysor_cortex.rds"
    shell:
        """
        Rscript {input.script} {input.tenx_baysor} tenx cortex
        """

rule baysor_cortex_tenx_rep2_seurat:
    input:
        script = "scripts/baysor_to_seurat.R",
        tenx_baysor = "data/baysor_tenx_rep2/segmentation.csv",
    output:
        "data/tenx_rep2_baysor_cortex.rds"
    shell:
        """
        Rscript {input.script} {input.tenx_baysor} tenx_rep2 cortex
        """

rule baysor_cortex_tenx_rep3_seurat:
    input:
        script = "scripts/baysor_to_seurat.R",
        tenx_baysor = "data/baysor_tenx_rep3/segmentation.csv",
    output:
        "data/tenx_rep3_baysor_cortex.rds"
    shell:
        """
        Rscript {input.script} {input.tenx_baysor} tenx_rep3 cortex
        """

rule baysor_cortex_vizgen_seurat:
    input:
        script = "scripts/baysor_to_seurat.R",
        vizgen_baysor = "data/baysor_vizgen/segmentation.csv",
    output:
        "data/vizgen_baysor_cortex.rds"
    shell:
        """
        Rscript {input.script} {input.vizgen_baysor} vizgen cortex
        """

rule baysor_cortex_vizgen_rep2_seurat:
    input:
        script = "scripts/baysor_to_seurat.R",
        vizgen_baysor = "data/baysor_vizgen_rep2/segmentation.csv",
    output:
        "data/vizgen_rep2_baysor_cortex.rds"
    shell:
        """
        Rscript {input.script} {input.vizgen_baysor} vizgen_rep2 cortex
        """

rule baysor_cortex_vizgen_rep3_seurat:
    input:
        script = "scripts/baysor_to_seurat.R",
        vizgen_baysor = "data/baysor_vizgen_rep3/segmentation.csv",
    output:
        "data/vizgen_rep3_baysor_cortex.rds"
    shell:
        """
        Rscript {input.script} {input.vizgen_baysor} vizgen_rep3 cortex
        """

rule baysor_cortex_merfish_seurat:
    input:
        script = "scripts/baysor_to_seurat.R",
        vizgen_baysor = "data/baysor_merfish/segmentation.csv",
    output:
        "data/merfish_baysor_cortex.rds"
    shell:
        """
        Rscript {input.script} {input.vizgen_baysor} merfish cortex
        """

rule baysor_cortex_merfish_rep2_seurat:
    input:
        script = "scripts/baysor_to_seurat.R",
        vizgen_baysor = "data/baysor_merfish_rep2/segmentation.csv",
    output:
        "data/merfish_rep2_baysor_cortex.rds"
    shell:
        """
        Rscript {input.script} {input.vizgen_baysor} merfish_rep2 cortex
        """

rule baysor_cortex_resolve_seurat:
    input:
        script = "scripts/baysor_to_seurat.R",
        resolve_baysor = "data/baysor_resolve/segmentation.csv",
    output:
        "data/resolve_baysor_cortex.rds"
    shell:
        """
        Rscript {input.script} {input.resolve_baysor} resolve cortex
        """

rule baysor_cortex_eelfish_seurat:
    input:
        script = "scripts/baysor_to_seurat.R",
        eelfish_baysor = "data/baysor_eelfish/segmentation.csv",
    output:
        "data/eelfish_baysor_cortex.rds"
    shell:
        """
        Rscript {input.script} {input.eelfish_baysor} eelfish cortex
        """

rule baysor_thalamus_tenx_seurat:
    input:
        script = "scripts/baysor_to_seurat.R",
        tenx_baysor = "data/baysor_tenx_thalamus/segmentation.csv",
    output:
        "data/tenx_baysor_thalamus.rds"
    shell:
        """
        Rscript {input.script} {input.tenx_baysor} tenx thalamus
        """

rule baysor_thalamus_vizgen_seurat:
    input:
        script = "scripts/baysor_to_seurat.R",
        vizgen_baysor = "data/baysor_vizgen_thalamus/segmentation.csv",
    output:
        "data/vizgen_baysor_thalamus.rds"
    shell:
        """
        Rscript {input.script} {input.vizgen_baysor} vizgen thalamus
        """

rule annotate_seurat_baysor_cortex_tenx:
    input:
        script = "scripts/annotate_seurat.R",
        reference = "data/linnarsson_fullbrain_reference.rds",
        tenx = "data/tenx_baysor_cortex.rds",
    output:
        obj = "data/tenx_baysor_cortex_annotated.rds"
    shell:
        """
        Rscript {input.script} {input.reference} {input.tenx} {output.obj}
        """

rule annotate_seurat_baysor_cortex_tenx_rep2:
    input:
        script = "scripts/annotate_seurat.R",
        reference = "data/linnarsson_fullbrain_reference.rds",
        tenx = "data/tenx_rep2_baysor_cortex.rds",
    output:
        obj = "data/tenx_rep2_baysor_cortex_annotated.rds"
    shell:
        """
        Rscript {input.script} {input.reference} {input.tenx} {output.obj}
        """

rule annotate_seurat_baysor_cortex_tenx_rep3:
    input:
        script = "scripts/annotate_seurat.R",
        reference = "data/linnarsson_fullbrain_reference.rds",
        tenx = "data/tenx_rep3_baysor_cortex.rds",
    output:
        obj = "data/tenx_rep3_baysor_cortex_annotated.rds"
    shell:
        """
        Rscript {input.script} {input.reference} {input.tenx} {output.obj}
        """

rule annotate_seurat_baysor_cortex_vizgen:
    input:
        script = "scripts/annotate_seurat.R",
        reference = "data/linnarsson_fullbrain_reference.rds",
        vizgen = "data/vizgen_baysor_cortex.rds",
    output:
        obj = "data/vizgen_baysor_cortex_annotated.rds"
    shell:
        """
        Rscript {input.script} {input.reference} {input.vizgen} {output.obj}
        """

rule annotate_seurat_baysor_cortex_vizgen_rep2:
    input:
        script = "scripts/annotate_seurat.R",
        reference = "data/linnarsson_fullbrain_reference.rds",
        vizgen = "data/vizgen_rep2_baysor_cortex.rds",
    output:
        obj = "data/vizgen_rep2_baysor_cortex_annotated.rds"
    shell:
        """
        Rscript {input.script} {input.reference} {input.vizgen} {output.obj}
        """

rule annotate_seurat_baysor_cortex_vizgen_rep3:
    input:
        script = "scripts/annotate_seurat.R",
        reference = "data/linnarsson_fullbrain_reference.rds",
        vizgen = "data/vizgen_rep3_baysor_cortex.rds",
    output:
        obj = "data/vizgen_rep3_baysor_cortex_annotated.rds"
    shell:
        """
        Rscript {input.script} {input.reference} {input.vizgen} {output.obj}
        """

rule annotate_seurat_baysor_cortex_merfish:
    input:
        script = "scripts/annotate_seurat.R",
        reference = "data/linnarsson_fullbrain_reference.rds",
        merfish = "data/merfish_baysor_cortex.rds",
    output:
        obj = "data/merfish_baysor_cortex_annotated.rds"
    shell:
        """
        Rscript {input.script} {input.reference} {input.merfish} {output.obj}
        """

rule annotate_seurat_baysor_cortex_merfish_rep2:
    input:
        script = "scripts/annotate_seurat.R",
        reference = "data/linnarsson_fullbrain_reference.rds",
        merfish = "data/merfish_rep2_baysor_cortex.rds",
    output:
        obj = "data/merfish_rep2_baysor_cortex_annotated.rds"
    shell:
        """
        Rscript {input.script} {input.reference} {input.merfish} {output.obj}
        """

rule annotate_seurat_baysor_cortex_resolve:
    input:
        script = "scripts/annotate_seurat.R",
        reference = "data/linnarsson_fullbrain_reference.rds",
        resolve = "data/resolve_baysor_cortex.rds",
    output:
        obj = "data/resolve_baysor_cortex_annotated.rds"
    shell:
        """
        Rscript {input.script} {input.reference} {input.resolve} {output.obj}
        """

rule annotate_seurat_baysor_cortex_eelfish:
    input:
        script = "scripts/annotate_seurat.R",
        reference = "data/linnarsson_fullbrain_reference.rds",
        eelfish = "data/eelfish_baysor_cortex.rds",
    output:
        obj = "data/eelfish_baysor_cortex_annotated.rds"
    shell:
        """
        Rscript {input.script} {input.reference} {input.eelfish} {output.obj}
        """

rule annotate_seurat_baysor_thalamus_tenx:
    input:
        script = "scripts/annotate_seurat.R",
        reference = "data/linnarsson_fullbrain_reference.rds",
        tenx = "data/tenx_baysor_thalamus.rds",
    output:
        obj = "data/tenx_baysor_thalamus_annotated.rds"
    shell:
        """
        Rscript {input.script} {input.reference} {input.tenx} {output.obj}
        """

rule annotate_seurat_baysor_thalamus_vizgen:
    input:
        script = "scripts/annotate_seurat.R",
        reference = "data/linnarsson_fullbrain_reference.rds",
        vizgen = "data/vizgen_baysor_thalamus.rds",
    output:
        obj = "data/vizgen_baysor_thalamus_annotated.rds"
    shell:
        """
        Rscript {input.script} {input.reference} {input.vizgen} {output.obj}
        """

rule export_baysor_annotations_and_markers_tenx:
    input:
        script = "scripts/export_baysor_annotations.R",
        seurat_object = "data/tenx_baysor_cortex_annotated.rds",
    output:
        annotations = "data/tenx_baysor_cortex_annotations.csv",
        markers = "data/tenx_baysor_cortex_markers.csv"
    shell:
        """
        Rscript {input.script} {input.seurat_object} {output.annotations} {output.markers}
        """

rule export_baysor_annotations_and_markers_tenx_rep2:
    input:
        script = "scripts/export_baysor_annotations.R",
        seurat_object = "data/tenx_rep2_baysor_cortex_annotated.rds",
    output:
        annotations = "data/tenx_rep2_baysor_cortex_annotations.csv",
        markers = "data/tenx_rep2_baysor_cortex_markers.csv"
    shell:
        """
        Rscript {input.script} {input.seurat_object} {output.annotations} {output.markers}
        """

rule export_baysor_annotations_and_markers_tenx_rep3:
    input:
        script = "scripts/export_baysor_annotations.R",
        seurat_object = "data/tenx_rep3_baysor_cortex_annotated.rds",
    output:
        annotations = "data/tenx_rep3_baysor_cortex_annotations.csv",
        markers = "data/tenx_rep3_baysor_cortex_markers.csv"
    shell:
        """
        Rscript {input.script} {input.seurat_object} {output.annotations} {output.markers}
        """

rule export_baysor_annotations_and_markers_vizgen:
    input:
        script = "scripts/export_baysor_annotations.R",
        seurat_object = "data/vizgen_baysor_cortex_annotated.rds",
    output:
        annotations = "data/vizgen_baysor_cortex_annotations.csv",
        markers = "data/vizgen_baysor_cortex_markers.csv"
    shell:
        """
        Rscript {input.script} {input.seurat_object} {output.annotations} {output.markers}
        """

rule export_baysor_annotations_and_markers_vizgen_rep2:
    input:
        script = "scripts/export_baysor_annotations.R",
        seurat_object = "data/vizgen_rep2_baysor_cortex_annotated.rds",
    output:
        annotations = "data/vizgen_rep2_baysor_cortex_annotations.csv",
        markers = "data/vizgen_rep2_baysor_cortex_markers.csv"
    shell:
        """
        Rscript {input.script} {input.seurat_object} {output.annotations} {output.markers}
        """

rule export_baysor_annotations_and_markers_vizgen_rep3:
    input:
        script = "scripts/export_baysor_annotations.R",
        seurat_object = "data/vizgen_rep3_baysor_cortex_annotated.rds",
    output:
        annotations = "data/vizgen_rep3_baysor_cortex_annotations.csv",
        markers = "data/vizgen_rep3_baysor_cortex_markers.csv"
    shell:
        """
        Rscript {input.script} {input.seurat_object} {output.annotations} {output.markers}
        """

rule export_baysor_annotations_and_markers_resolve:
    input:
        script = "scripts/export_baysor_annotations.R",
        seurat_object = "data/resolve_baysor_cortex_annotated.rds",
    output:
        annotations = "data/resolve_baysor_cortex_annotations.csv",
        markers = "data/resolve_baysor_cortex_markers.csv"
    shell:
        """
        Rscript {input.script} {input.seurat_object} {output.annotations} {output.markers}
        """

rule export_baysor_annotations_and_markers_eelfish:
    input:
        script = "scripts/export_baysor_annotations.R",
        seurat_object = "data/eelfish_baysor_cortex_annotated.rds",
    output:
        annotations = "data/eelfish_baysor_cortex_annotations.csv",
        markers = "data/eelfish_baysor_cortex_markers.csv"
    shell:
        """
        Rscript {input.script} {input.seurat_object} {output.annotations} {output.markers}
        """

rule export_baysor_annotations_and_markers_merfish:
    input:
        script = "scripts/export_baysor_annotations.R",
        seurat_object = "data/merfish_baysor_cortex_annotated.rds",
    output:
        annotations = "data/merfish_baysor_cortex_annotations.csv",
        markers = "data/merfish_baysor_cortex_markers.csv"
    shell:
        """
        Rscript {input.script} {input.seurat_object} {output.annotations} {output.markers}
        """

rule export_baysor_annotations_and_markers_merfish_rep2:
    input:
        script = "scripts/export_baysor_annotations.R",
        seurat_object = "data/merfish_rep2_baysor_cortex_annotated.rds",
    output:
        annotations = "data/merfish_rep2_baysor_cortex_annotations.csv",
        markers = "data/merfish_rep2_baysor_cortex_markers.csv"
    shell:
        """
        Rscript {input.script} {input.seurat_object} {output.annotations} {output.markers}
        """

rule export_baysor_annotations_and_markers_tenx_thalamus:
    input:
        script = "scripts/export_baysor_annotations.R",
        seurat_object = "data/tenx_baysor_thalamus_annotated.rds",
    output:
        annotations = "data/tenx_baysor_thalamus_annotations.csv",
        markers = "data/tenx_baysor_thalamus_markers.csv"
    shell:
        """
        Rscript {input.script} {input.seurat_object} {output.annotations} {output.markers}
        """

rule export_baysor_annotations_and_markers_vizgen_thalamus:
    input:
        script = "scripts/export_baysor_annotations.R",
        seurat_object = "data/vizgen_baysor_thalamus_annotated.rds",
    output:
        annotations = "data/vizgen_baysor_thalamus_annotations.csv",
        markers = "data/vizgen_baysor_thalamus_markers.csv"
    shell:
        """
        Rscript {input.script} {input.seurat_object} {output.annotations} {output.markers}
        """

rule figure3_exclusive_expression_plots:
    input:
        script = "scripts/figure3_exclusive_expression_plots.R",
        tenx = "data/tenx.rds",
        vizgen = "data/vizgen.rds",
    output:
        p1 = "data/figure3_exclusive_expression_Gad1_Slc17a7_tenx_default_segs_plot.png",
        p2 = "data/figure3_exclusive_expression_Gad1_Slc17a7_vizgen_default_segs_plot.png",
    shell:
        """
        Rscript {input.script} {input.tenx} {input.vizgen} {output.p1} {output.p2}
        """

rule figure3_mecr_vs_frac_mols_in_cells:
    input:
        script = "scripts/figure3_mecr_vs_frac_mols_in_cells.R",
        vizgen = "data/vizgen_rep3_annotated.rds",
        tenx = "data/tenx_rep3_annotated.rds",
        resolve = "data/resolve_annotated.rds",
        eelfish = "data/eelfish_annotated.rds",
        merfish = "data/merfish_annotated.rds",
        vizgen_mols = "data/vizgen_rep3_molecules.csv",
        tenx_mols = "data/tenx_rep3_molecules.csv",
        resolve_mols = "data/resolve_molecules.csv",
        eelfish_mols = "data/eelfish_molecules.csv",
        merfish_mols = "data/merfish_molecules.csv",
        linnarsson = "data/linnarsson_fullbrain_reference.rds",
    output:
        plot = "data/figure3_mecr_vs_frac_mols_in_cells.png",
    shell:
        """
        Rscript {input.script} {input.vizgen} {input.tenx} {input.resolve} {input.eelfish} {input.merfish} \
        {input.vizgen_mols} {input.tenx_mols} {input.resolve_mols} {input.eelfish_mols} {input.merfish_mols} \
        {input.linnarsson} {output.plot}
        """

rule figure3_mecr_counts_plot:
    input:
        script = "scripts/figure3_counts_curve.py",
        markers = "data/wilcox_class_markers_linnarsson.csv",
        mols_1 = "data/baysor_tenx_rep3/segmentation.csv",
        mols_2 = "data/baysor_vizgen_rep3/segmentation.csv",
        mols_3 = "data/baysor_resolve/segmentation.csv",
        mols_4 = "data/baysor_merfish/segmentation.csv",
        mols_5 = "data/baysor_eelfish/segmentation.csv",
        ann_1 = "data/tenx_rep3_baysor_cortex_annotations.csv",
        ann_2 = "data/vizgen_rep3_baysor_cortex_annotations.csv",
        ann_3 = "data/resolve_baysor_cortex_annotations.csv",
        ann_4 = "data/merfish_baysor_cortex_annotations.csv",
        ann_5 = "data/eelfish_baysor_cortex_annotations.csv",
    output:
        all_df = "data/figure3_mecr_plot_data.csv",
        mean_counts_plot = "data/figure3_mecr_mean_counts_plot.png",
        median_counts_plot = "data/figure3_mecr_median_counts_plot.png",
        mean_features_plot = "data/figure3_mecr_mean_features_plot.png",
        median_features_plot = "data/figure3_mecr_median_features_plot.png",
        threshold_mean_counts_plot = "data/figure3_mecr_threshold_mean_counts_plot.png",
    conda:
        "base"
    shell:
        """
        python {input.script} {input.mols_1} {input.mols_2} {input.mols_3} {input.mols_4} {input.mols_5} \
        {input.ann_1} {input.ann_2} {input.ann_3} {input.ann_4} {input.ann_5} \
        {input.markers} 10x Vizgen Resolve MERFISH EELFISH {output.all_df} \
        {output.mean_counts_plot} {output.median_counts_plot} {output.mean_features_plot} \
        {output.median_features_plot} {output.threshold_mean_counts_plot}
        """

rule figure4_dotplot:
    input:
        scripts = "scripts/figure4_dotplot.R",
        macosko = "/brahms/hartmana/spatial_sensitivity_comparison/mackosko_23_mouse_brain_atlas/TH_and_CTX_seurat.rds",
    output:
        violin = "data/figure4_scrna_violin.png",
        dotplot = "data/figure4_scrna_dotplot.png",
    shell:
        """
        Rscript {input.scripts} {input.macosko} {output.violin} {output.dotplot}
        """

rule figure4_volcanoplot:
    input:
        scripts = "scripts/figure4_volcanoplot.R",
        tenx = "data/tenx.rds",
        tenx_ann = "data/tenx_annotated.rds",
    output:
        plot = "data/figure4_tenx_volcanoplot.png",
    shell:
        """
        Rscript {input.scripts} {input.tenx} {input.tenx_ann} {output.plot}
        """

rule export_annotations:
    input:
        script = "scripts/figure4_export_annotations.R",
        obj = "data/tenx_rep2_annotated.rds",
    output:
        "data/tenx_rep2_annotations.csv",
    shell:
        """
        Rscript {input.script} {input.obj} {output}
        """

rule figure4_spatial_molecule_plots:
    input:
        script = "scripts/figure4_spatial_molecule_plots.py",
        annotations = "data/tenx_rep2_annotations.csv",
        mols = "data/tenx_molecules.csv",
        transcripts = "/brahms/hartmana/spatial_sensitivity_comparison/10x_xenium/fresh_frozen_mouse_brain_rep2/transcripts.csv.gz",
        boundaries = "/brahms/hartmana/spatial_sensitivity_comparison/10x_xenium/fresh_frozen_mouse_brain_rep2/cell_boundaries.csv.gz",
    output:
        plot1 = "data/figure4_Slc17a6_Aqp4_molecules.png",
        plot2 = "data/figure4_Satb2_Aqp4_molecules.png",
    conda:
        "segmentation-env"
    shell:
        """
        python {input.script} {input.mols} {input.transcripts} {input.annotations} \
         {input.boundaries} {output.plot1} {output.plot2}
        """

rule figure4_tenx_violin_plots:
    input:
        script = "scripts/figure4_tenx_violin_plots.R",
        tenx = "data/tenx.rds",
        tenx_ann = "data/tenx_annotated.rds",
    output:
        plot = "data/figure4_tenx_violinplot.png",
    shell:
        """
        Rscript {input.script} {input.tenx} {input.tenx_ann} {output.plot}
        """

rule supfig_baysor_violin_plots:
    input:
        script = "scripts/supfig_baysor_violin_plots.R",
        tenx_cortex = "data/tenx_baysor_cortex_annotated.rds",
        tenx_thalamus = "data/tenx_baysor_thalamus_annotated.rds",
        vizgen_cortex = "data/vizgen_baysor_cortex_annotated.rds",
        vizgen_thalamus = "data/vizgen_baysor_thalamus_annotated.rds",
    output:
        tenx_volcano = "data/supfig_tenx_volcanoplot.png",
        tenx_vln = "data/supfig_tenx_violinplot.png",
        vizgen_volcano = "data/supfig_vizgen_volcanoplot.png",
        vizgen_vln = "data/supfig_vizgen_violinplot.png",
    shell:
        """
        Rscript {input.script} {input.tenx_cortex} {input.tenx_thalamus} {input.vizgen_cortex} {input.vizgen_thalamus} \
        {output.tenx_volcano} {output.tenx_vln} {output.vizgen_volcano} {output.vizgen_vln}
        """

rule supfig_default_vizgen_plots:
    input:
        script = "scripts/supfig_default_vizgen_plots.R",
        vizgen = "data/vizgen.rds",
        vizgen_ann = "data/vizgen_annotated.rds",
    output:
        vln = "data/supfig_vizgen_default_violinplot.png",
        volcano = "data/supfig_vizgen_default_volcanoplot.png",
    shell:
        """
        Rscript {input.script} {input.vizgen} {input.vizgen_ann} {output.vln} {output.volcano}
        """

rule supfig_scrna_dotplots:
    input:
        script = "scripts/supfig_scrna_dotplots.R",
        macosko = "/brahms/hartmana/spatial_sensitivity_comparison/mackosko_23_mouse_brain_atlas/TH_and_CTX_seurat.rds",
        vizgen = "data/vizgen.rds",
        vizgen_ann = "data/vizgen_annotated.rds",
        vizgen_cortex = "data/vizgen_baysor_cortex_annotated.rds",
        vizgen_thalamus = "data/vizgen_baysor_thalamus_annotated.rds",
        tenx_cortex = "data/tenx_baysor_cortex_annotated.rds",
        tenx_thalamus = "data/tenx_baysor_thalamus_annotated.rds",
    output:
        vizgen_dotplot = "data/supfig_vizgen_default_dotplot.png",
        vizgen_baysor_dotplot = "data/supfig_vizgen_baysor_dotplot.png",
        tenx_baysor_dotplot = "data/supfig_tenx_baysor_dotplot.png",
    shell:
        """
        Rscript {input.script} {input.macosko} {input.vizgen} {input.vizgen_ann} {input.vizgen_cortex} \
        {input.vizgen_thalamus} {input.tenx_cortex} {input.tenx_thalamus} \
        {output.vizgen_dotplot} {output.vizgen_baysor_dotplot} {output.tenx_baysor_dotplot}
        """               

rule supplemental_mecr_counts_plot:
    input:
        script = "scripts/supplemental_counts_curves_plots.py",
        markers = "data/wilcox_class_markers_linnarsson.csv",
        mols_1 = "data/baysor_tenx_rep2/segmentation.csv",
        mols_2 = "data/baysor_vizgen_rep2/segmentation.csv",
        mols_3 = "data/baysor_resolve/segmentation.csv",
        mols_4 = "data/baysor_merfish/segmentation.csv",
        mols_5 = "data/baysor_eelfish/segmentation.csv",
        mols_6 = "data/baysor_vizgen/segmentation.csv",
        mols_7 = "data/baysor_tenx/segmentation.csv",
        mols_8 = "data/baysor_vizgen_rep3/segmentation.csv",
        mols_9 = "data/baysor_tenx_rep3/segmentation.csv",
        ann_1 = "data/tenx_rep2_baysor_cortex_annotations.csv",
        ann_2 = "data/vizgen_rep2_baysor_cortex_annotations.csv",
        ann_3 = "data/resolve_baysor_cortex_annotations.csv",
        ann_4 = "data/merfish_baysor_cortex_annotations.csv",
        ann_5 = "data/eelfish_baysor_cortex_annotations.csv",
        ann_6 = "data/vizgen_baysor_cortex_annotations.csv",
        ann_7 = "data/tenx_baysor_cortex_annotations.csv",
        ann_8 = "data/vizgen_rep3_baysor_cortex_annotations.csv",
        ann_9 = "data/tenx_rep3_baysor_cortex_annotations.csv",
    output:
        all_df = "data/supp_mecr_plot_data.csv",
        mean_counts_plot = "data/supp_mecr_mean_counts_plot.png",
        median_counts_plot = "data/supp_mecr_median_counts_plot.png",
        mean_features_plot = "data/supp_mecr_mean_features_plot.png",
        median_features_plot = "data/supp_mecr_median_features_plot.png",
        threshold_mean_counts_plot = "data/supp_mecr_threshold_mean_counts_plot.png",
    conda:
        "base"
    shell:
        """
        python {input.script} {input.mols_1} {input.mols_2} {input.mols_3} {input.mols_4} {input.mols_5} \
        {input.mols_6} {input.mols_7} {input.mols_8} {input.mols_9} \
        {input.ann_1} {input.ann_2} {input.ann_3} {input.ann_4} {input.ann_5} \
        {input.ann_6} {input.ann_7} {input.ann_8} {input.ann_9} \
        {input.markers} 10x Vizgen Resolve MERFISH EELFISH Vizgen_rep1 10x_rep1 Vizgen_rep3 10x_rep3 {output.all_df} \
        {output.mean_counts_plot} {output.median_counts_plot} {output.mean_features_plot} \
        {output.median_features_plot} {output.threshold_mean_counts_plot}
        """

rule supplemental_mecr_counts_plot_thalamus:
    input:
        script = "scripts/supplemental_counts_curves_plots.py",
        markers = "data/wilcox_class_markers_linnarsson.csv",
        mols_1 = "data/baysor_tenx_thalamus/segmentation.csv",
        mols_2 = "data/baysor_vizgen_thalamus/segmentation.csv",
        mols_3 = "data/baysor_resolve/segmentation.csv",
        mols_4 = "data/baysor_merfish/segmentation.csv",
        mols_5 = "data/baysor_eelfish/segmentation.csv",
        mols_6 = "data/baysor_vizgen/segmentation.csv",
        mols_7 = "data/baysor_tenx/segmentation.csv",
        mols_8 = "data/baysor_vizgen_rep3/segmentation.csv",
        mols_9 = "data/baysor_tenx_rep3/segmentation.csv",
        ann_1 = "data/tenx_baysor_thalamus_annotations.csv",
        ann_2 = "data/vizgen_baysor_thalamus_annotations.csv",
        ann_3 = "data/resolve_baysor_cortex_annotations.csv",
        ann_4 = "data/merfish_baysor_cortex_annotations.csv",
        ann_5 = "data/eelfish_baysor_cortex_annotations.csv",
        ann_6 = "data/vizgen_baysor_cortex_annotations.csv",
        ann_7 = "data/tenx_baysor_cortex_annotations.csv",
        ann_8 = "data/vizgen_rep3_baysor_cortex_annotations.csv",
        ann_9 = "data/tenx_rep3_baysor_cortex_annotations.csv",
    output:
        all_df = "data/supp_mecr_plot_thalamus_data.csv",
        mean_counts_plot = "data/supp_mecr_mean_counts_thalamus_plot.png",
        median_counts_plot = "data/supp_mecr_median_counts_thalamus_plot.png",
        mean_features_plot = "data/supp_mecr_mean_features_thalamus_plot.png",
        median_features_plot = "data/supp_mecr_median_features_thalamus_plot.png",
        threshold_mean_counts_plot = "data/supp_mecr_threshold_mean_counts_thalamus_plot.png",
    conda:
        "base"
    shell:
        """
        python {input.script} {input.mols_1} {input.mols_2} {input.mols_3} {input.mols_4} {input.mols_5} \
        {input.mols_6} {input.mols_7} {input.mols_8} {input.mols_9} \
        {input.ann_1} {input.ann_2} {input.ann_3} {input.ann_4} {input.ann_5} \
        {input.ann_6} {input.ann_7} {input.ann_8} {input.ann_9} \
        {input.markers} 10x Vizgen Resolve MERFISH EELFISH Vizgen_rep1 10x_rep1 Vizgen_rep3 10x_rep3 {output.all_df} \
        {output.mean_counts_plot} {output.median_counts_plot} {output.mean_features_plot} \
        {output.median_features_plot} {output.threshold_mean_counts_plot}
        """
