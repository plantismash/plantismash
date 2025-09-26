var gene_colors = [
  { label: "Cytochrome 450", color: "#D50000", members : ["plants/p450"] },
  { label: "Terpene synthase", color: "#35F800", members : ["plants/Terpene_synth", "plants/Terpene_synth_C", "plants/SQHop_cyclase_C", "plants/SQHop_cyclase_N", "plants/Lycopene_cycl"] },
  { label: "Copper amine oxidase", color: "#E9C63E", members : ["plants/Cu_amine_oxid"] },
  { label: "Pictet-Spengler enzyme (Bet v1)", color: "#735100", members : ["plants/Bet_v_1"] },
  { label: "Glycosyltransferase", color: "#F477A6", members : ["plants/Glycos_transf_1", "plants/Glycos_transf_2", "plants/Glyco_transf_28", "plants/UDPGT", "plants/UDPGT_2", "plants/A__Glycosidic_branch_elongating", "plants/G__Monoterpenoid_cyanogenic_glucoside"] },
  { label: "Ketosynthase", color: "#75FFD8", members : ["plants/Chal_sti_synt_C", "plants/Chal_sti_synt_N"] },
  { label: "Squalene epoxidase", color: "#206B14", members : ["plants/SE"] },
  { label: "Carboxylesterase", color: "#A221FF", members : ["plants/COesterase"] },
  { label: "Methyltransferase", color: "#CF0FFF", members : ["plants/Methyltransf_2", "plants/Methyltransf_3", "plants/Methyltransf_7", "plants/Methyltransf_11", "plants/cMT", "plants/nMT", "plants/oMT"] },
  { label: "BAHD acyltransferase", color: "#0003EC", members : ["plants/Transferase"] },
  { label: "Scl acyltransferase", color: "#2935E2", members : ["plants/Peptidase_S10"] },
  { label: "Dioxygenase", color: "#F5B2FF", members : ["plants/DIOX_N", "plants/2OG-FeII_Oxy"] },
  { label: "Amino oxidase", color: "#2C113A", members : ["plants/Amino_oxidase"] },
  { label: "Aminotransferase", color: "#E42E00", members : ["plants/Aminotran_1_2", "plants/Aminotran_3"] },
  { label: "Prenyltransferase", color: "#552288", members : ["plants/Prenyltrans", "plants/Prenyltransf", "plants/UbiA"] },
  { label: "Strictosidine synthase-like", color: "#dd71dd", members : ["plants/Str_synth"] },
  { label: "PRISE enzymes", color: "#13C600", members : ["plants/PRISE"] },
  { label: "Dirigent enzymes", color: "#606000", members : ["plants/Dirigent"] },
  { label: "Cellulose synthase-like", color: "#009115", members :  ["plants/Cellulose_synt"]},
  { label: "Lipoxygenase", color: "#071daf", members : ["plants/Lipoxygenase"] },
  { label: "Aldo-keto reductase", color: "#8299f6", members : ["plants/Aldo_ket_red"] },
  { label: "Polyprenyl synthetase", color: "#ffa83f", members : ["plants/polyprenyl_synt"] },
  { label: "Short-chain dehydrogenase/reductase", color: "#F245BD", members : ["plants/Epimerase", "plants/adh_short", "plants/adh_short_C2",'plants/3Beta_HSD'] },
  { label: "Oxidoreductase", color: "#40D885", members : ["plants/BBE", "plants/Oxidored_FMN", "plants/NAD_binding_1", "plants/GMC_oxred_N", "plants/GMC_oxred_C"] },
  { label: "Alcohol dehydrogenase", color: "#40D885", members : ["plants/ADH_N", "plants/ADH_N_2", "plants/ADH_zinc_N"] },
  { label: "CoA-ligase", color: "#8F400B", members : ["plants/AMP-binding","plants/ECH_2"] },
  { label: "Fatty acid desaturase", color: "#ef900b", members : ["plants/FA_desaturase", "plants/FA_desaturase_2"]},
  { label: "Fatty acid hydroxylase", color: "#ef900b", members : ["plants/FA_hydroxylase", "plants/CER1-like_C"]},
  { label: "pyridoxal synthase", color: "#a8eeab", members : ["plants/YjeF_N", "plants/Pyridox_oxidase", "plants/PNPOx_C"] },
  { label: "pyridoxal-dependent Amino acid decarboxylase", color: "#FFC300", members : ["plants/Orn_Arg_deC_N", "plants/Orn_DAP_Arg_deC", "plants/Pyridoxal_deC"] },
  { label: "Amino acid Dehydrogenase", color: "#FFC300", members : ["plants/PALP", "plants/Thr_dehydrat_C", "plants/E1_dh"] },
  { label: "Transporter", color: "#0abaef", members : ["plants/MatE", "plants/LTP_2", "plants/ABC2_membrane", "plants/ABC_tran"] },
  { label: "Hydrolase", color: "#1040ef", members : ["plants/Abhydrolase_3"] },
  { label: "Glycosyl hydrolase", color: "#1040ef", members : ["plants/Glyco_hydro_1"] },
  { label: "Peptide cyclases", color: "#05f5ce", members : ["plants/BURP"] }, 
  { label: "Male sterility protein", color: "#FF6E54", members : ["plants/NAD_binding_4"] },
  { label: "FAE1/Type III polyketide synthase-like protein", color: "#7C29F0", members : ["plants/FAE1_CUT1_RppA"] },
  { label: "Glycerol-3-phosphate acyltransferase RAM2-like, HAD-like domain", color: "#2EB67D", members : ["plants/HAD_RAM2_N"] }, 
  { label: "Chalcone Synthase-like", color: "#FFD700", members : ["plants/Chalcone", "plants/Chalcone_2", "plants/Chalcone_3"] },
  { label: "Tryptophan Synthase", color: "#C71585", members : ["plants/Trp_syntA"] },
  { label: "Histidine Biosynthesis", color: "#FF69B4", members : ["plants/His_biosynth"] },
  { label: "Shikimate Pathway Enzyme", color: "#BC8F8F", members : ["plants/DAHP_synth_1", "plants/DAHP_synth_2"] },
  { label: "Aromatic Ring Lyase", color: "#FF8C00", members : ["plants/Lyase_aromatic"] },
  { label: "Sterol Biosynthesis", color: "#8B4513", members : ["plants/ERG4_ERG24"] },
  { label: "Thiopurine Methyltransferase", color: "#9932CC", members : ["plants/TPMT"] },
  { label: "Hydroxymethylglutaryl-CoA Lyase", color: "#A52A2A", members : ["plants/HMGL-like"] },
  { label: "Acetyltransferase", color: "#A0522D", members : ["plants/Acetyltransf_1"] },
  { label: "Protein Dimerisation Domain", color: "#4682B4", members : ["plants/Dimerisation"] },
  { label: "Putative Pyridoxamine Oxidase", color: "#7FFFD4", members : ["plants/Putative_PNPOx"] },
  { label: "Pyridoxine Oxidase C-terminal", color: "#00CED1", members : ["plants/PNP_phzG_C"] }
];

function get_gene_color(orf) {
	var color = "gray";
	var newColor = false;
	if (orf.hasOwnProperty("domains") && orf["domains"].length > 0) {
		var domain = orf["domains"][0];
		for (var i = 0; i < gene_colors.length; i++) {
			var elm = gene_colors[i];
			if (elm.members.indexOf(domain) >= 0) {
				color = elm.color;
				newColor = true;
				break;
			}
		}
	}
	if (!newColor) {
		if (orf["type"] == "biosynthetic") {
			color = "#810e15";
		}
	}
	return color;
}

function get_legend_obj(orf) {
	var legobj = null;
	if (orf.hasOwnProperty("domains") && orf["domains"].length > 0) {
		var domain = orf["domains"][0];
		for (var i = 0; i < gene_colors.length; i++) {
			var elm = gene_colors[i];
			if (elm.members.indexOf(domain) >= 0) {
				legobj = elm;
				break;
			}
		}
	}
	return legobj;
}