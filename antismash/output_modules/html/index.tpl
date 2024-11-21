<!doctype html>
<html>
  <head>
    <title>antiSMASH results</title>
    <link rel="stylesheet" type="text/css" href="css/vis.min.css">
    <link rel="stylesheet" type="text/css" href="css/style.css">
    <link rel="stylesheet" type="text/css" href="css/datatable.css">
    <link rel="stylesheet" type="text/css" href="css/coexpress.d3.hiveplot.css">
    <link rel="stylesheet" type="text/css" href="css/modal.css">
    <meta charset="utf-8" />
  </head>
  <body>
    <div id="header">
      <div class="top-header">
        <img class="antismash-logo" src="images/plantismash.png" alt="plantiSMASH">
        <span class="antismash-title"><a href="http://plantismash.secondarymetabolites.org">Plant Secondary Metabolite Analysis</a><br>
            <span class="white">Version 1.0.0-beta<span id="antismash-version" style="display: none;"></span></span>
        </span>
        <div id="icons">
          <a href="http://plantismash.secondarymetabolites.org/"><img src="images/home.png" alt="home" title="Go to start page"></a>
          <a href="http://plantismash.secondarymetabolites.org/help.html"><img src="images/help.png" alt="help" title="Get help using plantiSMASH"></a>
          <a href="http://plantismash.secondarymetabolites.org/about.html"><img src="images/about.png" alt="about" title="About plantiSMASH"></a>
          <a href="#" id="download"><img src="images/download.png" alt="download" title="Download results"></a>
          <div id="downloadmenu">
            <ul id="downloadoptions">
            </ul>
          </div>
        </div>
      </div>
      <div id="buttons">
        <span id="cluster-type">Select Gene Cluster:</span>
        <ul id="clusterbuttons">
          <li class="clbutton"><a href="#">Overview</a></li>
        </ul>
      </div>
    </div>

    <!-- overview page -->
    <div class="page" id="overview">
      <h3>Identified secondary metabolite clusters<span id="truncated"></span></h3>
      <table id="cluster-overview">
        <thead>
          <tr>
            <th>Cluster</th>
            <th>Type</th>
            <th>From</th>
            <th>To</th>
            <th>Size (kb)</th>
            <th>Core domains</th>
            <th>Product/substrate predicted by subgroup</th>
            <th>Most similar known cluster</th>
            <th>MIBiG BGC-ID</th>
          </tr>
        </thead>
        <tbody>
        </tbody>
      </table>
    </div>

    <div id="footer">
      <div id="logos">
      	<table id="logo-table">
      		<tr>
      			<td>
        			<img src='images/wur-logo.png' />
        		</td>
        		<td>
        			<img src='images/unila-logo.png' />
        		</td>
        		<td>
        			<img src='images/cfb-logo.png' />
        		</td>
        		<td>
				      <img src='images/jic-logo.png' />
        		</td>
        	</tr>
        </table>
      </div>
      <div id="copyright">
        If you have found plantiSMASH useful, please <a href="http://plantismash.secondarymetabolites.org/about">cite us</a>.
      </div>
    </div>

    <script src="js/jquery.js"></script>
    <script src="js/purl.js"></script>
    <script src="js/d3.v2.js"></script>
    <script src="js/svgene.js"></script>
    <script src="js/inchlib-1.2.0.1-satria.js"></script>
    <script src="js/vis.min.js"></script>
    <script src="js/kinetic-v5.1.0.min.js"></script>
    <script src="js/plotly-latest.min.js"></script>
    <script src="js/jsdomain.js"></script>
    <script src="js/clusterblast.js"></script>
    <script src="js/datatable.js"></script>
    <script src="js/svg-pan-zoom.min.js"></script>
    <script src="js/coexpress.js"></script>
    <script src="js/coexpress.d3.hiveplot.js"></script>
    <script src="geneclusters.js"></script>
    <script src="js/gene_colors.js"></script>
    <script src="js/modal.js"></script>
    <script type="text/javascript">
function toggle_downloadmenu(event) {
    event.preventDefault();
    $("#downloadmenu").fadeToggle("fast", "linear");
}

function switch_to_cluster() {
    setTimeout(function() {
        var url = $.url();
        $(".page").hide();
        var anchor = url.data.attr.fragment;
        if (anchor == "") {
            anchor = "overview";
        }
        $("#" + anchor).show();
        if (geneclusters[anchor] !== undefined) {
            svgene.drawClusters(anchor+"-svg", [geneclusters[anchor]], 20, 700);
        }
        if ($("#" + anchor + "-details-svg").length > 0) {
            jsdomain.drawDomains(anchor+ "-details-svg", details_data[anchor], 40, 700);
        }
        $("#" + anchor + " .clusterblast-selector").change();
        // --- add zoom & scroll ---
        if (anchor !== "overview") {
            clusterSvgHeight =  parseInt($("#" + anchor + "-svg>svg").attr("height"));
            $("#" + anchor + "-svg>svg").attr("height", clusterSvgHeight + 42);
            $(".svgene-tooltip").css("margin-top", "2em");
            panZoom = svgPanZoom("#" + anchor + "-svg>svg", {
                panEnabled: true,
                zoomEnabled: true,
                dblClickZoomEnabled: false,
                maxZoom: 1.77,
                minZoom: 1,
                beforePan: function(oldPan, newPan){
                    var sizes = this.getSizes()
                    , leftLimit = -((sizes.viewBox.x + sizes.viewBox.width) * sizes.realZoom) + sizes.viewBox.width
                    , rightLimit = sizes.width - sizes.width - (sizes.viewBox.x * sizes.realZoom);

                    customPan = {};
                    customPan.x = Math.max(leftLimit, Math.min(rightLimit, newPan.x));
                    customPan.y = false;

                    return customPan
                }
            });
            $("#pz-btn-container").remove();
            $("#" + anchor + "-svg").append("<div id='pz-btn-container'></div>");
            $("#pz-btn-container").append("<button onClick='javascript:panZoom.panBy({x: 10, y: 0});'><</button>");
            $("#pz-btn-container").append("<button onClick='javascript:panZoom.zoomOut();'>-</button>");
            $("#pz-btn-container").append("<button onClick='javascript:panZoom.reset();'>reset</button>");
            $("#pz-btn-container").append("<button onClick='javascript:panZoom.zoomIn();'>+</button>");
            $("#pz-btn-container").append("<button onClick='javascript:panZoom.panBy({x: -10, y: 0});'>></button>");
        }
        // --- end of add zoom & scroll ---
        // --- apply colors of the genes ---
    		update_gene_colors(anchor);
    		// --- end of apply colors of the genes ---
        // -- coexpress show/hide options --
        $("#" + anchor + "-coexpress-rec_id").val($("#coexpress-rec_id").val());
        $("#" + anchor + "-coexpress-show_ego").prop("checked", $("#coexpress-show_ego").val() == "true");
        if ($("#coexpress-show_heatmap_fluctuation").val() == "true") {
          $("#" + anchor + "-coexpress-show_heatmap_fluctuation").prop("checked", true);
        } else {
          $("#" + anchor + "-coexpress-show_heatmap_intensity").prop("checked", true);
        }
        $("#" + anchor + "-coexpress-network_cutoff").val($("#coexpress-network_cutoff").val());
        // -- end of coexpress show/hide options --
        // --- update coexpress summary on first time page opened ---
        if ($("#" + anchor + "-coexpress-summary").length < 1) {
            draw_coexpress_summary(anchor, $("#" + anchor + "-coexpress-rec_id").val());
        }
        // --- end of update coexpress summary on first time page opened ---
    }, 1);
}

function update_gene_colors(anchor) {
	var gene_domains = {};
	$("#" + anchor + " .content .legend .legend-container:not(.lcSkip)").addClass("hidden");
	$("#" + anchor + "-svg div.svgene-tooltip").each(function(index) {
		var match = $(this).attr("id").match(/(\d+)-cluster(\d+)-(.+)-tooltip/);
		if (match && match.length == 4) {
			var orf = match[3];
			var locus_tag = "";
			$("#" + anchor + "-svg svg g text.svgene-locustag").each(function(index1) {
				var pat = new RegExp("(\\d+)-cluster(\\d+)-" + orf + "-label");
				if (pat.test($(this).attr("id"))) {
					locus_tag = $(this).text();
				}
			});
			$("#" + anchor + "-svg svg g polygon.svgene-orf").each(function(index2) {
				var pat = new RegExp("(\\d+)-cluster(\\d+)-" + orf + "-orf");
				if (pat.test($(this).attr("id"))) {
					var domain = "";
					for (var i = 0; i < geneclusters[anchor]["orfs"].length; i++) {
						if (geneclusters[anchor]["orfs"][i]["locus_tag"] == locus_tag) {
							if (geneclusters[anchor]["orfs"][i].hasOwnProperty("domains")) {
								domain = geneclusters[anchor]["orfs"][i]["domains"][0];
							}
							break;
						}
					}
					for (var i = 0; i < gene_colors.length; i++) {
						var elm = gene_colors[i];
						if (elm.members.indexOf(domain) >= 0) {
							$(this).css("fill", elm.color);
							$("#" + anchor + " .content .legend .legend-container.lc" + i).removeClass("hidden");
							break;
						}
					}
				}
			});
		}
	});
}

function update_legends() {
	$(".page .content .legend").each(function(index){
		$(this).html("");
		$(this).append("<h4>Legend:</h4>");
		var gcs = [];
		gene_colors.forEach(function(elm) { gcs.push(elm); });
		gcs.push({ label: "(Other) Biosynthetic Genes", color: "#810e15", members : [] });
		gcs.push({ label: "Other Genes", color: "gray", members : [] });
		for (var i = 0; i < gcs.length; i++) {
			var div = $("<div class='legend-container' style='width: 20em; float: left; overflow: hidden; margin-bottom: 1em;'/>");
			div.append("<div class='legend-field' style='float: left; border: 2px solid " + gcs[i].color + "; background-color: " + gcs[i].color + ";'/>");
			div.append("<div class='legend-label' style='float: left;'>" + gcs[i].label + "</div>");
			if (gene_colors.indexOf(gcs[i]) >= 0) {
				div.addClass("lc" + gene_colors.indexOf(gcs[i]));
			} else {
				div.addClass("lcSkip");
			}
			$(this).append(div);
		}
		$(this).append("<hr style='width: 100%; clear: both;' />");
		var borders = [
			["biosynthetic", "biosynthetic genes"],
			["transport", "transport-related genes"],
			["regulatory", "regulatory genes"],
			["other", "other genes"]
		];
		for (var i = 0; i < borders.length; i++) {
			var div = $("<div style='width: 20em; float: left; overflow: hidden; margin-bottom: 1em;'/>");
			div.append("<div class='legend-field legend-type-" + borders[i][0] + "' style='float: left;'/>");
			div.append("<div class='legend-label' style='float: left;'>" + borders[i][1] + "</div>");
			$(this).append(div);
		}
	});
}

function toggle_cluster_rules(ev) {
    ev.preventDefault();
    var id = $(this).attr('id').replace(/-header/, '');
    var rules = $('#' + id);
    if (rules.css('display') == "none") {
        $(this).text('Hide pHMM detection rules used');
    } else {
        $(this).text('Show pHMM detection rules used');
    }
    rules.fadeToggle("fast", "linear");
}

function map_type_to_desc(type) {
    switch(type) {
      case "nrps": return "NRPS";
      case "t1pks": return "Type I PKS";
      case "t2pks": return "Type II PKS";
      case "t3pks": return "Type III PKS";
      case "t4pks": return "Type IV PKS";
      default: return type;
    }
}

function copyToClipboard (text) {
    window.prompt ("Copy to clipboard: Ctrl+C, Enter", text);
}

function updateShowHideGenes(cluster_id) {
  var cluster = geneclusters[cluster_id];
  if (cluster !== "undefined") {
    var div = $("#" + cluster_id + ">.content>.showhide");
    var core_checked = div.find("div>input#showhide-core-" + cluster_id).is(":checked");
    var other_checked = div.find("div>input#showhide-other-" + cluster_id).is(":checked");
    $("#" + cluster_id + "-svg polygon.svgene-orf").css("display", (other_checked? "":"none"));
    $("#" + cluster_id + "-svg polygon.svgene-type-biosynthetic").css("display", (core_checked? "":"none"));
  }
}

function addShowHideGenes() {
  for (var cluster_id in geneclusters) {
    // --- add show / hide ---
    var div = $("<div class='showhide'>");
    div.append("<h4>Show:</h4>");
    div.append("<div></div>");
    div.find("div").append("<input type='checkbox' id='showhide-core-" + cluster_id + "' class='showhide-field' checked='checked' onChange='javascript:updateShowHideGenes(\"" + cluster_id + "\");'>biosynthetic genes</input>");
    div.find("div").append("<input type='checkbox' id='showhide-other-" + cluster_id + "' class='showhide-field' checked='checked' onChange='javascript:updateShowHideGenes(\"" + cluster_id + "\");'>other genes</input>");
    $("#" + cluster_id + ">.content>.legend").before(div);
    // --- end of add show / hide ---
  }
}

function applyHybridColoring() {
    $(".clbutton.hybrid").each(function(idx, elmt) {
        var hybrid_class = $(elmt).attr("class").match(/clbutton (.*) hybrid/);
        if (hybrid_class.length > 1) {
            if (true) { // this should be replaced with checking whether a specific style for the hybrid cluster not exist
                var bgcolor = "";
                var gdcolors = [];
                var classes = hybrid_class[1].split("-");
                for (var i = 0; i < classes.length; i++) {
                    if ($("." + classes[i]).length < 1) { $("body").append("<div style='display: none;' class='" + classes[i] + "'>"); }
                    var clcolor = $("." + classes[i]).css("background-color");
                    if (clcolor != undefined) {
                        var rgb = clcolor.match(/\d+/g);
                        gdcolors.push(Math.round((parseInt(rgb[0]) + parseInt(rgb[1]) + parseInt(rgb[2])) / 3));
                        if (bgcolor.length > 0) {
                            bgcolor += ", ";
                        }
                        bgcolor += clcolor;
                    }
                }
                var gdsum = 0;
                var ftcolor = "white";
                for (var i = 0; i < gdcolors.length; i++) {
                    gdsum += gdcolors[i];
                }
                if ((gdsum / gdcolors.length) > 150) {
                    ftcolor = "black";
                }
                $("." + hybrid_class[1] + " a").css("color", ftcolor);
                $("." + hybrid_class[1]).css("background", "white");
                $("." + hybrid_class[1]).css("background", "-webkit-linear-gradient(" + bgcolor + ")");
                $("." + hybrid_class[1]).css("background", "-o-linear-gradient(" + bgcolor + ")");
                $("." + hybrid_class[1]).css("background", "-moz-linear-gradient(" + bgcolor + ")");
                $("." + hybrid_class[1]).css("background", "linear-gradient(" + bgcolor + ")");
            }
        }
    })
}

function addCoExpressSummary() {
    for (var clusterid in geneclusters) {
        if (geneclusters.hasOwnProperty(clusterid)) {
            var coexpress = geneclusters[clusterid]["geo"];
            if (coexpress == undefined || coexpress.length < 1) { continue; }
            var gbdiv = $("<div class='coexpress'>");
            var options = "";
            for (var i = 0; i < geo_dataset_info.length; i++) {
                var geo_data = geo_dataset_info[i];
                options += "<option value='" + geo_data["id"] + "'>" + geo_data["id"] + " : " + geo_data["title"] + "</value>";
            }
            gbdiv.append("<h3>Coexpression Analysis</h3>");
            gbdiv.append("<h4>GEO Record : <select id='" + clusterid + "-coexpress-rec_id' onchange='javascript:draw_coexpress_summary(\"" + clusterid + "\", this.value);'>" + options + "</select></h4>");
            gbdiv.append("<div id='" + clusterid + "-coexpress-heatmap_header' style='margin-left: 1em;'><div><h4>Expression Heatmap</h4>Show : <input id='" + clusterid + "-coexpress-show_heatmap_fluctuation' name='" + clusterid + "-coexpress-show_heatmap' type='radio' checked='checked' onchange='javascript:draw_coexpress_summary(\"" + clusterid + "\", $(\"#" + clusterid + "-coexpress-rec_id\").val());'/> Expression fluctuation <input id='" + clusterid + "-coexpress-show_heatmap_intensity' name='" + clusterid + "-coexpress-show_heatmap' type='radio' onchange='javascript:draw_coexpress_summary(\"" + clusterid + "\", $(\"#" + clusterid + "-coexpress-rec_id\").val());'/> Expression intensity</div></div>");
            gbdiv.append("<div id='" + clusterid + "-coexpress-network_header' style='margin-left: 1em;'><div><h4>Correlation Network</h4>Show : <input id='" + clusterid + "-coexpress-show_ego' type='checkbox' onchange='javascript:draw_coexpress_summary(\"" + clusterid + "\", $(\"#" + clusterid + "-coexpress-rec_id\").val());'/> Ego Networks" + "<div>Network distance cutoff (0-200): <input id='" + clusterid + "-coexpress-network_cutoff' type='text' /> <button onclick='javascript:draw_correlation_network(\"" + clusterid + "\", $(\"#" + clusterid + "-coexpress-rec_id\").val(), $(\"#" + clusterid + "-coexpress-network_cutoff\").val());'>Update</button>" + "</div>");
            $("#" + clusterid + ">.content").append(gbdiv);
        }
    }
    $("body").append("<input type='hidden' id='coexpress-rec_id'>");
    $("body").append("<input type='hidden' id='coexpress-show_ego' value='true'>");
    $("body").append("<input type='hidden' id='coexpress-show_heatmap_fluctuation' value='true'>");
    $("body").append("<input type='hidden' id='coexpress-network_cutoff' value='50'>");
    if (typeof inter_cluster_data !== 'undefined') {
        // Draw hiveplot overview
        $("#overview").append("<br /><h3>Cluster expression relation overview<span id='truncated'></span></h3>");
        var options = "";
        for (var i = 0; i < geo_dataset_info.length; i++) {
            var geo_data = geo_dataset_info[i];
            options += "<option value='" + geo_data["id"] + "'>" + geo_data["id"] + " : " + geo_data["title"] + "</value>";
        }
        $("#overview").append("<h4>GEO Record : <select id='cluster_exp_hiveplot-selected' onchange='javascript:drawHivePlot(this.value); drawSignalPlot(this.value);'>" + options + "</select></h4>");
        $("#overview").append("<div id='cluster_exp_overview'><div id='cluster_exp_hiveplot'></div><div id='cluster_exp_hiveplot_legend'></div><div id='cluster_exp_signalplot'></div>");
        drawHivePlot($("#cluster_exp_hiveplot-selected").val());
    }
}

function addCdhitSummary() {
    var sibling_col_idx = 0;
    for (var clusterid in geneclusters) {
        if (geneclusters.hasOwnProperty(clusterid)) {
            var cluster = geneclusters[clusterid];
                if ($("#cluster-overview>thead th:contains('CD-HIT')").length === 0) {
                sibling_col_idx = $("#cluster-overview>thead th:contains('Most similar known cluster')").index();
                $("#cluster-overview>thead th:contains('Most similar known cluster')").before("<th>CD-HIT</th");
            }
            if (sibling_col_idx > 0) {
                if (!("cdhitclusters" in cluster)) {
                   cluster["cdhitclusters"] = [];
                }
                $("#cluster-overview>tbody>tr>td>a[href='#" + clusterid + "']").parents("tr").find("td").eq(sibling_col_idx).before("<td>" + cluster["cdhitclusters"].length + "</td>");
            }
        }
    }
}

function addTableSummary() {
  for (var anchor in geneclusters) {
    if (geneclusters.hasOwnProperty(anchor)) {
      var contStyle = "font-size: 80%; overflow: scroll; padding-left: 0.5em; margin-bottom: 0.5em;";
      var theTable = "<table style='margin-top: 1em; border: 1px solid black;'>";
      theTable += "<thead><tr style='border: 1px solid black;'><th style='width: 13em;'>Locus tag</th><th style='width: 15em;'>Functional annotation</th><th style='width: 5.5em;'>From</th><th style='width: 5.5em;'>To</th><th style='width: 5em;'>Strand</th><th style='width: 10em;'>Category</th><th style='width: 15em;'>Domains</th><th style='width: 15em;'>Subgroup</th></tr></thead><tbody>";
      for (var i in geneclusters[anchor]["orfs"]) {
        var orf = geneclusters[anchor]["orfs"][i];
        var locus_tag = orf["locus_tag"];
        var from = orf["start"];
        var to = orf["end"];
        var strand = "+";
        if (orf["strand"] < 0) {
          strand = "-";
        }
        var annot_name = "-";
        var subgroup = orf["subgroup"];
        var domains = "";
        var category = get_legend_obj(orf);
        if (category === null) {
          if (orf["type"] == "biosynthetic") {
            category = "(Other) biosynthetic genes";
          } else {
            category = "Other genes";
          }
        } else {
          category = category["label"];
        }
        var desc = orf["description"];
		    var match = desc.match(/<span class=\"svgene-tooltip-bold\">(.*)<\/span>/);
    		if (match) {
          if (match[1].length > 0) {
            annot_name = match[1];
          }
        }
        if (orf["domains"].length > 0) {
          var split = desc.split("<br>");
          for (var x in split) {
            match = split[x].match(/(.+) \(E-value: (.+), bitscore: (.+), seeds: (.+)\)/);
            if (match) {
              if (domains.length > 0) {
                domains += ", ";
              }
              domains += match[1] + " (E=" + match[2] + ")";
            }
          }
        }
        if (domains === "") {
          domains = orf['domain_present']
        }
        theTable += "<tr id='cl-summary-" + anchor +"-tab-" + locus_tag.replace(/(:|\.)/g, '-') + "' onmouseover='javascript: highlightGene(\"" + anchor + "\", \"" + locus_tag + "\");' onmouseout='javascript: deHighlightGene(\"" + anchor + "\", \"" + locus_tag + "\");'><td>" + locus_tag + "</td><td>" + annot_name + "</td><td style='text-align: right;'>" + from + "</td><td style='text-align: right;'>" + to + "</td><td style='text-align: center;'>" + strand + "</td><td>" + category + "</td><td>" + domains + "</td><td>" + subgroup + "</td></tr>"
      }
      theTable += "</tbody></table>";
      $("#" + anchor + ".page .description-container").after("<div id='cl-summary-" + anchor +"' style='" + contStyle + "'><div><b>Genes:</b> <button class='showhidebutton' onclick='javascript: showHideTableSummary(\"" + anchor + "\");'>hide</button></div><div class='cl-summary-tab'>" + theTable + "</div></div>");

      $("#cl-summary-" + anchor + " .cl-summary-tab tr td:last-child").on("click", function() {
          var locusTag = $(this).closest("tr").find("td:first").text();
          var subgroup = $(this).text(); // Get the subgroup text
          openSvgInWindow(locusTag, subgroup); // Pass subgroup as an additional parameter
      });
    }
  }
}

var svgWindow = null;

function openSvgInWindow(locusTag, subgroup) {
    var svgUrl = "subgroup/tree_svg/" + locusTag + ".svg";
    var legendUrl = "subgroup/tree_svg/" + locusTag + "_legend.svg";
    var windowTitle = locusTag + " - subgroup: " + subgroup;

    var content = `
        <html>
        <head><title>${windowTitle}</title></head>
        <body>
            <div><img src='${svgUrl}' alt='SVG'></div>
            <div style='position: fixed; top: 0; left: 0;'><img src='${legendUrl}' alt='SVG Legend'></div>
        </body>
        </html>
    `;

    if (svgWindow && !svgWindow.closed) {
        svgWindow.document.open();
        svgWindow.document.write(content);
        svgWindow.document.close();
        svgWindow.focus();
    } else {
        svgWindow = window.open("", "_blank", "width=600,height=800,left=" + (window.screen.width - 600) + ",top=0");
        svgWindow.document.write(content);
    }
}

function showHideTableSummary(anchor) {
  var tableDiv = $("#cl-summary-" + anchor + " .cl-summary-tab");
  if ($(tableDiv).hasClass("hidden")) {
    $(tableDiv).removeClass("hidden");
    $("#cl-summary-" + anchor + " .showhidebutton").text("hide");
  } else {
    $(tableDiv).addClass("hidden")
    $("#cl-summary-" + anchor + " .showhidebutton").text("show");
  }
}

function highlightGene(anchor, locus_tag) {
  $("[id*=-" + anchor.replace("-", "") + "-" + locus_tag.replace(/(:|\.)/g, '-') + "-orf]").each(function(idx, elm) {
    $(elm).css("stroke-width", "5");
    $(elm).mouseover();
  });
  $("#cl-summary-" + anchor +"-tab-" + locus_tag.replace(/(:|\.)/g, '-')).css("background-color", "#CCC");
}

function deHighlightGene(anchor, locus_tag) {
  $("[id*=-" + anchor.replace("-", "") + "-" + locus_tag.replace(/(:|\.)/g, '-') + "-orf]").each(function(idx, elm) {
    $(elm).css("stroke-width", "1");
    $(elm).mouseout();
  });
  $("#cl-summary-" + anchor +"-tab-" + locus_tag.replace(/(:|\.)/g, '-')).css("background-color", "#FFF");
}

function addGeneHiddenSummary() {
    var sibling_col_idx = 0;
    for (var clusterid in geneclusters) {
        if (geneclusters.hasOwnProperty(clusterid)) {
            var cluster = geneclusters[clusterid];
                if ($("#cluster-overview>thead th:contains('Genes')").length === 0) {
                sibling_col_idx = $("#cluster-overview>thead th:contains('Most similar known cluster')").index();
                $("#cluster-overview>thead th:contains('Most similar known cluster')").before("<th style='display: none;'>Genes</th");
            }
            if (sibling_col_idx > 0) {
				var orfs = "";
				for (var i in cluster["orfs"]) {
					var orf = cluster["orfs"][i];
					if (orfs != "") {
						orfs += ", ";
					}
					orfs += orf["locus_tag"];
				}
                $("#cluster-overview>tbody>tr>td>a[href='#" + clusterid + "']").parents("tr").find("td").eq(sibling_col_idx).before("<td style='display: none;'>" + orfs + "</td>");
            }
        }
    }
}

function updateClusterBlastLabel() {
	$(".clusterblast-locustag").each(function(idx, elm) {
		var cbId = $(elm).attr("id");
		var match = cbId.match(/(clusterblast\-)(\d+)_([^0-9]+)(\d+)_(\d+)_(\d+)(\-label)/);
		if (match) {
			var clId = match[4];
			var gId = parseInt(match[6]);
			if (gId < geneclusters["cluster-" + clId]["orfs"].length) {
				$(elm).text(geneclusters["cluster-" + clId]["orfs"][gId]["locus_tag"]);
			}
		}
	});
}

$(document).ready(function() {

    applyHybridColoring();
    update_legends();
    addCdhitSummary();
    addTableSummary();
    addGeneHiddenSummary();
    addCoExpressSummary();

    // Build datatables
    var cur_sep_text = "";
    $('#cluster-overview thead tr').find("th").eq(0).after("<th>Record</th>");
    $('#cluster-overview tbody tr').each(function(idx, elmt) {
        if ($(elmt).hasClass("separator-row")) {
            cur_sep_text = $(elmt).text();
            cur_sep_text = cur_sep_text.replace("The following clusters are", "This row is");
            cur_sep_text = cur_sep_text.replace(/:+$/,"");
        } else {
            $(elmt).attr("title", cur_sep_text);
        }
        $(elmt).find("td").eq(0).after("<td>" + cur_sep_text.split("from record ")[1] + "</td>");
    });
    $('#cluster-overview tbody tr.separator-row').remove();
    $('#cluster-overview tbody tr td:nth-of-type(1)').each(function(idx, elmt) {
        $(elmt).attr("data-order", $(elmt).find("a").text().split(" ")[1]);
        $(elmt).attr("nowrap", true);
    });
    $('#cluster-overview tbody tr td:nth-of-type(4)').each(function(idx, elmt) {
        $(elmt).attr("data-order", $(elmt).text());
    });
    $('#cluster-overview tbody tr td:nth-of-type(5)').each(function(idx, elmt) {
        $(elmt).attr("data-order", $(elmt).text());
    });
    $('#cluster-overview *').css("font-size", "small");
    $('#cluster-overview').DataTable({
        paging: false
    });

    $("#download").click(toggle_downloadmenu);

    addShowHideGenes();

    $(".clbutton").click(function() {
        /* Make sure that even if user missed the link and clicked the
           background we still have the correct anchor */
        var href = $(this).children().first().attr('href');

        if (href === undefined) {
            return;
        }
        window.location.href = href;

        switch_to_cluster();
    }).mouseover(function() {
        /* Set the select cluster label text to cluster type */
        var classes = $(this).attr('class').split(' ');
        if (classes.length < 2) {
          return;
        }
        if (classes[1] == 'separator') {
          return;
        }
        var cluster_type = map_type_to_desc(classes[1]);
        var label = $('#cluster-type');
        label.data("orig_text", label.text());
        label.text(cluster_type + ":");
    }).mouseout(function() {
        /* and reset the select cluster label text */
        var label = $('#cluster-type');
        label.text(label.data("orig_text"));
    });

    $('.clusterblast-selector').change(function() {
        var id = $(this).attr('id').replace('-select', '');
        var url = $(this).val();
        $.get(url, function(data) {
            $('#' + id + '-svg').html(data);
            clusterblast.init(id + '-svg');
            updateClusterBlastLabel();
        }, 'html');
        $('#' + id + '-download').off('click');
        $('#' + id + '-download').click(function () {
            var url = $("#" + id + "-select").val();
            window.open(url, '_blank');
        });
    });

    $('.cluster-rules-header').click(toggle_cluster_rules);

    switch_to_cluster();

});
    </script>

  </body>
</html>
