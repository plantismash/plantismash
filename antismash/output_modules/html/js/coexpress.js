function draw_coexpress_summary(anchor, rec_id) {
  if (geneclusters.hasOwnProperty(anchor) && geneclusters[anchor].hasOwnProperty("geo")) {
    var gbdiv = $("#" + anchor + ">.content>.coexpress");
    if ($("#" + anchor + "-coexpress-summary").length < 1) {
        gbdiv.append("<div id='" + anchor + "-coexpress-summary'>");
        var gsdiv = $("#" + anchor + "-coexpress-summary");
        $("#" + anchor + "-coexpress-heatmap_header").appendTo("#" + anchor + "-coexpress-summary");
        gsdiv.append("<div style='margin-left: 1em;'>" + "<span style='font-size: xx-small;'>Heatmap was reconstructed using complete-linkage hierarchical clustering algorithm. The number on the scale above the tree represents the Pearson Correlation distance, where 0 = 100% positively correlated (1.00 PCC) and 200 = 100% negatively correlated (-1.00 PCC).</span>" + "</div>");
        gsdiv.append("<div id='" + anchor + "-coexpress-summary-heatmap' style='margin-left: 1em;'>");
        $("#" + anchor + "-coexpress-network_header").appendTo("#" + anchor + "-coexpress-summary");
        gsdiv.append("<div style='margin-left: 1em; margin-top: 0.5em;'><span style='font-size: xx-small;'>(Loading network graph might take a while depending on the complexity of the network. Afterwards, you can scroll, drag and move the nodes and the screen around to give a better view)</span></div>");
        gsdiv.append("<div id='" + anchor + "-coexpress-summary-network' style='margin-left: 1em; margin-top: 1em; border: 1px solid #ccc; background-color: #fff;'>");
    }
  	var show_fluctuation = $("#" + anchor + "-coexpress-show_heatmap_fluctuation").is(":checked");
	// draw legends	
	var divLegend = $("div#" + anchor + ".page div.content div.legend").clone();
	$(divLegend).find("hr").remove();
	$(divLegend).children("div:not(.legend-container)").remove();
    $("#" + anchor + "-cxheatmap-legend").remove();
	$("#" + anchor + "-coexpress-summary-heatmap").after("<div class='legend' id='" + anchor + "-cxheatmap-legend'>" + $(divLegend).html() + "</div>");
    $("#" + anchor + "-cxnetwork-legend").remove();
	$("#" + anchor + "-coexpress-summary-network").after("<div class='legend' id='" + anchor + "-cxnetwork-legend'>" + $(divLegend).html() + "</div>");
	// draw graphs
    draw_heatmap(anchor, rec_id, show_fluctuation);
    draw_correlation_network(anchor, rec_id, $("#" + anchor + "-coexpress-network_cutoff").val());
  	$("#coexpress-rec_id").val(rec_id);
  	$("#coexpress-show_ego").val($("#" + anchor + "-coexpress-show_ego").is(":checked"));
  	$("#coexpress-show_heatmap_fluctuation").val(show_fluctuation);
  }
}

function draw_heatmap(anchor, rec_id, show_fluctuation = false) {
  if (!geneclusters[anchor]["geo"].hasOwnProperty(rec_id)) {
      $("#" + anchor + "-coexpress-summary-heatmap").html("<i><br />No coexpression data found for this cluster on " + rec_id + ", cannot draw heatmap.</i>");
      return;
  }
  var inchlib = new InCHlib({
      target: anchor + "-coexpress-summary-heatmap",
      metadata: false,
      metadata_colors: "RdGy",
      heatmap_colors: "BW",
      heatmap_font_color: "transparent",
      draw_row_ids: true,
      navigation_toggle: {
        color_scale: false,
        export_button: false,
        hint_button: false,
        filter_button: false,
		color_scale: true,
      },
      row_color_function: function(row_id) {
          var color = "gray";
          geneclusters[anchor]["orfs"].forEach(function(orf) {
            if (orf["locus_tag"] === row_id) {
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
            }
          });
          return color;
      }
  });
  if (show_fluctuation) {
	  inchlib.settings.independent_columns = false;
	  inchlib.settings.independent_rows = true;
	  inchlib.settings.heatmap_colors = "BW";
  } else {
	  inchlib.settings.independent_columns = true;
	  inchlib.settings.independent_rows = false;
	  inchlib.settings.heatmap_colors = "YlOrR";
  }
  inchlib.read_data(geneclusters[anchor]["geo"][rec_id]);
  set_heatmap_scales(rec_id, inchlib);
  inchlib.draw();
}

function set_heatmap_scales(rec_id, inchlib) {
	var ranges = {};
	for(i = 0; i < geo_dataset_info.length; i++) {
		var info = geo_dataset_info[i];
		if (info["id"] == rec_id) {
			ranges = info["ranges"];
			break;
		}
	}
	inchlib.settings.ranges = [];
	for(i = 0; i < inchlib.json.data.feature_names.length; i++) {
		inchlib.settings.ranges.push(ranges[inchlib.json.data.feature_names[i]]);
	}
}

function draw_correlation_network(anchor, rec_id, cutoff) {
  $("#coexpress-network_cutoff").val(cutoff);
  if (!geneclusters[anchor]["geo"].hasOwnProperty(rec_id)) {
      $("#" + anchor + "-coexpress-summary-network").html("<i>No coexpression data found for this cluster on " + rec_id + ", cannot draw network graph.</i>");
      return;
  }
  // create an array with nodes
  var node_array = []
  for (var i in geneclusters[anchor]["orfs"]) {
    var orf = geneclusters[anchor]["orfs"][i];
    var color = "gray";
    try {
      color = get_gene_color(orf);
    }
    catch(err) {
      if (orf["type"] === "biosynthetic") {
        color = "#810e15";
      }
    }
    var font = {strokeColor: "white", strokeWidth: 5};
    var gene_title = gene_id;
    node_array.push({id: orf["locus_tag"], label: orf["locus_tag"], title: gene_title, color: color, font: font, shape: "box"})
  }

  // create an array with edges
  var edge_array = []
  for (var i in geneclusters[anchor]["orfs"]) {
      var orf = geneclusters[anchor]["orfs"][i];
      if (orf.hasOwnProperty("geo")) {
          for (var j in orf["geo"]) {
              var geo = orf["geo"][j];
              if (geo["rec_id"] == rec_id) {
                  var dist = geo["dist"];
                  for (var gene_id in dist) {
                      var external_node = true;
                      for (var x in geneclusters[anchor]["orfs"]) {
                          if (geneclusters[anchor]["orfs"][x]["locus_tag"] == gene_id) {
                              external_node = false;
                              break;
                          }
                      }
                      if (external_node) {
            						  if (!$("#" + anchor + "-coexpress-show_ego").is(":checked")) {
            							  break;
            						  }
                      }
                      if (gene_id != orf["locus_tag"] && (dist[gene_id] < cutoff)) {
                          var edge_exist = false;
                          for (var j in edge_array) {
                              var fr = edge_array[j].from;
                              var to = edge_array[j].to;
                              if ((fr == orf["locus_tag"] && to == gene_id) || (fr == gene_id && to == orf["locus_tag"])) {
                                  edge_exist = true;
                                  break;
                              }
                          }
                          if (!edge_exist) {
                              var slength = 150;
                              if (external_node) {
                                  // node come from neighboring genes
                                  var node_exist = false;
                                  for (var x in node_array) {
                                      if (node_array[x]["id"] == gene_id) {
                                          node_exist = true;
                                          break;
                                      }
                                  }
                                  if (!node_exist) {
                  		      ncolor = "#ccc";
                                      var gene_title = gene_id;
                                      node_array.push({id: gene_id, label: gene_id, title: gene_title, color: ncolor, font: "10px arial black", shape: "ellipse", border: "1px black"});
                                  }
                                  slength = 250;
                              }
                              var pcc_value = 1 - (parseFloat(dist[gene_id]) / 100);
                              var edge_title = orf["locus_tag"] + " <-> " + gene_id + "; d: " + dist[gene_id] + " (PCC " + pcc_value.toFixed(2) + ")";
                              edge_array.push({from: orf["locus_tag"], to: gene_id, title: edge_title, color: "black", length: slength});
                          }
                      }
                  }
              }
          }
      }
  }

  // rename nodes in external clusters, set dashed edges for other remote nodes
  node_array.forEach(function(node) {
      for (var anc in geneclusters) {
        if (geneclusters.hasOwnProperty(anc)) {
          for (var x in geneclusters[anc]["orfs"]) {
              if (geneclusters[anc]["orfs"][x]["locus_tag"] == node["id"]) {
                  if (anc !== anchor) {
					var legobj = get_legend_obj(geneclusters[anc]["orfs"][x]);
					if (legobj !== null) {
						$("#" + anchor + "-cxnetwork-legend  div.legend-container").each(function(idx, elm) {
							if ($(elm).find(".legend-label").text() == legobj.label) {
								$(elm).removeClass("hidden");
							}
						});
					}
                    ncolor = get_gene_color(geneclusters[anc]["orfs"][x]);
                    node.label =  "[" + anc.split("-")[1] + "]\n" + node.label;
                    node.color = ncolor;
                    node.font =  {strokeColor: "white", strokeWidth: 5, color: "red"};
                  } else {
                    node.label =  "\n" + node.label + "\n";
                  }
                  return;
              }
          }
        }
      }
      for (var x in edge_array) {
          if ([edge_array[x]["from"], edge_array[x]["to"]].indexOf(node["id"]) > -1) {
              edge_array[x]["dashes"] = [5];
              edge_array[x]["color"] = "#bbb";
              edge_array[x]["length"] = 500;
          }
      }
  });

  // create a network
  var nodes = new vis.DataSet(node_array);
  var edges = new vis.DataSet(edge_array);
  var container = document.getElementById(anchor + "-coexpress-summary-network");
  var data = {
    nodes: nodes,
    edges: edges
  };
  var options = {
    height: "600px",
    edges: {
      smooth: {
        type:  "straightCross",
        roundness: 0
      }
    },
    physics: {
      enabled: true,
    	barnesHut: {
        springConstant: 0.01,
        centralGravity: 0
      }
    }
  };
  var network = new vis.Network(container, data, options);
}

function drawHivePlot(rec_id) {
    $("#cluster_exp_hiveplot").html("");
    for (var i in inter_cluster_data[rec_id]) { // fix link nodes order
      var sourceId = inter_cluster_data[rec_id][i].source.id;
      var targetId = inter_cluster_data[rec_id][i].target.id;
      var links = inter_cluster_data[rec_id][i].links;
      for (var j in links) {
        var tempLink = links[j];
        for (var x in geneclusters["cluster-" + sourceId]["orfs"]) {
          if (geneclusters["cluster-" + sourceId]["orfs"][x]["locus_tag"] == tempLink[1]) {
            links[j] = [tempLink[1], tempLink[0]];
            break;
          }
        }
      }
    }
    Hiveplot.draw(Hiveplot.fetchNodes(geneclusters), rec_id, "cluster_exp_hiveplot");
}

function drawSignalPlot(rec_id) {
    $("#cluster_exp_signalplot").html("");
    return;
    var i = 0;
    seq_records = Object.keys(geo_signals[rec_id]).sort();
    for (var i in seq_records) {
      var seq_record = seq_records[i];
      if (geo_signals[rec_id].hasOwnProperty(seq_record)) {
        $("#cluster_exp_signalplot").append("<h5 style=' width: 100%; text-align: center; clear: both;'>Signal plot: " + seq_record + "</h5><div style='width: 100%; margin-bottom: 4em; overflow-x: scroll;' id='cluster_exp_signalplot_" + i  + "'/>");
        var signal_data = geo_signals[rec_id][seq_record];
        var lines = [];
        lines.push({
          x: [],
          y: [],
          z: [],
          name: seq_record,
          mode: 'lines',
          type: 'scatter',
          marker: {
            color: "#cccccc"
          }
        });

        var drawingNewLine = false;
        for (var n = 0; n < signal_data.length; n++) {
          lines[0].x.push(signal_data[n][0]);
          lines[0].y.push(signal_data[n][1]);
          lines[0].z.push(signal_data[n][2]);
          if (signal_data[n][2] >= 0) {
            if (!drawingNewLine) {
              lines.push({
                x: [],
                y: [],
                name: 'Cluster ' + signal_data[n][2],
                mode: 'lines',
                type: 'scatter'
              });
              drawingNewLine = true;
            }
            lines[lines.length - 1].x.push(signal_data[n][0]);
            lines[lines.length - 1].y.push(signal_data[n][1]);
          } else {
            if (drawingNewLine) {
              drawingNewLine = false;
            }
          }
        }

        var layout = {
          title: "",
          height: 350,
          width: signal_data.length / 4,
          xaxis: {
            title: "Location"
          },
          yaxis: {
            title: "R",
            range: [0, 1]
          }
        };

        Plotly.newPlot("cluster_exp_signalplot_" + i, lines, layout);
        document.getElementById("cluster_exp_signalplot_" + i).on('plotly_click', function(data){
            var pointNumber = data.points[0].pointNumber;
            var allData = data.points[0].data;
            if (allData.z[pointNumber] > 0) {
                window.location.href =  "#cluster-" + allData.z[pointNumber];
                switch_to_cluster();
            } else {
                var chr_id = allData.name;
                var fr = parseInt(allData.x[Math.max(pointNumber - 3, 0)]);
                var to = parseInt(allData.x[Math.min(pointNumber + 3, allData.x.length - 1)]);
                window.open("http://www.ncbi.nlm.nih.gov/projects/sviewer/?Db=gene&DbFrom=protein&Cmd=Link&noslider=1&id=" + chr_id + "&from=" + fr + "&to=" + to);
            }
        });
        i++;
      }
    }
}
