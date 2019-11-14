/* Most codes adapted from http://bost.ocks.org/mike/hive/ & http://www.konstantinkashin.com/data/hiv-hiveplot.html */
var Hiveplot = {
  version : "1.0"
};

Hiveplot.draw = function(nodes, rec_id, containerId) {
  var edges = [];
  if (inter_cluster_data.hasOwnProperty(rec_id)) {
    edges = inter_cluster_data[rec_id];
  }
  function link() {
        // A shape generator for Hive links, based on a source and a target.
        // The source and target are defined in polar coordinates (angle and radius).
        // Ratio links can also be drawn by using a startRadius and endRadius.
        // This class is modeled after d3.svg.chord.
        var source = function (d) {
                return d.source;
            },
            target = function (d) {
                return d.target;
            },
            angle = function (d) {
                return d.angle;
            },
            startRadius = function (d) {
                return d.radius;
            }
            endRadius = startRadius,
            arcOffset = -Math.PI / 2;

        function link(d, i) {
            var s = node(source, this, d, i),
                t = node(target, this, d, i),
                x;
            if (t.a - s.a > Math.PI) s.a += 2 * Math.PI;
            // draw cubic bezier curves for nodes on different axes
            if (s.a != t.a) {
                var a1 = (s.a + ((t.a - s.a) / 3)),
                    a2 = (t.a - ((t.a - s.a) / 3)),
                    long = s.r1 > t.r1 ? long = 's' : long = 't';
                    rdiff = 50
                    rp1 = long === 's' ? rp1 = s.r1 + (rdiff * d.multiplier) - (rdiff*((d.maxpair+1)/2)) : rp1 = s.r1 - (rdiff * d.multiplier) + (rdiff*((d.maxpair+1)/2));
                    rp2 = long === 't' ? rp2 = t.r1 + (rdiff * d.multiplier) - (rdiff*((d.maxpair+1)/2)) : rp2 = t.r1 - (rdiff * d.multiplier) + (rdiff*((d.maxpair+1)/2));
                return s.r0 - s.r1 || t.r0 - t.r1
                ? "M" + Math.cos(s.a) * s.r0 + "," + Math.sin(s.a) * s.r0
                + "L" + Math.cos(s.a) * s.r1 + "," + Math.sin(s.a) * s.r1
                + "C" + Math.cos(a1) * rp1 + "," + Math.sin(a1) * rp1
                + " " + Math.cos(a2) * rp2 + "," + Math.sin(a2) * rp2
                + " " + Math.cos(t.a) * t.r1 + "," + Math.sin(t.a) * t.r1
                + "L" + Math.cos(t.a) * t.r0 + "," + Math.sin(t.a) * t.r0
                + "C" + Math.cos(a2) * rp2 + "," + Math.sin(a2) * rp2
                + " " + Math.cos(a1) * rp1 + "," + Math.sin(a1) * rp1
                + " " + Math.cos(s.a) * s.r0 + "," + Math.sin(s.a) * s.r0
                : "M" + Math.cos(s.a) * s.r0 + "," + Math.sin(s.a) * s.r0
                + "C" + Math.cos(a1) * rp1 + "," + Math.sin(a1) * rp1
                + " " + Math.cos(a2) * rp2 + "," + Math.sin(a2) * rp2
                + " " + Math.cos(t.a) * t.r1 + "," + Math.sin(t.a) * t.r1;
            }
            //draw quadratic bezier curves for nodes on same axis
            else {
                var maxopp = 300,       // How far (transversally to axis) should curves go
                    minopp = 100,           // How short should curves be minimum
                    adjdiff = 50,           // How distant should curves be to each other if they have the same source and target
                    distance = Math.abs(s.r0 - t.r0),
                    short = s.r0 < t.r0 ? short = s.r0 : short = t.r0;      //In theory source is always before target, but just in case we used the same model elsewhere
                    // adj model depends on number of links
                    adj = short + (distance/2),
                    adj = adj + (adjdiff*d.multiplier) - (adjdiff*((d.maxpair+1)/2)),
                    // opp model
                    opp = 305 - (((maxopp-minopp)/(400-10))*distance),
                    opp = opp > maxopp ? opp = maxopp : opp = opp;
                    opp = opp < minopp ? opp = minopp : opp = opp;
                    // Pythagorean trig
                    alpha = Math.atan(opp/adj),
                    hyp = opp/Math.sin(alpha),
                    alpha = d.source.type === "pos" ? alpha = alpha : -(Math.PI/2)-alpha;
                return "M" + Math.cos(s.a) * s.r0 + "," + Math.sin(s.a) * s.r0          // P0 = start
                + "Q" + Math.cos(alpha) * hyp + "," + Math.sin(alpha) * hyp             // P1 = control
                + " " + Math.cos(t.a) * t.r1 + "," + Math.sin(t.a) * t.r1;              // P2 = end
            }
        }

        function node(method, thiz, d, i) {
            var node = method.call(thiz, d, i),
                a = +(typeof angle === "function" ? angle.call(thiz, node, i) : angle) + arcOffset,
                r0 = +(typeof startRadius === "function" ? startRadius.call(thiz, node, i) : startRadius),
                r1 = (startRadius === endRadius ? r0 : +(typeof endRadius === "function" ? endRadius.call(thiz, node, i) : endRadius));
            return {
                r0: r0,
                r1: r1,
                a: a
            };
        }

        link.source = function (_) {
            if (!arguments.length) return source;
            source = _;
            return link;
        };

        link.target = function (_) {
            if (!arguments.length) return target;
            target = _;
            return link;
        };

        link.angle = function (_) {
            if (!arguments.length) return angle;
            angle = _;
            return link;
        };

        link.radius = function (_) {
            if (!arguments.length) return startRadius;
            startRadius = endRadius = _;
            return link;
        };

        link.startRadius = function (_) {
            if (!arguments.length) return startRadius;
            startRadius = _;
            return link;
        };

        link.endRadius = function (_) {
            if (!arguments.length) return endRadius;
            endRadius = _;
            return link;
        };

        return link;
  }

  function compare(a, b) {
      if (a.group < b.group)
          return -1;
      if (a.group > b.group)
          return 1;
      if (a.group == b.group) {
          if (a.index < b.index)
              return -1;
          if (a.index > b.index)
              return 1;
          return 0;
      }
  }

  function color(g) {
      //alert($("li.clbutton." + g).text());
      return $("li.clbutton." + g).css("background-color");
  }

  function linkColor(g) {
      if (g == '0') { return "#8e5e00"} //brown
      else if (g == 'x') { return "#6699ff"} //blue
      else { return "#66ffcc"} //green
  }

  function linkWidth(g) {
      if (g == '0') { return "2px"}
      else { return "2px"}
  }

  function linkOpacity(g) {
      if (g == '0') { return ".3"}
      else { return ".3"}
  }

  function degrees(radians) {
      return radians / Math.PI * 180 - 90;
  }

  function capitalizeFirstLetter(string) {
      return string.charAt(0).toUpperCase() + string.slice(1);
  }

  function linkClick(d) {
      if (d3.select(this).classed("static")) {
          svg.selectAll(".static").classed("static", false);
          svgLegend.classed("static", false);
      } else {
          svg.selectAll(".static").classed("static", false);
          svg.selectAll(".link").classed("static", function (p) {
              return p === d;
          });
          svg.selectAll("circle").classed("static", function (p) {
              return p === d.source.node || p === d.target.node;
          });
          svg.selectAll(".link").classed("active", function (p) {
              return p === d;
          });
          svg.selectAll("circle").classed("active", function (p) {
              return p === d.source.node || p === d.target.node;
          });
          svgLegend.classed("static", true);
      }
      svgLegend.classed("hp-hidden", false);
      displayLegendHighlighted(d.source.node);
  }

  function nodeClick(d) {
      if (d3.select(this).classed("static")) {
          svg.selectAll(".static").classed("static", false);
          svgLegend.classed("static", false);
      } else {
          svg.selectAll(".static").classed("static", false);
          svg.selectAll(".link").classed("static", function(p) { return p.source.node === d || p.target.node === d; });
          d3.select(this).classed("static", true);
          svg.selectAll(".link").classed("active", function(p) { return p.source.node === d || p.target.node === d; });
          svg.selectAll("circle").classed("active", function(p) { return p === d; });
          svgLegend.classed("static", true);
      }
      svgLegend.classed("hp-hidden", false);
      displayLegendHighlighted(d);
  }

  function linkMouseover(d) {
      var isStatic = false;
      svg.selectAll(".active.static").forEach(function (elm){ isStatic = elm.length > 0; });
      if (true) {//!isStatic) { //uncomment to disable highglighting when static is present
          svg.selectAll(".link").classed("active", function (p) {
              return p === d;
          });
          svg.selectAll("circle").classed("active", function (p) {
              return p === d.source.node || p === d.target.node;
          });
          if (!svgLegend.classed("static")) {
              svgLegend.classed("hp-hidden", false);
              displayLegendHighlighted(d.source.node);
          }
      }
  }

  function nodeMouseover(d) {
      var isStatic = false;
      svg.selectAll(".active.static").forEach(function (elm){ isStatic = elm.length > 0; });
      if (true) {//!isStatic) { //uncomment to disable highglighting when static is present
          // Highlight the node and connected links on mouseover.
          svg.selectAll(".link").classed("active", function(p) { return p.source.node === d || p.target.node === d; });
          svg.selectAll("circle").classed("active", function(p) { return p === d; });
          if (!svgLegend.classed("static")) {
              svgLegend.classed("hp-hidden", false);
              displayLegendHighlighted(d);
          }
      }
  }

  function mouseout() {
      // Clear any highlighted nodes or links.
      svg.selectAll(".active:not(.static)").classed("active", false);
      if (!svgLegend.classed("static")) {
          //svgLegend.classed("hp-hidden", true);
          displayLegendHighlighted();
      }
  }

  function displayLegendHighlighted(sourceNode) {
      if (sourceNode === undefined) {
          svgLegend.selectAll("*").remove();
          svgLegend.classed("standby", true);
          svgLegend.append("div")
              .style("height", "100%")
              .style("width", "100%")
              .style("line-height", (height - 250) + "px")
              .style("text-align", "center")
              .text("Select and click on links/nodes to show details.");
          return;
      }
      // populate rows
      var data = [];
      svgLegend.classed("standby", false);
      svg.selectAll(".link.active").each(function(d) {
          var source = d.source.node === sourceNode? d.source : d.target;
          var target = d.source.node === sourceNode? d.target : d.source;
          var source_genes = [];
          var target_genes = [];
          d.links.forEach(function(l) {
              var source_tupid = d.source.node === sourceNode? 0 : 1;
              var target_tupid = d.source.node === sourceNode? 1 : 0;
              if (source_genes.indexOf(l[source_tupid]) < 0) {
                  source_genes.push(l[source_tupid]);
              }
              if (target_genes.indexOf(l[target_tupid]) < 0) {
                  target_genes.push(l[target_tupid]);
              }
          });
          var row = {};
          row["source"] = source.node;
          row["target"] = target.node;
          row["score"] = "<span title='{significant inter-relations, genes from source, genes from target}'>";
          row["score"] += "{" + d.links.length + ", " + source_genes.length + ", " + target_genes.length + "}";
          row["score"] += "</span>";
          row["link"] = "<a href='javascript:Hiveplot.showDetail(\"" + rec_id + "\", " + source.id + ", " + target.id + ");'>show</a>";
          data.push(row);
      });
      if (data.length < 1) {
          var row = {};
          row["source"] = sourceNode;
          row["target"] = null;
          row["score"] = "";
          row["link"] = "";
          data.push(row);
      }
      // draw table
      var legActive = svgLegend.classed("hp-hidden", false);
      legActive.selectAll("*").remove();
      var table = legActive.append("table").style("width", "100%")
          .style("margin", "1em 0.5em");
      var thead = table.append('thead').style("border-bottom", "1px solid black");
  		var	tbody = table.append('tbody');
      thead.append('tr')
    		  .selectAll('th')
    		  .data(["Source", "Target", "Attributes <a href='javascript:Hiveplot.showScoreExplanation();'>[?]</a>", "Details"]).enter()
    		  .append('th')
          .style("width", "50px")
  		    .html(function (column) { return column; });
      var rows = tbody.selectAll('tr')
    		  .data(data)
    		  .enter()
    		  .append('tr');
      var cells = rows.selectAll('td')
    		  .data(function (row) {
    		    return ["source", "target", "score", "link"].map(function (column) {
              var tdcol = { column: column, value: row[column], class: "", css: "" };
              if (["source", "target"].indexOf(column) >= 0) {
                  if (row[column] === null) {
                      tdcol["value"] = "";
                      return tdcol;
                  }
                  tdcol["value"] = "<a href='#cluster-" + (row[column].index + 1) + "'>Cluster " + (row[column].index + 1) + "</a>";
                  tdcol["class"] = "clbutton " + row[column].clusterType;
                  tdcol["css"] = d3.selectAll(".clbutton." + row[column].clusterType).attr("style");
                  tdcol["css"] += "; " + d3.selectAll(".clbutton." + row[column].clusterType + " a").attr("style");
              }
              return tdcol;
    		    });})
    		  .enter()
    		  .append('td')
          .attr("class", function (d) { return d.class; })
          .attr("style", function (d) { return d.css; })
  		    .html(function (d) { return d.value; });
  }

  // Assign nodes by group
  var nodesByGroup = d3.nest().key(function (d) {
      return d.axisNumber;
  }).sortKeys(d3.ascending).entries(nodes);

  // sort (apply function defined above)
  nodesByGroup.forEach(function (group) {
      group.values = group.values.sort(compare) // return array in order of
  });

  // now assign plotIndex values
  nodesByGroup.forEach(function (group) {
      var lastGroup = group.values[0].group,
          count = 0;
      group.values.forEach(function (d, i) {
          if (d.group != lastGroup) lastGroup = d.group, count += 2;
          d.plotIndex = count++;
          nodes[d.index].plotIndex = d.plotIndex;
      });
      group.count = count - 1;
  });
  nodesByGroup[0].key = "neg"
  nodesByGroup[1].key = "pos"

  // convert node reference from edge objects into the real object
  var link_counts = {};
  edges.forEach(function(edge) {
      // source and target id (from inter_cluster_data) corresponds equal to cluster id number, however nodes array (from Hiveplot.fetchNodes)
      // is a zero-based array, so cluster1 is in nodes[0]. Thus, the -1 below (and in target).
      var sourceIdx = +edge.source.id -1;
      edge.source["type"] = nodes[sourceIdx].axisNumber == 0? "neg" : "pos";
      edge.source["node"] = nodes[sourceIdx];
      var targetIdx = +edge.target.id -1;
      edge.target["type"] = nodes[targetIdx].axisNumber == 0? "neg" : "pos";
      edge.target["node"] = nodes[targetIdx];
      var pair_idx = "" + Math.min(sourceIdx, targetIdx) + "," + Math.max(sourceIdx, targetIdx);
      if (!(pair_idx in link_counts)) {
          link_counts[pair_idx] = 0;
      }
      link_counts[pair_idx]++;
      edge.multiplier = link_counts[pair_idx];
      edge.highlight = '0';
  });

  edges.forEach(function(edge) {
      var sourceIdx = +edge.source.id -1;
          targetIdx = +edge.target.id -1;
          pair_idx = "" + Math.min(sourceIdx, targetIdx) + "," + Math.max(sourceIdx, targetIdx);
      edge.maxpair = link_counts[pair_idx]
  });

  // set variables
  var width = 400 + (4 * nodesByGroup[1].count),
      height = 425 + (4 * nodesByGroup[0].count),
      innerRadius = 40,
      outerRadius = 540,
      majorAngle = Math.PI / 4,
      minorAngle = Math.PI / 8;
  var angle = d3.scale.ordinal()
      .domain(["neg", "pos"])
      .range([0, Math.PI / 2]);
  var radius = d3.scale.linear()
      .range([innerRadius, 350]);
  radius.domain([0, d3.max([nodesByGroup[0].count, nodesByGroup[1].count])]);
  var formatNumber = d3.format(",d"),
      defaultInfo;

  // Start drawing the svg
  var svg = d3.select("#" + containerId).append("svg")
      .attr("style", "float: left;")
      .attr("class", "d3-hiveplot")
      .attr("width", width)
      .attr("height", height)
      .append("g")
      .attr("transform", "translate(160," + (250 + (4 * nodesByGroup[0].count)) + ")");

  // Draw the axes.
  svg.selectAll(".axis")
      .data(nodesByGroup)
      .enter().append("line")
      .attr("class", "axis")
      .attr("stroke", "rgb(0,0,0)")
      .attr("transform", function (d) {
          return "rotate(" + degrees(angle(d.key)) + ")";
      })
      .attr("x1", radius(-2))
      .attr("x2", function (d) {
          return radius(d.count + 2);
      });

  // Draw the links.
  svg.append("g")
      .attr("class", "links")
      .selectAll(".link")
      .data(edges)
      .enter().append("path")
      .attr("class", "link")
      .attr("d", link()
          .angle(function (d) {
              return angle(d.type);
          })
          .radius(function (d) {
              return radius(d.node.plotIndex);
          }))
      .attr("stroke", function (d) {
          return linkColor(d.highlight);
      })
      .attr("stroke-width", function (d) {
          return linkWidth(d.highlight);
      })
      .attr("stroke-opacity", function (d) {
          return linkOpacity(d.highlight)
      })
      .attr("fill", "transparent")
      .on("mouseover", linkMouseover)
      .on("mouseout", mouseout)
      .on("click", linkClick);

  // Draw the nodes
  svg.append("g")
      .attr("class", "nodes")
      .selectAll(".type")
      .data(nodesByGroup)
      .enter().append("g")
      .attr("class", function (d) {
          return d.key + "nodes";
      })
      .attr("transform", function (d) {
          return "rotate(" + degrees(angle(d.key)) + ")";
      })
      .selectAll("circle")
      .data(function (d) {
          return d.values;
      })
      .enter().append("circle")
      .attr("cx", function (d) {
          return radius(d.plotIndex);
      })
      .attr("r", 4)
      .attr("fill", function (d) {
          return color(d.clusterType);
      })
      .attr("stroke", "#000")
      .attr("stroke-width", "1px")
      .attr("title", function (d) {
          return "Cluster #" + d.index;
      })
      .on("mouseover", nodeMouseover)
      .on("mouseout", mouseout)
      .on("click", nodeClick)
      .classed("");

  // Legends
  var svgLegend =  d3.select("#" + containerId).append("div")
      .classed("d3-hiveplot-leg", true)
      .style("width", "500px")
      .style("height", (height - 250) + "px")
      .classed("hp-hidden", false);//true);
  displayLegendHighlighted();
};

Hiveplot.showDetail = function(rec_id, source_id, target_id) {
	function getGeneColor(gid, clid) {
		var color = "gray";
		var cat = "Other Genes";
	  geneclusters["cluster-" + clid]["orfs"].forEach(function(orf) {
		if (orf["locus_tag"] === gid) {
			var newColor = false;
			if (orf.hasOwnProperty("domains") && orf["domains"].length > 0) {
				var domain = orf["domains"][0];
				for (var i = 0; i < gene_colors.length; i++) {
					var elm = gene_colors[i];
					if (elm.members.indexOf(domain) >= 0) {
						color = elm.color;
						newColor = true;
						cat = elm.label;
						break;
					}
				}
			}
			if (!newColor) {
				if (orf["type"] == "biosynthetic") {
					color = "#810e15";
					cat = "(Other) Biosynthetic Genes";
				}
			}
		}
		});
		return [color, cat];
	}
	var title = "<span>" + "Cluster " + source_id + " (triangle)" + "</span> &times; " + "<span>" + "Cluster " + target_id + " (square)" + "</span>";
    openDialog("<div style='width: 100%; height: 650px; clear: both;'><h2 style='width: 100%; clear: both; text-align: center'>" + title + "</h2><div style='width: 20%; height: 650px; float: left;' id='hv-coexpress-network-legend'></div><div style='width: 70%; height: 650px; float: left;' id='hv-coexpress-network'></div></div>");
    // fetch nodes & links
    var genes = [[], []];
    var links = [];
	var gene_domains = {};
    for (var i in inter_cluster_data[rec_id]) {
      var inter_cluster = inter_cluster_data[rec_id][i];
      if ((inter_cluster.source.id == source_id && inter_cluster.target.id == target_id)
          || (inter_cluster.target.id == source_id && inter_cluster.source.id == target_id)) {
        var src = 0;
        var tgt = 1;
        if (inter_cluster.target.id == source_id && inter_cluster.source.id == target_id) {
          src = 1;
          tgt = 0
        }
        for (var j in inter_cluster.links) {
          var link = inter_cluster.links[j];
          if (genes[src].indexOf(link[0]) < 0) {
            genes[src].push(link[0]);
          }
          if (genes[tgt].indexOf(link[1]) < 0) {
            genes[tgt].push(link[1]);
          }
          links.push([link[0], link[1]]);
        }
        for (var j in inter_cluster.source.links) {
          var link = inter_cluster.source.links[j];
          if (genes[src].indexOf(link[0]) < 0) {
            genes[src].push(link[0]);
          }
          if (genes[src].indexOf(link[1]) < 0) {
            genes[src].push(link[1]);
          }
          links.push([link[0], link[1]]);
        }
        for (var j in inter_cluster.target.links) {
          var link = inter_cluster.target.links[j];
          if (genes[tgt].indexOf(link[0]) < 0) {
            genes[tgt].push(link[0]);
          }
          if (genes[tgt].indexOf(link[1]) < 0) {
            genes[tgt].push(link[1]);
          }
          links.push([link[0], link[1]]);
        }
      }
    }
  	// create a network
    var container = document.getElementById("hv-coexpress-network");
  	var node_array = [];
  	var edge_array = [];
    // build nodes and edges
	var categories = {};
    for (var i in genes[0]) {
      var gid = genes[0][i];
	  var color = "#";
	  var gc = getGeneColor(gid, source_id);
	  if (!categories.hasOwnProperty(gc[1])) {
		 categories[gc[1]] = gc[0];
	  }
      node_array.push({id: gid, label: gid, title: gid, shape: "triangle", color: gc[0]});
    }
    for (var i in genes[1]) {
      var gid = genes[1][i];
	  var gc = getGeneColor(gid, target_id);
	  if (!categories.hasOwnProperty(gc[1])) {
		 categories[gc[1]] = gc[0];
	  }
      node_array.push({id: gid, label: gid, title: gid, shape: "square", color: gc[0]});
    }
    for (var i in links) {
      var link = links[i];
      var dist = 400;
	  var color = "#000";
      if ((genes[0].indexOf(link[0]) > -1 && genes[0].indexOf(link[1]) > -1)
          || (genes[1].indexOf(link[0]) > -1 && genes[1].indexOf(link[1]) > -1)) {
        dist = 200;
		color = "#CCC";
      }
      edge_array.push({from: link[0], to: link[1], length: dist, color: color});
    }
  	var nodes = new vis.DataSet(node_array);
  	var edges = new vis.DataSet(edge_array);
  	var data = {
  		nodes: nodes,
  		edges: edges
  	};
	// draw legends
	$("#hv-coexpress-network-legend").html("<b>Legends:</b><br /><ul style='padding-left: 0; list-style: none;'></ul>");
	for (var cat in categories) {
		if (categories.hasOwnProperty(cat)) {
			$("#hv-coexpress-network-legend ul").append("<li><span class='legend-field' style='background-color: " + categories[cat] + ";'></span>" + cat + "</li>");
		}
	}
	// draw network
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
  				centralGravity: 0,
          avoidOverlap: 1
  			}
  		}
  	};
  	var network = new vis.Network(container, data, options);
    return;
}

Hiveplot.showScoreExplanation = function() {
    var score_explanation = "----------------------------------------------------------" + "\n";
    score_explanation += "Each link in the graph represents a co-expression network among genes in two different clusters. The given score/attributes describes the network the following way: {a, b, c}" + "\n";
    score_explanation += "a) Number of significant inter-cluster co-expressions" + "\n";
    score_explanation += "b) Genes from Source cluster" + "\n";
    score_explanation += "c) Genes from Target cluster" + "\n";
    score_explanation += "----------------------------------------------------------";
    alert(score_explanation);
}

Hiveplot.fetchNodes = function(geneclusters) {
    nodes = [];

    for (cluster_name in geneclusters) {
        node = {};
        cluster = geneclusters[cluster_name];
        node["index"] = cluster["idx"] - 1;
        node["cluster"] = cluster_name;
        node["clusterType"] = cluster["type"];
        node["group"] = cluster["type"];
        nodes.push(node);
    }

    // plot axis number
    var type_sizes = {};
        cur_type = 'foo'
        total_nodes=0;
        half_point=0;
    for (i in nodes) {
        //nodes[i]["axisNumber"] = (nodes[i]["index"] > (nodes.length / 2)) ? 1 : 0; // temporary
        cur_type = nodes[i]['clusterType']
        if (!(cur_type in type_sizes)){
            type_sizes[cur_type]=0;
        }
        type_sizes[cur_type]++;
        total_nodes++;
    }
    var type_sizes_sorted=sortProperties(type_sizes,true)
    half_point=total_nodes/2;
    total_nodes=0;

    for (i in type_sizes_sorted){
        cur_type=type_sizes_sorted[i][0];
        total_nodes = total_nodes + type_sizes_sorted[i][1];
        type_sizes[cur_type] = total_nodes <= half_point ? 0 : 1;       //here we recycle type_sizes into an object specifying the axis of each cluster_type
    }

    for (i in nodes) {
        cur_type = nodes[i]['clusterType'];
        nodes[i]["axisNumber"] = type_sizes[cur_type];
    }

    return nodes.sort(function(a, b) { return a["index"] - b["index"]; });
};


/**
 * Sort object properties (only own properties will be sorted).
 * @param {object} obj object to sort properties
 * @param {bool} isNumericSort true - sort object properties as numeric value, false - sort as string value.
 * @returns {Array} array of items in [[key,value],[key,value],...] format.
 */
function sortProperties(obj, isNumericSort)
{
    isNumericSort=isNumericSort || false; // by default text sort
    var sortable=[];
    for(var key in obj)
        if(obj.hasOwnProperty(key))
            sortable.push([key, obj[key]]);
    if(isNumericSort)
        sortable.sort(function(a, b)
        {
            return a[1]-b[1];
        });
    else
        sortable.sort(function(a, b)
        {
            var x=a[1].toLowerCase(),
                y=b[1].toLowerCase();
            return x<y ? -1 : x>y ? 1 : 0;
        });
    return sortable; // array in format [ [ key1, val1 ], [ key2, val2 ], ... ]
}
