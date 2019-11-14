/* Copyright 2012 Kai Blin. Licensed under the Apache License v2.0, see LICENSE file */

var svgene = {
    version: "0.1.5",
    label_height: 14,
    extra_label_width: 100,
    unique_id: 0
};

svgene.geneArrowPoints = function (orf, height, offset, border, scale) {
  var top_ = offset + svgene.label_height + border;
  var bottom = offset + svgene.label_height + height - border;
  var middle = offset + svgene.label_height + (height/2);
  if (orf.strand == 1) {
      var start = scale(orf.start);
      var box_end = Math.max(scale(orf.end) - (2*border), start);
      var point_end = scale(orf.end);
      points  = "" + start + "," + top_;
      points += " " + box_end + "," + top_;
      points += " " + point_end + "," + middle;
      points += " " + box_end + "," + bottom;
      points += " " + start + "," + bottom;
      return points;
  }
  if (orf.strand == -1) {
      var point_start = scale(orf.start);
      var end = scale(orf.end);
      var box_start = Math.min(scale(orf.start) + (2*border), end);
      points = "" + point_start + "," + middle;
      points += " " + box_start + "," + top_;
      points += " " + end + "," + top_;
      points += " " + end + "," + bottom;
      points += " " + box_start + "," + bottom;
      return points;
  }
};

svgene.drawOrderedClusterOrfs = function(cluster, chart, all_orfs, scale,
                                         i, idx, height, width,
                                         single_cluster_height, offset) {
  chart.append("line")
    .attr("x1", 0)
    .attr("y1", (single_cluster_height * i) + svgene.label_height + (height/2))
    .attr("x2", width)
    .attr("y2", (single_cluster_height * i) + svgene.label_height + (height/2))
    .attr("class", "svgene-line");
  chart.selectAll("polygon")
    .data(all_orfs)
  .enter().append("polygon")
    .attr("points", function(d) { return svgene.geneArrowPoints(d, height, (single_cluster_height * i), offset, scale); })
    .attr("class", function(d) { return "svgene-type-" + d.type + " svgene-orf"; })
    .attr("id", function(d) { return idx + "-cluster" + cluster.idx + "-" + svgene.tag_to_id(d.locus_tag) + "-orf"; })
    .attr("style", function(d) { if (d.color !== undefined) { return "fill:" + d.color; } })
  chart.selectAll("text")
    .data(all_orfs)
  .enter().append("text")
    .attr("x", function(d) { return scale(d.start); })
    .attr("y", (single_cluster_height * i) + svgene.label_height + offset/2)
    .attr("class", "svgene-locustag")
    .attr("id", function(d) { return idx + "-cluster" + cluster.idx + "-" + svgene.tag_to_id(d.locus_tag) + "-label"; })
    .text(function(d) { return d.locus_tag; });

};

svgene.drawUnorderedClusterOrfs = function(cluster, chart, all_orfs, scale,
                                           i, idx, height, width,
                                           single_cluster_height, offset) {
  chart.selectAll("rect")
    .data(all_orfs)
  .enter().append("rect")
    .attr("x", function(d) { return scale(d.start);})
    .attr("y", (single_cluster_height * i) + svgene.label_height + offset)
    .attr("height", height - (2 * offset))
    .attr("width", function(d) { return scale(d.end) - scale(d.start)})
    .attr("rx", 3)
    .attr("ry", 3)
    .attr("class", function(d) { return "svgene-type-" + d.type + " svgene-orf"; })
    .attr("id", function(d) { return idx + "-cluster" + cluster.idx + "-" + svgene.tag_to_id(d.locus_tag) + "-orf"; })
    .attr("style", function(d) { if (d.color !== undefined) { return "fill:" + d.color; } })
  chart.selectAll("text")
    .data(all_orfs)
  .enter().append("text")
    .attr("x", function(d) { return scale(d.start); })
    .attr("y", (single_cluster_height * i) + svgene.label_height + offset/2)
    .attr("class", "svgene-locustag")
    .attr("id", function(d) { return idx + "-cluster" + cluster.idx + "-" + svgene.tag_to_id(d.locus_tag) + "-label"; })
    .text(function(d) { return d.locus_tag; });
};

svgene.drawClusters = function(id, clusters, height, width) {
  var container = d3.select("#" + id);
  var single_cluster_height = height + svgene.label_height;
  container.selectAll("svg").remove();
  container.selectAll("div").remove();
  var chart = container.append("svg")
    .attr("height", single_cluster_height * clusters.length)
    .attr("width", width + svgene.extra_label_width);
  var all_orfs = [];

  for (i=0; i < clusters.length; i++) {
      var cluster = clusters[i];
      for (j=0; j < cluster.orfs.length; j++) {
          orf = cluster.orfs[j];
          if (orf.type == "biosynthetic") {
              all_orfs.push(orf)
          } else {
              all_orfs.splice(0, 0, orf)
          }
      }
      var idx = svgene.unique_id++;
      var offset = height/10;
      var x = d3.scale.linear()
        .domain([cluster.start, cluster.end])
        .range([0, width]);
      if (cluster.unordered) {
          svgene.drawUnorderedClusterOrfs(cluster, chart, all_orfs, x,
                                          i, idx, height, width,
                                          single_cluster_height, offset);
      } else {
          svgene.drawOrderedClusterOrfs(cluster, chart, all_orfs, x,
                                        i, idx, height, width,
                                        single_cluster_height, offset);
      }
      container.selectAll("div")
        .data(all_orfs)
      .enter().append("div")
        .attr("class", "svgene-tooltip")
        .attr("id", function(d) { return idx + "-cluster" + cluster.idx + "-" + svgene.tag_to_id(d.locus_tag) + "-tooltip"; })
        .html(function(d) { return d.description});
  }
  for (i=0; i < clusters.length; i++) {
      var cluster = clusters[i];
      if (cluster.label !== undefined) {
        chart.append("text")
            .text(cluster.label)
            .attr("class", "svgene-clusterlabel")
            .attr("x", function() { return width + svgene.extra_label_width - this.getComputedTextLength() - 5})
            .attr("y", function() { return (single_cluster_height * i) + svgene.label_height } )
            .attr("font-size", svgene.label_height);
      }
  }
  svgene.init();
};

svgene.tag_to_id = function(tag) {
    return tag.replace(/(:|\.)/g, '-').replace(/-orf/g, '_orf');
}


svgene.tooltip_handler = function(ev) {
    var id = $(this).attr("id").replace("-orf", "-tooltip");
    var tooltip = $("#"+id);

    if (svgene.active_tooltip) {
        svgene.active_tooltip.hide();
    }
    svgene.active_tooltip = tooltip;

    if (tooltip.css("display") == 'none') {
        var offset = $(this).offset();
        tooltip.css("left", offset.left + 10);
        var this_parent = $(this).parent();
        var top_offset = this_parent.height()/(this_parent.children('line').length * 2);
        tooltip.css("top", offset.top + top_offset);
        tooltip.show();
        tooltip.click(function(){$(this).hide()});
        var timeout = setTimeout(function(){ tooltip.slideUp("fast") }, 5000);
        tooltip.data("timeout", timeout);
        tooltip.mouseover(function() {
            clearTimeout(tooltip.data("timeout"));
        }).mouseout(function() {
            timeout = setTimeout(function(){ tooltip.slideUp("fast") }, 5000);
            tooltip.data("timeout", timeout);
        });
    } else {
        tooltip.hide();
    }
};

svgene.init = function() {
    $(".svgene-orf").mouseover(function(e) {
        var id = $(this).attr("id").replace("-orf", "-label");
        $("#"+id).show();
    }).mouseout(function(e) {
        var id = $(this).attr("id").replace("-orf", "-label");
        $("#"+id).hide();
    }).click(svgene.tooltip_handler);
    $(".svgene-textarea").click(function(event) {
        event.stopPropagation();
    });

};
