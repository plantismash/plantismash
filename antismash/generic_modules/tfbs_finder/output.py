from .tfbs_detection import TFBSFinderResults

def generate_details_div(cluster, seq_record, options, js_data, current_html=None):
    print("✅ [TFBS] generate_details_div called")  # For debug

    results = options.extrarecord[seq_record.id].extradata.get("TFBSFinderResults")
    if not results:
        print("⚠️  [TFBS] No results found in extrarecord")
        return current_html

    div = '<div class="tfbs-finder-details">\n'
    div += '<h3>TFBS Finder Summary</h3>\n'

    for record_id, hits in results.hits_by_record.items():
        if not hits:
            div += f"<p>No hits found for record {record_id}</p>\n"
            continue

        div += f"<p><b>Record:</b> {record_id} — {len(hits)} hit(s) detected</p>\n"
        div += "<ul>\n"
        for hit in hits:
            strand = "+" if hit.strand == 1 else "−"
            div += f"<li><b>{hit.name}</b> at <b>{hit.start}</b> ({strand} strand), " \
                   f"score: {hit.score:.1f}/{hit.max_score:.1f}, " \
                   f"confidence: {hit.confidence.name.title()}, " \
                   f"species: {hit.species}</li>\n"
        div += "</ul>\n"

    div += '\n</div>'

    from pyquery import PyQuery as pq
    if current_html is None:
        current_html = pq('<div>')
    current_html.append(pq(div))
    return current_html
