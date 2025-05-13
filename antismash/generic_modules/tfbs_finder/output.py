from .tfbs_detection import TFBSFinderResults
from .tfbs_detection import run_analyses

def generate_details_div(cluster, seq_record, options, js_data, current_html=None):
    """Generates a detail div to insert into the cluster page HTML."""
    results = options.extrarecord[seq_record.id].extradata.get("TFBSFinderResults")
    if not results:
        return current_html

    div = '<div class="tfbs-finder-details">\n'
    div += results.format_html()
    div += '\n</div>'

    from pyquery import PyQuery as pq
    if current_html is None:
        current_html = pq('<div>')
    current_html.append(pq(div))
    return current_html
