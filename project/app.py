# ============================================================
# app.py  —  Flask backend for RegSeqDB
# ============================================================
from flask import Flask, request, render_template
from regseqdb import RegSeqDB
from plotting import (
    plot_tf_affinity_vs_expression,
    plot_rnap_energy_vs_expression,
    plot_tf_affinity_vs_rnap_energy,
    plot_compare_tf_affinity,
    plot_compare_rnap_energy,
    plot_compare_tf_rnap,
    fig_to_json,
)

app = Flask(__name__)
app.config['APPLICATION_ROOT'] = '/students_26/Team12/project/app'

# ── DB credentials ────────────────────────────────────────────
DB_CREDENTIALS = dict(
    host     = "bioed-new.bu.edu",
    port     = 4253,
    database = "Team12",
    username = "mgdouq",
    password = "ha3563douq",
)

db = RegSeqDB()
db.connect(**DB_CREDENTIALS)


def get_db():
    """Return DB connection, reconnecting if it has gone away."""
    global db
    try:
        db.cursor.execute("SELECT 1")
    except Exception:
        db = RegSeqDB()
        db.connect(**DB_CREDENTIALS)
    return db


# ══════════════════════════════════════════════════════════════
# Static page routes
# ══════════════════════════════════════════════════════════════

@app.route('/')
def home():
    return render_template('defaultpage.html')

@app.route('/about')
def about():
    return render_template('about.html')

@app.route('/help')
def help_page():
    return render_template('help.html')
    
@app.route('/comparison')
def comparison():
    return render_template('comparison.html')

@app.route('/singlesearch')
def singlesearch():
    return render_template('singlesearch.html')


# ══════════════════════════════════════════════════════════════
# Single search route
# ══════════════════════════════════════════════════════════════

@app.route('/search', methods=['GET'])
def search():

    promoter  = request.args.get("promoter_name", "").strip()
    tf        = request.args.get("TF_name", "").strip()
    condition = request.args.get("Condition", "").strip()

    if not promoter or not tf or not condition:
        return render_template("singlesearch.html")

    try:
        data = get_db().get_promoter_expr_and_binding(
            promoter=promoter,
            condition=condition,
            tf=tf,
            include_rnap=True,
        )
        results  = data["results"]
        colnames = data["colnames"]
        rowcount = data["rowcount"]
    except ValueError as e:
        return render_template("singlesearch.html", error=str(e))

    graph_tf_expr   = None
    graph_rnap_expr = None
    graph_tf_rnap   = None

    if rowcount > 0:
        graph_tf_expr   = fig_to_json(plot_tf_affinity_vs_expression(results, colnames))
        graph_rnap_expr = fig_to_json(plot_rnap_energy_vs_expression(results, colnames))
        graph_tf_rnap   = fig_to_json(plot_tf_affinity_vs_rnap_energy(results, colnames))

    col_idx = {name: i for i, name in enumerate(colnames)}
    table_rows = []
    for row in results:
        num_dna    = row[col_idx["num_DNA"]]
        num_rna    = row[col_idx["num_RNA"]]
        energy     = row[col_idx.get("energy",   -1)]
        affinity   = row[col_idx.get("affinity", -1)]
        expression = (num_rna / num_dna) if num_dna > 0 else 0
        table_rows.append({
            "sID":        row[col_idx["sID"]],
            "seq":        row[col_idx["seq"]] if "seq" in col_idx else "—",
            "num_DNA":    num_dna,
            "num_RNA":    num_rna,
            "expression": expression,
            "affinity":   affinity if affinity is not None else 0,
            "energy":     energy   if energy   is not None else 0,
        })

    return render_template(
        "singlesearch.html",
        promoter_name   = promoter,
        tf_name         = tf,
        condition       = condition,
        rowcount        = rowcount,
        graph_tf_expr   = graph_tf_expr,
        graph_rnap_expr = graph_rnap_expr,
        graph_tf_rnap   = graph_tf_rnap,
        table_rows      = table_rows,
        sigma_factor    = "σ70",
    )


# ══════════════════════════════════════════════════════════════
# Comparison route — handles BOTH form display AND results
# ══════════════════════════════════════════════════════════════

@app.route('/compare', methods=['GET'])
def compare():

    promoter1  = request.args.get("promoter_name1", "").strip()
    promoter2  = request.args.get("promoter_name2", "").strip()
    tf1        = request.args.get("TF_name1", "").strip()
    tf2        = request.args.get("TF_name2", "").strip()
    condition  = request.args.get("Condition", "").strip()

    # If form not submitted, show empty comparison page
    if not promoter1 or not promoter2 or not tf1 or not tf2 or not condition:
        return render_template("comparison.html")

    error = None
    data1 = data2 = None

    try:
        data1 = get_db().get_promoter_expr_and_binding(
            promoter=promoter1, condition=condition, tf=tf1, include_rnap=True)
    except ValueError as e:
        error = f"Promoter 1 error: {e}"

    try:
        data2 = get_db().get_promoter_expr_and_binding(
            promoter=promoter2, condition=condition, tf=tf2, include_rnap=True)
    except ValueError as e:
        error = (error + f" | Promoter 2 error: {e}") if error else f"Promoter 2 error: {e}"

    if error:
        return render_template("comparison.html", error=error)

    results1, colnames1, rowcount1 = data1["results"], data1["colnames"], data1["rowcount"]
    results2, colnames2, rowcount2 = data2["results"], data2["colnames"], data2["rowcount"]

    graph_compare_tf   = None
    graph_compare_rnap = None
    graph_compare_both = None

    if rowcount1 > 0 and rowcount2 > 0:
        graph_compare_tf   = fig_to_json(plot_compare_tf_affinity(
            results1, colnames1, promoter1,
            results2, colnames2, promoter2,
        ))
        graph_compare_rnap = fig_to_json(plot_compare_rnap_energy(
            results1, colnames1, promoter1,
            results2, colnames2, promoter2,
        ))
        graph_compare_both = fig_to_json(plot_compare_tf_rnap(
            results1, colnames1, promoter1,
            results2, colnames2, promoter2,
        ))

    return render_template(
        "comparison.html",
        promoter1          = promoter1,
        promoter2          = promoter2,
        tf1                = tf1,
        tf2                = tf2,
        condition          = condition,
        rowcount1          = rowcount1,
        rowcount2          = rowcount2,
        graph_compare_tf   = graph_compare_tf,
        graph_compare_rnap = graph_compare_rnap,
        graph_compare_both = graph_compare_both,
    )


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5001, debug=True)
