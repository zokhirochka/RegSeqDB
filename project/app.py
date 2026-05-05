import sys
import os
import hmac

# Add the project directory to the path
sys.path.insert(0, '/var/www/html/students_26/Team12/project')

import json
import math
import decimal
from flask import Flask, render_template, request, session, redirect, url_for, jsonify
from regseqDB import RegSeqDB

app = Flask(__name__)

_auth_path = os.path.join(os.path.dirname(__file__), 'auth.json')
with open(_auth_path) as _f:
    _auth = json.load(_f)

app.secret_key = _auth['secret_key']
ADMIN_PASSWORD = _auth['admin_password']


@app.before_request
def require_login():
    if request.endpoint not in ('login', 'logout', 'static') and not session.get('authenticated'):
        return redirect(url_for('login'))


@app.route('/login', methods=['GET', 'POST'])
def login():
    error = None
    if request.method == 'POST':
        password = request.form.get('password', '')
        if ADMIN_PASSWORD and hmac.compare_digest(password, ADMIN_PASSWORD):
            session['authenticated'] = True
            return redirect(url_for('home'))
        error = 'Incorrect password.'
    return render_template('login.html', error=error)


@app.route('/logout')
def logout():
    session.pop('authenticated', None)
    return redirect(url_for('login'))

# ── Decimal JSON Encoder ──────────────────────────────────────────
class DecimalEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, decimal.Decimal):
            return float(obj)
        return super().default(obj)

# ── Database credentials ──────────────────────────────────────────
DB_CREDENTIALS = dict(
    host     = "bioed-new.bu.edu",
    port     = 4253,
    database = "Team12",
    username = "mgdouq",
    password = "ha3563douq",
)

# ── Connect once at startup ───────────────────────────────────────
db = RegSeqDB()
try:
    db.connect(**DB_CREDENTIALS)
    print("✓ Connected to database.")
except Exception as e:
    print(f"✗ Database connection failed: {e}")


# ── Routes ────────────────────────────────────────────────────────

@app.route('/')
def home():
    return render_template('defaultpage.html')


@app.route('/about')
def about():
    return render_template('about.html')


@app.route('/help')
def help_page():
    return render_template('help.html')


@app.route('/singlesearch')
def singlesearch():
    return render_template('singlesearch.html')


@app.route('/comparison')
def comparison():
    return render_template('comparison.html')


# ── Debug ─────────────────────────────────────────────────────────
@app.route('/debug-methods')
def debug_methods():
    import regseqDB, inspect
    methods = [m for m in dir(RegSeqDB()) if not m.startswith('_')]
    return f"<pre>File: {inspect.getfile(regseqDB)}\nMethods: {methods}</pre>"


# ── Single Search ─────────────────────────────────────────────────
@app.route('/search')
def search():
    promoter  = request.args.get('promoter_name', '').strip()
    tf        = request.args.get('TF_name', '').strip()
    condition = request.args.get('Condition', 'glucose').strip()

    if not promoter or not tf:
        return render_template('singlesearch.html')

    try:
        data = db.get_promoter_expr_and_binding(promoter, condition, tf)
    except ValueError as e:
        return render_template('singlesearch.html',
                               error=str(e),
                               promoter_name=promoter,
                               tf_name=tf,
                               condition=condition)

    results  = data['results']   # (sID, num_DNA, num_RNA, energy, affinity)
    rowcount = data['rowcount']

    plot_data  = []
    table_rows = []

    for row in results:
        sID, num_DNA, num_RNA, energy, affinity = row

        if num_DNA is None or num_RNA is None or num_DNA <= 0:
            continue

        plot_data.append({
            "sID":      sID,
            "dna":      float(num_DNA),
            "rna":      float(num_RNA),
            "energy":   float(energy)   if energy   is not None else None,
            "affinity": float(affinity) if affinity is not None else None,
        })

        expr = float(num_RNA) / float(num_DNA)
        table_rows.append({
            "sID":        sID,
            "seq":        "",
            "num_DNA":    num_DNA,
            "num_RNA":    num_RNA,
            "expression": expr,
            "affinity":   float(affinity) if affinity is not None else 0,
            "energy":     float(energy)   if energy   is not None else 0,
        })

    has_data   = len(plot_data) > 0
    graph_json = json.dumps(plot_data, cls=DecimalEncoder) if has_data else None

    # ── Binding coords for genomic locus ─────────────────────────
    locus = None
    try:
        coords = db.get_promoter_binding_coords(promoter, condition, tf)
        if coords['rowcount'] > 0:
            _, _, tss, seq, rnap_start, rnap_stop, tf_start, tf_stop = coords['results'][0]

            # Normalize all positions to a percentage-based canvas (10%–90%)
            all_pos   = [int(tss), int(rnap_start), int(rnap_stop), int(tf_start), int(tf_stop)]
            pos_min   = min(all_pos) - 20   # small padding
            pos_max   = max(all_pos) + 20
            pos_range = pos_max - pos_min

            def to_pct(val):
                return round(((int(val) - pos_min) / pos_range) * 80 + 10, 1)

            locus = {
                "tss":        to_pct(tss),
                "rnap_start": to_pct(rnap_start),
                "rnap_stop":  to_pct(rnap_stop),
                "tf_start":   to_pct(tf_start),
                "tf_stop":    to_pct(tf_stop),
            }
    except Exception:
        pass   # locus stays None; the template should handle missing locus gracefully

    return render_template('singlesearch.html',
        promoter_name  = promoter,
        tf_name        = tf,
        condition      = condition,
        rowcount       = rowcount,
        graph_json     = graph_json,
        graph_tf_expr  = has_data,
        graph_rnap_expr= has_data,
        graph_tf_rnap  = has_data,
        table_rows     = table_rows,
        locus          = locus,
    )


# ── Comparison ────────────────────────────────────────────────────
@app.route('/compare')
def compare():
    promoter   = request.args.get('promoter_name', '').strip()
    tf         = request.args.get('TF_name', '').strip()
    condition1 = request.args.get('Condition1', 'glucose').strip()
    condition2 = request.args.get('Condition2', 'xylose').strip()

    if not promoter or not tf:
        return render_template('comparison.html')

    # FIX: use db.get_condition_comparison() instead of a hardcoded duplicate
    # query. This keeps the SQL in one place and runs proper input validation.
    try:
        results_dict = db.get_condition_comparison(promoter, tf, condition1, condition2)
        results  = results_dict['results']
        rowcount = results_dict['rowcount']
    except ValueError as e:
        return render_template('comparison.html',
                               error=str(e),
                               promoter_name=promoter,
                               tf_name=tf,
                               condition1=condition1,
                               condition2=condition2)
    except Exception as e:
        return f"<pre>QUERY ERROR:\n{e}</pre>", 500

    try:
        plot_data = []
        for row in results:
            sID, dna1, rna1, dna2, rna2, energy, affinity = row
            dna1 = float(dna1) if dna1 is not None else 0
            rna1 = float(rna1) if rna1 is not None else 0
            dna2 = float(dna2) if dna2 is not None else 0
            rna2 = float(rna2) if rna2 is not None else 0
            if dna1 > 0 and rna1 > 0 and dna2 > 0 and rna2 > 0:
                plot_data.append({
                    "sID":      sID,
                    "dna1":     dna1,
                    "rna1":     rna1,
                    "dna2":     dna2,
                    "rna2":     rna2,
                    "energy":   float(energy)   if energy   is not None else None,
                    "affinity": float(affinity) if affinity is not None else None,
                })
        has_data = len(plot_data) > 0
        graph_payload = {"rows": plot_data, "cond1": condition1, "cond2": condition2}
        graph_json = json.dumps(graph_payload, cls=DecimalEncoder) if has_data else None
    except Exception as e:
        return f"<pre>DATA PROCESSING ERROR:\n{e}</pre>", 500

    table_rows = []
    for row in results:
        sID, dna1, rna1, dna2, rna2, energy, affinity = row
        d1 = float(dna1) if dna1 is not None else 0
        r1 = float(rna1) if rna1 is not None else 0
        d2 = float(dna2) if dna2 is not None else 0
        r2 = float(rna2) if rna2 is not None else 0
        expr1 = r1 / d1 if d1 > 0 else None
        expr2 = r2 / d2 if d2 > 0 else None
        logfc = math.log2(expr1 / expr2) if expr1 and expr2 else None
        table_rows.append({
            'sID':      sID,
            'dna1':     d1, 'rna1': r1,
            'dna2':     d2, 'rna2': r2,
            'expr1':    expr1,
            'expr2':    expr2,
            'logfc':    logfc,
            'affinity': float(affinity) if affinity is not None else None,
            'energy':   float(energy)   if energy   is not None else None,
        })

    locus = None
    try:
        coords = db.get_promoter_binding_coords(promoter, condition1, tf)
        if coords['rowcount'] > 0:
            _, _, tss, seq, rnap_start, rnap_stop, tf_start, tf_stop = coords['results'][0]
            all_pos   = [int(tss), int(rnap_start), int(rnap_stop), int(tf_start), int(tf_stop)]
            pos_min   = min(all_pos) - 20
            pos_max   = max(all_pos) + 20
            pos_range = pos_max - pos_min
            def to_pct(val):
                return round(((int(val) - pos_min) / pos_range) * 80 + 10, 1)
            locus = {
                'tss':        to_pct(tss),
                'rnap_start': to_pct(rnap_start),
                'rnap_stop':  to_pct(rnap_stop),
                'tf_start':   to_pct(tf_start),
                'tf_stop':    to_pct(tf_stop),
            }
    except Exception:
        pass

    return render_template('comparison.html',
        promoter_name=promoter, tf_name=tf,
        condition1=condition1, condition2=condition2,
        rowcount1=rowcount, rowcount2=rowcount,
        graph_tf_expr=has_data, graph_json=graph_json,
        table_rows=table_rows, locus=locus,
    )


# ── Autocomplete API ──────────────────────────────────────────────
@app.route('/api/promoters')
def api_promoters():
    q = request.args.get('q', '').strip()
    if not q:
        return jsonify([])
    try:
        res = db.search_promoter_names(q)
        return jsonify([row[0] for row in res['results']])
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/tfs')
def api_tfs():
    q = request.args.get('q', '').strip()
    if not q:
        return jsonify([])
    try:
        res = db.search_tf_names(q)
        return jsonify([row[0] for row in res['results']])
    except Exception as e:
        return jsonify({'error': str(e)}), 500

# ── Run ───────────────────────────────────────────────────────────
if __name__ == '__main__':
    app.run(debug=True)
