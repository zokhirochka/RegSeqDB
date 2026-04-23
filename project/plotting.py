"""
plotting.py
-----------
Visualization module for the RegSeqDB website project.
Based on Dr. James Galagan's lab data (BU Microbiology).
"""

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# ---------------------------------------------------------------------------
# Color palette
# ---------------------------------------------------------------------------
PALETTE       = ["#00bfff"]
COLOR_P1      = "#00bfff"   # blue  — promoter 1
COLOR_P2      = "#f97066"   # red   — promoter 2
BG_COLOR      = "#0f1117"
PAPER_COLOR   = "#1a1d27"
FONT_COLOR    = "#e8eaf0"
GRID_COLOR    = "#2e3147"

BASE_LAYOUT = dict(
    paper_bgcolor=PAPER_COLOR,
    plot_bgcolor=BG_COLOR,
    font=dict(color=FONT_COLOR, family="monospace", size=13),
    margin=dict(l=60, r=30, t=60, b=60),
    hoverlabel=dict(bgcolor="#2e3147", font_size=12),
    xaxis=dict(gridcolor=GRID_COLOR, zerolinecolor=GRID_COLOR),
    yaxis=dict(gridcolor=GRID_COLOR, zerolinecolor=GRID_COLOR),
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _to_df(results, colnames):
    return pd.DataFrame(results, columns=colnames)


def _compute_expression(df):
    df = df[df["num_DNA"] > 0].copy()
    df["expression"] = df["num_RNA"] / df["num_DNA"]
    df = df[df["expression"] > 0]
    return df


def _apply_log_y(fig, enabled):
    if enabled:
        fig.update_yaxes(type="log", gridcolor=GRID_COLOR, zerolinecolor=GRID_COLOR)
    return fig


# ---------------------------------------------------------------------------
# Plot 1 — TF Affinity vs. Expression
# ---------------------------------------------------------------------------
def plot_tf_affinity_vs_expression(results, colnames, title=None, log_scale=True):
    df = _to_df(results, colnames)
    df = _compute_expression(df)

    color_col = "tf_name" if "tf_name" in df.columns else None

    fig = px.scatter(
        df,
        x="affinity",
        y="expression",
        color=color_col,
        hover_data=["sID"] + (["tf_name"] if color_col else []),
        labels={
            "affinity":   "TF Binding Affinity",
            "expression": "Expression (RNA / DNA)",
            "tf_name":    "Transcription Factor",
        },
        color_discrete_sequence=PALETTE,
        opacity=0.75,
        trendline="ols",
        trendline_options=dict(log_y=log_scale),
    )

    fig.update_traces(marker=dict(size=7, line=dict(width=0.5, color="rgba(255,255,255,0.19)")))
    fig.update_layout(
        **BASE_LAYOUT,
        title=dict(
            text=title or "TF Binding Affinity vs. Promoter Expression",
            font=dict(size=16, color=FONT_COLOR),
        ),
    )
    _apply_log_y(fig, log_scale)
    return fig


# ---------------------------------------------------------------------------
# Plot 2 — RNAP Binding Energy vs. Expression
# ---------------------------------------------------------------------------
def plot_rnap_energy_vs_expression(results, colnames, title=None, log_scale=True):
    df = _to_df(results, colnames)
    df = _compute_expression(df)

    color_col = "sigma" if "sigma" in df.columns else None

    fig = px.scatter(
        df,
        x="energy",
        y="expression",
        color=color_col,
        hover_data=["sID"],
        labels={
            "energy":     "RNAP Binding Energy",
            "expression": "Expression (RNA / DNA)",
            "sigma":      "Sigma Factor",
        },
        color_discrete_sequence=PALETTE,
        opacity=0.75,
        trendline="ols",
        trendline_options=dict(log_y=log_scale),
    )

    fig.update_traces(marker=dict(size=7, symbol="diamond", line=dict(width=0.5, color="rgba(255,255,255,0.19)")))
    fig.update_layout(
        **BASE_LAYOUT,
        title=dict(
            text=title or "RNAP Binding Energy vs. Promoter Expression",
            font=dict(size=16, color=FONT_COLOR),
        ),
    )
    _apply_log_y(fig, log_scale)
    return fig


# ---------------------------------------------------------------------------
# Plot 3 — Expression by Condition
# ---------------------------------------------------------------------------
def plot_expression_by_condition(results, colnames, title=None, log_scale=True):
    df = _to_df(results, colnames)
    df = _compute_expression(df)

    fig = px.box(
        df,
        x="cond",
        y="expression",
        color="cond",
        points="all",
        hover_data=["sID"],
        labels={
            "cond":       "Experimental Condition",
            "expression": "Expression (RNA / DNA)",
        },
        color_discrete_sequence=PALETTE,
    )

    fig.update_traces(marker=dict(size=4, opacity=0.6))
    fig.update_layout(
        **BASE_LAYOUT,
        showlegend=False,
        title=dict(
            text=title or "Promoter Expression Across Conditions",
            font=dict(size=16, color=FONT_COLOR),
        ),
    )
    _apply_log_y(fig, log_scale)
    return fig


# ---------------------------------------------------------------------------
# Plot 4 — TF Affinity vs. RNAP Energy
# ---------------------------------------------------------------------------
def plot_tf_affinity_vs_rnap_energy(results, colnames, title=None):
    df = _to_df(results, colnames)

    color_col  = "tf_name" if "tf_name" in df.columns else None
    symbol_col = "sigma"   if "sigma"   in df.columns else None

    fig = px.scatter(
        df,
        x="affinity",
        y="energy",
        color=color_col,
        symbol=symbol_col,
        hover_data=["sID"],
        labels={
            "affinity": "TF Binding Affinity",
            "energy":   "RNAP Binding Energy",
            "tf_name":  "Transcription Factor",
            "sigma":    "Sigma Factor",
        },
        color_discrete_sequence=PALETTE,
        opacity=0.8,
    )

    fig.update_traces(marker=dict(size=8, line=dict(width=0.5, color="rgba(255,255,255,0.19)")))
    fig.update_layout(
        **BASE_LAYOUT,
        title=dict(
            text=title or "TF Affinity vs. RNAP Binding Energy",
            font=dict(size=16, color=FONT_COLOR),
        ),
    )
    return fig


# ---------------------------------------------------------------------------
# Plot C1 — Comparison: TF Affinity vs. Expression (two promoters overlaid)
# ---------------------------------------------------------------------------
def plot_compare_tf_affinity(
    results1, colnames1, label1,
    results2, colnames2, label2,
    log_scale=True,
):
    df1 = _compute_expression(_to_df(results1, colnames1))
    df2 = _compute_expression(_to_df(results2, colnames2))

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=df1["affinity"], y=df1["expression"],
        mode="markers",
        name=label1,
        marker=dict(color=COLOR_P1, size=6, opacity=0.7,
                    line=dict(width=0.4, color="rgba(255,255,255,0.2)")),
        hovertemplate=f"<b>{label1}</b><br>Affinity: %{{x:.4f}}<br>Expression: %{{y:.3f}}<extra></extra>",
    ))

    fig.add_trace(go.Scatter(
        x=df2["affinity"], y=df2["expression"],
        mode="markers",
        name=label2,
        marker=dict(color=COLOR_P2, size=6, opacity=0.7, symbol="diamond",
                    line=dict(width=0.4, color="rgba(255,255,255,0.2)")),
        hovertemplate=f"<b>{label2}</b><br>Affinity: %{{x:.4f}}<br>Expression: %{{y:.3f}}<extra></extra>",
    ))

    fig.update_layout(
        **BASE_LAYOUT,
        title=dict(
            text=f"TF Affinity vs. Expression — {label1} vs. {label2}",
            font=dict(size=16, color=FONT_COLOR),
        ),
        xaxis_title="TF Binding Affinity",
        yaxis_title="Expression (RNA / DNA)",
        legend=dict(bgcolor="rgba(0,0,0,0)", bordercolor=GRID_COLOR, borderwidth=1),
    )
    _apply_log_y(fig, log_scale)
    return fig


# ---------------------------------------------------------------------------
# Plot C2 — Comparison: RNAP Energy vs. Expression (two promoters overlaid)
# ---------------------------------------------------------------------------
def plot_compare_rnap_energy(
    results1, colnames1, label1,
    results2, colnames2, label2,
    log_scale=True,
):
    df1 = _compute_expression(_to_df(results1, colnames1))
    df2 = _compute_expression(_to_df(results2, colnames2))

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=df1["energy"], y=df1["expression"],
        mode="markers",
        name=label1,
        marker=dict(color=COLOR_P1, size=6, opacity=0.7,
                    line=dict(width=0.4, color="rgba(255,255,255,0.2)")),
        hovertemplate=f"<b>{label1}</b><br>Energy: %{{x:.4f}}<br>Expression: %{{y:.3f}}<extra></extra>",
    ))

    fig.add_trace(go.Scatter(
        x=df2["energy"], y=df2["expression"],
        mode="markers",
        name=label2,
        marker=dict(color=COLOR_P2, size=6, opacity=0.7, symbol="diamond",
                    line=dict(width=0.4, color="rgba(255,255,255,0.2)")),
        hovertemplate=f"<b>{label2}</b><br>Energy: %{{x:.4f}}<br>Expression: %{{y:.3f}}<extra></extra>",
    ))

    fig.update_layout(
        **BASE_LAYOUT,
        title=dict(
            text=f"RNAP Energy vs. Expression — {label1} vs. {label2}",
            font=dict(size=16, color=FONT_COLOR),
        ),
        xaxis_title="RNAP Binding Energy",
        yaxis_title="Expression (RNA / DNA)",
        legend=dict(bgcolor="rgba(0,0,0,0)", bordercolor=GRID_COLOR, borderwidth=1),
    )
    _apply_log_y(fig, log_scale)
    return fig


# ---------------------------------------------------------------------------
# Plot C3 — Comparison: TF Affinity vs. RNAP Energy (two promoters overlaid)
# ---------------------------------------------------------------------------
def plot_compare_tf_rnap(
    results1, colnames1, label1,
    results2, colnames2, label2,
):
    df1 = _to_df(results1, colnames1)
    df2 = _to_df(results2, colnames2)

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=df1["affinity"], y=df1["energy"],
        mode="markers",
        name=label1,
        marker=dict(color=COLOR_P1, size=6, opacity=0.7,
                    line=dict(width=0.4, color="rgba(255,255,255,0.2)")),
        hovertemplate=f"<b>{label1}</b><br>Affinity: %{{x:.4f}}<br>Energy: %{{y:.4f}}<extra></extra>",
    ))

    fig.add_trace(go.Scatter(
        x=df2["affinity"], y=df2["energy"],
        mode="markers",
        name=label2,
        marker=dict(color=COLOR_P2, size=6, opacity=0.7, symbol="diamond",
                    line=dict(width=0.4, color="rgba(255,255,255,0.2)")),
        hovertemplate=f"<b>{label2}</b><br>Affinity: %{{x:.4f}}<br>Energy: %{{y:.4f}}<extra></extra>",
    ))

    fig.update_layout(
        **BASE_LAYOUT,
        title=dict(
            text=f"TF Affinity vs. RNAP Energy — {label1} vs. {label2}",
            font=dict(size=16, color=FONT_COLOR),
        ),
        xaxis_title="TF Binding Affinity",
        yaxis_title="RNAP Binding Energy",
        legend=dict(bgcolor="rgba(0,0,0,0)", bordercolor=GRID_COLOR, borderwidth=1),
    )
    return fig


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------
def fig_to_json(fig):
    import json
    import plotly
    return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
