import sys
import numpy as np
import plotly.graph_objects as go
import math
import dash
from dash import dcc, html
from dash.dependencies import Input, Output

#lsof -i:8050 & kill -9 to kill previous jobs

pdb_file = next((a for a in sys.argv[1:] if a.endswith('.pdb')), "1.renumber.pdb")
non_pdb_args = [a for a in sys.argv[1:] if not a.endswith('.pdb')]

# Parse alternating file label pairs: file1 label1 file2 label2 ...
if non_pdb_args and len(non_pdb_args) % 2 == 0:
    filenames = non_pdb_args[0::2]
    labels    = non_pdb_args[1::2]
elif non_pdb_args:
    filenames = non_pdb_args
    labels    = [f"Dataset {i+1}" for i in range(len(filenames))]
else:
    filenames = ["correl_01_1OCB.dat", "correl_01_1OC5.dat", "correl_01_1OCJ.dat"]  # Max 4
    labels    = ["WT", "D405N", "D416A"]

max_arcs = 2500           # Max arcs drawn per region pair
min_corr_threshold = 0.4  # Only draw arcs above/below this |value|
index_correction = 91      # Converts 0-based matrix indices to PDB residue numbers (PDB starts at 91)
label_radius = 1.16      # Radius for residue-number labels (just outside the colored arc at 1.06)
radius = 1.0              # Base circle radius

showIndices = True       # Show residue index numbers on outer ring
showColorBar = False      # Show colorbar

regions = [
    (87, 105),(223, 239),(265, 281),(313, 331)]
region_names = ['Loop1','Loop2', 'Loop3', 'Loop4']
region_colors = ["#00a6a6", "#4be2e2", "#fb5607", "#4ffb7a"]

n_regions = len(regions)
angle_per_region = 360 / n_regions

# --- Load residue names from PDB ---
residue_names = {}  # key: resSeq (int) -> three-letter name e.g. {84: 'ALA', 85: 'GLY', ...}
try:
    with open(pdb_file) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                resname = line[17:20].strip()
                resseq  = int(line[22:26].strip())
                if resseq not in residue_names:
                    residue_names[resseq] = resname
except FileNotFoundError:
    pass  # labels will fall back to number-only if PDB not found


def residue_angles(region_idx):
    start, end = regions[region_idx]
    n_res = end - start + 1
    rel_pos = np.linspace(0, 1, n_res, endpoint=False) + 0.5 / n_res
    base_angle = angle_per_region * region_idx
    angles_deg = base_angle + rel_pos * angle_per_region
    return np.deg2rad(angles_deg)


def interpolate_color(corr):
    corr = max(-1, min(1, corr))

    if corr < 0:
        t = (corr + 1)
        r = int(t * 255)
        g = int(t * 255)
        b = 255
    else:
        t = corr
        r = 255
        g = int((1 - t) * 255)
        b = int((1 - t) * 255)

    return f'rgba({r},{g},{b},0.5)'

def generate_figure(data, title, show_colorbar=False):
    fig = go.Figure()
    region_angles_arr = [residue_angles(i) for i in range(n_regions)]
    drawn_pairs = set()
    shapes = []

    # Pre-scan to find which residues have arcs, and cache arc data
    residues_with_arcs = set()
    arc_cache = {}
    for i in range(n_regions):
        start_i, end_i = regions[i]
        for j in range(i + 1, n_regions):
            start_j, end_j = regions[j]
            submatrix = data[start_i:end_i + 1, start_j:end_j + 1]
            mask = np.abs(submatrix) >= min_corr_threshold
            if not np.any(mask):
                continue
            r_idx, c_idx = np.where(mask)
            corrs = submatrix[r_idx, c_idx]
            if len(corrs) > max_arcs:
                sorted_idx = np.argsort(np.abs(corrs))[::-1][:max_arcs]
                r_idx = r_idx[sorted_idx]
                c_idx = c_idx[sorted_idx]
                corrs = corrs[sorted_idx]
            arc_cache[(i, j)] = (r_idx, c_idx, corrs)
            for ri in r_idx:
                residues_with_arcs.add(start_i + ri + index_correction)
            for ci in c_idx:
                residues_with_arcs.add(start_j + ci + index_correction)

    for i in range(n_regions):
        start_i, end_i = regions[i]
        angles_i = region_angles_arr[i]
        indices = np.arange(start_i, end_i + 1) + index_correction

        x_lbl = label_radius * radius * np.cos(angles_i)
        y_lbl = label_radius * radius * np.sin(angles_i)

        fig.add_trace(go.Scatter(
            x=x_lbl,
            y=y_lbl,
            mode='text',
            text=[f"{residue_names.get(idx, '')}{idx}" if idx in residues_with_arcs else '' for idx in indices] if showIndices else None,
            textposition='top center',
            textfont=dict(size=7, color='black', family='Arial'),
            hoverinfo='none',
            showlegend=False
        ))

        start_angle = math.radians(angle_per_region * i)
        end_angle = math.radians(angle_per_region * (i + 1))
        arc_x = 1.06 * radius * np.cos(np.linspace(start_angle, end_angle, 30))
        arc_y = 1.06 * radius * np.sin(np.linspace(start_angle, end_angle, 30))

        fig.add_trace(go.Scatter(
            x=arc_x,
            y=arc_y,
            mode='lines',
            line=dict(width=6, color=region_colors[i]),
            showlegend=False,
            hoverinfo='skip'
        ))

        mid_angle = (start_angle + end_angle) / 2
        x_txt = 1.5 * radius * math.cos(mid_angle)
        y_txt = 1.5 * radius * math.sin(mid_angle)

        fig.add_trace(go.Scatter(
            x=[x_txt],
            y=[y_txt],
            mode='text',
            text=[region_names[i]],
            textfont=dict(size=18, color=region_colors[i]),
            showlegend=False
        ))

    for i in range(n_regions):
        start_i, end_i = regions[i]
        angles_i = region_angles_arr[i]

        for j in range(i + 1, n_regions):
            if (i, j) not in arc_cache:
                continue
            start_j, end_j = regions[j]
            angles_j = region_angles_arr[j]
            r_idx, c_idx, corrs = arc_cache[(i, j)]

            x1 = radius * np.cos(angles_i[r_idx])
            y1 = radius * np.sin(angles_i[r_idx])
            x2 = radius * np.cos(angles_j[c_idx])
            y2 = radius * np.sin(angles_j[c_idx])

            for xi1, yi1, xi2, yi2, corr, ri, ci in zip(x1, y1, x2, y2, corrs, r_idx, c_idx):
                res_i = start_i + ri + index_correction
                res_j = start_j + ci + index_correction

                pair_key = tuple(sorted((res_i, res_j)))
                if pair_key in drawn_pairs:
                    continue

                drawn_pairs.add(pair_key)

                color = interpolate_color(corr)
                width = 0.5 + abs(corr) * 3.5
                path = f'M {xi1},{yi1} Q 0,0 {xi2},{yi2}'

                shapes.append(dict(
                    type='path',
                    path=path,
                    line=dict(color=color, width=width),
                    layer='below',
                    name=f'{res_i},{res_j}'
                ))

                all_arcs[filenames.index(title)].append((i, j, res_i, res_j, corr))

    fig.update_layout(
        width=480,
        height=480,
        plot_bgcolor='white',
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        margin=dict(l=10, r=70, t=60, b=20),
        shapes=shapes
    )

    if show_colorbar:
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode='markers',
            marker=dict(
                colorscale=[[0, 'blue'], [0.5, 'white'], [1, 'red']],
                cmin=-1,
                cmax=1,
                colorbar=dict(
                    title="Correlation",
                    tickvals=[-1, -0.5, 0, 0.5, 1],
                    ticktext=["-1", "-0.5", "0", "0.5", "1"],
                    len=0.6,
                    x=1.2,
                ),
                showscale=True,
                size=0
            ),
            hoverinfo='none',
            showlegend=False
        ))

    fig.update_yaxes(scaleanchor='x', scaleratio=1)
    return fig


all_arcs = [[] for _ in range(len(filenames))]

figures = [
    generate_figure(
        np.loadtxt(fname, comments='#'),
        fname,
        show_colorbar=(showColorBar and i == len(filenames) - 1)
    )
    for i, fname in enumerate(filenames)
]

# --- Dash ---
app = dash.Dash(__name__)

app.layout = html.Div([
    html.H2("Filter Residues"),
    dcc.Input(
        id='residue-input',
        type='text',
        placeholder='e.g., 45,102,200',
        style={'width': '600px'}
    ),
    html.Div([
        html.Button('Highlight', id='submit-button', n_clicks=0),
        html.Button('Table', id='table-button', n_clicks=0, style={'marginLeft': '10px'})
    ]),
    html.Div("Leave empty and press highlight to reset.", style={'marginBottom': '10px'}),

    html.Div(id='graph-container', children=[
        html.Div([
            html.H4(labels[i], style={'textAlign': 'center'}),
            dcc.Graph(id=f'graph-{i}', figure=figures[i])
        ], style={'width': '32%', 'minWidth': '480px', 'margin': '5px'})
        for i in range(len(filenames))
    ], style={'display': 'flex', 'flexWrap': 'wrap', 'justifyContent': 'center'}),

    html.Div([
        html.H4("Highlighted Correlations Table"),
        html.Pre(id='table-output', style={
            'whiteSpace': 'pre-wrap',
            'fontFamily': 'monospace',
            'border': '1px solid #ccc',
            'padding': '10px',
            'marginTop': '20px',
            'maxHeight': '300px',
            'overflowY': 'auto',
            'width': '50%'
        })
    ])
])


@app.callback(
    [Output(f'graph-{i}', 'figure') for i in range(len(filenames))],
    Input('submit-button', 'n_clicks'),
    Input('residue-input', 'value')
)
def update_all_graphs(n_clicks, input_value):
    selected = set()
    if input_value:
        parts = input_value.split(',')
        for part in parts:
            part = part.strip()
            if '-' in part:
                try:
                    start, end = map(int, part.split('-'))
                    selected.update(str(i) for i in range(start, end + 1))
                except ValueError:
                    continue
            else:
                if part.isdigit():
                    selected.add(part)

    output_figs = []
    for fig in figures:
        new_fig = go.Figure(fig)
        for shape in new_fig.layout.shapes:
            if 'name' not in shape:
                continue
            res_pair = shape['name'].split(',')
            shape.visible = (
                True if not selected or any(r in selected for r in res_pair) else False
            )
        output_figs.append(new_fig)
    return output_figs


@app.callback(
    Output('table-output', 'children'),
    Input('table-button', 'n_clicks'),
    Input('residue-input', 'value'),
    prevent_initial_call=True
)
def print_table(n_clicks_table, input_value):
    selected = set()
    if input_value:
        parts = input_value.split(',')
        for part in parts:
            part = part.strip()
            if '-' in part:
                try:
                    start, end = map(int, part.split('-'))
                    selected.update(range(start, end + 1))
                except ValueError:
                    continue
            elif part.isdigit():
                selected.add(int(part))

    output = ""
    for idx, arcs in enumerate(all_arcs):
        filtered = [
            (region1, region2, res1, res2, corr)
            for region1, region2, res1, res2, corr in arcs
            if not selected or res1 in selected or res2 in selected
        ]

        if not filtered:
            continue

        filtered.sort(key=lambda x: x[4], reverse=True)

        output += f"\n=== {labels[idx]} ===\n"
        output += f"{'Region1':<10} {'Region2':<10} {'Residue1':<10} {'Residue2':<10} {'Correlation':>10}\n"

        for region1, region2, res1, res2, corr in filtered:
            output += f"{region_names[region1]:<10} {region_names[region2]:<10} {res1:<10} {res2:<10} {corr:10.3f}\n"

    return output or "No correlations matched the filter."


if __name__ == '__main__':
    app.run(debug=True)
