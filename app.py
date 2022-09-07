import pandas as pd
import numpy as np
import plotly.express as px
import plotly.figure_factory as ff
import dash
from dash import html
from dash import dcc
from dash.dependencies import Input, Output, State
from dash import dash_table
import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.colors import n_colors

#counts_url = "https://www.dropbox.com/s/pivx6ktd4zfj4jn/counts.csv?dl=1"
#seurat_counts = pd.read_table(counts_url,sep=',', header=(0))

TSNEdf = pd.read_csv('TSNEdf.csv')
TSNEdf['clusters'] = TSNEdf['clusters'].astype(str)
UMAPdf = pd.read_csv('UMAPdf.csv')
UMAPdf['clusters'] = UMAPdf['clusters'].astype(str)

app = dash.Dash(__name__,external_stylesheets=[dbc.themes.LUX])
server = app.server

app.layout = html.Div([
    dbc.Tabs(
        [dbc.Tab(label="UMAP/tSNE", tab_id="umap"),
         dbc.Tab(label="Clustering", tab_id="clustering"),
         dbc.Tab(label="Gene Search", tab_id="gene_search"),
         dbc.Tab(label="Annotation", tab_id="annotation"),
         dbc.Tab(label="Table and Code", tab_id="table")],
        id="tabs",
        active_tab="umap",
    ),
    html.Div(id="tab-content", className="p-4"),
])

@app.callback(
    Output("tab-content", "children"),
    [Input("tabs", "active_tab")])

def render_tab_content(active_tab):
    if active_tab == "umap":
        return html.Div([
            dbc.Row([
                dbc.Col([
                    dbc.RadioItems(
                        options=[
                            {"label": "UMAP", "value": 'UMAP'},
                            {"label": "tSNE", "value": 'tSNE'}],
                        value='UMAP',
                        id="radioitems-input_1",
                        inline=True
                    ),
                ],width={'size':2,"offset":0,'order':1},style={'text-align':'center'}), #,style={'padding-top' : 10}
                dbc.Col([
                    dbc.RadioItems(
                        options=[
                            {"label": "2D", "value": '2D'},
                            {"label": "3D", "value": '3D'}],
                        value='2D',
                        id="radioitems-input_2",
                        inline=True
                    ),
                ],width={'size':2,"offset":0,'order':1},style={'text-align':'center'}), #,style={'padding-top' : 10}
                dbc.Col([
                    dbc.Row([
                        dbc.Col([
                            html.H6('X_AXIS',style={'padding-top' : 10,'padding-right' : 2})
                        ],width={'size':4,"offset":0,'order':1}),
                        dbc.Col([
                            dcc.Dropdown(id="x_axis",
                                         options=[],
                                         value=[],
                                         multi=False,
                                         disabled=False,
                                         clearable=False,
                                         searchable=True)
                        ],width={'size':8,"offset":0,'order':1},style={'text-align':'center'})
                    ])
                ],width={'size':2,"offset":0,'order':1},style={'text-align':'center'}),
                dbc.Col([
                    dbc.Row([
                        dbc.Col([
                            html.H6('Y_AXIS',style={'padding-top' : 10,'padding-right' : 2})
                        ],width={'size':4,"offset":0,'order':1}),
                        dbc.Col([
                            dcc.Dropdown(id="y_axis",
                                         options=[],
                                         value=[],
                                         multi=False,
                                         disabled=False,
                                         clearable=False,
                                         searchable=True)
                        ],width={'size':8,"offset":0,'order':1},style={'text-align':'center'})
                    ])
                ],width={'size':2,"offset":0,'order':1},style={'text-align':'center'}),
                dbc.Col([
                    dbc.Row([
                        dbc.Col([
                            html.H6('Z_AXIS',style={'padding-top' : 10,'padding-right' : 2})
                        ],width={'size':4,"offset":0,'order':1}),
                        dbc.Col([
                            dcc.Dropdown(id="z_axis",
                                         options=[],
                                         value=[],
                                         multi=False,
                                         disabled=False,
                                         clearable=False,
                                         searchable=True)
                        ],width={'size':8,"offset":0,'order':1},style={'text-align':'center'})
                    ])
                ],width={'size':2,"offset":0,'order':1},style={'text-align':'center'}),
            ], className='p-2 align-items-stretch'),

            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                             dbc.Row([
                                dbc.Col([
                                    dbc.Button("SVG", size="sm",className="me-1",id='btn_1',color="secondary"),
                                    dcc.Download(id='download_1'),
                                    dbc.Button("HTML", size="sm",className="me-1",id='btn_2',color="secondary"),
                                    dcc.Download(id='download_2'),
                                    dbc.Button("csv", size="sm",className="me-1",id='btn_3',color="secondary"),
                                    dcc.Download(id='download_3')
                                ],width={'size':12,'offset':0,'order':1},style={'text-align':'right'}),
                            ]),

                            dbc.Row([
                                dbc.Col([
                                    html.Div(id='chart_title'),
                                ],width={'size':12,'offset':0,'order':1},style={'text-align':'center'}),
                            ]),
                            dbc.Row([
                                dbc.Col([
                                     dcc.Graph(id='scatter_chart',figure={},style={'height': '450px'},selectedData={'points': [{'hovertext':'X24_CTGACACAATGC'}]})
                                ],width={'size':12,'offset':0,'order':1}),
                            ]),
                        ])
                    ], className='h-100 text-left')
                ], xs=6),
                dbc.Col([
                    dbc.Row([
                        dbc.Col([
                            dbc.Card([
                                dbc.CardBody([
                                    dbc.Row([
                                        dbc.Col([
                                            dbc.Button("SVG", size="sm",className="me-1",id='btn_4',color="secondary"),
                                            dcc.Download(id='download_4'),
                                            dbc.Button("HTML", size="sm",className="me-1",id='btn_5',color="secondary"),
                                            dcc.Download(id='download_5')
                                        ],width={'size':12,'offset':0,'order':1},style={'text-align':'right'}),
                                    ]),
                                    dbc.Row([
                                        dbc.Col([
                                            html.Span(f"Sample Violin Plot for Samples",style={'text-align':'center'}),
                                            dcc.Graph(id='violin_chart',figure={},style={'height': '200px'})
                                        ],width={'size':12,'offset':0,'order':1},style={'text-align':'center'}),
                                    ]),
                                ])
                            ], className='h-100 text-left')
                        ])
                    ]),
                    dbc.Row([
                        dbc.Col([
                            dbc.Card([
                                dbc.CardBody([
                                    dbc.Row([
                                        dbc.Col([
                                            dbc.Button("SVG", size="sm",className="me-1",id='btn_6',color="secondary"),
                                            dcc.Download(id='download_6'),
                                            dbc.Button("HTML", size="sm",className="me-1",id='btn_7',color="secondary"),
                                            dcc.Download(id='download_7')
                                        ],width={'size':12,'offset':0,'order':1},style={'text-align':'right'}),
                                    ]),
                                    dbc.Row([
                                        dbc.Col([
                                            html.Div(id='chart_title_2'),
                                        ],width={'size':12,'offset':0,'order':1},style={'text-align':'center'}),
                                    ]),
                                    dbc.Row([
                                        dbc.Col([
                                            html.H6(f'Gene',style={'padding-top' : 10,'padding-right' : 2}),
                                        ],width={'size':3,'offset':0,'order':1},style={'text-align':'center'}),
                                        dbc.Col([
                                            dcc.Dropdown(id="gene_2",
                                                         options=[],
                                                         value=[],
                                                         multi=False,
                                                         disabled=False,
                                                         clearable=False,
                                                         searchable=True)
                                        ],width={'size':6,'offset':0,'order':1},style={'text-align':'center'}),
                                    ]),

                                    dbc.Row([
                                        dbc.Col([
                                            dcc.Loading(children=[dcc.Graph(id='violin_chart_2',figure={},style={'height': '200px'})],color='#119DFF',type='dot'),
                                        ],width={'size':12,'offset':0,'order':1},style={'text-align':'center'}),
                                    ]),

                                ])
                            ], className='h-100 text-left')
                        ])
                    ])
                ], xs=6),
            ], className='p-2 align-items-stretch'),
        ])

    elif active_tab == "clustering":
        return html.Div([
            dbc.Row([
                dbc.Col([
                    dbc.Row([
                        dbc.Col([
                            html.H6('Cluster', style={'padding-top': 10, 'padding-right': 2}),
                            dcc.Input(type="number", id="num_of_cluster", value=4, min=0,
                                      style={'width': '100%', 'text-align': 'center'})
                        ], width={'size': 3, "offset": 0, 'order': 1},
                            style={'text-align': 'center', 'display': 'flex'}),  # ,style={'padding-top' : 10}
                        dbc.Col([
                            html.H6('ID Numbers', style={'padding-top': 10, 'padding-right': 2}),
                        ], width={'size': 3, "offset": 0, 'order': 1}, style={'text-align': 'center'}),
                        dbc.Col([
                            dcc.Input(type="number", id="num_of_id", value=0, min=0,
                                      style={'width': '40%', 'text-align': 'center'}),
                            dcc.Input(type="number", id="num_of_id_2", value=15, min=0,
                                      style={'width': '40%', 'text-align': 'center'})
                        ], width={'size': 3, "offset": 0, 'order': 1},
                            style={'text-align': 'center', 'display': 'flex'}),  # ,style={'padding-top' : 10}
                        dbc.Col([
                            dbc.Button("Submit", size="sm", className="me-1", id='submit', color="primary")
                        ], width={'size': 3, "offset": 0, 'order': 1}, style={'text-align': 'right'}),
                        # ,style={'padding-top' : 10}
                    ])
                ], width={'size': 6, "offset": 0, 'order': 1}, style={'text-align': 'center'}),
                dbc.Col([
                    dbc.Row([
                        dbc.Col([
                            html.H6('GENE:', style={'padding-top': 10, 'padding-right': 2})
                        ], width={'size': 3, "offset": 0, 'order': 1}),
                        dbc.Col([
                            dcc.Dropdown(id="gene",
                                         options=[],
                                         value=[],
                                         multi=False,
                                         disabled=False,
                                         clearable=False,
                                         searchable=True)
                        ], width={'size': 9, "offset": 0, 'order': 1}, style={'text-align': 'center'})
                    ])
                ], width={'size': 6, "offset": 0, 'order': 1}, style={'text-align': 'center'}),
            ], className='p-2 align-items-stretch'),

            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            dbc.Row([
                                dbc.Col([
                                    dbc.Button("SVG", size="sm", className="me-1", id='btn_8', color="secondary"),
                                    dcc.Download(id='download_8'),
                                    dbc.Button("HTML", size="sm", className="me-1", id='btn_9', color="secondary"),
                                    dcc.Download(id='download_9'),
                                    dbc.Button("csv", size="sm", className="me-1", id='btn_10', color="secondary"),
                                    dcc.Download(id='download_10')
                                ], width={'size': 12, 'offset': 0, 'order': 1}, style={'text-align': 'right'}),
                            ]),
                            dbc.Row([
                                dbc.Col([
                                    html.Span(f"Sample Violin Plot for Samples", style={'text-align': 'center'}),
                                    dcc.Loading(
                                        children=[dcc.Graph(id='violin_chart_3', figure={}, style={'height': '350px'})],
                                        color='#119DFF', type='dot'),
                                ], width={'size': 12, 'offset': 0, 'order': 1}, style={'text-align': 'center'}),
                            ]),
                        ])
                    ], className='h-100 text-left')
                ], xs=6),
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            dbc.Row([
                                dbc.Col([
                                    dbc.Button("SVG", size="sm", className="me-1", id='btn_11', color="secondary"),
                                    dcc.Download(id='download_11'),
                                    dbc.Button("HTML", size="sm", className="me-1", id='btn_12', color="secondary"),
                                    dcc.Download(id='download_12')
                                ], width={'size': 12, 'offset': 0, 'order': 1}, style={'text-align': 'right'}),
                            ]),
                            dbc.Row([
                                dbc.Col([
                                    html.Div(id='chart_title_3'),
                                ], width={'size': 12, 'offset': 0, 'order': 1}, style={'text-align': 'center'}),
                            ]),
                            dbc.Row([
                                dbc.Col([
                                    dcc.Loading(
                                        children=[dcc.Graph(id='violin_chart_4', figure={}, style={'height': '350px'})],
                                        color='#119DFF', type='dot'),
                                ], width={'size': 12, 'offset': 0, 'order': 1}),
                            ]),
                        ])
                    ], className='h-100 text-left')
                ], xs=6),
            ], className='p-2 align-items-stretch'),

            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            dbc.Row([
                                dbc.Col([
                                    dbc.Button("SVG", size="sm", className="me-1", id='btn_13', color="secondary"),
                                    dcc.Download(id='download_13'),
                                    dbc.Button("HTML", size="sm", className="me-1", id='btn_14', color="secondary"),
                                    dcc.Download(id='download_14')
                                ], width={'size': 12, 'offset': 0, 'order': 1}, style={'text-align': 'right'}),
                            ]),
                            dbc.Row([
                                dbc.Col([
                                    html.Span(f"Number of Cells in Each Cluster", style={'text-align': 'center'}),
                                    dcc.Loading(
                                        children=[dcc.Graph(id='bar_chart', figure={}, style={'height': '300px'})],
                                        color='#119DFF', type='dot'),
                                ], width={'size': 12, 'offset': 0, 'order': 1}, style={'text-align': 'center'}),
                            ]),
                        ])
                    ], className='h-100 text-left')
                ], xs=12),
            ], className='p-2 align-items-stretch'),

        ])



@app.callback([Output('x_axis', 'options'),
               Output('x_axis', 'value'),
               Output('y_axis', 'options'),
               Output('y_axis', 'value'),
               Output('z_axis', 'options'),
               Output('z_axis', 'value')],
              [Input('radioitems-input_1', 'value')])

def update_axis(radioitems_input_1):
    if radioitems_input_1 == 'UMAP':
        value = [{'label': 'UMAP1', 'value': 'UMAP1'},
                 {'label': 'UMAP2', 'value': 'UMAP2'},
                 {'label': 'UMAP3', 'value': 'UMAP3'}]
        x = 'UMAP1'
        y = 'UMAP2'
        z = 'UMAP3'
        return value, x, value, y, value, z
    elif radioitems_input_1 == 'tSNE':
        value = [{'label': 'TSNE1', 'value': 'TSNE1'},
                 {'label': 'TSNE2', 'value': 'TSNE2'},
                 {'label': 'TSNE3', 'value': 'TSNE3'}]
        x = 'TSNE1'
        y = 'TSNE2'
        z = 'TSNE3'
        return value, x, value, y, value, z


@app.callback([Output('scatter_chart', 'figure'),
               Output('chart_title', 'children')],
              [Input('radioitems-input_1', 'value'),
               Input('radioitems-input_2', 'value'),
               Input('x_axis', 'value'),
               Input('y_axis', 'value'),
               Input('z_axis', 'value')
               ])
def update_scatter_chart(radioitems_input_1, radioitems_input_2, x_axis, y_axis, z_axis):
    if radioitems_input_1 == 'tSNE' and radioitems_input_2 == '2D':
        figTSNE = px.scatter(TSNEdf, x=TSNEdf[x_axis], y=TSNEdf[y_axis], color='clusters',
                             labels='Id',
                             hover_name='Id')
        figTSNE.update_layout(template='plotly_white',
                              margin=dict(l=0, r=0, t=0, b=0),
                              clickmode='event+select',
                              # dragmode='select',
                              hovermode='closest')
        # figTSNE.update_traces(marker_size = 3)
        return figTSNE, html.Span(f'2D TSNE', style={'text-align': 'center'})

    elif radioitems_input_1 == 'tSNE' and radioitems_input_2 == '3D':
        figTSNE3D = px.scatter_3d(TSNEdf, x=TSNEdf[x_axis], y=TSNEdf[y_axis], z=TSNEdf[z_axis],
                                  color='clusters',
                                  hover_name='Id')
        figTSNE3D.update_layout(template='plotly_white',
                                legend_traceorder="normal",
                                margin=dict(l=0, r=0, t=0, b=0),
                                clickmode='event+select',
                                # dragmode='select',
                                hovermode='closest')
        figTSNE3D.update_traces(marker_size=2)
        return figTSNE3D, html.Span(f'3D TSNE', style={'text-align': 'center'})

    elif radioitems_input_1 == 'UMAP' and radioitems_input_2 == '2D':
        figUMAP = px.scatter(UMAPdf, x=UMAPdf[x_axis], y=UMAPdf[y_axis], color='clusters',
                             labels='Id',
                             hover_name='Id')
        figUMAP.update_layout(template='plotly_white',
                              margin=dict(l=0, r=0, t=0, b=0),
                              clickmode='event+select',
                              # dragmode='select',
                              hovermode='closest')
        # figUMAP.update_traces(marker_size = 3)
        return figUMAP, html.Span(f'2D UMAP', style={'text-align': 'center'})

    elif radioitems_input_1 == 'UMAP' and radioitems_input_2 == '3D':
        figUMAP3D = px.scatter_3d(UMAPdf, x=UMAPdf[x_axis], y=UMAPdf[y_axis], z=UMAPdf[z_axis],
                                  color='clusters',
                                  hover_name='Id')
        figUMAP3D.update_layout(template='plotly_white',
                                legend_traceorder="normal",
                                margin=dict(l=0, r=0, t=0, b=0),
                                clickmode='event+select',
                                # dragmode='select',
                                hovermode='closest')
        figUMAP3D.update_traces(marker_size=2)
        return figUMAP3D, html.Span(f'3D UMAP', style={'text-align': 'center'})

@app.callback(Output('bar_chart', 'figure'),
             [Input('gene','value')])

def update_violin_chart(gene):
    Cell_count_clus = UMAPdf[['Id','clusters']].groupby(['clusters'],as_index = False).count()
    bar_chart = px.bar(Cell_count_clus, x='clusters', y='Id').update_xaxes(categoryorder="total descending")
    bar_chart.update_layout(template='plotly_white',margin=dict(l=0,r=0,t=0,b=0))
    return bar_chart


if __name__ == '__main__':
    app.run_server(debug=True,host="0.0.0.0",port=8080)