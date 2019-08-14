'''The morph-tool command line launcher'''
from __future__ import print_function

import logging
import os
from itertools import chain
from time import time
from pathlib2 import Path

import dash
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
from neurom import COLS, NeuriteType, iter_neurites, iter_segments, load_neuron

from neuroc.axon_shrinker.shrink import cut_and_graft, cut_axon_end

L = logging.getLogger(__name__)


TREE_COLOR = {NeuriteType.basal_dendrite: 'red',
              NeuriteType.apical_dendrite: 'purple',
              NeuriteType.axon: 'blue',
              NeuriteType.soma: 'black',
              NeuriteType.undefined: 'green'}

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = html.Div(children=[

    dcc.Checklist(
        options=[
            {'label': 'Shrink', 'value': 'shrink'},
            {'label': 'Cut-end', 'value': 'cut-end'},
        ],
        value=['shrink', 'cut-end'],
        id='options',
    ),

    html.Button('Reset clicks', id='reset'),
    html.Button('Apply', id='apply'),
    html.Button('Next morphology', id='next'),
    dcc.RadioItems(
        options=[
            {'label': 'XY', 'value': 'xy'},
            {'label': 'ZY', 'value': 'zy'},
        ],
        value='xy',
        id='plane'
    ),
    html.Pre(id='dummy'),
    html.Div([], id='explanation'),
    html.Pre(id='click-data-x', hidden=True),
    html.Div(id='click-data-y', hidden=True),
    dcc.Input(
        placeholder='Filename',
        type='text',
        value='',
        id='neuron',
        style={'width': '100%'},
    ),

    html.Div(
        [
            dcc.Graph(id='graph'),
        ],
    ),

]
)


def _make_trace(neuron, plane):
    '''Create the trace to be plotted'''
    for neurite in iter_neurites(neuron):
        segments = list(iter_segments(neurite))

        segs = [(s[0][COLS.XYZ], s[1][COLS.XYZ]) for s in segments]

        coords = dict(x=list(chain.from_iterable((p1[0], p2[0], None) for p1, p2 in segs)),
                      y=list(chain.from_iterable((p1[1], p2[1], None) for p1, p2 in segs)),
                      z=list(chain.from_iterable((p1[2], p2[2], None) for p1, p2 in segs)))

        color = TREE_COLOR.get(neurite.root_node.type, 'black')
        if plane.lower() == '3d':
            plot_fun = go.Scatter3d
        else:
            plot_fun = go.Scatter
            coords = dict(x=coords[plane[0]], y=coords[plane[1]])
        yield plot_fun(
            line=dict(color=color, width=2),
            mode='lines',
            **coords
        )


def get_figure(neuron, plane, x, y):
    '''Create the figure'''
    points = [{
        'x': x or [],
        'y': y or [],
        'mode': 'markers',
        'color': 'red',
        'marker': {'size': 10}
    }]

    return {
        'data': list(_make_trace(neuron, plane)) + points,
        'layout': {
            'clickmode': 'event+select'
        }
    }


OUTPUT_FOLDER = None
INPUT_FOLDER = None
FILENAMES = None


def set_output_folder(output):
    '''Globally set the output folder'''
    global OUTPUT_FOLDER  # pylint: disable=global-statement
    OUTPUT_FOLDER = output


def set_input_folder(output):
    '''Globally set the input folder'''
    global INPUT_FOLDER, FILENAMES  # pylint: disable=global-statement
    INPUT_FOLDER = output

    def is_allowed(filename):
        s = filename.split('.')
        if len(s) < 2:
            return False
        return s[-1].lower() in {'asc', 'h5', 'swc'}
    FILENAMES = (os.path.join(INPUT_FOLDER, f) for f in os.listdir(INPUT_FOLDER) if is_allowed(f))


@app.callback(
    Output('neuron', 'value'),
    [
        Input('next', 'n_clicks')
    ],
    [State('neuron', 'value')]
)
def next_file(_, filename):
    '''Get next morphology from input folder'''
    filename = next(FILENAMES)
    return filename


@app.callback(
    Output('graph', 'figure'),
    [
        Input('click-data-x', 'children'),
        Input('click-data-y', 'children'),
        Input('neuron', 'value'),
        Input('plane', 'value')]
)
def display(x, y, filename, plane):
    '''Retrigger the graph display'''
    if not filename:
        return None
    figure = get_figure(load_neuron(filename), plane, x, y)
    figure['layout']['height'] = 1000

    return figure


@app.callback(
    Output('explanation', 'children'),
    [Input('options', 'value')]
)
def explain(options):
    '''Explain what to do depending on the configuration of the 'cut-end' and 'shrink' checkboxes'''
    if not options:
        return 'Please activate at least one of: shrink, cut end'

    if len(options) == 2:
        return ('Click 3 times and press "Apply". The first 2 clicks will define the region to be'
                ' shrinked (ie. removed), the 3rd one will define the start of piece of axon to be'
                ' thrown away')

    if options[0] == 'shrink':
        return 'Click twice to define the region to be shrinked (ie. removed),'

    if options[0] == 'cut-end':
        return 'Click once to indicate of piece of axon to be thrown away'
    return 'Should never be visible'


@app.callback(
    Output('dummy', 'children'),
    [Input('apply', 'n_clicks')],
    [State('click-data-x', 'children'),
     State('click-data-y', 'children'),
     State('neuron', 'value'),
     State('options', 'value')]
)
def apply(n_clicks, _, y, filename, options):
    '''Produces the shrinked neuron'''
    if not n_clicks:
        return ''
    path = Path(filename)
    name_no_ext = path.stem
    ext = path.suffix
    data = y
    if options == ['cut-end']:
        neuron = cut_axon_end(filename, data[0])
        path = os.path.join(OUTPUT_FOLDER, '{}-end-cut.{}'.format(name_no_ext, ext))
        neuron.write(path)
        return path

    if options == ['shrink'] and data and len(data) >= 2:
        upward = (data[1] > data[0])
        cut_neuron, _ = cut_and_graft(filename, upward, data[0], data[1], 0)
        path = os.path.join(OUTPUT_FOLDER, '{}-shrinked.{}'.format(name_no_ext, ext))
        cut_neuron.write(path)
        return path

    if set(options) == {'shrink', 'cut-end'} and data and len(data) >= 3:
        upward = (data[1] > data[0])
        shrinked_neuron, _ = cut_and_graft(filename, upward, data[0], data[1], 0)
        neuron = cut_axon_end(shrinked_neuron, data[2])
        path = os.path.join(OUTPUT_FOLDER, '{}-shrinked-and-end-end.{}'.format(name_no_ext, ext))
        neuron.write(path)
        return path

    return filename


@app.callback(
    Output('click-data-y', 'children'),
    [Input('graph', 'clickData'),
     Input('reset', 'n_clicks_timestamp')],
    [State('click-data-y', 'children')]
)
def display_click_data_y(clickData, reset_time, pos_s):
    '''Display click data'''
    if not clickData:
        return None

    if reset_time:
        new_click = (time() - (reset_time / 1000.)) < 1
        if new_click:
            return []

    pos = clickData['points'][0]
    pos = pos['y']
    if not pos_s:
        return [pos]
    elif len(pos_s) == 3:
        return pos_s
    else:
        return pos_s + [pos]


@app.callback(
    Output('click-data-x', 'children'),
    [Input('graph', 'clickData'),
     Input('reset', 'n_clicks_timestamp'),
     Input('next', 'n_clicks_timestamp')],
    [State('click-data-x', 'children')]
)
def display_click_data_x(clickData, reset_time, next_time, pos_s):
    '''Display click data'''
    if not clickData:
        return None

    if next_time:
        new_click = (time() - (next_time / 1000.)) < 0.1
        if new_click:
            return []

    if reset_time:
        new_click = (time() - (reset_time / 1000.)) < 0.1
        if new_click:
            return []

    pos = clickData['points'][0]
    pos = pos['x']
    if not pos_s:
        return [pos]
    elif len(pos_s) == 3:
        return pos_s
    else:
        return pos_s + [pos]
