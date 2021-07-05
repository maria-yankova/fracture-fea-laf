from plotly import graph_objs as go

def plot_nodes(plot_x, plot_y, plot_z=None, indices=None, xrange=None, yrange=None,
               figsize=[400, 400], dticks=None, mode='markers+text',
              markersize=5, showlegend=False, show_labels=True, threed=False):
    """
    Parameters
    ----------
    show_labels : bool, optional
        If True, show node labels.
        
    """
    labels = ['node labels', 'element labels']
    labels = len(plot_x) * ['h']
    data = []
    if threed:
        if type(plot_x) == list:
            i = 0
            for plx, ply, plz in zip(plot_x, plot_y, plot_z):
                if indices == None:
                    text = None
                else:
                    text = indices[i]
                trace = {
                    'type': 'scatter3d',
                    'x': plx,
                    'y': ply,
                    'z': plz,
                     'mode': mode,
                    'marker': {'size': markersize},
                     'name': labels[i],                
                }
                if show_labels:
                    trace.update({'text': text, 'type': 'scatter3d'})

                data.append(trace)

                i += 1
        else:    
            if indices is not None:
                indices = indices
                print('here')
                print('indices shape: ', indices.shape)
                print('plot_x shape: ', plot_x.shape)
                data.append({
                     'type': 'scatter3d',
                        'x': plot_x,
                        'y': plot_y,
                        'z': plot_z,
                         'mode': mode,
                        'marker': {'size': 5},
                        'text': indices,
                        'name': labels[0],
                    })
            else:
                data.append({
                     'type': 'scatter3d',
                        'x': plot_x,
                        'y': plot_y,
                         'z': plot_z,
                        'text': None,
                         'mode': 'markers',
                        'marker': {'size': 5},
                        'name': labels[0],
                    })
    else:
        if type(plot_x) == list:
            i = 0
            for plx, ply in zip(plot_x,plot_y):
                if indices == None:
                    text = None
                else:
                    text = indices[i]
                trace = {
                    'type': 'scattergl',
                    'x': plx,
                    'y': ply,
                     'mode': mode,
                    'marker': {'size': markersize},
                    'name': labels[i],           
                }
                if show_labels:
                    trace.update({'text': text, 'type': 'scatter'})

                data.append(trace)

                i += 1
        else:    
            if indices is not None:
                indices = indices
                data.append({
                        'x': plot_x,
                        'y': plot_y,
                         'mode': mode,
                        'marker': {'size': 5},
                        'text': indices,
                         'name': labels[0],
                    })
            else:
                data.append({
                        'x': plot_x,
                        'y': plot_y,
                        'text': None,
                         'mode': 'markers',
                        'marker': {'size': 5},
                         'name': labels[0],
                    })

    layout = {
        'width': figsize[0],
        'height': figsize[1],
        'xaxis': {
            'zeroline': False,
            'showline': True,
            'mirror': 'all',
            'title': r'x (mm)',
            'titlefont':{'size': 13},
        },
        'yaxis': {
            'zeroline': False,
            'showline': True,
            'mirror': 'all',
            'title': r'y (mm)',
            'titlefont':{'size': 13},
            'scaleanchor': 'x',
            'scaleratio': 1,
        },
        # 'annotations':[
        # {
        #     'text':  r'$\LARGE{x}$',
        #     'x': 0.5,
        #     'y': -0.1,
        #     'xref': 'paper',
        #     'yref': 'paper',
        #     'showarrow': False,
        #     'font': {'size': 11},
        # },
        #     {
        #     'text':  r'$\LARGE{y}$',
        #     'x': -0.1,
        #     'y': 0.5,
        #     'xref': 'paper',
        #     'yref': 'paper',
        #     'showarrow': False,
        #     'font': {'size': 11},
        #     'textangle': -90
        # },
        # ],
        'legend':{
            'x': 0.6,
            'y': 0.7,
        },
         'margin': {
                    'l': 100,
                    'r': 20,
                    'b': 80,
                    't': 40,
         },

    }
    if xrange:
        layout['xaxis']['range'] = xrange
    if yrange:
        layout['yaxis']['range'] = yrange
    if dticks:
        layout['xaxis']['dtick'] = dticks[0]
        layout['yaxis']['dtick'] = dticks[1]
    layout['showlegend'] = showlegend
    f = go.FigureWidget(data, layout)
    return f
