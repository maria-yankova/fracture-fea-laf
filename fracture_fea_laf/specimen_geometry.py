def standard_specimen(spec_type, dimensions='2D', fraction='half', aw_ratio=0.5):
    """
    Returns the dimensions of standard fracture specimens.

    Parameters
    ----------
    spec_type : string
        Geometry of the specimen, one of 'ct-1t' | 'charpy-senb-0.5' | 'charpy-senb-0.2'. More geometries can be added.
    dimensions : string
        One of '2D' | '3D'
    fraction : string
        One of 'quarter' | 'half' | 'full'
    aw_ratio : float
        Fixed at 0.5 for CT and variable for SENB.
    
    Returns
    -------
    if 2D:
        a/w, W, A, C, D, E, F
    if 3D:h
        a/w, W, A, B, C, D, E, F
    """
    
    specimens_dims = {
    'ct-1t':{
            'a/w': 0.5,
            'W': 50,
            'A': 62.5,
            'B': 25,
            'C': 12.5,
            'D': 23.5,
            'E': 60,
            'F': 37.5
        },
    'charpy-senb-0.5':{
            'a/w': 0.5,
            'W': 10,
            'L': 55,
            'B': 10,
    },
    'charpy-senb-0.1':{
            'a/w': 0.1,
            'W': 10,
            'L': 55,
            'B': 10,
    },
      'charpy-senb-0.2':{
            'a/w': 0.2,
            'W': 10,
            'L': 55,
            'B': 10,
    },
    }
    specimen = specimens_dims[spec_type]


    if spec_type=='ct-1t':
        if fraction=='half':
            specimen['E'] *= 0.5
        elif fraction=='quarter':
            specimen['E'] *= 0.5
            specimen['B'] *= 0.5
    elif spec_type[:4]=='char':
        if fraction=='half':
            specimen['L'] *= 0.5
        elif fraction=='quarter':
            specimen['L'] *= 0.5
            specimen['B'] *= 0.5
    
    if dimensions=='2D':
        specimen.pop('B')
        
    return specimens_dims[spec_type]
    