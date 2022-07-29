def standard_specimen(spec_type, dimensions='2D', fraction='half', aw_ratio=0.5):
    """
    Returns the dimensions of standard fracture specimens.

    Parameters
    ----------
    spec_type : string
        Geometry of the specimen, one of 'ct-1t' | 'ct-1t-jacobs' | 'ct-2t' |'charpy-senb-0.5' | 'charpy-senb-0.2'. More geometries can be added, where the beginning og the geometry name should be kept the same, e.g. 'ct' for compact tension. Compact tension (ct) specimen geometry parameters as in Fig. 5 (a) [1]. Single-edge notched bend (senb) specimen geometry parameters as follows: length 'L', width 'W' and thickness 'B'. In both geometries the 'a/w' parameter refers to the ratio of the crack length to specimen width. 
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

    References
    ----------
    [1] Heerens J, Hellmann D. Development of the Euro fracture toughness dataset. Eng Fract Mech 2002;69:421–49. https://doi.org/10.1016/S0013-7944(01)00067-4.
    """
    specimens_dims = {
    'ct-0.5t':{
            'a/w': 0.5,
            'W': 25,
            'A': 31.25,
            'B': 12.5,
            'C': 6.25,
            'D': 10,
            'E': 30,
            'F': 18.75
        },
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
    'ct-1t-prime':{
        'a/w': 0.5,
            'W': 50,
            'A': 62.5,
            'B': 25,
            'C': 12.5,
            'D': 22.5,
            'E': 60,
            'F': 35.5
        },
    'ct-2t':{
            'a/w': 0.5,
            'W': 100,
            'A': 125.0,
            'B': 50,
            'C': 25.0,
            'D': 51.0,
            'E': 120,
            'F': 75.0,
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


    if spec_type[:2]=='ct':
        if fraction=='half':
            specimen['E'] *= 0.5
        elif fraction=='quarter':
            specimen['E'] *= 0.5
            specimen['B'] *= 0.5
    elif spec_type[:11]=='charpy-senb':
        if fraction=='half':
            specimen['L'] *= 0.5
        elif fraction=='quarter':
            specimen['L'] *= 0.5
            specimen['B'] *= 0.5
    
    if dimensions=='2D':
        specimen.pop('B')
        
    return specimens_dims[spec_type]
    