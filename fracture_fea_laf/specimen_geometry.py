def standard_specimen(spec_type, dimensions='2D', fraction='half', aw_ratio=0.5):
    """
    Returns the dimensions of standard fracture specimens.

    Parameters
    ----------
    spec_type
        CT-1T/SENB
    dimensions
        2D/3D
    fraction
        quarter/half/full
    aw_ratio
        fix at 0.5 for CT
    
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
    'senb-0.5':{
            'a/w': 0.5,
            'W': 10,
            'L': 55,
            'B': 10,
    },
    'senb-0.1':{
            'a/w': 0.1,
            'W': 10,
            'L': 55,
            'B': 10,
    },
      'senb-0.2':{
            'a/w': 0.2,
            'W': 10,
            'L': 55,
            'B': 10,
    },
    }
    specimen = specimens_dims[spec_type]

    if dimensions=='2D':
        specimen.pop('B')
        
    if spec_type=='ct-1t':
        if fraction=='half':
            specimen['E'] *= 0.5
        elif fraction=='quarter':
            specimen['E'] *= 0.5
            specimen['B'] *= 0.5
    elif spec_type[:4]=='senb':

        if fraction=='half':
            specimen['L'] *= 0.5
        elif fraction=='quarter':
            specimen['L'] *= 0.5
            specimen['B'] *= 0.5

    return specimens_dims[spec_type]
    