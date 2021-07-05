import numpy as np
from more_itertools import pairwise
from fracture_fea_laf import utils

def triangle_area(tri):
    x1, y1, x2, y2, x3, y3 = tri[0][0], tri[0][1], tri[1][0], tri[1][1], tri[2][0], tri[2][1]
    return abs(0.5 * (((x2-x1)*(y3-y1))-((x3-x1)*(y2-y1))))

def donut_polar_grid_quad(inner_r, outer_r, num_rays, num_arcs, quad=1,     include_borders=True):
    """
    Create a donut polar grid of nodes for one quadrant.
    
    Parameters
    ----------
    inner_r : float
        Size of inner radius of donut
    outer_r : float
        Size of outer radius of donut
    num_rays : int
        Number of rays within a quadrant of the circle (splits angles)
    num_arcs : int
        Number of arcs within a quadrant of the circle (splits radius)
    quad : int
        Circle quadrant to be populated, accepted values are 1, 2, 3 or 4 (default is 1)
    include_borders : bool
        Include border nodes. Default value is True.    
    
    Outputs
    -------
    tuple (ndarray, ndarray), shape=(N,)
        Polar coordinates distance-angle for nodes.
        
    
    """
    if quad == 1:
        ang_lim_low = 0.0
        ang_lim_upp = np.pi / 2 + 0.1
        
    r = np.linspace(inner_r, outer_r, num_arcs )
    theta = np.linspace(0, 2 * np.pi, (num_rays - 1) * 4 +1)

    radius_matrix, theta_matrix = np.meshgrid(r,theta)
    quat_idx = np.intersect1d(np.where(theta_matrix.flatten() >= ang_lim_low)[0],
               np.where(theta_matrix.flatten() <= ang_lim_upp)[0])
    
    
    return theta_matrix.flatten()[quat_idx], radius_matrix.flatten()[quat_idx]


def order_nodes_elem(elem_nodes_all, all_nodes, all_nodes_labs):
    """
    Order the nodes of an element in anticlockwise direction, consistent with the element definition in Abaqus.
    
    Parameters:
    elem_nodes_all : array of (M, n)
        The nodes that form part of M elements with n corners each.
    all_nodes : array of (N, 2)
        The coordinates of all N nodes in the mesh.
    all_nodes_labs : array of (N,)
        The labels of all nodes in the mesh.
    
    """
    elem_nodes_all_ordered = []
    for el_idx, elem_nodes in enumerate(elem_nodes_all):
        
        
        # # sort nodes in growing order
        # node_el_sorted = np.sort(elem_nodes)
        
        # indices of nodes in node list
        node_idx = np.where(np.isin(all_nodes_labs, elem_nodes))[0]
       
        if len(node_idx)==3:
            vals, counts = np.unique(elem_nodes, return_counts=True)
            rep_idx = np.where(counts == 2)[0]
            node_idx = np.insert(node_idx, rep_idx, node_idx[np.where(counts == 2)])

            
        node_el_sorted_by_node_list = all_nodes_labs[node_idx]    
        # coords of nodes in the order they are in the all nodes list?
        node_coord_sorted = all_nodes[node_idx]
        # order nodes in correct order
        
        node_ord = utils.order_coplanar_points(
            np.hstack((node_coord_sorted, np.zeros((len(node_coord_sorted),1)))).T,
            np.array([[0,0,1]]).T
        )
        
        elem_nodes_all_ordered.append(node_el_sorted_by_node_list[node_ord])
        
    return np.concatenate(elem_nodes_all_ordered).reshape(-1, 4)



def cells_from_nodes(nodes_array):
    """
    Find the cell vertices in a regular grid defined by an array of nodes
    
    """
    c = 0
    cells = []

    for ri, row in enumerate(nodes_array):
        count = 0
        if ri != len(nodes_array)-1:
            cells.append([])
        for item1, item2 in pairwise(row):
            if ri == 0:
                cells[ri].append([item1, item2])
            elif ri == len(nodes_array)-1:
                cells[ri-1][count].append(item1)
                cells[ri-1][count].append(item2)
    
            else:
                cells[ri].append([item1, item2])
                cells[ri-1][count].append(item1)
                cells[ri-1][count].append(item2)
            count += 1
    
    return np.array(cells)


def make_donut_mesh(crack_tip_radius_microns, fine_mesh_element_length=0.025, fan_box_num_side_elements=3, fan_box_width=5):
    """
    
    Parameters
    ----------
    crack_tip_radius_microns : float
        Inner radius in micrometers
    
    fan_box_num_side_elements : integer
     number of arcs (TO CALCULATE LATER) = 3 
    fan_box_width : integer
       aproximate width of fan mesh in multiples of r 
    ret_rin_nodes : bool
        Return the inner nodes compromising the inner radius in anti-clockwise direction. 
       
    Returns
    -------
    
    
    """    
    rad_mm = crack_tip_radius_microns * 1e-3
    
    # width of fan mesh in multiples of elem side
    fbox_side_mult_elem = int((fan_box_width*rad_mm) / fine_mesh_element_length)
    fbox_side_mult_elem = round(fbox_side_mult_elem/3)*3
    
    # side size of fan mesh in mm
    fbox_side_mm = fbox_side_mult_elem * fine_mesh_element_length
    
    # number of rays
    rn = 2 * fbox_side_mult_elem + 1
    # rn += 1
    thetaout, rout  = donut_polar_grid_quad(rad_mm, fbox_side_mm,
                                            rn, fan_box_num_side_elements)
    
    xout, yout = utils.polar2cart_2D(rout, thetaout)
    
    # Move nodes at ray ends to the sides of a square
    ##################################################
    # coords from 0 to fan_box_num_side_elements of side square box in multiples of ctip radii
    # half coordinates of the end of fan rays crossing square
    end_ray_half_coords =  np.linspace(0, fbox_side_mm, int((rn-1)/2+1))
    # print('end_ray_half_coords: ', end_ray_half_coords)

    end_ray_repeat_coords = (fbox_side_mm * np.ones(rn))[None]
    end_ray_coords = np.concatenate((
        np.concatenate((end_ray_half_coords[:-1],
                    np.flip(end_ray_half_coords)))[None],
            end_ray_repeat_coords
        ), axis=0)

    # The second half of array indices: 
    # 1st half + mid element + next element - 1 for zero based
    # int((rn-1)/2) + 1 + 1 - 1 
    flip_idx = int((rn-1)/2) + 1 - 1
    end_ray_coords[:, :flip_idx] = np.flip(end_ray_coords[:, :flip_idx], axis=0)

    # Array of nodes as (Rows, Cols, (x, y)_node)
    print('rn: ', rn)
    # change from rn to rn + 1
    don_reshaped = np.concatenate((xout[None].T, yout[None].T), axis=1).reshape(rn, fan_box_num_side_elements, 2)
    don_reshaped[:,-1,:] = end_ray_coords.T
    don_nodes_flattened = don_reshaped.reshape(rn*fan_box_num_side_elements, 2)
    don_node_labs = np.arange(1, fan_box_num_side_elements * rn + 1).reshape(rn, fan_box_num_side_elements)
    don_node_labs_flattened = don_node_labs.flatten()
    
    don_cells = cells_from_nodes(don_node_labs)
    don_cells_flattened = don_cells.reshape(-1, 4)
    

    don_cells_flattened = order_nodes_elem(don_cells_flattened, don_nodes_flattened,don_node_labs_flattened)

    don_cell_centres = np.mean(
        cells_from_nodes(don_reshaped), axis=2
    )
    don_cell_centres_flattened = don_cell_centres.reshape((fan_box_num_side_elements-1) * (rn-1), 2)

    don_cell_labs = np.arange(1, (fan_box_num_side_elements-1) * (rn-1) + 1).reshape(fan_box_num_side_elements-1, rn-1)
    don_cell_labs_flattened = don_cell_labs.flatten()

    out = {
        'nodes_flattened': don_nodes_flattened,
        'node_labels_flattened': don_node_labs_flattened,
        'cells_flattened': don_cells_flattened,
        'cell_centres_flattened': don_cell_centres_flattened,
        'cell_labels_flattened': don_cell_labs_flattened,
    }

    return out


def make_fine_plus_donut_mesh(crack_tip_radius_microns, fine_mesh_length=0.2, fine_mesh_element_length = 0.025, fan_box_num_side_elements=4,  fan_box_width=5, ret_crack_definition=True, size_behind_crack=0.25):
    """
    Build refined mesh = 4 rectilinear + 1 donut mesh at crack tip
    ''''''''''''''''''''''''''''''''''''''
    '               '    '               '
    '        5      '  4 '        3      '
    '               '    '               '
    '               '    '               '
    '               '    '               '
    '               '    '               '
    ''''''''''''''''''''''''''''''''''''''
    '        6      ' d-1'        2      '
    ''''''''''''''''' '  '               '
                      ''''''''''''''''''''
    
    Parameters
    ----------
    size_behind_crack : float
        Size of the mesh behind the crack as fraction of fine_mesh_length.
    
    """
    rad_mm = crack_tip_radius_microns * 1e-3
    # number of elements on side of fine mesh
    num_elem_side = int(fine_mesh_length / fine_mesh_element_length)
    line_elem =  num_elem_side - 1 #  substract corner element
    line_elem = round(line_elem / 3) * 3
    num_elem_side = line_elem + 1
    ref_mesh_size = [2 * num_elem_side * fine_mesh_element_length, num_elem_side * fine_mesh_element_length]
    
    #refined mesh nodes and cells
    ref_mesh_nodes_all = []
    ref_mesh_labs_all = []
    ref_mesh_cells_all = []
    ref_mesh_cell_centres_all = []
    ref_mesh_cell_labs_all = []
    
    # donut mesh nodes and cells    
    donut_mesh = make_donut_mesh(
        crack_tip_radius_microns,
        fine_mesh_element_length=fine_mesh_element_length, fan_box_num_side_elements=fan_box_num_side_elements, 
        fan_box_width=fan_box_width,
    )
        
    rn = int(donut_mesh['nodes_flattened'].shape[0] / fan_box_num_side_elements)
    don_mesh_size = [
            (donut_mesh['nodes_flattened'].reshape(rn, fan_box_num_side_elements, 2))[0][-1][0],
            (donut_mesh['nodes_flattened'].reshape(rn, fan_box_num_side_elements, 2))[0][-1][0]
    ]
    if ret_crack_definition:
        crack_line = donut_mesh['node_labels_flattened'].reshape(rn, fan_box_num_side_elements)[:,0]
        crack_front = donut_mesh['node_labels_flattened'].reshape(rn, fan_box_num_side_elements)[0]
    
    ref_mesh_nodes_all.append(donut_mesh['nodes_flattened'])
    ref_mesh_labs_all.append(donut_mesh['node_labels_flattened'])
    ref_mesh_cells_all.append(donut_mesh['cells_flattened'])
    ref_mesh_cell_centres_all.append(donut_mesh['cell_centres_flattened'])
    ref_mesh_cell_labs_all.append(donut_mesh['cell_labels_flattened'])


    rect_mesh_shifts = [
        [don_mesh_size[0], 0],
        [don_mesh_size[0], don_mesh_size[1]],
        [0, don_mesh_size[1]],
        [-ref_mesh_size[0] * size_behind_crack, don_mesh_size[1]],
        [-ref_mesh_size[0] * size_behind_crack, 0],
    ]
    rect_mesh_sizes = [
        [ref_mesh_size[0] / 2 - don_mesh_size[0], don_mesh_size[1]],
        [ref_mesh_size[0] / 2 - don_mesh_size[0], ref_mesh_size[1] - don_mesh_size[1]],
        [don_mesh_size[0], ref_mesh_size[1] - don_mesh_size[1]],
        [ref_mesh_size[0] * size_behind_crack, ref_mesh_size[1] - don_mesh_size[1]],
        [ref_mesh_size[0] * size_behind_crack, don_mesh_size[1] - rad_mm]
    ]

    for i in range(5):
        mesh_size = rect_mesh_sizes[i]
        mesh_grid = [int(mesh_size[0]/fine_mesh_element_length) + 1,
                     int(mesh_size[1]/fine_mesh_element_length) + 1]

        w = np.linspace(0, mesh_size[0], mesh_grid[0])  
        h = np.linspace(0, mesh_size[1], mesh_grid[1])

        if i == 4:
            h = list(donut_mesh['nodes_flattened'].reshape(rn, fan_box_num_side_elements, 2)[-1, :, 1])  # DON'T have this
            mesh_grid[1] = len(h)
        meshw, meshh = np.meshgrid(w, h)

        # shift meshgrid 
        meshw += rect_mesh_shifts[i][0]
        meshh += rect_mesh_shifts[i][1]

        # x, y coordinates of nodes
        nodes_coord = np.concatenate([meshw.flatten()[None],
                                      meshh.flatten()[None]], axis=0).T
        node_labs = np.zeros_like(meshw.flatten())
        reshaped = nodes_coord.reshape(mesh_grid[1], mesh_grid[0], 2)
        
        for ni, n in enumerate(nodes_coord):
            # find which nodes already been created 
            cond = np.isclose(n, np.concatenate(ref_mesh_nodes_all), 0.0001)
            found_idx = np.where(np.all(cond, axis=1))        
            if found_idx[0].size>0:
                node_labs[ni] = np.concatenate(ref_mesh_labs_all)[found_idx]

        # Apply mask over repeated border nodes
        mesh_mask = np.zeros_like(meshw)
        if i == 0:
            mesh_mask[:,0] = 1
        elif i == 1:
            mesh_mask[0,:] = 1
        elif i == 2:
            mesh_mask[0,:] = 1
            mesh_mask[:,-1] = 1
        elif i == 3:
            mesh_mask[:,-1] = 1
        elif i == 4:
            mesh_mask[-1,:] = 1
            mesh_mask[:,-1] = 1

        meshw, meshh = [np.ma.masked_array(i, mask=mesh_mask) for i in np.meshgrid(w, h)]
        # shift meshgrid 
        meshw += rect_mesh_shifts[i][0]
        meshh += rect_mesh_shifts[i][1]
        nodes_coord_masked = np.concatenate([meshw.compressed()[None],
                                      meshh.compressed()[None]], axis=0).T

        new_node_idx = np.where(node_labs==0)
        node_labs[new_node_idx] = np.arange(ref_mesh_labs_all[-1][-1]+1,
                                            ref_mesh_labs_all[-1][-1]+1 + len(meshw.compressed()))
        node_labs_grid = node_labs.reshape(meshw.shape)
        if i==0:
            # print('node_labs_grid[0,1:]: ', node_labs_grid[0,1:])
            crack_front = np.concatenate((crack_front, node_labs_grid[0,1:]))
            # print('front: ', front)
        elif i==4:
            crack_lip = node_labs_grid[0,:-1]
            # print('crack_lip: ', crack_lip)
        ref_mesh_nodes_all.append(nodes_coord_masked)
        ref_mesh_labs_all.append(node_labs[new_node_idx])

        # CHANGE!
        cells = (cells_from_nodes(node_labs_grid)).reshape(-1, 4)
        cells = order_nodes_elem(cells, np.concatenate(ref_mesh_nodes_all),
                                 np.concatenate(ref_mesh_labs_all))
        cell_centres = np.mean(
            cells_from_nodes(reshaped), axis=2
        )
        cell_centres_flattened = cell_centres.reshape((mesh_grid[0]-1) * (mesh_grid[1]-1), 2)
        cell_labs = np.arange(ref_mesh_cell_labs_all[-1][-1]+1,
                                            ref_mesh_cell_labs_all[-1][-1]+1 + (mesh_grid[0]-1) * (mesh_grid[1]-1))

        ref_mesh_cells_all.append(cells.astype('int'))
        ref_mesh_cell_centres_all.append(cell_centres_flattened)
        ref_mesh_cell_labs_all.append(cell_labs)
    crack_nodes = {
        'crack_lip': crack_lip.astype('int'),
        'crack_line': crack_line.astype('int'),
        'crack_front': crack_front.astype('int'),
    }

    out = {
        'node_coordinates': ref_mesh_nodes_all,
        'node_labels': ref_mesh_labs_all,
        'element_nodes': ref_mesh_cells_all,
        'element_centre_coordinates': ref_mesh_cell_centres_all,
        'element_labels': ref_mesh_cell_labs_all,
    }
    if ret_crack_definition:
        out['crack_definition'] = crack_nodes

    return out
