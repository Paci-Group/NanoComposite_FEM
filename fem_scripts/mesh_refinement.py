import numpy as np
import os
from sklearn.neighbors import KNeighborsRegressor


def scalar_field_from_parsed(filepath):
    """
    Load a scalar field from a file in Gmsh parsed data format.
    
    The field should be a Scalar Field defined on the elements of
    a tetrahedral mesh. The file should be an ASCII .pos file.
    
    Parameters
    ----------
    
    filepath : str
        The path to the .pos file
        
    Returns
    -------
    
    tets : numpy.ndarray of np.float64
        The coordinates of each of the N tetrahedral elements 4 nodes in a N x 4 x 3 array. 
        
    scalar_vals : numpy.ndarray of np.float64
        The scalar values defined at each node in a N x 4 array. 
    """
    tets = []
    scalar_vals = []
    with open(filepath, 'r') as f:
        for line in f.readlines():
            try:
                line = line[3:-3]  # get rid of SS
                line = line.replace('{', '')  # make it easier to split
                pts, vals = line.split(')')
                tets.append(pts.split(','))
                scalar_vals.append(vals.split(','))
            except ValueError:
                pass
    tets = np.array(tets, dtype=np.float64)
    scalar_vals = np.array(scalar_vals, dtype=np.float64)
    return tets, scalar_vals


def pad_point_errors_by_box(points, scalars, xlo, xhi, ylo, yhi, zlo, zhi, pad=1.1):
    """
    Periodically repeat points and associated scalar data in 3D.
    
    Box must be centered at the origin.
    
    Parameters
    ----------
    points : np.ndarray
        N x 3 array of 3D coordinates
    scalars : np.ndarray
        N x 0 array of scalar data associated with the coordinates
    xlo : float
        The minimum x value of the periodic box
    xhi : float
        The maximum x value of the periodic box
    ylo : float
        The minimum y value of the periodic box
    yhi : float
        The maximum y value of the periodic box
    zlo : float
        The minimum z value of the periodic box
    zhi : float
        The maximum z value of the periodic box
    pad : float, default=0.1
        The proportion of each box side to periodically repeat and pad the cell with
        
    Returns
    -------
    (padded_points, padded_scalars) : tuple
        np.ndarrays of the resulting points and scalars that have been padded with repeats. 
    
    """
    
    # TODO: Clean this up. Right now, I build a potentially very large 
    # array and then shrink it back down. Should just add what we need 
    # to to start
    
    add_points = []
    add_scalars = []
    
    x_sel = [xlo, 0, xhi]
    y_sel = [ylo, 0, yhi]
    z_sel = [zlo, 0, zhi]
    
    for x_trans in [-1,0,1]:
        for y_trans in [-1,0,1]:
            for z_trans in [-1,0,1]:
                if (x_trans, y_trans, z_trans) == (0,0,0):
                    continue
                translation_vector = np.array([x_trans * (xhi-xlo), y_trans * (yhi-ylo), z_trans * (zhi-zlo)])
                add_points.append(points + translation_vector)
                add_scalars.append(scalars)
    
    padded_points = np.concatenate([points] + add_points, axis=0) 
    padded_scalars = np.concatenate([scalars] + add_scalars, axis=0)
    
    min_point = np.array([xlo, ylo, zlo]) * pad
    max_point = np.array([xhi, yhi, zhi]) * pad

    inside = np.all(np.logical_and(min_point <= padded_points, padded_points <= max_point), axis=1)
    
    return padded_points[inside], padded_scalars[inside]
                
    
    
#     # left add to right
#     add_points.append((points[points[:,0] < xlo + (xhi-xlo)*pad] + np.array([0, 0, xhi-xlo])).reshape(-1, 3))
#     add_scalars.append((scalars[points[:,0] < xlo + (xhi-xlo)*pad]))
#     # right add to left
#     add_points.append((points[points[:,0] > xhi - (xhi-xlo)*pad] - np.array([0, 0, xhi-xlo])).reshape(-1, 3))
#     add_scalars.append((scalars[points[:,0] > xhi - (xhi-xlo)*pad]))
    
#     # back add to front
#     add_points.append((points[points[:,1] < ylo + (yhi-ylo)*pad] + np.array([0, yhi-ylo, 0])).reshape(-1, 3))
#     add_scalars.append((scalars[points[:,1] < ylo + (yhi-ylo)*pad]))
#     # front add to back
#     add_points.append((points[points[:,1] > yhi - (yhi-ylo)*pad] - np.array([0, yhi-ylo, 0])).reshape(-1, 3))
#     add_scalars.append((scalars[points[:,1] > yhi - (yhi-ylo)*pad]))
    
#     # bottom add to top
#     add_points.append((points[points[:,2] < zlo + (zhi-zlo)*pad] + np.array([0, 0, zhi-zlo])).reshape(-1, 3))
#     add_scalars.append((scalars[points[:,2] < zlo + (zhi-zlo)*pad]))
#     # top add to botto
#     add_points.append((points[points[:,2] > zhi - (zhi-zlo)*pad] - np.array([0, 0, zhi-zlo])).reshape(-1, 3))
#     add_scalars.append((scalars[points[:,2] > zhi - (zhi-zlo)*pad]))
    
#     return np.concatenate([points] + add_points, axis=0), np.concatenate([scalars] + add_scalars, axis=0)


def mean_tet_edge_approx(vols):
    """
    Edge length of regular tetrahedra with given volumes
    
    V = a^3 / (6 * sqrt(2))
    
    Parameters
    ----------
    vols : numpy.ndarray
        The volumes of an array of tetrahedral elements
        
    Returns
    -------
    
    numpy.ndarray
        The side lengths
    """ 
    return (6 * np.sqrt(2) * vols)**(1/3)


def max_tet_edge(tets):
    """
    Compute the maximum edge length of tetrahedral elements
    
    Parameters
    ----------
    tets : numpy.ndarray
        An N x 4 x 3 array of the coordinates of the 4 nodes 
        of N tetrahedral elements
        
    Returns
    -------
    
    numpy.ndarray
        N x 0 array of the maximum side lengths of each element
    """ 
    a = np.sum((tets[:, 0, :] - tets[:, 1, :])**2, 1)**0.5
    b = np.sum((tets[:, 0, :] - tets[:, 2, :])**2, 1)**0.5
    c = np.sum((tets[:, 0, :] - tets[:, 3, :])**2, 1)**0.5
    d = np.sum((tets[:, 1, :] - tets[:, 2, :])**2, 1)**0.5
    e = np.sum((tets[:, 1, :] - tets[:, 3, :])**2, 1)**0.5
    f = np.sum((tets[:, 2, :] - tets[:, 3, :])**2, 1)**0.5
    return np.max(np.concatenate([a.reshape(-1,1), b.reshape(-1,1), c.reshape(-1,1), 
                                  d.reshape(-1,1), e.reshape(-1,1), f.reshape(-1,1)], axis=1), axis=1)


def mean_tet_edge(tets):
    """
    Compute the mean edge length of tetrahedral elements
    
    Parameters
    ----------
    tets : numpy.ndarray
        An N x 4 x 3 array of the coordinates of the 4 nodes 
        of N tetrahedral elements
        
    Returns
    -------
    
    numpy.ndarray
        N x 0 array of the mean side lengths of each element
    """
    a = np.sum((tets[:, 0, :] - tets[:, 1, :])**2, 1)**0.5
    b = np.sum((tets[:, 0, :] - tets[:, 2, :])**2, 1)**0.5
    c = np.sum((tets[:, 0, :] - tets[:, 3, :])**2, 1)**0.5
    d = np.sum((tets[:, 1, :] - tets[:, 2, :])**2, 1)**0.5
    e = np.sum((tets[:, 1, :] - tets[:, 3, :])**2, 1)**0.5
    f = np.sum((tets[:, 2, :] - tets[:, 3, :])**2, 1)**0.5
    return np.mean(np.concatenate([a.reshape(-1,1), b.reshape(-1,1), c.reshape(-1,1), 
                                  d.reshape(-1,1), e.reshape(-1,1), f.reshape(-1,1)], axis=1), axis=1)


def get_knn(loc, locs, k=10):
    """
    Find the indices of the k nearest neighbors of a given coordinate.
    
    Parameters
    ----------
    loc : numpy.ndarray
        A given coordinate composed of n parameters. n x 0 array
        
    locs : numpy.ndarray
        N x n array of coordinates
        
    k : int, default=10
        The number of nearest neighbors
        
    Returns
    -------
    
    numpy.ndarray
        k x 0 array of indices of the nearest neigbors of `loc` in `locs`
    """
    rs = locs - loc
    rnorm = np.linalg.norm(rs, axis=1)
    return rnorm.argsort()[:k+1]


def smooth_knn(data, locs, k=10):
    """
    Perform k nearest neighbors smoothing of a dataset.
    
    Return a new dataset where each value is the mean of 
    the former value and its k nearest neighbors, as determined
    by the Euclidean distance between the location of the current
    data point and others in locs.
    
    Parameters
    ----------
    data : numpy.ndarray
        A dataset of N scalar values. N x 0 array.
        
    locs : numpy.ndarray
        N x n array of coordinates in n-dimensional space
        
    k : int, default=10
        The number of nearest neighbors
        
    Returns
    -------
    
    numpy.ndarray
        N x 0 array that is the smoothed version of `data`.
    """
    print("getting unique...")
    ulocs, idx = np.unique(locs, return_index=True, axis=0)
    _, idx_inv = np.unique(locs, return_inverse=True, axis=0)
    data = data.copy()
    udata = data[idx]
    print("Starting KNN Smoothing...")
    for i in range(len(udata)):
        print(f"Smoothing element {i+1}", end='\r')
        indices_knn = get_knn(ulocs[i], ulocs, k)
        udata[i] = np.mean(udata[indices_knn])
    return udata[idx_inv]


def smooth_knn_fast(data, locs, k=10):
    """
    Perform k nearest neighbors smoothing of a dataset.
    
    Return a new dataset where each value is the mean of 
    the former value and its k nearest neighbors, as determined
    by the Euclidean distance between the location of the current
    data point and others in locs.
    
    This function uses scipy.neighbors.KNeighborsRegressor 
    for much improved performance over `mesh_refinement.smooth_knn`.
    
    Parameters
    ----------
    data : numpy.ndarray
        A dataset of N scalar values. N x 0 array.
        
    locs : numpy.ndarray
        N x n array of coordinates in n-dimensional space
        
    k : int, default=10
        The number of nearest neighbors
        
    Returns
    -------
    
    numpy.ndarray
        N x 0 array that is the smoothed version of `data`.
    """
    knn = KNeighborsRegressor(n_neighbors=k+1)
    knn.fit(locs, data)
    return knn.predict(locs)
   
    
def compute_size_field(prev_lengths, err, N):
    """
    Compute a new size field given a current size field and corresponding error estimates
    
    This function is adapted from the official Gmsh demo found in their GitLab Repo at
    https://gitlab.onelab.info/gmsh/gmsh/-/blob/master/demos/api/adapt_mesh.py.
    
    Parameters
    ----------
    prev_lengths : numpy.ndarray
        An N x 0 array of the current sizes of N elements
        
    err : numpy.ndarray
        N x 0 array of error estimates or error proxies for each element
        
    N : int
        The target number of elements in the refined mesh
        
    Returns
    -------
    
    numpy.ndarray
        N x 0 array, the new size field.
    """
    a = 2.0
    d = 2.0
    fact = (a**((2. + a)/(1. + a)) + a**(1./(1. + a))) * np.sum(err**(2./(1. + a)))
    ri = err**(2./(2.*(1 + a))) * a**(1./(d*(1. + a))) * ((1. + a)*N/fact)**(1./d)
    new_lengths = prev_lengths / ri
    
    # replace infinite (0 error) lengths with their previous values
    inf_idx = np.isinf(new_lengths)
    new_lengths[inf_idx] = prev_lengths[inf_idx]
    
    return new_lengths


def compute_new_scaling(prev_lengths, err, N):
    """
    Compute the scaling to apply to the current size field given error estimates
    
    This function is adapted from the official Gmsh demo found in their GitLab Repo at
    https://gitlab.onelab.info/gmsh/gmsh/-/blob/master/demos/api/adapt_mesh.py.
    
    Parameters
    ----------
    prev_lengths : numpy.ndarray
        An N x 0 array of the current sizes of N elements
        
    err : numpy.ndarray
        N x 0 array of error estimates or error proxies for each element
        
    N : int
        The target number of elements in the refined mesh
        
    Returns
    -------
    
    numpy.ndarray
        N x 0 array that is the factors by which to divide `prev_lengths`.
    """
    a = 2.0
    d = 2.0
    fact = (a**((2. + a)/(1. + a)) + a**(1./(1. + a))) * np.sum(err**(2./(1. + a)))
    ri = err**(2./(2.*(1 + a))) * a**(1./(d*(1. + a))) * ((1. + a)*N/fact)**(1./d)
    return ri

def compute_new_scaling_target(err, target):
    err = err / target 
    return (err)**(1/3)

def write_sf_to_pos(sf, tets, viewname='sf', outdir='.'):
    """
    Write a size field in Gmsh ASCII parsed format (.pos).
    
    Parameters
    ----------
    sf : numpy.ndarray
        N x 0 array. Scalar field to write to file
        
    tets : numpy.ndarray
        N x 4 x 3 array of the 3D coordinates of each of the 4 
        nodes of N tetrahedral elements
        
    viewname : str, default="sf"
        The name of the resulting Gmsh PostView
        
    outdir : str, default="."
        The output directory of the PostView file
      
    Returns
    -------
    None
    """
    with open(os.path.join(outdir, viewname), 'w+') as f:
        f.write(f'View "{viewname}" ' + '{\n')
        for tet, s in zip(tets, sf):
            f.write(f"SS({','.join([str(t) for t in tet.flatten()])})" + "{" + f"{s},{s},{s},{s}" + "};\n")
        f.write("};")
        
        
def results_stats(tets, err_proxy_1, err_proxy_2=None):
    mean_tet_edges = mean_tet_edge(tets)
    mesh_stats = f"""
Mesh Stats
----------
Number of Tetrahedra: {len(tets)}
Mean Tet Side Length: {mean_tet_edges.mean()}
Std. Tet Side Length: {mean_tet_edges.std()}
"""

    err_proxy_1_stats = f"""
Error Proxy #1 Stats
--------------------
Mean Err. Proxy: {err_proxy_1.mean()}
Std. Err. Proxy: {err_proxy_1.std()}
Min. Err. Proxy: {err_proxy_1.min()}
Max. Err. Proxy: {err_proxy_1.max()}
"""

    stats = [mesh_stats, err_proxy_1_stats]

    if err_proxy_2 is not None:
        err_proxy_2_stats = f"""
Error Proxy #2 Stats
--------------------
Mean Err. Proxy: {err_proxy_2.mean()}
Std. Err. Proxy: {err_proxy_2.std()}
Min. Err. Proxy: {err_proxy_2.min()}
Max. Err. Proxy: {err_proxy_2.max()}
"""
        stats.append(err_proxy_2_stats)
    
    return '\n\n'.join(stats)
    
  
    