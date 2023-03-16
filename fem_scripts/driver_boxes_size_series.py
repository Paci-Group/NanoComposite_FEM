#!/Users/bhenders/opt/miniconda3/envs/onelab/bin/python

# This script shows how to call getdp/gmsh from python with direct access to the
# onelab database. (The same basic principle can be used to create a python-based
# optimizer driving onelab clients.)

# You should run the script by opening it with Gmsh: either interactively (with
# 'File->Open') or in batch mode (with 'gmsh driver.py -')
# import the onelab python module
import onelab
import sys
sys.path.append('path/to/curr/dir')
import numpy as np
from scipy.optimize import least_squares, minimize
import json
import mesh_refinement as mr

# create a new onelab client
client = onelab.client(__file__)

# get Gmsh location from Gmsh options
mygmsh = client.getString('General.ExecutableFileName')

# Either force GetDP location, or also get it from Gmsh options: adapt the
# following line if you want to force the value:
mygetdp = '/Applications/Gmsh.app/Contents/MacOS/getdp'

client.sendInfo('Will use gmsh={0} and getdp={1}'.format(mygmsh, mygetdp))
client.sendInfo('PYTHONPATH={0}'.format(sys.path))

def get_energy_per_tet(results='results'):
    all_tets, etet = mr.scalar_field_from_parsed(f'{results}/etet.pos')
    return all_tets, etet.min(), etet.max(), etet.mean(), etet.std()

def get_inc_vol(results='results'):
    all_tets, vols = mr.scalar_field_from_parsed(f'{results}/vol.pos')
    return vols[:,0].sum()

# create a onelab variable for the model name
composite = client.defineString('Composite', value='box')

# get model file names with correct path
composite_geo = client.getPath(composite + '.geo')
composite_msh = client.getPath(composite + '.msh')
composite_pro = client.getPath(composite + '.pro')
results = client.getPath('results_cuboids')
summary = results + '/box_size_series.csv'

def generate_mesh(t, R, S, rf, per2pi=20, lcoarse=1.5):
    client.setNumber('Geometry/Parameters/lcoarse', value=lcoarse)
    client.addNumberChoice('Geometry/Parameters/lcoarse', lcoarse)
    client.setNumber('Geometry/Parameters/per2pi', value=per2pi)
    client.addNumberChoice('Geometry/Parameters/per2pi', per2pi)

    # Change the sphere radius
    client.setNumber('Geometry/Box/Sx', value=S[0])
    client.addNumberChoice('Geometry/Box/Sx', S[0])
    client.setNumber('Geometry/Box/Sy', value=S[1])
    client.addNumberChoice('Geometry/Box/Sy', S[1])
    client.setNumber('Geometry/Box/Sz', value=S[2])
    client.addNumberChoice('Geometry/Box/Sz', S[2])
    client.setNumber('Geometry/Inclusion/Rx', value=R[0])
    client.addNumberChoice('Geometry/Inclusion/Rx', R[0])
    client.setNumber('Geometry/Inclusion/Ry', value=R[1])
    client.addNumberChoice('Geometry/Inclusion/Ry', R[1])
    client.setNumber('Geometry/Inclusion/Rz', value=R[2])
    client.addNumberChoice('Geometry/Inclusion/Rz', R[2])
    client.setNumber("Geometry/Inclusion/Rf", value=rf)
    client.addNumberChoice("Geometry/Inclusion/Rf", rf)
    client.setNumber('Geometry/Inclusion/tinterface', value=t)
    client.addNumberChoice('Geometry/Inclusion/tinterface', t)

    # run gmsh as a subclient
    client.runSubClient('myGmsh', mygmsh + ' ' + composite_geo + ' -3 -optimize_netgen -v 2')


def get_eps_effective(eps_inc, eps_int, eps_mat, t, R, S, rf, target=1E-15):
    client.sendInfo(f'Trying to update file {composite_pro}')
    # Update epsilon epsilon values
    lines = []
    with open(composite_pro, 'r') as f:
        for line in f.readlines():
            if line.startswith("  epsr[Region[3]] ="):
                line = f"  epsr[Region[3]] = {eps_inc};\n"
            elif line.startswith("  epsr[Region[2]] ="):
                line = f"  epsr[Region[2]] = {eps_int};\n"
            elif line.startswith("  epsr[Region[1]] ="):
                line = f"  epsr[Region[1]] = {eps_mat};\n"
            lines.append(line)

    with open(composite_pro, 'w') as f:
        for line in lines:
            f.write(line)

    client.sendInfo(f'Updated file {composite_pro}')


    # unpack side lengths
    Sx, Sy, Sz = S

    # unpack inclusion Lengths
    Rx, Ry, Rz = R

    # volume loading
    c = Rx * 2
    S = Sx
    vol = (c - 2*rf)**3 + 6*rf*(c - 2*rf)**2 + 3*np.pi*rf**2*(c - 2*rf) + 4/3 * np.pi * rf**3
    vl = vol / S**3

    client.sendInfo(f'Beginning Calculation for VL~{vl}% -- Side Lengths: {Sx} {Sy} {Sz}, Rs: {Rx} {Ry} {Rz}, Interface: {t}, Fillet Radius: {rf}')

    # Model info needed to convert energy to a dielectric constant
    d = Sz / 2
    w = Sx / 2
    l = Sy / 2
    A = l * w
    V = 0.5
    eps0 = 8.85418782E-12

    # Start with these parameters
    per2pi = 20
    lcoarse = 1.5

    # target max energy per tet
    max_e = target * 5
    nit = 0

    while True:   # calculate solution until energy criterion met
        # run getdp as a subclient
        client.runSubClient('myGetDP', mygetdp + ' ' + composite_pro +
                            ' -msh ' + composite_msh +
                            ' -solve Electrostatics_v -v 2' +
                            ' -pos Electrostatics_v')

        client.setString('Gmsh/Action', value='refresh')

        # retrieve the energy and convert to dielectric constant
        energy = client.getNumber('Output/Energy [J]')
        epsr = 2 * energy * d / (A * eps0 * V ** 2)
        client.setNumber('Output/Dielectric Constant', value=epsr)
        client.addNumberChoice('Output/Dielectric Constant', epsr)

        dtot = client.getNumber("Output/Total Dipole Moment")
        dh = client.getNumber("Output/Matrix Dipole Moment Ver. 2")
        di = client.getNumber("Output/Inclusion Dipole Moment Ver. 2")
        dint = client.getNumber("Output/Interface Dipole Moment Ver. 2")

        # get the energy per tet
        tets, min_e, max_e, mean_e, std_e = get_energy_per_tet(results)
        v_inc = get_inc_vol(results)
        vl_actual = v_inc / (Sx * Sy * Sz / 8)

        client.setNumber('Output/Minimum Tet Energy [J]', value=min_e)
        client.addNumberChoice('Output/Minimum Tet Energy [J]', min_e)
        client.setNumber('Output/Maximum Tet Energy [J]', value=max_e)
        client.addNumberChoice('Output/Maximum Tet Energy [J]', max_e)
        client.setNumber('Output/Mean Tet Energy [J]', value=mean_e)
        client.addNumberChoice('Output/Mean Tet Energy [J]', mean_e)
        client.setNumber('Output/Std Tet Energy [J]', value=std_e)
        client.addNumberChoice('Output/Std Tet Energy [J]', std_e)
        client.setNumber('Output/Number of Tets', value=len(tets))
        client.addNumberChoice('Output/Number of Tets', len(tets))

        client.sendInfo(f'ITERATION {nit + 1}: EPS={epsr:.6f}')

        # If energy criterion is hit, don't re-mesh
        if max_e < target or nit > 10:
            with open(summary, 'a') as f:
                f.write(f"{Sx},{Sy},{Sz},{rf},{Rx},{Ry},{Rz},{vl},{vl_actual},{eps_inc},{eps_int},{eps_mat},{t},{lcoarse},{per2pi},{len(tets)},{target},{mean_e},{std_e},{max_e},{min_e},{energy},{epsr},{V},{di},{dint},{dh},{dtot}\n")
            return epsr, False    # False indicates no re-mesh
            break

        # else, refine mesh
        nit += 1
        per2pi += 20
        lcoarse *= 0.7

        generate_mesh(t, R, S, rf, per2pi=per2pi, lcoarse=lcoarse)

    return epsr, True    # True indicates that a re-mesh was required

def get_residuals(params, eps_mat, vls, Ss, cpmd_values, target):
    eps_inc, eps_int, t = params
    epsilons = np.ones(len(vls))
    try:
        for i, vl in enumerate(vls):
            epsilons[i] *= get_eps_effective(eps_inc, eps_int, eps_mat, t, vl, Ss[i], target)

        return epsilons - cpmd_values

    except Exception as e:
        client.sendInfo('get_residuals: Exited with error: {0}'.format(e))

def callback(xk):
    client.sendInfo(f'eps_inc: {xk[0]}, eps_int: {xk[1]}, t: {xk[2]}')

def get_filleted_cube(vol, rf):
    """Return the side length for a cube with fillet radius rf and volume vol"""
    c = 0.088189*(-976.3 * rf**3 + ((729*vol - 976.3 * rf**3)**2 - 1344607*rf**6)**(1/2) + 729 * vol)**(1/3) + 9.7337*rf**2 / ((-976.3 * rf**3 + ((729*vol - 976.3 * rf**3)**2 - 1344607*rf**6)**(1/2) + 729 * vol)**(1/3))
    return c

def get_objective(params, eps_mat, vls, Ss, cpmd_values, target):
    eps_inc, eps_int, t = params
    epsilons = np.ones(len(vls))
    try:
        for i, vl in enumerate(vls):
            client.sendInfo(f'eps_inc: {eps_inc}, eps_int: {eps_int}, eps_mat: {eps_mat}, t_int: {t}')

            # fillet radius (use VdW radius for silver)
            rf = 1.72
            vol = (vl/ 100) * Ss[i]**3
            R = get_filleted_cube(vol, rf) / 2  # get sidelength of inclusion
            epsilons[i] *= get_eps_effective(eps_inc, eps_int, eps_mat, t, [R]*3, [Ss[i]]*3, rf, target)

        objective = np.abs((epsilons - cpmd_values)).sum()
        client.sendInfo(f'Objective: {objective}')
        return objective

    except Exception as e:
        client.sendInfo('get_objective: Exited with error: {0}'.format(e))

# we're done if we don't do the actual calculation
if client.action == 'check':
   exit(0)

#######################
# The Optimization loop
#######################
try:
    with open(summary, 'w+') as f:
        f.write("Sx,Sy,Sz,Rf,Rx,Ry,Rz,vl,vl_fem,eps_inc,eps_int,eps_mat,t_interface,lcoarse,per2pi,Nel,target E,mean E,Std E,max E,min E,Energy,Eps,V,Di,Dint,Dh,Dtot\n")

    target=5E-10
    eps_incs = np.array([191])
    eps_ints = np.array([3.04])
    eps_mat = 3.04
    Ss = np.loadtxt('/Users/bhenders/Desktop/MgO_Ag_Project/Figures/Epsilon_infty_FEM/sidelengths.csv', delimiter=',', dtype=float)
    ts = np.array([1.])
    lcoarse=0.5
    per2pi = 50
    # fillet radius (use VdW radius for silver)
    rf = 1.72
    R = 2.83489575  # sidelength of inclusion

    client.sendInfo(f'Beginning Calculations...')

    for S in Ss:
        for t, eps_inc, eps_int in zip(ts, eps_incs, eps_ints):
            # remesh only when t or vl changes
            generate_mesh(t, [R]*3, S, rf, per2pi=50, lcoarse=0.5)
            # client.sendInfo(f'eps_inc: {eps_inc}, eps_int: {eps_int}, eps_mat: {eps_mat}, t_int: {t}')
            client.sendInfo(f'eps_inc: {eps_inc}, eps_int: {eps_int}, eps_mat: {eps_mat}, t_int: {t}')
            _, remesh = get_eps_effective(eps_inc, eps_int, eps_mat, t, [R]*3, S, rf, target)

            if remesh:  # regenerate the base mesh
                generate_mesh(t, [R]*3, S, rf, per2pi=per2pi, lcoarse=lcoarse)


    client.sendInfo(f'Finished Calculations.')
except Exception as e:
    client.sendInfo('Opt Loop: Exited with error: {0}'.format(e))
