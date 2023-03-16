#!/Users/bhenders/opt/miniconda3/envs/onelab/bin/python

# This script shows how to call getdp/gmsh from python with direct access to the
# onelab database. (The same basic principle can be used to create a python-based
# optimizer driving onelab clients.)

# You should run the script by opening it with Gmsh: either interactively (with
# 'File->Open') or in batch mode (with 'gmsh driver.py -')

# import the onelab python module
import onelab
import sys
sys.path.append('/path/to/curr/dir')
import numpy as np
from scipy.optimize import least_squares, minimize
import json

# create a new onelab client
client = onelab.client(__file__)

# get Gmsh location from Gmsh options
mygmsh = client.getString('General.ExecutableFileName')

# Either force GetDP location, or also get it from Gmsh options: adapt the
# following line if you want to force the value:
mygetdp = '/Applications/Gmsh.app/Contents/MacOS/getdp'

client.sendInfo('Will use gmsh={0} and getdp={1}'.format(mygmsh, mygetdp))
client.sendInfo('PYTHONPATH={0}'.format(sys.path))

try:
    import mesh_refinement as mr
except Exception as e:
    client.sendInfo('Exited with error: {0}'.format(e))


def get_energy_per_tet(results='results'):
    all_tets, etet = mr.scalar_field_from_parsed(f'{results}/etet.pos')
    return all_tets, etet.min(), etet.max(), etet.mean(), etet.std()

results = client.getPath('results')

R = 1
summary = results + '/Optimization_run.csv'

def get_eps_effective(composite, eps_inc, eps_int, eps_mat, t_interface, vl, S, target=1E-15):

    # get model file names with correct path
    # composite_geo = client.getPath(composite + '.geo')
    composite_msh = client.getPath(composite + '.msh')
    composite_pro = client.getPath(composite + '.pro')

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
    # Side Length for prescribed VL and Radius 1
    R = (S**3 * 3/4 / np.pi * (vl / 100))**(1/3)
    # S = (4/3 * np.pi * R**3 / (vl / 100))**(1/3)
    client.sendInfo(f'Beginning Calculation for VL={vl}% -- Side Length: {S}, R: {R}, Interface: {t_interface}')

    # Model info needed to convert energy to a dielectric constant
    d = S / 2
    w = S / 2
    l = S / 2
    A = l * w
    V = 0.5
    eps0 = 8.85418782E-12

    # Start with these parameters
    per2pi = 30
    lcoarse = 1

    # target max energy per tet
    max_e = target * 5
    nit = 0

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

    # get the energy per tet
    tets, min_e, max_e, mean_e, std_e = get_energy_per_tet(results)

    dtot = client.getNumber("Output/Total Dipole Moment")
    dh = client.getNumber("Output/Matrix Dipole Moment Ver. 2")
    di = client.getNumber("Output/Inclusion Dipole Moment Ver. 2")
    dint = client.getNumber("Output/Interface Dipole Moment Ver. 2")

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

    with open(summary, 'a') as f:
        f.write(f"{S},{R},{vl},{eps_inc},{eps_int},{eps_mat},{t_interface},{lcoarse},{per2pi},{len(tets)},{target},{mean_e},{std_e},{max_e},{min_e},{energy},{epsr}\n")

    return epsr, di / dtot

def get_residuals(params, eps_mat, t_interface, composites, vls, Ss, cpmd_values, target):
    eps_inc, eps_int = params
    epsilons = np.ones(len(vls))
    try:
        for i, vl in enumerate(vls):
            epsilons[i] *= get_eps_effective(composites[i], eps_inc, eps_int, eps_mat, t_interface, vl, Ss[i], target)

        return epsilons - cpmd_values

    except Exception as e:
        client.sendInfo('get_residuals: Exited with error: {0}'.format(e))

def callback(xk):
    client.sendInfo(f'eps_inc: {xk[0]}, eps_int: {xk[1]}')

def get_objective(params, eps_mat, t_interface, composites, vls, Ss, cpmd_values, target):
    eps_inc, eps_int = params
    epsilons = np.ones(len(vls))
    try:
        di_dtots = []
        for i, vl in enumerate(vls):
            client.sendInfo(f'eps_inc: {eps_inc}, eps_int: {eps_int}, eps_mat: {eps_mat}, t_interface: {t_interface}')
            eps, di_dtot = get_eps_effective(composites[i], eps_inc, eps_int, eps_mat, t_interface, vl, Ss[i], target)
            epsilons[i] *= eps
            di_dtots.append(di_dtot)
        client.sendInfo(f'Di / Dtot: {di_dtots}')
        objective = np.abs((epsilons - cpmd_values) / cpmd_values).sum() + np.abs((di_dtots[0]-0.0962)/0.0962)
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
        f.write("S,R,vl,eps_inc,eps_int,eps_mat,t_interface,lcoarse,per2pi,Nel,target E,mean E,Std E,max E,min E,Energy,Eps\n")
    target=5E-10
    eps_inc = 60
    eps_int = 40
    eps_mat = 9.26
    t_interface = 1
    composites = [client.defineString('Composite_6', value='ellipsoid_6'), client.defineString('Composite_8', value='ellipsoid_8')]

    initial_guess = np.array([eps_inc, eps_int])
    bounds = (np.array([9.26, 10000]), np.array([9.26, 10000]))
    vls = np.array([6.754800892, 2.975746799])
    Ss = np.array([12.926, 17.077])
    cpmd_values = np.array([11.94, 10.303])

    client.sendInfo(f'Beginning Optimization Loop...')
    result = minimize(get_objective, initial_guess, callback=callback, bounds=bounds, args=(eps_mat, t_interface, composites, vls, Ss, cpmd_values, target),
                      method='Nelder-Mead')
    client.sendInfo(f'Finished Optimization Loop. Dumping Results')

    np.savetxt("optimize.csv", result)

    # with open(results + 'optimize.json', 'w+') as f:
    #     f.write(result_string)

    client.sendInfo(f'Dumped Results to {results + "optimize.json"}')
except Exception as e:
    client.sendInfo('Opt Loop: Exited with error: {0}'.format(e))
