import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

'''
     (C) Patrick Lestrange 2017
     
     A post-processing script to grab the PHF rotated densities
     in a CQ output file and calculate the Euclidean distance 
     between them
'''

def grab_matrices(cq_file,nbasis):
    
    # determine the total number of grid points
    search_file = open(cq_file, "r")
    for line in search_file:
        if 'LebedevTrap Grid' in line:
            contents = re.split('\s+|\(|\)|\,',line)
            ngridT = int(contents[5])*int(contents[6])
    search_file.close()

    read_den = False
    angles, ang_dict = [], {}
    quaterions = ['Scalar', 'Mz', 'My', 'Mx']
    factor, ang_count, den_count, quat_ind = 1, -1, 0, 0
    den_arrays = np.zeros([ngridT,4,nbasis*nbasis],dtype=complex)
    search_file = open(cq_file, "r")
    for line in search_file:

        # collect the angles in a list of dictionaries
        if 'Alpha' in line:
            ang_dict = {}
            contents = re.split('\s+',line)
            ang_dict['alpha'] = float(contents[2])
        if 'Beta' in line:
            contents = re.split('\s+',line)
            ang_dict['beta'] = float(contents[2])
        if 'Gamma' in line:
            contents = re.split('\s+',line)
            ang_dict['gamma'] = float(contents[2])
            angles.append(ang_dict)
            ang_count += 1

        # determine what kind of density we're grabbing
        if 'Rot Density' in line:
            contents = re.split('\s+|\(|\)',line)
            quat_ind = quaterions.index(contents[3])
            if 'Re' in line:
                factor = 1
            else:
                factor = 1j
 
        # Grab densities from the output file
        if '-----' in line:
            read_den = not read_den
        if read_den:
            contents = re.split('\s+|--+',line)
            for element in contents:
                if 'e+' in element or 'e-' in element:
                    den_arrays[ang_count,quat_ind,den_count] += factor*float(element)
                    den_count += 1
        if not read_den:
            den_count = 0 

    search_file.close()

    # redimension the rot densities into matrices
    rot_dens = np.zeros([ngridT,4,nbasis,nbasis],dtype=complex)
    for g in range(ngridT):
        for q in range(4):
            ij = 0
            for i in range(nbasis):
                for j in range(nbasis):
                    rot_dens[g,q,i,j] = den_arrays[g,q,ij]
                    ij += 1                        

    # make list of angle labels
    ang_list = []
    for g in range(ngridT):
        string = "(%.2f,%.2f,%.2f)" % (angles[g]['alpha'], angles[g]['beta'], angles[g]['gamma'])
        ang_list.append(string)
        
#   print "rot den = \n", rot_dens[1,1].real
#   print angles[1]['alpha']

    return ngridT, angles, ang_list, rot_dens

def euc_distance(rot_dens,ngridT):

    # calculate the euclidean distance between all rotated densities
    min_dist = []
    dist = np.zeros([4,ngridT,ngridT])
    for q in range(4):
        for g in range(ngridT):
            dist_min = 1000.
            dist_ind = g
            for h in range(g):
#               adiff = angles[g]['alpha'] - angles[h]['alpha']
#               bdiff = angles[g]['beta']  - angles[h]['beta']
#               gdiff = angles[g]['gamma'] - angles[h]['gamma']
                dist[q,g,h] = np.linalg.norm(rot_dens[g,q]-rot_dens[h,q],'fro')
                if abs(dist[q,g,h]) < dist_min and g!=h:
                    dist_min = dist[q,g,h]
                    dist_ind = h
            if (g,dist_ind) not in min_dist and (dist_ind,g) not in min_dist:
                min_dist.append((g,dist_ind)) 
        dist[q] /= np.max(dist[q])

    return dist
   
def track_zeros():

    # determine the number of zeros in the original ordering
    ntotal, nzeros = 0, 0
    for g in range(ngridT-1):
        ntotal += 1
        dist = np.linalg.norm(rot_dens[g,0]-rot_dens[g+1,0],'fro')
        if abs(dist) < 1e-6:
            nzeros += 1
    print "NGrid = %d, NZeros = %d" % (ngridT, nzeros)
    print ntotal

    for i in range(16,30):
        for j in range(16,i):
            print "i = ", i, " angle = ",angles[i]
            print "j = ", j, " angle = ",angles[j]
            print np.linalg.norm(rot_dens[i,0]-rot_dens[j,0],'fro')

    # keep track of the min distances for each angle
    for angle in angles:
        angle['closest'] = []
    for i, j in min_dist:
        angles[i]['closest'].append(j)
        angles[j]['closest'].append(i)
    for g, angle in enumerate(angles):
#       print len(angle['closest'])
        if len(angle['closest']) >= 4:
            print "Angle", g, " = ", angle       
            for i in angle['closest']:
                print angles[i]  
            print "\n"

def plot_euc_heatmap(dist):

    fig = plt.figure()

    plt.subplot(221)
    plt.imshow(dist[0], cmap='Blues', interpolation='nearest')
    plt.title("Scalar P")
    plt.colorbar()

    plt.subplot(222)
    plt.imshow(dist[1], cmap='Blues', interpolation='nearest')
    plt.title("Mz P")
    plt.colorbar()

    plt.subplot(223)
    plt.imshow(dist[2], cmap='Blues', interpolation='nearest')
    plt.title("My P")
    plt.colorbar()

    plt.subplot(224)
    plt.imshow(dist[3], cmap='Blues', interpolation='nearest')
    plt.title("Mx P")
    plt.colorbar()

    plt.tight_layout()
    plt.show()

def arc_lengths(angles,ngridT):

#   arcs = np.zeros([ngridT,ngridT])
#   for g, ang1 in enumerate(angles):
#       arcs[g,g]  =  np.cos(ang1['alpha'])*np.cos(ang1['alpha'])
#       arcs[g,g] += (np.sin(ang1['alpha'])*np.cos(ang1['gamma'])*
#                     np.sin(ang1['alpha'])*np.cos(ang1['gamma']))
#       arcs[g,g] += (np.sin(ang1['alpha'])*np.sin(ang1['gamma'])*np.cos(ang1['beta'])*
#                     np.sin(ang1['alpha'])*np.sin(ang1['gamma'])*np.cos(ang1['beta']))
#       arcs[g,g] += (np.sin(ang1['alpha'])*np.sin(ang1['gamma'])*np.sin(ang1['beta'])*
#                     np.sin(ang1['alpha'])*np.sin(ang1['gamma'])*np.sin(ang1['beta']))
#       arcs[g,g] = np.sqrt(arcs[g,g])
#       print "vector length = ", arcs[g,g]

    arcs = np.zeros([ngridT,ngridT])
    for g, ang1 in enumerate(angles):
#       for h, ang2 in enumerate(angles):
        for h in range(g):
            arcs[g,h]  =  np.cos(ang1['alpha'])*np.cos(angles[h]['alpha'])
            arcs[g,h] += (np.sin(ang1['alpha'])*np.cos(ang1['gamma'])*
                          np.sin(angles[h]['alpha'])*np.cos(angles[h]['gamma']))
            arcs[g,h] += (np.sin(ang1['alpha'])*np.sin(ang1['gamma'])*np.cos(ang1['beta'])*
                          np.sin(angles[h]['alpha'])*np.sin(angles[h]['gamma'])*np.cos(angles[h]['beta']))
            arcs[g,h] += (np.sin(ang1['alpha'])*np.sin(ang1['gamma'])*np.sin(ang1['beta'])*
                          np.sin(angles[h]['alpha'])*np.sin(angles[h]['gamma'])*np.sin(angles[h]['beta']))
            # necessary due to small inaccuracies in the angles, but these points are actually -1
            if arcs[g,h] < -1: 
#               print "arc length out of range: ", arcs[g,h]
                arcs[g,h] = -1
            arcs[g,h]  =  np.arccos(arcs[g,h])
#   arcs /= np.max(arcs)

    return arcs

def plot_arc_heatmap(dist,arcs):

    fig = plt.figure()

    plt.subplot(131)
    plt.imshow(arcs/np.max(arcs), cmap='Blues', interpolation='nearest')
    plt.title("Arc Lengths")
    plt.colorbar()

    plt.subplot(132)
    plt.imshow(dist[0], cmap='Blues', interpolation='nearest')
    plt.title("Scalar P")
    plt.colorbar()

    plt.subplot(133)
    plt.imshow(arcs-dist[0], cmap='bwr', interpolation='nearest')
    plt.title("Difference")
    plt.colorbar()

    plt.tight_layout()
    plt.show()

def angle_differences(angles,ngridT):

    diffs = np.zeros([3,ngridT,ngridT])
    for i, key in enumerate(angles[0]):
        for g in range(ngridT):
            for h in range(g):
                diffs[i,g,h] = angles[g][key] - angles[h][key]
#       diffs[i] /= np.max(diffs[i])

    return diffs

def plot_angle_differences(dist,diffs):

    fig = plt.figure()

    plt.subplot(331)
    plt.imshow(dist[0], cmap='Blues', interpolation='nearest')
    plt.title("Scalar P")
    plt.colorbar()

    plt.subplot(332)
    plt.imshow(diffs[0], cmap='bwr', interpolation='nearest')
    plt.title("Alpha")
    plt.colorbar()

    plt.subplot(333)
    plt.imshow(diffs[1], cmap='bwr', interpolation='nearest')
    plt.title("Beta")
    plt.colorbar()
    
    plt.subplot(334)
    plt.imshow(diffs[2], cmap='bwr', interpolation='nearest')
    plt.title("Gamma")
    plt.colorbar()

    sum_diffs = diffs[0]+diffs[1]+diffs[2]

    plt.subplot(335)
    plt.imshow(sum_diffs, cmap='bwr', interpolation='nearest')
    plt.title("Sum")
    plt.colorbar()

    plt.subplot(336)
    plt.imshow((sum_diffs/np.max(sum_diffs))-dist[0], cmap='bwr', interpolation='nearest')
    plt.title("Normed Sum-P")
    plt.colorbar()

    plt.tight_layout()
    plt.show()

def q_inv(q):

    norm = 0.
    for i in range(4):
        norm += q[i]**2

    q_inv = q
    for i in range(1,4):
        q_inv[i] *= -1
    q_inv /= norm

    return q_inv

def q_mult(a,b):

    ab = np.zeros([4])

    ab[0] = a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3]
    ab[1] = a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2]
    ab[2] = a[0]*b[2] + a[2]*b[0] + a[3]*b[1] - a[1]*b[3]
    ab[3] = a[0]*b[3] + a[3]*b[0] + a[1]*b[2] - a[2]*b[1]

    return ab

def q_slerp(v0,v1,t):

    dot = 0.
    for i in range(4):
        dot += v0[i]*v1[i]

    # if they're too close, linearly interpolate
    dot_thresh = 0.9995
    if dot > dot_thresh:
        slerp = v0 * t*(v1 - v0)
        slerp /= np.linalg.norm(slerp)

    # Find slerp
    else: 
        if dot < 0:
            v1 = -v1
            dot = -dot
            
        # clamp dot product to right range
        max(min(dot, 1.), -1)
        theta_0 = np.arccos(dot)
        theta = theta_0 * t
     
        v2 = v1 - v0*dot
        v2 /= np.linalg.norm(v2)
    
        slerp = v0*np.cos(theta) + v2*np.sin(theta)

    return slerp

def quaternion_rep(angles,ngridT):

    quats     = np.zeros([ngridT,4])
    q_angles  = np.zeros([ngridT,ngridT])
    slerp_len = np.zeros([ngridT,ngridT])

    # transform euler angles to unit quaternions
    for g, angle in enumerate(angles):
        quats[g,0] = ( np.cos(angle['beta']/2) *
                       np.cos((angle['alpha']+angle['gamma'])/2) )
        quats[g,1] = (-np.sin(angle['beta']/2) *
                       np.sin((angle['alpha']-angle['gamma'])/2) )
        quats[g,2] = ( np.sin(angle['beta']/2) *
                       np.cos((angle['alpha']-angle['gamma'])/2) )
        quats[g,3] = ( np.cos(angle['beta']/2) *
                       np.sin((angle['alpha']+angle['gamma'])/2) )

    # evaluate angles between all unit quaternions ensuring we're going the short way
    for g in range(ngridT):
        for h in range(g):
#           q_angles[g,h] = 2*np.arccos(np.abs(np.dot(quats[g],quats[h])))
            q_angles[g,h] = np.arccos(2*np.dot(quats[g],quats[h])**2 - 1)

    # integrate over rotation to get arc length
    nslerp = 10
    for g in range(ngridT):
        for h in range(g):
            # smoothly iterpolate between two quaternions (slerp)
            ang = q_angles[g,h]
            for i in range(nslerp-1):
                t = float(i)/float(nslerp)
                slerp1 = q_slerp(quats[g],quats[h],t)
                t = float(i+1)/float(nslerp)
                slerp2 = q_slerp(quats[g],quats[h],t)
                slerp_len[g,h] += np.linalg.norm(slerp1-slerp2)
                
    return quats, q_angles, slerp_len

def plot_quaternion_angles(dist,q_angles,slerp_len):

    fig = plt.figure()

    plt.subplot(131)
    plt.imshow(dist[0], cmap='Blues', interpolation='nearest')
    plt.title("Scalar P")
    plt.colorbar()

    plt.subplot(132)
    plt.imshow(2*q_angles, cmap='Blues', interpolation='nearest')
    plt.title("Quaternion Angles")
    plt.colorbar()

    plt.subplot(133)
    plt.imshow(slerp_len, cmap='Blues', interpolation='nearest')
    plt.title("Arc Lengths")
    plt.colorbar()

    plt.tight_layout()
    plt.show()
 
if __name__ == '__main__':
 
    # grab rotated densities from the output file 
#   cq_file = 'o2-631g-read.out'
#   nbasis  = 18
    cq_file = 'h3-sto3g-read-dense.out' 
    nbasis  = 3
    ngridT, angles, ang_list, rot_dens = grab_matrices(cq_file,nbasis)

    # calculate euclidean distance between each rotated density matrix
    dist = euc_distance(rot_dens,ngridT)
#   plot_euc_heatmap(dist)

    # calculate the differences between each rotation angle
#   diffs = angle_differences(angles,ngridT)
#   plot_angle_differences(dist,diffs)

    # calculate the angles and slerp lengths between the quaternion 
    # representation of the rotations
    quats, q_angles, slerp_len = quaternion_rep(angles,ngridT)
    plot_quaternion_angles(dist,q_angles,slerp_len)

    # calculate the arc lengths between each point
#   arcs = arc_lengths(angles,ngridT)
#   plot_arc_heatmap(dist,arcs)


