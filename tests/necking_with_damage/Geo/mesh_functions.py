import numpy as np
import matplotlib.pyplot as plt
from matplotlib import path
from numpy import linspace

class foreground():
    # Contains coor, areas, numpts, dx,dy,dz and right and left boundaries for areas
    def __init__(self, coor, vols, xyx_numpts_vec, dx_vec, dy_vec, dh_vec, xr, xl, yr, yl, hr, hl):
        self.coor = coor
        self.vols = vols
        self.xyz_numpts_vec = xyx_numpts_vec
        self.dx_vec = dx_vec
        self.dy_vec = dy_vec
        self.dh_vec = dh_vec
        self.xr = xr
        self.xl = xl
        self.yr = yr
        self.yl = yl
        self.hr = hr
        self.hl = hl

def generate_unif_particles(origin, xyzmax_vec, xyz_numpts_vec):
    # Generate areas and positions from input origin, 3 orthogonal max displacements in x, y and z as [xmax, ymax, zmax], and a vector of the number of points in x, y and z to generate as np.array([numptsx, numptsy, numptsz])

    xyz_totalnumpts_vec = xyz_numpts_vec
    node_num=xyz_totalnumpts_vec[0]*xyz_totalnumpts_vec[1]*xyz_totalnumpts_vec[2]
    x=np.zeros(xyz_totalnumpts_vec[0])
    y=np.zeros(xyz_totalnumpts_vec[1])
    h=np.zeros(xyz_totalnumpts_vec[2])
    
    dx_vec=np.zeros(node_num)
    dy_vec=np.zeros(node_num)
    dh_vec=np.zeros(node_num)
    
    xr_vec=np.zeros(node_num)
    xl_vec=np.zeros(node_num)
    yr_vec=np.zeros(node_num)
    yl_vec=np.zeros(node_num)
    hr_vec=np.zeros(node_num)
    hl_vec=np.zeros(node_num)
    
    coor=np.zeros((int(node_num),4))
    print('Generating Uniform Particles...')
    print("%s Nodes" %node_num)
    vol=np.zeros(node_num)
    dx=(xyzmax_vec[0]-origin[0])/(xyz_numpts_vec[0])
    dy=(xyzmax_vec[1]-origin[1])/(xyz_numpts_vec[1])
    if xyz_numpts_vec[2]>1: dh=(np.double(xyzmax_vec[2])-origin[2])/(np.double(xyz_numpts_vec[2]))
    if xyz_numpts_vec[2]==1: dh=0
    node=0
    for i in range(xyz_totalnumpts_vec[2]):
        if i>0: h[i]=h[i-1]+dh
        if i==0: h[i]=origin[2]+(0.5)*dh
        hr=h[i]+dh/2.0
        hl=h[i]-dh/2.0
        for j in range(xyz_totalnumpts_vec[0]):
            if j>0: x[j]=x[j-1]+dx
            if j==0: x[j]=origin[0]+(0.5)*dx
            xr=x[j]+dx/2.0
            xl=x[j]-dx/2.0
            for k in range(xyz_totalnumpts_vec[1]):
                if k>0: y[k]=y[k-1]+dy
                if k==0: y[k]=origin[1]+(0.5)*dy
                yr=y[k]+dy/2.0
                yl=y[k]-dy/2.0
                vol[node]=(xr-xl)*(yr-yl)
                if hr!=hl: vol[node]=vol[node]*(hr-hl)
                coor[node,:]=[node+1,x[j],y[k],h[i]]
                xr_vec[node]=xr
                xl_vec[node]=xl
                yr_vec[node]=yr
                yl_vec[node]=yl
                hr_vec[node]=hr
                hl_vec[node]=hl
                dx_vec[node]=xr-xl
                dy_vec[node]=yr-yl
                dh_vec[node]=hr-hl
                node=node+1
    coor=np.vstack(([node,0,0,0], coor))

    G=foreground(coor, vol, xyz_numpts_vec, dx_vec, dy_vec, dh_vec, xr_vec, xl_vec, yr_vec, yl_vec, hr_vec, hl_vec)
    print("Complete")
    return(G)

def vis_particles2d(G):
    # Generate Scatter plot of saved point cloud data
    plt.scatter(G.coor[1:,1], G.coor[1:,2])
    plt.show()

def save_geometry(G):
    # outputs Geometry.dat including input coordinates and areas of the PD nodes
    print("Saving input files...")
    coor=G.coor
    vol=G.vols
    dx_vec=np.transpose(np.array([G.dx_vec]))
    dy_vec=np.transpose(np.array([G.dy_vec]))
    dz_vec=np.transpose(np.array([G.dh_vec]))
    nodes=coor.shape[0]-1

    with open("Geometry.dat", "w") as f:
        f.write("%6d \t\t\t x \t\t\t\t\t  y\t\t\t\t\t\t z \t\t\t\t\tarea \n" % (nodes) )
    data=np.hstack((coor[1:nodes+1,:], np.transpose(np.array([vol[:nodes]]))))
    with open("Geometry.dat", "ab") as f:
        np.savetxt(f, np.vstack((data)), fmt='%6i %e %e %e %e')
    print("Save Complete, %d nodes"%nodes)


def subt_rect_dom_fg(G, p1, p2, p3, p4):
#Subtracts rectangle primitive domain from particles domain given 4 corner point vectors as nparray p1 p2 p3 p4
    xyz=G.coor[1:,1:-1]
    temp=G.coor[1:,:]
    p = path.Path([p1,p2,p3,p4])
    bools=p.contains_points(xyz)
    bools=np.invert(p.contains_points(xyz))
    nodenum=np.sum(bools)
    G.coor=temp[bools]
    for i in range(G.coor.shape[0]):
        G.coor[i,0]=i+1
    G.coor=np.vstack((np.array([nodenum,0,0,0]),G.coor))
    G.vols=G.vols[bools]
    G.dx_vec = G.dx_vec[bools]
    G.dy_vec = G.dy_vec[bools]
    G.dh_vec = G.dh_vec[bools]
    G.xr = G.xr[bools]
    G.xl = G.xl[bools]
    G.yr = G.yr[bools]
    G.yl = G.yl[bools]
    G.hr = G.hr[bools]
    G.hl = G.hl[bools]
    return(G)

def in_cube(xyz, xmin, xmax, ymin, ymax, zmin, zmax):
    out = np.logical_and(np.logical_and(np.logical_and(xyz[:, 0]>=xmin, xyz[:, 0]<=xmax), np.logical_and(xyz[:, 1]>=ymin, xyz[:, 1]<=ymax)), np.logical_and(xyz[:, 2]>=zmin, xyz[:, 2]<=zmax))
    return(out)

def subt_cubic_dom_fg(G, xmin, xmax, ymin, ymax, zmin, zmax):
    #Subtracts cubic primitive domain from particles domain given 4 corner point vectors as nparray p1 p2 p3 p4
    xyz=G.coor[1:,1:]
    temp=G.coor[1:,:]
    bools=in_cube(xyz, xmin, xmax, ymin, ymax, zmin, zmax)
    bools=np.invert(bools)
    nodenum=np.sum(bools)
    G.coor=temp[bools]
    for i in range(G.coor.shape[0]):
        G.coor[i,0]=i+1
    G.coor=np.vstack((np.array([nodenum,0,0,0]),G.coor))
    G.vols=G.vols[bools]
    G.dx_vec = G.dx_vec[bools]
    G.dy_vec = G.dy_vec[bools]
    G.dh_vec = G.dh_vec[bools]
    G.xr = G.xr[bools]
    G.xl = G.xl[bools]
    G.yr = G.yr[bools]
    G.yl = G.yl[bools]
    G.hr = G.hr[bools]
    G.hl = G.hl[bools]
    return(G)


def in_circle(xyz, center, radius):
    # Determine if point(s) is(are) within the 2d circular region defined by center and radius
    out = np.sqrt(np.sum(np.multiply(np.array(xyz)-center,np.array(xyz)-center), axis=1))>=radius
    return(out)
    
def subt_circular_domain(G, center, radius):
    # Subtracts cirucular domain from particles given center and radius
    xyz=G.coor[1:,1:]
    temp=G.coor[1:,:]
    bools=in_circle(xyz, center, radius)
    nodenum=np.sum(bools)
    G.coor=temp[bools]
    for i in range(G.coor.shape[0]):
        G.coor[i,0]=i+1
    G.coor=np.vstack((np.array([nodenum,0,0,0]),G.coor))
    G.vols=G.vols[bools]
    G.dx_vec = G.dx_vec[bools]
    G.dy_vec = G.dy_vec[bools]
    G.dh_vec = G.dh_vec[bools]
    G.xr = G.xr[bools]
    G.xl = G.xl[bools]
    G.yr = G.yr[bools]
    G.yl = G.yl[bools]
    G.hr = G.hr[bools]
    G.hl = G.hl[bools]
    return(G)
    
