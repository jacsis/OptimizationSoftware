import numpy as np
from models.mesh import Vertex, Triangle, Edge
from models.parameters import eps, epsH, alpha
from scipy.spatial import cKDTree
import random
from numpy.linalg import lstsq



def A_tot(Tarr):
    AT=0
    for tri in Tarr:
        AT += tri.areav
    return AT

def functional_C_Tot(Tarr, lamb, Qtot):
    FC = 0
    for tri1 in Tarr:
        FC0 = 0
        for tri2 in Tarr:
            if np.allclose(tri1.center ,tri2.center):
                continue
            FC0 += tri1.areav * tri2.scalar * tri2.areav / np.linalg.norm(tri1.center - tri2.center) 
        FC0 += lamb * tri1.areav
    FC += FC0**2
    AC=0
    for tri1 in Tarr:
        AC += tri1.scalar * tri1.areav
    AC += - Qtot
    return FC + AC**2




def functional_C0(Tarr, lamb, Qtot):
    FC = 0
    for tri1 in Tarr:
        FC0 = 0
        for tri2 in Tarr:
            if np.allclose(tri1.center ,tri2.center):
                continue
            dist = np.linalg.norm(tri1.center - tri2.center)
            FC0 += tri1.areav * tri2.scalar * tri2.areav / dist
        FC0 += lamb * tri1.areav
        FC += FC0  
    return FC


def functional_C(Tarr, lamb, Qtot):
    AC=0
    for tri1 in Tarr:
        AC += tri1.scalar * tri1.areav
    AC += - Qtot
    return functional_C0(Tarr, lamb, Qtot)**2 + AC**2

# def Matrix_pot(Earr, Tarr):
#     n=len(Tarr)
#     rand_EC = random.sample(Earr, n)
#     P_Tc = np.array([t.center for t in Tarr])
#     P_Ec = np.array([e.mp for e in rand_EC])

#         # Compute all pairwise distances: result is (n, n)
#     diff = P_Tc[:, np.newaxis, :] - P_Ec[np.newaxis, :, :]  # shape (n, n, 3)
#     dists_T = 1/np.linalg.norm(diff, axis=2)  # shape (n, n)
#     dists = dists_T.T
#     areas = np.array([[Tarr[i].areav] for i in range(n)])
#     dists_areas = dists * areas
#     new_column = np.array([[-1] for i in range(n)])
#     dists_areas = np.hstack((dists_areas, new_column))
#     new_row = areas
#     new_row = np.append(new_row, 0.)
#     matrix = np.vstack((dists_areas, new_row))

#     return matrix

def Matrix_pot(Earr, Tarr):
    P_Tc = np.array([t.center for t in Tarr])
    P_Ec = np.array([e.mp for e in Earr])

        # Compute all pairwise distances: result is (n, n)
    diff = P_Tc[:, np.newaxis, :] - P_Ec[np.newaxis, :, :]  # shape (n, n, 3)
    dists_T = 1/np.linalg.norm(diff, axis=2)  # shape (n, n)
    dists = dists_T.T
    areas = np.array([[Tarr[i].areav] for i in range(len(Tarr))])
    dists_areas = dists * areas.T
    new_column = np.array([[-1] for i in range(len(Earr))])
    dists_areas = np.hstack((dists_areas, new_column))
    new_row = areas
    new_row = np.append(new_row, 0.)
    matrix = np.vstack((dists_areas, new_row))

    return matrix

# def Matrix_pot(Earr, Tarr):
#     import random
#     import numpy as np

#     n = len(Tarr)
#     rand_EC = random.sample(Earr, n)

#     P_Tc = np.array([t.center for t in Tarr])     # Triangle centers
#     P_Ec = np.array([e.mp for e in rand_EC])      # Edge midpoints

#     # Compute pairwise distances: shape (n, n)
#     diff = P_Tc[:, np.newaxis, :] - P_Ec[np.newaxis, :, :]
#     dists = 1 / np.linalg.norm(diff, axis=2)

#     # Get triangle areas (shape (n, 1))
#     areas = np.array([t.area for t in Tarr]).reshape(-1, 1)

#     # Multiply each row of dists by corresponding triangle area
#     dists *= areas  # shape (n, n)

#     # Append column of -1s for Lagrange multiplier
#     new_column = np.full((n, 1), -1.0)
#     dists = np.hstack((dists, new_column))  # shape (n, n+1)

#     # Append row: triangle areas + 0 (for Lagrange term)
#     new_row = np.append(areas.flatten(), 0.0).reshape(1, -1)  # shape (1, n+1)

#     matrix = np.vstack((dists, new_row))  # shape (n+1, n+1)

#     return matrix


# def charg_dens_solv(Earr, Tarr, Qtot):
#     matrix = Matrix_pot(Earr, Tarr)
#     n = matrix.shape[0]
#     b = np.array([[0] for i in range(n-1)])
#     b = np.append(b, Qtot)
#     sol = np.linalg.solve(matrix, b)
#     rho = [sol[i] for i in range(len(sol)-1)]
#     print(Qtot)
#     #print(matrix)
#     #print(sol[len(sol)-1])
#     for i in range(len(rho)):
#         Tarr[i].scalar = rho[i]

def charg_dens_solv(Earr, Tarr, Qtot):
    matrix = Matrix_pot(Earr, Tarr)
    m = len(Earr)
    n = len(Tarr)

    b = np.zeros((m, 1))          # zero potentials at edge midpoints
    b = np.vstack((b, [[Qtot]]))  # add total charge constraint
    sol , _, _, _ = lstsq(matrix, b, rcond=None)
    rho = [sol[i] for i in range(len(sol)-1)]
    #print(matrix)
    #print(sol[len(sol)-1])
    for i in range(len(rho)):
        Tarr[i].scalar = rho[i]
    print(matrix)
    print(rho)
    print(b)
    

# # Optimized functional_C00 with k-d tree for efficient distance calculations
# def functional_C00(t, Tarr, lamb, Qtot, kdtree):
#     FC = 0
#     # Query the k-d tree for distances to other triangles
#     _, indices = kdtree.query(t.center, k=len(Tarr))
    
#     for idx in indices:
#         tri = Tarr[idx]
#         if np.allclose(tri.center, t.center):
#             continue
#         dist = np.linalg.norm(tri.center - t.center)  # Can be optimized further
#         FC += tri.areav * tri.scalar / dist
#     FC += lamb
#     return FC

# # Optimized Charge_Sum using vectorized approach
# def Charge_Sum(Tarr):
#     return np.sum([tri.areav * tri.scalar for tri in Tarr])

# # Optimized grad_Charge_r using vectorization and k-d tree
# def grad_Charge_r(t, Tarr, lamb, Qtot, kdtree):
#     F = 2 * t.areav * (Charge_Sum(Tarr) - Qtot)
#     _, indices = kdtree.query(t.center, k=len(Tarr))  # Query for all neighbors
#     for idx in indices:
#         tri = Tarr[idx]
#         if np.allclose(tri.center, t.center):
#             continue
#         dist = np.linalg.norm(tri.center - t.center)
#         c0 = functional_C00(tri, Tarr, lamb, Qtot, kdtree)
#         F += 2 * t.areav / dist * c0
#     return F

# # Optimized grad_Charge_l using batch processing and k-d tree
# def grad_Charge_l(Tarr, lamb, Qtot, kdtree):
#     FC = 0
#     # Query all pairs of triangles efficiently using k-d tree
#     for tri1 in Tarr:
#         FC0 = 0
#         _, indices = kdtree.query(tri1.center, k=len(Tarr))  # Nearest neighbors
#         for idx in indices:
#             tri2 = Tarr[idx]
#             if np.allclose(tri2.center, tri1.center):
#                 continue
#             dist = np.linalg.norm(tri1.center - tri2.center)
#             FC0 += tri2.scalar * tri2.areav / dist
#         FC0 += lamb
#         FC += 2 * FC0
#     return FC

# # Minimization using the optimized gradient calculations
# def Minimization(Tarr, lamb, Qtot, alpha=0.01):
#     # Build a k-d tree to speed up the distance calculations
#     points = np.array([tri.center for tri in Tarr])
#     kdtree = cKDTree(points)
    
#     # Iterate over the triangles and update the charge distribution
#     for tri in Tarr:
#         tri.scalar -= alpha * grad_Charge_r(tri, Tarr, lamb, Qtot, kdtree)
#     lamb -= alpha * grad_Charge_l(Tarr, lamb, Qtot, kdtree)
    
#     return lamb



def functional_C00(t : Triangle,Tarr, lamb, Qtot):
    FC = 0
    for tri in Tarr:
        if np.allclose(tri.center , t.center):
                continue
        dist = np.linalg.norm(tri.center - t.center)
        FC += tri.areav * tri.scalar  / dist
    FC += lamb   
    return FC

def Charge_Sum(Tarr):
    CS = 0
    for tri in Tarr:
        CS += tri.areav * tri.scalar
    return CS
        
def grad_Charge_r(t : Triangle, Tarr, lamb, Qtot):
    F = 2 * t.areav * (Charge_Sum(Tarr)-Qtot)
    for tri in Tarr:
        c0 = functional_C00(tri, Tarr, lamb, Qtot)
        if np.allclose(tri.center,t.center):
            continue
        dist = np.linalg.norm(tri.center - t.center)
        F += 2 * t.areav / dist * c0 
    return F

def grad_Charge_l(Tarr, lamb, Qtot):
    FC = 0
    for tri1 in Tarr:
        FC0 = 0
        for tri2 in Tarr:
            if np.allclose(tri2.center ,tri1.center):
                continue
            dist = np.linalg.norm(tri1.center - tri2.center)
            FC0 += tri2.scalar * tri2.areav / dist
        FC0 += lamb
        FC += 2 * FC0
    return FC



def Minimization(Tarr, lamb, Qtot):
    for tri in Tarr:
        tri.scalar = tri.scalar - alpha * grad_Charge_r(tri, Tarr, lamb, Qtot)
    lamb = lamb - alpha * grad_Charge_l(Tarr, lamb, Qtot)
    return lamb

def vertex_coords_to_array(Varr):
    return np.concatenate([v.coords for v in Varr])

def array_to_vertex_coords(arr, Varr):
    for i, v in enumerate(Varr):
        v.coords = arr[3*i : 3*(i+1)]

def Grad_H(Varr):
    Varr_none = []
    for v in Varr:
        if v.bound == None:
            Varr_none.append(v)
    n = 3 * len(Varr_none)
    G = np.zeros((n), dtype=np.float64)
    Varr_or = [v.coords.copy() for v in Varr_none]

    for i in range(n):
        ei = np.zeros(3); ei[np.mod(i,3)] = eps

        Varr_none[np.floor_divide(i,3)].coords = Varr_or[np.floor_divide(i,3)] + ei
        fp = functional(Varr_none[np.floor_divide(i,3)])
        Varr_none[np.floor_divide(i,3)].coords = Varr_or[np.floor_divide(i,3)]
        f0 = functional(Varr_none[np.floor_divide(i,3)])

        G[i] = (fp - f0) / (eps)

    return G

def Hessian(Varr):
    Varr_none = [v for v in Varr if v.bound is None]
    n = 3 * len(Varr_none)
    H = np.zeros((n, n), dtype=np.float64)
    epsH2 = epsH ** 2
    Varr_or = [v.coords.copy() for v in Varr_none]

    for i in range(n):
        vi, ki = divmod(i, 3)
        ei = np.zeros(3)
        ei[ki] = epsH

        for j in range(i, n):  # Only upper triangle
            vj, kj = divmod(j, 3)
            ej = np.zeros(3)
            ej[kj] = epsH

            # Save original coords
            Varr_none[vi].coords = Varr_or[vi].copy()
            Varr_none[vj].coords = Varr_or[vj].copy()

            if vi == vj:
                # Same vertex: combine ei and ej into one vector
                Varr_none[vi].coords = Varr_or[vi] + ei + ej
                fpp = functional_2(Varr_none[vi], Varr_none[vi])

                Varr_none[vi].coords = Varr_or[vi] + ei - ej
                fpm = functional_2(Varr_none[vi], Varr_none[vi])

                Varr_none[vi].coords = Varr_or[vi] - ei + ej
                fmp = functional_2(Varr_none[vi], Varr_none[vi])

                Varr_none[vi].coords = Varr_or[vi] - ei - ej
                fmm = functional_2(Varr_none[vi], Varr_none[vi])
            else:
                # Different vertices: perturb separately
                Varr_none[vi].coords = Varr_or[vi] + ei
                Varr_none[vj].coords = Varr_or[vj] + ej
                fpp = functional_2(Varr_none[vi], Varr_none[vj])

                Varr_none[vi].coords = Varr_or[vi] + ei
                Varr_none[vj].coords = Varr_or[vj] - ej
                fpm = functional_2(Varr_none[vi], Varr_none[vj])

                Varr_none[vi].coords = Varr_or[vi] - ei
                Varr_none[vj].coords = Varr_or[vj] + ej
                fmp = functional_2(Varr_none[vi], Varr_none[vj])

                Varr_none[vi].coords = Varr_or[vi] - ei
                Varr_none[vj].coords = Varr_or[vj] - ej
                fmm = functional_2(Varr_none[vi], Varr_none[vj])

            H_ij = (fpp - fpm - fmp + fmm) / (4 * epsH2)
            H[i, j] = H_ij
            if i != j:
                H[j, i] = H_ij  # Fill symmetric entry

            # Restore coordinates
            Varr_none[vi].coords = Varr_or[vi]
            if vi != vj:
                Varr_none[vj].coords = Varr_or[vj]

    return H


def G_Min(Varr):
    G = Grad_H(Varr)
    s = - alpha * Grad_H(Varr)

    i=0
    for v in Varr:
        if v.bound == None:
            v.coords += s[3 * i:3 * i + 3]
            i=i+1

def H_Min(Varr):
    H = Hessian(Varr)
    G = Grad_H(Varr)
    s = -np.linalg.solve(H, G)  # safer and faster than inv @ G

    i=0
    for v in Varr:
        if v.bound == None:
            v.coords += s[3 * i:3 * i + 3]
            i=i+1
