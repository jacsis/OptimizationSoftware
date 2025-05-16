import numpy as np
from models.mesh import Vertex, B_Cond
from models.parameters import eps, epsH, alpha



def functional(v : Vertex):
    A=0
    for tri in v.t:
        A = A + tri.area()
    return A

def functional_2(v1 : Vertex, v2 : Vertex):
    A = 0
    T_tot =  list(set(v1.t + v2.t))
    for tri in T_tot:
        A = A + tri.area()
    return A





        
def grad(v: Vertex, k):
    veps = eps * np.eye(3)[k]  # directional perturbation
    f0 = functional(v)
    perturbed_coords = v.coords + veps
    original_coords = v.coords.copy()

    v.coords = perturbed_coords
    f_eps = functional(v)
    v.coords = original_coords

    return (f_eps - f0) / eps



def Minimization(Varr):
    for ver in Varr:
        if ver.bound is None:
            g = np.array([grad(ver, k) for k in range(3)])
            ver.coords -= alpha * g

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
