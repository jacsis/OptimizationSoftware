import numpy as np
import plotly.graph_objects as go
from vedo import Plotter, Mesh
from skimage import measure
import copy

class Vertex:
    def __init__(self, id, coords):
        self.id = id
        self.coords = np.array(coords,np.float64)  # Store coordinates as a NumPy array
        self.e = []  # To store the edges' ids that contain the vertex
        self.t = []  # To store the triangles' ids that contain the vertex
        self.bound = None # To be changed to store the type of boundary condition, if any
    
    def copy(self):
        return copy.deepcopy(self)

class Edge:
    def __init__(self, id, v1: Vertex, v2: Vertex):
        self.id = id
        self.v1 = v1
        self.v2 = v2
        self.t = []  # To store the triangles' ids that contain the edge
        self.mp = [] # To store the middle point of the edge
        
    def len(self):
        l = np.linalg.norm(self.v1.coords - self.v2.coords)
        return float(l)

def area_v(v1 : Vertex, v2 : Vertex, v3 : Vertex ):
    AB = v2.coords - v1.coords
    AC = v3.coords - v1.coords
    cross_prod = np.cross(AB, AC)
    area = 0.5 * np.linalg.norm(cross_prod)
    return float(area)

class Triangle:
    def __init__(self, id, v1: Vertex, v2: Vertex, v3: Vertex):
        self.id = id
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3
        self.center = (v1.coords+v2.coords+v3.coords)/3 
        self.areav = area_v(self.v1,self.v2,self.v3)
        self.n = None  # Will store the normal vector after calculation
        self.scalar = None # Will store the attribute after calculation
        
        

    #def per(self):
    #    return sum(edge.len() for edge in self.edges)

    def area(self):
        AB = self.v2.coords - self.v1.coords
        AC = self.v3.coords - self.v1.coords
        cross_prod = np.cross(AB, AC)
        area = 0.5 * np.linalg.norm(cross_prod)
        return float(area)
    
    def get_vertices(self):
        return [self.v1,self.v2,self.v3]



    def normal(self):
        AB = self.v2.coords - self.v1.coords
        AC = self.v3.coords - self.v2.coords
        cross_prod = np.cross(AB, AC)
        norm = np.linalg.norm(cross_prod)
        if norm == 0:
            raise ValueError("Degenerate triangle with zero area")
        self.n = cross_prod / norm
        return self.n
    
class Quadr:
    def __init__(self, id, v1: Vertex, v2: Vertex, v3: Vertex, v4: Vertex):
        self.id = id
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3
        self.v4 = v4

class Mesh:
    def __init__(self, id):
        self.id = id
        self.verts = []
        self.edges= []
        self.triangles= []

    def area(self):
        A=0
        for t in self.triangles:
            A = A + t.area()




class B_Cond:
    def __init__(self, id, shape, type):
        self.id = id
        self.shape = shape
        self.type = type

        thetas = np.linspace(0, 2 * np.pi, 10001)  # Sample θ from 0 to 2π
        self.shape_p  = np.array([[i(theta) for i in self.shape ] for theta in thetas])  # Shape (10001, 2)

    def shape_f(self ,x):
        return [f(x) for f in self.shape]
    


def Vmove(coo0, b : B_Cond):
    
    dist = np.linalg.norm(b.shape_p - coo0, axis=1) # Compute distances to coo0

    # Find closest point
    index = np.argmin(dist)
    m = b.shape_p[index]
    return m

def Triangulation(Varr, Qarr):
    Tarr0 = []
    Varr0 = Varr.copy()
    lV = len(Varr0)
    next_vid = lV + 1
    next_tid = 1

    for q in Qarr:
        # Compute quad center
        qm = Vertex(next_vid, (q.v1.coords + q.v2.coords + q.v3.coords + q.v4.coords) / 4)
        Varr0.append(qm)
        next_vid += 1

        # Create four triangles
        Tarr0.append(Triangle(next_tid, q.v2, q.v1, qm)); next_tid += 1
        Tarr0.append(Triangle(next_tid, q.v4, q.v2, qm)); next_tid += 1
        Tarr0.append(Triangle(next_tid, q.v3, q.v4, qm)); next_tid += 1
        Tarr0.append(Triangle(next_tid, q.v3, qm, q.v1)); next_tid += 1

    # Re-index vertices by ID
    Varr00 = [None] * len(Varr0)
    for v in Varr0:
        Varr00[v.id - 1] = v
        v.t = []

    # Assign triangles to vertices
    for tri in Tarr0:
        tri.v1.t.append(tri)
        tri.v2.t.append(tri)
        tri.v3.t.append(tri)

    return Varr00, Tarr0


def Refine(Varr, Tarr):
    import copy

    Varr0 = copy.copy(Varr)  # Shallow copy is enough here
    Tarr0 = []
    edge_midpoints = {}
    
    lV = len(Varr0)
    next_vid = lV + 1  # Start new Vertex IDs from here
    next_tid = 1       # Start new Triangle IDs from 1

    def get_or_create_midpoint(v1, v2):
        nonlocal next_vid
        key = tuple(sorted((v1.id, v2.id)))
        if key in edge_midpoints:
            return edge_midpoints[key]
        
        midpoint_coords = (v1.coords + v2.coords) / 2
        vm = Vertex(next_vid, midpoint_coords)
        next_vid += 1

        # Boundary condition interpolation
        if v1.bound != v2.bound:
            vm.bound = None
        else:
            vm.bound = v1.bound
            if vm.bound is not None and vm.bound.type == "d":
                vm.coords = Vmove(vm.coords, vm.bound)

        Varr0.append(vm)
        edge_midpoints[key] = vm
        return vm

    for t in Tarr:
        m12 = get_or_create_midpoint(t.v1, t.v2)
        m23 = get_or_create_midpoint(t.v2, t.v3)
        m13 = get_or_create_midpoint(t.v1, t.v3)

        # Create 4 new triangles
        Tarr0.append(Triangle(next_tid, m12, t.v1, m13)); next_tid += 1
        Tarr0.append(Triangle(next_tid, t.v2, m12, m23)); next_tid += 1
        Tarr0.append(Triangle(next_tid, m12, m13, m23)); next_tid += 1
        Tarr0.append(Triangle(next_tid, m23, m13, t.v3)); next_tid += 1

    # Reset triangle references in vertices
    for v in Varr0:
        v.t = []

    for tri in Tarr0:
        tri.v1.t.append(tri)
        tri.v2.t.append(tri)
        tri.v3.t.append(tri)

    return Varr0, Tarr0


def mesh_generator(B_Cond_Arr):
    n_ini=4
    Vb1 = [ Vertex(i+1, B_Cond_Arr[0].shape_f(i*2*np.pi/n_ini) ) for i in range(n_ini)]
    for i in Vb1:
        i.bound = B_Cond_Arr[0]
    Vb2 = [ Vertex(i+1+len(Vb1), B_Cond_Arr[1].shape_f(i*2*np.pi/n_ini) ) for i in range(n_ini)]
    for i in Vb2:
        i.bound = B_Cond_Arr[1]
    Vint = [ Vertex(i+1+len(Vb1)+len(Vb2), (Vb1[i].coords + Vb2[i].coords)/2  ) for i in range(n_ini)]
    Varr = Vb1 + Vb2 + Vint

    Tarr=[]
    Tb1 = [ Triangle(i+1+len(Tarr), Vb1[i], Vb1[i+1],Vint[i]) for i in range(n_ini-1)]
    Tb1 = Tb1 + [ Triangle(1+len(Tb1), Vb1[len(Vint)-1], Vb1[0],Vint[len(Vint)-1])]
    Tarr = Tb1
    Tb2 = [ Triangle(i+1+len(Tarr), Vb2[i], Vb2[i+1],Vint[i]) for i in range(n_ini-1)]
    Tb2 = Tb2 + [ Triangle(1+len(Tarr)+len(Tb2), Vb2[len(Vint)-1], Vb2[0],Vint[len(Vint)-1])]
    Tarr = Tarr + Tb2
    Tint = [ Triangle(i+1+len(Tarr), Vb1[i+1], Vint[i+1], Vint[i]) for i in range(n_ini-1)]
    Tint = Tint + [ Triangle(1+len(Tarr)+len(Tint), Vb1[0], Vint[0], Vint[len(Vint)-1]) ]
    Tint = Tint + [ Triangle(i+1+len(Tarr)+len(Tint), Vb2[i+1], Vint[i+1], Vint[i]) for i in range(n_ini-1)] 
    Tint = Tint + [ Triangle(1+len(Tarr)+len(Tint), Vb2[0], Vint[0], Vint[len(Vint)-1]) ]

    Tarr = Tb1 + Tb2 + Tint
    Varr, Tarr = Refine(Varr,Tarr)
    return Varr, Tarr

from vedo import IcoSphere

def mesh_func_generator(func):
    # Create 3D grid
    #x = np.linspace(-2.5, 2.5, 10)
    #y = np.linspace(-2.5, 2.5, 10)
    #z = np.linspace(-3.5, 2.5, 10)
    X, Y, Z = np.ogrid[-1.:1.:20j, -1.:1.:20j, -1.:1.:20j]
    #X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    volume = func(X,Y,Z)
    verts, faces, _, _ = measure.marching_cubes(volume, level=0)
    sf =  IcoSphere(subdivisions=4)
    verts, faces = sf.vertices,sf.cells

    Varr = [Vertex(i+1, coord) for i, coord in enumerate(verts)] # Convert to Vertex objects

    Tarr = []   
    Earr = []  
    unique_edges = []                                  # Convert to Triangle objects
    for i, (i1, i2, i3) in enumerate(faces):
        v1 = Varr[i1]
        v2 = Varr[i2]
        v3 = Varr[i3]
        Tarr.append(Triangle(i+1, v1, v2, v3))
    for tri in Tarr:                
        tri.v1.t.append(tri)
        tri.v2.t.append(tri)
        tri.v3.t.append(tri)
        if {tri.v1,tri.v2} in unique_edges:
            pos = unique_edges.index({tri.v1,tri.v2})
            Earr[pos].t.append(tri.id)
        else:
            unique_edges.append({tri.v1,tri.v2})
            Earr.append(Edge(len(Earr)+1,tri.v1,tri.v2))
            Earr[len(Earr)-1].t.append(tri.id)
        if {tri.v1,tri.v3} in unique_edges:
            pos = unique_edges.index({tri.v1,tri.v3})
            Earr[pos].t.append(tri.id)
        else:
            unique_edges.append({tri.v1,tri.v3})
            Earr.append(Edge(len(Earr)+1,tri.v1,tri.v3))
            Earr[len(Earr)-1].t.append(tri.id)
        if {tri.v2,tri.v3} in unique_edges:
            pos = unique_edges.index({tri.v2,tri.v3})
            Earr[pos].t.append(tri.id)
        else:
            unique_edges.append({tri.v2,tri.v3})
            Earr.append(Edge(len(Earr)+1,tri.v2,tri.v3))
            Earr[len(Earr)-1].t.append(tri.id)
        
    for e in Earr:
        e.mp = (e.v1.coords + e.v2.coords)/2
    return Varr, Tarr, Earr



# def mesh_func_generator(func):
#     # Create 3D grid
#     x = np.linspace(-2.5, 2.5, 20)
#     y = np.linspace(-2.5, 2.5, 20)
#     z = np.linspace(-3.5, 2.5, 20)
#     X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
#     volume = func(X, Y, Z)
#     verts, faces, _, _ = measure.marching_cubes(volume, level=0.5)

#     # Deduplicate vertices based on coordinates
#     unique_verts = {}
#     Varr = []
#     for coord in verts:
#         key = tuple(np.round(coord, decimals=8))
#         if key not in unique_verts:
#             vertex = Vertex(len(Varr) + 1, coord)
#             unique_verts[key] = vertex
#             Varr.append(vertex)

#     # Map original indices to deduplicated Vertex objects
#     index_map = {i: unique_verts[tuple(np.round(coord, 8))] for i, coord in enumerate(verts)}

#     # Deduplicate triangles based on vertex ids
#     Tarr = []
#     unique_tris = set()
#     for i1, i2, i3 in faces:
#         v1 = index_map[i1]
#         v2 = index_map[i2]
#         v3 = index_map[i3]

#         # Use sorted coordinate tuples as key, rounded for numerical stability
#         coords = sorted([tuple(np.round(v1.coords, 8)),
#                          tuple(np.round(v2.coords, 8)),
#                          tuple(np.round(v3.coords, 8))])
#         tri_key = tuple(coords)

#         if tri_key not in unique_tris:
#             tri = Triangle(len(Tarr) + 1, v1, v2, v3)
#             Tarr.append(tri)
#             unique_tris.add(tri_key)

#     # Optional: link triangles to vertices
#     for tri in Tarr:
#         tri.v1.t.append(tri)
#         tri.v2.t.append(tri)
#         tri.v3.t.append(tri)

#     return Varr, Tarr




def plot_triangles(triangles):
    # Step 1: Gather all unique vertices
    vertex_dict = {}
    index = 0
    for tri in triangles:
        for v in [tri.v1, tri.v2, tri.v3]:
            if v.id not in vertex_dict:
                vertex_dict[v.id] = (index, v.coords)
                index += 1

    # Step 2: Create vertex coordinate arrays
    sorted_vertices = sorted(vertex_dict.items(), key=lambda item: item[1][0])  # sort by index
    coords = np.array([v[1][1] for v in sorted_vertices])
    x, y, z = coords[:, 0], coords[:, 1], coords[:, 2]

    # Step 3: Extract triangle indices
    i, j, k = [], [], []
    id_to_index = {v_id: idx for v_id, (idx, _) in vertex_dict.items()}
    for tri in triangles:
        i.append(id_to_index[tri.v1.id])
        j.append(id_to_index[tri.v2.id])
        k.append(id_to_index[tri.v3.id])

    # Step 4: Plot using Mesh3d
    fig = go.Figure(data=[go.Mesh3d(
        x=x, y=y, z=z,
        i=i, j=j, k=k,
        color='lightblue', opacity=0.5,
        flatshading=True
    )])

    #step 5: Normal vectors

    # Step 6: Add triangle edges as lines
    edge_x = []
    edge_y = []
    edge_z = []

    for tri in triangles:
        verts = [tri.v1.coords, tri.v2.coords, tri.v3.coords, tri.v1.coords]  # loop back to first vertex
        for a, b in zip(verts[:-1], verts[1:]):
            edge_x += [a[0], b[0], None]
            edge_y += [a[1], b[1], None]
            edge_z += [a[2], b[2], None]

    fig.add_trace(go.Scatter3d(
        x=edge_x, y=edge_y, z=edge_z,
        mode='lines',
        line=dict(color='grey', width=2),
        showlegend=False
    ))

    fig.update_layout(scene=dict(
        xaxis_title='X',
        yaxis_title='Y',
        zaxis_title='Z'
    ))
    fig.show()

def plot_triangles_canvas(canvas, triangles):
    # Step 1: Gather all unique vertices
    vertex_dict = {}
    index = 0
    for tri in triangles:
        for v in [tri.v1, tri.v2, tri.v3]:
            if v.id not in vertex_dict:
                vertex_dict[v.id] = (index, v.coords)
                index += 1

    # Step 2: Create vertex coordinate arrays
    sorted_vertices = sorted(vertex_dict.items(), key=lambda item: item[1][0])  # sort by index
    coords = np.array([v[1][1] for v in sorted_vertices])
    x, y, z = coords[:, 0], coords[:, 1], coords[:, 2]

    # Step 3: Extract triangle indices
    i, j, k = [], [], []
    id_to_index = {v_id: idx for v_id, (idx, _) in vertex_dict.items()}
    for tri in triangles:
        i.append(id_to_index[tri.v1.id])
        j.append(id_to_index[tri.v2.id])
        k.append(id_to_index[tri.v3.id])

    # Step 4: Plot using the canvas
    ax = canvas.axes  # Access the axes of the canvas
    ax.clear()  # Clear the previous plot

    # Plot the mesh triangles
    ax.plot_trisurf(x, y, z, triangles=list(zip(i, j, k)), color='lightblue', alpha=0.5)

    # Plot triangle edges
    for tri in triangles:
        verts = [tri.v1.coords, tri.v2.coords, tri.v3.coords, tri.v1.coords]  # loop back to first vertex
        verts = np.array(verts)
        ax.plot(verts[:, 0], verts[:, 1], verts[:, 2], color='k', lw=2)

    # Set axis labels (optional)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    canvas.draw()  # Redraw the canvas





def plot_triangles_vedo(triangles, qt_widget=None):
    vertex_dict = {}
    index = 0
    for tri in triangles:
        for v in [tri.v1, tri.v2, tri.v3]:
            if v.id not in vertex_dict:
                vertex_dict[v.id] = (index, v.coords)
                index += 1

    sorted_vertices = sorted(vertex_dict.items(), key=lambda item: item[1][0])
    points = np.array([v[1][1] for v in sorted_vertices])

    faces = []
    id_to_index = {v_id: idx for v_id, (idx, _) in vertex_dict.items()}
    scalars_per_face = []

    for tri in triangles:
        i1 = id_to_index[tri.v1.id]
        i2 = id_to_index[tri.v2.id]
        i3 = id_to_index[tri.v3.id]
        faces.append([i1, i2, i3])
        scalars_per_face.append(tri.scalar)

    faces = np.array(faces)
    faces_with_size = np.c_[np.full(len(faces), 3), faces]

    mesh2 = Mesh([points, faces_with_size])
    mesh2.celldata["scalar"] = scalars_per_face  # bind scalar per triangle
    mesh2.cmap("viridis", scalars_per_face).add_scalarbar()
    mesh2.lighting("plastic")

    plt = Plotter(qt_widget=qt_widget)
    plt.show(mesh2, axes=1)

    return plt



def plot_triangles_vedo_old(triangles, qt_widget=None):
    vertex_dict = {}
    index = 0
    for tri in triangles:
        for v in [tri.v1, tri.v2, tri.v3]:
            if v.id not in vertex_dict:
                vertex_dict[v.id] = (index, v.coords)
                index += 1

    sorted_vertices = sorted(vertex_dict.items(), key=lambda item: item[1][0])
    points = np.array([v[1][1] for v in sorted_vertices])

    faces = []
    id_to_index = {v_id: idx for v_id, (idx, _) in vertex_dict.items()}
    for tri in triangles:
        i1 = id_to_index[tri.v1.id]
        i2 = id_to_index[tri.v2.id]
        i3 = id_to_index[tri.v3.id]
        faces.append([i1, i2, i3])

    faces = np.array(faces)
    faces_with_size = np.c_[np.full(len(faces), 3), faces]

    mesh2 = Mesh([points, faces_with_size])

    mesh2.color("lightblue").alpha(0.5).lighting("plastic")

    # Create a Plotter and display the sphere mesh
    plt = Plotter()
    plt.show(mesh2, axes=1)

    return plt


  
