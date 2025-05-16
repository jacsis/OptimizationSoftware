from PyQt5.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QLineEdit, QPushButton, QLabel, QSizePolicy
import vedo
from vedo import Plotter, Mesh, Sphere  # Import Mesh from vedo
from models.mesh import mesh_func_generator, Triangle  # Ensure Triangle is imported correctly
from models.electrostatic import Minimization, H_Min, A_tot, Minimization, charg_dens_solv
from .task_base import Task
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import numpy as np

class VedoWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.layout = QVBoxLayout()
        self.setLayout(self.layout)

        # Set up the VTK widget for rendering
        self.vtkWidget = QVTKRenderWindowInteractor(self)
        self.layout.addWidget(self.vtkWidget)

        # Create a Vedo Plotter using the VTK widget
        self.plt = Plotter(qt_widget=self.vtkWidget, bg="white")

        # Create a default sphere mesh
        self.mesh = Sphere(r=1, c="lightgray")  # Create a sphere with radius 1 and color lightgray
        self.plt.show(self.mesh)

    def clear_plot(self):
        if self.plt:
            self.plt.clear()

    def show_mesh(self, mesh_data):
        if not mesh_data:
            return

        if isinstance(mesh_data, list) and all(isinstance(t, Triangle) for t in mesh_data):
            # Extract vertices and faces from list of Triangle objects
            vertices = []
            faces = []
            vertex_index = 0
            vertex_map = {}
            scalars = []  # List to store the scalar values for each triangle

            for tri in mesh_data:
                v = [vertices.coords for vertices in tri.get_vertices()]  # Get vertices for the triangle
                face = []
                for vert in v:
                    key = tuple(vert)
                    if key not in vertex_map:
                        vertex_map[key] = vertex_index
                        vertices.append(vert)
                        vertex_index += 1
                    face.append(vertex_map[key])
                faces.append(face)
                scalars.append(tri.scalar)  # Collect the scalar value for each triangle

            vertices = np.array(vertices)
            faces = np.array(faces)

            # Normalize the scalar values to the range [0, 1]
            scalar_min = min(scalars)
            scalar_max = max(scalars)+0.0000000000001
            scalars_normalized = [(scalar - scalar_min+0.0000000000001) / (scalar_max - scalar_min) for scalar in scalars]
            

            # Apply colormap based on the normalized scalar values
            #colormap = 'viridis'  # You can change this to any colormap you like (e.g., 'plasma', 'inferno', etc.)
            #colors = [vedo.colors.color_map(scalar, colormap, vmin=0, vmax=1) for scalar in scalars_normalized]
            #colors = [[int(round(scalar*255)),0,0,255] for scalar in scalars_normalized]
            colors = [[int(np.round(scalar*255)), 255-int(np.round(scalar*255)), 0, 255] for scalar in scalars_normalized]
            

            # Create the mesh with vertices and faces, and apply the colors
            vedo_mesh = Mesh([vertices, faces])
            vedo_mesh.cellcolors=colors

            # Show the mesh in the plotter
            self.plt.show(vedo_mesh, resetcam=True)

        else:
            self.plt.show(mesh_data, resetcam=True)

    # def show_mesh(self, mesh_data):
    #     if not mesh_data:
    #         return

    #     if isinstance(mesh_data, list) and all(isinstance(t, Triangle) for t in mesh_data):
    #         # Extract vertices and faces from list of Triangle objects
    #         vertices = []
    #         faces = []
    #         vertex_index = 0
    #         vertex_map = {}

    #         for tri in mesh_data:
    #             v = [vertices.coords for vertices in tri.get_vertices()]  # Should return list or array of 3 vertices
    #             face = []
    #             for vert in v:
    #                 key = tuple(vert)
    #                 if key not in vertex_map:
    #                     vertex_map[key] = vertex_index
    #                     vertices.append(vert)
    #                     vertex_index += 1
    #                 face.append(vertex_map[key])
    #             faces.append(face)

    #         vertices = np.array(vertices)
    #         faces = np.array(faces)
    #         vedo_mesh = Mesh([vertices, faces]).c("lightblue").alpha(0.5)
    #         self.plt.show(vedo_mesh, resetcam=True)

    #     else:
    #         self.plt.show(mesh_data, resetcam=True)

lamb = 1.
class Charge_Minimizer(Task):
    def get_widget(self, parent=None):
        widget = QWidget(parent)
        main_layout = QVBoxLayout(widget)

        main_layout.addWidget(QLabel("Insert the function defining the surface in the implicit form (e.g. x**2 + y**2 + z**2 -1)"))

        funcs = []
        mesh_data = {"V": [], "T": [], "E" :[]}
        Qtot = 1.
        

        # Function input row
        row1 = QHBoxLayout()
        store_imp_fun = QPushButton("Store Function")
        store_imp_fun.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        store_imp_fun_input = QLineEdit()
        store_imp_fun.clicked.connect(lambda: store_function_impl_list(store_imp_fun_input, funcs))
        row1.addWidget(store_imp_fun)
        row1.addWidget(store_imp_fun_input)
        main_layout.addLayout(row1)

        # Chage input row 
        row2 = QHBoxLayout()
        main_layout.addWidget(QLabel("Insert total charge"))
        Qtot_input = QLineEdit()
        row2.addWidget(Qtot_input)
        main_layout.addLayout(row2)
       

        # Mesh control row
        row3 = QHBoxLayout()
        b_gen = QPushButton("Generate mesh")
        b_gen.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        clearm_button = QPushButton("Clear mesh")
        clearm_button.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        row3.addWidget(b_gen)
        row3.addWidget(clearm_button)
        main_layout.addLayout(row3)

        # Vedo plot area (using the VedoWidget class)
        self.plotter_widget = VedoWidget()  # Create VedoWidget instance
        main_layout.addWidget(self.plotter_widget)

        # Solver for charge
        row4 = QHBoxLayout()
        b_solv = QPushButton("Solve for charge")
        b_solv.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        b_hess = QPushButton("Minimization: hessian")
        b_hess.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        row4.addWidget(b_solv)
        row4.addWidget(b_hess)
        main_layout.addLayout(row4)

        # # Minimization control row
        # row4 = QHBoxLayout()
        # b_grad = QPushButton("Minimization: gradient")
        # b_grad.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        # b_hess = QPushButton("Minimization: hessian")
        # b_hess.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        # row4.addWidget(b_grad)
        # row4.addWidget(b_hess)
        # main_layout.addLayout(row4)

        # Result label
        result_label = QLabel("Result: ")
        main_layout.addWidget(result_label)

        # Functionality
        def generate_func_mesh():
            if not funcs:
                print("No function stored.")
                return
            V, T, E = mesh_func_generator(funcs[0])
            mesh_data["V"], mesh_data["T"], mesh_data["E"] = V, T, E
            input_text = Qtot_input.text().strip()
            try:
                Qtot = float(input_text)
            except ValueError:
                print("Invalid total charge input.")
                return
            dens = Qtot / A_tot(mesh_data["T"])
            for tri in mesh_data["T"]:
                tri.scalar = dens
            self.plotter_widget.clear_plot()
            self.plotter_widget.show_mesh(T)
          
        def charg_solver():
            if not mesh_data["T"]:
                return
            charg_dens_solv(mesh_data["E"],mesh_data["T"], Qtot)
            c=0
            a=0
            for t in mesh_data["T"]:
                c+=t.scalar * t.areav
                a+=t.areav
                print(t.areav)
            print(a)
            self.plotter_widget.show_mesh(mesh_data["T"])  

        def grad_mesh():
            if not mesh_data["T"]:
                return
            global lamb  # if you're using lamb as a global variable
            lamb = Minimization(mesh_data["T"], lamb, Qtot)
            print((lamb, mesh_data["T"][0].scalar))
            self.plotter_widget.show_mesh(mesh_data["T"])
            

        def hess_mesh():
            if not mesh_data["V"]:
                return
            H_Min(mesh_data["V"])
            self.plotter_widget.show_mesh(mesh_data["T"])

        def clear_func():
            funcs.clear()
            mesh_data["V"].clear()
            mesh_data["T"].clear()
            self.plotter_widget.clear_plot()

        def store_function_impl_list(entr, funcs):
            input_text = entr.text().strip()
            try:
                func = eval(f"lambda x, y, z: {input_text}", {"np": np, "x": 0, "y": 0, "z": 0})
                funcs.clear()
                funcs.append(func)
                print("Function stored.")
            except Exception as e:
                print(f"Error parsing input: {e}")

        # Connect buttons
        b_gen.clicked.connect(generate_func_mesh)
        clearm_button.clicked.connect(clear_func)
        b_solv.clicked.connect(charg_solver)
        b_hess.clicked.connect(hess_mesh)

        return widget
