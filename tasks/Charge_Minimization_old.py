from PyQt5.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QLineEdit, QPushButton, QLabel, QSizePolicy
import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from models.mesh import mesh_func_generator, plot_triangles_canvas
from models.surf_min import Minimization, H_Min
from .task_base import Task
from vedo import Plotter

class MplCanvas(FigureCanvas):
    def __init__(self, parent=None):
        fig = Figure()
        self.axes = fig.add_subplot(111, projection='3d')
        super().__init__(fig)
        self.setParent(parent)

class VedoCanvas:
    def __init__(self, parent=None):
        self.plotter = Plotter(qtWidget=parent, bg="white", size=(800, 600))
    
    def clear(self):
        self.plotter.clear()

    def show_mesh(self, mesh):
        self.plotter.clear()
        self.plotter.show(mesh, axes=1, interactive=False)

    def add_actor(self, actor):
        self.plotter.show(actor, resetcam=False, interactive=False)





class Charge_Minimizer(Task):
    def get_widget(self, parent=None):
        widget = QWidget(parent)
        main_layout = QVBoxLayout(widget)

        main_layout.addWidget(QLabel("Insert the function defining the surface in the implicit form (e.g. x**2 + y**2 + z**2 -1)"))


        funcs = []
        mesh_data = {"V": [], "T": []}


        # Bound 1 row
        row1 = QHBoxLayout()
        store_imp_fun = QPushButton("Store Function")
        store_imp_fun.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        store_imp_fun_input = QLineEdit()
        store_imp_fun.clicked.connect(lambda: store_function_impl_list(store_imp_fun_input, funcs))
        row1.addWidget(store_imp_fun)
        row1.addWidget(store_imp_fun_input)
        main_layout.addLayout(row1)


        # Mesh control row
        row3 = QHBoxLayout()
        b_gen = QPushButton("Generate mesh")
        b_gen.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        clearm_button = QPushButton("Clear mesh")
        clearm_button.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        b_refine = QPushButton("Refine mesh")
        b_refine.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        row3.addWidget(b_gen)
        row3.addWidget(clearm_button)
        row3.addWidget(b_refine)
        main_layout.addLayout(row3)

        # Plot area
        self.canvas = MplCanvas(widget)
        main_layout.addWidget(self.canvas)

       

        # Minimization control row
        row4 = QHBoxLayout()
        b_grad = QPushButton("Minimization: gradient")
        b_grad.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        b_hess = QPushButton("Minimization: hessian")
        b_hess.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        row4.addWidget(b_grad)
        row4.addWidget(b_hess)
        main_layout.addLayout(row4)

        # Back or result row
        result_label = QLabel("Result: ")
        main_layout.addWidget(result_label)

        # Functionality
        def generate_func_mesh():
            V, T = mesh_func_generator(funcs[0])
            mesh_data["V"], mesh_data["T"] = V, T
            plot_triangles_canvas(self.canvas, T)

        # def refine_mesh():
        #     V_new, T_new = Refine(mesh_data["V"], mesh_data["T"])
        #     mesh_data["V"], mesh_data["T"] = V_new, T_new
        #     plot_triangles_canvas(self.canvas, T_new)

        def grad_mesh():
            Minimization(mesh_data["V"])
            plot_triangles_canvas(self.canvas, mesh_data["T"])

        def hess_mesh():
            H_Min(mesh_data["V"])
            plot_triangles_canvas(self.canvas, mesh_data["T"])

        def clear_func():
            funcs.clear()

        def store_function_impl_list(entr, funcs):
            input_text = entr.text().strip()
            try:
                func = eval(f"lambda x, y, z: {input_text}", {"np": np, "x": 0, "y": 0, "z": 0})
                funcs.clear()        # Clear the existing list
                funcs.append(func)   # Add the new lambda to it
            except Exception as e:
                print(f"Error parsing input: {e}")

        # Connect buttons
        b_gen.clicked.connect(generate_func_mesh)
        #b_refine.clicked.connect(refine_mesh)
        clearm_button.clicked.connect(clear_func)
        b_grad.clicked.connect(grad_mesh)
        b_hess.clicked.connect(hess_mesh)

        return widget
