from PyQt5.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QLineEdit, QPushButton, QLabel, QSizePolicy, QFrame
import numpy as np
from models.mesh import mesh_func_generator, plot_triangles_vedo
from models.surf_min import Minimization, H_Min
from .task_base import Task

class Charge_Minimizer(Task):
    def get_widget(self, parent=None):
        widget = QWidget(parent)
        main_layout = QVBoxLayout(widget)

        main_layout.addWidget(QLabel("Insert the function defining the surface in the implicit form (e.g. x**2 + y**2 + z**2 -1)"))

        funcs = []
        mesh_data = {"V": [], "T": []}

        # Function input row
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
        row3.addWidget(b_gen)
        row3.addWidget(clearm_button)
        main_layout.addLayout(row3)

        # Vedo plot area (QFrame placeholder)
        self.plot_frame = QFrame()
        self.plot_frame.setMinimumSize(800, 600)
        main_layout.addWidget(self.plot_frame)
        self.plotter = None  # keep reference to the Plotter instance

        # Minimization control row
        row4 = QHBoxLayout()
        b_grad = QPushButton("Minimization: gradient")
        b_grad.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        b_hess = QPushButton("Minimization: hessian")
        b_hess.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        row4.addWidget(b_grad)
        row4.addWidget(b_hess)
        main_layout.addLayout(row4)

        # Result label
        result_label = QLabel("Result: ")
        main_layout.addWidget(result_label)

        # Functionality
        def generate_func_mesh():
            if not funcs:
                print("No function stored.")
                return
            V, T = mesh_func_generator(funcs[0])
            mesh_data["V"], mesh_data["T"] = V, T
            self.plotter = plot_triangles_vedo(T, qt_widget=self.plot_frame)

        def grad_mesh():
            if not mesh_data["V"]:
                return
            Minimization(mesh_data["V"])
            self.plotter = plot_triangles_vedo(mesh_data["T"], qt_widget=self.plot_frame)

        def hess_mesh():
            if not mesh_data["V"]:
                return
            H_Min(mesh_data["V"])
            self.plotter = plot_triangles_vedo(mesh_data["T"], qt_widget=self.plot_frame)

        def clear_func():
            funcs.clear()
            mesh_data["V"].clear()
            mesh_data["T"].clear()
            if self.plotter:
                self.plotter.clear()

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
        b_grad.clicked.connect(grad_mesh)
        b_hess.clicked.connect(hess_mesh)

        return widget
