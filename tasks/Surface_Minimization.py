from PyQt5.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QLineEdit, QPushButton, QLabel, QSizePolicy
import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from models.mesh import B_Cond, mesh_generator, Refine, plot_triangles_canvas
from models.surf_min import Minimization, H_Min
from .task_base import Task

class MplCanvas(FigureCanvas):
    def __init__(self, parent=None):
        fig = Figure()
        self.axes = fig.add_subplot(111, projection='3d')
        super().__init__(fig)
        self.setParent(parent)

class Surface_Minimizer(Task):
    def get_widget(self, parent=None):
        widget = QWidget(parent)
        main_layout = QVBoxLayout(widget)

        main_layout.addWidget(QLabel("Insert the boundary constraints as parametric functions (e.g. [np.cos(s),np.sin(s), z ])"))

        bounds = []
        bounds_list = []
        B_Cond_Arr = []
        mesh_data = {"V": [], "T": []}


        # Bound 1 row
        row1 = QHBoxLayout()
        store_button1 = QPushButton("Store Bound 1")
        store_button1.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        bound_input1 = QLineEdit()
        bounds.append(bound_input1)
        store_button1.clicked.connect(lambda: store_function_list(bound_input1, bounds_list))
        row1.addWidget(store_button1)
        row1.addWidget(bound_input1)
        main_layout.addLayout(row1)

        # Bound 2 row
        row2 = QHBoxLayout()
        store_button2 = QPushButton("Store Bound 2")
        store_button2.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        bound_input2 = QLineEdit()
        bounds.append(bound_input2)
        store_button2.clicked.connect(lambda: store_function_list(bound_input2, bounds_list))
        row2.addWidget(store_button2)
        row2.addWidget(bound_input2)
        main_layout.addLayout(row2)

        # Mesh control row
        row3 = QHBoxLayout()
        b_gen = QPushButton("Generate mesh")
        b_gen.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        clearb_button = QPushButton("Clear bounds")
        clearb_button.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        b_refine = QPushButton("Refine mesh")
        b_refine.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        row3.addWidget(b_gen)
        row3.addWidget(clearb_button)
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
        def generate_mesh():
            V, T = mesh_generator(B_Cond_Arr)
            mesh_data["V"], mesh_data["T"] = V, T
            plot_triangles_canvas(self.canvas, T)

        def refine_mesh():
            V_new, T_new = Refine(mesh_data["V"], mesh_data["T"])
            mesh_data["V"], mesh_data["T"] = V_new, T_new
            plot_triangles_canvas(self.canvas, T_new)

        def grad_mesh():
            Minimization(mesh_data["V"])
            plot_triangles_canvas(self.canvas, mesh_data["T"])

        def hess_mesh():
            H_Min(mesh_data["V"])
            plot_triangles_canvas(self.canvas, mesh_data["T"])

        def clear_bound():
            B_Cond_Arr.clear()

        def store_function_list(entr, list_expr):
            list_expr.clear()
            input_text = entr.text().strip()
            if not input_text.startswith("["):
                input_text = f"[{input_text}]"
            try:
                raw_items = input_text.strip("[]").split(",")
                cleaned_items = [item.strip() for item in raw_items if item.strip()]
                for expr in cleaned_items:
                    list_expr.append(expr)
                funcs = [eval(f"lambda s: {e}", {"np": np, "s": 0}) for e in list_expr]
                B_Cond_Arr.append(B_Cond(len(B_Cond_Arr)+1, funcs, "d"))
            except Exception as e:
                print(f"Error parsing input: {e}")

        # Connect buttons
        b_gen.clicked.connect(generate_mesh)
        b_refine.clicked.connect(refine_mesh)
        clearb_button.clicked.connect(clear_bound)
        b_grad.clicked.connect(grad_mesh)
        b_hess.clicked.connect(hess_mesh)

        return widget
