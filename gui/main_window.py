from PyQt5.QtWidgets import QWidget, QVBoxLayout, QStackedLayout, QPushButton
from tasks.Surface_Minimization import Surface_Minimizer
from tasks.Charge_Minimization import Charge_Minimizer

class MainWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Task Launcher")
        self.resize(800, 850)

        self.stack = QStackedLayout()
        self.setLayout(self.stack)

        self.task_choice_page = QWidget()
        task_choice_layout = QVBoxLayout(self.task_choice_page)

        # Buttons for tasks
        calc_btn = QPushButton("Surface Minimization")
        reverse_btn = QPushButton("Charge Optimization")

        calc_btn.clicked.connect(lambda: self.load_task_page("Surf_Minimization"))
        reverse_btn.clicked.connect(lambda: self.load_task_page("Charge_Minimization"))

        task_choice_layout.addWidget(calc_btn)
        task_choice_layout.addWidget(reverse_btn)

        # Pages dictionary
        self.pages = {
            "TaskChoice": self.task_choice_page,
            "Surf_Minimization": self.create_task_page(Surface_Minimizer()),
            "Charge_Minimization": self.create_task_page(Charge_Minimizer())
        }

        # Add all pages to stack
        for page in self.pages.values():
            self.stack.addWidget(page)

        self.stack.setCurrentWidget(self.task_choice_page)

    def create_task_page(self, task_obj):
        """Creates a page for the given task with a back button"""
        page = QWidget()
        layout = QVBoxLayout(page)

        task_widget = task_obj.get_widget(self)
        back_button = QPushButton("Back")
        back_button.clicked.connect(lambda: self.stack.setCurrentWidget(self.pages["TaskChoice"]))

        layout.addWidget(task_widget)
        layout.addWidget(back_button)

        return page

    def load_task_page(self, task_name):
        self.stack.setCurrentWidget(self.pages[task_name])
