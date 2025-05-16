from PyQt5.QtWidgets import QApplication, QVBoxLayout, QWidget
from vedo import Plotter, Sphere  # Import Sphere from vedo

from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import sys

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

        # Create a sphere mesh and display it
        self.mesh = Sphere(r=1, c="lightgray")  # Create a sphere with radius 1 and color lightgray
        self.plt.show(self.mesh)

    def load_mesh(self, filepath):
        # Optionally: Load a mesh from file (if needed in the future)
        pass

if __name__ == "__main__":

    # Create the Qt application
    app = QApplication(sys.argv)

    # Create and show the main widget
    viewer = VedoWidget()
    viewer.resize(800, 600)
    viewer.show()

    # Start the Qt event loop
    sys.exit(app.exec_())
