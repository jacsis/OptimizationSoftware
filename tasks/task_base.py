from abc import ABC, abstractmethod
from PyQt5.QtWidgets import QWidget

class Task(ABC):
    @abstractmethod
    def get_widget(self, parent=None) -> QWidget:
        """Returns a QWidget for the task."""
        pass