__author__ = 'wilsonp'
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QApplication, QWidget

from ui_calculatorform import Ui_CalculatorForm


class CalculatorForm(QWidget):
    def __init__(self, parent=None):
        super(CalculatorForm, self).__init__(parent)

        self.ui = Ui_CalculatorForm()

        self.ui.setupUi(self)

    @pyqtSlot(int)
    def on_inputSpinBox1_valueChanged(self, value):
        self.ui.outputWidget.setText(str(value + self.ui.inputSpinBox2.value()))

    @pyqtSlot(int)
    def on_inputSpinBox2_valueChanged(self, value):
        self.ui.outputWidget.setText(str(value + self.ui.inputSpinBox1.value()))


if __name__ == '__main__':
    import sys

    app = QApplication(sys.argv)
    calculator = CalculatorForm()
    calculator.show()
    sys.exit(app.exec_())