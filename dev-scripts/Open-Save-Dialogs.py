__author__ = 'wilsonp'

import sys


from PyQt5.QtWidgets import (QApplication, QCheckBox, QDialog,
        QErrorMessage, QFileDialog, QFrame, QGridLayout,
        QLabel, QMessageBox, QPushButton)


class Dialog(QDialog):

    def __init__(self, parent=None):
        super(Dialog, self).__init__(parent)

        self.openFilesPath = ''

        self.errorMessageDialog = QErrorMessage(self)

        frameStyle = QFrame.Sunken | QFrame.Panel

        self.openFileNameLabel = QLabel()
        self.openFileNameLabel.setFrameStyle(frameStyle)
        self.openFileNameButton = QPushButton("Open()")

        self.saveFileNameLabel = QLabel()
        self.saveFileNameLabel.setFrameStyle(frameStyle)
        self.saveFileNameButton = QPushButton("Save()")


        self.openFileNameButton.clicked.connect(self.setOpenFileName)
        self.saveFileNameButton.clicked.connect(self.setSaveFileName)

        self.native = QCheckBox()
        self.native.setText("Use native file dialog.")
        self.native.setChecked(True)
        if sys.platform not in ("win32", "darwin"):
            self.native.hide()

        layout = QGridLayout()
        layout.setColumnStretch(1, 1)
        layout.setColumnMinimumWidth(1, 250)

        # layout.addWidget(self.directoryLabel, 0, 0)
        layout.addWidget(self.openFileNameButton, 1, 0)
        layout.addWidget(self.openFileNameLabel, 1, 1)
        layout.addWidget(self.saveFileNameButton, 2, 0)
        layout.addWidget(self.saveFileNameLabel, 2, 1)

        layout.addWidget(self.native, 15, 0)
        self.setLayout(layout)

        self.setWindowTitle("Open Save File")



    def setOpenFileName(self):
        options = QFileDialog.Options()
        if not self.native.isChecked():
            options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,
                "QFileDialog.getOpenFileName()", self.openFileNameLabel.text(),
                "All Files (*);;Text Files (*.txt)", options=options)
        if fileName:
            self.openFileNameLabel.setText(fileName)


    def setSaveFileName(self):
        options = QFileDialog.Options()
        if not self.native.isChecked():
            options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self,
                "QFileDialog.getSaveFileName()",
                self.saveFileNameLabel.text(),
                "All Files (*);;Text Files (*.txt)", options=options)
        if fileName:
            self.saveFileNameLabel.setText(fileName)




if __name__ == '__main__':
    app = QApplication(sys.argv)
    dialog = Dialog()
    dialog.show()
    sys.exit(app.exec_())