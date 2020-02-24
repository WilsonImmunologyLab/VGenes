__author__ = 'wilsonp'



import sys, os, re


from PyQt5.QtWidgets import (QApplication, QCheckBox, QDialog,
        QErrorMessage, QFileDialog, QFrame, QGridLayout,
        QLabel, QMessageBox, QPushButton)


class Dialog(QDialog):

    def __init__(self, parent=None):
        super(Dialog, self).__init__(parent)

        self.openFilesPath = ''

        self.errorMessageDialog = QErrorMessage(self)

        frameStyle = QFrame.Sunken | QFrame.Panel

        self.DescriptionLabel = QLabel()
        self.DescriptionLabel.setFrameStyle(frameStyle)
        self.DescriptionLabel.setText('This utility will convert IMGT FASTA files to NCBI FASTA files')


        self.openFileNameLabel = QLabel()
        self.openFileNameLabel.setFrameStyle(frameStyle)
        self.openFileNameButton = QPushButton("Open (IMGT FASTA filename):")

        self.saveFileNameLabel = QLabel()
        self.saveFileNameLabel.setFrameStyle(frameStyle)
        self.saveFileNameButton = QPushButton("Save as (NCBI FASTA filename):")

        self.initiateButton = QPushButton("Initiate conversion")

        self.openFileNameButton.clicked.connect(self.setOpenFileName)
        self.saveFileNameButton.clicked.connect(self.setSaveFileName)
        self.initiateButton.clicked.connect(self.setInitiateConversion)

        self.native = QCheckBox()
        self.native.setText("Use native file dialog.")
        self.native.setChecked(True)
        if sys.platform not in ("win32", "darwin"):
            self.native.hide()

        layout = QGridLayout()
        layout.setColumnStretch(1, 1)
        layout.setColumnMinimumWidth(1, 250)

        layout.addWidget(self.DescriptionLabel, 0, 0)
        layout.addWidget(self.openFileNameButton, 1, 0)
        layout.addWidget(self.openFileNameLabel, 1, 1)
        layout.addWidget(self.saveFileNameButton, 2, 0)
        layout.addWidget(self.saveFileNameLabel, 2, 1)
        layout.addWidget(self.initiateButton, 3, 0)

        layout.addWidget(self.native, 4, 0)
        self.setLayout(layout)

        self.setWindowTitle("IMGT to NCBI FASTA converter")



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

    def setInitiateConversion(self):
        print('converted')

        Ofilenamed = self.openFileNameLabel.text()
        SfileNamed = self.saveFileNameLabel.text()

        # (dirname, filename) = os.path.split(SfileNamed)
        # print(dirname)
        # print(filename)

        IMGTFile  = open(Ofilenamed, 'r')
        NCBIfile = open(SfileNamed, 'w')

        NCBILine = ''
        StartLine = True
        for IMGTlines in IMGTFile:

            if IMGTlines[0] == '>':
                if StartLine == False:
                    NCBIfile.write('\n')

                SIMGTLines = IMGTlines.split('|')
                NCBILine = SIMGTLines[0] + '-' + SIMGTLines[1] + '\n'
                NCBILine = NCBILine.replace('/', '-')
                print(NCBILine)
                NCBIfile.write(NCBILine)
                StartLine = False


            else:
                NCBILine = IMGTlines
                NCBILine2 = ''
                NCBILine3 = ''
                SIMGTLines = NCBILine.split('.')  #  top get rid of .
                for seqparts in SIMGTLines:
                    NCBILine2 += seqparts
                NCBILine3 = NCBILine2.replace('\n', '')

                NCBIfile.write(NCBILine3)
                print(NCBILine3)

        NCBIfile.write('\n')



        IMGTFile.close()
        NCBIfile.close()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    dialog = Dialog()
    dialog.show()
    sys.exit(app.exec_())