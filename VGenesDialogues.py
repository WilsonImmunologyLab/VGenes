__author__ = 'wilsonp'
from PyQt5 import QtWidgets, QtGui, QtCore
import os
global LastFileName
LastFileName = ''
def newFile(self):
            global LastFileName
            if LastFileName != '':
                workingdir, filename = os.path.split(LastFileName)
                os.chdir(workingdir)

            options = QtWidgets.QFileDialog.Options()

            # options |= QtWidgets.QFileDialog.DontUseNativeDialog
            DBFilename, _ = QtWidgets.QFileDialog.getSaveFileName(self, "New Database", "New database", "All Files (*);;SqLLite database Files (*.db)", options=options)
            LastFileName = DBFilename


            return DBFilename

def saveFile(self, typeSave):
    import time
    options = QtWidgets.QFileDialog.Options()
    global LastFileName
    if LastFileName != '':
        if len(LastFileName)==1:
            workingdir, filename = os.path.split(LastFileName)
            os.chdir(workingdir)
        elif len(LastFileName)>1:
            workingdir, filename = os.path.split(LastFileName[0])
            os.chdir(workingdir)


    if typeSave == 'db':
        queryIs = "Database " + time.strftime('%c')
        FileTypes  = "VGenes database Files (*.vdb);;All Files (*)"
    elif typeSave == 'Nucleotide':
        queryIs = "FASTA " + time.strftime('%c')
        FileTypes  = "Nucleotide Files (*.nt);;FASTA Files (*.fasta);;Text Files (*.txt);;All Files (*)"
    elif typeSave == 'FASTA':
        queryIs = "FASTA " + time.strftime('%c')
        FileTypes  = "FASTA Files (*.fasta);;Nucleotide Files (*.nt);;Text Files (*.txt);;All Files (*)"
    elif typeSave == 'fastq':
        # queryIs = "What do you want to name the merged fastq file?" + time.strftime('%c')
        filenamed = ''
        # LastName = LastFileName[0]
        workingdir, filename = os.path.split(LastFileName[0])
        cutto = len(filename)-6
        filename = filename[:cutto]
        filenamed = filename + '_x_'
        workingdir, filename = os.path.split(LastFileName[1])
        cutto = len(filename)-6
        filename = filename[:cutto]
        filenamed = filenamed + filename
        queryIs = filenamed + '_Merged'
        FileTypes  = "fastq Files (*.fastq);;fq Files (*.fq);;All Files (*)"

    elif typeSave == 'seq':
        queryIs = "Sequence " + time.strftime('%c')
        FileTypes  = "All Files (*);;Seq Files (*.seq);;Text Files (*.txt)"
    elif typeSave == 'csv' or 'CSV':
        queryIs = "CSV " + time.strftime('%c')
        FileTypes  = "comma separated values (*.csv);;All Files (*)"
    elif typeSave == 'text':

        queryIs = "Text: " + time.strftime('%c')
        FileTypes  = "text (*.txt);;All Files (*)"


    # options |= QtWidgets.QFileDialog.DontUseNativeDialog
    Filename, _ = QtWidgets.QFileDialog.getSaveFileName(self, queryIs, queryIs, FileTypes, options=options)


    if Filename:
        LastFileName = Filename

        return Filename
def openFile(self, typeOpen):
    global LastFileName
    queryIs = ''
    FileTypes = ''
    if LastFileName != '':
        workingdir, filename = os.path.split(LastFileName)
        os.chdir(workingdir)


    if typeOpen == 'db':
        queryIs = "Open Database"
        FileTypes  = "VGenes database Files (*.vdb);;"
    elif typeOpen == 'Nucleotide':
        queryIs = "Open nucleotide file"
        FileTypes = "Nucleotide Files (*.nt);;FASTA Files (*.fasta);;Text Files (*.txt);;All Files (*)"
    elif typeOpen == 'FASTA':
        queryIs = "FASTA "
        FileTypes = "FASTA Files (*.fasta);;Nucleotide Files (*.nt);;Text Files (*.txt);;All Files (*)"
    elif typeOpen == 'seq':
        queryIs = "Open a single sequence file"
        FileTypes  = "All Files (*);;Seq Files (*.seq);;Text Files (*.txt)"
    elif typeOpen == 'CSV':
        queryIs = "Open a Comma Separated Values file"
        FileTypes  = "All Files (*);;CSV Files (*.csv);;Text Files (*.txt)"
    elif typeOpen == 'All':
        queryIs = "All file types"
        FileTypes  = "All Files (*);;"

    options = QtWidgets.QFileDialog.Options()

    # options |= QtWidgets.QFileDialog.DontUseNativeDialog
    fileName, _ = QtWidgets.QFileDialog.getOpenFileName(self,
            queryIs,
            queryIs,
            FileTypes, options=options)

    if fileName:
        LastFileName = fileName
        return fileName


def openFiles(self, typeOpen):
    global LastFileName
    if LastFileName != '' and len(LastFileName) == 1:
        workingdir, filename = os.path.split(LastFileName)
        os.chdir(workingdir)

    if typeOpen == 'db':
        queryIs = "Open Database"
        FileTypes  = "VGenes database Files (*.vdb);;"
    elif typeOpen == 'Nucleotide':
        queryIs = "Open nucleotide file"
        FileTypes = "Nucleotide Files (*.nt);;FASTA Files (*.fasta);;Text Files (*.txt);;All Files (*)"
    elif typeOpen == 'FASTA':
        queryIs = "FASTA "
        FileTypes = "FASTA Files (*.fasta);;Nucleotide Files (*.nt);;Text Files (*.txt);;All Files (*)"
    elif typeOpen == 'seq':
        queryIs = "Open a single sequence file"
        FileTypes  = "All Files (*);;Seq Files (*.seq);;Text Files (*.txt)"
    elif typeOpen == 'CSV':
        queryIs = "Open a Comma Separated Values file"
        FileTypes  = "All Files (*);;CSV Files (*.csv);;Text Files (*.txt)"
    elif typeOpen == 'All':
        queryIs = "All file types"
        FileTypes  = "All Files (*);;"

    # queryIs = "Open single sequence files"
    # FileTypes  = "All Files (*);;Seq Files (*.seq);;Text Files (*.txt)"


    options = QtWidgets.QFileDialog.Options()

    # options |= QtWidgets.QFileDialog.DontUseNativeDialog
    fileName, _ = QtWidgets.QFileDialog.getOpenFileNames(self,
            queryIs,
            queryIs,
            FileTypes, options=options)

    if fileName:
        LastFileName = fileName[0]
        return fileName

def openfastq(self):
    global LastFileName
    if LastFileName != '' and len(LastFileName) == 1:
        workingdir, filename = os.path.split(LastFileName)
        os.chdir(workingdir)

    queryIs = "Choose two paired fastq files to be merged"
    FileTypes  = "fastq Files (*.fastq);;fq Files (*.fq);;Text Files (*.txt);;All Files (*)"


    options = QtWidgets.QFileDialog.Options()

    # options |= QtWidgets.QFileDialog.DontUseNativeDialog
    fileName, _ = QtWidgets.QFileDialog.getOpenFileNames(self,
            queryIs,
            queryIs,
            FileTypes, options=options)

    if fileName:
        LastFileName = fileName
        return fileName

def progressDialogue(self):
    response = QtWidgets.QProgressDialog()


def questionMessage(self, question, buttons):

    type = ''

    if buttons == 'YN':
        type = QtWidgets.QMessageBox.No | QtWidgets.QMessageBox.Yes
    elif buttons == "YNC":
        type  = QtWidgets.QMessageBox.Cancel | QtWidgets.QMessageBox.No | QtWidgets.QMessageBox.Yes # , QtWidgets.QMessageBox.Yes
    elif buttons == "YC":
        type  = QtWidgets.QMessageBox.Cancel | QtWidgets.QMessageBox.Yes # , QtWidgets.QMessageBox.Yes
    elif buttons == "OKC":
        type  = QtWidgets.QMessageBox.Cancel | QtWidgets.QMessageBox.Ok # | QtWidgets.QMessageBox.Yes # , QtWidgets.QMessageBox.Yes
    elif buttons == 'OK':
        type  = QtWidgets.QMessageBox.Ok #  | QtWidgets.QMessageBox.Yes # , QtWidgets.QMessageBox.Yes
    elif buttons == "YNCA":
        type  = QtWidgets.QMessageBox.Cancel | QtWidgets.QMessageBox.No | QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.NoToAll | QtWidgets.QMessageBox.YesToAll# , QtWidgets.QMessageBox.Yes
    # QtWidgets.QMessageBox.setDefaultButton(self, QtWidgets.QMessageBox.Yes)

    reply = QtWidgets.QMessageBox.question(self, "QMessageBox.question()", question, type)





    if reply == QtWidgets.QMessageBox.Yes:
        return 'Yes'
    elif reply == QtWidgets.QMessageBox.No:
        return 'No'
    elif reply == QtWidgets.QMessageBox.Ok:
        return 'OK'
    elif reply == QtWidgets.QMessageBox.NoToAll:
        return 'NoToAll'
    elif reply == QtWidgets.QMessageBox.YesToAll:
        return 'YesToAll'

    else:
        return 'Cancel'

def informationMessage(self, statement, buttons):

    #
    # reply = QtWidgets.QMessageBox.information(self,
    #         "QMessageBox.information()", statement)
    # if reply == QtWidgets.QMessageBox.Ok:
    #     self.informationLabel.setText("OK")
    # else:
    #     self.informationLabel.setText("Escape")
    # #
    #
    type = ''

    if buttons == 'YN':
        type = QtWidgets.QMessageBox.No | QtWidgets.QMessageBox.Yes
    elif buttons == "YNC":
        type  = QtWidgets.QMessageBox.Cancel | QtWidgets.QMessageBox.No | QtWidgets.QMessageBox.Yes # , QtWidgets.QMessageBox.Yes
    elif buttons == "OKC":
        type  = QtWidgets.QMessageBox.Cancel | QtWidgets.QMessageBox.Ok # | QtWidgets.QMessageBox.Yes # , QtWidgets.QMessageBox.Yes
    elif buttons == 'OK':
        type  = QtWidgets.QMessageBox.Ok #  | QtWidgets.QMessageBox.Yes # , QtWidgets.QMessageBox.Yes


    # QtWidgets.QMessageBox.setDefaultButton(self, QtWidgets.QMessageBox.Yes)

    reply = QtWidgets.QMessageBox.information(self, "QMessageBox.information()", statement, type)





    if reply == QtWidgets.QMessageBox.Yes:
        return 'Yes'
    elif reply == QtWidgets.QMessageBox.No:
        return 'No'
    elif reply == QtWidgets.QMessageBox.Ok:
        return 'OK'

    elif reply == QtWidgets.QMessageBox.Cancel:
        return 'Cancel'

def getItemDial(self, queryIs, items):
    # dialogue to generate combo box, OK, cancel

    item, ok = QtWidgets.QInputDialog.getItem(self, "QInputDialog.getItem()",
            queryIs, items, 0, False)
    if ok and item:
        return item

def setText(self, QueryIS, DefaultText):
    text, ok = QtWidgets.QInputDialog.getText(self, "Input",
            QueryIS, QtWidgets.QLineEdit.Normal, DefaultText) #QtCore.QDir.home().dirName())
            # QueryIS, DefaultText, QtCore.QDir.home().dirName())
    if ok and text != '':
        return text
    else:
        return 'Cancelled Action'





def setItem(self, items, title):
    # items = ("Spring", "Summer", "Fall", "Winter")

    item, ok = QtWidgets.QInputDialog.getItem(self, "QInputDialog.getItem()",
            title, items, 0, False)
    if ok and item:
        return item
    else:
        return'Cancel'



























