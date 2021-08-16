# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'rename_dialog.ui'
#
# Created by: PyQt5 UI code generator 5.13.0
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_RenameDialog(object):
    def setupUi(self, RenameDialog):
        RenameDialog.setObjectName("RenameDialog")
        RenameDialog.resize(539, 524)
        self.gridLayout_2 = QtWidgets.QGridLayout(RenameDialog)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.gridLayoutRename = QtWidgets.QGridLayout()
        self.gridLayoutRename.setObjectName("gridLayoutRename")
        self.gridLayout_2.addLayout(self.gridLayoutRename, 1, 0, 1, 3)
        self.pushButtonClearList = QtWidgets.QPushButton(RenameDialog)
        self.pushButtonClearList.setObjectName("pushButtonClearList")
        self.gridLayout_2.addWidget(self.pushButtonClearList, 2, 2, 1, 1)
        self.pushButtonEdit = QtWidgets.QPushButton(RenameDialog)
        self.pushButtonEdit.setMaximumSize(QtCore.QSize(150, 16777215))
        self.pushButtonEdit.setObjectName("pushButtonEdit")
        self.gridLayout_2.addWidget(self.pushButtonEdit, 4, 2, 1, 1)
        self.radioButtonHidePath = QtWidgets.QRadioButton(RenameDialog)
        self.radioButtonHidePath.setChecked(True)
        self.radioButtonHidePath.setObjectName("radioButtonHidePath")
        self.gridLayout_2.addWidget(self.radioButtonHidePath, 2, 0, 1, 1)
        self.frame = QtWidgets.QFrame(RenameDialog)
        self.frame.setMaximumSize(QtCore.QSize(16777215, 50))
        self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame.setObjectName("frame")
        self.gridLayout = QtWidgets.QGridLayout(self.frame)
        self.gridLayout.setObjectName("gridLayout")
        self.label_4 = QtWidgets.QLabel(self.frame)
        self.label_4.setObjectName("label_4")
        self.gridLayout.addWidget(self.label_4, 1, 2, 1, 2)
        self.spinBox = QtWidgets.QSpinBox(self.frame)
        self.spinBox.setObjectName("spinBox")
        self.gridLayout.addWidget(self.spinBox, 1, 1, 1, 1)
        self.label_3 = QtWidgets.QLabel(self.frame)
        self.label_3.setObjectName("label_3")
        self.gridLayout.addWidget(self.label_3, 1, 0, 1, 1)
        self.radioButtonBack = QtWidgets.QRadioButton(self.frame)
        self.radioButtonBack.setObjectName("radioButtonBack")
        self.gridLayout.addWidget(self.radioButtonBack, 1, 5, 1, 1)
        self.radioButtonFront = QtWidgets.QRadioButton(self.frame)
        self.radioButtonFront.setChecked(True)
        self.radioButtonFront.setObjectName("radioButtonFront")
        self.gridLayout.addWidget(self.radioButtonFront, 1, 4, 1, 1)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem, 1, 6, 1, 1)
        self.gridLayout_2.addWidget(self.frame, 3, 0, 1, 3)
        self.pushButtonUndo = QtWidgets.QPushButton(RenameDialog)
        self.pushButtonUndo.setEnabled(False)
        self.pushButtonUndo.setObjectName("pushButtonUndo")
        self.gridLayout_2.addWidget(self.pushButtonUndo, 4, 0, 1, 1)
        self.label = QtWidgets.QLabel(RenameDialog)
        self.label.setMaximumSize(QtCore.QSize(16777215, 20))
        font = QtGui.QFont()
        font.setPointSize(15)
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.gridLayout_2.addWidget(self.label, 0, 0, 1, 3)
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout_2.addItem(spacerItem1, 4, 1, 1, 1)

        self.retranslateUi(RenameDialog)
        QtCore.QMetaObject.connectSlotsByName(RenameDialog)

    def retranslateUi(self, RenameDialog):
        _translate = QtCore.QCoreApplication.translate
        RenameDialog.setWindowTitle(_translate("RenameDialog", "Edit file names in batches"))
        self.pushButtonClearList.setText(_translate("RenameDialog", "Clear list"))
        self.pushButtonEdit.setText(_translate("RenameDialog", "Edit file name"))
        self.radioButtonHidePath.setText(_translate("RenameDialog", "hide file path"))
        self.label_4.setText(_translate("RenameDialog", "Characters from"))
        self.label_3.setText(_translate("RenameDialog", "Remove"))
        self.radioButtonBack.setText(_translate("RenameDialog", "back"))
        self.radioButtonFront.setText(_translate("RenameDialog", "front"))
        self.pushButtonUndo.setText(_translate("RenameDialog", "Undo last edits"))
        self.label.setText(_translate("RenameDialog", "Drag and drop files below:"))
