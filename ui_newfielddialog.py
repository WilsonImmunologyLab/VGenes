# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'newfielddialog.ui'
#
# Created by: PyQt5 UI code generator 5.13.0
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_NewFieldDialog(object):
    def setupUi(self, NewFieldDialog):
        NewFieldDialog.setObjectName("NewFieldDialog")
        NewFieldDialog.resize(581, 287)
        self.gridLayout_2 = QtWidgets.QGridLayout(NewFieldDialog)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.comboBoxFrom = QtWidgets.QComboBox(NewFieldDialog)
        self.comboBoxFrom.setObjectName("comboBoxFrom")
        self.gridLayout_2.addWidget(self.comboBoxFrom, 4, 0, 1, 3)
        self.label_2 = QtWidgets.QLabel(NewFieldDialog)
        self.label_2.setObjectName("label_2")
        self.gridLayout_2.addWidget(self.label_2, 3, 0, 1, 1)
        self.labelTip = QtWidgets.QLabel(NewFieldDialog)
        self.labelTip.setObjectName("labelTip")
        self.gridLayout_2.addWidget(self.labelTip, 6, 0, 1, 2)
        self.DisplayTip = QtWidgets.QPushButton(NewFieldDialog)
        self.DisplayTip.setObjectName("DisplayTip")
        self.gridLayout_2.addWidget(self.DisplayTip, 6, 2, 1, 1)
        self.radioButton = QtWidgets.QRadioButton(NewFieldDialog)
        self.radioButton.setObjectName("radioButton")
        self.gridLayout_2.addWidget(self.radioButton, 3, 2, 1, 1)
        self.lineEdit = QtWidgets.QLineEdit(NewFieldDialog)
        self.lineEdit.setObjectName("lineEdit")
        self.gridLayout_2.addWidget(self.lineEdit, 2, 0, 1, 3)
        self.gridLayoutFig = QtWidgets.QGridLayout()
        self.gridLayoutFig.setObjectName("gridLayoutFig")
        self.gridLayout_2.addLayout(self.gridLayoutFig, 7, 0, 1, 3)
        self.label = QtWidgets.QLabel(NewFieldDialog)
        self.label.setObjectName("label")
        self.gridLayout_2.addWidget(self.label, 1, 0, 1, 1)
        spacerItem = QtWidgets.QSpacerItem(233, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout_2.addItem(spacerItem, 9, 0, 1, 1)
        self.OK = QtWidgets.QPushButton(NewFieldDialog)
        self.OK.setMinimumSize(QtCore.QSize(100, 0))
        self.OK.setMaximumSize(QtCore.QSize(100, 16777215))
        self.OK.setDefault(True)
        self.OK.setObjectName("OK")
        self.gridLayout_2.addWidget(self.OK, 9, 2, 1, 1)
        self.frame = QtWidgets.QFrame(NewFieldDialog)
        self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame.setObjectName("frame")
        self.gridLayout = QtWidgets.QGridLayout(self.frame)
        self.gridLayout.setObjectName("gridLayout")
        self.label_3 = QtWidgets.QLabel(self.frame)
        font = QtGui.QFont()
        font.setPointSize(15)
        self.label_3.setFont(font)
        self.label_3.setAlignment(QtCore.Qt.AlignCenter)
        self.label_3.setObjectName("label_3")
        self.gridLayout.addWidget(self.label_3, 0, 0, 1, 1)
        self.gridLayout_2.addWidget(self.frame, 0, 0, 1, 3)
        self.Cancel = QtWidgets.QPushButton(NewFieldDialog)
        self.Cancel.setMinimumSize(QtCore.QSize(100, 0))
        self.Cancel.setMaximumSize(QtCore.QSize(100, 16777215))
        self.Cancel.setObjectName("Cancel")
        self.gridLayout_2.addWidget(self.Cancel, 9, 1, 1, 1)
        spacerItem1 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout_2.addItem(spacerItem1, 8, 0, 1, 1)

        self.retranslateUi(NewFieldDialog)
        QtCore.QMetaObject.connectSlotsByName(NewFieldDialog)

    def retranslateUi(self, NewFieldDialog):
        _translate = QtCore.QCoreApplication.translate
        NewFieldDialog.setWindowTitle(_translate("NewFieldDialog", "Dialog"))
        self.label_2.setText(_translate("NewFieldDialog", "Copy value from (optional): "))
        self.labelTip.setText(_translate("NewFieldDialog", "Current field has more than 30 distint values! Hide figure by default!"))
        self.DisplayTip.setText(_translate("NewFieldDialog", "Display anyway"))
        self.radioButton.setText(_translate("NewFieldDialog", "Numerical"))
        self.label.setText(_translate("NewFieldDialog", "Field name:"))
        self.OK.setText(_translate("NewFieldDialog", "OK"))
        self.label_3.setText(_translate("NewFieldDialog", "Add new field to your VGene DB"))
        self.Cancel.setText(_translate("NewFieldDialog", "Cancel"))
