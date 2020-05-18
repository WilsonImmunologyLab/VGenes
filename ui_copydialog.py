# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'copydialog.ui'
#
# Created by: PyQt5 UI code generator 5.13.0
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_CopyDialog(object):
    def setupUi(self, CopyDialog):
        CopyDialog.setObjectName("CopyDialog")
        CopyDialog.resize(572, 295)
        self.gridLayout_2 = QtWidgets.QGridLayout(CopyDialog)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.labelTip = QtWidgets.QLabel(CopyDialog)
        self.labelTip.setObjectName("labelTip")
        self.gridLayout_2.addWidget(self.labelTip, 1, 0, 1, 3)
        self.DisplayTip = QtWidgets.QPushButton(CopyDialog)
        self.DisplayTip.setObjectName("DisplayTip")
        self.gridLayout_2.addWidget(self.DisplayTip, 1, 3, 1, 2)
        self.gridLayoutFig = QtWidgets.QGridLayout()
        self.gridLayoutFig.setObjectName("gridLayoutFig")
        self.gridLayout_2.addLayout(self.gridLayoutFig, 2, 0, 2, 5)
        self.comboBoxFrom = QtWidgets.QComboBox(CopyDialog)
        self.comboBoxFrom.setObjectName("comboBoxFrom")
        self.gridLayout_2.addWidget(self.comboBoxFrom, 5, 0, 1, 5)
        self.label_2 = QtWidgets.QLabel(CopyDialog)
        self.label_2.setObjectName("label_2")
        self.gridLayout_2.addWidget(self.label_2, 6, 0, 1, 1)
        self.comboBoxTo = QtWidgets.QComboBox(CopyDialog)
        self.comboBoxTo.setObjectName("comboBoxTo")
        self.gridLayout_2.addWidget(self.comboBoxTo, 7, 0, 1, 5)
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout_2.addItem(spacerItem, 8, 1, 1, 1)
        spacerItem1 = QtWidgets.QSpacerItem(265, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout_2.addItem(spacerItem1, 9, 0, 1, 2)
        self.Cancel = QtWidgets.QPushButton(CopyDialog)
        self.Cancel.setObjectName("Cancel")
        self.gridLayout_2.addWidget(self.Cancel, 9, 2, 1, 1)
        self.OK = QtWidgets.QPushButton(CopyDialog)
        self.OK.setDefault(True)
        self.OK.setObjectName("OK")
        self.gridLayout_2.addWidget(self.OK, 9, 3, 1, 2)
        self.frame = QtWidgets.QFrame(CopyDialog)
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
        self.gridLayout_2.addWidget(self.frame, 0, 0, 1, 5)
        self.label = QtWidgets.QLabel(CopyDialog)
        self.label.setObjectName("label")
        self.gridLayout_2.addWidget(self.label, 4, 0, 1, 3)
        self.radioButton = QtWidgets.QRadioButton(CopyDialog)
        self.radioButton.setObjectName("radioButton")
        self.gridLayout_2.addWidget(self.radioButton, 4, 3, 1, 2)

        self.retranslateUi(CopyDialog)
        QtCore.QMetaObject.connectSlotsByName(CopyDialog)

    def retranslateUi(self, CopyDialog):
        _translate = QtCore.QCoreApplication.translate
        CopyDialog.setWindowTitle(_translate("CopyDialog", "Dialog"))
        self.labelTip.setText(_translate("CopyDialog", "Current field has more than 30 distint values! Hide figure by default!"))
        self.DisplayTip.setText(_translate("CopyDialog", "Display anyway"))
        self.label_2.setText(_translate("CopyDialog", "Copy value to:"))
        self.Cancel.setText(_translate("CopyDialog", "Cancel"))
        self.OK.setText(_translate("CopyDialog", "OK"))
        self.label_3.setText(_translate("CopyDialog", "Copy value from one column to another"))
        self.label.setText(_translate("CopyDialog", "Copy value from:"))
        self.radioButton.setText(_translate("CopyDialog", "Numerical"))
