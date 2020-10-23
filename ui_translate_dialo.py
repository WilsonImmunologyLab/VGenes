# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'translate_dialog.ui'
#
# Created by: PyQt5 UI code generator 5.13.0
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Translate_Dialog(object):
    def setupUi(self, Translate_Dialog):
        Translate_Dialog.setObjectName("Translate_Dialog")
        Translate_Dialog.resize(1458, 918)
        self.gridLayout = QtWidgets.QGridLayout(Translate_Dialog)
        self.gridLayout.setObjectName("gridLayout")
        self.label = QtWidgets.QLabel(Translate_Dialog)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem, 0, 1, 1, 1)
        self.pushButton = QtWidgets.QPushButton(Translate_Dialog)
        self.pushButton.setObjectName("pushButton")
        self.gridLayout.addWidget(self.pushButton, 0, 2, 1, 1)
        self.textEditNT = QtWidgets.QTextEdit(Translate_Dialog)
        self.textEditNT.setObjectName("textEditNT")
        self.gridLayout.addWidget(self.textEditNT, 1, 0, 1, 3)
        self.label_2 = QtWidgets.QLabel(Translate_Dialog)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 2, 0, 1, 1)
        self.textEdit1 = QtWidgets.QTextEdit(Translate_Dialog)
        self.textEdit1.setObjectName("textEdit1")
        self.gridLayout.addWidget(self.textEdit1, 3, 0, 1, 3)
        self.label_3 = QtWidgets.QLabel(Translate_Dialog)
        self.label_3.setObjectName("label_3")
        self.gridLayout.addWidget(self.label_3, 4, 0, 1, 1)
        self.textEdit2 = QtWidgets.QTextEdit(Translate_Dialog)
        self.textEdit2.setObjectName("textEdit2")
        self.gridLayout.addWidget(self.textEdit2, 5, 0, 1, 3)
        self.label_4 = QtWidgets.QLabel(Translate_Dialog)
        self.label_4.setObjectName("label_4")
        self.gridLayout.addWidget(self.label_4, 6, 0, 1, 1)
        self.textEdit3 = QtWidgets.QTextEdit(Translate_Dialog)
        self.textEdit3.setObjectName("textEdit3")
        self.gridLayout.addWidget(self.textEdit3, 7, 0, 1, 3)
        self.label_5 = QtWidgets.QLabel(Translate_Dialog)
        self.label_5.setObjectName("label_5")
        self.gridLayout.addWidget(self.label_5, 8, 0, 1, 1)
        self.textEdit4 = QtWidgets.QTextEdit(Translate_Dialog)
        self.textEdit4.setObjectName("textEdit4")
        self.gridLayout.addWidget(self.textEdit4, 9, 0, 1, 3)
        self.label_6 = QtWidgets.QLabel(Translate_Dialog)
        self.label_6.setObjectName("label_6")
        self.gridLayout.addWidget(self.label_6, 10, 0, 1, 1)
        self.textEdit5 = QtWidgets.QTextEdit(Translate_Dialog)
        self.textEdit5.setObjectName("textEdit5")
        self.gridLayout.addWidget(self.textEdit5, 11, 0, 1, 3)
        self.label_7 = QtWidgets.QLabel(Translate_Dialog)
        self.label_7.setObjectName("label_7")
        self.gridLayout.addWidget(self.label_7, 12, 0, 1, 1)
        self.textEdit6 = QtWidgets.QTextEdit(Translate_Dialog)
        self.textEdit6.setObjectName("textEdit6")
        self.gridLayout.addWidget(self.textEdit6, 13, 0, 1, 3)

        self.retranslateUi(Translate_Dialog)
        QtCore.QMetaObject.connectSlotsByName(Translate_Dialog)

    def retranslateUi(self, Translate_Dialog):
        _translate = QtCore.QCoreApplication.translate
        Translate_Dialog.setWindowTitle(_translate("Translate_Dialog", "Dialog"))
        self.label.setText(_translate("Translate_Dialog", "DNA or RNA sequence"))
        self.pushButton.setText(_translate("Translate_Dialog", "Translate"))
        self.label_2.setText(_translate("Translate_Dialog", "5\'->3\' Frame 1"))
        self.label_3.setText(_translate("Translate_Dialog", "5\'->3\' Frame 2"))
        self.label_4.setText(_translate("Translate_Dialog", "5\'->3\' Frame 3"))
        self.label_5.setText(_translate("Translate_Dialog", "3\'->5\' Frame 1"))
        self.label_6.setText(_translate("Translate_Dialog", "3\'->5\' Frame 2"))
        self.label_7.setText(_translate("Translate_Dialog", "3\'->5\' Frame 3"))
