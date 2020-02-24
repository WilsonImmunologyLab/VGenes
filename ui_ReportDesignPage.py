# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ReportDesignPage.ui'
#
# Created: Thu Dec 31 13:58:29 2015
#      by: PyQt5 UI code generator 5.3.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_WizardPage(object):
    def setupUi(self, WizardPage):
        WizardPage.setObjectName("WizardPage")
        WizardPage.resize(510, 300)
        self.listWidget = QtWidgets.QListWidget(WizardPage)
        self.listWidget.setGeometry(QtCore.QRect(130, 50, 131, 192))
        self.listWidget.setTabKeyNavigation(True)
        self.listWidget.setObjectName("listWidget")
        self.listWidget_2 = QtWidgets.QListWidget(WizardPage)
        self.listWidget_2.setGeometry(QtCore.QRect(310, 50, 131, 192))
        self.listWidget_2.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.listWidget_2.setTabKeyNavigation(True)
        self.listWidget_2.setDragEnabled(True)
        self.listWidget_2.setDragDropMode(QtWidgets.QAbstractItemView.InternalMove)
        self.listWidget_2.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.listWidget_2.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectColumns)
        self.listWidget_2.setObjectName("listWidget_2")
        self.pushButton = QtWidgets.QPushButton(WizardPage)
        self.pushButton.setGeometry(QtCore.QRect(370, 250, 114, 32))
        self.pushButton.setObjectName("pushButton")
        self.toolButton_2 = QtWidgets.QToolButton(WizardPage)
        self.toolButton_2.setGeometry(QtCore.QRect(270, 100, 31, 41))
        self.toolButton_2.setObjectName("toolButton_2")
        self.toolButton_3 = QtWidgets.QToolButton(WizardPage)
        self.toolButton_3.setGeometry(QtCore.QRect(270, 50, 31, 41))
        self.toolButton_3.setObjectName("toolButton_3")
        self.toolButton_4 = QtWidgets.QToolButton(WizardPage)
        self.toolButton_4.setGeometry(QtCore.QRect(270, 150, 31, 41))
        self.toolButton_4.setObjectName("toolButton_4")
        self.toolButton_5 = QtWidgets.QToolButton(WizardPage)
        self.toolButton_5.setGeometry(QtCore.QRect(270, 200, 31, 41))
        self.toolButton_5.setObjectName("toolButton_5")
        self.label = QtWidgets.QLabel(WizardPage)
        self.label.setGeometry(QtCore.QRect(180, 30, 62, 16))
        self.label.setObjectName("label")
        self.label_2 = QtWidgets.QLabel(WizardPage)
        self.label_2.setGeometry(QtCore.QRect(360, 30, 62, 16))
        self.label_2.setObjectName("label_2")
        self.pushButton_2 = QtWidgets.QPushButton(WizardPage)
        self.pushButton_2.setGeometry(QtCore.QRect(250, 250, 114, 32))
        self.pushButton_2.setObjectName("pushButton_2")

        self.retranslateUi(WizardPage)
        QtCore.QMetaObject.connectSlotsByName(WizardPage)

    def retranslateUi(self, WizardPage):
        _translate = QtCore.QCoreApplication.translate
        WizardPage.setWindowTitle(_translate("WizardPage", "WizardPage"))
        self.pushButton.setText(_translate("WizardPage", "Continue"))
        self.toolButton_2.setText(_translate("WizardPage", ">>"))
        self.toolButton_3.setText(_translate("WizardPage", ">"))
        self.toolButton_4.setText(_translate("WizardPage", "<"))
        self.toolButton_5.setText(_translate("WizardPage", "<<"))
        self.label.setText(_translate("WizardPage", "Fields"))
        self.label_2.setText(_translate("WizardPage", "Report"))
        self.pushButton_2.setText(_translate("WizardPage", "Cancel"))

