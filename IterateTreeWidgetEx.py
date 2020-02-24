__author__ = 'PCW-MacBookProRet'
# not needed for python 3
import sip
sip.setapi('QString', 2)

from collections import defaultdict
from PyQt5 import QtGui, QtCore, QtWidgets

class Window(QtWidgets.QWidget):
    def __init__(self):
        QtWidgets.QWidget.__init__(self)
        self.tree = QtWidgets.QTreeWidget(self)
        self.tree.header().hide()
        for index in range(5):
            parent = QtWidgets.QTreeWidgetItem(self.tree, ['NUS2K%s' % index])
            if index % 3:
                parent.setCheckState(0, QtCore.Qt.PartiallyChecked)
            else:
                parent.setCheckState(0, QtCore.Qt.Unchecked)
            features = 'Loader Reports Logging'.split()
            for count, item in enumerate(features):
                child = QtWidgets.QTreeWidgetItem(parent, [item])
                if index % 3 and count % 3:
                    child.setCheckState(0, QtCore.Qt.Checked)
                else:
                    child.setCheckState(0, QtCore.Qt.Unchecked)
            parent.setExpanded(True)
        self.button = QtWidgets.QPushButton('Export', self)
        self.button.clicked.connect(self.handleExport)
        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(self.tree)
        layout.addWidget(self.button)

    def handleExport(self):
        mapping = self.exportTree()
        for machine, features in mapping.items():
            print('%s:' % machine)
            for feature in features:
                print('  %s' % feature)

    def exportTree(self):
        mapping = defaultdict(list)
        root = self.tree.invisibleRootItem()
        for index in range(root.childCount()):
            parent = root.child(index)
            if parent.checkState(0) == QtCore.Qt.PartiallyChecked:
                features = mapping[parent.text(0)]
                for row in range(parent.childCount()):
                    child = parent.child(row)
                    if child.checkState(0) == QtCore.Qt.Checked:
                        features.append(child.text(0))
        return mapping

if __name__ == '__main__':

    import sys
    app = QtWidgets.QApplication(sys.argv)
    window = Window()
    window.setGeometry(800, 300, 300, 300)
    window.show()
    sys.exit(app.exec_())