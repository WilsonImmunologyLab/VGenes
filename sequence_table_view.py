from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QColor, QPainter, QStandardItem, QStandardItemModel
from PyQt5.QtWidgets import QTableView

ORIGINAL_VALUE_ROLE = Qt.UserRole + 1


class SequenceTableItem(QStandardItem):
    def __init__(self, text='', original_value=None):
        super().__init__(text)
        if original_value is None:
            original_value = text
        self.set_original_value(original_value)

    def original_value(self):
        return self.data(ORIGINAL_VALUE_ROLE)

    def set_original_value(self, value):
        self.setData(value, ORIGINAL_VALUE_ROLE)


class SequenceTableModel(QStandardItemModel):
    pass


class SequenceTableView(QTableView):
    itemChanged = pyqtSignal(object)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._model = SequenceTableModel(self)
        self._model.itemChanged.connect(self.itemChanged.emit)
        self.setModel(self._model)
        self._placeholder_text = ''

    def model(self):
        return self._model

    def setRowCount(self, rows):
        self._model.setRowCount(rows)

    def setColumnCount(self, columns):
        self._model.setColumnCount(columns)

    def rowCount(self):
        return self._model.rowCount()

    def columnCount(self):
        return self._model.columnCount()

    def setHorizontalHeaderLabels(self, labels):
        self._model.setHorizontalHeaderLabels(labels)

    def setItem(self, row, column, item):
        self._model.setItem(row, column, self._coerce_item(item))

    def item(self, row, column):
        return self._model.item(row, column)

    def selectedItems(self):
        return [self._model.itemFromIndex(index) for index in self.selectedIndexes()]

    def indexFromItem(self, item):
        return item.index()

    def setCurrentCell(self, row, column):
        index = self._model.index(row, column)
        if index.isValid():
            self.setCurrentIndex(index)

    def setPlaceholderText(self, text):
        self._placeholder_text = text
        self.viewport().update()

    def paintEvent(self, event):
        super().paintEvent(event)
        if self._model.rowCount() != 0 or not self._placeholder_text:
            return

        painter = QPainter(self.viewport())
        painter.setRenderHint(QPainter.TextAntialiasing)
        painter.setPen(QColor('#7a6c52'))
        rect = self.viewport().rect().adjusted(24, 24, -24, -24)
        painter.drawText(rect, Qt.AlignCenter | Qt.TextWordWrap, self._placeholder_text)
        painter.end()

    def _coerce_item(self, item):
        if isinstance(item, QStandardItem):
            return item

        text = ''
        try:
            text = item.text()
        except:
            text = str(item)

        standard_item = SequenceTableItem(text)
        try:
            standard_item.setFlags(item.flags())
        except:
            pass
        try:
            standard_item.setCheckState(item.checkState())
        except:
            pass
        try:
            standard_item.set_original_value(item.original_value())
        except:
            standard_item.set_original_value(text)
        return standard_item
