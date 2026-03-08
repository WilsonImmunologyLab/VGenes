from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QAbstractItemView

from sequence_table_view import SequenceTableItem


class SequenceTableAdapter:
    """Centralize main sequence table state and row operations."""

    def __init__(self, table_widget):
        self.table = table_widget
        self._page_size = 20
        self._edit_tag = False
        self._fields = []
        self._db_fields = []
        self._row_names = []
        self._row_index_by_name = {}

    def initialize_state(self, page_size=20):
        self._page_size = page_size
        self._edit_tag = False
        self._fields = []
        self._db_fields = []
        self._row_names = []
        self._row_index_by_name = {}

    def clear(self):
        self.table.setColumnCount(0)
        self.table.setRowCount(0)
        self._fields = []
        self._db_fields = []
        self._row_names = []
        self._row_index_by_name = {}

    def page_size(self):
        return self._page_size

    def set_page_size(self, page_size):
        self._page_size = page_size

    def edit_tag(self):
        return self._edit_tag

    def set_edit_tag(self, edit_tag):
        self._edit_tag = edit_tag

    def header_label(self, column):
        return self._fields[column]

    def db_field(self, column):
        return self._db_fields[column]

    def row_count(self):
        return self.table.rowCount()

    def column_count(self):
        return self.table.columnCount()

    def edit_triggers(self):
        return self.table.editTriggers()

    def set_edit_triggers(self, triggers):
        self.table.setEditTriggers(triggers)

    def row_name(self, row):
        try:
            return self._row_names[row]
        except:
            item = self.table.item(row, 1)
            if item is None:
                return ''
            return item.text()

    def update_row_name(self, row, name):
        old_name = ''
        try:
            old_name = self._row_names[row]
        except:
            pass
        if old_name in self._row_index_by_name:
            del self._row_index_by_name[old_name]
        if row < len(self._row_names):
            self._row_names[row] = name
        self._row_index_by_name[name] = row

    def check_item(self, row):
        return self.table.item(row, 0)

    def item(self, row, column):
        return self.table.item(row, column)

    def cell_text(self, row, column, default=''):
        item = self.item(row, column)
        if item is None:
            return default
        return item.text()

    def set_cell_text(self, row, column, text):
        item = self.item(row, column)
        if item is not None:
            item.setText(text)

    def is_checked(self, row):
        item = self.check_item(row)
        if item is None:
            return False
        return item.checkState() == Qt.Checked

    def set_checked(self, row, checked):
        item = self.check_item(row)
        if item is None:
            return
        target_state = Qt.Checked if checked else Qt.Unchecked
        if item.checkState() != target_state:
            item.setCheckState(target_state)

    def selected_rows(self):
        selection_model = self.table.selectionModel()
        if selection_model is None:
            return []

        if self.table.selectionBehavior() == QAbstractItemView.SelectRows:
            return sorted(index.row() for index in selection_model.selectedRows())

        rows = {index.row() for index in selection_model.selectedIndexes()}
        return sorted(rows)

    def selected_items(self):
        return self.table.selectedItems()

    def selected_names(self):
        return [self.row_name(row) for row in self.selected_rows()]

    def set_current_name(self, name, column=1):
        row = self._row_index_by_name.get(name)
        if row is None:
            return False
        self.table.setCurrentCell(row, column)
        return True

    def begin_reload(self):
        self.table.blockSignals(True)
        self.table.setUpdatesEnabled(False)
        self.table.setSortingEnabled(False)

    def end_reload(self):
        self.table.blockSignals(False)
        self.table.setUpdatesEnabled(True)
        self.table.setSortingEnabled(False)

    def configure(self, db_fields, display_fields, row_names):
        self.table.setRowCount(len(row_names))
        self.table.setColumnCount(len(db_fields))
        self.table.setHorizontalHeaderLabels(display_fields)
        self._fields = list(display_fields)
        self._db_fields = list(db_fields)
        self._row_names = list(row_names)
        self._row_index_by_name = {name: index for index, name in enumerate(self._row_names)}

    def apply_selection_behavior(self, select_rows):
        self.table.setSelectionMode(QAbstractItemView.ExtendedSelection)
        if select_rows:
            self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        else:
            self.table.setSelectionBehavior(QAbstractItemView.SelectItems)

    def populate_rows(self, data_rows, checked_records):
        checked_names = set(str(name) for name in checked_records)
        for row_index, row_data in enumerate(data_rows):
            check_item = SequenceTableItem()
            check_item.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled | Qt.ItemIsSelectable)
            if str(row_data[0]) in checked_names:
                check_item.setCheckState(Qt.Checked)
            else:
                check_item.setCheckState(Qt.Unchecked)
            self.table.setItem(row_index, 0, check_item)

            for col_index, value in enumerate(row_data):
                unit = SequenceTableItem(str(value), value)
                self.table.setItem(row_index, col_index + 1, unit)

    def apply_edit_mode(self):
        if self._edit_tag:
            self.table.setEditTriggers(QtWidgets.QAbstractItemView.DoubleClicked)
        else:
            self.table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)

    def set_sort_indicator(self, column_index, sort_order):
        self.table.horizontalHeader().setSortIndicatorShown(True)
        self.table.horizontalHeader().setSortIndicator(column_index, sort_order)

    def resize_primary_columns(self, check_width=10, name_width=None):
        self.table.horizontalHeader().resizeSection(0, check_width)
        if name_width is not None and self.table.columnCount() > 1:
            self.table.horizontalHeader().resizeSection(1, name_width)
