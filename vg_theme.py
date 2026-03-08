def db_tab_stylesheet():
    return """
QTableView#SeqTable {
    background: #f6f4ee;
    alternate-background-color: #ebe6d8;
    border: 1px solid #c7bea8;
    border-radius: 10px;
    gridline-color: #d9d0bd;
    selection-background-color: #0f6c73;
    selection-color: #fbfaf6;
    color: #2d261b;
    padding: 4px;
}
QTableView#SeqTable::item {
    padding: 6px 8px;
    border-bottom: 1px solid #e2dbc9;
}
QHeaderView::section {
    background: #d8ceb8;
    color: #2f2618;
    border: none;
    border-right: 1px solid #c1b69f;
    border-bottom: 1px solid #b6aa91;
    padding: 8px 6px;
    font-weight: 600;
}
QPushButton#ShowTable,
QPushButton#pushButtonRefresh,
QPushButton#pushButtonRefresh1,
QPushButton#pushButtonFirstPage,
QPushButton#pushButtonPreviousPage,
QPushButton#pushButtonNextPage,
QPushButton#pushButtonLastPage,
QPushButton#pushButtonJumpTo {
    background: #e6dcc8;
    color: #31281b;
    border: 1px solid #c3b79d;
    border-radius: 8px;
    padding: 6px 12px;
}
QPushButton#ShowTable:hover,
QPushButton#pushButtonRefresh:hover,
QPushButton#pushButtonRefresh1:hover,
QPushButton#pushButtonFirstPage:hover,
QPushButton#pushButtonPreviousPage:hover,
QPushButton#pushButtonNextPage:hover,
QPushButton#pushButtonLastPage:hover,
QPushButton#pushButtonJumpTo:hover {
    background: #f0e8d7;
}
QPushButton#ShowTable:pressed,
QPushButton#pushButtonRefresh:pressed,
QPushButton#pushButtonRefresh1:pressed,
QPushButton#pushButtonFirstPage:pressed,
QPushButton#pushButtonPreviousPage:pressed,
QPushButton#pushButtonNextPage:pressed,
QPushButton#pushButtonLastPage:pressed,
QPushButton#pushButtonJumpTo:pressed {
    background: #d7cab1;
}
QLineEdit#lineEditPageNumber,
QSpinBox#spinBoxPageSize {
    background: #fbfaf6;
    color: #2d261b;
    border: 1px solid #c7bea8;
    border-radius: 6px;
    padding: 4px 8px;
}
QCheckBox#checkBoxAll,
QCheckBox#checkBoxAll1 {
    color: #3d3322;
    spacing: 8px;
}
QCheckBox#checkBoxAll::indicator,
QCheckBox#checkBoxAll1::indicator {
    width: 16px;
    height: 16px;
}
QCheckBox#checkBoxAll::indicator:unchecked,
QCheckBox#checkBoxAll1::indicator:unchecked {
    background: #fbfaf6;
    border: 1px solid #bfae8a;
    border-radius: 4px;
}
QCheckBox#checkBoxAll::indicator:checked,
QCheckBox#checkBoxAll1::indicator:checked {
    background: #0f6c73;
    border: 1px solid #0f6c73;
    border-radius: 4px;
}
QLabel#dbTableStatusLabel {
    color: #5b4b33;
    background: #efe7d6;
    border: 1px solid #d1c5ac;
    border-radius: 8px;
    padding: 6px 10px;
}
"""
