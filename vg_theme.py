from dataclasses import dataclass

from PyQt5.QtGui import QColor, QPalette


@dataclass(frozen=True)
class ThemeSpec:
    name: str
    window: str
    window_text: str
    central_start: str
    central_end: str
    pane: str
    pane_border: str
    tab_idle_bg: str
    tab_idle_text: str
    tab_active_bg: str
    tab_active_text: str
    panel_bg: str
    panel_border: str
    group_title: str
    tree_bg: str
    tree_alt_bg: str
    tree_hover_bg: str
    tree_border: str
    field_bg: str
    field_border: str
    field_text: str
    button_bg: str
    button_text: str
    button_border: str
    button_hover_bg: str
    button_pressed_bg: str
    button_flat_hover: str
    item_bg: str
    item_border: str
    status_bg: str
    menu_bg: str
    menu_hover_bg: str
    toolbar_bg: str
    toolbar_border: str
    toolbar_hover_bg: str
    toolbar_checked_bg: str
    accent: str
    accent_text: str
    table_bg: str
    table_alt_bg: str
    table_border: str
    table_grid: str
    table_text: str
    table_row_border: str
    header_bg: str
    header_text: str
    header_border: str
    check_border: str
    db_status_bg: str
    db_status_text: str
    danger: str


LIGHT_THEME = ThemeSpec(
    name="light",
    window="#efe8da",
    window_text="#2d261b",
    central_start="#f6f1e7",
    central_end="#ece3d2",
    pane="#f8f4eb",
    pane_border="#d0c3aa",
    tab_idle_bg="#ded3bf",
    tab_idle_text="#3a2f22",
    tab_active_bg="#f8f4eb",
    tab_active_text="#1f4f58",
    panel_bg="rgba(255, 252, 246, 0.88)",
    panel_border="#d7cbb4",
    group_title="#58472f",
    tree_bg="#fbf8f2",
    tree_alt_bg="#f0e9dd",
    tree_hover_bg="#e5dcc8",
    tree_border="#cdbfa6",
    field_bg="#fffdf9",
    field_border="#ccbea4",
    field_text="#2c2418",
    button_bg="#eadfc9",
    button_text="#302719",
    button_border="#cabda2",
    button_hover_bg="#f4ecd9",
    button_pressed_bg="#d9c9ac",
    button_flat_hover="rgba(225, 214, 190, 0.7)",
    item_bg="#fbf8f2",
    item_border="#d4c7af",
    status_bg="#e6dcc9",
    menu_bg="#f8f4eb",
    menu_hover_bg="#0f6c73",
    toolbar_bg="#d9cfbc",
    toolbar_border="#b9ad94",
    toolbar_hover_bg="#efe6d5",
    toolbar_checked_bg="#cfc3aa",
    accent="#0f6c73",
    accent_text="#fbfaf6",
    table_bg="#f6f4ee",
    table_alt_bg="#ebe6d8",
    table_border="#c7bea8",
    table_grid="#d9d0bd",
    table_text="#2d261b",
    table_row_border="#e2dbc9",
    header_bg="#d8ceb8",
    header_text="#2f2618",
    header_border="#c1b69f",
    check_border="#bfae8a",
    db_status_bg="#efe7d6",
    db_status_text="#5b4b33",
    danger="#8b2d1f",
)


DARK_THEME = ThemeSpec(
    name="dark",
    window="#171b20",
    window_text="#ecf1f6",
    central_start="#1c232a",
    central_end="#12171d",
    pane="#202932",
    pane_border="#33404d",
    tab_idle_bg="#283440",
    tab_idle_text="#c7d2de",
    tab_active_bg="#202932",
    tab_active_text="#89d4d9",
    panel_bg="rgba(31, 40, 50, 0.96)",
    panel_border="#384655",
    group_title="#d6dde6",
    tree_bg="#1d252d",
    tree_alt_bg="#242f39",
    tree_hover_bg="#30404e",
    tree_border="#394756",
    field_bg="#11181f",
    field_border="#415061",
    field_text="#edf2f7",
    button_bg="#2c3945",
    button_text="#edf2f7",
    button_border="#42515e",
    button_hover_bg="#374755",
    button_pressed_bg="#23303c",
    button_flat_hover="rgba(87, 110, 130, 0.45)",
    item_bg="#1d252d",
    item_border="#344252",
    status_bg="#202831",
    menu_bg="#222b34",
    menu_hover_bg="#0f6c73",
    toolbar_bg="#2a323a",
    toolbar_border="#414c57",
    toolbar_hover_bg="#36414c",
    toolbar_checked_bg="#233640",
    accent="#0f6c73",
    accent_text="#f7fbfc",
    table_bg="#1b2229",
    table_alt_bg="#222c35",
    table_border="#374452",
    table_grid="#303c49",
    table_text="#edf2f7",
    table_row_border="#2f3a46",
    header_bg="#283441",
    header_text="#edf2f7",
    header_border="#3a4858",
    check_border="#667789",
    db_status_bg="#232d37",
    db_status_text="#d8e1ea",
    danger="#ff7a6b",
)


OCEAN_THEME = ThemeSpec(
    name="ocean",
    window="#dbe8ee",
    window_text="#132a35",
    central_start="#edf5f8",
    central_end="#d6e5ec",
    pane="#f4f9fb",
    pane_border="#a8c0cc",
    tab_idle_bg="#bfd4de",
    tab_idle_text="#1a3948",
    tab_active_bg="#f4f9fb",
    tab_active_text="#0e6773",
    panel_bg="rgba(247, 251, 252, 0.93)",
    panel_border="#b7ccd6",
    group_title="#345665",
    tree_bg="#f7fbfc",
    tree_alt_bg="#e7f0f4",
    tree_hover_bg="#d7e7ee",
    tree_border="#adc4cf",
    field_bg="#ffffff",
    field_border="#b4cad4",
    field_text="#17303b",
    button_bg="#d7e7ee",
    button_text="#18323e",
    button_border="#afc6d0",
    button_hover_bg="#e6f1f5",
    button_pressed_bg="#c6dbe4",
    button_flat_hover="rgba(173, 198, 208, 0.45)",
    item_bg="#f7fbfc",
    item_border="#bfd1d9",
    status_bg="#dbe8ee",
    menu_bg="#f4f9fb",
    menu_hover_bg="#0e6773",
    toolbar_bg="#c9dbe3",
    toolbar_border="#a7bcc6",
    toolbar_hover_bg="#e7f1f5",
    toolbar_checked_bg="#b9d0d9",
    accent="#0e6773",
    accent_text="#f7fbfc",
    table_bg="#f4f9fb",
    table_alt_bg="#e6f0f4",
    table_border="#b3c9d2",
    table_grid="#ccdbe1",
    table_text="#17303b",
    table_row_border="#dbe7ec",
    header_bg="#d3e2e9",
    header_text="#18303a",
    header_border="#b5c9d1",
    check_border="#98b5c0",
    db_status_bg="#e9f2f6",
    db_status_text="#30505d",
    danger="#a33b2e",
)


FOREST_THEME = ThemeSpec(
    name="forest",
    window="#e1e8df",
    window_text="#213126",
    central_start="#f0f4ed",
    central_end="#d9e3d6",
    pane="#f7faf5",
    pane_border="#b7c6b2",
    tab_idle_bg="#d1dccb",
    tab_idle_text="#294032",
    tab_active_bg="#f7faf5",
    tab_active_text="#2f6c54",
    panel_bg="rgba(250, 252, 248, 0.92)",
    panel_border="#c2d0bd",
    group_title="#486251",
    tree_bg="#fbfcf9",
    tree_alt_bg="#edf2e8",
    tree_hover_bg="#dce8d8",
    tree_border="#becdb8",
    field_bg="#ffffff",
    field_border="#c1d0bb",
    field_text="#25372b",
    button_bg="#dce8d7",
    button_text="#24362b",
    button_border="#b8c8b3",
    button_hover_bg="#e8f0e4",
    button_pressed_bg="#cdddc8",
    button_flat_hover="rgba(184, 200, 179, 0.5)",
    item_bg="#fbfcf9",
    item_border="#c8d4c4",
    status_bg="#dfe8db",
    menu_bg="#f7faf5",
    menu_hover_bg="#2f6c54",
    toolbar_bg="#cedbc8",
    toolbar_border="#aebea8",
    toolbar_hover_bg="#e8f0e4",
    toolbar_checked_bg="#c0cfb9",
    accent="#2f6c54",
    accent_text="#f6faf8",
    table_bg="#f7faf5",
    table_alt_bg="#e9f0e5",
    table_border="#bdcab8",
    table_grid="#d6ded2",
    table_text="#25372b",
    table_row_border="#e0e8dd",
    header_bg="#d7e1d2",
    header_text="#26382d",
    header_border="#bcc9b7",
    check_border="#9db099",
    db_status_bg="#edf3ea",
    db_status_text="#4b5f51",
    danger="#9b3d33",
)


THEMES = {
    LIGHT_THEME.name: LIGHT_THEME,
    DARK_THEME.name: DARK_THEME,
    OCEAN_THEME.name: OCEAN_THEME,
    FOREST_THEME.name: FOREST_THEME,
}

THEME_LABELS = {
    "light": "Light Studio",
    "dark": "Dark Slate",
    "ocean": "Ocean Mist",
    "forest": "Forest Paper",
}


def theme_spec(theme_name="light"):
    return THEMES.get(theme_name, LIGHT_THEME)


def theme_choices():
    return [(name, THEME_LABELS.get(name, name.title())) for name in THEMES]


def build_palette(theme_name="light"):
    theme = theme_spec(theme_name)
    palette = QPalette()
    palette.setColor(QPalette.Window, QColor(theme.window))
    palette.setColor(QPalette.WindowText, QColor(theme.window_text))
    palette.setColor(QPalette.Base, QColor(theme.field_bg))
    palette.setColor(QPalette.AlternateBase, QColor(theme.tree_alt_bg))
    palette.setColor(QPalette.ToolTipBase, QColor(theme.field_bg))
    palette.setColor(QPalette.ToolTipText, QColor(theme.field_text))
    palette.setColor(QPalette.Text, QColor(theme.field_text))
    palette.setColor(QPalette.Button, QColor(theme.button_bg))
    palette.setColor(QPalette.ButtonText, QColor(theme.button_text))
    palette.setColor(QPalette.BrightText, QColor(theme.danger))
    palette.setColor(QPalette.Highlight, QColor(theme.accent))
    palette.setColor(QPalette.HighlightedText, QColor(theme.accent_text))
    palette.setColor(QPalette.Link, QColor(theme.accent))
    return palette


def main_window_stylesheet(theme_name="light"):
    theme = theme_spec(theme_name)
    return f"""
QMainWindow {{
    background: {theme.window};
    color: {theme.window_text};
}}
QWidget#centralwidget {{
    background: qlineargradient(
        x1: 0, y1: 0, x2: 1, y2: 1,
        stop: 0 {theme.central_start},
        stop: 1 {theme.central_end}
    );
}}
QWidget {{
    color: {theme.window_text};
}}
QTabWidget::pane {{
    border: 1px solid {theme.pane_border};
    background: {theme.pane};
    border-radius: 10px;
    top: -2px;
}}
QTabWidget {{
    background: transparent;
}}
QTabWidget > QWidget {{
    background: {theme.pane};
    color: {theme.window_text};
}}
QTabBar {{
    background: transparent;
    font-size: 15px;
    font-weight: 600;
}}
QTabBar::tab-bar {{
    alignment: left;
}}
QTabBar::tab {{
    background: {theme.tab_idle_bg};
    color: {theme.tab_idle_text};
    border: 1px solid {theme.pane_border};
    border-bottom: none;
    border-top-left-radius: 8px;
    border-top-right-radius: 8px;
    padding: 8px 14px;
    margin-right: 4px;
    min-width: 110px;
}}
QTabBar::tab:selected {{
    background: {theme.tab_active_bg};
    color: {theme.tab_active_text};
    margin-bottom: -2px;
}}
QTabBar::tab:!selected {{
    margin-top: 3px;
}}
QFrame#frame_3,
QFrame#frame_4,
QGroupBox#groupBox_2,
QLabel#label_Name {{
    background: {theme.panel_bg};
    border: 1px solid {theme.panel_border};
    border-radius: 10px;
}}
QGroupBox#groupBox_2 {{
    margin-top: 12px;
    padding-top: 12px;
    font-weight: 600;
}}
QGroupBox,
QCheckBox,
QRadioButton {{
    color: {theme.window_text};
}}
QGroupBox#groupBox_2::title {{
    subcontrol-origin: margin;
    left: 10px;
    padding: 0 6px;
    color: {theme.group_title};
}}
QTreeWidget#treeWidget {{
    background: {theme.tree_bg};
    alternate-background-color: {theme.tree_alt_bg};
    border: 1px solid {theme.tree_border};
    border-radius: 10px;
    padding: 4px;
}}
QTreeWidget#treeWidget::item:selected {{
    background: {theme.accent};
    color: {theme.accent_text};
}}
QTreeWidget#treeWidget::item:hover {{
    background: {theme.tree_hover_bg};
    color: {theme.window_text};
}}
QAbstractScrollArea,
QListWidget,
QTableWidget,
QTextBrowser,
QScrollArea {{
    background: {theme.item_bg};
    color: {theme.window_text};
    border: 1px solid {theme.item_border};
    border-radius: 8px;
}}
QAbstractItemView {{
    selection-background-color: {theme.accent};
    selection-color: {theme.accent_text};
}}
QPushButton {{
    background: {theme.button_bg};
    color: {theme.button_text};
    border: 1px solid {theme.button_border};
    border-radius: 2px;
    padding: 0px 12px;
}}
QPushButton:hover {{
    background: {theme.button_hover_bg};
}}
QPushButton:pressed {{
    background: {theme.button_pressed_bg};
}}
QPushButton[flat="true"] {{
    background: transparent;
    border-color: transparent;
    padding: 6px 8px;
}}
QPushButton[flat="true"]:hover {{
    background: {theme.button_flat_hover};
    border-color: {theme.button_border};
}}
QLineEdit,
QComboBox,
QSpinBox,
QPlainTextEdit,
QTextEdit {{
    background: {theme.field_bg};
    color: {theme.field_text};
    border: 1px solid {theme.field_border};
    border-radius: 8px;
    padding: 5px 8px;
}}
QComboBox {{
    padding-right: 28px;
}}
QComboBox::drop-down {{
    subcontrol-origin: padding;
    subcontrol-position: top right;
    width: 26px;
    border-left: 1px solid {theme.field_border};
    background: transparent;
    border-top-right-radius: 2px;
    border-bottom-right-radius: 2px;
}}
QComboBox::down-arrow {{
    width: 0px;
    height: 0px;
    border-left: 5px solid transparent;
    border-right: 5px solid transparent;
    border-top: 7px solid {theme.field_text};
    margin-right: 8px;
}}
QComboBox QAbstractItemView {{
    background: {theme.menu_bg};
    color: {theme.window_text};
    selection-background-color: {theme.accent};
    selection-color: {theme.accent_text};
    border: 1px solid {theme.pane_border};
}}
QLabel {{
    color: {theme.window_text};
}}
QMenuBar,
QStatusBar {{
    background: {theme.status_bg};
    color: {theme.window_text};
}}
QMenuBar::item:selected,
QMenu::item:selected {{
    background: {theme.menu_hover_bg};
    color: {theme.accent_text};
}}
QMenu {{
    background: {theme.menu_bg};
    color: {theme.window_text};
    border: 1px solid {theme.pane_border};
}}
QToolBar {{
    background: {theme.toolbar_bg};
    color: {theme.window_text};
    border: none;
    border-bottom: 1px solid {theme.toolbar_border};
    spacing: 4px;
    padding: 3px 6px;
}}
QToolBar::separator {{
    background: {theme.toolbar_border};
    width: 1px;
    margin: 4px 6px;
}}
QToolButton {{
    background: transparent;
    color: {theme.window_text};
    border: 1px solid transparent;
    border-radius: 2px;
    padding: 2px 4px;
}}
QToolButton:hover {{
    background: {theme.toolbar_hover_bg};
    border-color: {theme.toolbar_border};
}}
QToolButton:checked,
QToolButton:pressed {{
    background: {theme.toolbar_checked_bg};
    border-color: {theme.toolbar_border};
}}
QLabel#tableTabInfoLabel {{
    color: {theme.group_title};
    background: {theme.panel_bg};
    border: 1px solid {theme.panel_border};
    border-radius: 2px;
    padding: 8px 10px;
}}
"""


def db_tab_stylesheet(theme_name="light"):
    theme = theme_spec(theme_name)
    return f"""
QTableView#SeqTable {{
    background: {theme.table_bg};
    alternate-background-color: {theme.table_alt_bg};
    border: 1px solid {theme.table_border};
    border-radius: 2px;
    gridline-color: {theme.table_grid};
    selection-background-color: {theme.accent};
    selection-color: {theme.accent_text};
    color: {theme.table_text};
    padding: 4px;
}}
QTableView#SeqTable::item {{
    padding: 6px 8px;
    border-bottom: 1px solid {theme.table_row_border};
}}
QHeaderView::section {{
    background: {theme.header_bg};
    color: {theme.header_text};
    border: none;
    border-right: 1px solid {theme.header_border};
    border-bottom: 1px solid {theme.header_border};
    padding: 8px 6px;
    font-weight: 600;
}}
QPushButton#ShowTable,
QPushButton#pushButtonRefresh,
QPushButton#pushButtonRefresh1,
QPushButton#pushButtonFirstPage,
QPushButton#pushButtonPreviousPage,
QPushButton#pushButtonNextPage,
QPushButton#pushButtonLastPage,
QPushButton#pushButtonJumpTo {{
    background: {theme.button_bg};
    color: {theme.button_text};
    border: 1px solid {theme.button_border};
    border-radius: 8px;
    padding: 0px 12px;
}}
QPushButton#ShowTable:hover,
QPushButton#pushButtonRefresh:hover,
QPushButton#pushButtonRefresh1:hover,
QPushButton#pushButtonFirstPage:hover,
QPushButton#pushButtonPreviousPage:hover,
QPushButton#pushButtonNextPage:hover,
QPushButton#pushButtonLastPage:hover,
QPushButton#pushButtonJumpTo:hover {{
    background: {theme.button_hover_bg};
}}
QPushButton#ShowTable:pressed,
QPushButton#pushButtonRefresh:pressed,
QPushButton#pushButtonRefresh1:pressed,
QPushButton#pushButtonFirstPage:pressed,
QPushButton#pushButtonPreviousPage:pressed,
QPushButton#pushButtonNextPage:pressed,
QPushButton#pushButtonLastPage:pressed,
QPushButton#pushButtonJumpTo:pressed {{
    background: {theme.button_pressed_bg};
}}
QLineEdit#lineEditPageNumber,
QSpinBox#spinBoxPageSize {{
    background: {theme.field_bg};
    color: {theme.field_text};
    border: 1px solid {theme.field_border};
    border-radius: 2px;
    padding: 4px 8px;
}}
QCheckBox#checkBoxAll,
QCheckBox#checkBoxAll1 {{
    color: {theme.window_text};
    spacing: 8px;
}}
QCheckBox#checkBoxAll::indicator,
QCheckBox#checkBoxAll1::indicator {{
    width: 16px;
    height: 16px;
}}
QCheckBox#checkBoxAll::indicator:unchecked,
QCheckBox#checkBoxAll1::indicator:unchecked {{
    background: {theme.field_bg};
    border: 1px solid {theme.check_border};
    border-radius: 2px;
}}
QCheckBox#checkBoxAll::indicator:checked,
QCheckBox#checkBoxAll1::indicator:checked {{
    background: {theme.accent};
    border: 1px solid {theme.accent};
    border-radius: 2px;
}}
QLabel#dbTableStatusLabel {{
    color: {theme.db_status_text};
    background: {theme.db_status_bg};
    border: 1px solid {theme.table_border};
    border-radius: 2px;
    padding: 6px 10px;
}}
"""
