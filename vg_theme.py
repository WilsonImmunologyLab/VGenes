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


CERULEAN_THEME = ThemeSpec(
    name="cerulean",
    window="#eef5f9",
    window_text="#243746",
    central_start="#f7fbfd",
    central_end="#e3eff6",
    pane="#ffffff",
    pane_border="#b9d1df",
    tab_idle_bg="#d7e7f1",
    tab_idle_text="#25526f",
    tab_active_bg="#ffffff",
    tab_active_text="#2a7db8",
    panel_bg="rgba(255, 255, 255, 0.94)",
    panel_border="#c7dbe7",
    group_title="#2f668c",
    tree_bg="#ffffff",
    tree_alt_bg="#eef6fb",
    tree_hover_bg="#dcecf7",
    tree_border="#bfd5e3",
    field_bg="#ffffff",
    field_border="#b9d1df",
    field_text="#274055",
    button_bg="#2fa4e7",
    button_text="#ffffff",
    button_border="#2389c5",
    button_hover_bg="#48afea",
    button_pressed_bg="#1f7db6",
    button_flat_hover="rgba(47, 164, 231, 0.16)",
    item_bg="#ffffff",
    item_border="#c7dbe7",
    status_bg="#dcebf5",
    menu_bg="#ffffff",
    menu_hover_bg="#2fa4e7",
    toolbar_bg="#e6f1f8",
    toolbar_border="#b7cedc",
    toolbar_hover_bg="#f4f9fc",
    toolbar_checked_bg="#d4e7f3",
    accent="#2fa4e7",
    accent_text="#ffffff",
    table_bg="#ffffff",
    table_alt_bg="#eff7fb",
    table_border="#b9d1df",
    table_grid="#d9e8f1",
    table_text="#274055",
    table_row_border="#e3eef5",
    header_bg="#d8e8f2",
    header_text="#23425a",
    header_border="#bfd4e1",
    check_border="#9bbad0",
    db_status_bg="#eef6fb",
    db_status_text="#4c6880",
    danger="#c1453d",
)


BRITE_THEME = ThemeSpec(
    name="brite",
    window="#fff7eb",
    window_text="#231815",
    central_start="#fffdf8",
    central_end="#ffeccf",
    pane="#fffaf0",
    pane_border="#231815",
    tab_idle_bg="#ffd44d",
    tab_idle_text="#231815",
    tab_active_bg="#ffffff",
    tab_active_text="#111111",
    panel_bg="rgba(255, 250, 240, 0.96)",
    panel_border="#231815",
    group_title="#111111",
    tree_bg="#ffffff",
    tree_alt_bg="#fff3d8",
    tree_hover_bg="#ffe18a",
    tree_border="#231815",
    field_bg="#ffffff",
    field_border="#231815",
    field_text="#111111",
    button_bg="#ff6f61",
    button_text="#111111",
    button_border="#231815",
    button_hover_bg="#ff866f",
    button_pressed_bg="#ef5d50",
    button_flat_hover="rgba(255, 111, 97, 0.18)",
    item_bg="#ffffff",
    item_border="#231815",
    status_bg="#ffe8b6",
    menu_bg="#fffaf0",
    menu_hover_bg="#ff6f61",
    toolbar_bg="#ffe083",
    toolbar_border="#231815",
    toolbar_hover_bg="#fff0b8",
    toolbar_checked_bg="#ffd05d",
    accent="#ff6f61",
    accent_text="#111111",
    table_bg="#ffffff",
    table_alt_bg="#fff4dc",
    table_border="#231815",
    table_grid="#eacb83",
    table_text="#111111",
    table_row_border="#f0d79d",
    header_bg="#ffd44d",
    header_text="#111111",
    header_border="#231815",
    check_border="#231815",
    db_status_bg="#ffe7b4",
    db_status_text="#231815",
    danger="#c13c2e",
)


COSMO_THEME = ThemeSpec(
    name="cosmo",
    window="#f7fafc",
    window_text="#23313d",
    central_start="#ffffff",
    central_end="#ecf2f7",
    pane="#ffffff",
    pane_border="#d4dee7",
    tab_idle_bg="#dbe6ef",
    tab_idle_text="#31475b",
    tab_active_bg="#ffffff",
    tab_active_text="#2780e3",
    panel_bg="rgba(255, 255, 255, 0.95)",
    panel_border="#d4dee7",
    group_title="#355166",
    tree_bg="#ffffff",
    tree_alt_bg="#f1f6fa",
    tree_hover_bg="#e3edf6",
    tree_border="#d4dee7",
    field_bg="#ffffff",
    field_border="#cbd7e2",
    field_text="#243746",
    button_bg="#2780e3",
    button_text="#ffffff",
    button_border="#1d68bc",
    button_hover_bg="#3a8bea",
    button_pressed_bg="#1f70ca",
    button_flat_hover="rgba(39, 128, 227, 0.14)",
    item_bg="#ffffff",
    item_border="#d8e1e8",
    status_bg="#edf3f8",
    menu_bg="#ffffff",
    menu_hover_bg="#2780e3",
    toolbar_bg="#eef4f8",
    toolbar_border="#c8d5df",
    toolbar_hover_bg="#ffffff",
    toolbar_checked_bg="#dde8f2",
    accent="#2780e3",
    accent_text="#ffffff",
    table_bg="#ffffff",
    table_alt_bg="#f4f8fb",
    table_border="#d4dee7",
    table_grid="#e2e9ef",
    table_text="#243746",
    table_row_border="#e9eff4",
    header_bg="#e5edf3",
    header_text="#274055",
    header_border="#ccd8e2",
    check_border="#a8b9c9",
    db_status_bg="#f0f5f9",
    db_status_text="#486173",
    danger="#c13f31",
)


SANDSTONE_THEME = ThemeSpec(
    name="sandstone",
    window="#f5f1e8",
    window_text="#3e372f",
    central_start="#fbf9f4",
    central_end="#efe8db",
    pane="#fdfcf9",
    pane_border="#d7ccbd",
    tab_idle_bg="#e4d8c3",
    tab_idle_text="#5b5146",
    tab_active_bg="#fdfcf9",
    tab_active_text="#325d88",
    panel_bg="rgba(253, 252, 249, 0.94)",
    panel_border="#ddd2c3",
    group_title="#6b5f53",
    tree_bg="#fdfcf9",
    tree_alt_bg="#f5efe4",
    tree_hover_bg="#ebe2d0",
    tree_border="#d7ccbd",
    field_bg="#ffffff",
    field_border="#d0c5b6",
    field_text="#433b33",
    button_bg="#325d88",
    button_text="#ffffff",
    button_border="#294f74",
    button_hover_bg="#3b6794",
    button_pressed_bg="#294f74",
    button_flat_hover="rgba(50, 93, 136, 0.14)",
    item_bg="#fdfcf9",
    item_border="#ddd3c6",
    status_bg="#ece4d8",
    menu_bg="#fdfcf9",
    menu_hover_bg="#325d88",
    toolbar_bg="#e8dece",
    toolbar_border="#cdbfae",
    toolbar_hover_bg="#f5efe4",
    toolbar_checked_bg="#ddd0bc",
    accent="#325d88",
    accent_text="#ffffff",
    table_bg="#ffffff",
    table_alt_bg="#f7f2ea",
    table_border="#d7ccbd",
    table_grid="#e5ddd0",
    table_text="#433b33",
    table_row_border="#ebe3d8",
    header_bg="#e8dece",
    header_text="#4a4138",
    header_border="#d0c4b4",
    check_border="#b9ab98",
    db_status_bg="#f3ede3",
    db_status_text="#62574c",
    danger="#b94a48",
)


UNITED_THEME = ThemeSpec(
    name="united",
    window="#f8f4ef",
    window_text="#3a3026",
    central_start="#fffdfa",
    central_end="#f2e7da",
    pane="#ffffff",
    pane_border="#decfc0",
    tab_idle_bg="#eadbc8",
    tab_idle_text="#5b4838",
    tab_active_bg="#ffffff",
    tab_active_text="#e95420",
    panel_bg="rgba(255, 255, 255, 0.94)",
    panel_border="#e2d3c5",
    group_title="#715846",
    tree_bg="#ffffff",
    tree_alt_bg="#faf4ed",
    tree_hover_bg="#f4e4d3",
    tree_border="#decfc0",
    field_bg="#ffffff",
    field_border="#d8c9ba",
    field_text="#3f342a",
    button_bg="#e95420",
    button_text="#ffffff",
    button_border="#c9461c",
    button_hover_bg="#ef6736",
    button_pressed_bg="#cd4a1d",
    button_flat_hover="rgba(233, 84, 32, 0.12)",
    item_bg="#ffffff",
    item_border="#e4d6c8",
    status_bg="#f2e6d8",
    menu_bg="#ffffff",
    menu_hover_bg="#e95420",
    toolbar_bg="#f0e2d2",
    toolbar_border="#d8c6b6",
    toolbar_hover_bg="#fbf2e8",
    toolbar_checked_bg="#ead7c5",
    accent="#e95420",
    accent_text="#ffffff",
    table_bg="#ffffff",
    table_alt_bg="#fbf5ee",
    table_border="#decfc0",
    table_grid="#ede1d5",
    table_text="#3f342a",
    table_row_border="#efe5db",
    header_bg="#eadbc8",
    header_text="#493b2f",
    header_border="#dcccbc",
    check_border="#bda895",
    db_status_bg="#f6ece1",
    db_status_text="#6a5542",
    danger="#c73e1c",
)


THEMES = {
    LIGHT_THEME.name: LIGHT_THEME,
    DARK_THEME.name: DARK_THEME,
    OCEAN_THEME.name: OCEAN_THEME,
    FOREST_THEME.name: FOREST_THEME,
    CERULEAN_THEME.name: CERULEAN_THEME,
    BRITE_THEME.name: BRITE_THEME,
    COSMO_THEME.name: COSMO_THEME,
    SANDSTONE_THEME.name: SANDSTONE_THEME,
    UNITED_THEME.name: UNITED_THEME,
}

THEME_LABELS = {
    "light": "Light Studio",
    "dark": "Dark Slate",
    "ocean": "Ocean Mist",
    "forest": "Forest Paper",
    "cerulean": "Cerulean Sky",
    "brite": "Brite Pop",
    "cosmo": "Cosmo Blue",
    "sandstone": "Sandstone Calm",
    "united": "United Ember",
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
