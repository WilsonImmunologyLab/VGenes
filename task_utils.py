import math
import traceback

from PyQt5.QtCore import QObject, QRunnable, Qt, pyqtSignal, pyqtSlot

import VGenesSQL


class TaskSignals(QObject):
    started = pyqtSignal(str)
    progress = pyqtSignal(int, str)
    result = pyqtSignal(object)
    error = pyqtSignal(str, str)
    finished = pyqtSignal()


class FunctionTask(QRunnable):
    def __init__(self, fn, *args, task_name="Task", **kwargs):
        super().__init__()
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.task_name = task_name
        self.signals = TaskSignals()

    @pyqtSlot()
    def run(self):
        try:
            self.signals.started.emit(self.task_name)
            result = self.fn(*self.args, emit_progress=self.signals.progress.emit, **self.kwargs)
        except Exception as exc:
            self.signals.error.emit(str(exc), traceback.format_exc())
        else:
            self.signals.result.emit(result)
        finally:
            self.signals.finished.emit()


def build_seq_table_order_clause(sort_field, sort_order, default_sort_fields):
    order_fields = []
    if sort_field:
        direction = "ASC" if sort_order == Qt.AscendingOrder else "DESC"
        order_fields.append(f"{sort_field} {direction}")

    for field in default_sort_fields:
        if field in ("", None, "None"):
            continue
        if any(entry.split()[0] == field for entry in order_fields):
            continue
        order_fields.append(field)

    if not order_fields:
        order_fields.append("SeqName")
    return " ORDER BY " + ", ".join(order_fields)


def fetch_seq_table_page(
    db_filename,
    page_size,
    current_page_index,
    sort_field,
    sort_order,
    default_sort_fields,
    emit_progress=None,
):
    if emit_progress is None:
        emit_progress = lambda pct, label: None

    if not db_filename or db_filename == "none":
        return {"status": "no_db"}

    emit_progress(5, "Counting records ...")
    total_records = VGenesSQL.RunSQL(db_filename, "SELECT COUNT(*) FROM vgenesDB")[0][0]
    total_pages = max(1, math.ceil(total_records / page_size))

    if total_records == 0:
        return {
            "status": "empty",
            "total_records": 0,
            "total_pages": total_pages,
            "page_size": page_size,
            "current_page": 1,
        }

    emit_progress(15, "Loading display fields ...")
    header_rows = VGenesSQL.RunSQL(
        db_filename,
        'SELECT Field,FieldNickName FROM fieldsname WHERE display = "yes" ORDER BY display_priority,ID LIMIT 0,200',
    )
    current_field_list = [row[0] for row in header_rows]
    current_nickname_list = [row[1] for row in header_rows]

    fields = ",".join(current_field_list)
    if fields == "":
        fields = "*"

    emit_progress(40, "Loading page records ...")
    record_limit_statement = f" LIMIT {current_page_index * page_size},{page_size}"
    sql_statement = "select " + fields + " from vgenesDB"
    sql_statement += build_seq_table_order_clause(sort_field, sort_order, default_sort_fields)
    sql_statement += record_limit_statement
    data_rows = VGenesSQL.RunSQL(db_filename, sql_statement)

    emit_progress(90, "Preparing table view ...")
    return {
        "status": "ok",
        "total_records": total_records,
        "total_pages": total_pages,
        "page_size": page_size,
        "current_page": current_page_index + 1,
        "field_list": current_field_list,
        "nickname_list": current_nickname_list,
        "rows": data_rows,
        "row_names": [str(row[0]) for row in data_rows],
    }
