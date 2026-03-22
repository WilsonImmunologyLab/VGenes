# VGenes Development Log

## Purpose

This file tracks the current modernization goals for VGenes, the rationale behind them, and the progress made while refining the project.

## Current Assessment

VGenes is a substantial PyQt desktop application with strong domain functionality, but the current codebase has several structural constraints:

- The main application logic is concentrated in `VGenesMain.py`, which is very large and tightly coupled.
- The data model is stored in SQLite, but database access is spread across the codebase and often uses string-built SQL.
- The UI is feature-rich but relies heavily on `QTableWidget`, modal dialogs, and inconsistent styling.
- Packaging and dependency management are incomplete.
- Some code depends on legacy Biopython APIs that complicate dependency upgrades.

## Project Goals

### Goal 1: Improve Data Load and Storage

Keep SQLite as the primary editable store for now, but make it safer, faster, and easier to maintain.

Target outcomes:

- Centralize database access behind a small service or repository layer.
- Replace string-built SQL with parameterized queries.
- Add indexes for common filtering and sorting paths.
- Reduce dependence on `os.chdir()` for database and file operations.
- Evaluate a hybrid model:
  - SQLite for editable application state.
  - Parquet for large analytical exports and cached report datasets.
  - DuckDB only if reporting/query workloads later outgrow SQLite.

### Goal 2: Make the UI More Modern

Improve the interface without rewriting the entire application.

Target outcomes:

- Move high-volume tables from `QTableWidget` to `QTableView` with models.
- Create a small, centralized visual system for colors, spacing, and widget states.
- Reduce dialog sprawl by favoring persistent panels or clearer task flows.
- Improve progress feedback, empty states, error states, and visual hierarchy.
- Keep the domain-heavy workflows intact while making the product feel more polished.

### Goal 3: Update Core Packages Safely

Modernize dependencies in a controlled way, without breaking biological sequence workflows.

Target outcomes:

- Add a real dependency manifest and pin a known-good environment.
- Audit imports that depend on removed or deprecated Biopython APIs.
- Upgrade packages in stages rather than all at once.
- Separate "safe package updates" from "refactor-required updates".

## Recommended Technical Direction

### Data Storage Recommendation

Do not replace SQLite immediately.

Why:

- It fits a desktop app well.
- It supports local, portable, queryable storage.
- The current problem is mostly schema/design and access patterns, not the database engine itself.

Recommended direction:

1. Keep SQLite as the source of truth.
2. Add migration/version handling for schema changes.
3. Normalize a few repeated or metadata-heavy areas where it improves maintainability.
4. Add Parquet export paths for large analysis outputs.
5. Revisit DuckDB later only for heavy read-only analytics.

### UI Recommendation

Do not attempt a full UI rewrite first.

Recommended direction:

1. Improve the existing PyQt UI architecture.
2. Replace the most expensive tables with model/view implementations.
3. Add centralized theming and reusable styled widgets.
4. Modernize the main workflow screens before polishing secondary dialogs.

### Package Upgrade Recommendation

Upgrade in phases.

Safe-to-start areas:

- dependency manifest
- packaging metadata
- version pinning
- import auditing

Refactor-first areas:

- legacy Biopython usage such as `Bio.Alphabet`, `Bio._py3k`, and `pairwise2`
- code tightly coupled to generated Qt UI code and old widget patterns

## Phased Roadmap

### Phase 1: Baseline and Safety

- Create and maintain this development log.
- Add a dependency manifest.
- Document current runtime assumptions.
- Audit data access entry points.
- Identify the highest-traffic UI tables and dialogs.

### Phase 2: Data Layer Hardening

- Introduce a dedicated database module boundary.
- Replace string-concatenated SQL in the most critical paths.
- Add indexes and measure load/filter performance.
- Remove `os.chdir()` from DB and file access flows where practical.
- Add schema version checks.

### Phase 3: UI Foundation

- Introduce a central theme/style module.
- Convert the main sequence table to model/view architecture.
- Standardize status/progress/error messaging.
- Simplify repeated dialog patterns.

### Phase 4: Dependency Modernization

- Create an installable/pinned environment definition.
- Refactor Biopython-dependent modules away from removed APIs.
- Upgrade secondary packages once compatibility is clear.

### Phase 5: Packaging and Developer Experience

- Replace placeholder packaging metadata.
- Document local setup, build steps, and tool dependencies.
- Clarify what is source, generated code, vendor code, and runtime data.

## Immediate Next Steps

The most practical next slice of work is:

1. create a dependency manifest and environment baseline
2. map the current database read/write hotspots
3. identify which table should be migrated first from `QTableWidget` to `QTableView`
4. start a small database access wrapper instead of editing SQL inline everywhere

## Dependency Baseline

The project now has an explicit dependency baseline in:

- `requirements.txt` for the core desktop runtime
- `requirements-optional.txt` for specialized workflows

Current guidance:

- Prefer a dedicated virtual environment rather than the checked-in `venv/`.
- Treat Biopython as pinned until legacy imports are removed from the codebase.
- Keep optional scientific and immunology-specific tooling separate from the
  base GUI install.

## Progress Log

### 2026-03-10

- Added a small reusable async task layer in `task_utils.py` with:
  - `TaskSignals`
  - `FunctionTask`
  - a background fetch helper for paged DB-table loading
- Migrated the Database-tab table load path in `VGenesMain.py` off the old
  `QThread` UI-mutation pattern and onto:
  - background data fetch on `QThreadPool`
  - main-thread table application
  - token-guarded progress dialog lifecycle
- Kept the old `LoadTable_thread` class in place for now as a legacy reference,
  but the active DB-table loading path no longer depends on it.

### 2026-03-08

- Reviewed repository structure and core application files.
- Identified `VGenesMain.py` as the central architectural bottleneck.
- Evaluated whether SQLite should be replaced and concluded that it should be retained for now.
- Identified a hybrid storage direction: SQLite for editing, Parquet for analysis/export, DuckDB only if needed later.
- Identified UI modernization priorities around model/view tables, centralized theming, and workflow cleanup.
- Identified package upgrade constraints, especially legacy Biopython usage.
- Created `develop.md` to track goals, roadmap, and progress.
- Added `requirements.txt` with a conservative core runtime baseline.
- Added `requirements-optional.txt` for specialized workflows such as Change-O,
  pRESTO, clustering, and Mongo-backed utilities.
- Added `db_utils.py` as a small shared SQLite helper for connection lifecycle
  and identifier validation.
- Refactored key `VGenesSQL.py` paths to stop depending on `os.chdir()` for
  connection setup.
- Parameterized direct field updates in `VGenesSQL.py` for `UpdateField`,
  `UpdateFieldTable`, and `UpdateFieldbySeqName`.
- Moved selected utility functions such as `RunSQL`, `RunUpdateSQL`,
  `ImportDB`, and `DumpDB` onto the shared connection helper.
- Added helper-backed batch edit primitives in `VGenesSQL.py` for:
  - fetching record IDs by sequence names
  - reading a single field by record ID
  - updating a field where it matches an old value
- Migrated the batch field update and bulk replace flows in
  `VGenesMain.py` away from inline update SQL and onto the safer helpers.
- Added a structured multi-field update helper in `VGenesSQL.py` for saving
  edited records by sequence name.
- Migrated `on_btnSaveChange_clicked` in `VGenesMain.py` away from a large
  hand-built `UPDATE vgenesDB SET ...` string to the helper-based update path.
- Migrated additional repeated update flows in `VGenesMain.py` to helper-backed
  operations, including:
  - range/mapping batch updates
  - shared clone annotation writes
  - pair rename barcode writes
  - Change-O clonal pool imports
  - clone ranking and duplicate-marking writes
  - isotype writeback
  - subgroup append updates
- Migrated another batch of dynamic field updates to helper-backed operations,
  including:
  - barcode modification threads updating dynamic fields by `ID`
  - selection-based mark/update actions
  - whole-column field copy operations
  - clonal pool reset-by-selection updates
  - reanalysis writeback of many fields per record
- Added helper-backed anchored multi-field updates for CSV/table import style
  workflows that update records by a chosen anchor field, with optional HC/LC
  filtering.
- The remaining `UPDATE vgenesDB SET ...` occurrences in `VGenesMain.py` are
  now in commented legacy code blocks rather than active paths.
- Began the UI/data presentation refactor on the main sequence table by
  replacing per-row `QCheckBox` cell widgets with native checkable table items.
- Added a small internal helper layer around `SeqTable` for row-name lookup,
  checked-state handling, and selected-row handling.
- Added SQL-backed persistent sort behavior for the main sequence table so
  paging and tree-to-table navigation stay consistent with the active sort.
- Extracted the main sequence table state and row-population logic into
  `seq_table_adapter.py`, reducing direct `QTableWidget` manipulation spread
  across `VGenesMain.py` and setting up the next model/view migration step.
- Moved the main table's selection, current-row focus, cell text access, and
  edit-trigger handling further behind `seq_table_adapter.py`, so navigation
  and edit rollback paths now depend less on direct widget item APIs.
- Replaced the main sequence grid at runtime with a model-backed
  `SequenceTableView` wrapper, so the `SeqTable` screen now runs on Qt's
  model/view architecture while preserving the existing app-level behavior.
- Added a dedicated DB-tab theme module and applied it to the main sequence
  grid, pager controls, and selection controls so the modernized table has a
  more coherent visual surface without restyling the whole application.
- Added DB-tab status feedback and an empty-state placeholder for the main
  sequence grid so users now get explicit messages when the table is hidden,
  loading, empty, or successfully populated.
- Moved main-table paging/edit/header metadata fully into
  `seq_table_adapter.py` and added indexed row-name lookup there, reducing
  direct widget-state coupling and avoiding repeated linear row scans when
  focusing a sequence in the current page.
- Replaced the main table's old ad hoc `last_name` item attribute with an
  explicit Qt model data role in `sequence_table_view.py`, making row metadata
  storage cleaner and more native to the model/view implementation.
- Fixed duplicate-name edit handling by making
  `VGenesSQL.UpdateFieldbySeqName` raise database errors instead of silently
  swallowing them, so the main-table warning/rollback path now executes.
- Hardened the `fieldsname` display editor by replacing its per-row raw
  `UPDATE fieldsname ...` statements with a helper-backed batch save path and
  removing a duplicate table reload inside that dialog.
- Moved the alter dialog's `fieldsname` insert/delete operations onto helper
  functions in `VGenesSQL.py` and fixed its fixed-field deletion guard to read
  the actual type text instead of comparing the widget object itself.
- Migrated the remaining active field-creation import paths to the same
  `fieldsname` helper functions, so CSV/VDB-driven metadata creation now uses
  the same insert logic as the alter dialog instead of assembling raw insert
  SQL in multiple places.
- Replaced the placeholder packaging metadata with a real `pyproject.toml`,
  reduced `setup.py` to a setuptools compatibility shim, expanded the README
  with developer setup guidance, and tightened `.gitignore` for common Python
  build/cache artifacts.
- Started the Biopython-compatibility refactor by removing direct
  `Bio.Alphabet` and `Bio._py3k` imports from the checked-in local sequence
  compatibility modules (`Seq.py` and `CodonTable.py`), redirecting them to
  the local `Alphabet`/`IUPAC` implementations instead.
- Removed active `pairwise2` usage by switching the main alignment workflow to
  `Bio.Align.PairwiseAligner` and simplifying `Clonify.py`'s equal-length
  distance path to direct match counting, while also fixing the most immediate
  Python 2 syntax blockers in `Clonify.py` so it now parses under Python 3.
- Updated the dependency baseline to `biopython==1.86` in
  `requirements.txt`, reflecting the removal of the main direct
  `Bio.Alphabet`, `Bio._py3k`, and active `pairwise2` blockers from the
  maintained code paths.
- Validated the upgraded dependency stack in the local environment:
  - installed the updated core requirements
  - installed the optional workflow dependencies
  - corrected `requirements-optional.txt` so `fastcluster` stays compatible
    with the pinned `numpy<2` baseline
  - added a `pyqtgraph` template compatibility fallback in
    `PyQtGraphPlotItem.py`
  - broke a `VGenesMain`/`VReports` circular import with lazy UI-class imports
  - confirmed `import VGenesMain` succeeds under the updated dependency set
- The current local Anaconda environment still emits duplicate Qt-library
  warnings because it mixes Conda Qt libraries with pip-installed PyQt wheel
  runtimes; a clean project virtual environment is still the recommended
  runtime target.
- Added small runtime-hygiene fixes in `VGenesMain.py` so the app now creates
  its `Temp/` directory before writing logs and uses that directory as
  `MPLCONFIGDIR`, avoiding Matplotlib cache warnings when the user home cache
  path is unavailable or unwritable.
- Performed a conservative repository reorganization pass:
  - moved sample/demo `.vdb` databases out of the root into `samples/vdb/`
  - moved archival docs/spreadsheets/build notes into `reference/`
  - moved one-off examples and legacy utility scripts into `dev-scripts/`
  - intentionally left runtime-coupled files such as `BackUP.vdb`,
    `UpdateRecord.nt`, root icon assets, and active application modules in
    place to avoid breaking import/spec/runtime paths during this pass
- Updated the main table load/check/edit wiring to use the lighter checkable
  item approach instead of embedded checkbox widgets.
- Added persistent main-table sort state and moved page loading order to SQL,
  so sorting can be applied consistently across paged data instead of only
  sorting the currently loaded widget rows in memory.
- This improves table loading overhead and simplifies a future migration from
  `QTableWidget` toward a fuller model/view approach.
- Reached the first manual-test checkpoint for database hardening. Remaining
  raw SQL updates are now concentrated in narrower specialized workflows rather
  than the most common editing paths.
- Integrated the active Change-O path into the application:
  - added `changeo_runner.py` to run IgBlast, `MakeDb.py`, and
    `DefineClones.py` end-to-end inside VGenes
  - wired the Change-O task onto the shared `FunctionTask` / progress-dialog
    async path
  - updated the Change-O dialog so users can run the internal pipeline
    directly while still keeping an external result-import fallback
  - used a default `DefineClones.py --dist 0.15` threshold as the initial
    built-in Change-O clone distance setting
- Split Change-O UI into two layers:
  - preserved the legacy `ChangeODialog` for the original manual/export flow
  - added a new `ChangeORunDialog` as the default entry point for integrated
    Change-O execution and parameter tuning
- Migrated the active HC/LC barcode-pairing workflow onto the shared async
  task model:
  - added `run_hclc_pairing(...)` in `task_utils.py`
  - replaced the active `HCLC_pair_thread` launch paths for
    `match HC/LC` and `pair-and-jump` with `FunctionTask` + shared progress
    dialog handling
  - corrected the worker to use the current checkout's `VGenesSQL.RunSQL(...)`
    API instead of a newer `run_sql_query(...)` helper that is not present in
    this branch
- Migrated the active "pair to new DB" HC/LC export flow onto the shared
  async task model with `run_hclc_export_db(...)` and a main-thread result
  adapter that preserves the existing `HCLC_finish_process(...)` behavior
- Migrated the active "record to new DB" export flow onto the shared async
  task model with `run_copy_records_db(...)`, reusing the same
  `HCLC_finish_process(...)` post-run UI path
- Fixed HC/LC async check-state synchronization so post-task updates now
  deduplicate `CheckedRecords` and run the normal tree/table sync path again,
  restoring correct "next checked" navigation after `match HC/LC`
- Migrated the active alignment-HTML workflow onto the shared async task
  model:
  - added `run_alignment_html_task(...)` in `VGenesMain.py`
  - replaced the active alignment HTML launch paths for the main checked-list,
    HC table, LC table, and clone alignment HTML actions with
    `FunctionTask` + shared progress handling
- Hardened `AlignSequencesHTMLBCR(...)` against non-numeric ruler coordinates
  such as `NA`, preventing clone-page alignment HTML generation from crashing
  on incomplete region metadata
- Migrated the active tree-HTML workflow onto the shared async task model:
  - added `run_tree_html_task(...)` in `VGenesMain.py`
  - replaced the active tree HTML launch paths for HC, LC, selected records,
    and clone-tree generation with `FunctionTask` + shared progress handling
- Migrated the active sequence-logo workflow onto the shared async task
  model:
  - added `run_seq_logo_task(...)` in `VGenesMain.py`
  - replaced the active main-logo and clone-logo launch paths with
    `FunctionTask` + shared progress handling
- Migrated the active sequence-similarity heatmap workflow onto the shared
  async task model:
  - added `run_seq_similarity_task(...)` in `VGenesMain.py`
  - replaced the active checked-record and clone similarity launch paths with
    `FunctionTask` + shared progress handling
- Migrated the active Pattern Search workflow onto the shared async task
  model:
  - added `run_pattern_search_task(...)` in `VGenesMain.py`
  - replaced the active `PatternThread` launch path in `SearchPattern(...)`
    with `FunctionTask` + shared progress handling
  - kept the existing pattern-result dialog rendering on the main thread
- Removed the Pattern Search IgBlast dependency from the active workflow:
  - the search now derives AA regions directly from stored V/D/J, FR/CDR, and
    junction annotations in the database
  - this preserves the async/progress flow while avoiding a full IgBlast run
    for each search
- Migrated the active Cookie/representative sampling workflow onto the shared
  async task model:
  - added `run_cookie_sampling_task(...)` in `VGenesMain.py`
  - added a shared `launchCookieSamplingTask(...)` path
  - replaced the four active `CookieThread` launch sites with
    `FunctionTask` + shared progress handling
  - kept the existing `cookieRes(...)` table rendering path on the main thread
- Fixed two sampling regressions:
  - moved the new Cookie async launcher onto `SamplingDialog`, which is where
    the representative sampling workflow actually runs
  - replaced proportional stratified size rounding with an exact allocation
    helper so requested totals are preserved much more closely instead of
    collapsing to per-level minimums
- Fixed Cookie sampling environment/messaging issues:
  - corrected the large-sample warning text so PF mode no longer claims it is
    running "without prime factor"
  - added `pyclustering` to `requirements-optional.txt`
  - installed `pyclustering` locally and added a clearer runtime error message
    if that optional dependency is missing
- Migrated the active protein-similarity workflow onto the shared async task
  model:
  - added `run_protein_similarity_task(...)` in `VGenesMain.py`
  - replaced the dialog's `protein_slimlar_thread` launch path with
    `FunctionTask` + main-thread result/error handling
  - kept `ShowProteinSimilarResults(...)` as the main-thread result renderer
- Migrated the active CSV/VDB merge-import workflows onto the shared async
  task model:
  - added `run_csv_import_task(...)` and `run_vdb_import_task(...)` in
    `VGenesMain.py`
  - replaced the active `CSV_thread` / `VDB_thread` launch paths in
    `ImportDataDialogue` with `FunctionTask` + main-thread result/error
    handling
- Migrated the active IgBlast/IMGT parser launch paths in `ImportDataDialogue`
  onto the shared async task model:
  - added `run_igblast_parser_task(...)`, `run_igblast_results_task(...)`,
    and `run_imgt_parser_task(...)` in `VGenesMain.py`
  - replaced the active `WorkThread`, `WorkThread1`, and
    `WorkThreadIMGTparser` launch paths used by the 10X/FASTA/SEQ import,
    IgBlast-results import, and IMGT import flows
  - kept the existing `multi_callback(...)` / `multiIMGT_callback(...)`
    database-insert logic intact by routing parsed results back through the
    same callbacks on the main thread

## Notes

The working tree already contains unrelated local modifications. Future code changes should be scoped carefully and recorded here as they are completed.
