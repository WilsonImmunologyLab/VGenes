<p align="center">
  <img src="/docs/img/logo.png"  width="200">
</p>

# VGenes
VGenes, an integrated graphical platform for multimodal analyses of B cell receptor repertoire

Instructions, documentation, and tutorials can be found at:

https://wilsonimmunologylab.github.io/VGenes/

Source coee of VGenes is hosted on GitHub, you can view and clone the repository at

https://github.com/WilsonImmunologyLab/VGenes

<img src="/docs/img/GraphicalAbstract.png"  width="800">

VGenes enables multimodal analysis of BCR repertoire data along with their transcriptome expression, surface protein expression, and antigen probe binding, by taking advantage of novel multimodal single-cell sequencing technics. By developing a graphical user interface (GUI), VGenes highly reduces the learning cost and allows users to efficiently and conveniently manage, analyze, edit their BCR sequences. In addition, VGenes seamlessly connects BCR sequences with antibody cloning and characterization, integrates a comprehensive collection of BCR-specific functions, aiming to help users better understand B cell repertoire and efficiently select candidates of neutralizing antibody from massive B cell populations. The use of VGenes will significantly facilitate B cell research and antibody candidate selection.

<img src="/docs/img/Overall.png"  width="800">

## Developer Setup

VGenes currently ships as a desktop application source tree rather than a clean
importable Python package. The recommended development workflow is:

1. Create a dedicated virtual environment with Python 3.10 or 3.11.
2. Install the core runtime dependencies:

```bash
python3 -m pip install -r requirements.txt
```

3. Install optional workflow-specific tools only if needed:

```bash
python3 -m pip install -r requirements-optional.txt
```

4. Launch the desktop application directly:

```bash
python3 VGenesMain.py
```

Notes:

- Biopython is pinned in `requirements.txt` to a tested modern release rather
  than left floating, because sequence/alignment compatibility still matters
  for the desktop workflows.
- The application expects bundled runtime assets such as `Tools/`, `IgBlast/`,
  `Resources/`, `Js/`, and `Temp/` to exist relative to the project root.
- Packaging metadata is now defined in `pyproject.toml`, but the runtime is
  still organized as a source-first desktop app rather than a polished library.
- Prefer a clean project virtual environment over a mixed Anaconda + `pip`
  base environment. During dependency validation, the mixed environment loaded
  duplicate Qt libraries and emitted warnings even though imports completed.

## Build Metadata

Basic build metadata now lives in `pyproject.toml`. For packaging work, treat
the current state as a reproducible dependency/install baseline, not a finished
distribution story.
