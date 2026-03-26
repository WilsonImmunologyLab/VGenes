# -*- mode: python ; coding: utf-8 -*-
import sys
sys.setrecursionlimit(5000)

block_cipher = None


def _is_qml_payload(entry):
    text = " ".join(str(part) for part in entry) if isinstance(entry, (tuple, list)) else str(entry)
    normalized = text.replace("\\", "/")
    return "/Qt5/qml/" in normalized or normalized.endswith("/Qt5/qml") or "QtQml/WorkerScript.2" in normalized

added_files = [
             ('./Js/*','Js'),
             ('./Conf/path_setting.txt','Conf'),
             ('./Conf/RecentPaths.vtx','Conf'),
             ('./Temp/ErLog.txt','Temp'),
             ('./Temp/ErLog2.txt','Temp'),
             ('./Tools/raxml','Tools'),
             ('./Tools/clustalo','Tools'),
             ('./Tools/muscle','Tools'),
             ('./Tools/makeblastdb','Tools'),
             ('./Data','Data'),
             ('./IgBlast','IgBlast'),
             ]

a = Analysis(['VGenesMain.py'],
             pathex=['.'],
             binaries=[],
             datas=added_files,
             hiddenimports=['scipy.special.cython_special','cmath'],
             hookspath=['hooks'],
             runtime_hooks=[],
             excludes=[
                 'PyQt5.QtQml',
                 'PyQt5.QtQuick',
                 'PyQt5.QtQuickWidgets',
             ],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)
a.datas = [entry for entry in a.datas if not _is_qml_payload(entry)]
a.binaries = [entry for entry in a.binaries if not _is_qml_payload(entry)]
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          [],
          exclude_binaries=True,
          name='VGene',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          console=False )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               upx_exclude=[],
               name='VGenes')
app = BUNDLE(coll,
             name='VGenes.app',
             icon='mAB.icns',
             bundle_identifier=None,
             info_plist={
              'NSHumanReadableCopyright':"Copyright @ 2023, Wilson Lab, All Rights Reserved",
              'NSHighResolutionCapable': 'True'
             })
