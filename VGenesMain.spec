# -*- mode: python ; coding: utf-8 -*-
import sys
sys.setrecursionlimit(5000)

block_cipher = None

added_files = [
             ('/Users/lel4003/Documents/Projects/VGene/VGenes/Js/*','Js'),
             ('/Users/lel4003/Documents/Projects/VGene/VGenes/Conf/path_setting.txt','Conf'),
             ('/Users/lel4003/Documents/Projects/VGene/VGenes/Conf/RecentPaths.vtx','Conf'),
             ('/Users/lel4003/Documents/Projects/VGene/VGenes/Temp/ErLog.txt','Temp'),
             ('/Users/lel4003/Documents/Projects/VGene/VGenes/Temp/ErLog2.txt','Temp'),
             ('/Users/lel4003/Documents/Projects/VGene/VGenes/Tools/raxml','Tools'),
             ('/Users/lel4003/Documents/Projects/VGene/VGenes/Tools/clustalo','Tools'),
             ('/Users/lel4003/Documents/Projects/VGene/VGenes/Tools/muscle','Tools'),
             ('/Users/lel4003/Documents/Projects/VGene/VGenes/Tools/makeblastdb','Tools'),
             ('/Users/lel4003/Documents/Projects/VGene/VGenes/Data','Data'),
             ('/Users/lel4003/Documents/Projects/VGene/VGenes/IgBlast','IgBlast'),
             ]

a = Analysis(['VGenesMain.py'],
             pathex=['/Users/lel4003/Documents/Projects/VGene/VGenes'],
             binaries=[],
             datas=added_files,
             hiddenimports=['scipy.special.cython_special','cmath'],
             hookspath=['hooks'],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)
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
              'NSHumanReadableCopyright':"Copyright @ 2021, Wilson Lab, All Rights Reserved",
              'NSHighResolutionCapable': 'True'
             })