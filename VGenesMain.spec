# -*- mode: python ; coding: utf-8 -*-

block_cipher = None

added_files = [
             ('/Users/leil/Documents/Projects/VGene/Resources/Js/echarts.min.js','Js'),
             ('/Users/leil/Documents/Projects/VGene/Resources/Js/jquery.js','Js')
             ]

a = Analysis(['VGenesMain.py'],
             pathex=['/Users/leil/Documents/Projects/VGene/VGenes'],
             binaries=[],
             datas=[],
             hiddenimports=[],
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
          a.binaries,
          a.zipfiles,
          a.datas,
          [],
          name='VGenes',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          upx_exclude=[],
          runtime_tmpdir=None,
          console=False )
app = BUNDLE(exe,
             name='VGenes.app',
             icon=None,
             bundle_identifier=None)
