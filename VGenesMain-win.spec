# -*- mode: python ; coding: utf-8 -*-

block_cipher = None

added_files = [
             (r'C:\Users\leili\Documents\GitHub\VGenes\Data','Data'),
             (r'C:\Users\leili\Documents\GitHub\VGenes\Conf\*','Conf'),
             (r'C:\Users\leili\Documents\GitHub\VGenes\IgBlast','IgBlast'),
             (r'C:\Users\leili\Documents\GitHub\VGenes\Tools\*','Tools'),
             (r'C:\Users\leili\Documents\GitHub\VGenes\Js\*','Js'),
             (r'C:\Users\leili\Documents\GitHub\VGenes\Temp\ErLog.txt','Temp'),
             (r'C:\Users\leili\Documents\GitHub\VGenes\Temp\ErLog.txt','Temp'),
             (r'C:\Users\leili\Documents\GitHub\VGenes\qtweb-resources','.')
             ]

a = Analysis(['VGenesMain.py'],
             pathex=['C:\\Users\\leili\\Documents\\GitHub\\VGenes'],
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
          name='VGenes',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          console=False,
          icon='mAB.ico')
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               upx_exclude=[],
               name='VGenes')
