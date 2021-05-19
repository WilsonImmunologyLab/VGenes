# -*- mode: python ; coding: utf-8 -*-
import sys
sys.setrecursionlimit(5000)

block_cipher = None

added_files = [
             ('/Users/leil/Documents/Projects/VGene/VGenes/Js/echarts.js','Js'),
             ('/Users/leil/Documents/Projects/VGene/VGenes/Js/jquery.js','Js'),
             ('/Users/leil/Documents/Projects/VGene/VGenes/Js/bootstrap-theme.min.css','Js'),
             ('/Users/leil/Documents/Projects/VGene/VGenes/Js/bootstrap.min.css','Js'),
             ('/Users/leil/Documents/Projects/VGene/VGenes/Js/bootstrap.min.js','Js'),
             ('/Users/leil/Documents/Projects/VGene/VGenes/Js/d3.js','Js'),
             ('/Users/leil/Documents/Projects/VGene/VGenes/Js/phylotree.css','Js'),
             ('/Users/leil/Documents/Projects/VGene/VGenes/Js/phylotree.js','Js'),
             ('/Users/leil/Documents/Projects/VGene/VGenes/Js/underscore-min.js','Js'),
             ('/Users/leil/Documents/Projects/VGene/VGenes/Js/underscore-min.map','Js'),
             ('/Users/leil/Documents/Projects/VGene/VGenes/Conf/path_setting.txt','Conf'),
             ('/Users/leil/Documents/Projects/VGene/VGenes/Conf/RecentPaths.vtx','Conf'),
             ('/Users/leil/Documents/Projects/VGene/VGenes/Temp/ErLog.txt','Temp'),
             ('/Users/leil/Documents/Projects/VGene/VGenes/Temp/ErLog2.txt','Temp'),
             ('/Users/leil/Documents/Projects/VGene/VGenes/Tools/raxml','Tools'),
             ('/Users/leil/Documents/Projects/VGene/VGenes/Tools/clustalo','Tools'),
             ('/Users/leil/Documents/Projects/VGene/VGenes/Tools/muscle','Tools'),
             ('/Users/leil/Documents/Projects/VGene/VGenes/Tools/makeblastdb','Tools'),
             ('/Users/leil/Documents/Projects/VGene/VGenes/Data/template.html','Data'),
             ('/Users/leil/Documents/Projects/VGene/VGenes/Data/template_raxml_tree.html','Data'),
             ('/Users/leil/Documents/Projects/VGene/VGenes/Data/template_newick_tree.html','Data'),
             ('/Users/leil/Documents/Projects/VGene/VGenes/Data/GibsonConnectors.txt','Data'),
             ('/Users/leil/Documents/Projects/VGene/VGenes/Data/HtConnectors.txt','Data'),
             ('/Users/leil/Documents/Projects/VGene/VGenes/Data/AbVecConnectors.txt','Data'),
             ('/Users/leil/Documents/Projects/VGene/VGenes/Data/template_tree.html','Data'),
             ('/Users/leil/Documents/Projects/VGene/VGenes/Data/VDJGenes.db','Data'),
             ('/Users/leil/Documents/Projects/VGene/VGenes/IgBlast','IgBlast'),
             ]

a = Analysis(['VGenesMain.py'],
             pathex=['/Users/leil/Documents/Projects/VGene/VGenes'],
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
               name='VGene')
app = BUNDLE(coll,
             name='VGene.app',
             icon='mAB.icns',
             bundle_identifier=None,
             info_plist={
              'NSHumanReadableCopyright':"Copyright @ 2021, Wilson Lab, All Rights Reserved",
              'NSHighResolutionCapable': 'True'
             })