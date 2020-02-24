__author__ = 'wilsonp'
def ModuleFind(self):
    from modulefinder import ModuleFinder
    from VGenesDialogues import openFile
    import os


    Pathname = openFile(self, 'CSV')

    workingdir, filename = os.path.split(Pathname)


    os.chdir(workingdir)
    # import filename

    finder = ModuleFinder()
    finder.run_script(filename)

    Doc = 'Loaded modules:\n'
    for name, mod in finder.modules.items():
        # Doc('%s: ' % name, end=''\n)
        Doc += (','.join(list(mod.globalnames.keys())[:3]))
        Doc += '\n'
        # Doc
    Doc +=('-'*50)
    # Doc +=('Modules not imported:')
    # Doc +=('\n'.join(finder.badmodules.keys()))

