import py3Dmol #visualize sample protein 
view = py3Dmol.view(query='pdb:1kzy')
view.setStyle({'cartoon':{'color':'spectrum'}})
view
