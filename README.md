# DSSP_IMPL
Secondary structures assignment using DSSP as reference
Installation
Requirements
•	Python 3.6
•	Biopython package
pip install biopython
Clone the repository
git clone https://github.com/kabhel/DSSP_implt.git
cd DSSP_implt
chmod +x script.py

hydratation: 
this algorithm does not hydrate proteins, so it does not add hydrogen to the level of pdb structures. for this the site http://molprobity.biochem.duke.edu/ is recommended for hydration. all you have to do is enter the site and upload the pdb file, then after clicking on upload the site will indicate several options for the file, you just have to indicate that you want to add hydrogen to have the pdb hydrated.
Cet algorithme n'hydrate pas les protéines, donc il n'ajoute pas les hydrogène au niveau des structure pdb. Pour cela le site http://molprobity.biochem.duke.edu/ est recommandée pour l'hydratation. Il suffit de rentre dans le site et d'y déposer le fichier pdb puis après avoir cliqué sur upload le site vous indiquera plusieurs options pour le fichier il suffit d'indiquer qu’on veut ajouter des hydrogènes pour avoir le pdb hydraté.
Run the program
Examples
Python script.py 1bta.pdb A
The first argument is the protein in question and the second argument is the string in question. The first argument must be a pdb file for the algorithm to work.
Le premier argument est la protéine en question et le deuxième argument est la chaîne en question. Le premier argument doit être un fichier pdb pourque l'algorithme fonctionne.

Reference
W. Kabsch and C.Sander, Biopolymers 22 (1983) 2577-2637
