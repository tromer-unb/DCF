
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# Carregar estrutura do arquivo CIF
structure = Structure.from_file("input.cif")

# Obter célula primitiva
analyzer = SpacegroupAnalyzer(structure, symprec=1e-5)
primitive = analyzer.find_primitive()

# Salvar como novo arquivo CIF
primitive.to(filename="celula_minima.cif")
