import basis_set_exchange as bse

basisSet = "aug-cc-pCVTZ"
fmt = "cfour"

basis = bse.api.get_basis(basisSet, elements=[17], fmt=None)
print(basis)
