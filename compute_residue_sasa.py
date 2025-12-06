def compute_residue_sasa(pdb_file, WATER_NAMES, ACCESS_THRESHOLD):
    import freesasa

    structure = freesasa.Structure(pdb_file)
    result = freesasa.calc(structure)

    res_dict = result.residueAreas()
    out = []

    for chain, chain_dict in res_dict.items():
        for resnum, area_obj in chain_dict.items():
            resname = area_obj.residueType

            if resname.upper() in WATER_NAMES:
                continue

            sasa = area_obj.total
            accessible = "Yes" if sasa >= ACCESS_THRESHOLD else "No"
            out.append([chain, str(resnum), resname, sasa, accessible])
    return out
