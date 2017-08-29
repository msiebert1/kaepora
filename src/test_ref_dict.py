def build_bsnip_ref_dict():
	with open('../data/info_files/bsnip_references.txt') as f:
		txt = f.readlines()
		ref_dict = {}
		for line in txt:
			info = line.split()
			for el in info:
				if '.flm' in el:
					ref_dict[el] = info[-1]

	return ref_dict