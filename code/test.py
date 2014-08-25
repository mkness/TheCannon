name = 'self_tags'
file_in2 = open(name+'.pickle', 'r') 
params, icovs_params = pickle.load(file_in2)
file_in2.close()

fn = 'coeffs.pickle'
fd = open(fn, "r")
dataall, metaall, labels, offsets, coeffs, covs, scatters = pickle.load(fd) 
fd.close()
