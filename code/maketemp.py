J,H,K = loadtxt("starsin_new_all_ordered.txt", usecols = (13,14,15), unpack =1) 
JmK = J-K
JmKo = JmK 
a1 = 0.5816
a2 = 0.9134
a3 = -0.1443
eqn = [a1 + a2*X + a3*X**2 for X in JmKo]
temp_alonso = [5040/a for a in eqn]
savetxt("Temperature_Alonso.txt", temp_alonso, fmt = "%s") 
