import numpy as np 
dark_type = np.dtype([
    ('mass', '>f4'),
    ('r', '>f4', (3,)),
    ('v', '>f4', (3,)),
    ('eps', '>f4'),
    ('phi', '>f4')
])



# (i) 
print("Minimum coordinates (x, y, z):", np.min(dark['r'], axis=0))
print("Maximum coordinates (x, y, z):", np.max(dark['r'], axis=0))

# (ii)
print("Total mass:", np.sum(dark['mass']))