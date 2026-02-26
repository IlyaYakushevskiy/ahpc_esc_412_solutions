#!/usr/bin/env python
import numpy as np
from sys import argv, exit

header_type = np.dtype([
  ('time', '>f8'),('N', '>i4'), ('Dims', '>i4'),
  ('Ngas', '>i4'), ('Ndark', '>i4'), ('Nstar', '>i4'), 
  ('pad', '>i4')])
gas_type  = np.dtype([
  ('mass','>f4'), 
  ('r', '>f4', (3,)),  ('v', '>f4', (3,)),
  ('rho','>f4'), ('temp','>f4'), ('hsmooth','>f4'), 
  ('metals','>f4'), ('phi','>f4')])

dark_type = np.dtype([
    ('mass', '>f4'),
    ('r', '>f4', (3,)),
    ('v', '>f4', (3,)),
    ('eps', '>f4'),
    ('phi', '>f4')
])


star_type  = np.dtype([
  ('mass','>f4'), 
  ('r', '>f4', (3,)),  ('v', '>f4', (3,)),
  ('metals','>f4'), ('tform','>f4'), 
  ('eps','>f4'), ('phi','>f4')])

file = "Task01/B100.00100" 

with open(file ,'rb') as tipsy:
  header = np.fromfile(tipsy,dtype=header_type,count=1)
  header=header[0] # there is only a single row
  gas  = np.fromfile(tipsy,dtype=gas_type,count=header['Ngas'])
  dark = np.fromfile(tipsy,dtype=dark_type,count=header['Ndark'])
  star = np.fromfile(tipsy,dtype=star_type,count=header['Nstar'])

print(header)
print(dark)
print(dark['r'])
print(dark['r'][0])
print(dark['r'][:,0])

# (i) 
print("Minimum coordinates (x, y, z):", np.min(dark['r'], axis=0))
print("Maximum coordinates (x, y, z):", np.max(dark['r'], axis=0))

# (ii)
print("Total mass:", np.sum(dark['mass']))