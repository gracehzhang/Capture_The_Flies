'''
Created on Jul 22, 2015

@author: linzhang
'''
import jug
jug.init('foldx_montecarlo_jug_v2.py', 'foldx_montecarlo_jug_v2.jugdata')
import foldx_montecarlo_jug_v2
count = 0
total = 0
for i in jug.task.value(foldx_montecarlo_jug_v2.rand):
    if i >= jug.task.value(foldx_montecarlo_jug_v2.obs):
        count += 1
    total += i
r_value = count / float(len(foldx_montecarlo_jug_v2.rand))
print jug.task.value(foldx_montecarlo_jug_v2.obs)
print jug.task.value(len(foldx_montecarlo_jug_v2.rand))
print total
print total / len(foldx_montecarlo_jug_v2.rand)
print "r_value: " + str(r_value)